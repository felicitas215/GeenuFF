import os
import copy
import itertools
import bisect
import math
import logging
import shutil
from collections import Counter, defaultdict
from dataclasses import dataclass

import intervaltree
from abc import ABC, abstractmethod
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from dustdas import gffhelper, fastahelper
from geenuff.base import orm
from geenuff.base import types
from geenuff.base import helpers
from geenuff.base.helpers import (get_strand_direction, strand_or_none, get_geenuff_start_end,
                                  has_inframe_stop_codon, spliced_cds_sequence, START_CODON,
                                  STOP_CODONS)

logger = logging.getLogger(__name__)

# the GFF types of each kind, as sets, so that every line's type is looked up once
KNOWN_TYPES = frozenset(x.value for x in types.AllKnownGFFFeatures)
IGNORABLE_TYPES = frozenset(x.value for x in types.IgnorableGFFFeatures)
GENE_TYPES = frozenset(x.value for x in types.SuperLocusAll)
TRANSCRIPT_TYPES = frozenset(x.value for x in types.TranscriptLevel)
# where an exon, CDS or UTR line goes among its transcript's lines, by type
PIECE_KEYS = {x.value: key for level, key in [(types.ExonLevel, 'exons'), (types.CDSLevel, 'cds'),
                                              (types.UTRLevel, 'utrs')] for x in level}


# how much of a locus its own errors take away, lower being better (see
# GFFErrorHandling._locus_severity); used to rank candidates when an overlap is resolved
NOT_MASKED = 0       # nothing wrong with it
FLANK_MASKED = 1     # an edge is unknown, e.g. a missing UTR, but the coding sequence is correct
MASKED_OUTRIGHT = 2  # the coding sequence itself is masked, so keeping the locus recovers nothing


@dataclass(frozen=True)
class Finding:
    """One error of a coding transcript (see GFFErrorHandling._find_errors).

    extent is what the error masks: '5p' or '3p' the flank on that side of the locus, from
    handlers[0] outward; 'whole' the locus and both flanks; or an explicit (start, end) range,
    which is empty for an error that is only recorded. side is the end of every handler that is
    not where the feature really ends ('5p', '3p' or 'both'), or None."""
    error_type: str
    extent: str | tuple[int, int]
    side: str | None = None
    handlers: tuple['FeatureImporter', ...] = ()


# what each orm.SuperLocus.excluded_from_export value means, for the import summary
UNEXPORTED_REASONS = {
    types.UNPLACEABLE_STRAND: 'their transcript or its exon/CDS lines are not all on one '
                              'definite strand',
    types.UNPLACEABLE_COORDINATES: 'a line of theirs runs backwards, its start past its end',
    types.OUTSIDE_SEQUENCE: 'a line of theirs starts before or ends past their sequence',
    types.OVERLAP_DROPPED: 'dropped so an overlapping partner could be kept whole (see below)',
}

# why a GFF line was left out while the lines are grouped into loci by their ID and Parent
# attributes (see OrganizedGFFEntries.load_organized_entries), for the import summary
DROPPED_LINE_REASONS = {
    'no_id': 'gene and transcript lines without an ID, which no line can name as its parent',
    'shared_id': 'gene and transcript lines sharing their ID with another such line, so that '
                 'which one a child names is undecidable',
    'several_genes': 'transcript lines naming more than one gene as their parent',
    'parent_not_gene': 'transcript lines naming another transcript as their parent',
    'missing_parent_unplaceable': 'transcript lines naming a parent that matches no line, where the '
                                  'transcripts naming it lie on different sequences or strands',
    'no_parent': 'transcript, exon and UTR lines naming no parent (CDS lines naming none are masked)',
    'unknown_parent': 'exon and UTR lines naming a parent that matches no line (CDS lines doing so '
                      'are masked)',
    'gene_parented_duplicate': 'exon, UTR and CDS lines whose parent is a gene with transcripts, '
                               'duplicating a line of one of those transcripts',
    'gene_parented': 'exon, UTR and CDS lines whose parent is a gene with transcripts, duplicating no '
                     'line of them, and exon and UTR lines whose parent is a gene without transcripts '
                     '(CDS lines doing so are masked)',
    'floating_cds_unplaceable': 'CDS lines to be masked together that lie on different sequences or '
                                'strands, or on another sequence or strand than the gene they name',
    'parent_dropped': 'lines whose parent was left out',
    'other_sequence': 'lines on another sequence than their parent',
}


class ImportStatistics(object):
    """Aggregate counts for one add_genome() call. A single summary logged at the end of
    the import gives a full picture of what was found (see log_summary), instead of one
    log line per locus/transcript/error, which floods the log without being any easier to
    get an overview from."""

    def __init__(self) -> None:
        self.total_super_loci: int = 0
        self.coding_super_loci: int = 0  # with at least one transcript with CDS lines
        self.total_transcripts: int = 0
        self.total_coding_transcripts: int = 0
        # transcript.longest == True, one per exported gene; a gene left out of exports is not
        # counted, nothing about it having been checked either (see clean_and_insert)
        self.longest_transcripts: int = 0
        self.longest_error_free_transcripts: int = 0  # of the above, the ones with no error at all
        self.empty_super_loci: int = 0  # a gene with no transcripts at all
        # genes inferred for a Parent ID of transcripts that names no line in the file, one per
        # such ID, shared by every transcript naming it
        self.genes_inferred_for_missing_parents: int = 0
        # CDS lines without a transcript line are masked by a transcript made for them, in a gene
        # made for them unless they name a gene without transcripts (see
        # OrganizedGFFEntries._place_floating_cds)
        self.floating_cds_transcripts: int = 0
        self.floating_cds_genes: int = 0
        # transcripts with CDS lines but no exon lines, their exons built from their CDS and UTR lines
        self.transcripts_with_exons_built: int = 0
        # CDS missing only their stop codon, as in GFF3 converted from GTF, extended by the stop codon
        # after them; GTF itself is not read
        self.stop_codons_recovered: int = 0
        # GFF lines left out while grouping, keyed by a DROPPED_LINE_REASONS key and then by GFF
        # type; a line with several parents counts once per parent it is left out for
        self.dropped_lines: defaultdict[str, defaultdict[str, int]] = defaultdict(lambda: defaultdict(int))
        self.unstranded_super_loci: int = 0  # e.g. NCBI's '?' strand for trans-spliced genes
        # kept in full but never exported, counted per reason and keyed by the value stored in
        # orm.SuperLocus.excluded_from_export
        self.unexported_super_loci: defaultdict[str, int] = defaultdict(int)
        # transcripts left out entirely, no line of theirs having a storable range, or one lying
        # partly outside their sequence with no CDS line on it; a coding one is clipped instead and
        # masked whole, counted among the errors as beyond_sequence_edge
        self.unstorable_transcripts_dropped: int = 0
        self.transcripts_outside_sequence_dropped: int = 0
        # every pair of super loci sharing genomic range, of whatever kind, as recorded in the
        # super_locus_overlap table; these mask nothing, see _record_overlap_pairs
        self.overlap_pairs_recorded: int = 0
        self.super_loci_in_overlap_pairs: int = 0
        # loci whose exported transcript overlaps another's, counted once each however many partners it has
        self.super_loci_overlapping_exported: int = 0
        # of the pairs where both genes are exported, how each was settled; the four always sum
        # to that number (see GFFErrorHandling._compute_overlap_masks)
        self.overlap_pairs_resolved: int = 0  # one locus kept whole, the other dropped
        self.overlap_pairs_both_masked_outright: int = 0  # neither worth keeping, both masked whole
        self.overlap_pairs_nested: int = 0  # one inside the other, never decided, both masked whole
        self.overlap_pairs_in_chains: int = 0  # a partner overlaps a third locus, so not decided
        self.overlap_loci_dropped: int = 0
        self.overlap_loci_in_chains: int = 0  # each counted once however many pairs it is in
        # keyed by types.Errors value, counted per transcript from the error types detected, so
        # an error masking nothing still shows up here (see docs/spec_vs_gff.md)
        self.errors: defaultdict[str, int] = defaultdict(int)
        self.unrecognized_feature_types: defaultdict[str, int] = defaultdict(int)  # keyed by the raw, unknown GFF type

    @staticmethod
    def _per_type(counts: dict[str, int], text: str) -> tuple[int, str]:
        """One summary entry for a tally kept per GFF feature type: the total, with the types
        themselves in brackets, rather than the same sentence once per type."""
        if counts:
            text += ' (' + ', '.join(f'{name} {n}' for name, n in sorted(counts.items())) + ')'
        return sum(counts.values()), text

    def log_summary(self, species: str) -> None:
        """Logs one summary of the import, as sections of counted lines with the counts in a
        column, so a reader can scan down them rather than through a paragraph."""
        overlap_pairs = (self.overlap_pairs_resolved + self.overlap_pairs_both_masked_outright
                         + self.overlap_pairs_nested + self.overlap_pairs_in_chains)
        # every transcript of an undecided isolated pair is in that one pair only, two per pair
        masked_whole = (2 * self.overlap_pairs_both_masked_outright + 2 * self.overlap_pairs_nested
                        + self.overlap_loci_in_chains)
        sections = [
            ('loci', [
                (self.total_super_loci, 'genes imported, the inferred and made ones below included'),
                (self.coding_super_loci, 'of them with a coding transcript'),
                (self.empty_super_loci, 'of them with no transcript under them'),
                (self.genes_inferred_for_missing_parents, 'genes inferred for a Parent ID of '
                                                          'transcripts naming no line'),
                (self.floating_cds_genes, 'genes made to mask CDS lines naming no parent or a parent '
                                          'matching no line'),
                (self.unstranded_super_loci, "gene lines on neither strand, e.g. '.' or NCBI's '?'"),
            ] + [(count, f'genes kept but never exported: {UNEXPORTED_REASONS[reason]}')
                 for reason, count in sorted(self.unexported_super_loci.items())]),
            ('transcripts', [
                (self.total_transcripts, 'transcripts imported, the made ones below included'),
                (self.total_coding_transcripts, 'of them with a CDS'),
                (self.floating_cds_transcripts, 'made to mask CDS lines without a transcript line, '
                                                'masked whole as floating_cds'),
                (self.transcripts_with_exons_built, 'with exons built from their CDS and UTR lines, '
                                                    'having no exon lines'),
                (self.stop_codons_recovered, 'with the stop codon after their CDS added to it, the '
                                             'CDS missing only that, as in GFF3 converted from GTF '
                                             '(GTF itself is not read)'),
                (self.longest_transcripts, 'selected for export, one per locus that will be exported'),
                (self.longest_error_free_transcripts, 'of those selected with no error at all'),
                (self.unstorable_transcripts_dropped, 'dropped, no line of theirs having a '
                                                      'storable range'),
                (self.transcripts_outside_sequence_dropped, 'dropped, a line of theirs starting '
                                                            'before or ending past their sequence '
                                                            'and no CDS line on it'),
            ]),
            ('transcripts overlapping another transcript that will be exported on the same strand', [
                (self.super_loci_overlapping_exported, f'transcripts overlap at least one other, in '
                                                       f'{overlap_pairs} pairs, settled as follows:'),
                (self.overlap_pairs_resolved, '  kept whole, their partner dropped from the export'),
                (self.overlap_loci_dropped, '  dropped from the export, their region masked on the '
                                            'transcript kept'),
                (masked_whole, '  masked whole with both flanks, the sum of:'),
                (2 * self.overlap_pairs_both_masked_outright, '    in a pair where both are masked '
                                                              'outright for their own errors'),
                (2 * self.overlap_pairs_nested, '    in a pair where one is nested inside the other'
                                                'on the same strand'),
                (self.overlap_loci_in_chains, '    in a chain, overlapping more than one other transcript'),
            ]),
        ]
        if self.dropped_lines:
            sections.append(('GFF lines left out while grouping lines into genes by ID and Parent',
                             [self._per_type(self.dropped_lines[reason], text)
                              for reason, text in DROPPED_LINE_REASONS.items()
                              if reason in self.dropped_lines]))
        # the overlap masks are reported per outcome in the overlap section above instead
        errors = [(count, error_type) for error_type, count in sorted(self.errors.items())
                  if error_type != types.SL_OVERLAP_ERROR]
        if errors:
            sections.append(('errors, counted once per transcript they occur in, however often '
                             'they occur in that one', errors))
        if self.unrecognized_feature_types:
            sections.append(('lines skipped for a feature type GeenuFF has no use for',
                             [(count, feature_type) for feature_type, count
                              in sorted(self.unrecognized_feature_types.items())]))

        width = max(len(str(count)) for _, entries in sections for count, _ in entries)
        lines = [f'import summary for "{species}":']
        for heading, entries in sections:
            lines.append(f'  {heading}:')
            # leading spaces of a text indent its whole line, the count included, so that the
            # parts of a total stand below it
            lines += [f'    {" " * (len(text) - len(text.lstrip(" ")))}{count:>{width}}  {text.lstrip(" ")}'
                      for count, text in entries]
        for note in self._summary_notes():
            lines.append(f'  note: {note}')
        logger.info('\n'.join(lines))

    def _summary_notes(self) -> list[str]:
        """Explanations for the parts of the summary above that a count alone does not convey."""
        notes = [f'errors are counted as found, so one whose masked range came out empty, such as '
                 f'a missing UTR with no room to mask it in, still counts (see docs/spec_vs_gff.md)',
                 f'the super_locus_overlap table separately records {self.overlap_pairs_recorded} '
                 f'pairs over {self.super_loci_in_overlap_pairs} coding genes, counting every '
                 f'coding transcript a gene has rather than only the exported one; it masks '
                 f'nothing and is there for consumers other than Helixer']
        return notes


# core queue prep
class InsertionQueue(helpers.QueueController):
    def __init__(self, session, engine):
        super().__init__(session, engine)
        self.super_locus = helpers.CoreQueue(orm.SuperLocus.__table__.insert())
        self.transcript = helpers.CoreQueue(orm.Transcript.__table__.insert())
        self.transcript_error = helpers.CoreQueue(orm.TranscriptError.__table__.insert())
        self.transcript_piece = helpers.CoreQueue(orm.TranscriptPiece.__table__.insert())
        self.protein = helpers.CoreQueue(orm.Protein.__table__.insert())
        self.feature = helpers.CoreQueue(orm.Feature.__table__.insert())
        self.association_transcript_piece_to_feature = helpers.CoreQueue(
            orm.association_transcript_piece_to_feature.insert())
        self.association_protein_to_feature = helpers.CoreQueue(
            orm.association_protein_to_feature.insert())
        self.association_transcript_to_protein = helpers.CoreQueue(
            orm.association_transcript_to_protein.insert())
        self.super_locus_overlap = helpers.CoreQueue(orm.super_locus_overlap.insert())

        # super_locus_overlap references super_locus twice, so it follows it here
        self.ordered_queues = [
            self.super_locus, self.transcript, self.transcript_error, self.transcript_piece,
            self.protein, self.feature, self.association_transcript_piece_to_feature,
            self.association_protein_to_feature, self.association_transcript_to_protein,
            self.super_locus_overlap
        ]


class InsertCounterHolder(object):
    """provides incrementing unique integers to be used as primary keys for bulk inserts"""
    feature = helpers.Counter(orm.Feature)
    protein = helpers.Counter(orm.Protein)
    transcript = helpers.Counter(orm.Transcript)
    super_locus = helpers.Counter(orm.SuperLocus)
    transcript_piece = helpers.Counter(orm.TranscriptPiece)
    genome = helpers.Counter(orm.Genome)


class OrganizedGeenuffImporterGroup(object):
    """Stores the handler objects for a super locus in an organised fashion.
    The format is similar to the one of OrganizedGFFEntryGroup, but it stores objects
    according to the Geenuff way of saving genomic annotations. This format can then
    be checked for errors and changed accordingly before being inserted into the db.

    The importers are organised in the following way:

    importers = {
        'super_locus' = super_locus_importer,
        'transcripts': [
            {
                'transcript': transcript_importer,
                'transcript_piece: transcript_piece_importer,
                'transcript_feature': transcript_importer,
                'protein': protein_importer,
                'cds': cds_importer,
                'introns': [intron_importer1, intron_importer2, ...],
                'masks': [],  # geenuff_mask importers, filled in by GFFErrorHandling.resolve_errors
                'detected_error_types': set(),  # every error type found, stored as orm.TranscriptError
                'mask_whole': None,  # or the error type it is masked whole for, nothing else checked
                'coding_length': int,  # coding transcripts only, see _set_longest_transript
                'findings': [...],  # coding transcripts only, see GFFErrorHandling._find_errors
            },
            ...
        ],
    }
    """

    def __init__(self, organized_gff_entries, coord, controller):
        self.coord = coord
        self.controller = controller
        self.importers = {'transcripts': []}
        try:
            self._parse_gff_entries(organized_gff_entries)
        except Exception as e:
            logger.error(f'Error originally raised while parsing the following entries: '
                         f'{organized_gff_entries}')
            raise e

    def _drop_unstorable_transcript(self, sl_i, t_id, what, gff_start, gff_end):
        """Leaves a transcript out entirely because one of its lines has no range it could be
        stored as, its start running past its end. Trans-spliced organellar genes wrapping a
        circular molecule's origin are written this way; nothing can be made of them here, and
        the locus they belong to is kept out of exports."""
        sl_i.excluded_from_export = types.UNPLACEABLE_COORDINATES
        self.controller.stats.unstorable_transcripts_dropped += 1
        logger.warning(f"dropping transcript '{t_id}': its {what} runs from {gff_start} to "
                       f"{gff_end}, i.e. backwards, so it has no storable range; super locus "
                       f"'{sl_i.given_name}' will not be exported")

    @staticmethod
    def _cds_extended_by_a_codon(t_entries, is_plus_strand):
        """Copies of a transcript's CDS lines extended by the 3 bases following them along its
        exons, possibly across an intron, or None where fewer are left. Used to recover the stop
        codon GTF leaves out of the CDS, in GFF3 converted from GTF; GTF itself is not read."""
        cds = t_entries['cds']
        if is_plus_strand:
            following = (p for exon in t_entries['exons']
                         for p in range(max(exon.start, cds[-1].end + 1), exon.end + 1))
        else:
            following = (p for exon in reversed(t_entries['exons'])
                         for p in range(min(exon.end, cds[0].start - 1), exon.start - 1, -1))
        following = list(itertools.islice(following, 3))
        # fewer bases are left in the transcript, or overlapping exons repeat one
        if len(set(following)) < 3:
            return None

        pieces = [copy.copy(c) for c in cds]
        # outward from the CDS, each base extends the piece it touches or, past an intron, starts one
        for p in following:
            touching = next((c for c in pieces if (c.end + 1 if is_plus_strand else c.start - 1) == p), None)
            if touching is None:
                touching = copy.copy(cds[-1] if is_plus_strand else cds[0])
                touching.start = touching.end = p
                pieces.append(touching)
            elif is_plus_strand:
                touching.end = p
            else:
                touching.start = p
        return sorted(pieces, key=lambda c: c.start)

    def _clip_to_sequence(self, t, t_entries):
        """Clips a transcript reaching past an edge of its sequence to the part on it, leaving out
        its exon and CDS lines wholly beyond. Returns False, changing nothing, where no CDS line
        of it lies on the sequence."""
        length = self.coord.length
        on_sequence = [line for line in t_entries['cds'] if line.end >= 1 and line.start <= length]
        if not on_sequence:
            return False
        t_entries['cds'] = on_sequence
        t_entries['exons'] = [line for line in t_entries['exons'] if line.end >= 1 and line.start <= length]
        for line in [t] + t_entries['exons'] + t_entries['cds']:
            line.start, line.end = max(line.start, 1), min(line.end, length)
        return True

    def _drop_transcript_outside_sequence(self, sl_i, t_id, line):
        """Leaves a transcript with no CDS line on its sequence out entirely because one of its
        lines starts before its sequence or ends past it, which an annotation of another assembly
        version, or a gene crossing the origin of a circular molecule written with an end past
        its length, both produce. Such a line cannot be checked or labelled against the sequence,
        and one starting before it cannot be stored at all; the locus it belongs to is kept out
        of exports."""
        sl_i.excluded_from_export = types.OUTSIDE_SEQUENCE
        self.controller.stats.transcripts_outside_sequence_dropped += 1
        logger.warning(f"dropping transcript '{t_id}': its {line.type} runs from {line.start} to "
                       f"{line.end}, outside its sequence of {self.coord.length} bp; super locus "
                       f"'{sl_i.given_name}' will not be exported")

    def _parse_gff_entries(self, entries):
        """Changes the GFF format into the GeenuFF format. Does all the parsing."""
        sl = entries['super_locus']
        stats = self.controller.stats
        stats.total_super_loci += 1
        if any(t_entries['cds'] for t_entries in entries['transcripts'].values()):
            stats.coding_super_loci += 1
        try:
            sl_is_plus_strand = get_strand_direction(sl)
        except ValueError:
            # e.g. NCBI marks trans-spliced genes with strand '?': this pipeline has no way
            # to represent a gene that isn't on a single +/- strand, so the super locus is
            # saved on its own, with no transcripts/features under it. 'transcript.longest'
            # therefore stays False for everything here, which keeps it out of every
            # longest-only-filtered export query (see GeenuffExportController._genome_query)
            # without needing a dedicated error/mask type, i.e. left alone entirely
            stats.unstranded_super_loci += 1
            logger.debug(f"super locus '{sl.get_ID()}' has strand '{sl.strand}' (not '+' "
                         f"or '-'), most likely a trans-spliced gene; saving it without "
                         f"any transcripts so it is never exported or masked")
            self.importers['super_locus'] = SuperLocusImporter(entry_type=sl.type,
                                                               given_name=sl.get_ID(),
                                                               coord=self.coord,
                                                               is_plus_strand=None,
                                                               start=sl.start,
                                                               end=sl.end,
                                                               controller=self.controller)
            return

        sl_start, sl_end = get_geenuff_start_end(sl.start, sl.end, sl_is_plus_strand)
        sl_i = self.importers['super_locus'] = SuperLocusImporter(entry_type=sl.type,
                                                                  given_name=sl.get_ID(),
                                                                  coord=self.coord,
                                                                  is_plus_strand=sl_is_plus_strand,
                                                                  start=sl_start,
                                                                  end=sl_end,
                                                                  controller=self.controller)
        for t, t_entries in entries['transcripts'].items():
            stats.total_transcripts += 1
            stats.floating_cds_transcripts += t_entries['floating']
            stats.transcripts_with_exons_built += t_entries['exons_built']
            # masked whole without further checks, see GFFErrorHandling._find_errors
            t_importers = {'masks': [], 'detected_error_types': set(),
                           'mask_whole': types.FLOATING_CDS if t_entries['floating'] else None}
            t_id = t.get_ID()
            if t.start > t.end:
                self._drop_unstorable_transcript(sl_i, t_id, 'transcript', t.start, t.end)
                continue
            # GFF coordinates are 1-based and inclusive, so a line lies within its sequence from
            # 1 up to the sequence length
            outside = [line for line in [t] + t_entries['exons'] + t_entries['cds']
                       if line.start < 1 or line.end > self.coord.length]
            if outside:
                if not self._clip_to_sequence(t, t_entries):
                    self._drop_transcript_outside_sequence(sl_i, t_id, outside[0])
                    continue
                t_importers['mask_whole'] = types.BEYOND_SEQUENCE_EDGE
            t_is_plus_strand = strand_or_none(t)
            # a transcript on neither strand cannot be placed, but its features are still kept,
            # on the gene's strand, rather than dropped; where it codes, the locus is marked
            # below and left out of exports rather than masked, there being no one range to mask
            unusable_strand = t_is_plus_strand is None
            if unusable_strand:
                t_is_plus_strand = sl_is_plus_strand
                logger.debug(f"transcript '{t_id}' has strand '{t.strand}' (not '+' or '-'), most "
                             f"likely trans-spliced; keeping its features on the strand of super "
                             f"locus '{sl.get_ID()}'")
            # create transcript handler
            t_i = TranscriptImporter(entry_type=t.type,
                                     given_name=t_id,
                                     super_locus_id=sl_i.id,
                                     controller=self.controller)
            # create transcript piece handler
            tp_i = TranscriptPieceImporter(given_name=t_id,
                                           transcript_id=t_i.id,
                                           position=0,
                                           controller=self.controller)
            # create transcript feature handler
            tf_i = FeatureImporter(self.coord,
                                   t_is_plus_strand,
                                   types.GEENUFF_TRANSCRIPT,
                                   given_name=t_id,
                                   score=t.score,
                                   source=t.source,
                                   controller=self.controller)
            tf_i.set_start_end_from_gff(t.start, t.end)

            # insert everything so far into the dict
            t_importers['transcript'] = t_i
            t_importers['transcript_piece'] = tp_i
            t_importers['transcript_feature'] = tf_i

            # if it is not a non-coding gene or something like that
            if t_entries['cds']:
                stats.total_coding_transcripts += 1
                # a coding transcript that is not on one definite strand, or whose CDS lines
                # disagree with it, cannot be analysed; the locus is left out of exports while
                # its features stay in the database (see docs/trans_splicing.md for why nothing
                # is masked in its place either)
                if unusable_strand or {strand_or_none(x) for x in t_entries['cds']} != {t_is_plus_strand}:
                    sl_i.excluded_from_export = types.UNPLACEABLE_STRAND
                    logger.debug(f"coding transcript '{t_id}' is not on a single definite strand; "
                                 f"super locus '{sl.get_ID()}' will not be exported")
                # create protein handler
                # floating CDS lines masked together can name several proteins or none
                protein_id = t_id if t_entries['floating'] else self._get_protein_id_from_cds_list(t_entries['cds'])
                p_i = ProteinImporter(given_name=protein_id,
                                      super_locus_id=sl_i.id,
                                      controller=self.controller)
                # create coding features from exon limits
                phase_5p = t_entries['cds'][0 if t_is_plus_strand else -1].phase
                # spliced once here, from the raw per-piece CDS list we already have on hand,
                # to check for a truncated (non-codon-multiple) length, a premature stop codon,
                # or a missing start/stop codon, none of which the file's own annotation is
                # trusted to rule out. This also makes the start/stop codon check safe when the
                # terminal CDS piece is shorter than 3bp (the codon is split across a splice
                # junction right at the transcript's edge): reading 3 contiguous genomic bases
                # from the boundary, as the old per-feature check did, would read into the
                # intron and report a false MISSING_START/STOP_CODON in that case.
                cds_seq = spliced_cds_sequence(self.coord.sequence, t_entries['cds'], t_is_plus_strand)
                is_truncated = len(cds_seq) % 3 != 0
                has_inframe_stop = has_inframe_stop_codon(cds_seq)
                has_stop_codon = cds_seq[-3:] in STOP_CODONS
                # GTF leaves the stop codon out of the CDS, and GFF3 converted from GTF can keep it
                # that way (GTF itself is not read). Where it alone is missing, the CDS being
                # whole codons without any stop codon, and the next codon is one, it is added:
                # translation ends there whatever the file says. The other checks above hold for
                # the extended CDS as well
                if not (t_importers['mask_whole'] or is_truncated or has_inframe_stop or has_stop_codon):
                    extended = self._cds_extended_by_a_codon(t_entries, t_is_plus_strand)
                    if extended is not None:
                        extended_seq = spliced_cds_sequence(self.coord.sequence, extended, t_is_plus_strand)
                        if extended_seq[-3:] in STOP_CODONS:
                            t_entries['cds'], cds_seq, has_stop_codon = extended, extended_seq, True
                            stats.stop_codons_recovered += 1
                # the lines are sorted by start, so any overlap shows between neighbours
                cds_lines = t_entries['cds']
                overlapping_cds = any(b.start <= a.end for a, b in zip(cds_lines, cds_lines[1:]))
                cds_i = FeatureImporter(self.coord,
                                        t_is_plus_strand,
                                        types.GEENUFF_CDS,
                                        phase_5p=phase_5p,
                                        # the file's starting phase isn't trusted; the first
                                        # CDS piece always starts a fresh codon, i.e. phase 0
                                        # -> same value in both phase conventions, Helixer later
                                        # computes the phases from the starting phase
                                        phase=0,
                                        is_truncated=is_truncated,
                                        has_inframe_stop=has_inframe_stop,
                                        has_start_codon=cds_seq[:3] == START_CODON,
                                        has_stop_codon=has_stop_codon,
                                        has_overlapping_pieces=overlapping_cds,
                                        score=t.score,
                                        source=t.source,
                                        controller=self.controller)
                # the next two lines are normally enough to define cds start & end
                gff_cds_start = t_entries['cds'][0].start
                gff_cds_end = t_entries['cds'][-1].end
                # however, we have to handle partial gene models that can end in / have hanging introns
                # for a hanging intron the 'exon start' doesn't line up with the transcript start (same for ends)
                gff_exon_start = t_entries['exons'][0].start
                gff_exon_end = t_entries['exons'][-1].end
                if gff_exon_start == gff_cds_start != t.start:
                    # hanging intron at start, CDS feature will be extended to end of transcript
                    # this will be wrong if the start codon exactly aligned w/ exon start and the final exon is
                    # non-coding, but that is rare, and this implementation is more cautious / conservative
                    # (AKA: will create more error masks to reflect the ambiguity)
                    gff_cds_start = t.start
                if gff_exon_end == gff_cds_end != t.end:
                    # hanging intron at end
                    gff_cds_end = t.end

                if gff_cds_start > gff_cds_end:
                    self._drop_unstorable_transcript(sl_i, t_id, 'CDS', gff_cds_start, gff_cds_end)
                    continue
                cds_i.set_start_end_from_gff(gff_cds_start, gff_cds_end)

                # insert everything so far into the dict
                t_importers['protein'] = p_i
                t_importers['cds'] = cds_i

                # create all the introns by taking the transcript, and subtracting the exons
                # (this handles literal sequence-edge cases more reliably than the gaps between exons)
                # what remains is introns and gets FeatureImport setup and inserted into the previously created list
                introns = []
                itree = intervaltree.IntervalTree()
                etree = intervaltree.IntervalTree()
                # bc interval tree only operates w/ start < end
                inv_t_start, inv_t_end = sorted([tf_i.start, tf_i.end])
                itree[inv_t_start:inv_t_end] = t_is_plus_strand
                exons = t_entries['exons']
                for exon in exons:
                    # an exon on another strand than its transcript, or on none at all, cannot be
                    # chopped out of the span below; without its exons that span is not intronic
                    # either, so the whole locus is given up on rather than left with what would
                    # look like one intron covering the entire transcript
                    e_is_plus_strand = strand_or_none(exon)
                    if e_is_plus_strand != t_is_plus_strand:
                        sl_i.excluded_from_export = types.UNPLACEABLE_STRAND
                        itree.clear()
                        break
                    e_start, e_end = get_geenuff_start_end(exon.start, exon.end, e_is_plus_strand)
                    inv_e_start, inv_e_end = sorted([e_start, e_end])
                    itree.chop(inv_e_start, inv_e_end)
                    # also check for and mark overlapping exons (for now with a backward 'intron' # todo refactor)
                    overlapping = etree[inv_e_start:inv_e_end]
                    if overlapping:
                        overlapper_begin = min([o.begin for o in overlapping])
                        overlapper_end = max([o.end for o in overlapping])
                        if len(overlapping) != 1:
                            logger.warning('handling overlaps of >1 exon... (masking as if unioned), but this is a '
                                           'weird enough sort of error that you should really check what is going on '
                                           'if you read this (around {} {}-{})'.format(self.coord, overlapper_begin,
                                                                                       overlapper_end))

                        ovlp_end = max(overlapper_begin, inv_e_start)
                        ovlp_start = min(overlapper_end, inv_e_end)
                        logger.debug('found an overlap {}-{}'.format(ovlp_end, ovlp_start))
                        if not e_is_plus_strand:
                            ovlp_start, ovlp_end = ovlp_end, ovlp_start
                        # insert a dummy 'backwards' intron, which will later be turned into an overlap error
                        intron_err = FeatureImporter(self.coord,
                                                     is_plus_strand=e_is_plus_strand,
                                                     feature_type=types.GEENUFF_INTRON,
                                                     start=ovlp_start,
                                                     end=ovlp_end,
                                                     score=t.score,
                                                     source=t.source,
                                                     controller=self.controller)
                        introns.append(intron_err)
                    etree[inv_e_start:inv_e_end] = e_is_plus_strand

                for interval in itree:
                    i_start, i_end = interval.begin, interval.end
                    if not interval.data:  # if minus strand (data = is_plus_strand)
                        i_start, i_end = i_end, i_start
                    intron_i = FeatureImporter(self.coord,
                                               is_plus_strand=interval.data,
                                               feature_type=types.GEENUFF_INTRON,
                                               start=i_start,
                                               end=i_end,
                                               score=t.score,
                                               source=t.source,
                                               controller=self.controller)
                    introns.append(intron_i)

                introns = sorted(introns, key=lambda x: x.start)
                if t_is_plus_strand:
                    t_importers['introns'] = introns
                else:
                    t_importers['introns'] = introns[::-1]
            self.importers['transcripts'].append(t_importers)
        self._set_longest_transript()

    def _set_longest_transript(self):
        """Looks for the transcript with the longest spliced CDS (see _coding_length) and sets the
        'longest' parameters in all the TranscriptImporters. In case of a tie, the first
        transcript found will be set as longest."""
        max_exon_len = -1
        longest_importer = None
        for t in self.importers['transcripts']:
            if 'cds' in t:
                # kept to rank the locus against an overlapping one by the same measure (see
                # GFFErrorHandling._decide_overlap_pair)
                t['coding_length'] = self._coding_length(t['cds'], t['introns'])
                if t['coding_length'] > max_exon_len:
                    max_exon_len = t['coding_length']
                    longest_importer = t['transcript']
        for t in self.importers['transcripts']:
            if t['transcript'] is longest_importer:
                t['transcript'].longest = True
            else:
                t['transcript'].longest = False

    @staticmethod
    def _coding_length(cds, introns):
        """Spliced length of a CDS: its span less the parts of the transcript's introns within it."""
        c_lo, c_hi = sorted((cds.start, cds.end))
        length = c_hi - c_lo
        for intron in introns:
            lo, hi = sorted((intron.start, intron.end))
            length -= max(0, min(c_hi, hi) - max(c_lo, lo))
        return length

    @staticmethod
    def _get_protein_id_from_cds_entry(cds_entry):
        # check if anything is labeled as protein_id
        protein_id = cds_entry.attrib_filter(tag='protein_id')
        # failing that, try and get parent ID (presumably transcript, maybe gene)
        if not protein_id:
            protein_id = cds_entry.get_Parent()
        # hopefully take single hit
        if len(protein_id) == 1:
            protein_id = protein_id[0]
            if isinstance(protein_id, gffhelper.GFFAttribute):
                protein_id = protein_id.value
                assert len(protein_id) == 1
                protein_id = protein_id[0]
        # or handle other cases
        elif len(protein_id) == 0:
            protein_id = None
        else:
            raise ValueError('indeterminate single protein id {}'.format(protein_id))
        return protein_id

    @staticmethod
    def _get_protein_id_from_cds_list(cds_entry_list):
        """Returns the protein id of a list of cds gff entries. If multiple ids or no id at all
        are found, an error is raised."""
        protein_ids = set()
        for cds_entry in cds_entry_list:
            protein_id = OrganizedGeenuffImporterGroup._get_protein_id_from_cds_entry(cds_entry)
            if protein_id:
                protein_ids.add(protein_id)
        if len(protein_ids) != 1:
            raise ValueError('No protein_id or more than one protein_ids for one transcript')
        return protein_ids.pop()


class OrganizedGFFEntryGroup(object):
    """Holds the entries of one super locus as grouped by OrganizedGFFEntries, and returns the
    corresponding OrganizedGeenuffImporterGroup. Does not perform error checking, which happens
    later.

    The entries are organised in the following way, exons, CDS and UTRs sorted by start:

    entries = {
        'super_locus' = super_locus_entry,
        'transcripts' = {
            transcript_entry1: {
                'exons': [ordered_exon_entry1, ordered_exon_entry2, ...],
                'cds': [ordered_cds_entry1, ordered_cds_entry2, ...],
                'utrs': [ordered_utr_entry1, ordered_utr_entry2, ...],
                'exons_built': False,  # True where the exons were built from CDS and UTR lines
                'floating': False,  # True for one made only to mask CDS lines without a transcript line
            },
            ...
        }
    }
    """

    def __init__(self, entries, fasta_importer, controller):
        self.controller = controller
        self.entries = entries
        self.coord = fasta_importer.gffid_to_coords[entries['super_locus'].seqid]

    def get_geenuff_importers(self):
        geenuff_importer_group = OrganizedGeenuffImporterGroup(self.entries, self.coord,
                                                               self.controller)
        return geenuff_importer_group.importers


class OrganizedGFFEntries(object):
    """Groups the gff entries coming from gffhelper into genes by their ID and Parent attributes,
    whatever the order of the lines in the file. Also does some basic gff value cleanup.
    The entries are organized in the following way, one entry group per gene (see
    OrganizedGFFEntryGroup):

    organized_entries = {
        'seqid1': [entry_group_gene1, entry_group_gene2, ...],
        'seqid2': [...],
        ...
    }

    Lines that cannot be placed are left out and counted (see DROPPED_LINE_REASONS), together
    with every line below them.
    """

    def __init__(self, gff_file, stats=None):
        self.gff_file = gff_file
        self.organized_entries = {}
        # only present when constructed as part of a full ImportController.add_gff() run;
        # standalone/test use (no controller) gets its own throwaway counter
        self.stats = stats if stats is not None else ImportStatistics()

    def load_organized_entries(self):
        genes, transcripts, pieces = [], [], []
        for entry in self._useful_gff_entries():
            if entry.type in GENE_TYPES:
                genes.append(entry)
            elif entry.type in TRANSCRIPT_TYPES:
                transcripts.append(entry)
            elif entry.type in PIECE_KEYS:
                pieces.append((entry, PIECE_KEYS[entry.type]))
            else:
                logger.debug(f'ignoring {entry.type} at {entry.seqid}:{entry.start}-{entry.end}')

        # the IDs of every gene or transcript line left out, so that the lines below it go too
        dropped_ids = set()
        genes_by_id = self._lines_by_unique_id(genes, transcripts, dropped_ids)
        transcripts_by_id = self._lines_by_unique_id(transcripts, genes, dropped_ids)
        groups = {gene_id: {'super_locus': gene, 'transcripts': {}}
                  for gene_id, gene in genes_by_id.items()}
        self._place_transcripts(transcripts_by_id, groups, dropped_ids)
        self._place_pieces(pieces, groups, dropped_ids)

        self.organized_entries = {}
        for group in groups.values():
            for t, t_entries in group['transcripts'].items():
                for key in ['exons', 'cds', 'utrs']:
                    t_entries[key].sort(key=lambda e: e.start)
                if t_entries['cds'] and not t_entries['exons']:
                    t_entries['exons'] = self._built_exons(t, t_entries)
                    t_entries['exons_built'] = True
            self.organized_entries.setdefault(group['super_locus'].seqid, []).append(group)

    @staticmethod
    def _new_t_entries(floating=False):
        """The lines of one transcript; floating is True for one made only to mask CDS lines
        without a transcript line, see _place_floating_cds."""
        return {'exons': [], 'cds': [], 'utrs': [], 'exons_built': False, 'floating': floating}

    def _line(self, seqid, feature_type, start, end, strand, attributes):
        """A GFF line made by GeenuFF itself, for a gene, transcript or exon the file lacks."""
        line = gffhelper.GFFObject(f'{seqid}\tGeenuFF\t{feature_type}\t{start}\t{end}\t.\t'
                                   f'{strand or "."}\t.\t{attributes}')
        self._clean_entry(line)
        return line

    def _built_exons(self, t, t_entries):
        """Exons for a transcript with CDS lines but no exon lines, built from its CDS and UTR lines:
        pieces that touch, such as a UTR ending right before the CDS starts, form one exon.
        Overlapping pieces stay apart, so that the overlap is found as an error."""
        exons = []
        for piece in sorted(t_entries['cds'] + t_entries['utrs'], key=lambda e: e.start):
            if exons and piece.start == exons[-1].end + 1:
                exons[-1].end = piece.end
            else:
                exons.append(self._line(piece.seqid, 'exon', piece.start, piece.end, piece.strand,
                                        f'Parent={t.get_ID()}'))
        return exons

    @staticmethod
    def _id_of(entry):
        """The line's ID, or None where it has none, for which gffhelper's get_ID raises an error."""
        return next((attribute.value[0] for attribute in entry.attributes if attribute.tag == 'ID'),
                    None)

    def _drop(self, entry, reason):
        self.stats.dropped_lines[reason][entry.type] += 1
        logger.debug(f'leaving out {entry.type} at {entry.seqid}:{entry.start}-{entry.end}: '
                     f'{DROPPED_LINE_REASONS[reason]}')

    def _lines_by_unique_id(self, lines, other_parent_lines, dropped_ids):
        """{ID: line} of the gene or transcript lines whose ID no other gene or transcript line
        uses. A line without an ID, which no line can name as its parent, or with an ID another
        line shares, which makes what a child names undecidable, is left out."""
        counts = Counter(self._id_of(line) for line in lines + other_parent_lines)
        by_id = {}
        for line in lines:
            line_id = self._id_of(line)
            if line_id is None:
                self._drop(line, 'no_id')
            elif counts[line_id] > 1:
                self._drop(line, 'shared_id')
                dropped_ids.add(line_id)
            else:
                by_id[line_id] = line
        return by_id

    def _place_transcripts(self, transcripts_by_id, groups, dropped_ids):
        """Puts every transcript under the one gene it names. The transcripts naming one Parent ID
        that matches no line share a gene inferred for them, spanning them all (see _inferred_gene),
        unless they lie on different sequences or strands, which leaves them all out."""
        named_missing_parent = defaultdict(list)
        for t_id, t in transcripts_by_id.items():
            parents = t.get_Parent() or []
            if len(parents) == 1 and parents[0] in groups:
                reason = 'other_sequence' if t.seqid != groups[parents[0]]['super_locus'].seqid else None
            elif not parents:
                reason = 'no_parent'
            elif len(parents) > 1:
                reason = 'several_genes'
            elif parents[0] in dropped_ids:
                reason = 'parent_dropped'
            elif parents[0] in transcripts_by_id:
                reason = 'parent_not_gene'
            else:
                named_missing_parent[parents[0]].append(t)
                continue
            if reason is None:
                groups[parents[0]]['transcripts'][t] = self._new_t_entries()
            else:
                self._drop(t, reason)
                dropped_ids.add(t_id)

        for parent_id, members in named_missing_parent.items():
            if len({(t.seqid, t.strand) for t in members}) > 1:
                for t in members:
                    self._drop(t, 'missing_parent_unplaceable')
                    dropped_ids.add(t.get_ID())
                continue
            groups[parent_id] = {'super_locus': self._inferred_gene(parent_id, members),
                                 'transcripts': {t: self._new_t_entries() for t in members}}
            self.stats.genes_inferred_for_missing_parents += 1

    def _inferred_gene(self, gene_id, transcripts):
        """A gene line for transcripts on one sequence and strand that the file gives none,
        spanning them all."""
        return self._line(transcripts[0].seqid, 'gene', min(t.start for t in transcripts),
                          max(t.end for t in transcripts), transcripts[0].strand, f'ID={gene_id}')

    def _place_pieces(self, pieces, groups, dropped_ids):
        """Puts every exon, CDS and UTR line under each transcript it names, several being allowed.
        An exon or UTR line naming a gene is left out, whether or not it duplicates a line of one
        of that gene's transcripts: it could equally belong to an isoform the file gives no
        transcript line of its own. CDS lines without a transcript line are masked instead (see
        _place_floating_cds), unless they name a gene that has transcripts."""
        transcripts = {t.get_ID(): (t, t_entries) for group in groups.values()
                       for t, t_entries in group['transcripts'].items()}
        gene_parented, floating = [], defaultdict(list)
        for piece, key in pieces:
            parents = piece.get_Parent() or []
            if not parents:
                if key == 'cds':
                    name = self._id_of(piece) or f'{piece.seqid}:{piece.start}-{piece.end}'
                    floating[('no_parent', name)].append(piece)
                else:
                    self._drop(piece, 'no_parent')
            for parent_id in parents:
                if parent_id in transcripts:
                    t, t_entries = transcripts[parent_id]
                    if piece.seqid != t.seqid:
                        self._drop(piece, 'other_sequence')
                    else:
                        t_entries[key].append(piece)
                elif parent_id in groups:
                    if key == 'cds' and not groups[parent_id]['transcripts']:
                        floating[('gene', parent_id)].append(piece)
                    else:
                        gene_parented.append((piece, key, parent_id))
                elif parent_id in dropped_ids:
                    self._drop(piece, 'parent_dropped')
                elif key == 'cds':
                    floating[('unknown_parent', parent_id)].append(piece)
                else:
                    self._drop(piece, 'unknown_parent')

        # judged only once every transcript line is placed, so the order of lines does not matter
        for piece, key, gene_id in gene_parented:
            position = (piece.start, piece.end)
            duplicate = any(position == (other.start, other.end)
                            for t_entries in groups[gene_id]['transcripts'].values()
                            for other in t_entries[key])
            self._drop(piece, 'gene_parented_duplicate' if duplicate else 'gene_parented')
        self._place_floating_cds(floating, groups)

    def _place_floating_cds(self, floating, groups):
        """CDS lines without a transcript line: naming a gene without transcripts or a parent that
        matches no line, grouped by that parent's ID; or naming no parent, grouped only where they
        share an identical ID (e.g. NCBI), each on its own otherwise. What they code is unknown,
        but some gene was annotated there, so each group
        gets a transcript spanning it that is masked whole (see GFFErrorHandling._find_errors),
        under the gene it names or a gene made for it. A group is left out if its lines lie on
        different sequences or strands, or on another sequence or strand than the gene they name."""
        for (kind, name), lines in floating.items():
            gene = groups[name]['super_locus'] if kind == 'gene' else None
            places = {(line.seqid, line.strand) for line in lines + ([gene] if gene else [])}
            if len(places) > 1:
                for line in lines:
                    self._drop(line, 'floating_cds_unplaceable')
                continue
            gene_id = name if gene else f'floating_gene_{name}'
            t = self._line(lines[0].seqid, 'mRNA', min(line.start for line in lines),
                           max(line.end for line in lines), lines[0].strand,
                           f'ID=floating_{name};Parent={gene_id}')
            if gene is None:
                groups[gene_id] = {'super_locus': self._inferred_gene(gene_id, [t]), 'transcripts': {}}
                self.stats.floating_cds_genes += 1
            t_entries = groups[gene_id]['transcripts'][t] = self._new_t_entries(floating=True)
            t_entries['cds'] = lines

    def _useful_gff_entries(self):
        reader = self._gff_gen()
        for entry in reader:
            if entry.type not in IGNORABLE_TYPES:
                yield entry

    def _gff_gen(self):
        reader = gffhelper.read_gff_file(self.gff_file)
        for entry in reader:
            if entry.type not in KNOWN_TYPES:
                # an unrecognized Sequence Ontology term shouldn't take down the whole
                # import; skip just this line and keep going
                self.stats.unrecognized_feature_types[entry.type] += 1
                # debug -> leads to overprinting in some files if set to warning
                logger.debug(f'skipping line with unrecognized feature type '
                             f'"{entry.type}" at {entry.seqid}:{entry.start}-{entry.end} '
                             f'(not a Sequence Ontology term GeenuFF knows how to handle)')
                continue
            self._clean_entry(entry)
            yield entry

    @staticmethod
    def _clean_entry(entry):
        # always present and integers
        entry.start = int(entry.start)
        entry.end = int(entry.end)
        # clean up score
        if entry.score == '.':
            entry.score = None
        else:
            entry.score = float(entry.score)

        # clean up phase
        if entry.phase == '.':
            entry.phase = None
        else:
            entry.phase = int(entry.phase)
        assert entry.phase in [None, 0, 1, 2]

        # clean up strand
        if entry.strand == '.' or entry.strand == '?':  # new ? for test purposes
            entry.strand = None
        else:
            assert entry.strand in ['+', '-']


class GFFErrorHandling(object):
    """Deals with error detection and handling of the input features. Does the handling
    in the space of GeenuFF importers.
    Assumes all super locus handler groups to be ordered 5p to 3p and of one strand.
    Works with a list of OrganizedGeenuffImporterGroup, which correspond to a list of
    super loci, and looks for errors, recording them per transcript and merging what they mask
    into its geenuff_mask features (see resolve_errors).
    """

    def __init__(self, geenuff_importer_groups, controller):
        self.groups = geenuff_importer_groups
        self.controller = controller
        self._overlap_masks = []
        self._extents = {}
        self._starts, self._ends = [], []
        if self.groups:
            self.is_plus_strand = self.groups[0]['super_locus'].is_plus_strand
            self.coord = self.groups[0]['super_locus'].coord
            # 5p to 3p by the selected transcript, the gene line of a locus without one
            self.groups.sort(key=self._sort_start, reverse=not self.is_plus_strand)
            # found before overlaps are settled, which grade each locus by its errors
            for group in self.groups:
                for transcript in group['transcripts']:
                    if 'cds' in transcript:
                        transcript['findings'] = self._find_errors(transcript)
            self._compute_overlap_masks()

    @staticmethod
    def _sort_start(group):
        transcript = GFFErrorHandling._selected_transcript(group)
        if transcript is None:
            return group['super_locus'].start
        return transcript['transcript_feature'].start

    def _compute_overlap_masks(self):
        """Settles what happens to every pair of super loci on this strand that share genomic
        range, in two records serving two different consumers.

        The unfiltered record goes to the super_locus_overlap table: one row per pair, for every
        kind of locus, including ones with no transcript at all. It is the annotation's geometry
        as given, it masks nothing, and no consumer that reads geenuff_mask features sees it.

        The second record decides what an export is handed. Two loci cannot share a base: one
        label per base is written, so whichever is written last silently wins. Of an isolated
        crossing pair the better locus is kept whole and its partner is dropped from the export,
        leaving one coherent gene rather than two with a hole through them (see
        _decide_overlap_pair). Where both are masked outright, or nested inside the other on the
        same strand, both genes are masked over their whole length and the flank on both sides.
        """
        # multiplying by sign turns both strands into ascending coordinates running 5p to 3p, so
        # everything below can compare and sort without asking which strand it is on
        sign = 1 if self.is_plus_strand else -1
        self._record_overlap_pairs(sign)

        # the span each locus would occupy in an export, keyed by its index in self.groups. A
        # locus with nothing to export, non-coding or already excluded, is simply left out, which
        # is what keeps it from being anyone's overlap partner further down
        extents = {}
        for i, group in enumerate(self.groups):
            transcript = self._exported_transcript(group)
            if transcript is not None:
                tf = transcript['transcript_feature']
                extents[i] = (sign * tf.start, sign * tf.end)
        # every mask measures its flank over these spans, toward the closest of their edges (see
        # _buffered_span); a locus dropped below stays among them, being an annotated gene all the
        # same, which a gene's true end is not assumed to run through
        self._extents = extents
        self._starts = sorted(lo for lo, _ in extents.values())
        self._ends = sorted(hi for _, hi in extents.values())
        # graded once per locus, however many pairs it is in
        severity = {i: self._locus_severity(self.groups[i]) for i in extents}

        # _overlapping_pairs wants (lo, hi, key) triples, hence appending the index to each span.
        # pairs is used as an ordered set, the shared range it could hold being recomputed in
        # _decide_overlap_pair from the extents; partners counts how many other loci each one
        # runs into, which is what distinguishes an isolated pair from a chain
        pairs, partners = {}, defaultdict(set)
        for i, j, _, _ in self._overlapping_pairs([e + (i,) for i, e in extents.items()]):
            pairs[(min(i, j), max(i, j))] = None
            partners[i].add(j)
            partners[j].add(i)

        stats = self.controller.stats
        # every locus with at least one partner, counted once however many partners it has
        stats.super_loci_overlapping_exported += len(partners)
        # ranges to mask, gathered per locus and merged at the end; a locus in several pairs
        # collects an entry from each of them
        masks = [[] for _ in self.groups]
        # what each decision came to, collected here because the masks it implies cannot be
        # worked out until every pair has been decided (see the second loop)
        # a locus in two undecided pairs lands in a set twice over and is masked once
        resolved, masked_whole, in_chains = [], set(), set()
        for i, j in pairs:
            # only an isolated pair is decided; in a chain a dropped locus can sit between two
            # kept ones, and the one masked region recorded per locus cannot express that. Where
            # nothing is decided, both genes are masked end to end
            if len(partners[i]) > 1 or len(partners[j]) > 1:
                in_chains.update((i, j))
                stats.overlap_pairs_in_chains += 1
            elif self._is_nested(extents[i], extents[j]):
                masked_whole.update((i, j))
                stats.overlap_pairs_nested += 1
            else:
                keeper, dropped, pieces = self._decide_overlap_pair(i, j, extents, severity)
                if keeper is None:
                    masked_whole.update((i, j))
                    stats.overlap_pairs_both_masked_outright += 1
                else:
                    self.groups[dropped]['super_locus'].excluded_from_export = types.OVERLAP_DROPPED
                    resolved.append((keeper, dropped, pieces))
                    stats.overlap_pairs_resolved += 1
                    stats.overlap_loci_dropped += 1
        stats.overlap_loci_in_chains += len(in_chains)
        masked_whole |= in_chains

        # an unresolved locus is masked like any gene with an error masking it whole, its flank
        # included on both sides
        for i in masked_whole:
            masks[i].append(self._buffered_span(i))
        for keeper, dropped, pieces in resolved:
            # an erroneous dropped locus has untrustworthy outer boundaries like any other erroneous
            # gene, so the piece sticking out past the kept locus gets the flank on its outer side
            erroneous = severity[dropped] != NOT_MASKED
            buffered_lo, buffered_hi = self._buffered_span(dropped)
            for (lo, hi), side in pieces:
                if erroneous and side == -1:
                    lo = buffered_lo
                elif erroneous and side == 1:
                    hi = buffered_hi
                # the mask goes on the locus that is kept, the dropped one's own features never
                # reaching an export to carry it (see orm.SuperLocus.excluded_from_export)
                masks[keeper].append((lo, hi))
        # kept sign-normalised, to be merged with the ranges of the locus' own errors (see
        # resolve_errors)
        self._overlap_masks = masks

    def _buffered_span(self, i):
        """The span of locus i's selected transcript with the flank on both sides, sign-normalised.

        Each flank reaches into the gap toward the closest edge of another locus' selected
        transcript, i.e. the largest end at or before the span's start and the smallest start at
        or after its end, or toward the end of the sequence where there is none. Taking the
        closest edge rather than the next locus in sort order matters where a long transcript
        encloses a shorter one: the shorter one need not end closest to whatever follows the long
        one. A transcript overlapping locus i is not a neighbour, sharing sequence with it rather
        than bounding it, so the flank is measured past it."""
        lo, hi = self._extents[i]
        lower, upper = self._sequence_bounds()
        j = bisect.bisect_right(self._ends, lo)
        previous_end = self._ends[j - 1] if j > 0 else lower
        j = bisect.bisect_left(self._starts, hi)
        next_start = self._starts[j] if j < len(self._starts) else upper
        return lo - self._flank_reach(lo - previous_end), hi + self._flank_reach(next_start - hi)

    @staticmethod
    def _flank_reach(gap):
        """How far a flank reaches into a gap, min(gap // 2, int(sqrt(gap)) * 10). Nothing at all
        where there is no gap, as between touching loci."""
        if gap <= 0:
            return 0
        return min(gap // 2, int(math.sqrt(gap)) * 10)

    @staticmethod
    def _is_nested(span_a, span_b):
        """Whether one span lies wholly inside the other, which on one strand is more often an
        annotation mistake than two real genes, so such a pair is never decided. Identical spans
        are not nested, there being no inner locus."""
        (a_lo, a_hi), (b_lo, b_hi) = span_a, span_b
        if (a_lo, a_hi) == (b_lo, b_hi):
            return False
        return a_lo <= b_lo and b_hi <= a_hi or b_lo <= a_lo and a_hi <= b_hi

    def _decide_overlap_pair(self, i, j, extents, severity):
        """Works out which locus of a crossing pair is kept, and what is masked in place of the
        other. Returns (keeper, dropped, pieces) or three Nones where neither is kept.

        The keeper is the less damaged locus by severity ({index: _locus_severity}), on a tie the
        one with the longer spliced CDS, the measure that picks the longest isoform, and on a tie
        of both the 5'-most one, whichever of the two has coding sequence in the range they share.
        A locus masked outright for its own errors is never kept, keeping it recovering nothing.

        The pieces are (range, side) pairs: what the dropped locus covers beyond the kept one,
        side being the direction the piece points away from the keeper, past which it may run on
        into the flank (see _buffered_span). Identical spans need no mask, the kept locus covering
        the same bases.
        """
        options = [k for k in (i, j) if severity[k] != MASKED_OUTRIGHT]
        if not options:
            return None, None, None
        def rank(k):
            # groups are ordered 5p to 3p, so the lower index is the 5'-most locus
            return severity[k], -self._selected_transcript(self.groups[k])['coding_length'], k

        keeper = min(options, key=rank)
        dropped = j if keeper == i else i
        (keeper_lo, keeper_hi), (dropped_lo, dropped_hi) = extents[keeper], extents[dropped]
        pieces = []
        if dropped_lo < keeper_lo:
            pieces.append(((dropped_lo, keeper_lo), -1))
        if dropped_hi > keeper_hi:
            pieces.append(((keeper_hi, dropped_hi), 1))
        return keeper, dropped, pieces

    @staticmethod
    def _locus_severity(group):
        """How much of a locus its own errors already take away, lower being better. Keeping a
        locus whose coding sequence is masked anyway recovers nothing, while one whose errors
        only mask a flank still has its coding sequence labelled (see _find_errors)."""
        extents = [f.extent for f in GFFErrorHandling._selected_transcript(group)['findings']]
        if any(e == 'whole' or isinstance(e, tuple) and e[0] != e[1] for e in extents):
            return MASKED_OUTRIGHT
        if any(e in ('5p', '3p') for e in extents):
            return FLANK_MASKED
        return NOT_MASKED

    @staticmethod
    def _coding_extent(group, sign):
        """Everything a locus' coding transcripts reach over, as a sign-normalised ascending
        range. Taken from the transcripts rather than the gene line, which can be far wider or
        narrower than the transcripts hanging off it."""
        bounds = [(sign * t['transcript_feature'].start, sign * t['transcript_feature'].end)
                  for t in group['transcripts']
                  if 'cds' in t and t.get('transcript_feature') is not None]
        return min(lo for lo, _ in bounds), max(hi for _, hi in bounds)

    def _record_overlap_pairs(self, sign):
        """Queues one super_locus_overlap row per pair of coding super loci sharing genomic range.
        Measured over every coding transcript a locus has, not just the exported one, so this
        record covers every overlap the masking sweep below can find and more.

        Loci with no CDS anywhere under them take no part: a bare gene line whose exons name the
        gene as their parent, or a purely non-coding gene, is not a gene anything here can be said
        about, and such records outnumber the real ones several times over."""
        spans = [self._coding_extent(g, sign) + (i,) for i, g in enumerate(self.groups)
                 if any('cds' in t for t in g['transcripts'])]
        involved = set()
        for i, j, lo, hi in self._overlapping_pairs(spans):
            self.controller.insertion_queues.super_locus_overlap.queue.append({
                'super_locus_id': self.groups[i]['super_locus'].id,
                'partner_id': self.groups[j]['super_locus'].id,
                'start': sign * lo,
                'end': sign * hi,
                'is_plus_strand': self.is_plus_strand,
                'coordinate_id': self.coord.id,
            })
            self.controller.stats.overlap_pairs_recorded += 1
            involved.update((i, j))
        self.controller.stats.super_loci_in_overlap_pairs += len(involved)

    @staticmethod
    def _selected_transcript(group):
        """The locus' longest coding transcript, whether or not the locus is exported, or None
        for a locus with no coding transcript at all."""
        for transcript in group['transcripts']:
            if 'cds' in transcript and transcript['transcript'].longest:
                return transcript
        return None

    @staticmethod
    def _exported_transcript(group):
        """The transcript exported for this super locus, or None where it is not exported. A
        dropped locus still has a selected transcript, it just reaches no consumer."""
        if group['super_locus'].excluded_from_export is not None:
            return None
        return GFFErrorHandling._selected_transcript(group)

    @staticmethod
    def _overlapping_pairs(spans):
        """Sweeps (lo, hi, key) spans, yielding (key_a, key_b, lo, hi) for every pair sharing
        range, with (lo, hi) the shared range itself. Keeping every span that still reaches
        past the current start ('active') is what finds every partner: comparing each span only
        against its immediate predecessor in sort order misses both the partner that lies
        downstream (a span overlapped only from the 3p side is never the second argument of such
        a comparison, so it is never examined at all) and the non-adjacent partner, e.g.
        A(1-1000), B(50-100) nested in A, C(500-600) also inside A but not touching B, where
        comparing C only to B finds nothing."""
        active = []
        for lo, hi, key in sorted(spans):
            # drop the spans that end at or before this one starts; directly adjacent is fine
            active = [(a_hi, a_key) for a_hi, a_key in active if a_hi > lo]
            for a_hi, a_key in active:
                # spans are swept in ascending start order, so the shared range always starts
                # at the current span's own start
                yield a_key, key, lo, min(hi, a_hi)
            active.append((hi, key))

    @staticmethod
    def _merge_ranges(ranges):
        """Merges ascending, half-open (lo, hi) ranges into a minimal, sorted list of ranges."""
        merged = []
        for lo, hi in sorted(ranges):
            if merged and lo <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], hi)
            else:
                merged.append([lo, hi])
        return [tuple(r) for r in merged]

    def _find_errors(self, transcript):
        """Every error of a coding transcript, as a list of Finding."""
        cds, tf, introns = transcript['cds'], transcript['transcript_feature'], transcript['introns']
        # made only to mask CDS lines without a transcript line, or cut off by its sequence's edge:
        # nothing else about it is known
        if transcript['mask_whole']:
            return [Finding(transcript['mask_whole'], 'whole', 'both', (cds, tf))]
        proper = [x for x in introns if not self._is_backwards(x)]
        # overlapping exon lines (a backwards intron standing in for them) or overlapping CDS lines
        # most likely come from a wrong annotation, so the locus is masked whole together with
        # the flank on both sides. Nothing else is checked: the spliced CDS sequence the checks
        # below rely on repeats the overlapping bases and is wrong
        overlaps = []
        if len(proper) < len(introns):
            overlaps.append(types.OVERLAPPING_EXONS)
        if cds.has_overlapping_pieces:
            overlaps.append(types.OVERLAPPING_CDS)
        if overlaps:
            return [Finding(error_type, 'whole', 'both', (cds, tf)) for error_type in overlaps]

        found = []
        # a missing UTR leaves only where the transcript ends unknown, the CDS being sound
        if cds.start == tf.start:
            found.append(Finding(types.MISSING_UTR_5P, '5p', '5p', (cds, tf)))
        if cds.end == tf.end:
            found.append(Finding(types.MISSING_UTR_3P, '3p', '3p', (cds, tf)))

        # every error below means the CDS boundaries cannot be trusted, so the locus is masked
        # whole together with the flank on both sides. Codons and the reading frame are checked
        # once in _parse_gff_entries, on the spliced CDS sequence (see there for why)
        if not cds.has_start_codon:
            found.append(Finding(types.MISSING_START_CODON, 'whole', '5p', (cds,)))
        if not cds.has_stop_codon:
            found.append(Finding(types.MISSING_STOP_CODON, 'whole', '3p', (cds,)))
        if cds.is_truncated:  # not a multiple of 3
            found.append(Finding(types.TRUNCATED_CDS, 'whole'))
        if cds.has_inframe_stop:
            found.append(Finding(types.INFRAME_STOP_CODON, 'whole'))
        # cannot be spliced; an intron the transcript starts or ends in is only partly inside it,
        # so its length says nothing, and it is a truncated intron instead (see below)
        min_length = self.controller.config['min_intron_length']
        complete = [x for x in proper if x.start != tf.start and x.end != tf.end]
        if any(abs(x.end - x.start) < min_length for x in complete):
            found.append(Finding(types.TOO_SHORT_INTRON, 'whole'))
        # the transcript starting or ending inside an intron, its outermost exon missing, so where
        # the transcript really starts or ends is unknown as well
        for intron in proper:
            if intron.start == tf.start:
                found.append(Finding(types.TRUNCATED_INTRON, 'whole', '5p', (intron, tf)))
            if intron.end == tf.end:
                found.append(Finding(types.TRUNCATED_INTRON, 'whole', '3p', (intron, tf)))

        # a starting phase in the file other than 0 is only recorded, the importer setting the
        # phase itself whatever the file says
        if cds.phase_5p != 0:
            found.append(Finding(types.WRONG_PHASE_5P, (cds.start, cds.start)))
        return found

    def _is_backwards(self, feature):
        """Whether a feature runs against the strand, as the stand-in for overlapping exons does."""
        if self.is_plus_strand:
            return feature.end < feature.start
        return feature.end > feature.start  # GeenuFF convention: - strand start > end

    def resolve_errors(self):
        """For each coding transcript of an exported locus: marks the ends its errors make unknown,
        records which errors it has, and merges every range they mask, together with the ranges
        given up to overlapping loci, into its geenuff_mask features."""
        for i, group in enumerate(self.groups):
            sign = 1 if self.is_plus_strand else -1
            reason = group['super_locus'].excluded_from_export
            # a locus dropped for an overlap was checked like any other, so its errors are recorded,
            # but it gets no mask, its features reaching no export, and no count in the import
            # statistics (see clean_and_insert)
            if reason == types.OVERLAP_DROPPED:
                for transcript in group['transcripts']:
                    for finding in transcript.get('findings', []):
                        transcript['detected_error_types'].add(finding.error_type)
                continue
            # a locus excluded for any other reason is not analysed at all, and gets no mask: where
            # its features belong is exactly what is unknown about it, so any mask would be
            # guesswork, either covering the wrong strand or covering sequence that is fine
            if reason is not None:
                continue

            # the case of no transcript for a super locus
            if not group['transcripts']:
                self.controller.stats.empty_super_loci += 1
                logger.debug('{} is a gene without any transcripts; it has no features and will '
                             'never be exported or masked'.format(group['super_locus'].given_name))
            # the selected transcript with the flank on both sides, which every mask reaches to
            buffered = self._buffered_span(i) if i in self._extents else None
            for transcript in group['transcripts']:
                if 'cds' not in transcript:
                    continue
                ranges = []
                for finding in transcript['findings']:
                    for handler in finding.handlers:
                        if finding.side in ('5p', 'both'):
                            handler.start_is_biological_start = False
                        if finding.side in ('3p', 'both'):
                            handler.end_is_biological_end = False
                    # recorded whether or not it masks anything, as a wrong starting phase does not
                    transcript['detected_error_types'].add(finding.error_type)
                    ranges.append(self._finding_range(finding, buffered, sign))

                # the case of this super locus sharing genomic range with another one. Only the
                # exported transcript is masked for it, the ranges having been worked out from that
                # transcript's own span (see _compute_overlap_masks)
                if transcript['transcript'].longest and self._overlap_masks[i]:
                    transcript['detected_error_types'].add(types.SL_OVERLAP_ERROR)
                    ranges += self._overlap_masks[i]

                # a range can come out empty, e.g. a flank measured from the selected transcript
                # for an isoform reaching past it
                transcript['masks'] = [FeatureImporter(self.coord, self.is_plus_strand, types.GEENUFF_MASK,
                                                       start=sign * lo, end=sign * hi,
                                                       controller=self.controller)
                                       for lo, hi in self._merge_ranges(r for r in ranges if r[0] < r[1])]

                # the backwards introns standing in for overlapping exons are not saved, the
                # OVERLAPPING_EXONS error being descriptive enough
                transcript['introns'][:] = [x for x in transcript['introns'] if not self._is_backwards(x)]

    @staticmethod
    def _finding_range(finding, buffered, sign):
        """The sign-normalised range a finding masks, given the buffered span of its locus. A '5p'
        or '3p' extent runs from the finding's first handler out into the flank on that side,
        'whole' covers the buffered span, and a (start, end) extent is its own range."""
        if isinstance(finding.extent, tuple):
            start, end = finding.extent
            return sign * start, sign * end
        lo, hi = buffered
        if finding.extent == '5p':
            return lo, sign * finding.handlers[0].start
        if finding.extent == '3p':
            return sign * finding.handlers[0].end, hi
        return lo, hi

    def _sequence_bounds(self):
        """The sequence as a sign-normalised, half open range. On the minus strand the bases
        length - 1 down to 0 become -(length - 1) up to 0, so the exclusive end is 1."""
        if self.is_plus_strand:
            return 0, self.coord.length
        return -(self.coord.length - 1), 1


##### main flow control #####
class ImportController(object):
    def __init__(self, database_path, config={}, replace_db=False):
        self.database_path = database_path
        self.latest_genome = None
        self.stats = ImportStatistics()
        self._mk_session(replace_db)
        # queues for adding to db
        self.insertion_queues = InsertionQueue(session=self.session, engine=self.engine)
        self.config = {'min_intron_length': 20}  # the default config
        self.config.update(config)

    def _mk_session(self, replace_db):
        if os.path.exists(self.database_path):
            if replace_db:
                os.remove(self.database_path)
                logger.info('removed existing database at {}'.format(self.database_path))
            else:
                logger.error('database already existing at {} and --replace-db not set'.format(
                    self.database_path))
                exit()
        self.engine = create_engine(helpers.full_db_path(self.database_path), echo=False)
        orm.Base.metadata.create_all(self.engine)
        self.session = sessionmaker(bind=self.engine)()

    def make_genome(self, genome_args=None):
        if genome_args is None:
            genome_args = {}
        genome = orm.Genome(**genome_args)
        self.latest_fasta_importer = FastaImporter(genome)
        self.session.add(genome)
        self.session.commit()

    def run_analyze(self):
        # run ANALYZE; on the db for hopefully more performant queries
        logger.info('Running ANALYZE on the database')
        with self.engine.connect() as con:
            con.execute('ANALYZE;')

    def add_genome(self, fasta_path, gff_path, genome_args=None, clean_gff=True):
        if genome_args is None:
            genome_args = {}
        species = genome_args.get('species', 'unnamed genome')
        logger.info(f'Starting to add genome: {species}')
        logger.info(f'FASTA path: {fasta_path}')
        logger.info(f'GFF path: {gff_path}')

        self.clean_tmp_data()
        self.stats = ImportStatistics()
        self.add_sequences(fasta_path, genome_args)
        try:
            self.add_gff(gff_path, clean=clean_gff)
            self.run_analyze()
        except Exception:
            self.session.close()
            # the path may be given as an SQLAlchemy URL, and an in-memory database has no file
            file_path = self.database_path.removeprefix('sqlite:///')
            if os.path.isfile(file_path):
                part_path = f'{file_path}.partial'
                shutil.move(file_path, part_path)
                logger.error(f'Aborting due to error, attempt so far saved at {part_path} '
                             f'for debugging purposes')
            else:
                logger.error('Aborting due to error')
            raise
        self.stats.log_summary(species)

    def add_sequences(self, seq_path, genome_args=None):
        if genome_args is None:
            genome_args = {}
        if self.latest_genome is None:
            self.make_genome(genome_args)

        logger.info('Starting to add sequences from the FASTA file')
        self.latest_fasta_importer.add_sequences(seq_path)
        self.session.commit()
        logger.info(f'Added {len(self.latest_fasta_importer.genome.coordinates)} sequences')

    def clean_tmp_data(self):
        self.latest_genome = None
        self.latest_super_loci = []

    def add_gff(self, gff_file, clean=True):
        def insert_importer_groups(self, groups):
            """Initiates the calling of the add_to_queue() function of the importers
            in the correct order. Also initiates the insert of the many2many rows.
            """
            for group in groups:
                group['super_locus'].add_to_queue()
                # insert all features as well as transcript and protein related entries
                for transcript in group['transcripts']:
                    # make shortcuts
                    tp = transcript['transcript_piece']
                    tf = transcript['transcript_feature']
                    # add transcript handler that are always present
                    transcript['transcript'].add_to_queue()
                    tp.add_to_queue()
                    tf.add_to_queue()
                    tf.insert_feature_piece_association(tp.id)
                    # if coding transcript
                    if 'protein' in transcript:
                        transcript['protein'].add_to_queue()
                        transcript['protein'].insert_transcript_protein_association(transcript['transcript'].id)
                        transcript['cds'].insert_feature_protein_association(transcript['protein'].id)
                        transcript['cds'].add_to_queue()
                        transcript['cds'].insert_feature_piece_association(tp.id)
                    # if there are introns
                    if 'introns' in transcript:
                        for intron in transcript['introns']:
                            intron.add_to_queue()
                            intron.insert_feature_piece_association(tp.id)
                    # insert the masks, and which errors they are for
                    for mask in transcript['masks']:
                        mask.add_to_queue()
                        mask.insert_feature_piece_association(tp.id)
                    for error_type in transcript['detected_error_types']:
                        self.insertion_queues.transcript_error.queue.append(
                            {'transcript_id': transcript['transcript'].id, 'type': error_type})

        def clean_and_insert(self, groups, clean, is_final_coord):
            # a super locus whose own strand is neither '+' nor '-' is kept apart rather than
            # falling into the minus bucket, where it would otherwise be free to sort first and
            # give GFFErrorHandling a strand of None to work every neighbour comparison against
            plus = [g for g in groups if g['super_locus'].is_plus_strand is True]
            minus = [g for g in groups if g['super_locus'].is_plus_strand is False]
            unstranded = [g for g in groups if g['super_locus'].is_plus_strand is None]
            if clean:
                # check and correct for errors
                # do so for each strand seperately
                # all changes should be made by reference
                GFFErrorHandling(plus, self).resolve_errors()
                # reverse order on minus strand
                GFFErrorHandling(minus[::-1], self).resolve_errors()
            # tally the final, selected transcripts (one per coding locus) now that error
            # resolution is done, so 'error-free' reflects every check that ran
            for group in plus + minus + unstranded:
                # an unexported locus was never checked for errors either, so counting its
                # transcript here would report a gene nothing was ever looked at as error-free
                reason = group['super_locus'].excluded_from_export
                if reason is not None:
                    self.stats.unexported_super_loci[reason] += 1
                    continue
                for transcript in group['transcripts']:
                    # counted from the error types detected, not from what ended up masked: an
                    # error whose range comes out empty (e.g. a missing UTR of a gene nested in
                    # another, which has no unclaimed sequence to extend the mask into) is still a
                    # real finding about the transcript, and dropping it from the count would
                    # understate the error rate
                    if transcript['transcript'].longest:
                        self.stats.longest_transcripts += 1
                        if not transcript['detected_error_types']:
                            self.stats.longest_error_free_transcripts += 1
                    # each error type is counted once per transcript it occurs in, no matter
                    # how many times it occurs within that one transcript
                    for error_type in transcript['detected_error_types']:
                        self.stats.errors[error_type] += 1
            # insert importers
            insert_importer_groups(self, plus)
            insert_importer_groups(self, minus)
            insert_importer_groups(self, unstranded)
            if is_final_coord or self.insertion_queues.total_size() > 10000:  # todo: change to 1000, RAM issue fix maybe
                self.insertion_queues.execute_so_far()

        assert self.latest_fasta_importer is not None, 'No recent genome found'
        logger.info('Starting to parse the GFF file')
        self.latest_fasta_importer.mk_mapper(gff_file)
        gff_organizer = OrganizedGFFEntries(gff_file, self.stats)
        gff_organizer.load_organized_entries()

        organized_gff_entries = gff_organizer.organized_entries
        n_organized_gff_entries = len(organized_gff_entries)
        geenuff_importer_groups = []
        for i, seqid in enumerate(organized_gff_entries.keys()):
            for entry_group in organized_gff_entries[seqid]:
                organized_entries = OrganizedGFFEntryGroup(entry_group, self.latest_fasta_importer,
                                                           self)
                geenuff_importer_groups.append(organized_entries.get_geenuff_importers())
            # never do error checking across fasta sequence borders
            is_final_coord = (i == (n_organized_gff_entries - 1))
            clean_and_insert(self, geenuff_importer_groups, clean, is_final_coord)
            logger.info(f'Finished importing features from {len(geenuff_importer_groups)} super loci '
                        f'from coordinate with seqid {seqid} ({i + 1}/{n_organized_gff_entries})')
            geenuff_importer_groups = []


class Insertable(ABC):
    @abstractmethod
    def add_to_queue(self):
        pass


class FastaImporter(object):
    def __init__(self, genome):
        self.genome = genome
        self.mapper = None
        self._coords_by_seqid = None
        self._gffid_to_coords = None
        self._gff_seq_ids = None

    @property
    def gffid_to_coords(self):
        if not self._gffid_to_coords:
            self._gffid_to_coords = {}
            for gffid in self._gff_seq_ids:
                fa_id = self.mapper(gffid)
                x = self.coords_by_seqid[fa_id]
                self._gffid_to_coords[gffid] = x
        return self._gffid_to_coords

    @property
    def coords_by_seqid(self):
        if not self._coords_by_seqid:
            self._coords_by_seqid = {c.seqid: c for c in self.genome.coordinates}
        return self._coords_by_seqid

    def mk_mapper(self, gff_file=None):
        fa_ids = [e.seqid for e in self.genome.coordinates]
        if gff_file is not None:  # allow setup without ado when we know IDs match exactly
            self._gff_seq_ids = helpers.get_seqids_from_gff(gff_file)
        else:
            self._gff_seq_ids = fa_ids
        mapper, is_forward = helpers.two_way_key_match(fa_ids, self._gff_seq_ids)
        self.mapper = mapper

        if not is_forward:
            raise NotImplementedError("Still need to implement backward match if fasta IDs "
                                      "are subset of gff IDs")

    def add_sequences(self, seq_file):
        # todo, parallelize sequence & annotation format, then import directly from ~Slice
        for seqid, seq in self.parse_fasta(seq_file):
            coord = orm.Coordinate(sequence=seq,
                                   length=len(seq),
                                   seqid=seqid,
                                   sha1=helpers.sequence_hash(seq),
                                   genome=self.genome)
            logger.info(f'Added coordinate object for FASTA sequence with seqid {seqid} to the queue')

    def parse_fasta(self, seq_file, id_delim=' '):
        fp = fastahelper.FastaParser()
        for fasta_header, seq in fp.read_fasta(seq_file):
            seq = seq.upper()  # this may perform poorly
            seqid = fasta_header.split(id_delim)[0]
            yield seqid, seq


class SuperLocusImporter(Insertable):
    def __init__(self,
                 entry_type,
                 given_name,
                 controller,
                 coord=None,
                 is_plus_strand=None,
                 start=-1,
                 end=-1,
                 excluded_from_export=None):
        self.id = InsertCounterHolder.super_locus()
        self.entry_type = entry_type
        self.given_name = given_name
        self.controller = controller
        # not neccessary for insert but helpful for certain error cases
        self.coord = coord
        self.is_plus_strand = is_plus_strand
        self.start = start
        self.end = end
        # a reason string from types, or None; see orm.SuperLocus.excluded_from_export
        self.excluded_from_export = excluded_from_export

    def add_to_queue(self):
        to_add = {'type': self.entry_type, 'given_name': self.given_name, 'id': self.id,
                  'excluded_from_export': self.excluded_from_export}
        self.controller.insertion_queues.super_locus.queue.append(to_add)

    def __repr__(self):
        params = {'id': self.id, 'type': self.entry_type, 'given_name': self.given_name}
        return helpers.get_repr('SuperLocusImporter', params)


class FeatureImporter(Insertable):
    def __init__(self,
                 coord,
                 is_plus_strand,
                 feature_type,
                 controller,
                 start=-1,
                 end=-1,
                 given_name=None,
                 phase_5p=0,
                 phase=None,
                 is_truncated=False,
                 has_inframe_stop=False,
                 has_start_codon=True,
                 has_stop_codon=True,
                 has_overlapping_pieces=False,
                 score=None,
                 source=None):
        self.id = InsertCounterHolder.feature()
        self.coord = coord
        self.given_name = given_name
        self.is_plus_strand = is_plus_strand
        self.feature_type = feature_type
        # start/end may have to be adapted to geenuff
        self.start = start
        self.end = end
        self.phase_5p = phase_5p  # the file's own value, only used for the WRONG_PHASE_5P check
        # the phase actually saved to the db; defaults to phase_5p for every non-CDS feature
        # type (where phase is irrelevant and stays 0), but a CDS passes its own recomputed
        # value here instead of trusting the file's phase_5p (see _parse_gff_entries)
        self.phase = phase if phase is not None else phase_5p
        self.is_truncated = is_truncated  # only used for the TRUNCATED_CDS check
        self.has_inframe_stop = has_inframe_stop  # only used for the INFRAME_STOP_CODON check
        self.has_start_codon = has_start_codon  # only used for the MISSING_START_CODON check
        self.has_stop_codon = has_stop_codon  # only used for the MISSING_STOP_CODON check
        self.has_overlapping_pieces = has_overlapping_pieces  # only used for the OVERLAPPING_CDS check
        self.score = score
        self.source = source
        self.start_is_biological_start = True
        self.end_is_biological_end = True
        self.controller = controller

    def add_to_queue(self):
        feature = {
            'id': self.id,
            'type': self.feature_type,
            'given_name': self.given_name,
            'coordinate_id': self.coord.id,
            'is_plus_strand': self.is_plus_strand,
            'score': self.score,
            'source': self.source,
            'phase': self.phase,
            'start': self.start,
            'end': self.end,
            'start_is_biological_start': self.start_is_biological_start,
            'end_is_biological_end': self.end_is_biological_end,
        }
        self.controller.insertion_queues.feature.queue.append(feature)

    def insert_feature_piece_association(self, transcript_piece_id):
        features2pieces = {
            'feature_id': self.id,
            'transcript_piece_id': transcript_piece_id,
        }
        self.controller.insertion_queues.association_transcript_piece_to_feature.\
            queue.append(features2pieces)

    def insert_feature_protein_association(self, protein_id):
        features2protein = {
            'feature_id': self.id,
            'protein_id': protein_id,
        }
        self.controller.insertion_queues.association_protein_to_feature.\
            queue.append(features2protein)

    def set_start_end_from_gff(self, gff_start, gff_end):
        self.start, self.end = get_geenuff_start_end(gff_start, gff_end, self.is_plus_strand)

    def pos_cmp_key(self):
        sortable_start = self.start
        sortable_end = self.end
        if not self.is_plus_strand:
            sortable_start *= -1
            sortable_end *= -1
        return self.coord.seqid, self.is_plus_strand, sortable_start, sortable_end, self.feature_type

    def __repr__(self):
        params = {
            'id': self.id,
            'coord_id': self.coord.id,
            'type': self.feature_type,
            'is_plus_strand': self.is_plus_strand,
            'phase': self.phase_5p,
        }
        if self.given_name:
            params['given_name'] = self.given_name
        return helpers.get_repr('FeatureImporter', params, str(self.start) + '--' + str(self.end))


class TranscriptImporter(Insertable):
    def __init__(self, entry_type, given_name, super_locus_id, controller, longest=False):
        self.id = InsertCounterHolder.transcript()
        self.entry_type = entry_type
        self.given_name = given_name
        self.super_locus_id = super_locus_id
        self.controller = controller
        self.longest = longest

    def add_to_queue(self):
        transcript = self._get_params_dict()
        self.controller.insertion_queues.transcript.queue.append(transcript)

    def _get_params_dict(self):
        d = {
            'id': self.id,
            'type': self.entry_type,
            'given_name': self.given_name,
            'super_locus_id': self.super_locus_id,
            'longest': self.longest,
        }
        return d

    def __repr__(self):
        return helpers.get_repr('TranscriptImporter', self._get_params_dict())


class TranscriptPieceImporter(Insertable):
    def __init__(self, given_name, transcript_id, position, controller):
        self.id = InsertCounterHolder.transcript_piece()
        self.given_name = given_name
        self.transcript_id = transcript_id
        self.position = position
        self.controller = controller

    def add_to_queue(self):
        transcript_piece = self._get_params_dict()
        self.controller.insertion_queues.transcript_piece.queue.append(transcript_piece)

    def _get_params_dict(self):
        d = {
            'id': self.id,
            'given_name': self.given_name,
            'transcript_id': self.transcript_id,
            'position': self.position,
        }
        return d

    def __repr__(self):
        return helpers.get_repr('TranscriptPieceImporter', self._get_params_dict())


class ProteinImporter(Insertable):
    def __init__(self, given_name, super_locus_id, controller):
        self.id = InsertCounterHolder.protein()
        self.given_name = given_name
        self.super_locus_id = super_locus_id
        self.controller = controller

    def add_to_queue(self):
        protein = self._get_params_dict()
        self.controller.insertion_queues.protein.queue.append(protein)

    def _get_params_dict(self):
        d = {
            'id': self.id,
            'given_name': self.given_name,
            'super_locus_id': self.super_locus_id,
        }
        return d

    def insert_transcript_protein_association(self, transcript_id):
        transcript2protein = {
            'transcript_id': transcript_id,
            'protein_id': self.id,
        }

        self.controller.insertion_queues.association_transcript_to_protein.\
            queue.append(transcript2protein)

    def __repr__(self):
        return helpers.get_repr('ProteinImporter', self._get_params_dict())
