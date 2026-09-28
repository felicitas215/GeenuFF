import os
import math
import logging
import shutil
from collections import defaultdict

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
                                  STOP_CODONS, in_enum_values)

logger = logging.getLogger(__name__)


class GFFValidityError(Exception):
    pass


# how much of a locus its own errors take away, lower being better (see
# GFFErrorHandling._locus_severity); used to rank candidates when an overlap is resolved
NOT_MASKED = 0       # nothing wrong with it
FLANK_MASKED = 1     # an edge is unknown, e.g. a missing UTR, but the coding sequence is correct
MASKED_OUTRIGHT = 2  # the coding sequence itself is masked, so keeping the locus recovers nothing


# what each orm.SuperLocus.excluded_from_export value means, for the import summary
UNEXPORTED_REASONS = {
    types.UNPLACEABLE_STRAND: 'their transcript or its exon/CDS lines are not all on one '
                              'definite strand',
    types.UNPLACEABLE_COORDINATES: 'a line of theirs runs backwards, its start past its end',
    types.OVERLAP_DROPPED: 'dropped so an overlapping partner could be kept whole (see below)',
}


class ImportStatistics(object):
    """Aggregate counts for one add_genome() call. A single summary logged at the end of
    the import gives a full picture of what was found (see log_summary), instead of one
    log line per locus/transcript/error, which floods the log without being any easier to
    get an overview from."""

    def __init__(self) -> None:
        self.total_super_loci: int = 0
        self.total_transcripts: int = 0
        self.total_coding_transcripts: int = 0
        # transcript.longest == True, one per exported gene; a gene left out of exports is not
        # counted, nothing about it having been checked either (see clean_and_insert)
        self.longest_transcripts: int = 0
        self.longest_error_free_transcripts: int = 0  # of the above, the ones with no error at all
        self.empty_super_loci: int = 0  # a gene with no transcripts at all
        self.unstranded_super_loci: int = 0  # e.g. NCBI's '?' strand for trans-spliced genes
        # kept in full but never exported, counted per reason and keyed by the value stored in
        # orm.SuperLocus.excluded_from_export
        self.unexported_super_loci: defaultdict[str, int] = defaultdict(int)
        # transcripts left out entirely, no line of theirs having a storable range
        self.unstorable_transcripts_dropped: int = 0
        self.discontinuous_gene_ids_reused: int = 0
        self.gene_parented_duplicates_dropped: int = 0  # see OrganizedGFFEntryGroup._is_parented_to_super_locus
        self.gene_parented_features_reparented: int = 0  # same, but attached to the sole transcript instead
        self.gene_parented_ambiguous_dropped: int = 0  # same, but position matched no transcript
        self.backwards_errors_removed: int = 0  # from overlapping super loci, see _remove_backwards_errors
        # every pair of super loci sharing genomic range, of whatever kind, as recorded in the
        # super_locus_overlap table; these mask nothing, see _record_overlap_pairs
        self.overlap_pairs_recorded: int = 0
        self.super_loci_in_overlap_pairs: int = 0
        # loci whose exported transcript overlaps another's, counted once each however many partners it has
        self.super_loci_overlapping_exported: int = 0
        # of the pairs where both genes are exported, how each was settled; the three always sum
        # to that number (see GFFErrorHandling._compute_overlap_masks)
        self.overlap_pairs_resolved: int = 0  # one locus kept whole, the other dropped
        self.overlap_pairs_refused: int = 0  # neither could be kept, both masked over the shared range
        self.overlap_pairs_in_chains: int = 0  # a partner overlaps a third locus, so not decided
        self.overlap_loci_dropped: int = 0
        # keyed by types.Errors value, counted per transcript from the error types detected
        # rather than from the error features inserted, so an error left with a zero-length
        # range to mask still shows up here (see docs/spec_vs_gff.md)
        self.errors: defaultdict[str, int] = defaultdict(int)
        self.unrecognized_feature_types: defaultdict[str, int] = defaultdict(int)  # keyed by the raw, unknown GFF type

    def log_summary(self, species: str) -> None:
        """Logs one summary of the import, as sections of counted lines with the counts in a
        column, so a reader can scan down them rather than through a paragraph."""
        overlap_pairs = (self.overlap_pairs_resolved + self.overlap_pairs_refused
                         + self.overlap_pairs_in_chains)
        sections = [
            ('loci', [
                (self.total_super_loci, 'gene lines in the file'),
                (self.empty_super_loci, 'of them with no transcript under them'),
                (self.discontinuous_gene_ids_reused, 'gene IDs used by more than one gene line'),
                (self.unstranded_super_loci, "gene lines on neither strand, e.g. '.' or NCBI's '?'"),
            ] + [(count, f'genes kept but never exported: {UNEXPORTED_REASONS[reason]}')
                 for reason, count in sorted(self.unexported_super_loci.items())]),
            ('transcripts', [
                (self.total_transcripts, 'transcripts in the file'),
                (self.total_coding_transcripts, 'of them with a CDS'),
                (self.longest_transcripts, 'selected for export, one per locus that will be exported'),
                (self.longest_error_free_transcripts, 'of those selected with no error at all'),
                (self.unstorable_transcripts_dropped, 'dropped, no line of theirs having a '
                                                      'storable range'),
            ]),
            ('transcripts overlapping another transcript that will be exported on the same strand', [
                (self.super_loci_overlapping_exported, 'transcripts overlap at least one other'),
                (overlap_pairs, 'pairs between them, each settled one of three ways:'),
                (self.overlap_pairs_resolved, 'pairs where one transcript was kept whole and its '
                                              'partner dropped from the export'),
                (self.overlap_pairs_refused, 'pairs where neither transcript could be kept, so both '
                                             'are masked over the sequence they share'),
                (self.overlap_pairs_in_chains, 'pairs left alone, one of the two overlapping a '
                                               'further transcript as well'),
            ]),
            ('exon/CDS lines naming their gene as the parent instead of a transcript', [
                (self.gene_parented_duplicates_dropped, 'dropped as duplicates of a transcript\'s '
                                                        'own lines'),
                (self.gene_parented_features_reparented, 'given to the gene\'s sole transcript'),
                (self.gene_parented_ambiguous_dropped, 'dropped, matching no transcript'),
            ]),
        ]
        if self.errors:
            sections.append(('errors, counted once per transcript they occur in, however often '
                             'they occur in that one',
                             [(count, error_type)
                              for error_type, count in sorted(self.errors.items())]))
        if self.unrecognized_feature_types:
            sections.append(('lines skipped for a feature type GeenuFF has no use for',
                             [(count, feature_type) for feature_type, count
                              in sorted(self.unrecognized_feature_types.items())]))

        width = max(len(str(count)) for _, entries in sections for count, _ in entries)
        lines = [f'import summary for "{species}":']
        for heading, entries in sections:
            lines.append(f'  {heading}:')
            lines += [f'    {count:>{width}}  {text}' for count, text in entries]
        for note in self._summary_notes():
            lines.append(f'  note: {note}')
        logger.info('\n'.join(lines))

    def _summary_notes(self) -> list[str]:
        """Explanations for the parts of the summary above that a count alone does not convey."""
        notes = []
        if types.SL_OVERLAP_ERROR in self.errors:
            notes.append(f"'{types.SL_OVERLAP_ERROR}' marks sequence claimed by two transcripts at "
                         f'once: on a transcript that was kept, the region its dropped partner used to '
                         f'cover; on a pair where neither could be kept, the sequence they share')
        notes.append(f'errors are counted as found, so one whose masked range came out empty, '
                     f'such as a missing UTR with no room to mask it in, still counts '
                     f'(see docs/spec_vs_gff.md)')
        if self.backwards_errors_removed:
            notes.append(f'{self.backwards_errors_removed} error masks were discarded for coming '
                         f'out with their start and end the wrong way round')
        notes.append(f'the super_locus_overlap table separately records {self.overlap_pairs_recorded} '
                     f'pairs over {self.super_loci_in_overlap_pairs} coding genes, counting every '
                     f'coding transcript a gene has rather than only the exported one; it masks '
                     f'nothing and is there for consumers other than Helixer')
        return notes


# core queue prep
class InsertionQueue(helpers.QueueController):
    def __init__(self, session, engine):
        super().__init__(session, engine)
        self.super_locus = helpers.CoreQueue(orm.SuperLocus.__table__.insert())
        self.transcript = helpers.CoreQueue(orm.Transcript.__table__.insert())
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
            self.super_locus, self.transcript, self.transcript_piece, self.protein, self.feature,
            self.association_transcript_piece_to_feature, self.association_protein_to_feature,
            self.association_transcript_to_protein, self.super_locus_overlap
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
                'errors': [],  # errors are filled in later
                'detected_error_types': set(),  # every error type found, see _add_error
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

    def _parse_gff_entries(self, entries):
        """Changes the GFF format into the GeenuFF format. Does all the parsing."""
        sl = entries['super_locus']
        stats = self.controller.stats
        stats.total_super_loci += 1
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
            t_importers = {'errors': [], 'detected_error_types': set()}
            # check for multi inheritance and throw NotImplementedError if found
            if t.get_Parent() is None:
                raise GFFValidityError(f"transcript level feature without Parent found {t}, attributes: {t.attributes}")
            if len(t.get_Parent()) > 1:
                raise NotImplementedError
            t_id = t.get_ID()
            if t.start > t.end:
                self._drop_unstorable_transcript(sl_i, t_id, 'transcript', t.start, t.end)
                continue
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
                protein_id = self._get_protein_id_from_cds_list(t_entries['cds'])
                p_i = ProteinImporter(given_name=protein_id,
                                      super_locus_id=sl_i.id,
                                      controller=self.controller)
                # create coding features from exon limits
                if t_is_plus_strand:
                    phase_5p = t_entries['cds'][0].phase
                    phase_3p = t_entries['cds'][-1].phase
                else:
                    phase_5p = t_entries['cds'][-1].phase
                    phase_3p = t_entries['cds'][0].phase
                # spliced once here, from the raw per-piece CDS list we already have on hand,
                # to check for a truncated (non-codon-multiple) length, a premature stop codon,
                # or a missing start/stop codon, none of which the file's own annotation is
                # trusted to rule out. This also makes the start/stop codon check safe when the
                # terminal CDS piece is shorter than 3bp (the codon is split across a splice
                # junction right at the transcript's edge): reading 3 contiguous genomic bases
                # from the boundary, as the old per-feature check did, would read into the
                # intron and report a false MISSING_START/STOP_CODON in that case.
                cds_seq = spliced_cds_sequence(self.coord.sequence, t_entries['cds'], t_is_plus_strand)
                cds_i = FeatureImporter(self.coord,
                                        t_is_plus_strand,
                                        types.GEENUFF_CDS,
                                        phase_5p=phase_5p,
                                        phase_3p=phase_3p,
                                        # the file's starting phase isn't trusted; the first
                                        # CDS piece always starts a fresh codon, i.e. phase 0
                                        # -> same value in both phase conventions, Helixer later
                                        # computes the phases from the starting phase
                                        phase=0,
                                        is_truncated=len(cds_seq) % 3 != 0,
                                        has_inframe_stop=has_inframe_stop_codon(cds_seq),
                                        has_start_codon=cds_seq[:3] == START_CODON,
                                        has_stop_codon=cds_seq[-3:] in STOP_CODONS,
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
        """Looks for the transcript with the longest exon length (cds_length - sum(intron_lengths))
        and sets the 'longest' parameters in all the TranscriptImporters. In case of a tie, the
        first transcript found will be set as longest."""
        def filter_coding_introns(cds_i, introns):
            """Filters the non-coding introns out of introns"""
            cds_range = sorted([cds_i.start, cds_i.end])
            coding_introns = []
            for intron in introns:
                intron_range = sorted([intron.start, intron.end])
                if min(cds_range[1], intron_range[1]) - max(cds_range[0], intron_range[0]) > 0:
                    coding_introns.append(intron)
            return coding_introns

        max_exon_len = -1
        longest_importer = None
        for t in self.importers['transcripts']:
            if 'cds' in t:
                cds_len = abs(t['cds'].start - t['cds'].end)
                coding_introns = filter_coding_introns(t['cds'], t['introns'])
                intron_lengths = sum([abs(i.start - i.end) for i in coding_introns])
                exon_len = cds_len - intron_lengths
                if exon_len > max_exon_len:
                    max_exon_len = exon_len
                    longest_importer = t['transcript']
        for t in self.importers['transcripts']:
            if t['transcript'] is longest_importer:
                t['transcript'].longest = True
            else:
                t['transcript'].longest = False

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
    """Takes an entry group (all entries of one super locus) and stores the entries
    in an orderly fashion. Can then return a corresponding OrganizedGeenuffImporterGroup.
    Does not perform error checking, which happens later.

    The entries are organized in the following way:

    entries = {
        'super_locus' = super_locus_entry,
        'transcripts' = {
            transcript_entry1: {
                'exons': [ordered_exon_entry1, ordered_exon_entry2, ...],
                'cds': [ordered_cds_entry1, ordered_cds_entry2, ...]
            },
            transcript_entry2: {
                'exons': [ordered_exon_entry1, ordered_exon_entry2, ...],
                'cds': [ordered_cds_entry1, ordered_cds_entry2, ...]
            },
            ...
        }
    }
    """

    def __init__(self, gff_entry_group, fasta_importer, controller):
        self.fasta_importer = fasta_importer
        self.controller = controller
        self.entries = {'transcripts': {}}
        self.coord = None
        self.add_gff_entry_group(gff_entry_group)

    def add_gff_entry_group(self, entries):
        latest_transcript = None
        for entry in list(entries):
            if in_enum_values(entry.type, types.SuperLocusAll):
                assert 'super_locus' not in self.entries
                self.entries['super_locus'] = entry
            elif in_enum_values(entry.type, types.TranscriptLevel):
                self.entries['transcripts'][entry] = {'exons': [], 'cds': []}
                latest_transcript = entry
            elif latest_transcript is not None:
                is_exon = in_enum_values(entry.type, types.ExonLevel)
                is_cds = in_enum_values(entry.type, types.CDSLevel)
                key = 'exons' if is_exon else 'cds'
                if (is_exon or is_cds) and self._is_parented_to_super_locus(entry):
                    # some GFF3 sources (e.g. NCBI/EMBL) redundantly echo a transcript's
                    # already-nested exon/CDS a second time as a standalone feature parented
                    # directly to the gene; grouping is otherwise positional (not Parent-id
                    # aware), so without this check such a line would get glued onto whatever
                    # transcript happens to be "latest", corrupting it with a bogus overlap
                    if self._matches_existing_transcript_feature(entry, key):
                        # a duplicate echo of something already attached: drop it, keeping
                        # the copy that's already correctly attached
                        self.controller.stats.gene_parented_duplicates_dropped += 1
                        logger.debug(f"skipping {entry.type} parented directly to the gene "
                                     f"'{self.entries['super_locus'].get_ID()}' instead of a "
                                     f"transcript (a redundant echo of an already-nested "
                                     f"feature)")
                    elif len(self.entries['transcripts']) == 1:
                        # not a duplicate, but this gene has exactly one transcript declared
                        # so far, so Parent=gene can only mean this one transcript; unlike
                        # the ambiguous case below, there is no other candidate to confuse it
                        # with, so it's safe to attach directly (mirrors how positional
                        # grouping already handles this correctly for a single-transcript
                        # gene when the line isn't gene-parented at all)
                        self.controller.stats.gene_parented_features_reparented += 1
                        logger.debug(f"reparenting {entry.type} from gene "
                                     f"'{self.entries['super_locus'].get_ID()}' to its sole "
                                     f"transcript '{latest_transcript.get_ID()}'")
                        self.entries['transcripts'][latest_transcript][key].append(entry)
                    else:
                        # neither a known duplicate nor safely attributable to a sole
                        # transcript: genuinely ambiguous (which of several transcripts, if
                        # any, does it belong to?), so it's dropped, but flagged louder
                        self.controller.stats.gene_parented_ambiguous_dropped += 1
                        logger.warning(f"skipping {entry.type} parented directly to the gene "
                                       f"'{self.entries['super_locus'].get_ID()}' instead of a "
                                       f"transcript, and its position matches no transcript "
                                       f"already collected for this gene; check this gene by "
                                       f"hand")
                elif is_exon:
                    self.entries['transcripts'][latest_transcript]['exons'].append(entry)
                elif is_cds:
                    self.entries['transcripts'][latest_transcript]['cds'].append(entry)
                else:
                    logger.warning(f'Found unexpected entry type: {entry.type}')
            else:
                # is overprinting
                logger.debug(f'Ignoring {entry.type} without transcript found in {entry.seqid}: '
                             f'{entries[0].attribute}')
                # todo: add counter?

        # set the coordinate
        self.coord = self.fasta_importer.gffid_to_coords[self.entries['super_locus'].seqid]

        # order exon and cds lists by start value (disregard strand for now)
        for _, value_dict in self.entries['transcripts'].items():
            for key in ['exons', 'cds']:
                value_dict[key].sort(key=lambda e: e.start)

    def _is_parented_to_super_locus(self, entry):
        """True if entry's Parent= names the gene itself rather than any specific
        transcript"""
        parent_ids = entry.get_Parent()
        return bool(parent_ids) and self.entries['super_locus'].get_ID() in parent_ids

    def _matches_existing_transcript_feature(self, entry, key):
        """True if (entry.start, entry.end) already belongs to some transcript of this
        gene collected so far, under 'exons' or 'cds' (whichever key is given). Order
        dependent: a gene-parented duplicate seen before the transcript it echoes has been
        fully collected will not match yet."""
        position = (entry.start, entry.end)
        return any(position == (existing.start, existing.end)
                   for t_entries in self.entries['transcripts'].values()
                   for existing in t_entries[key])

    def get_geenuff_importers(self):
        geenuff_importer_group = OrganizedGeenuffImporterGroup(self.entries, self.coord,
                                                               self.controller)
        return geenuff_importer_group.importers


class OrganizedGFFEntries(object):
    """Structures the gff entries coming from gffhelper by seqid and gene. Also does some
    basic gff value cleanup.
    The entries are organized in the following way:

    organized_entries = {
        'seqid1': [
            [gff_entry1_gene1, gff_entry2_gene1, ...],
            [gff_entry1_gene2, gff_entry2_gene2, ...],
        ],
        'seqid2': [
            [gff_entry1_gene1, gff_entry2_gene1, ...],
            [gff_entry1_gene2, gff_entry2_gene2, ...],
        ],
        ...
    }
    """

    def __init__(self, gff_file, stats=None):
        self.gff_file = gff_file
        self.organized_entries = {}
        # only present when constructed as part of a full ImportController.add_gff() run;
        # standalone/test use (no controller) gets its own throwaway counter
        self.stats = stats if stats is not None else ImportStatistics()

    def load_organized_entries(self):
        self.organized_entries = {}
        gene_level = [x.value for x in types.SuperLocusAll]

        reader = self._useful_gff_entries()
        first = next(reader, None)

        if first is not None:
            seqid = first.seqid
            gene_group = [first]
            # GFF3 allows one ID to be split across several 'gene' lines (a discontinuous
            # feature); these are kept as separate super locus records below rather than
            # merged into one, since merging would make any real gene lying between the two
            # occurrences falsely look nested/overlapping in _compute_overlap_masks
            seen_gene_ids = {first.get_ID()} - {None}
            self.organized_entries[seqid] = []

            for entry in reader:
                if entry.type in gene_level:
                    entry_id = entry.get_ID()
                    if entry_id is not None:
                        if entry_id in seen_gene_ids:
                            self.stats.discontinuous_gene_ids_reused += 1
                            logger.debug(f"'gene' ID '{entry_id}' reused at "
                                         f"{entry.seqid}:{entry.start}-{entry.end}; kept "
                                         f"as a separate super locus record")
                        seen_gene_ids.add(entry_id)
                    self.organized_entries[seqid].append(gene_group)
                    gene_group = [entry]
                    if entry.seqid != seqid:
                        self.organized_entries[entry.seqid] = []
                        seqid = entry.seqid
                else:
                    gene_group.append(entry)
            self.organized_entries[seqid].append(gene_group)

    def _useful_gff_entries(self):
        skipable = [x.value for x in types.IgnorableGFFFeatures]
        reader = self._gff_gen()
        for entry in reader:
            if entry.type not in skipable:
                yield entry

    def _gff_gen(self):
        known = [x.value for x in types.AllKnownGFFFeatures]
        reader = gffhelper.read_gff_file(self.gff_file)
        for entry in reader:
            if entry.type not in known:
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
    super loci, and looks for errors. Error features may be inserted and importers be
    removed when deemed necessary.
    """

    def __init__(self, geenuff_importer_groups, controller):
        self.groups = geenuff_importer_groups
        self.controller = controller
        self._overlap_masks = []
        if self.groups:
            self.is_plus_strand = self.groups[0]['super_locus'].is_plus_strand
            self.coord = self.groups[0]['super_locus'].coord
            # make sure self.groups is sorted correctly
            self.groups.sort(key=lambda g: g['super_locus'].start, reverse=not self.is_plus_strand)
            self._compute_overlap_masks()

    def _compute_overlap_masks(self):
        """Settles what happens to every pair of super loci on this strand that share genomic
        range, in two records serving two different consumers.

        The unfiltered record goes to the super_locus_overlap table: one row per pair, for every
        kind of locus, including ones with no transcript at all. It is the annotation's geometry
        as given, it masks nothing, and no consumer that reads error features sees it.

        The second record decides what an export is handed. Two loci cannot share a base: one
        label per base is written, so whichever is written last silently wins. Where one of the
        two can be kept without mislabelling the other's coding sequence, it is kept whole and
        its partner is dropped from the export, leaving one coherent gene rather than two with a
        hole through them (see _decide_overlap_pair). Where neither can be kept, both genes are
        masked over their whole length.
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
        resolved, masked_whole = [], set()
        for i, j in pairs:
            keeper = dropped = mask = None
            # only an isolated pair is decided; in a chain a dropped locus can sit between two
            # kept ones, and the one masked region recorded per locus cannot express that
            isolated = len(partners[i]) == 1 and len(partners[j]) == 1
            if isolated:
                keeper, dropped, mask = self._decide_overlap_pair(i, j, extents, sign)
            if keeper is None:
                # nothing could be saved here, so both genes are masked end to end; a locus in
                # two undecided pairs lands in the set twice over and is masked once
                masked_whole.update((i, j))
                if isolated:
                    stats.overlap_pairs_refused += 1
                else:
                    stats.overlap_pairs_in_chains += 1
                continue
            # marking it here, inside the loop, is what makes the dropped locus invisible to
            # _exported_transcript, and so to the neighbour search the second loop relies on
            self.groups[dropped]['super_locus'].excluded_from_export = types.OVERLAP_DROPPED
            resolved.append((keeper, dropped, mask))
            stats.overlap_pairs_resolved += 1
            stats.overlap_loci_dropped += 1

        # masks are extended only now that every drop is settled, so a border is measured toward
        # a gene that really is exported
        for i in masked_whole:
            masks[i].append(self._extend_unknown_ends(i, extents[i], extents))
        for keeper, dropped, mask in resolved:
            nested = (extents[keeper][0] <= extents[dropped][0]
                      and extents[dropped][1] <= extents[keeper][1])
            if not nested:
                # a nested gene's hole is interior, both its sides being the kept gene's own
                # intron or UTR, so there is nothing out there to extend into. Crossing, the
                # dropped gene sticks out on whichever side it starts or ends past the kept one,
                # and that is the only side its own end can be unknown on
                outward = 1 if extents[dropped][0] > extents[keeper][0] else -1
                mask = self._extend_unknown_ends(dropped, mask, extents, sides=(outward,))
            # the mask goes on the locus that is kept, the dropped one's own features never
            # reaching an export to carry it (see orm.SuperLocus.excluded_from_export)
            if mask[0] < mask[1]:  # two loci with identical spans leave nothing to mask
                masks[keeper].append(mask)
        # back out of sign-normalised coordinates into the strand-oriented ones the features use,
        # merging first so a locus masked from several pairs ends up with the fewest ranges
        self._overlap_masks = [[(sign * lo, sign * hi) for lo, hi in self._merge_ranges(ranges)]
                               for ranges in masks]

    def _extend_unknown_ends(self, i, mask, extents, sides=(-1, 1)):
        """Runs a mask on past the gene it covers, on whichever of the given sides that gene's
        end is not where the gene really ended: a truncated CDS, a missing start/stop codon or an
        unannotated UTR all leave it unknown how much further the gene ran, so the sequence out
        there cannot be taught as intergenic either. Reaches as far as any other error mask
        would, the border with the next exported gene."""
        lo, hi = mask
        for direction in sides:
            if self._end_is_known(self.groups[i], direction):
                continue
            limit = self._outward_limit(i, direction, extents)
            if direction > 0:
                hi += self._border_offset(limit - hi)
            else:
                lo -= self._border_offset(lo - limit)
        return lo, hi

    def _outward_limit(self, i, direction, extents):
        """How far out a mask on this gene may reach: the next exported gene that way, or the end
        of the sequence. A gene overlapping it gives no room at all, the gap being negative.

        Neighbours are taken in sort order, which for a gene nested inside another is not the
        nearest one by coordinate, so a mask can reach into a gene further along. That only
        happens between genes that overlap, and those are masked whole anyway, so the masked
        sequence comes out the same."""
        neighbour = self._exported_neighbour(i, direction)
        if neighbour is not None:
            return extents[neighbour][0] if direction > 0 else extents[neighbour][1]
        lower, upper = (0, self.coord.length) if self.is_plus_strand else (-self.coord.length, 0)
        return upper if direction > 0 else lower

    def _end_is_known(self, group, direction):
        """Whether a gene's annotated end on one side is where the gene really ends. Sign-normalised
        coordinates run 5p to 3p, so direction 1 asks about its 3p end and -1 about its 5p one."""
        transcript = self._selected_transcript(group)
        cds, tf = transcript['cds'], transcript['transcript_feature']
        if cds.is_truncated or cds.has_inframe_stop:
            return False  # the reading frame is wrong, so neither end can be trusted
        if direction > 0:
            return cds.end != tf.end and cds.has_stop_codon
        return cds.start != tf.start and cds.has_start_codon

    def _decide_overlap_pair(self, i, j, extents, sign):
        """Works out whether one locus of an overlapping pair can be kept, and what has to be
        masked in place of the other. Returns (keeper, dropped, mask range) or three Nones.

        Keeping one means its labels cover the range the two share, so the pair is only decided
        where that stays true of the sequence underneath: the shared range must not hold the
        dropped locus' coding sequence, which would then read as UTR or intron. With coding
        sequence from both there is nothing to choose and the pair is left alone.

        How much is masked depends on the shape. A nested locus is masked over its whole span,
        which is interior to the other and so leaves both of that one's ends visible. A crossing
        partner is masked only over the part sticking out past the locus that is kept; identical
        spans need no mask at all, the kept locus covering the same bases. Either way the mask
        may run further out still, see _extend_past_unknown_end.
        """
        spans = {k: extents[k] for k in (i, j)}
        overlap = (max(spans[i][0], spans[j][0]), min(spans[i][1], spans[j][1]))
        codes = {k: self._codes_within(self.groups[k], overlap, sign) for k in (i, j)}
        if codes[i] and codes[j]:
            return None, None, None

        outer, inner = None, None
        if spans[i][0] <= spans[j][0] and spans[j][1] <= spans[i][1]:
            outer, inner = i, j
        elif spans[j][0] <= spans[i][0] and spans[i][1] <= spans[j][1]:
            outer, inner = j, i
        if outer is not None:
            # only the outer locus can be kept here: keeping the inner instead would mean masking
            # the whole outer one to save the small gene inside it
            if self._locus_severity(self.groups[outer]) == MASKED_OUTRIGHT:
                return None, None, None
            return outer, inner, spans[inner]

        # crossing: whoever owns coding sequence in the shared range has to be the one kept
        forced = i if codes[i] else (j if codes[j] else None)
        options = [k for k in (i, j) if (forced is None or k == forced)
                   and self._locus_severity(self.groups[k]) != MASKED_OUTRIGHT]
        if not options:
            return None, None, None
        keeper = min(options, key=lambda k: (self._locus_severity(self.groups[k]),
                                             -self._coding_length(self.groups[k]), k))
        dropped = j if keeper == i else i
        # only the part sticking out past the kept locus is masked, the shared range keeping the
        # kept locus' labels, which are true for it
        mask = ((spans[dropped][0], spans[keeper][0]) if spans[dropped][0] < spans[keeper][0]
                else (spans[keeper][1], spans[dropped][1]))
        return keeper, dropped, mask

    def _codes_within(self, group, span, sign):
        """Whether any exonic coding base of this locus' selected transcript falls in the range.
        The cds feature spans the introns inside it, so those are taken back out."""
        transcript = self._selected_transcript(group)
        cds = transcript['cds']
        lo = max(sign * cds.start, span[0])
        hi = min(sign * cds.end, span[1])
        if lo >= hi:
            return False
        coding = hi - lo
        for intron in transcript['introns']:
            coding -= max(0, min(hi, sign * intron.end) - max(lo, sign * intron.start))
        return coding > 0

    def _coding_length(self, group):
        """Spliced coding length of the selected transcript, the same measure that picks the
        longest isoform, so one rule ranks both."""
        transcript = self._selected_transcript(group)
        cds = transcript['cds']
        length = abs(cds.end - cds.start)
        for intron in transcript['introns']:
            lo, hi = sorted((intron.start, intron.end))
            c_lo, c_hi = sorted((cds.start, cds.end))
            length -= max(0, min(c_hi, hi) - max(c_lo, lo))
        return length

    @staticmethod
    def _locus_severity(group):
        """How much of a locus its own errors already take away, lower being better. Keeping a
        locus whose coding sequence is masked anyway recovers nothing, while one that only lacks
        a UTR or a codon at its edge still has its coding sequence labelled: the mask for those
        covers the flank, not the body (see _add_overlapping_error)."""
        transcript = GFFErrorHandling._selected_transcript(group)
        cds, tf = transcript['cds'], transcript['transcript_feature']
        if cds.is_truncated or cds.has_inframe_stop:
            return MASKED_OUTRIGHT
        if (cds.start == tf.start or cds.end == tf.end
                or not cds.has_start_codon or not cds.has_stop_codon):
            return FLANK_MASKED
        return NOT_MASKED

    @staticmethod
    def _coding_extent(group, sign):
        """Everything a locus' coding transcripts reach over, as a sign-normalised ascending
        range. Taken from the transcripts rather than the gene line, which for one line of a
        discontinuous feature can be far narrower than the transcripts hanging off it."""
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

    def _3p_cds_start(self, transcript):
        """returns the start of the 3p most cds feature"""
        cds = transcript['cds']
        start = cds.start
        # introns are ordered by coordinate with no respect to strand
        intron_ends = [x.end for x in transcript["introns"]]
        if self.is_plus_strand:
            i_ends_within = [i for i in intron_ends if cds.start < i < cds.end]
            if i_ends_within:
                start = max(i_ends_within)
        else:
            i_ends_within = [i for i in intron_ends if cds.end < i < cds.start]
            if i_ends_within:
                start = min(i_ends_within)
        return start

    def resolve_errors(self):
        for i, group in enumerate(self.groups):
            # a locus that is not exported at all is not analysed either, and gets no mask: where
            # its features belong is exactly what is unknown about it, so any mask would be
            # guesswork, either covering the wrong strand or covering sequence that is fine
            if group['super_locus'].excluded_from_export is not None:
                continue

            # the case of no transcript for a super locus
            if not group['transcripts']:
                self.controller.stats.empty_super_loci += 1
                logger.debug('{} is a gene without any transcripts; it has no features and will '
                             'never be exported or masked'.format(group['super_locus'].given_name))
            # other cases
            for transcript in group['transcripts']:
                # if coding transcript
                if 'cds' in transcript:
                    cds = transcript['cds']
                    introns = transcript['introns']
                    tf = transcript['transcript_feature']

                    # the case of missing of implicit UTR ranges
                    # the solution is similar to the one above
                    if cds.start == tf.start:
                        self._add_overlapping_error(i, transcript, cds, '5p', types.MISSING_UTR_5P,
                                                    mark_other_handlers=[tf])
                    if cds.end == tf.end:
                        self._add_overlapping_error(i, transcript, cds, '3p', types.MISSING_UTR_3P,
                                                    mark_other_handlers=[tf])

                    # the case of missing start/stop codon; already computed once in
                    # _parse_gff_entries, from the spliced CDS sequence (see there for why)
                    if not cds.has_start_codon:
                        self._add_overlapping_error(i, transcript, cds, '5p', types.MISSING_START_CODON)
                    if not cds.has_stop_codon:
                        self._add_overlapping_error(i, transcript, cds, '3p', types.MISSING_STOP_CODON)

                    # the case of a truncated (not a multiple of 3) or a premature-stop-codon
                    # containing CDS; both were already computed once in _parse_gff_entries
                    if cds.is_truncated:
                        self._add_error(i, transcript, cds.start, cds.end, self.is_plus_strand,
                                        types.TRUNCATED_CDS)
                    if cds.has_inframe_stop:
                        self._add_error(i, transcript, cds.start, cds.end, self.is_plus_strand,
                                        types.INFRAME_STOP_CODON)

                    # the case of wrong 5p phase
                    if cds.phase_5p != 0:
                        self._add_overlapping_error(i, transcript, cds, '5p', types.WRONG_PHASE_5P)

                    # the case of this super locus sharing genomic range with another one. Only
                    # the exported transcript carries the mask, the range having been worked out
                    # from that transcript's own span, so this error type occurs once per gene
                    # (see _compute_overlap_masks)
                    if transcript['transcript'].longest:
                        for error_start, error_end in self._overlap_masks[i]:
                            self._add_error(i, transcript, error_start, error_end,
                                            self.is_plus_strand, types.SL_OVERLAP_ERROR)

                    if introns:
                        # the case of wrong 3p phase
                        len_3p_exon = abs(cds.end - self._3p_cds_start(transcript))
                        if cds.phase_3p != len_3p_exon % 3:
                            self._add_overlapping_error(i, transcript, cds, '3p',
                                                        types.MISMATCHED_PHASE_3P)

                    faulty_introns = []
                    for j, intron in enumerate(introns):
                        # the case of overlapping exons
                        if ((tf.is_plus_strand and intron.end < intron.start)
                                or (not self.is_plus_strand and intron.end > intron.start)):
                            # mark the overlapping cds regions as errors
                            if j > 0:
                                error_start = introns[j - 1].end
                            else:
                                error_start = tf.start
                            if j < len(introns) - 1:
                                error_end = introns[j + 1].start
                            else:
                                error_end = tf.end
                            self._add_error(i, transcript, error_start, error_end,
                                            self.is_plus_strand, types.OVERLAPPING_EXONS)
                            faulty_introns.append(intron)
                        # the case of a too short intron
                        # todo put the minimum length in a config somewhere
                        elif abs(intron.end - intron.start) < self.controller.config['min_intron_length']:
                            self._add_error(i, transcript, intron.start, intron.end,
                                            self.is_plus_strand, types.TOO_SHORT_INTRON)
                    # do not save faulty introns, the error should be descriptive enough
                    for intron in faulty_introns:
                        introns.remove(intron)

                    # finally, introns can be partial (although this normally happens at a sequence end)
                    for intron in transcript['introns']:
                        if intron.start == tf.start:
                            self._add_overlapping_error(i, transcript, intron, '5p', types.TRUNCATED_INTRON)
                        if intron.end == tf.end:
                            self._add_overlapping_error(i, transcript, intron, '3p', types.TRUNCATED_INTRON)
        # remove all errors that are in the wrong order (caused by overlapping super loci)
        # these can only be removed now as they were needed for further processing
        self._remove_backwards_errors()

    def _remove_backwards_errors(self):
        def is_backwards(e):
            return ((self.is_plus_strand and e.end < e.start)
                    or (not self.is_plus_strand and e.end > e.start))

        for group in self.groups:
            n_removed = 0
            for transcript in group['transcripts']:
                full_len = len(transcript['errors'])
                transcript['errors'] = [e for e in transcript['errors'] if not is_backwards(e)]
                n_removed += full_len - len(transcript['errors'])
            if n_removed > 0:
                self.controller.stats.backwards_errors_removed += n_removed
                msg = ('removed {count} backwards error(s) from overlapping super loci: '
                       'seqid: {seqid}, {geneid}').format(count=n_removed,
                                                          seqid=self.coord.seqid,
                                                          geneid=group['super_locus'].given_name)
                logger.debug(msg)

    def _add_error(self, i, transcript_g, start, end, is_plus_strand, error_type):
        error_i = FeatureImporter(self.coord,
                                  is_plus_strand,
                                  error_type,
                                  start=start,
                                  end=end,
                                  controller=self.controller)
        transcript_g['errors'].append(error_i)
        # a set, so a transcript hit twice by the same error type is still only counted once
        # when the types are tallied at the end of the import (see clean_and_insert)
        transcript_g['detected_error_types'].add(error_type)
        strand_str = 'plus' if is_plus_strand else 'minus'
        logger.debug(f'marked as erroneous: seqid: {self.coord.seqid}, {start}--{end}:'
                     f'{self.groups[i]["super_locus"].given_name}, on {strand_str} strand, '
                     f'with type: {error_type}')

    def _add_overlapping_error(self, i, transcript_g, handler, direction, error_type,
                               mark_other_handlers=None):
        """Constructs an error features that overlaps halfway to the next super locus
        in the given direction from the given handler if possible. Otherwise, mark until the end.
        If the direction is 'whole', the handler parameter is ignored.

        Also sets handler.start_is_biological_start=False (or the end) if necessary
        """
        if mark_other_handlers is None:
            mark_other_handlers = []

        assert direction in ['5p', '3p', 'whole']
        coord = self.groups[i]['super_locus'].coord
        # the error type is recorded as detected here, before any of the extent handling below
        # can shorten the range to nothing or drop it, so the import statistics report what was
        # found rather than what survived (see clean_and_insert)
        transcript_g['detected_error_types'].add(error_type)

        # set correct upstream error starting point
        if direction in ['5p', 'whole']:
            previous = self._exported_neighbour(i, -1)
            if previous is not None:
                anchor_5p = self._error_border_mark(self.groups[previous]['super_locus'],
                                                    self.groups[i]['super_locus'])
            else:
                if self.is_plus_strand:
                    anchor_5p = 0
                else:
                    anchor_5p = coord.length

        # set correct downstream error end point
        if direction in ['3p', 'whole']:
            following = self._exported_neighbour(i, 1)
            if following is not None:
                anchor_3p = self._error_border_mark(self.groups[i]['super_locus'],
                                                    self.groups[following]['super_locus'])
            else:
                if self.is_plus_strand:
                    anchor_3p = coord.length
                else:
                    anchor_3p = -1

        if direction == '5p':
            error_5p = anchor_5p
            error_3p = handler.start
            for h in [handler] + mark_other_handlers:
                if isinstance(h, FeatureImporter):
                    h.start_is_biological_start = False

        elif direction == '3p':
            error_5p = handler.end
            error_3p = anchor_3p
            for h in [handler] + mark_other_handlers:
                if isinstance(h, FeatureImporter):
                    h.end_is_biological_end = False

        elif direction == 'whole':
            error_5p = anchor_5p
            error_3p = anchor_3p

        # _error_border_mark's dist <= 0 case (the anchor's own neighbour already extends
        # into or past it, e.g. an overlapping predecessor) can place the computed anchor
        # past the handler's own, exact boundary; without clamping, that silently produces
        # a backwards region that _remove_backwards_errors later discards outright, losing
        # the error entirely instead of just its precise extent. Clamp the computed anchor
        # back to the reliable, handler-based side, not the other way around, so a
        # degenerate case collapses to a zero-length marker right at the handler's own
        # boundary instead of ballooning out to wherever the bad anchor landed. Such a marker
        # masks nothing, which is the correct amount when there is no unclaimed sequence left
        # to mask, while still recording the error (see docs/spec_vs_gff.md).
        sign = 1 if self.is_plus_strand else -1
        if sign * error_3p < sign * error_5p:
            if direction == '5p':
                error_5p = error_3p
            else:
                error_3p = error_5p

        if not self._zero_len_coords_at_sequence_edge(error_5p, error_3p, direction, coord):
            self._add_error(i, transcript_g, error_5p, error_3p, self.is_plus_strand, error_type)

    def _zero_len_coords_at_sequence_edge(self, error_5p, error_3p, direction, coordinate):
        """Check if error 5p-3p is of zero length due to hitting start or end of sequence"""
        out = False
        if self.is_plus_strand:
            if direction == '5p':
                if error_5p == error_3p == 0:
                    out = True
            elif direction == '3p':
                if error_5p == error_3p == coordinate.length:
                    out = True
        else:
            if direction == '5p':
                if error_5p == error_3p == coordinate.length - 1:
                    out = True
            elif direction == '3p':
                if error_5p == error_3p == -1:
                    out = True
        return out

    def _exported_neighbour(self, i, step):
        """The index of the nearest exported gene to one side, or None. A mask that runs part of
        the way to the next gene is measuring toward whatever a consumer will be handed, so
        unexported loci, the dropped partner of a resolved overlap among them, are skipped rather
        than cutting the mask short at sequence nothing is written for."""
        j = i + step
        while 0 <= j < len(self.groups):
            if self._exported_transcript(self.groups[j]) is not None:
                return j
            j += step
        return None

    @staticmethod
    def _border_offset(gap):
        """How far a mask reaches into a gap of unclaimed sequence,
        min(gap / 2, sqrt(gap) * 10). Nothing at all where there is no gap, as between a nested
        gene and the one around it."""
        if gap <= 0:
            return 0
        return min(gap // 2, int(math.sqrt(gap)) * 10)

    def _error_border_mark(self, sl, sl_next):
        """The point between two super loci that an error mask reaching toward the next gene
        stops at."""
        if self.is_plus_strand:
            return sl.end + self._border_offset(sl_next.start - sl.end)
        return sl.end - self._border_offset(sl.end - sl_next.start)


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
        except Exception as e:
            self.session.close()
            part_path = f'{self.database_path}.partial'
            shutil.move(self.database_path, part_path)
            logger.error(f'Aborting due to error, attempt so far saved at {part_path} '
                         f'for debugging purposes')
            raise e
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
                    # insert the errors
                    for error in transcript['errors']:
                        error.add_to_queue()
                        error.insert_feature_piece_association(tp.id)

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
                    # counted from the error types detected, not from the error features that
                    # ended up being inserted: an error whose masked range collapses to nothing
                    # (e.g. a missing UTR of a gene nested in another, which has no unclaimed
                    # sequence to extend the mask into) is still a real finding about the
                    # transcript, and dropping it from the count would understate the error rate
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
                 phase_3p=0,
                 phase=None,
                 is_truncated=False,
                 has_inframe_stop=False,
                 has_start_codon=True,
                 has_stop_codon=True,
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
        self.phase_3p = phase_3p  # only used for error checking
        # the phase actually saved to the db; defaults to phase_5p for every non-CDS feature
        # type (where phase is irrelevant and stays 0), but a CDS passes its own recomputed
        # value here instead of trusting the file's phase_5p (see _parse_gff_entries)
        self.phase = phase if phase is not None else phase_5p
        self.is_truncated = is_truncated  # only used for the TRUNCATED_CDS check
        self.has_inframe_stop = has_inframe_stop  # only used for the INFRAME_STOP_CODON check
        self.has_start_codon = has_start_codon  # only used for the MISSING_START_CODON check
        self.has_stop_codon = has_stop_codon  # only used for the MISSING_STOP_CODON check
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
