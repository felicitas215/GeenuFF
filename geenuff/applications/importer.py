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
KNOWN_TYPES = types.KNOWN_GFF_TYPES
IGNORABLE_TYPES = types.IGNORABLE_GFF_TYPES
GENE_TYPES = frozenset(x.value for x in types.SuperLocusAll)
TRANSCRIPT_TYPES = frozenset(x.value for x in types.TranscriptLevel)
# where an exon, CDS or UTR line goes among its transcript's lines, by type
PIECE_KEYS = {x.value: key for level, key in [(types.ExonLevel, 'exons'), (types.CDSLevel, 'cds'),
                                              (types.UTRLevel, 'utrs')] for x in level}
# above these shares of coding transcripts with CDS lines outside their exon lines, or with
# overlapping exon lines, the import warns
CDS_OUTSIDE_EXONS_WARN_SHARE = 0.05
OVERLAPPING_EXONS_WARN_SHARE = 0.05


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
    types.UNPLACEABLE_STRAND: 'their gene, transcript or its exon/CDS lines are not all on one '
                              'definite strand',
    types.UNPLACEABLE_COORDINATES: 'a line of theirs runs backwards, its start past its end',
    types.OUTSIDE_SEQUENCE: 'a line of theirs starts before or ends past their sequence',
    types.OVERLAP_DROPPED: 'dropped so an overlapping partner could be kept whole',
}


class MissingGeneLines(Exception):
    """Coding transcript lines name no Parent, as gene predictors write them without gene lines."""


# why a GFF line was left out while the lines are grouped into loci by their ID and Parent
# attributes (see OrganizedGFFEntries.load_organized_entries), for the import summary, where the
# GFF types of the lines stand beside it
DROPPED_LINE_REASONS = {
    'no_id': 'without an ID',
    'shared_id': 'sharing an ID with another gene or transcript line',
    'no_gene': 'naming no gene, without CDS lines',
    'several_genes': 'naming several genes',
    'parent_not_gene': 'naming a transcript as parent',
    'missing_parent_unplaceable': 'naming a missing parent, on different sequences or strands',
    'no_parent': 'naming no parent',
    'unknown_parent': 'naming a parent that matches no line',
    'gene_parented_duplicate': 'naming a gene with transcripts, duplicating their lines',
    'gene_parented': "naming a gene, duplicating none of its transcripts' lines",
    'gene_parented_cds': 'CDS naming a gene with transcripts',
    'floating_cds_unplaceable': 'CDS without a transcript line, spread over sequences or strands, not masked either',
    'other_sequence': 'on another sequence than their parent',
}


# what the masks of an exported transcript leave of it, for the import summary: nothing masked;
# only sequence beside it masked, for a missing UTR, for an overlapping gene dropped in its favor,
# or both; or the transcript itself masked in full, for floating CDS lines, for errors of its own,
# or for an overlap left unresolved alone (see exported_outcome)
EXPORTED_OUTCOMES = {
    'unmasked': 'nothing masked',
    'utr_flank': 'the flank of a missing UTR masked beside the gene',
    'overlap': 'the overhang of an overlapping gene dropped in their favor masked beside the gene',
    'utr_flank_and_overlap': 'both masked beside the gene, the two masks possibly overlapping',
    'floating': 'masked in full, created to mask CDS lines without a transcript line',
    'own_errors': 'masked in full for errors of their own',
    'nested': 'masked in full for an overlap with a gene inside or around them, left unresolved',
    'chain': 'masked in full for overlaps with more than one other gene, left unresolved',
}
MASKED_IN_FULL_OUTCOMES = ('floating', 'own_errors', 'nested', 'chain')
# the errors masking nothing or only sequence beside the transcript; every other one masks it in
# full (see GFFErrorHandling._find_errors)
NOT_MASKING_IN_FULL = frozenset({types.MISSING_UTR_5P, types.MISSING_UTR_3P, types.WRONG_PHASE_5P,
                                 types.SL_OVERLAP_ERROR})


def error_severity(errors) -> str:
    """What a transcript's own errors take away, from its error types, as the overlap decision
    grades it (see GFFErrorHandling._locus_severity): 'masked_outright', 'flank_masked' or 'none'.
    The masks of an overlap are not the transcript's own."""
    if errors - NOT_MASKING_IN_FULL:
        return 'masked_outright'
    if errors & {types.MISSING_UTR_5P, types.MISSING_UTR_3P}:
        return 'flank_masked'
    return 'none'


def overlap_length(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> int:
    """The bases two lists of sorted, merged, half-open (lo, hi) ranges share."""
    shared, i, j = 0, 0, 0
    while i < len(a) and j < len(b):
        shared += max(0, min(a[i][1], b[j][1]) - max(a[i][0], b[j][0]))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return shared


def flank_masked(side: str, span: tuple[int, int], masks: list[tuple[int, int]], is_plus_strand: bool) -> bool:
    """Whether a mask touches a transcript's ascending span on its 5p or 3p side, as the flank of a
    missing UTR does unless it came out empty for lack of room beside the gene."""
    if (side == '5p') == is_plus_strand:
        return any(lo < span[0] <= hi for lo, hi in masks)
    return any(lo <= span[1] < hi for lo, hi in masks)


def ascending_range(start: int, end: int, is_plus_strand: bool) -> tuple[int, int]:
    """The half-open (lo, hi) bases of a feature, its start and end running downwards on the minus
    strand."""
    return (start, end) if is_plus_strand else (end + 1, start + 1)


def exported_outcome(errors, transcript_feature, masks, nested: bool) -> str:
    """The EXPORTED_OUTCOMES key of an exported transcript, from its error types, transcript
    feature and geenuff_mask features; nested says whether an overlap left unresolved has one gene
    inside the other rather than a chain of more than two."""
    if not masks:
        return 'unmasked'
    if helpers.is_masked_whole(transcript_feature, masks):
        if types.FLOATING_CDS in errors:
            return 'floating'
        if errors - NOT_MASKING_IN_FULL:
            return 'own_errors'
        return 'nested' if nested else 'chain'
    # a mask not covering the transcript lies beside it, for a missing UTR or a gene dropped in its favor
    utr = bool(errors & {types.MISSING_UTR_5P, types.MISSING_UTR_3P})
    overlap = types.SL_OVERLAP_ERROR in errors
    if utr and overlap:
        return 'utr_flank_and_overlap'
    return 'overlap' if overlap else 'utr_flank'


# the steps of the import summary's gene funnel, from all genes imported down to those giving
# labels; each gene ends at one of them (see funnel_step)
EMPTY, LEFT_OUT, NON_CODING, NOT_EXPORTED, EXPORTED = 'empty', 'left_out', 'non_coding', 'not_exported', 'exported'
# the headers of the import summary's steps, in order (see ImportStatistics.summary_steps); the
# counts of the first two only exist while the GFF3 is read
(GFF_LINES, GFF_GENE_LINES, GENES_IMPORTED, GENES_WITH_TRANSCRIPTS, CODING_GENES, EXPORTED_CODING_GENES,
 GENES_GIVING_LABELS) = SUMMARY_STEPS = ('GFF lines', 'GFF gene lines', 'Genes imported', 'Genes with transcripts',
                                        'Coding genes', 'Exported coding genes', 'Genes giving labels')
READING_ONLY_STEPS = SUMMARY_STEPS[:2]
# below this share of the coding transcripts read giving labels, the import warns of little
# training data and names the scripts preparing a GFF3 file for import
LOW_TRAINING_DATA_SHARE = 0.1
PREPARING_SCRIPTS = ('add_gene_lines.py for transcripts without a gene or naming themselves as Parent',
                     'add_transcript_lines.py for genes whose exon and CDS lines name the gene itself '
                     '(only use on single isoform files!)',
                     'rebuild_exons.py for files listing exon lines for the UTRs only',
                     'rename_gff_seqids.py for GFF3 sequence IDs differing from the FASTA ones')
# why all transcripts of a gene were left out on import, keyed as in left_out_reason
LEFT_OUT_REASONS = {
    'no_strand': "on neither strand, e.g. NCBI's '?'",
    types.OUTSIDE_SEQUENCE: 'with a line outside their sequence',
    types.UNPLACEABLE_COORDINATES: 'with a line running backwards',
}


def plural(n: int, noun: str) -> str:
    """'1 transcript', '2 transcripts'."""
    return f'{n} {noun}{"" if n == 1 else "s"}'


def funnel_step(n_imported: int, n_dropped: int, n_coding: int, excluded: str | None) -> str:
    """The funnel step a gene ends at, from its transcripts imported and left out on import, its
    coding transcripts imported and why it is excluded from export, if it is."""
    if not n_imported:
        return LEFT_OUT if n_dropped else EMPTY
    if not n_coding:
        return NON_CODING
    return NOT_EXPORTED if excluded else EXPORTED


def left_out_reason(excluded: str | None) -> str:
    """The LEFT_OUT_REASONS key of a gene whose transcripts were all left out on import: the one
    stored as excluded_from_export, or 'no_strand' for a gene on neither strand, which the
    database does not record apart."""
    return excluded if excluded in (types.OUTSIDE_SEQUENCE, types.UNPLACEABLE_COORDINATES) else 'no_strand'


def format_summary(title: str, steps: list[tuple[str, int | None, str, list[tuple[str, int, str, int]], list[str]]],
                   notes: list[str]) -> str:
    """A summary as a funnel: one header per step with its count, below it what was subtracted on
    the way to the next step, so the counts add up from the GFF lines down to the training labels.

    steps holds (header, count or None, bracketed detail or '', [(sign, count, text, level)], extra
    lines); sign is '-', '+' or ' ' for a breakdown, and level indents an entry further. An entry
    with a count of 0 is left out, and so is a step without count or entries. The counts of a step
    are right-aligned to its longest one."""
    steps = [(header, count, detail, [e for e in entries if e[1]], extra)
             for header, count, detail, entries, extra in steps]
    steps = [step for step in steps if step[1] is not None or step[3]]
    lines = [f'{title}:']
    for header, count, detail, entries, extra in steps:
        width = max((len(str(n)) for _, n, _, _ in entries), default=1)
        lines.append(header + (f': {count}' if count is not None else ':') + (f' ({detail})' if detail else ''))
        lines += [f'  {"    " * level}{sign} {n:>{width}}  {text}' for sign, n, text, level in entries]
        lines += [f'  {line}' for line in extra]
    lines += [f'note: {note}' for note in notes]
    return '\n'.join(lines)


class ImportStatistics(object):
    """Aggregate counts for one add_genome() call. A single summary logged at the end of
    the import gives a full picture of what was found (see log_summary), instead of one
    log line per locus/transcript/error, which floods the log without being any easier to
    get an overview from."""

    def __init__(self) -> None:
        self.gff_lines_read: int = 0  # every feature line of the GFF3 file
        # lines of an SO type GeenuFF has no use for, e.g. region or intron, skipped silently, keyed by type
        self.ignored_feature_types: defaultdict[str, int] = defaultdict(int)
        self.gene_lines_read: int = 0
        self.total_super_loci: int = 0
        # the genes at each step of the summary's funnel (see funnel_step); the exported ones are
        # longest_transcripts. Genes whose transcripts were all left out on import are keyed by
        # left_out_reason, coding genes not exported by orm.SuperLocus.excluded_from_export
        self.empty_super_loci: int = 0
        self.genes_left_out_on_import: defaultdict[str, int] = defaultdict(int)
        self.non_coding_genes: int = 0
        # of those, the ones with CDS lines naming the gene, left out as gene_parented_cds
        self.non_coding_genes_with_gene_parented_cds: int = 0
        self.coding_genes_not_exported: defaultdict[str, int] = defaultdict(int)
        self.coding_transcripts_imported: int = 0  # of the coding genes, exported or not
        self.total_transcripts: int = 0  # every transcript read, imported or dropped
        self.imported_transcripts: int = 0  # of the above, the ones stored in the database
        self.total_coding_transcripts: int = 0
        # transcript.longest == True, one per exported gene; a gene left out of exports is not
        # counted, nothing about it having been checked either (see clean_and_insert)
        self.longest_transcripts: int = 0
        self.longest_error_free_transcripts: int = 0  # of the above, the ones with no error at all
        # of the above, by what their masks leave of them, keyed by an EXPORTED_OUTCOMES key
        self.exported_outcomes: Counter[str] = Counter()
        # the masks of the above, as ascending (lo, hi) ranges keyed by (sequence, is_plus_strand), each
        # masking only its own strand, and the bases of the genome's sequences on one strand
        self.mask_ranges: defaultdict[tuple, list[tuple[int, int]]] = defaultdict(list)
        # the spans of the above, keyed the same way, for the share of exported gene bases masked
        self.gene_spans: defaultdict[tuple, list[tuple[int, int]]] = defaultdict(list)
        self.genome_bases: int = 0
        # coding transcripts of exported genes other than the selected one, never exported
        self.unselected_coding_transcripts: int = 0
        # transcripts that CDS lines name, left out while grouping, keyed by the DROPPED_LINE_REASONS
        # key the first line of their chain was left out for, e.g. shared_id for a transcript whose
        # gene line shares its ID with another gene line
        self.coding_transcripts_dropped: defaultdict[str, int] = defaultdict(int)
        # of those, the ones below a gene line left out, and the IDs of those genes, gene lines sharing
        # one ID counting as one, keyed the same way
        self.coding_transcripts_dropped_under_genes: defaultdict[str, int] = defaultdict(int)
        self.coding_gene_ids_dropped: defaultdict[str, set[str]] = defaultdict(set)
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
        # coding transcripts with CDS lines outside their exon lines, warned about above a share of
        # CDS_OUTSIDE_EXONS_WARN_SHARE
        self.cds_outside_exons_transcripts: int = 0
        # coding transcripts with exon lines overlapping each other, warned about above a share of
        # OVERLAPPING_EXONS_WARN_SHARE
        self.overlapping_exons_transcripts: int = 0
        # CDS missing only their stop codon, as in GFF3 converted from GTF, extended by the stop codon
        # after them; GTF itself is not read
        self.stop_codons_recovered: int = 0
        # of those, the ones past the end of their exons, the last exon and the transcript extended
        self.stop_codons_recovered_past_exons: int = 0
        # GFF lines left out while grouping, keyed by a DROPPED_LINE_REASONS key and then by GFF
        # type; a line with several parents counts once per parent it is left out for. The lines
        # below a line left out are counted apart, under the reason the first line of their chain
        # was left out for
        self.dropped_lines: defaultdict[str, defaultdict[str, int]] = defaultdict(lambda: defaultdict(int))
        self.dropped_lines_below: defaultdict[str, defaultdict[str, int]] = defaultdict(lambda: defaultdict(int))
        # the IDs of genes with CDS lines naming them that have transcripts, and how many of those CDS
        # lines duplicate a CDS line of the gene's transcripts
        self.gene_parented_cds_genes: set[str] = set()
        self.gene_parented_cds_duplicates: int = 0
        # kept in full but never exported, counted per reason and keyed by the value stored in
        # orm.SuperLocus.excluded_from_export
        self.unexported_super_loci: defaultdict[str, int] = defaultdict(int)
        # transcripts left out entirely, no line of theirs having a storable range, or one lying
        # partly outside their sequence with no CDS line on it; a coding one is clipped instead and
        # masked whole, counted among the errors as beyond_sequence_edge
        self.unstorable_transcripts_dropped: int = 0
        self.transcripts_outside_sequence_dropped: int = 0
        self.unstranded_transcripts_dropped: int = 0  # the transcripts of genes on neither strand
        # every pair of super loci sharing genomic range, of whatever kind, as recorded in the
        # super_locus_overlap table; these mask nothing, see _record_overlap_pairs
        self.overlap_pairs_recorded: int = 0
        self.super_loci_in_overlap_pairs: int = 0
        # of the pairs where both genes are exported, how each was settled; the four always sum
        # to that number (see GFFErrorHandling._compute_overlap_masks)
        self.overlap_pairs_resolved: int = 0  # one locus kept whole, the other dropped
        self.overlap_pairs_both_masked_outright: int = 0  # neither worth keeping, both masked whole
        self.overlap_pairs_nested: int = 0  # one inside the other, never decided, both masked whole
        self.overlap_pairs_in_chains: int = 0  # a partner overlaps a third locus, so not decided
        self.overlap_loci_dropped: int = 0
        self.overlap_loci_dropped_labeled: int = 0  # of those, the ones not masked whole for their own errors
        self.overlap_loci_in_chains: int = 0  # each counted once however many pairs it is in
        # keyed by types.Errors value, counted per transcript selected for export from the error
        # types detected, so an error masking nothing still shows up here (see docs/spec_vs_gff.md)
        self.errors: defaultdict[str, int] = defaultdict(int)
        # the errors masking a transcript in full, counted over the transcripts masked in full
        self.errors_masking_in_full: Counter[str] = Counter()
        # the exported transcripts by their missing UTRs, 'none', '5p', '3p' or 'both'; the same for
        # the ones giving labels, with the missing UTRs whose flank mask came out empty and those
        # also with the overhang of an overlapping gene dropped in their favor masked
        self.utrs_missing: Counter[str] = Counter()
        self.labeled_utrs: Counter[str] = Counter()
        self.labeled_empty_utr_flanks: Counter[str] = Counter()
        self.labeled_with_overhang: int = 0
        self.unrecognized_feature_types: defaultdict[str, int] = defaultdict(int)  # keyed by the raw, unknown GFF type

    def add_exported(self, errors: set[str], outcome: str, key: tuple, span: tuple[int, int],
                     masks: list[tuple[int, int]]) -> None:
        """Counts one exported transcript: its error types and EXPORTED_OUTCOMES key, its missing
        UTRs and, where it gives labels, whether their flanks were masked; key is (sequence,
        is_plus_strand), span and masks are its ascending, half-open (lo, hi) ranges."""
        self.longest_transcripts += 1
        self.longest_error_free_transcripts += not errors
        for error_type in errors:
            self.errors[error_type] += 1
        self.exported_outcomes[outcome] += 1
        self.mask_ranges[key] += masks
        self.gene_spans[key].append(span)
        missing = [side for side, error_type in (('5p', types.MISSING_UTR_5P), ('3p', types.MISSING_UTR_3P))
                   if error_type in errors]
        utrs = 'both' if len(missing) == 2 else missing[0] if missing else 'none'
        # a transcript made to mask floating CDS lines is never checked for its UTRs
        if outcome != 'floating':
            self.utrs_missing[utrs] += 1
        if outcome in MASKED_IN_FULL_OUTCOMES:
            self.errors_masking_in_full.update(errors - NOT_MASKING_IN_FULL)
            return
        self.labeled_utrs[utrs] += 1
        self.labeled_with_overhang += types.SL_OVERLAP_ERROR in errors
        self.labeled_empty_utr_flanks[utrs] += sum(not flank_masked(side, span, masks, key[1]) for side in missing)

    @staticmethod
    def _per_type(counts: dict[str, int], n_shown: int = 5) -> str:
        """'type n, type n, ... and k more types', the most frequent first."""
        top = sorted(counts.items(), key=lambda item: -item[1])
        more = len(top) - n_shown
        return ', '.join(f'{name} {n}' for name, n in top[:n_shown]) + (f' and {more} more types' if more > 0 else '')

    def _dropped_lines_entry(self, reason: str, prefix: str = '', show_types: bool = True) -> tuple[int, str]:
        """One summary entry for the lines left out for one reason: their total and why, the GFF
        types (given show_types) and the lines below them left out with them in brackets, and the
        coding transcripts lost, being potential training data, where there are any."""
        counts, below = self.dropped_lines[reason], self.dropped_lines_below.get(reason)
        details = [self._per_type(counts)] if show_types else []
        if below:
            details.append(f'with {plural(sum(below.values()), "line")} below them: {self._per_type(below)}')
        if reason == 'gene_parented_cds':
            details.append(f'on {plural(len(self.gene_parented_cds_genes), "gene")}, '
                           f'{self.gene_parented_cds_duplicates} of them duplicating a CDS line of a transcript of '
                           f'their gene')
        text = prefix + DROPPED_LINE_REASONS[reason] + (f' ({"; ".join(details)})' if details else '')
        lost = self.coding_transcripts_dropped.get(reason, 0)
        under_genes = self.coding_transcripts_dropped_under_genes.get(reason, 0)
        gene_ids = plural(len(self.coding_gene_ids_dropped.get(reason, ())), 'gene ID')
        if lost:
            text += f', {plural(lost, "coding transcript")} lost'
        if under_genes:
            text += f', from {gene_ids}' if under_genes == lost else f', {under_genes} of them from {gene_ids}'
        return sum(counts.values()), text

    def labeled_genes(self) -> int:
        """The exported genes giving labels, i.e. not masked in full."""
        return self.longest_transcripts - sum(self.exported_outcomes[outcome] for outcome in MASKED_IN_FULL_OUTCOMES)

    def no_training_data_reason(self) -> str:
        """The largest loss of coding transcripts or genes in the summary's funnel, or the lack of
        CDS lines where nothing coding was lost."""
        losses = [(n, f'{plural(n, "coding transcript")} left out while grouping lines into genes, '
                      f'{DROPPED_LINE_REASONS[reason]}') for reason, n in self.coding_transcripts_dropped.items()]
        losses += [(n, f'{plural(n, "gene")} with all transcripts left out on import, {LEFT_OUT_REASONS[reason]}')
                   for reason, n in self.genes_left_out_on_import.items()]
        losses += [(n, f'{plural(n, "coding gene")} not exported, {UNEXPORTED_REASONS[reason]}')
                   for reason, n in self.coding_genes_not_exported.items()]
        top_errors = ', '.join(error_type for error_type, _ in self.errors_masking_in_full.most_common(3))
        for outcome in MASKED_IN_FULL_OUTCOMES:
            n = self.exported_outcomes[outcome]
            text = f'{plural(n, "exported gene")} {EXPORTED_OUTCOMES[outcome]}'
            losses.append((n, text + (f', most often {top_errors}' if outcome == 'own_errors' and top_errors else '')))
        n, reason = max(losses, key=lambda loss: loss[0], default=(0, ''))
        return f'the largest loss is {reason}' if n else 'the GFF3 file has no CDS lines'

    def log_summary(self, species: str) -> None:
        """Logs one summary of the import (see format_summary), and at DEBUG level the complete
        lists of the GFF types skipped, of which the summary names the most frequent only."""
        for label, counts in (('SO types GeenuFF does not use', self.ignored_feature_types),
                              ('unknown types', self.unrecognized_feature_types)):
            if counts:
                logger.debug(f'GFF lines skipped, {label}: {self._per_type(counts, len(counts))}')
        logger.info(format_summary(f'import summary for "{species}"', self.summary_steps(),
                                   self.summary_notes()))

    @staticmethod
    def _error_lines(heading: str, errors: Counter[str]) -> list[str]:
        """A heading and below it the error types, the most frequent first."""
        if not errors:
            return []
        width = len(str(max(errors.values())))
        return [heading] + [f'  {n:>{width}}  {error_type}' for error_type, n in errors.most_common()]

    def summary_steps(self, reading_counts: bool = True
                      ) -> list[tuple[str, int | None, str, list[tuple[str, int, str, int]], list[str]]]:
        """The steps of the import summary's gene funnel (see format_summary), from the GFF lines
        down to the genes giving labels; each step's count less what is subtracted below it is the
        next step's count. reading_counts is False where the counts that only exist while the GFF3
        is read are unknown, as for a summary rebuilt from the database: the GFF line steps are left
        out then."""
        steps = []
        if reading_counts:
            n_dropped = sum(sum(c.values()) for c in self.dropped_lines.values())
            n_dropped += sum(sum(c.values()) for c in self.dropped_lines_below.values())
            reasons = [reason for reason in DROPPED_LINE_REASONS if reason in self.dropped_lines]
            # a reason that left out gene lines alone is given once, under the gene lines, with what
            # was left out below them and the coding transcripts lost through them
            gene_only = [reason for reason in reasons if set(self.dropped_lines[reason]) <= GENE_TYPES]
            n_coding_lost = sum(n for reason, n in self.coding_transcripts_dropped.items() if reason not in gene_only)
            n_under_genes = sum(sum(self.dropped_lines[reason].values())
                                + sum(self.dropped_lines_below.get(reason, {}).values()) for reason in gene_only)
            dropped_text = 'left out while grouping lines into genes by ID and Parent'
            if n_under_genes:
                dropped_text += f', {n_under_genes} of them gene lines and the lines below them (see {GFF_GENE_LINES})'
            if n_coding_lost:
                dropped_text += f', {plural(n_coding_lost, "coding transcript")} lost'
            steps.append((GFF_LINES, self.gff_lines_read, '', [
                ('-', sum(self.ignored_feature_types.values()),
                 f'skipped, SO types GeenuFF does not use ({self._per_type(self.ignored_feature_types)})', 0),
                ('-', sum(self.unrecognized_feature_types.values()),
                 f'skipped, unknown types ({self._per_type(self.unrecognized_feature_types)})', 0),
                ('-', n_dropped, dropped_text, 0),
            ] + [(' ', *self._dropped_lines_entry(reason), 1) for reason in reasons if reason not in gene_only],
                          []))
            steps.append((GFF_GENE_LINES, self.gene_lines_read, '', [
                ('-', *self._dropped_lines_entry(reason, prefix='left out, ', show_types=False), 0)
                if reason in gene_only else
                ('-', sum(n for gff_type, n in self.dropped_lines[reason].items() if gff_type in GENE_TYPES),
                 f'left out, {DROPPED_LINE_REASONS[reason]} (see {GFF_LINES})', 0)
                for reason in reasons
            ] + [
                ('+', self.genes_inferred_for_missing_parents,
                 'created for transcripts naming a Parent that matches no line', 0),
                ('+', self.floating_cds_genes, 'created to mask CDS lines without a transcript line', 0),
            ], []))

        n_with_transcripts = (self.total_super_loci - self.empty_super_loci
                              - sum(self.genes_left_out_on_import.values()))
        n_exported = self.longest_transcripts
        n_transcripts_left_out = self.total_transcripts - self.imported_transcripts
        steps.append((GENES_IMPORTED, self.total_super_loci, '', [
            ('-', self.empty_super_loci, 'without a transcript under them', 0),
        ] + [('-', n, f'{LEFT_OUT_REASONS[reason]}: all their transcripts left out', 0)
             for reason, n in sorted(self.genes_left_out_on_import.items())], []))

        transcripts_detail = plural(self.imported_transcripts, 'transcript')
        if n_transcripts_left_out:
            transcripts_detail += f', {n_transcripts_left_out} more left out on import'
        non_coding_text = 'non-coding, no CDS line(s)'
        if self.non_coding_genes_with_gene_parented_cds:
            non_coding_text += (f', {self.non_coding_genes_with_gene_parented_cds} of them with CDS lines naming the '
                                f'gene instead (use add_transcript_lines.py to fix)')
        steps.append((GENES_WITH_TRANSCRIPTS, n_with_transcripts, transcripts_detail,
                      [('-', self.non_coding_genes, non_coding_text, 0)], []))

        coding_detail = [plural(self.coding_transcripts_imported, 'coding transcript')]
        if self.transcripts_with_exons_built:
            coding_detail.append(f'exons built from CDS & UTRs for {self.transcripts_with_exons_built}')
        if self.stop_codons_recovered:
            coding_detail.append(f'stop codon found & added to {self.stop_codons_recovered}')
        not_exported = [('-', self.coding_genes_not_exported.get(types.OVERLAP_DROPPED, 0),
                         f'not exported, dropped in favor of an overlapping gene '
                         f'({self.overlap_loci_dropped_labeled} would have been labeled in full)', 0)]
        not_exported += [('-', n, f'not exported, {UNEXPORTED_REASONS[reason]}', 0)
                         for reason, n in sorted(self.coding_genes_not_exported.items())
                         if reason != types.OVERLAP_DROPPED]
        steps.append((CODING_GENES, n_with_transcripts - self.non_coding_genes, '; '.join(coding_detail),
                      not_exported, []))

        exported_detail = 'one transcript each, the longest coding one'
        if self.unselected_coding_transcripts:
            exported_detail += (f'; {plural(self.unselected_coding_transcripts, "further coding isoform")} '
                                f'not exported')
        utrs = self.utrs_missing
        utr_lines = [(utrs['none'], 'both UTRs annotated'), (utrs['5p'], "only the 5' UTR missing"),
                     (utrs['3p'], "only the 3' UTR missing"), (utrs['both'], 'both UTRs missing')]
        utr_width = len(str(max(n for n, _ in utr_lines)))
        n_floating = self.exported_outcomes['floating']
        not_counted = f'; the {n_floating} created to mask CDS lines not counted' if n_floating else ''
        utr_block = [f"UTRs of these, each gene in one line (5' UTR missing in "
                     f"{self.errors.get(types.MISSING_UTR_5P, 0)}, 3' UTR missing in "
                     f"{self.errors.get(types.MISSING_UTR_3P, 0)}; can overlap with other "
                     f"errors{not_counted}):"] * bool(sum(utrs.values()))
        utr_block += [f'  {n:>{utr_width}}  {text}' for n, text in utr_lines if n]
        n_wrong_phase = self.errors.get(types.WRONG_PHASE_5P, 0)
        steps.append((EXPORTED_CODING_GENES, n_exported, exported_detail,
                      [('-', self.exported_outcomes[outcome], EXPORTED_OUTCOMES[outcome], 0)
                       for outcome in MASKED_IN_FULL_OUTCOMES],
                      self._error_lines('errors masking them in full, each type counted once per transcript '
                                        '(can overlap):',
                                        self.errors_masking_in_full)
                      + utr_block
                      + [f'first CDS phase in the file not 0, masking nothing: {n_wrong_phase}'] * bool(n_wrong_phase)))

        # overlapping masks, and the spans of overlapping exported genes, on one strand count once
        merged_masks = {key: GFFErrorHandling._merge_ranges(ranges) for key, ranges in self.mask_ranges.items()}
        merged_spans = {key: GFFErrorHandling._merge_ranges(ranges) for key, ranges in self.gene_spans.items()}
        n_masked_bases = sum(hi - lo for ranges in merged_masks.values() for lo, hi in ranges)
        n_gene_bases = sum(hi - lo for ranges in merged_spans.values() for lo, hi in ranges)
        n_masked_gene_bases = sum(overlap_length(spans, merged_masks.get(key, []))
                                  for key, spans in merged_spans.items())
        both_strands = 2 * self.genome_bases
        n_labeled = self.labeled_genes()
        labeled, empty = self.labeled_utrs, self.labeled_empty_utr_flanks
        n_one_empty = empty['5p'] + empty['3p']
        one_detail = f"{labeled['5p']} 5', {labeled['3p']} 3'" + (f'; {n_one_empty} of these masks empty, no room '
                                                                  f'beside the gene' if n_one_empty else '')
        both_detail = f' ({empty["both"]} of these masks empty, no room beside the gene)' if empty['both'] else ''
        steps.append((GENES_GIVING_LABELS, n_labeled, 'each labeled in full, masks only beside it', [
            (' ', labeled['none'], 'labeled in full, no UTR flank masked', 0),
            (' ', labeled['both'], f'the flanks of both UTRs masked beside the gene{both_detail}', 0),
            (' ', labeled['5p'] + labeled['3p'], f'the flank of one UTR masked beside the gene ({one_detail})', 0),
        ], [f'{self.labeled_with_overhang} of them also with the overhang of an overlapping gene dropped in their '
            f'favor masked beside them'] * bool(self.labeled_with_overhang)
                      + [f'{n_masked_bases} of the genome\'s {both_strands} bases on both strands masked '
                         f'({100 * n_masked_bases / both_strands if both_strands else 0:.2f}%), each mask on its '
                         f'gene\'s strand only; Helixer gets no labels there',
                         f'{n_masked_gene_bases} of the {n_gene_bases} bases of the exported genes masked '
                         f'({100 * n_masked_gene_bases / n_gene_bases if n_gene_bases else 0:.2f}%), union of '
                         f'overlapping genes spans used; missing UTRs are of unknown length, not counted']))
        return steps

    def summary_notes(self) -> list[str]:
        """Explanations for the parts of the summary above that a count alone does not convey."""
        overlap_pairs = (self.overlap_pairs_resolved + self.overlap_pairs_both_masked_outright
                         + self.overlap_pairs_nested + self.overlap_pairs_in_chains)
        notes = []
        if overlap_pairs:
            notes.append(f'{overlap_pairs} pairs of exported genes overlap on one strand: '
                         f'{self.overlap_pairs_resolved} settled by dropping one gene, the rest masked in full '
                         f'(docs/overlap_masking.md)')
        if self.errors:
            notes.append('errors count as found, also where their mask came out empty, e.g. a missing UTR '
                         'without room beside the gene')
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
    """Stores the handler objects for a super locus in an organized fashion.
    The format is similar to the one of OrganizedGFFEntryGroup, but it stores objects
    according to the Geenuff way of saving genomic annotations. This format can then
    be checked for errors and changed accordingly before being inserted into the db.

    The importers are organized in the following way:

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

    def _drop_unstorable_transcript(self, sl_i, t_id, what, gff_start, gff_end, coding):
        """Leaves a transcript out entirely because one of its lines has no range it could be
        stored as, its start running past its end. Trans-spliced organellar genes wrapping a
        circular molecule's origin are written this way; nothing can be made of them here, and
        the locus they belong to is kept out of exports."""
        sl_i.excluded_from_export = types.UNPLACEABLE_COORDINATES
        sl_i.n_dropped_transcripts += 1
        sl_i.n_dropped_coding_transcripts += coding
        self.controller.stats.unstorable_transcripts_dropped += 1
        logger.warning(f"dropping transcript '{t_id}': its {what} runs from {gff_start} to "
                       f"{gff_end}, i.e. backwards, so it has no storable range; super locus "
                       f"'{sl_i.given_name}' will not be exported")

    @staticmethod
    def _cds_extended_by_a_codon(t_entries, is_plus_strand, seq_length):
        """Copies of a transcript's CDS lines extended by the 3 bases following them along its
        exons, possibly across an intron, and directly past its last exon where fewer are left in
        the exons, with how many of the 3 lie past them; None where the sequence ends first. Used
        to recover the stop codon GTF leaves out of the CDS, in GFF3 converted from GTF, which can
        also lack the exon part holding it; GTF itself is not read."""
        cds, exons = t_entries['cds'], t_entries['exons']
        if is_plus_strand:
            in_exons = (p for exon in exons for p in range(max(exon.start, cds[-1].end + 1), exon.end + 1))
            past_exons = range(max(line.end for line in exons + cds) + 1, seq_length + 1)
        else:
            in_exons = (p for exon in reversed(exons) for p in range(min(exon.end, cds[0].start - 1),
                                                                    exon.start - 1, -1))
            past_exons = range(min(line.start for line in exons + cds) - 1, 0, -1)
        following = list(itertools.islice(in_exons, 3))
        # overlapping exons repeat a base
        if len(set(following)) < len(following):
            return None
        n_past = 3 - len(following)
        following += itertools.islice(past_exons, n_past)
        # the sequence ends first
        if len(following) < 3:
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
        return sorted(pieces, key=lambda c: c.start), n_past

    @staticmethod
    def _extend_to_cds(t, t_entries, is_plus_strand):
        """Extends a transcript line and a copy of its outermost exon at the 3' end to the end of
        its CDS, for a stop codon added past them."""
        exons = t_entries['exons']
        if is_plus_strand:
            last = copy.copy(exons[-1])
            last.end = t_entries['cds'][-1].end
            t_entries['exons'] = exons[:-1] + [last]
            t.end = max(t.end, last.end)
        else:
            first = copy.copy(exons[0])
            first.start = t_entries['cds'][0].start
            t_entries['exons'] = [first] + exons[1:]
            t.start = min(t.start, first.start)

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

    def _drop_transcript_outside_sequence(self, sl_i, t_id, line, coding):
        """Leaves a transcript with no CDS line on its sequence out entirely because one of its
        lines starts before its sequence or ends past it, which an annotation of another assembly
        version, or a gene crossing the origin of a circular molecule written with an end past
        its length, both produce. Such a line cannot be checked or labeled against the sequence,
        and one starting before it cannot be stored at all; the locus it belongs to is kept out
        of exports."""
        sl_i.excluded_from_export = types.OUTSIDE_SEQUENCE
        sl_i.n_dropped_transcripts += 1
        sl_i.n_dropped_coding_transcripts += coding
        self.controller.stats.transcripts_outside_sequence_dropped += 1
        logger.warning(f"dropping transcript '{t_id}': its {line.type} runs from {line.start} to "
                       f"{line.end}, outside its sequence of {self.coord.length} bp; super locus "
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
            logger.debug(f"super locus '{sl.get_ID()}' has strand '{sl.strand}' (not '+' "
                         f"or '-'), most likely a trans-spliced gene; saving it without "
                         f"any transcripts so it is never exported or masked")
            n_coding = sum(bool(t_entries['cds']) for t_entries in entries['transcripts'].values())
            sl_i = self.importers['super_locus'] = SuperLocusImporter(
                entry_type=sl.type,
                given_name=sl.get_ID(),
                coord=self.coord,
                is_plus_strand=None,
                start=sl.start,
                end=sl.end,
                excluded_from_export=types.UNPLACEABLE_STRAND if n_coding else None,
                controller=self.controller)
            sl_i.n_dropped_transcripts = len(entries['transcripts'])
            sl_i.n_dropped_coding_transcripts = n_coding
            stats.total_transcripts += len(entries['transcripts'])
            stats.unstranded_transcripts_dropped += len(entries['transcripts'])
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
            # masked whole without further checks, see GFFErrorHandling._find_errors
            t_importers = {'masks': [], 'detected_error_types': set(),
                           'mask_whole': types.FLOATING_CDS if t_entries['floating'] else None}
            t_id = t.get_ID()
            if t.start > t.end:
                self._drop_unstorable_transcript(sl_i, t_id, 'transcript', t.start, t.end, bool(t_entries['cds']))
                continue
            # GFF coordinates are 1-based and inclusive, so a line lies within its sequence from
            # 1 up to the sequence length
            outside = [line for line in [t] + t_entries['exons'] + t_entries['cds']
                       if line.start < 1 or line.end > self.coord.length]
            if outside:
                if not self._clip_to_sequence(t, t_entries):
                    self._drop_transcript_outside_sequence(sl_i, t_id, outside[0], bool(t_entries['cds']))
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
                # a coding transcript that is not on one definite strand, or whose gene or CDS
                # lines disagree with it, cannot be analyzed; the locus is left out of exports
                # while its features stay in the database (see docs/trans_splicing.md for why
                # nothing is masked in its place either)
                if (unusable_strand or t_is_plus_strand != sl_is_plus_strand
                        or {strand_or_none(x) for x in t_entries['cds']} != {t_is_plus_strand}):
                    sl_i.excluded_from_export = types.UNPLACEABLE_STRAND
                    logger.debug(f"coding transcript '{t_id}' is not on one definite strand, that of its "
                                 f"gene; super locus '{sl.get_ID()}' will not be exported")
                # create protein handler
                # floating CDS lines masked together can name several proteins or none
                protein_id = t_id if t_entries['floating'] else self._protein_id(t_entries['cds']) or t_id
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
                # translation ends there whatever the file says. Where the exons end before it,
                # the last exon and the transcript are extended to hold it. The other checks above
                # hold for the extended CDS as well
                if not (t_importers['mask_whole'] or is_truncated or has_inframe_stop or has_stop_codon):
                    extended = self._cds_extended_by_a_codon(t_entries, t_is_plus_strand, self.coord.length)
                    if extended is not None:
                        extended_cds, n_past_exons = extended
                        extended_seq = spliced_cds_sequence(self.coord.sequence, extended_cds, t_is_plus_strand)
                        if extended_seq[-3:] in STOP_CODONS:
                            t_entries['cds'], cds_seq, has_stop_codon = extended_cds, extended_seq, True
                            stats.stop_codons_recovered += 1
                            if n_past_exons:
                                self._extend_to_cds(t, t_entries, t_is_plus_strand)
                                tf_i.set_start_end_from_gff(t.start, t.end)
                                stats.stop_codons_recovered_past_exons += 1
                # the lines are sorted by start, so any overlap shows between neighbors
                cds_lines = t_entries['cds']
                overlapping_cds = any(b.start <= a.end for a, b in zip(cds_lines, cds_lines[1:]))
                cds_outside_exons = self._cds_outside_exons(cds_lines, t_entries['exons'])
                stats.cds_outside_exons_transcripts += cds_outside_exons
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
                                        is_outside_exons=cds_outside_exons,
                                        spliced_length=len(cds_seq),
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
                    self._drop_unstorable_transcript(sl_i, t_id, 'CDS', gff_cds_start, gff_cds_end, True)
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
                has_overlapping_exons = False
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
                        has_overlapping_exons = True
                        logger.debug(f'{t_id} on {exon.seqid}: exon {exon.start}-{exon.end} overlaps '
                                     f'{plural(len(overlapping), "other exon")} of its transcript')
                        overlapper_begin = min([o.begin for o in overlapping])
                        overlapper_end = max([o.end for o in overlapping])
                        ovlp_end = max(overlapper_begin, inv_e_start)
                        ovlp_start = min(overlapper_end, inv_e_end)
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
                stats.overlapping_exons_transcripts += has_overlapping_exons

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
            stats.imported_transcripts += 1
            stats.floating_cds_transcripts += t_entries['floating']
            stats.transcripts_with_exons_built += t_entries['exons_built']
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
    def _cds_outside_exons(cds_lines, exon_lines):
        """Whether a CDS line reaches outside the exon lines, touching or overlapping exons counting as
        one; both lists are GFF lines sorted by start."""
        covered = []
        for exon in exon_lines:
            if covered and exon.start <= covered[-1][1] + 1:
                covered[-1][1] = max(covered[-1][1], exon.end)
            else:
                covered.append([exon.start, exon.end])
        return any(not any(lo <= cds.start and cds.end <= hi for lo, hi in covered) for cds in cds_lines)

    @staticmethod
    def _coding_length(cds, introns):
        """Spliced length of a CDS: its span less the parts of the transcript's introns within it.
        The backwards introns standing in for overlapping exons are no introns and are skipped."""
        c_lo, c_hi = sorted((cds.start, cds.end))
        length = c_hi - c_lo
        for intron in introns:
            if (intron.end < intron.start) if intron.is_plus_strand else (intron.end > intron.start):
                continue
            lo, hi = sorted((intron.start, intron.end))
            length -= max(0, min(c_hi, hi) - max(c_lo, lo))
        return length

    @staticmethod
    def _protein_id(cds_entry_list):
        """The protein ID of a transcript's CDS lines: their protein_id attribute, else their
        Parent, several joined by commas, so a CDS line shared by transcripts gives each the same;
        None without either. Only stored, never used to decide anything."""
        protein_ids = set()
        for cds_entry in cds_entry_list:
            attributes = cds_entry.attrib_filter(tag='protein_id')
            protein_ids.update(attributes[0].value if attributes else cds_entry.get_Parent() or [])
        return ','.join(sorted(protein_ids)) or None


class OrganizedGFFEntryGroup(object):
    """Holds the entries of one super locus as grouped by OrganizedGFFEntries, and returns the
    corresponding OrganizedGeenuffImporterGroup. Does not perform error checking, which happens
    later.

    The entries are organized in the following way, exons, CDS and UTRs sorted by start:

    entries = {
        'super_locus' = super_locus_entry,
        'transcripts' = {
            transcript_entry1: {
                'exons': [ordered_exon_entry1, ordered_exon_entry2, ...],
                'cds': [ordered_cds_entry1, ordered_cds_entry2, ...],
                'utrs': [ordered_utr_entry1, ordered_utr_entry2, ...],
                'exons_built': False,  # True where the exons were built from CDS and UTR lines, not floating
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
        self.stats.gene_lines_read += len(genes)

        # {ID of a gene or transcript line left out: the reason the first line of its chain was left
        # out for}, so that the lines below it go too, and what they lose is put down to that reason
        dropped_ids = {}
        genes_by_id = self._lines_by_unique_id(genes, transcripts, dropped_ids)
        transcripts_by_id = self._lines_by_unique_id(transcripts, genes, dropped_ids)
        groups = {gene_id: {'super_locus': gene, 'transcripts': {}}
                  for gene_id, gene in genes_by_id.items()}
        cds_parents = {parent for piece, key in pieces if key == 'cds' for parent in piece.get_Parent() or []}
        self._place_transcripts(transcripts_by_id, groups, dropped_ids, cds_parents)
        self._place_pieces(pieces, groups, dropped_ids)

        self.organized_entries = {}
        for group in groups.values():
            for t, t_entries in group['transcripts'].items():
                for key in ['exons', 'cds', 'utrs']:
                    t_entries[key].sort(key=lambda e: e.start)
                if t_entries['cds'] and not t_entries['exons']:
                    t_entries['exons'] = self._built_exons(t, t_entries)
                    # floating CDS lines get exons too, but are counted as floating only
                    t_entries['exons_built'] = not t_entries['floating']
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

    def _drop_below(self, entry, reason):
        """Leaves out a line whose parent was left out, counted under the reason the first line of
        its chain was left out for."""
        self.stats.dropped_lines_below[reason][entry.type] += 1
        logger.debug(f'leaving out {entry.type} at {entry.seqid}:{entry.start}-{entry.end}, its '
                     f'parent being left out: {DROPPED_LINE_REASONS[reason]}')

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
                dropped_ids[line_id] = ('shared_id', line_id if line.type in GENE_TYPES else None)
            else:
                by_id[line_id] = line
        return by_id

    def _place_transcripts(self, transcripts_by_id, groups, dropped_ids, cds_parents):
        """Puts every transcript under the one gene it names. The transcripts naming one Parent ID
        that matches no line share a gene inferred for them, spanning them all (see _inferred_gene),
        unless they lie on different sequences or strands, which leaves them all out. A transcript
        naming no Parent at all is left out where no CDS line names it (cds_parents), giving no
        training data either way.

        Raises:
            MissingGeneLines: if coding transcripts name no Parent at all, which gene predictors
                write without gene lines; add_gene_lines.py adds them
        """
        parentless = [t_id for t_id, t in transcripts_by_id.items() if not t.get_Parent()]
        for t_id in [t_id for t_id in parentless if t_id not in cds_parents]:
            self._drop(transcripts_by_id.pop(t_id), 'no_gene')
            dropped_ids[t_id] = ('no_gene', None)
        parentless = [t_id for t_id in parentless if t_id in cds_parents]
        if parentless:
            # where they lie tells a few strays, e.g. on an organelle sequence, from a whole file without gene lines
            by_seqid = Counter(transcripts_by_id[t_id].seqid for t_id in parentless).most_common()
            on = ', '.join(f'{seqid}: {n}' for seqid, n in by_seqid[:3]) + (', ...' if len(by_seqid) > 3 else '')
            n_coding = sum(t_id in cds_parents for t_id in transcripts_by_id)
            share = len(parentless) / n_coding
            raise MissingGeneLines(
                f'{len(parentless)} of {n_coding} coding transcript lines '
                f'({f"{share:.1%}" if share >= 0.001 else "<0.1%"}) name no Parent gene or "Parent=" attribute, on '
                f'{len(by_seqid)} sequence{"s" if len(by_seqid) > 1 else ""} ({on}), e.g. '
                f'{", ".join(parentless[:3])}, as gene predictors such as '
                f'AUGUSTUS, BRAKER and Tiberius write them without gene lines. Add the gene lines with '
                f'add_gene_lines.py and import the GFF3 file it writes.')
        named_missing_parent = defaultdict(list)
        for t_id, t in transcripts_by_id.items():
            parents = t.get_Parent() or []
            if len(parents) == 1 and parents[0] in groups:
                reason = 'other_sequence' if t.seqid != groups[parents[0]]['super_locus'].seqid else None
            elif len(parents) > 1:
                reason = 'several_genes'
            elif parents[0] in dropped_ids:
                self._drop_below(t, dropped_ids[parents[0]][0])
                dropped_ids[t_id] = dropped_ids[parents[0]]
                continue
            elif parents[0] in transcripts_by_id:
                reason = 'parent_not_gene'
            else:
                named_missing_parent[parents[0]].append(t)
                continue
            if reason is None:
                groups[parents[0]]['transcripts'][t] = self._new_t_entries()
            else:
                self._drop(t, reason)
                dropped_ids[t_id] = (reason, None)

        for parent_id, members in named_missing_parent.items():
            if len({(t.seqid, t.strand) for t in members}) > 1:
                for t in members:
                    self._drop(t, 'missing_parent_unplaceable')
                    dropped_ids[t.get_ID()] = ('missing_parent_unplaceable', None)
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
        # {ID of a left out line that CDS lines name: (the reason its chain was left out for, the
        # gene ID the chain starts from, or None where it starts from a transcript line)}
        coding_dropped = {}
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
                    self._drop_below(piece, dropped_ids[parent_id][0])
                    if key == 'cds':
                        coding_dropped[parent_id] = dropped_ids[parent_id]
                elif key == 'cds':
                    floating[('unknown_parent', parent_id)].append(piece)
                else:
                    self._drop(piece, 'unknown_parent')
        for reason, gene_id in coding_dropped.values():
            self.stats.coding_transcripts_dropped[reason] += 1
            if gene_id is not None:
                self.stats.coding_gene_ids_dropped[reason].add(gene_id)
                self.stats.coding_transcripts_dropped_under_genes[reason] += 1

        # judged only once every transcript line is placed, so the order of lines does not matter
        genes_with_cds = set()
        for piece, key, gene_id in gene_parented:
            position = (piece.start, piece.end)
            duplicate = any(position == (other.start, other.end)
                            for t_entries in groups[gene_id]['transcripts'].values()
                            for other in t_entries[key])
            if key == 'cds':
                # counted per gene, its CDS lines together perhaps being an isoform
                self._drop(piece, 'gene_parented_cds')
                genes_with_cds.add(gene_id)
                self.stats.gene_parented_cds_duplicates += duplicate
            else:
                self._drop(piece, 'gene_parented_duplicate' if duplicate else 'gene_parented')
        self.stats.gene_parented_cds_genes |= genes_with_cds
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
            # a CDS ID can equal the parent ID of other floating CDS lines; both go under one gene
            if gene is None and gene_id not in groups:
                groups[gene_id] = {'super_locus': self._inferred_gene(gene_id, [t]), 'transcripts': {}}
                self.stats.floating_cds_genes += 1
            t_entries = groups[gene_id]['transcripts'][t] = self._new_t_entries(floating=True)
            t_entries['cds'] = lines

    def _useful_gff_entries(self):
        reader = self._gff_gen()
        for entry in reader:
            if entry.type in IGNORABLE_TYPES:
                self.stats.ignored_feature_types[entry.type] += 1
            else:
                yield entry

    def _gff_gen(self):
        for line in helpers.gff_feature_lines(self.gff_file):
            self.stats.gff_lines_read += 1
            entry = gffhelper.GFFObject(gffline=line)
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
                for k in (i, j):
                    self.groups[k]['unresolved_overlap'] = 'nested'
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
                    # given up although it would have given labels, i.e. what the overlap cost
                    stats.overlap_loci_dropped_labeled += severity[dropped] != MASKED_OUTRIGHT
        stats.overlap_loci_in_chains += len(in_chains)
        for k in in_chains:
            self.groups[k]['unresolved_overlap'] = 'chain'
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
        # kept sign-normalized, to be merged with the ranges of the locus' own errors (see
        # resolve_errors)
        self._overlap_masks = masks

    def _buffered_span(self, i):
        """The span of locus i's selected transcript with the flank on both sides, sign-normalized.

        Each flank reaches into the gap toward the closest edge of another locus' selected
        transcript, i.e. the largest end at or before the span's start and the smallest start at
        or after its end, or toward the end of the sequence where there is none. Taking the
        closest edge rather than the next locus in sort order matters where a long transcript
        encloses a shorter one: the shorter one need not end closest to whatever follows the long
        one. A transcript overlapping locus i is not a neighbor, sharing sequence with it rather
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
        only mask a flank still has its coding sequence labeled (see _find_errors)."""
        extents = [f.extent for f in GFFErrorHandling._selected_transcript(group)['findings']]
        if any(e == 'whole' or isinstance(e, tuple) and e[0] != e[1] for e in extents):
            return MASKED_OUTRIGHT
        if any(e in ('5p', '3p') for e in extents):
            return FLANK_MASKED
        return NOT_MASKED

    @staticmethod
    def _coding_extent(group, sign):
        """Everything a locus' coding transcripts reach over, as a sign-normalized ascending
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
        # overlapping exon lines (a backwards intron standing in for them), overlapping CDS lines or
        # CDS lines outside the exons (e.g. exon lines for the UTRs only) most likely come from a
        # wrong annotation, so the locus is masked whole together with the flank on both sides.
        # Nothing else is checked: the spliced CDS sequence the checks below rely on repeats the
        # overlapping bases and is wrong, and a CDS outside the exons lies in an intron
        overlaps = []
        if len(proper) < len(introns):
            overlaps.append(types.OVERLAPPING_EXONS)
        if cds.has_overlapping_pieces:
            overlaps.append(types.OVERLAPPING_CDS)
        if cds.is_outside_exons:
            overlaps.append(types.CDS_OUTSIDE_EXONS)
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
        # more likely a fragment or an artifact than a gene; Helixer drops predicted genes with less
        # coding sequence by default, too (HelixerPost's min_coding_length of 60)
        if cds.spliced_length < self.controller.config['min_cds_length']:
            found.append(Finding(types.TOO_SHORT_CDS, 'whole'))
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

    @staticmethod
    def _mark_unknown_ends(finding):
        """Sets the biological start or end flags of the features an error finding makes unknown to False."""
        for handler in finding.handlers:
            if finding.side in ('5p', 'both'):
                handler.start_is_biological_start = False
            if finding.side in ('3p', 'both'):
                handler.end_is_biological_end = False

    def resolve_errors(self):
        """For each coding transcript of an exported locus: marks the ends its errors make unknown,
        records which errors it has, and merges every range they mask, together with the ranges
        given up to overlapping loci, into its geenuff_mask features."""
        for i, group in enumerate(self.groups):
            sign = 1 if self.is_plus_strand else -1
            reason = group['super_locus'].excluded_from_export
            # a locus dropped for an overlap was checked like any other, so its errors and unknown
            # ends are recorded, but it gets no mask, its features reaching no export, and no count
            # in the import statistics (see clean_and_insert)
            if reason == types.OVERLAP_DROPPED:
                for transcript in group['transcripts']:
                    for finding in transcript.get('findings', []):
                        self._mark_unknown_ends(finding)
                        transcript['detected_error_types'].add(finding.error_type)
                continue
            # a locus excluded for any other reason is not analyzed at all, and gets no mask: where
            # its features belong is exactly what is unknown about it, so any mask would be
            # guesswork, either covering the wrong strand or covering sequence that is fine
            if reason is not None:
                continue

            # the case of no transcript for a super locus
            if not group['transcripts']:
                logger.debug('{} is a gene without any transcripts; it has no features and will '
                             'never be exported or masked'.format(group['super_locus'].given_name))
            # the selected transcript with the flank on both sides, which every mask reaches to
            buffered = self._buffered_span(i) if i in self._extents else None
            for transcript in group['transcripts']:
                if 'cds' not in transcript:
                    continue
                ranges = []
                for finding in transcript['findings']:
                    self._mark_unknown_ends(finding)
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
        """The sign-normalized range a finding masks, given the buffered span of its locus. A '5p'
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
        """The sequence as a sign-normalized, half open range. On the minus strand the bases
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
        self.config = {'min_intron_length': 20, 'min_cds_length': 60}  # the default config
        self.config.update(config)

    def _mk_session(self, replace_db):
        file_path = self.database_path.removeprefix('sqlite:///')
        if os.path.exists(file_path):
            if not replace_db:
                raise FileExistsError(f'database already existing at {file_path} and --replace-db not set')
            os.remove(file_path)
            logger.info(f'removed existing database at {file_path}')
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
        self.stats.genome_bases = sum(c.length for c in self.latest_fasta_importer.genome.coordinates)
        try:
            self.add_gff(gff_path, clean=clean_gff)
            self.run_analyze()
        except Exception as e:
            self.session.close()
            # the path may be given as an SQLAlchemy URL, and an in-memory database has no file
            file_path = self.database_path.removeprefix('sqlite:///')
            if os.path.isfile(file_path):
                part_path = f'{file_path}.partial'
                shutil.move(file_path, part_path)
                logger.error(f'Aborting due to error, attempt so far saved at {part_path} '
                             f'for debugging purposes: {type(e).__name__}: {e}')
            else:
                logger.error(f'Aborting due to error: {type(e).__name__}: {e}')
            raise
        self.stats.log_summary(species)
        n_total = self.stats.total_coding_transcripts
        if self.stats.cds_outside_exons_transcripts > CDS_OUTSIDE_EXONS_WARN_SHARE * n_total:
            share = self.stats.cds_outside_exons_transcripts / n_total
            logger.warning(f'{self.stats.cds_outside_exons_transcripts} of {n_total} coding transcripts '
                           f'({share:.1%}) have CDS lines outside their exon lines and are masked in full '
                           f'(cds_outside_exons). The file may list exon lines for the UTRs only; check it against '
                           f'its source, and see rebuild_exons.py in docs/scripts.md.')
        if self.stats.overlapping_exons_transcripts > OVERLAPPING_EXONS_WARN_SHARE * n_total:
            share = self.stats.overlapping_exons_transcripts / n_total
            logger.warning(f'{self.stats.overlapping_exons_transcripts} of {n_total} coding transcripts '
                           f'({share:.1%}) have overlapping exon lines and are masked in full (overlapping_exons). '
                           f'The exon lines may be written wrongly, e.g. nested; where the CDS lines are intact, '
                           f'drop the exon lines: GeenuFF then builds the exons from the CDS and UTR lines.')
        n_labeled = self.stats.labeled_genes()
        n_coding = n_total + sum(self.stats.coding_transcripts_dropped.values())
        # this then also shows if the file has 0 coding transcripts (i.e. no CDS lines at all)
        if not n_labeled or n_labeled < LOW_TRAINING_DATA_SHARE * n_coding:
            logger.warning(f'{"NO" if not n_labeled else "LITTLE"} TRAINING DATA for "{species}": {n_labeled} genes '
                           f'giving labels from {n_coding} coding transcripts in the GFF3 file; '
                           f'{self.stats.no_training_data_reason()} (see the summary above). Scripts that can '
                           f'help (GeenuFF/docs/scripts.md): {"; ".join(PREPARING_SCRIPTS)}.')

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
            # give GFFErrorHandling a strand of None to work every neighbor comparison against
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
                sl = group['super_locus']
                reason = sl.excluded_from_export
                coding = [t for t in group['transcripts'] if 'cds' in t]
                step = funnel_step(len(group['transcripts']), sl.n_dropped_transcripts, len(coding), reason)
                if step == EMPTY:
                    self.stats.empty_super_loci += 1
                elif step == LEFT_OUT:
                    self.stats.genes_left_out_on_import[left_out_reason(reason)] += 1
                elif step == NON_CODING:
                    self.stats.non_coding_genes += 1
                    self.stats.non_coding_genes_with_gene_parented_cds += (sl.given_name
                                                                           in self.stats.gene_parented_cds_genes)
                else:
                    self.stats.coding_transcripts_imported += len(coding)
                if reason is not None:
                    self.stats.unexported_super_loci[reason] += 1
                    if step == NOT_EXPORTED:
                        self.stats.coding_genes_not_exported[reason] += 1
                    continue
                for transcript in coding:
                    if not transcript['transcript'].longest:
                        self.stats.unselected_coding_transcripts += 1
                        continue
                    # counted from the error types detected, not from what ended up masked: an
                    # error whose range comes out empty (e.g. a missing UTR of a gene nested in
                    # another, which has no unclaimed sequence to extend the mask into) is still a
                    # real finding about the transcript. Each type is counted once per transcript
                    errors, tf = transcript['detected_error_types'], transcript['transcript_feature']
                    self.stats.add_exported(
                        errors, exported_outcome(errors, tf, transcript['masks'],
                                                 group.get('unresolved_overlap') == 'nested'),
                        (tf.coord.seqid, tf.is_plus_strand), ascending_range(tf.start, tf.end, tf.is_plus_strand),
                        [ascending_range(m.start, m.end, m.is_plus_strand) for m in transcript['masks']])
            # insert importers
            insert_importer_groups(self, plus)
            insert_importer_groups(self, minus)
            insert_importer_groups(self, unstranded)
            if is_final_coord or self.insertion_queues.total_size() > 10000:  # todo: change to 1000, RAM issue fix maybe
                self.insertion_queues.execute_so_far()

        assert self.latest_fasta_importer is not None, 'No recent genome found'
        logger.info('Starting to parse the GFF file')
        self.latest_fasta_importer.match_gff_seqids(gff_file)
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
        self.gffid_to_coords = None

    def match_gff_seqids(self, gff_file):
        """Maps each sequence ID of the GFF's feature lines to its FASTA sequence's Coordinate,
        raising helpers.NonMatchableIDs for IDs without exactly one (see helpers.match_seqids)."""
        coords_by_seqid = {c.seqid: c for c in self.genome.coordinates}
        seqid_pairs = helpers.match_seqids(coords_by_seqid, helpers.get_seqids_from_gff(gff_file))
        self.gffid_to_coords = {gff_id: coords_by_seqid[fa_id] for gff_id, fa_id in seqid_pairs.items()}

    def add_sequences(self, seq_file):
        # todo, parallelize sequence & annotation format, then import directly from ~Slice
        for seqid, seq in self.parse_fasta(seq_file):
            coord = orm.Coordinate(sequence=seq,
                                   length=len(seq),
                                   seqid=seqid,
                                   sha1=helpers.sequence_hash(seq),
                                   genome=self.genome)
            logger.info(f'Added coordinate object for FASTA sequence with seqid {seqid} to the queue')

    def parse_fasta(self, seq_file):
        fp = fastahelper.FastaParser()
        for fasta_header, seq in fp.read_fasta(seq_file):
            seq = seq.upper()  # this may perform poorly
            seqid = helpers.fasta_header_id(fasta_header)
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
        # transcripts of it left out of the database, and the coding ones among them
        self.n_dropped_transcripts = 0
        self.n_dropped_coding_transcripts = 0

    def add_to_queue(self):
        to_add = {'type': self.entry_type, 'given_name': self.given_name, 'id': self.id,
                  'excluded_from_export': self.excluded_from_export,
                  'dropped_transcripts': self.n_dropped_transcripts,
                  'dropped_coding_transcripts': self.n_dropped_coding_transcripts}
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
                 is_outside_exons=False,
                 spliced_length=None,
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
        self.is_outside_exons = is_outside_exons  # only used for the CDS_OUTSIDE_EXONS check
        self.spliced_length = spliced_length  # only used for the TOO_SHORT_CDS check
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
