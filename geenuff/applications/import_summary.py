"""The import summary: the counts an import collects (ImportStatistics) and the funnel they are
logged as, from the GFF lines down to the genes giving labels (format_summary). Shared by the
importer and scripts/summarize_geenuff_db.py, which rebuilds the summary from a database."""
import logging
from collections import Counter, defaultdict

from geenuff.base import helpers, types
from geenuff.base.helpers import merge_ranges, plural

logger = logging.getLogger(__name__)

GENE_TYPES = types.GENE_TYPES

# what each orm.SuperLocus.excluded_from_export value means, for the import summary
UNEXPORTED_REASONS = {
    types.UNPLACEABLE_STRAND: 'their gene, transcript or its exon/CDS lines are not all on one '
                              'definite strand',
    types.UNPLACEABLE_COORDINATES: 'a line of theirs runs backwards, its start past its end',
    types.OUTSIDE_SEQUENCE: 'a line of theirs starts before or ends past their sequence',
    types.OVERLAP_DROPPED: 'dropped so an overlapping partner could be kept whole',
}


# why a GFF line was left out while the lines are grouped into loci by their ID and Parent
# attributes (see importer.OrganizedGFFEntries.load_organized_entries), for the import summary,
# where the GFF types of the lines stand beside it
DROPPED_LINE_REASONS = {
    'no_id': 'without an ID',
    'shared_id': 'sharing an ID with another gene/transcript line',
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
# full (see importer.GFFErrorHandling._find_errors)
NOT_MASKING_IN_FULL = frozenset({types.MISSING_UTR_5P, types.MISSING_UTR_3P, types.WRONG_PHASE_5P,
                                 types.SL_OVERLAP_ERROR})
# the errors that mask a transcript in full before anything else, its UTRs among it, is checked
# (see importer.GFFErrorHandling._find_errors)
UTRS_NOT_CHECKED = frozenset({types.FLOATING_CDS, types.BEYOND_SEQUENCE_EDGE, types.OVERLAPPING_EXONS,
                              types.OVERLAPPING_CDS, types.CDS_OUTSIDE_EXONS})


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
# why all transcripts of a gene were left out on import, keyed as in left_out_reason
LEFT_OUT_REASONS = {
    'no_strand': 'on neither strand',
    types.OUTSIDE_SEQUENCE: 'with a line outside their sequence',
    types.UNPLACEABLE_COORDINATES: 'with a line running backwards',
}


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
    lines); sign is '-' (loose lines/genes/features), '+' (gain genes/features; very, very rare) or ' '
    for a breakdown. Level indents an entry further. An entry with a count of 0 is left out, and so
    is a step without count or entries. The counts of a step are right-aligned to its longest one,
    and a pipe line separates the steps."""
    steps = [(header, count, detail, [e for e in entries if e[1]], extra)
             for header, count, detail, entries, extra in steps]
    steps = [step for step in steps if step[1] is not None or step[3]]
    lines = [f'{title}:', '-' * 100]
    for i, (header, count, detail, entries, extra) in enumerate(steps):
        if i:
            lines.append(f'{" " * 8}|')
        width = max((len(str(n)) for _, n, _, _ in entries), default=1)
        lines.append(header + (f': {count}' if count is not None else ':') + (f' ({detail})' if detail else ''))
        lines += [f'  {"    " * level}{sign} {n:>{width}}  {text}' for sign, n, text, level in entries]
        lines += [f'  {line}' if line else '' for line in extra]
    lines += ['-' * 100] + [f'note: {note}' for note in notes]
    return '\n'.join(lines)


class ImportStatistics(object):
    """Aggregate counts for one add_genome() call in importer.py. A single summary logged at the
    end of the import gives a full picture of what was found (see log_summary)."""

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
        # counted, nothing about it having been checked either (see importer.clean_and_insert)
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
        # importer.OrganizedGFFEntries._place_floating_cds)
        self.floating_cds_transcripts: int = 0
        self.floating_cds_genes: int = 0
        # transcripts with CDS lines but no exon lines, their exons built from their CDS and UTR lines
        self.transcripts_with_exons_built: int = 0
        # coding transcripts with CDS lines outside their exon lines, warned about above a share of
        # importer.CDS_OUTSIDE_EXONS_WARN_SHARE
        self.cds_outside_exons_transcripts: int = 0
        # coding transcripts with exon lines overlapping each other, warned about above a share of
        # importer.OVERLAPPING_EXONS_WARN_SHARE
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
        self.transcripts_outside_sequence_examples: list[str] = []  # the first few, for the warning after import
        self.unstranded_transcripts_dropped: int = 0  # the transcripts of genes on neither strand
        # every pair of super loci sharing genomic range, of whatever kind, as recorded in the
        # super_locus_overlap table; these mask nothing, see importer.GFFErrorHandling._record_overlap_pairs
        self.overlap_pairs_recorded: int = 0
        self.super_loci_in_overlap_pairs: int = 0
        # of the pairs where both genes are exported, how each was settled; the four always sum
        # to that number (see importer.GFFErrorHandling._compute_overlap_masks)
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
        # the exported transcripts by their missing UTRs, 'none', '5p', '3p' or 'both', and those never
        # checked for them (UTRS_NOT_CHECKED); the same for the ones giving labels, with the missing
        # UTRs whose flank mask came out empty and those also with the overhang of an overlapping
        # gene dropped in their favor masked
        self.utrs_missing: Counter[str] = Counter()
        self.utrs_not_checked: int = 0
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
        if errors & UTRS_NOT_CHECKED:
            self.utrs_not_checked += 1
        else:
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
            details.append(f'{plural(sum(below.values()), "line")} below them: {self._per_type(below)}')
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

    def genes_with_transcripts(self) -> int:
        """The imported genes with a transcript under them, coding or not."""
        return self.total_super_loci - self.empty_super_loci - sum(self.genes_left_out_on_import.values())

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

        n_with_transcripts = self.genes_with_transcripts()
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
                                f'gene instead (use add_transcript_lines.py to fix if possible)')
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
                     (utrs['3p'], "only the 3' UTR missing"), (utrs['both'], 'both UTRs missing'),
                     (self.utrs_not_checked, 'not checked, errors making UTR checks impossible')]
        utr_width = len(str(max(n for n, _ in utr_lines)))
        utr_block = ['', f"UTRs of these (5' UTR missing in "
                         f"{self.errors.get(types.MISSING_UTR_5P, 0)}, 3' UTR missing in "
                         f"{self.errors.get(types.MISSING_UTR_3P, 0)}; can overlap with other "
                         f"errors):"] * bool(n_exported)
        utr_block += [f'  {n:>{utr_width}}  {text}' for n, text in utr_lines if n]
        n_wrong_phase = self.errors.get(types.WRONG_PHASE_5P, 0)
        error_lines = self._error_lines('errors masking them in full, each type counted once per transcript '
                                        '(can overlap):', self.errors_masking_in_full)
        steps.append((EXPORTED_CODING_GENES, n_exported, exported_detail,
                      [('-', self.exported_outcomes[outcome], EXPORTED_OUTCOMES[outcome], 0)
                       for outcome in MASKED_IN_FULL_OUTCOMES],
                      [''] * bool(error_lines) + error_lines
                      + utr_block
                      + ['', f'first CDS phase in the file not 0; corrected, masks nothing: {n_wrong_phase}']
                      * bool(n_wrong_phase)))

        # overlapping masks, and the spans of overlapping exported genes, on one strand count once
        merged_masks = {key: merge_ranges(ranges) for key, ranges in self.mask_ranges.items()}
        merged_spans = {key: merge_ranges(ranges) for key, ranges in self.gene_spans.items()}
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
                      + ['', f'{n_masked_bases} of the genome\'s {both_strands} bases on both strands masked '
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
