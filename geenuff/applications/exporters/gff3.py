import logging
from collections import defaultdict
from typing import TextIO

from geenuff.applications.exporter import GeenuffExportController, RangeMaker
from geenuff.base.orm import Transcript, SuperLocus, Feature
from geenuff.base.helpers import geenuff_to_gff_start_end, is_masked_whole, GFF_PHASE_FROM_NORMAL
from geenuff.base import types

logger = logging.getLogger(__name__)

# why a gene never reaches this file, phrased for the export's own log
UNWRITTEN_REASONS = {
    types.UNPLACEABLE_STRAND: 'their features are not all on one definite strand',
    types.UNPLACEABLE_COORDINATES: 'a line of theirs runs backwards, its start past its end',
    types.OUTSIDE_SEQUENCE: 'a line of theirs starts before or ends past their sequence',
    types.OVERLAP_DROPPED: 'they gave way to an overlapping gene that was kept instead',
}


class FilteredGff3ExportController(GeenuffExportController):
    """Writes a plain GFF3 file of one transcript per gene, the longest coding one, for comparing
    a Helixer prediction against the reference it was trained on.

    By default, it writes exactly the transcripts the h5 export gives labels from: those of genes
    reaching that export that are not masked whole. One with only sequence beside it masked, for a
    missing UTR or for an overlapping gene dropped in its favour, is labelled in full and written.

    With include_erroneous it writes every gene that can be written at all, whatever is wrong with
    it, which is the set to compare against when the question is what Helixer predicted per gene
    rather than how it did on sound ones. Only genes that cannot be written are left out then, and
    a gene dropped for overlapping another is not one of them: nothing is wrong with it beyond
    sharing sequence, which a GFF3 holds without trouble (see types.unrepresentable_reasons).
    The mRNA line of an erroneous transcript names its error types in a geenuff_errors attribute,
    and that of a gene kept out of the h5 export the reason in a geenuff_excluded attribute."""

    def write_filtered_gff3(self, file_out: str | None, include_erroneous: bool = False) -> None:
        handle_out = self._as_file_handle(file_out)
        handle_out.write('##gff-version 3\n')

        n_written = 0
        n_masked_whole = 0
        n_overlap_dropped = 0
        skipped = defaultdict(int)  # keyed by SuperLocus.excluded_from_export
        # the super locus is selected alongside so its reason and name come from the same query
        rows = (self.session.query(Transcript, SuperLocus)
                .join(SuperLocus, Transcript.super_locus_id == SuperLocus.id)
                .filter(Transcript.longest.is_(True))
                .order_by(SuperLocus.id))
        for transcript, super_locus in rows:
            reason = super_locus.excluded_from_export
            if reason is not None and (not include_erroneous
                                       or reason in types.unrepresentable_reasons):
                skipped[reason] += 1
                continue
            n_overlap_dropped += reason == types.OVERLAP_DROPPED
            features = [f for piece in transcript.transcript_pieces for f in piece.features]
            tx_feature = next(f for f in features if f.type.value == types.GEENUFF_TRANSCRIPT)
            masks = [f for f in features if f.type.value == types.GEENUFF_MASK]
            if is_masked_whole(tx_feature, masks):
                n_masked_whole += 1
                if not include_erroneous:
                    continue
            self._write_transcript(handle_out, transcript)
            n_written += 1

        if file_out is not None:
            handle_out.close()
        self._log_summary(n_written, n_masked_whole, n_overlap_dropped, skipped, include_erroneous)

    @staticmethod
    def _log_summary(n_written: int, n_masked_whole: int, n_overlap_dropped: int,
                     skipped: dict[str, int], include_erroneous: bool) -> None:
        if include_erroneous:
            logger.info(f'Wrote {n_written} transcripts to GFF3, one per gene, {n_masked_whole} of '
                        f'them masked whole and {n_overlap_dropped} dropped from the h5 export for '
                        f'an overlap')
        else:
            logger.info(f'Wrote {n_written} transcripts to GFF3, the ones the h5 export gives labels '
                        f'from ({n_masked_whole} masked whole left out)')
        for reason, count in sorted(skipped.items()):
            logger.info(f'  {count} genes left out, {UNWRITTEN_REASONS[reason]}')

    def _write_transcript(self, handle_out: TextIO, transcript: Transcript) -> None:
        range_maker = RangeMaker(transcript)
        features_by_type: dict[str, Feature] = {}
        for feature, _ in range_maker.feature_piece_pairs():
            features_by_type[feature.type.value] = feature
        tx_feature = features_by_type[types.GEENUFF_TRANSCRIPT]

        seqid = tx_feature.coordinate.seqid
        is_plus_strand = tx_feature.is_plus_strand
        strand = '+' if is_plus_strand else '-'
        source = tx_feature.source or 'GeenuFF'
        score = '.' if tx_feature.score is None else tx_feature.score

        gene_id = transcript.super_locus.given_name or f'gene{transcript.super_locus_id}'
        mrna_id = transcript.given_name or f'mRNA{transcript.id}'

        # what is wrong with the transcript, absent where nothing is
        mrna_attributes = f'ID={mrna_id};Parent={gene_id}'
        errors = sorted(e.type.value for e in transcript.errors)
        if errors:
            mrna_attributes += f';geenuff_errors={",".join(errors)}'
        if transcript.super_locus.excluded_from_export is not None:
            mrna_attributes += f';geenuff_excluded={transcript.super_locus.excluded_from_export}'

        tx_start, tx_end = geenuff_to_gff_start_end(tx_feature.start, tx_feature.end, is_plus_strand)
        handle_out.write(self._gff_line(seqid, source, 'gene', tx_start, tx_end, score, strand,
                                        '.', f'ID={gene_id}'))
        handle_out.write(self._gff_line(seqid, source, 'mRNA', tx_start, tx_end, score, strand,
                                        '.', mrna_attributes))

        exons = [group.ranges[0] for group in range_maker.exonic_ranges()]
        for i, exon in enumerate(exons, start=1):
            e_start, e_end = geenuff_to_gff_start_end(exon.start, exon.end, is_plus_strand)
            handle_out.write(self._gff_line(seqid, source, 'exon', e_start, e_end, score, strand,
                                            '.', f'ID={mrna_id}.exon{i};Parent={mrna_id}'))

        cds_pieces = [group.ranges[0] for group in range_maker.cds_exonic_ranges()]
        cds_feature = features_by_type[types.GEENUFF_CDS]
        cumulative_len = 0
        starting_phase = GFF_PHASE_FROM_NORMAL[cds_feature.phase]
        for cds in cds_pieces:
            c_start, c_end = geenuff_to_gff_start_end(cds.start, cds.end, is_plus_strand)
            phase = GFF_PHASE_FROM_NORMAL[(starting_phase + cumulative_len) % 3]
            handle_out.write(self._gff_line(seqid, source, 'CDS', c_start, c_end, score, strand,
                                            phase, f'ID={mrna_id}.cds;Parent={mrna_id}'))
            cumulative_len += abs(cds.end - cds.start)

    @staticmethod
    def _gff_line(seqid: str, source: str, feature_type: str, start: int, end: int,
                  score: float | str, strand: str, phase: int | str, attributes: str) -> str:
        cols = [seqid, source, feature_type, start, end, score, strand, phase, attributes]
        return '\t'.join(str(c) for c in cols) + '\n'
