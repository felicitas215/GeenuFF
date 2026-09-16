import logging
from typing import TextIO

from geenuff.applications.exporter import GeenuffExportController, RangeMaker
from geenuff.base.orm import Transcript, SuperLocus, Feature
from geenuff.base.helpers import geenuff_to_gff_start_end, GFF_PHASE_FROM_NORMAL
from geenuff.base import types

logger = logging.getLogger(__name__)


class FilteredGff3ExportController(GeenuffExportController):
    """Writes a plain GFF3 file containing exactly the transcripts Helixer's h5 export
    would use for training: the longest (transcript.longest) transcript per super locus,
    and only where that transcript carries no error feature at all. This is a stricter
    filter than the h5 export itself, which still includes and merely masks erroneous
    transcripts. Useful for comparing the Helixer predictions to the filtered reference,
    i.e. structurally sound transcripts."""

    def write_filtered_gff3(self, file_out: str | None) -> None:
        handle_out = self._as_file_handle(file_out)
        handle_out.write('##gff-version 3\n')

        n_written = 0
        n_skipped_erroneous = 0
        transcripts = (self.session.query(Transcript)
                       .join(SuperLocus, Transcript.super_locus_id == SuperLocus.id)
                       .filter(Transcript.longest.is_(True))
                       .order_by(SuperLocus.id))
        for transcript in transcripts:
            if not self._is_error_free(transcript):
                n_skipped_erroneous += 1
                continue
            self._write_transcript(handle_out, transcript)
            n_written += 1

        if file_out is not None:
            handle_out.close()
        logger.info(f'Wrote {n_written} error-free transcripts to GFF3 '
                    f'({n_skipped_erroneous} longest-but-erroneous transcripts skipped)')

    @staticmethod
    def _is_error_free(transcript: Transcript) -> bool:
        for piece in transcript.transcript_pieces:
            for feature in piece.features:
                if feature.type.value in types.geenuff_error_type_values:
                    return False
        return True

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

        tx_start, tx_end = geenuff_to_gff_start_end(tx_feature.start, tx_feature.end, is_plus_strand)
        handle_out.write(self._gff_line(seqid, source, 'gene', tx_start, tx_end, score, strand,
                                        '.', f'ID={gene_id}'))
        handle_out.write(self._gff_line(seqid, source, 'mRNA', tx_start, tx_end, score, strand,
                                        '.', f'ID={mrna_id};Parent={gene_id}'))

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
