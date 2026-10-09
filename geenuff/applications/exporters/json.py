import json
from abc import abstractmethod

from geenuff.applications.exporter import GeenuffExportController
from geenuff.applications.importer import error_severity
from geenuff.base import types
from geenuff.base.handlers import SuperLocusHandlerBase, TranscriptHandlerBase, CoordinateHandlerBase, \
    FeatureHandlerBase
from geenuff.base.helpers import is_masked_whole
from geenuff.base.orm import Coordinate, Transcript


class ToJsonable(object):
    """The query window is an ascending, half-open range [start, end) on one strand of one
    sequence, whichever the strand."""

    def __init__(self):
        self.jsonable_keys = ["id", "given_name"]

    @abstractmethod
    def is_fully_contained(self, coordinate, start, end, is_plus_strand):
        pass

    @abstractmethod
    def overlaps(self, coordinate, start, end, is_plus_strand):
        pass

    def pre_to_jsonable(self, data):
        out = {}
        for key in self.jsonable_keys:
            out[key] = data.__getattribute__(key)
        return out

    @abstractmethod
    def to_jsonable(self, data, coordinate, start, end, is_plus_strand):
        pass


# then multi inheritance of above and Handlers
class FeatureJsonable(FeatureHandlerBase, ToJsonable):
    def __init__(self, data=None):
        FeatureHandlerBase.__init__(self, data)
        ToJsonable.__init__(self)
        self.jsonable_keys += ["start", "start_is_biological_start", "end", "end_is_biological_end",
                               "is_plus_strand", "score", "source", "phase"]

    def to_jsonable(self, data, coordinate, start, end, is_plus_strand, transcript=None):
        assert transcript is not None, "needed to get the protein id"
        out = self.pre_to_jsonable(self.data)
        out['type'] = data.type.value
        out['is_fully_contained'] = self.is_fully_contained(coordinate, start, end, is_plus_strand)
        out['overlaps'] = self.overlaps(coordinate, start, end, is_plus_strand)
        out['protein_id'] = self.protein_ids(transcript)
        return out

    def protein_ids(self,  transcript):
        assert isinstance(transcript, Transcript)  # just for pycharm hints for now...
        p_ids = list(set([x.given_name for x in transcript.proteins]).intersection(
            set([x.given_name for x in self.data.proteins])
        ))
        l_proteins = len(p_ids)
        if l_proteins == 0:
            return None
        elif l_proteins == 1:
            return p_ids[0]
        else:
            # todo, actually, I don't think this should be able to happen, but double check
            raise NotImplementedError("what to do with {} proteins?".format(len(p_ids)))

    def _ascending_span(self):
        """The bases the feature covers, as an ascending half-open range; on the minus strand a
        GeenuFF start lies above its end (see docs/spec_vs_gff.md)."""
        if self.data.is_plus_strand:
            return self.data.start, self.data.end
        return self.data.end + 1, self.data.start + 1

    def _on_queried_strand(self, coordinate, is_plus_strand):
        return coordinate.id == self.data.coordinate_id and is_plus_strand == self.data.is_plus_strand

    def is_fully_contained(self, coordinate, start, end, is_plus_strand):
        lo, hi = self._ascending_span()
        return self._on_queried_strand(coordinate, is_plus_strand) and start <= lo and hi <= end

    def overlaps(self, coordinate, start, end, is_plus_strand):
        lo, hi = self._ascending_span()
        return self._on_queried_strand(coordinate, is_plus_strand) and lo < end and start < hi


class TranscriptJsonable(TranscriptHandlerBase, ToJsonable):
    def __init__(self, data=None):
        FeatureHandlerBase.__init__(self, data)
        ToJsonable.__init__(self)
        self.feature_handlers = self._mk_feature_handlers()

    def _mk_feature_handlers(self):
        return [FeatureJsonable(x) for x in self.sorted_features()]

    def sorted_features(self):
        out = []
        for piece in self.sorted_pieces:
            # todo, access the sort_cmp_keys for consistency with code base?
            if piece.features[0].is_plus_strand:
                features = sorted(piece.features, key=lambda f: (f.start, f.end))
            else:
                features = sorted(piece.features, key=lambda f: (f.start * -1, f.end * 1))
            out += features
        return out

    def _gene_model_handlers(self):
        """The feature handlers of the gene model itself, a geenuff_mask lying beside it as well"""
        return [fh for fh in self.feature_handlers if fh.data.type.value != types.GEENUFF_MASK]

    # todo, use pre-calculated jsons not fresh method call
    def overlaps(self, coordinate, start, end, is_plus_strand):
        return any(fh.overlaps(coordinate, start, end, is_plus_strand) for fh in self._gene_model_handlers())

    def is_fully_contained(self, coordinate, start, end, is_plus_strand):
        return all(fh.is_fully_contained(coordinate, start, end, is_plus_strand)
                   for fh in self._gene_model_handlers())

    def to_jsonable(self, data, coordinate, start, end, is_plus_strand):
        features = [fh.data for fh in self.feature_handlers]
        tx_feature = next((f for f in features if f.type.value == types.GEENUFF_TRANSCRIPT), None)
        masks = [f for f in features if f.type.value == types.GEENUFF_MASK]
        errors = {e.type.value for e in self.data.errors}

        out = self.pre_to_jsonable(self.data)
        out["type"] = self.data.type.value
        out["is_fully_contained"] = self.is_fully_contained(coordinate, start, end, is_plus_strand)
        out["overlaps"] = self.overlaps(coordinate, start, end, is_plus_strand)
        # the one transcript per gene an export selects, and whether its gene reaches the export
        out["selected_for_export"] = bool(self.data.longest)
        out["exported"] = out["selected_for_export"] and self.data.super_locus.excluded_from_export is None
        out["errors"] = sorted(errors)
        out["error_severity"] = error_severity(errors)
        # what the errors mask, merged into its geenuff_mask features, also listed among the features
        out["masks"] = [[m.start, m.end] for m in masks]
        out["masked_in_full"] = tx_feature is not None and is_masked_whole(tx_feature, masks)
        out["features"] = [fh.to_jsonable(fh.data, coordinate, start, end, is_plus_strand, self.data)
                           for fh in self.feature_handlers]
        return out


class SuperLocusJsonable(SuperLocusHandlerBase, ToJsonable):
    def __init__(self, data=None, longest=False):
        FeatureHandlerBase.__init__(self, data)
        ToJsonable.__init__(self)
        self.longest = longest
        self.transcript_handlers = self._mk_transcript_handlers()

    def _mk_transcript_handlers(self):
        """Every isoform, or with longest only the transcript selected for export, as the h5
        export uses it (see OrganizedGeenuffImporterGroup._set_longest_transript)."""
        return [TranscriptJsonable(x) for x in self.data.transcripts if not self.longest or x.longest]

    def is_fully_contained(self, coordinate, start, end, is_plus_strand):
        return all(th.is_fully_contained(coordinate, start, end, is_plus_strand) for th in self.transcript_handlers)

    def overlaps(self, coordinate, start, end, is_plus_strand):
        return any(th.overlaps(coordinate, start, end, is_plus_strand) for th in self.transcript_handlers)

    def to_jsonable(self, data, coordinate, start, end, is_plus_strand):
        out = self.pre_to_jsonable(self.data)
        out['type'] = self.data.type.value
        out['is_fully_contained'] = self.is_fully_contained(coordinate, start, end, is_plus_strand)
        out['overlaps'] = self.overlaps(coordinate, start, end, is_plus_strand)
        # why the gene is left out of exports, None where it is not
        out['excluded_from_export'] = self.data.excluded_from_export
        out['exported'] = (self.data.excluded_from_export is None
                           and any(t.longest for t in self.data.transcripts))
        out['transcripts'] = [th.to_jsonable(th, coordinate, start, end, is_plus_strand)
                              for th in self.transcript_handlers]
        return out


class CoordinateJsonable(CoordinateHandlerBase):
    def to_jsonable(self, start, end):
        return {'id': self.data.id,
                'seqid': self.data.seqid,
                'sequence': self.data.sequence[start:end],
                'start': start,
                'end': end}


class JsonExportController(GeenuffExportController):
    # todo, filter coordinate by start, end
    #  from orm obj (or join res) to json

    def coordinate_range_to_jsonable(self, species, seqid, start, end, is_plus_strand):
        """The genes with a transcript on one strand of [start, end) of a sequence, as written by
        query_and_write; end None reaches to the end of the sequence."""
        out = []
        for coordinate in self.session.query(Coordinate).filter(Coordinate.seqid == seqid).all():
            if coordinate.genome.species == species:
                if end is None:  # if end is not specified, take whole sequence
                    end = coordinate.length
                ch = CoordinateJsonable(coordinate)
                res = {'coordinate_piece': ch.to_jsonable(start, end),
                       'super_loci': []}
                for sl, sl_coordinate_seqid in self.genome_query(return_super_loci=True):
                    if sl_coordinate_seqid == seqid:
                        slh = SuperLocusJsonable(sl, longest=self.longest)
                        if slh.overlaps(coordinate, start, end, is_plus_strand):
                            res['super_loci'].append(slh.to_jsonable(slh.data, coordinate, start, end,
                                                                     is_plus_strand))
                out.append(res)
        return out

    def coordinate_range_to_json(self, species, seqid, start, end, is_plus_strand):
        return json.dumps(self.coordinate_range_to_jsonable(species, seqid, start, end, is_plus_strand))

    def query_and_write(self, species, seqid, start, end, is_plus_strand, file_out, pretty=False):
        jsonable = self.coordinate_range_to_jsonable(species, seqid, start, end, is_plus_strand)
        if pretty:
            dumps = json.dumps(jsonable, indent=2)
        else:
            dumps = json.dumps(jsonable)

        handle_out = self._as_file_handle(file_out)
        handle_out.write(dumps)
        handle_out.close()
