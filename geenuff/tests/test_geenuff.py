import argparse
import os
import logging
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pytest
from sqlalchemy import create_engine
from sqlalchemy import func
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import IntegrityError
import sqlalchemy

from geenuff.base import orm
from geenuff.base import types
from geenuff.base import helpers
from geenuff.base.orm import (Genome, Feature, Coordinate, Transcript, TranscriptPiece, SuperLocus,
                              Protein, TranscriptError)
from geenuff.base.handlers import SuperLocusHandlerBase, TranscriptHandlerBase
from geenuff.applications.importer import ImportController, MissingGeneLines, OrganizedGFFEntries, error_severity
from geenuff.applications.import_summary import GENES_GIVING_LABELS, SUMMARY_STEPS
from geenuff.applications.exporter import GeenuffExportController
from geenuff.applications.exporters.gff3 import FilteredGff3ExportController
from geenuff.parsing import gffhelper


@pytest.fixture(scope="session", autouse=True)
def prepare(request):
    if not os.getcwd().endswith('GeenuFF/geenuff'):
        pytest.exit('Tests need to be run from GeenuFF/geenuff directory')


### Helper functions ###
def mk_memory_session(db_path='sqlite:///:memory:'):
    engine = create_engine(db_path, echo=False)
    orm.Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    return Session()


def setup_data_handler(handler_type, data_type, **kwargs):
    data = data_type(**kwargs)
    handler = handler_type()
    handler.add_data(data)
    return data, handler


def run_script(name: str, *args: str) -> None:
    subprocess.run([sys.executable, str(Path(__file__).parents[2] / 'scripts' / name), *args], check=True)


def feature_lines(path: Path) -> list[list[str]]:
    return [line.split('\t') for line in path.read_text().splitlines() if line and not line.startswith('#')]


# Do not implement as __eq__ to not change __hash__ behavior
def matching_super_loci(s1, s2):
    if not isinstance(s1, s2.__class__):
        return False
    if s1.id is not None and s2.id is not None:
        if s1.id != s2.id:
            return False
    if s1.given_name != s2.given_name or s1.type != s2.type:
        return False
    return True


def matching_features(f1, f2):
    if not isinstance(f1, f2.__class__):
        return False
    if f1.id is not None and f2.id is not None:
        if f1.id != f2.id:
            return False
    if (f1.type != f2.type or f1.given_name != f2.given_name or f1.start != f2.start
            or f1.end != f2.end or f1.start_is_biological_start != f2.start_is_biological_start
            or f1.end_is_biological_end != f2.end_is_biological_end
            or f1.is_plus_strand != f2.is_plus_strand or f1.phase != f2.phase
            or f1.coordinate.id != f2.coordinate.id):
        return False
    return True


def matching_proteins(t1, t2):
    if not isinstance(t1, t2.__class__):
        return False
    if t1.id is not None and t2.id is not None:
        if t1.id != t2.id:
            return False
    if t1.given_name != t2.given_name or t1.super_locus.id != t2.super_locus.id:
        return False
    return True


def matching_transcripts(t1, t2):
    if not isinstance(t1, t2.__class__):
        return False
    if t1.id is not None and t2.id is not None:
        if t1.id != t2.id:
            return False
    if (t1.given_name != t2.given_name or t1.super_locus.id != t2.super_locus.id
            or t1.type != t2.type):
        return False
    return True


def orm_object_in_list(obj, obj_list):
    for o in obj_list[:]:  # make a copy at each iteration so we avoid weird errors
        if ((isinstance(o, Feature) and matching_features(o, obj))
                or (isinstance(o, Protein) and matching_proteins(o, obj))
                or (isinstance(o, Transcript) and matching_transcripts(o, obj))):
            obj_list.remove(o)
            return True
    return False


### The actual tests ###
def test_annogenome2coordinate_relation():
    """Check if everything is consistent when we add a Genome and a Coordinate
    to the db. Also check for correct deletion behavior.
    """
    sess = mk_memory_session()
    g = Genome(species='Athaliana', version='1.2', acquired_from='Phytozome12')
    coord = Coordinate(seqid='abc', genome=g)
    assert g is coord.genome
    # actually put everything in db
    sess.add(coord)
    sess.commit()
    # check primary keys were assigned
    assert g.id == 1
    assert coord.id == 1
    # check we can access coordinates from g
    coord_q = g.coordinates[0]
    assert coord is coord_q
    assert g is coord.genome
    assert g.id == coord_q.genome_id
    # check we get logical behavior on deletion
    sess.delete(coord)
    sess.commit()
    assert len(g.coordinates) == 0
    print(coord.genome)
    sess.delete(g)
    sess.commit()
    with pytest.raises(sqlalchemy.exc.InvalidRequestError):
        sess.add(coord)


def test_coordinate_constraints():
    """Check the coordinate constraints"""
    sess = mk_memory_session()
    g = Genome()

    # should be ok
    coors = Coordinate(length=30, seqid='abc', genome=g)
    coors2 = Coordinate(length=400, seqid='abcd', genome=g)
    sess.add_all([coors, coors2])
    sess.commit()

    # start/end constraints
    coors_bad1 = Coordinate(length=-12, seqid='abc', genome=g)
    coors_bad2 = Coordinate(genome=g)
    with pytest.raises(IntegrityError):
        sess.add(coors_bad1)  # start below 1
        sess.commit()
    sess.rollback()
    with pytest.raises(IntegrityError):
        sess.add(coors_bad2)  # no seqid
        sess.commit()


def test_coordinate_insert():
    """Test what happens when we insert two coordinates"""
    sess = mk_memory_session()
    g = Genome()
    coords = Coordinate(length=10, seqid='abc', genome=g)
    coords2 = Coordinate(length=11, seqid='def', genome=g)
    sl = SuperLocus()
    f0 = Feature(coordinate=coords)
    f1 = Feature(coordinate=coords2)
    # should be ok
    sess.add_all([g, sl, coords, coords2, f0, f1])
    assert f0.coordinate.length == 10
    assert f1.coordinate.length == 11


def test_many2many_with_features():
    """Test the many2many tables association_transcript_piece_to_feature and
    association_protein_to_feature
    """
    sl = SuperLocus()
    # one transcript, multiple proteins
    piece0 = TranscriptPiece()
    slated0 = Protein(super_locus=sl)
    slated1 = Protein(super_locus=sl)
    # features representing alternative start codon for proteins on one transcript
    feat0_tss = Feature(transcript_pieces=[piece0])
    feat2_stop = Feature(proteins=[slated0, slated1])
    feat3_start = Feature(proteins=[slated0])
    # test multi features per protein worked
    assert len(slated0.features) == 2
    # test mutli protein per feature worked
    assert len(feat2_stop.proteins) == 2
    assert len(feat3_start.proteins) == 1
    assert len(feat0_tss.proteins) == 0


def test_feature_has_its_things():
    """Test if properties of Feature table are correct and constraints are enforced"""
    sess = mk_memory_session()
    # test feature with nothing much set
    g = Genome()
    c = Coordinate(length=30, seqid='abc', genome=g)
    f = Feature(coordinate=c,
                start=1,
                end=30,
                start_is_biological_start=True,
                end_is_biological_end=True,
                is_plus_strand=True)
    sess.add_all([f, c])
    sess.commit()

    assert (f.source, f.score) == (None, None)
    # test feature with
    f1 = Feature(coordinate=c,
                 start=3,
                 end=-1,
                 start_is_biological_start=True,
                 end_is_biological_end=True,
                 is_plus_strand=False)
    assert not f1.is_plus_strand
    assert f1.start == 3
    assert f1.end == -1
    sess.add(f1)
    sess.commit()

    # test too low of start / end coordinates raise an error
    f_should_fail = Feature(coordinate=c,
                            start=-5,
                            end=10,
                            start_is_biological_start=True,
                            end_is_biological_end=True,
                            is_plus_strand=False)
    sess.add(f_should_fail)
    with pytest.raises(sqlalchemy.exc.IntegrityError):
        sess.commit()
    sess.rollback()
    f_should_fail = Feature(coordinate=c,
                            start=5,
                            end=-2,
                            start_is_biological_start=True,
                            end_is_biological_end=True,
                            is_plus_strand=False)
    sess.add(f_should_fail)
    with pytest.raises(sqlalchemy.exc.IntegrityError):
        sess.commit()
    sess.rollback()

    # test wrong class as parameter
    with pytest.raises(KeyError):
        f2 = Feature(coordinate=f)

    f2 = Feature(coordinate=c,
                 start=1,
                 end=3,
                 start_is_biological_start=True,
                 end_is_biological_end=True,
                 is_plus_strand=-1)  # note that 0, and 1 are accepted
    sess.add(f2)
    with pytest.raises(sqlalchemy.exc.StatementError):
        sess.commit()
    sess.rollback()

    # check 'phase is NULL or (not start_is_biological_start or phase = 0)'
    f = Feature(coordinate=c,
                start=1,
                end=3,
                start_is_biological_start=True,
                end_is_biological_end=True,
                is_plus_strand=True,
                phase=1)
    sess.add(f)
    with pytest.raises(sqlalchemy.exc.IntegrityError):
        sess.commit()
    sess.rollback()


def test_partially_remove_coordinate():
    """Add two coordinates to an annotated genome and test if everything ends
    up valid.
    """
    sess = mk_memory_session()
    g = Genome()
    place_holder = Genome()
    coord0 = Coordinate(length=30, seqid='abc', genome=g)
    coord1 = Coordinate(length=330, seqid='def', genome=g)
    sess.add_all([g, coord0, coord1])
    sess.commit()
    assert len(g.coordinates) == 2
    g.coordinates.remove(coord0)
    coord0.genome = place_holder  # else we'll fail the not NULL constraint
    sess.commit()
    # removed from g
    assert len(g.coordinates) == 1
    # but still in table
    assert len(sess.query(Coordinate).all()) == 2


def test_import_intron_at_seq_end():
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/intron_at_end.fa', 'testdata/intron_at_end.gff3', clean_gff=True)
    # one gene model on the + strand
    features = controller.session.query(Feature).filter(Feature.coordinate_id == 1).all()
    # here we have a partial gene model, that runs of the end (or + strand start) fo the
    # sequence in the middle of an intron
    # because the transcript continues to the end, this is in fact an unambiguous intron,
    # albeit with end is biological end being false
    # should produce
    # geenuff_transcript 559 - -1 (not/not biological start/end)
    # geenuff_cds 559 - -1 (y/not biological start/end)
    # geenuff_intron 49 - -1 (y/not biolical start/end)
    # errors missing_utr_5p and missing_utr_3p, the CDS spanning the whole transcript, and
    # truncated_intron and missing_stop_codon, the gene being partial and its CDS ending without a
    # stop codon: the whole gene and the flank on both sides is masked, the 3' one having no room
    # left before the sequence start, taking in the 5' flank of missing_utr_5p, so one
    # geenuff_mask 879 - -1
    assert len(features) == 4
    transcript = [f for f in features if f.type.value == types.GEENUFF_TRANSCRIPT][0]
    cds = [f for f in features if f.type.value == types.GEENUFF_CDS][0]
    intron = [f for f in features if f.type.value == types.GEENUFF_INTRON][0]
    mask = [f for f in features if f.type.value == types.GEENUFF_MASK][0]

    # coordinates
    assert (transcript.start, transcript.end) == (559, -1)
    assert (cds.start, cds.end) == (559, -1)
    assert (intron.start, intron.end) == (49, -1)
    # int(sqrt(1040)) * 10 = 320 of the 1040bp (560-1599) between the gene and the end of the
    # sequence, no gene lying that way to bound it
    assert (mask.start, mask.end) == (879, -1)
    t = controller.session.query(Transcript).one()
    assert error_types(t) == {types.MISSING_UTR_5P, types.MISSING_UTR_3P, types.TRUNCATED_INTRON,
                              types.MISSING_STOP_CODON}

    # biological start / ends marked correctly
    assert (transcript.start_is_biological_start, transcript.end_is_biological_end) == (False, False)
    assert (cds.start_is_biological_start, cds.end_is_biological_end) == (False, False)
    # note that the cds.start is presumably the start codon here, but it is ultimately ambiguous, so
    # conservatively we mark start_is_biological_start as False
    assert (intron.start_is_biological_start, intron.end_is_biological_end) == (True, False)


@pytest.fixture(scope='module')
def hanging_intron_session():
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/hanging_intron.fa', 'testdata/hanging_intron.gff3', clean_gff=True)
    return controller.session


def error_types(transcript):
    """The error types recorded for a transcript."""
    return {e.type.value for e in transcript.errors}


def transcript_features(session, given_name):
    """{feature type: sorted (start, end, start_is_biological_start, end_is_biological_end)} of
    one transcript, with its error types under 'errors'."""
    transcript = session.query(Transcript).filter_by(given_name=given_name).one()
    features = defaultdict(list)
    for f in transcript.transcript_pieces[0].features:
        features[f.type.value].append((f.start, f.end, f.start_is_biological_start,
                                       f.end_is_biological_end))
    return {feature_type: sorted(found) for feature_type, found in features.items()} | {
        'errors': error_types(transcript)}


def test_a_transcript_starting_in_an_intron_is_masked_whole(hanging_intron_session):
    """The mRNA starts 11bp before its first exon, so the transcript starts in an intron
    (233-244) and where it really starts is unknown: the gene (233-566) is masked whole with
    116 of the 233bp to the sequence start and 200 of the 434bp to its end. The CDS starts inside
    the first exon and is sound. See testdata/hanging_intron.gff3."""
    assert transcript_features(hanging_intron_session, 'rnaHang5') == {
        types.GEENUFF_TRANSCRIPT: [(233, 566, False, True)],
        types.GEENUFF_CDS: [(259, 501, True, True)],
        types.GEENUFF_INTRON: [(233, 244, False, True), (279, 350, True, True)],
        types.GEENUFF_MASK: [(117, 766, True, True)],
        'errors': {types.TRUNCATED_INTRON},
    }


def test_a_cds_flush_with_a_first_exon_after_a_hanging_intron_runs_to_the_transcript_start(
        hanging_intron_session):
    """As above, but the CDS starts at the first exon's start. Where the CDS really starts is then
    unknown too, so it is run on through the hanging intron to the transcript start, which leaves
    no 5' UTR: missing_utr_5p masks the 5' flank. See testdata/hanging_intron.gff3."""
    assert transcript_features(hanging_intron_session, 'rnaHang5Flush') == {
        types.GEENUFF_TRANSCRIPT: [(233, 566, False, True)],
        types.GEENUFF_CDS: [(233, 501, False, True)],
        types.GEENUFF_INTRON: [(233, 244, False, True), (279, 350, True, True)],
        types.GEENUFF_MASK: [(117, 766, True, True)],
        'errors': {types.MISSING_UTR_5P, types.TRUNCATED_INTRON},
    }


def test_a_cds_outside_its_exons_is_masked_whole(caplog):
    """Exon lines covering the UTRs only leave the CDS in an intron: cds_outside_exons masks the
    gene (100-700) whole with 50 of the 100bp to the sequence start and 150 of the 300bp to its end,
    while one exon spanning UTRs and CDS, or exon lines touching each other around the CDS, cover it.
    The import warns, as more than 5% of the coding transcripts have the error. See
    testdata/cds_outside_exons.gff3."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/cds_outside_exons.fa', 'testdata/cds_outside_exons.gff3', clean_gff=True)
    session = controller.session
    assert errors_of(session, 'rnaUtrOnly') == {types.CDS_OUTSIDE_EXONS}
    assert masks_of(session, 'rnaUtrOnly') == [(50, 850)]
    for transcript in ('rnaOneExon', 'rnaControl'):
        assert errors_of(session, transcript) == set()
        assert masks_of(session, transcript) == []
    assert any('1 of 3 coding transcripts' in r.message and r.levelname == 'WARNING' for r in caplog.records)


def test_rebuild_exons_covers_the_cds_of_utr_only_exon_lines(tmp_path):
    """rebuild_exons.py replaces the UTR-only exon lines of rnaUtrOnly with one exon spanning its
    UTR and CDS lines, as rnaOneExon has, which then import without an error, and leaves the other
    transcripts as they are."""
    out = tmp_path / 'rebuilt.gff3'
    run_script('rebuild_exons.py', '--gff', 'testdata/cds_outside_exons.gff3', '-o', str(out))
    exons = [columns[3:5] + [columns[8]] for columns in feature_lines(out) if columns[2] == 'exon']
    assert exons == [['101', '700', 'Parent=rnaUtrOnly'],
                     ['101', '700', 'ID=exon-rnaOneExon;Parent=rnaOneExon'],
                     ['101', '200', 'ID=exon-rnaControl-1;Parent=rnaControl'],
                     ['201', '599', 'ID=exon-rnaControl-2;Parent=rnaControl'],
                     ['600', '700', 'ID=exon-rnaControl-3;Parent=rnaControl']]

    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/cds_outside_exons.fa', str(out), clean_gff=True)
    for transcript in ('rnaUtrOnly', 'rnaOneExon', 'rnaControl'):
        assert errors_of(controller.session, transcript) == set()
        assert masks_of(controller.session, transcript) == []


def test_a_cds_flush_with_a_last_exon_before_a_hanging_intron_runs_to_the_transcript_end(
        hanging_intron_session):
    """The mRNA ends 46bp after its last exon (intron 654-700), and the CDS ends at that exon's
    end, so the CDS is run on to the transcript end and missing_utr_3p masks the 3' flank. The
    gene (300-700) is masked whole with 150 of the 300bp on either side.
    See testdata/hanging_intron.gff3."""
    assert transcript_features(hanging_intron_session, 'rnaHang3Flush') == {
        types.GEENUFF_TRANSCRIPT: [(300, 700, True, False)],
        types.GEENUFF_CDS: [(320, 700, True, False)],
        types.GEENUFF_INTRON: [(400, 470, True, True), (654, 700, True, False)],
        types.GEENUFF_MASK: [(150, 850, True, True)],
        'errors': {types.MISSING_UTR_3P, types.TRUNCATED_INTRON},
    }


# section: api
def test_transcript_piece_unique_constraints():
    """Add transcript pieces in valid and invalid configurations and test for
    valid outcomes.
    """
    sess = mk_memory_session()
    sl = SuperLocus()
    transcript0 = Transcript(super_locus=sl)
    transcript1 = Transcript(super_locus=sl)

    # test if same position for different transcript_id goes through
    piece_tr0_pos0 = TranscriptPiece(transcript=transcript0, position=0)
    piece_tr1_pos0 = TranscriptPiece(transcript=transcript1, position=0)
    sess.add_all([transcript0, piece_tr0_pos0, piece_tr1_pos0])
    sess.commit()

    # same transcibed_id but different position
    piece_tr0_pos1 = TranscriptPiece(transcript=transcript0, position=1)
    sess.add(piece_tr0_pos1)
    sess.commit()

    # test if unique constraint works
    piece_tr0_pos1_2nd = TranscriptPiece(transcript=transcript0, position=1)
    sess.add(piece_tr0_pos1_2nd)
    with pytest.raises(IntegrityError):
        sess.commit()


def test_order_pieces():
    """Add transcript pieces that consists of one or many features to the db and test
    if the order for the transcript pieces and the features is returned according to the
    position property instead of db insertion order.
    """
    sess = mk_memory_session()
    g = Genome(species='Athaliana', version='1.2', acquired_from='Phytozome12')
    coor = Coordinate(seqid='a', length=1000, genome=g)
    sess.add_all([g, coor])
    sess.commit()
    # setup one transcript handler with pieces
    sl, sl_h = setup_data_handler(SuperLocusHandlerBase, SuperLocus)
    t, t_h = setup_data_handler(TranscriptHandlerBase, Transcript, super_locus=sl)
    # insert in wrong order
    piece1 = TranscriptPiece(position=1)
    piece0 = TranscriptPiece(position=0)
    piece2 = TranscriptPiece(position=2)
    t.transcript_pieces = [piece0, piece1, piece2]
    sess.add_all([t, piece1, piece0, piece2])
    sess.commit()
    # see if they can be ordered as expected overall
    op = t_h.sorted_pieces
    print([piece0, piece1, piece2], 'expected')
    print(op, 'sorted')
    assert op == [piece0, piece1, piece2]


def test_fasta_import():
    """Import and test coordinate information from fasta files"""

    def import_fasta(path):
        controller = ImportController(database_path='sqlite:///:memory:')
        controller.add_sequences(path)
        return controller

    # test import of multiple sequences from one file
    controller = import_fasta('testdata/basic_sequences.fa')
    coords = controller.session.query(Coordinate).all()
    assert len(coords) == 5
    assert coords[0].seqid == '1'
    assert coords[0].length == 405
    assert coords[0].sha1 == 'dc6f3ba2b0c08f7d08053837b810f86cbaa06f38'
    assert coords[0].sequence == 'N' * 405
    assert coords[1].seqid == 'abc'
    assert coords[1].length == 808
    assert coords[1].sequence == 'AAGGCCTT' * 101
    assert coords[2].seqid == 'test123'
    assert coords[2].length == 100
    assert coords[2].sequence == 'A' * 100


def test_a_sequence_stored_in_chunks_reads_back_whole_and_in_ranges(tmp_path, monkeypatch):
    """A sequence longer than SEQUENCE_CHUNK_SIZE is stored in several sequence_chunk rows and reads
    back unchanged from a new session, whole and in ranges across the chunk borders."""
    monkeypatch.setattr(orm, 'SEQUENCE_CHUNK_SIZE', 10)
    db_path = f'sqlite:///{tmp_path / "chunks.sqlite3"}'
    sequence = 'ACGTTGCAAC' 'GGTTAACCAT' 'TTGGC'
    session = mk_memory_session(db_path)
    session.add(Coordinate(seqid='a', length=len(sequence), sequence=sequence, genome=Genome()))
    session.commit()

    coordinate = mk_memory_session(db_path).query(Coordinate).one()
    assert [chunk.start for chunk in coordinate.chunks] == [0, 10, 20]
    assert coordinate.sequence == sequence
    for start, end in ((0, 25), (8, 22), (10, 20), (19, 21), (22, 30)):
        assert coordinate.sequence_range(start, end) == sequence[start:end]


def test_dummyloci_errors():
    """The errors recorded for each transcript of dummyloci{.gff|.fa}, and the ranges they mask,
    merged into one geenuff_mask feature where they meet. Transcripts not listed have neither."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff', clean_gff=True)

    expected = {
        # test case 1 - see gff file for more documentation. Every mask reaches int(sqrt(1199)) * 10
        # = 340 of the 1199bp gap to gene_no_ATG: gene_empty and gene_non_coding lie between the two
        # but have no CDS, so they do not bound a mask (see _buffered_span)
        # x1's CDS ends where the transcript ends, and its 21bp of CDS (11-21, 111-120) are below
        # min_cds_length, masking the whole gene with both flanks as for y1 below, which takes in
        # the 3' flank from the CDS end
        'x1': ({types.MISSING_UTR_3P, types.TOO_SHORT_CDS}, [(0, 740)]),
        # y1's spliced CDS (11-21, 111-120, 201-301) is 122bp (not a multiple of 3) and contains
        # a premature stop codon; neither was designed on purpose, this is just what falls out of
        # the arbitrary dummy CDS boundaries chosen to test the other errors above. A wrong reading
        # frame masks the whole gene (0-400) and the flank on both sides: nothing lies before it
        'y1': ({types.TRUNCATED_CDS, types.INFRAME_STOP_CODON}, [(0, 740)]),
        # z1's CDS is the whole transcript, has no start codon, and its single 10bp piece (111-120)
        # is not a multiple of 3 and too short, masking the whole gene, which takes in both missing
        # UTR flanks
        'z1': ({types.MISSING_UTR_5P, types.MISSING_UTR_3P, types.MISSING_START_CODON,
                types.TRUNCATED_CDS, types.TOO_SHORT_CDS}, [(0, 740)]),
        # test case 4 (test cases 2 and 3 are without errors): the whole gene (1599-1800) with 340
        # of the same 1199bp gap as test case 1 before it and none of the 1bp left to the sequence
        # end after it, 1 // 2 being 0
        'y4': ({types.MISSING_START_CODON}, [(1259, 1800)]),
        # test case 5: x5's spliced CDS (40-151, 152-182) is 143bp, not a multiple of 3, masking the
        # gene 0-300 and 249 // 2 = 124 of the 249bp gap to test case 6
        'x5': ({types.TRUNCATED_CDS}, [(0, 424)]),
        # test case 6: x6's CDS starts at 524, 25bp before its exon and transcript, which masks it
        # whole with the same range as y6 below and leaves its wrong starting phase unchecked
        'x6': ({types.CDS_OUTSIDE_EXONS}, [(425, 1060)]),
        # y6's 4bp intron cannot be spliced and its spliced CDS (525-575, 580-600, 700-725) is
        # 98bp, masking the gene 549-750 with 124 of the 249bp gap before it, leaving base 424
        # between test case 5's mask and this one, and int(sqrt(1005)) * 10 = 310 of the 1005bp
        # left to the sequence end after it
        'y6': ({types.TOO_SHORT_INTRON, types.TRUNCATED_CDS}, [(425, 1060)]),
        # test case 7, 199 // 2 = 99 of the 199bp gap (1548-1350) to test case 8
        'x7': ({types.MISSING_UTR_5P}, [(1448, 1349)]),
        # test case 8: x8's overlapping exon and CDS lines mask the whole gene (1749-1549) with 2 of
        # the 5bp (1750-1754) to the end of the sequence before it and 99 of the 199bp (1548-1350)
        # to test case 7 after it. Its CDS is not checked further, the spliced sequence repeating
        # the overlapping bases, so its missing start codon is not reported
        'x8': ({types.OVERLAPPING_EXONS, types.OVERLAPPING_CDS}, [(1751, 1449)]),
    }
    found = {}
    for t in controller.session.query(Transcript):
        errors, masks = error_types(t), feature_ranges(t, types.GEENUFF_MASK)
        if errors or masks:
            found[t.given_name] = (errors, masks)
    assert found == expected


def test_case_1():
    """Confirm the existence of all features of test case 1 of dummyloci.gff except
    for mask features, which are tested in test_dummyloci_errors().
    Does not test the exact ids or exactly matching object relationships."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff', clean_gff=True)
    query = controller.session.query

    sl = query(SuperLocus).filter(SuperLocus.given_name == 'gene0').one()
    sl_h = SuperLocusHandlerBase(sl)

    super_locus = SuperLocus(given_name='gene0', type=types.SuperLocusAll.gene)
    assert matching_super_loci(sl, super_locus)

    coords = query(Coordinate).all()

    # confirm exisistence of all objects where things could go wrong
    # above db level and one piece has to exist for the sl_h.features query to work
    sl_h_features_wo_masks = [f for f in sl_h.features if f.type.value != types.GEENUFF_MASK]
    sl_objects = list(sl_h_features_wo_masks) + sl_h.data.transcripts + sl_h.data.proteins

    # first transcript
    transcript = Transcript(given_name='x1', type=types.TranscriptLevel.mRNA, super_locus=sl)
    assert orm_object_in_list(transcript, sl_objects)

    protein = Protein(given_name='x1.p', super_locus=sl)
    assert orm_object_in_list(protein, sl_objects)

    feature = Feature(given_name='x1',
                      type=types.GeenuffFeature.geenuff_transcript,
                      start=0,
                      end=120,
                      start_is_biological_start=True,
                      end_is_biological_end=False,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.GeenuffFeature.geenuff_cds,
                      start=10,
                      end=120,
                      start_is_biological_start=True,
                      end_is_biological_end=False,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.GeenuffFeature.geenuff_intron,
                      start=21,
                      end=110,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)

    # second transcript
    transcript = Transcript(given_name='y1', type=types.TranscriptLevel.mRNA, super_locus=sl)
    assert orm_object_in_list(transcript, sl_objects)

    protein = Protein(given_name='y1.p', super_locus=sl)
    assert orm_object_in_list(protein, sl_objects)

    feature = Feature(given_name='y1',
                      type=types.GeenuffFeature.geenuff_transcript,
                      start=0,
                      end=400,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.GeenuffFeature.geenuff_cds,
                      start=10,
                      end=301,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.GeenuffFeature.geenuff_intron,
                      start=21,
                      end=110,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.GeenuffFeature.geenuff_intron,
                      start=120,
                      end=200,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)

    # third transcript
    transcript = Transcript(given_name='z1', type=types.TranscriptLevel.mRNA, super_locus=sl)
    assert orm_object_in_list(transcript, sl_objects)

    protein = Protein(given_name='z1.p', super_locus=sl)
    assert orm_object_in_list(protein, sl_objects)

    feature = Feature(given_name='z1',
                      type=types.GeenuffFeature.geenuff_transcript,
                      start=110,
                      end=120,
                      start_is_biological_start=False,
                      end_is_biological_end=False,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.GeenuffFeature.geenuff_cds,
                      start=110,
                      end=120,
                      start_is_biological_start=False,
                      end_is_biological_end=False,
                      is_plus_strand=True,
                      phase=0,
                      coordinate=coords[0])
    assert orm_object_in_list(feature, sl_objects)

    # test if we have no extra objects
    assert not sl_objects


def test_case_8():
    """Analogous to test_case_1()"""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff', clean_gff=True)
    query = controller.session.query

    sl = query(SuperLocus).\
        filter(SuperLocus.given_name == 'gene_overlapping_exons_missing_start').one()
    sl_h = SuperLocusHandlerBase(sl)

    super_locus = SuperLocus(given_name='gene_overlapping_exons_missing_start',
                             type=types.SuperLocusAll.gene)
    assert matching_super_loci(sl, super_locus)

    coords = query(Coordinate).all()

    sl_h_features_wo_masks = [f for f in sl_h.features if f.type.value != types.GEENUFF_MASK]
    sl_objects = list(sl_h_features_wo_masks) + sl_h.data.transcripts + sl_h.data.proteins

    # first transcript
    transcript = Transcript(given_name='x8', type=types.TranscriptLevel.mRNA, super_locus=sl)
    assert orm_object_in_list(transcript, sl_objects)

    protein = Protein(given_name='x8.p', super_locus=sl)
    assert orm_object_in_list(protein, sl_objects)

    # x8's overlapping exon and CDS lines leave no boundary of it trustworthy
    feature = Feature(given_name='x8',
                      type=types.GeenuffFeature.geenuff_transcript,
                      start=1749,
                      end=1548,
                      start_is_biological_start=False,
                      end_is_biological_end=False,
                      is_plus_strand=False,
                      phase=0,
                      coordinate=coords[1])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.GeenuffFeature.geenuff_cds,
                      start=1724,
                      end=1573,
                      start_is_biological_start=False,
                      end_is_biological_end=False,
                      is_plus_strand=False,
                      phase=0,
                      coordinate=coords[1])
    assert orm_object_in_list(feature, sl_objects), (
        '{} not found in \n{}'.format(feature, '\n'.join([str(o) for o in sl_objects if isinstance(o, Feature)])))
    feature = Feature(given_name=None,
                      type=types.GeenuffFeature.geenuff_intron,
                      start=1718,
                      end=1649,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=False,
                      phase=0,
                      coordinate=coords[1])
    assert orm_object_in_list(feature, sl_objects)

    # second transcript
    transcript = Transcript(given_name='y8', type=types.TranscriptLevel.mRNA, super_locus=sl)
    assert orm_object_in_list(transcript, sl_objects)

    protein = Protein(given_name='y8.p', super_locus=sl)
    assert orm_object_in_list(protein, sl_objects)

    feature = Feature(given_name='y8',
                      type=types.GeenuffFeature.geenuff_transcript,
                      start=1749,
                      end=1548,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=False,
                      phase=0,
                      coordinate=coords[1])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.GeenuffFeature.geenuff_cds,
                      start=1729,
                      end=1669,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=False,
                      phase=0,
                      coordinate=coords[1])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.GeenuffFeature.geenuff_intron,
                      start=1658,
                      end=1599,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=False,
                      phase=0,
                      coordinate=coords[1])
    assert orm_object_in_list(feature, sl_objects)

    # test if we have no extra objects
    assert not sl_objects


def test_case_9_a_gene_without_errors_on_the_minus_strand():
    """Test case 9 of dummyloci.gff, a gene without errors on the minus strand, keeps its transcript,
    CDS and intron as given, with no error and no mask, and is selected for export. Minus strand
    features run from their 0-based start down to their exclusive end."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff', clean_gff=True)
    x9 = transcripts_by_name(controller)['x9']

    assert feature_ranges(x9, types.GEENUFF_TRANSCRIPT) == [(899, 699)]
    assert feature_ranges(x9, types.GEENUFF_CDS) == [(870, 729)]
    assert feature_ranges(x9, types.GEENUFF_INTRON) == [(840, 759)]
    assert not error_types(x9)
    assert not feature_ranges(x9, types.GEENUFF_MASK)
    assert x9.longest and x9.super_locus.excluded_from_export is None


def test_non_coding_intron():
    """checks non coding introns are handled correctly (cause no error)"""
    # test data has been simplified from an augustus run that previously resulted in erroneous masks
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/exonexonCDS.fa', 'testdata/exonexonCDS.gff3', clean_gff=True)
    # one gene model on the + strand
    features = controller.session.query(Feature).filter(Feature.coordinate_id == 1).all()
    assert len(features) == 3
    transcript = [f for f in features if f.type.value == types.GEENUFF_TRANSCRIPT][0]
    cds = [f for f in features if f.type.value == types.GEENUFF_CDS][0]
    intron = [f for f in features if f.type.value == types.GEENUFF_INTRON][0]
    assert (transcript.start, transcript.end) == (922, 4056)
    assert (cds.start, cds.end) == (2081, 3695)
    assert (intron.start, intron.end) == (1415, 2041)
    # one gene model on the - strand
    features = controller.session.query(Feature).filter(Feature.coordinate_id == 2).all()
    assert len(features) == 3
    transcript = [f for f in features if f.type.value == types.GEENUFF_TRANSCRIPT][0]
    cds = [f for f in features if f.type.value == types.GEENUFF_CDS][0]
    intron = [f for f in features if f.type.value == types.GEENUFF_INTRON][0]
    assert (transcript.start, transcript.end) == (3459, 334)
    assert (cds.start, cds.end) == (2194, 478)
    assert (intron.start, intron.end) == (3310, 2195)


def test_trans_spliced_gene_strand_not_crashing_and_excluded_from_export():
    """A '?' strand (NCBI's convention for trans-spliced genes) must not crash the importer:
    the gene is still saved as its own super locus, but with no transcripts under it, so
    'transcript.longest' never becomes True for it, and it is excluded from the export query
    that both h5 export and masking rely on (see GeenuffExportController._genome_query)."""
    db_path = 'testdata/trans_spliced.sqlite3'
    if os.path.exists(db_path):
        os.remove(db_path)
    try:
        controller = ImportController(database_path=db_path)
        controller.add_genome('testdata/trans_spliced.fa', 'testdata/trans_spliced.gff3',
                              clean_gff=True)

        # the trans-spliced gene is still saved as its own super locus record...
        sl_names = {sl.given_name for sl in controller.session.query(SuperLocus).all()}
        assert sl_names == {'transspliced1', 'gene2'}
        # ...but with no transcripts under it
        ts_sl = controller.session.query(SuperLocus).filter_by(given_name='transspliced1').one()
        assert ts_sl.transcripts == []

        exporter = GeenuffExportController(db_path, longest=True)
        coord_features = exporter.genome_query(longest_only=True)
        given_names = {f.given_name for fs in coord_features.values() for f in fs}
        assert 'rna_ts1' not in given_names
        assert 'cds_ts1' not in given_names
        assert 'rna2' in given_names
    finally:
        if os.path.exists(db_path):
            os.remove(db_path)


def test_filtered_gff3_export_writes_only_longest_labelled_transcripts(tmp_path):
    """FilteredGff3ExportController.write_filtered_gff3 must reproduce exactly the gene
    models old Helixer's h5 export gives labels from (longest transcript per locus, not masked
    whole). The coordinates/phases it writes back out must exactly round-trip the original GFF3
    input."""
    db_path = str(tmp_path / 'filtered_gff3.sqlite3')
    out_path = str(tmp_path / 'filtered.gff3')
    controller = ImportController(database_path=db_path)
    controller.add_genome('testdata/trans_spliced.fa', 'testdata/trans_spliced.gff3',
                          clean_gff=True)

    exporter = FilteredGff3ExportController(db_path)
    exporter.write_filtered_gff3(out_path)

    with open(out_path) as f:
        lines = [line.rstrip('\n') for line in f if not line.startswith('#')]

    # only gene2's transcript is coding and labeled; the trans-spliced locus has no transcript
    # at all and is absent entirely
    feature_types = [line.split('\t')[2] for line in lines]
    assert feature_types == ['gene', 'mRNA', 'exon', 'CDS']
    for line in lines:
        assert 'gene2' in line or 'rna2' in line

    cds_line = [line for line in lines if line.split('\t')[2] == 'CDS'][0]
    cols = cds_line.split('\t')
    # exactly the original input coordinates and phase (1030-1128, phase 0)
    assert (cols[0], cols[3], cols[4], cols[6], cols[7]) == ('NC_TEST.2', '1030', '1128', '+', '0')

    gene_line = [line for line in lines if line.split('\t')[2] == 'gene'][0]
    cols = gene_line.split('\t')
    assert (cols[3], cols[4]) == ('1000', '1200')


def test_filtered_gff3_export_with_include_erroneous_writes_all_it_can(tmp_path, caplog):
    """include_erroneous writes one transcript per gene whatever is wrong with it, so a
    prediction can be compared gene by gene rather than only against sound ones. The line it
    draws is representability, not correctness: a gene dropped from the h5 export for
    overlapping another is written, nothing being wrong with it beyond sharing sequence, which a
    GFF3 holds without trouble, while a gene whose features cannot be placed on one strand is
    not. Genes left out are logged with the reason, per types.unrepresentable_reasons."""
    def genes_written(fasta, gff3, include_erroneous):
        stem = f'{gff3}_{include_erroneous}'
        db_path = str(tmp_path / f'{stem}.sqlite3')
        out_path = str(tmp_path / f'{stem}.gff3')
        controller = ImportController(database_path=db_path)
        controller.add_genome(f'testdata/{fasta}.fa', f'testdata/{gff3}.gff3', clean_gff=True)
        FilteredGff3ExportController(db_path).write_filtered_gff3(
            out_path, include_erroneous=include_erroneous)
        with open(out_path) as handle:
            return {line.split('ID=')[1].strip() for line in handle
                    if line.split('\t')[2:3] == ['gene']}

    # the default export writes the genes the h5 export gives labels from: geneCrossCoder and
    # geneCleanNoCds, kept in their overlaps and with only the flanks of missing UTRs and the
    # overhang of the gene dropped in their favor masked; the others are masked whole or dropped.
    # With include_erroneous all eight come through, geneCrossGivesWay and geneTruncatedCoder among
    # them although both were dropped from the h5 export for overlapping a gene that was kept
    assert genes_written('overlapping_loci_pairs', 'overlapping_loci_pairs', False) == {
        'geneCrossCoder', 'geneCleanNoCds'}
    assert genes_written('overlapping_loci_pairs', 'overlapping_loci_pairs', True) == {
        'geneCrossCoder', 'geneCrossGivesWay', 'geneOuter', 'geneInner', 'geneBothTruncatedLeft',
        'geneBothTruncatedRight', 'geneCleanNoCds', 'geneTruncatedCoder'}

    # a gene that cannot be placed on one strand stays out even so, and is reported as left out
    caplog.set_level(logging.INFO)
    assert genes_written('unplaceable_strand', 'unplaceable_strand', True) == {'geneOK'}
    assert '2 genes left out, their features are not all on one definite strand' in caplog.text


def test_exon_lines_naming_a_gene_are_left_out_wherever_they_stand():
    """Exon lines naming their gene as the parent instead of a transcript, written above the
    gene's transcript lines as NCBI does. Lines are grouped by Parent, not by position, so where
    they stand does not matter. Each is left out, the one duplicating a line of the gene's
    transcript counted as such, and none is attached to a transcript, not even to a sole one: it
    could equally belong to an isoform the file gives no transcript line of its own, so attaching
    it could merge two proteins into one. See OrganizedGFFEntries._place_pieces; the four cases
    are in the gff3."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/overlapping_loci_chain.fa',
                          'testdata/features_above_transcript.gff3', clean_gff=True)

    # geneEchoAbove's echo; geneSoleTranscript's, geneAmbiguous's and geneNoTranscript's lines
    assert controller.stats.dropped_lines == {'gene_parented_duplicate': {'exon': 1},
                                              'gene_parented': {'exon': 3}}

    # the sole-transcript gene keeps only the exon nested under its transcript (800-900), so its
    # one intron runs from the transcript's start. Had the line above been attached, the intron
    # would instead lie between the two exons, at 700-799
    transcript = controller.session.query(Transcript).filter_by(
        given_name='rnaSoleTranscript').one()
    introns = sorted((f.start, f.end) for p in transcript.transcript_pieces for f in p.features
                     if f.type.value == types.GEENUFF_INTRON)
    assert introns == [(599, 799)]


def test_unrecognized_feature_type_is_skipped_not_fatal():
    """A GFF3 line whose feature type isn't a Sequence Ontology term GeenuFF knows (e.g. a
    stray GenBank-style 'misc_feature') must not crash the whole import: that one line is
    skipped and counted, while the rest of the file still imports normally."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/trans_spliced.fa', 'testdata/unrecognized_feature_type.gff3',
                          clean_gff=True)

    assert controller.stats.unrecognized_feature_types == {'misc_feature': 1}
    # the rest of the file (a normal, valid gene) still imported fine
    sl = controller.session.query(SuperLocus).filter_by(given_name='gene2').one()
    assert len(sl.transcripts) == 1


def test_exon_and_cds_lines_naming_a_gene_are_left_out_not_glued_onto_a_transcript():
    """Some GFF3 sources (e.g. NCBI/EMBL) redundantly echo a transcript's exon/CDS lines a second
    time, naming the gene as the parent. gene1 is a real pattern found in NCBI's TAIR12
    Arabidopsis thaliana annotation (gene HEC1/AT5G67060): its three echoes duplicate lines of
    its transcripts. gene2's gene-parented exon matches neither of its transcripts, and gene3's
    exon and CDS match its only transcript neither. All are left out, none attached to a
    transcript, so no transcript gets a bogus overlap."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/gene_parented_duplicate_exons.fa',
                          'testdata/gene_parented_duplicate_exons.gff3', clean_gff=True)

    assert controller.stats.dropped_lines == {'gene_parented_duplicate': {'exon': 3},
                                              'gene_parented': {'exon': 2},
                                              'gene_parented_cds': {'CDS': 1}}
    # gene3's CDS line, counted per gene, duplicating no CDS line of its transcript
    assert (len(controller.stats.gene_parented_cds_genes), controller.stats.gene_parented_cds_duplicates) == (1, 0)
    found = {e.type.value for e in controller.session.query(TranscriptError)}
    assert types.OVERLAPPING_EXONS not in found

    # rna2's own two real exons are untouched by the dropped gene-parented duplicates
    rna2 = controller.session.query(Transcript).filter_by(given_name='rna2').one()
    introns = [(f.start, f.end) for f in rna2.transcript_pieces[0].features
               if f.type.value == types.GEENUFF_INTRON]
    assert introns == [(108, 129)]  # the one real intron between rna2's two exons

    # rna3's own single real exon is untouched by the dropped, unmatched ambiguous one
    rna3 = controller.session.query(Transcript).filter_by(given_name='rna3').one()
    features = [(f.start, f.end) for p in rna3.transcript_pieces for f in p.features]
    assert features == [(299, 350)]

    # rna5 is gene3's sole transcript and keeps only its own exon (500-520), so the transcript,
    # running on to 560, ends in an intron
    rna5 = controller.session.query(Transcript).filter_by(given_name='rna5').one()
    introns5 = [(f.start, f.end) for f in rna5.transcript_pieces[0].features
                if f.type.value == types.GEENUFF_INTRON]
    assert introns5 == [(520, 560)]


def test_overlapping_loci_in_a_chain_are_each_masked_whole():
    """A pair that cannot be settled by keeping one of the two leaves both genes useless: masking
    only what they share would leave each with a hole through it, teaching no terminus and no
    continuity, so both are masked over their whole length instead. Every locus here overlaps two
    others, so every pair is part of a chain and none of them is even attempted
    (see GFFErrorHandling._compute_overlap_masks). GFF coordinates below are 1-based and
    inclusive, the masks are the GeenuFF 0-based, half-open equivalent.

    geneA (100-999) is a wide container. geneB (200-298) is fully nested inside it and geneC
    (900-1100) overlaps geneA's tail, but geneB and geneC do not touch each other. This covers
    the two ways an overlap can be missed when each locus is only compared to its immediate
    predecessor in sort order: geneA is never anyone's successor, so it would never be examined
    at all, and geneC's only predecessor is geneB, where there is nothing to find.

    geneF (2400-5000), geneG (3900-7000) and geneH (3950-3960) add the case of one locus
    overlapping two partners at once, of a locus nested in two of them, and of a tiny nested
    locus sorting between two wide ones that overlap each other.

    geneD_te_like is a bare 'gene' line with no mRNA/CDS at all (e.g. a transposable-element
    annotation), fully containing the real coding geneE. Since geneD_te_like never writes any
    CDS/exon/intron label of its own, there is no genuine per-base conflict, and geneE must NOT
    be masked just for sitting inside its span; only a CDS-containing transcript is eligible
    to create or receive an overlap mask, matching every other error check in this class."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/overlapping_loci_chain.fa',
                          'testdata/overlapping_loci_chain.gff3', clean_gff=True)
    session = controller.session

    # each is masked over its own whole span (geneA 99-999, geneB 199-298, geneC 899-1100), not
    # over what it shares. Every one of them has its CDS flush against its transcript on both
    # ends, and each mask also runs on into the flank on both sides, toward the closest edge of a
    # gene not overlapping it, which takes in the flanks masked for the missing UTRs: geneA back
    # to 50, 49 of the 99bp to the sequence start, and on to 1089, 90 of the 180bp to geneE, past
    # geneC overlapping it. geneB, enclosed by geneA, reaches back to 100 and on to 538,
    # int(sqrt(601)) * 10 = 240 toward geneC
    for name, mask in [('rnaA', (50, 1089)), ('rnaB', (100, 538)), ('rnaC', (659, 1139)),
                       # the same for the second chain: geneF, geneG and geneH overlap one another,
                       # so each reaches back toward geneE's end (1250) and on toward geneI's start
                       ('rnaF', (2069, 5440)), ('rnaG', (3389, 7004)), ('rnaH', (3439, 4510))]:
        assert types.SL_OVERLAP_ERROR in errors_of(session, name)
        assert masks_of(session, name) == [mask]
    # geneE sits inside a transcript-less 'gene' record: no genuine conflict. geneI starts after
    # geneG ends; directly adjacent or apart is not an overlap
    for name in ['rnaE', 'rnaI']:
        assert types.SL_OVERLAP_ERROR not in errors_of(session, name)


def test_features_on_an_unplaceable_strand_are_kept_but_never_exported(tmp_path):
    """A GFF3 strand of '?' (NCBI's marker for a trans-spliced feature) used to raise out of the
    import as soon as it appeared on an mRNA or CDS line, only the gene line being guarded. Such a
    locus is kept in the database in full, but excluded from every export rather than masked:
    where its features belong is precisely what is unknown, so a mask would either cover the wrong
    strand or cover sequence that is perfectly fine."""
    db_path = str(tmp_path / 'unplaceable_strand.sqlite3')
    controller = ImportController(database_path=db_path)
    controller.add_genome('testdata/unplaceable_strand.fa',
                          'testdata/unplaceable_strand.gff3', clean_gff=True)

    def features_of(given_name):
        transcript = controller.session.query(Transcript).filter_by(given_name=given_name).one()
        return [f for p in transcript.transcript_pieces for f in p.features]

    def excluded(given_name):
        return controller.session.query(SuperLocus).filter_by(given_name=given_name).one() \
                         .excluded_from_export

    for gene, rna in (('geneTransSpliced', 'rnaTransSpliced'),
                      ('geneCdsStrandClash', 'rnaCdsStrandClash')):
        features = features_of(rna)
        # nothing is dropped: the transcript and its CDS are in the database either way
        assert {types.GEENUFF_TRANSCRIPT, types.GEENUFF_CDS} <= {f.type.value for f in features}
        # and nothing is invented for them either, neither a mask nor an intron standing in for
        # the exons that could not be placed, nor an error
        assert not [f for f in features if f.type.value in (types.GEENUFF_MASK, types.GEENUFF_INTRON)]
        assert not errors_of(controller.session, rna)
        assert excluded(gene) == types.UNPLACEABLE_STRAND

    # only the ordinary gene is exported: the two above are excluded outright, and the non-coding
    # one is never selected, no transcript without a CDS being exported to Helixer
    exported = GeenuffExportController(db_path).genome_query(longest_only=True)
    exported_names = {f.given_name for features in exported.values() for f in features}
    assert 'rnaOK' in exported_names
    assert not {'rnaTransSpliced', 'rnaCdsStrandClash', 'rnaNonCoding'} & exported_names

    # the non-coding transcript is still kept in the database, and is not excluded: it needs no
    # marking, having nothing an export would have selected in the first place
    assert types.GEENUFF_TRANSCRIPT in {f.type.value for f in features_of('rnaNonCoding')}
    assert excluded('geneNonCoding') is None
    assert excluded('geneOK') is None
    assert controller.stats.unexported_super_loci == {types.UNPLACEABLE_STRAND: 2}


def test_an_overlapping_pair_keeps_one_locus_whole_where_it_can():
    """Two loci cannot share a base in an export, one label per base being written. Of an
    isolated crossing pair the less damaged locus, then the one with the longer CDS, is kept whole
    and its partner is dropped from exports, wherever either has coding sequence. Where both are
    masked outright, or one lies inside the other, neither is kept and both are masked over their
    whole length.
    See GFFErrorHandling._decide_overlap_pair; the four pairs are documented in the gff3."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/overlapping_loci_pairs.fa',
                          'testdata/overlapping_loci_pairs.gff3', clean_gff=True)

    session = controller.session

    def excluded(given_name):
        return controller.session.query(SuperLocus).filter_by(given_name=given_name).one() \
                         .excluded_from_export

    # every gene here gets super_loci_overlap_error, the kept one of a resolved pair as well
    for name in ['rnaCrossCoder', 'rnaOuter', 'rnaInner', 'rnaBothTruncatedLeft',
                 'rnaBothTruncatedRight', 'rnaCleanNoCds']:
        assert types.SL_OVERLAP_ERROR in errors_of(session, name)

    # crossing: the longer CDS is kept and only the partner's overhang is masked; the kept locus
    # itself carries no mask over the range they share
    assert excluded('geneCrossCoder') is None
    assert excluded('geneCrossGivesWay') == types.OVERLAP_DROPPED
    # geneCrossGivesWay lacks its 3' UTR, so the mask runs on past its own end (652) into the
    # flank, to 725, 73 of the 147bp to geneOuter. geneCrossCoder (99-399) lacks both UTRs: its 5'
    # flank 50-99, 49 of the 99bp to the sequence start, and its 3' flank, 200 of the 400bp to
    # geneOuter, which the overhang takes in
    assert masks_of(session, 'rnaCrossCoder') == [(50, 99), (399, 725)]
    # the dropped gene's own error is recorded, but masks nothing, its features reaching no export
    assert errors_of(session, 'rnaCrossGivesWay') == {types.MISSING_UTR_3P}
    assert masks_of(session, 'rnaCrossGivesWay') == []

    # nested: never decided, so both are masked over their whole length and the flank on both
    # sides, measured past the partner, taking in the flanks of their missing UTRs. geneOuter
    # (799-1399) from 726, 73 of the 147bp back to the dropped geneCrossGivesWay, to 1499, 100 of
    # the 200bp to geneBothTruncatedLeft; geneInner (999-1200) from 826, 173 of the 347bp back to
    # geneCrossGivesWay, to 1390, int(sqrt(399)) * 10 = 190 of the 399bp to geneBothTruncatedLeft
    assert excluded('geneOuter') is None and excluded('geneInner') is None
    assert masks_of(session, 'rnaOuter') == [(726, 1499)]
    assert masks_of(session, 'rnaInner') == [(826, 1390)]

    # both are masked outright, so neither can be kept and both are masked over their whole length
    # and the flank on both sides, measured past the partner overlapping it: geneBothTruncatedLeft
    # (1599-1899) from 1499, 100 of the 200bp to geneOuter, to 2119, 220 of the 500bp to
    # geneCleanNoCds; geneBothTruncatedRight (1849-2149) from 1639, 210 of the 450bp back to
    # geneOuter, to 2274, 125 of the 250bp to geneCleanNoCds
    assert excluded('geneBothTruncatedLeft') is None and excluded('geneBothTruncatedRight') is None
    assert masks_of(session, 'rnaBothTruncatedLeft') == [(1499, 2119)]
    assert masks_of(session, 'rnaBothTruncatedRight') == [(1639, 2274)]

    # the clean locus is kept although only the truncated one has CDS in the range they share. The
    # dropped geneTruncatedCoder is erroneous, so the mask over its overhang (2699-2950) runs on
    # int(sqrt(4350)) * 10 = 650 into the flank toward the end of the sequence. geneCleanNoCds
    # (2399-2699) lacks its 5' UTR: its 5' flank, 125 of the 250bp back to geneBothTruncatedRight
    assert excluded('geneCleanNoCds') is None
    assert excluded('geneTruncatedCoder') == types.OVERLAP_DROPPED
    assert masks_of(session, 'rnaCleanNoCds') == [(2274, 2399), (2699, 3600)]
    # geneTruncatedCoder's CDS is its whole transcript, 301bp and without codons
    assert errors_of(session, 'rnaTruncatedCoder') == {
        types.MISSING_UTR_5P, types.MISSING_UTR_3P, types.MISSING_START_CODON,
        types.MISSING_STOP_CODON, types.TRUNCATED_CDS}
    # its errors leave both ends unknown, although it is dropped and not masked
    features = transcript_features(session, 'rnaTruncatedCoder')
    for feature_type in (types.GEENUFF_TRANSCRIPT, types.GEENUFF_CDS):
        assert [f[2:] for f in features[feature_type]] == [(False, False)]
    # the import statistics count exported genes only: the two of pair 3, not geneTruncatedCoder
    assert controller.stats.errors[types.TRUNCATED_CDS] == 2

    assert controller.stats.overlap_pairs_resolved == 2
    assert controller.stats.overlap_loci_dropped == 2
    assert controller.stats.overlap_pairs_both_masked_outright == 1
    assert controller.stats.overlap_pairs_nested == 1
    assert controller.stats.overlap_pairs_in_chains == 0


@pytest.fixture(scope='module')
def overlap_grading_session():
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/overlap_grading.fa', 'testdata/overlap_grading.gff3',
                          clean_gff=True)
    return controller.session


def excluded_from_export(session, given_name):
    return session.query(SuperLocus).filter_by(given_name=given_name).one().excluded_from_export


def masks_of(session, given_name):
    """The geenuff_mask ranges of one transcript."""
    transcript = session.query(Transcript).filter_by(given_name=given_name).one()
    return feature_ranges(transcript, types.GEENUFF_MASK)


def errors_of(session, given_name):
    """The error types recorded for one transcript."""
    return error_types(session.query(Transcript).filter_by(given_name=given_name).one())


def test_a_too_short_intron_keeps_a_locus_from_being_kept_in_an_overlap(overlap_grading_session):
    """A locus with an intron too short to be spliced is masked outright, so its clean partner is
    kept although its CDS is the shorter one. The dropped locus' overhang (200-1200) is masked and
    runs on 100 of the 200bp to the sequence start, a too short intron leaving neither end of the
    locus trustworthy. See testdata/overlap_grading.gff3."""
    session = overlap_grading_session
    assert excluded_from_export(session, 'geneShortIntron') == types.OVERLAP_DROPPED
    assert excluded_from_export(session, 'geneClean') is None
    assert errors_of(session, 'rnaClean') == {types.SL_OVERLAP_ERROR}
    assert masks_of(session, 'rnaClean') == [(100, 1200)]


def test_a_wrong_starting_phase_does_not_count_against_a_locus(overlap_grading_session):
    """A wrong starting phase in the file changes no label, so it masks nothing and leaves the
    locus as undamaged as its clean partner, the longer CDS then deciding which one is kept.
    See testdata/overlap_grading.gff3."""
    session = overlap_grading_session
    assert excluded_from_export(session, 'genePhaseLong') is None
    assert excluded_from_export(session, 'geneCleanShort') == types.OVERLAP_DROPPED
    assert errors_of(session, 'rnaPhaseLong') == {types.WRONG_PHASE_5P, types.SL_OVERLAP_ERROR}
    # what the clean, dropped geneCleanShort covers beyond it is masked, but nothing for the phase
    assert masks_of(session, 'rnaPhaseLong') == [(1400, 2200)]


def test_overlapping_coding_locus_pairs_are_recorded_without_masking():
    """The super_locus_overlap table records every pair of coding super loci sharing genomic
    range, over all of their coding transcripts rather than only the one each exports, and masks
    nothing. Loci with no CDS anywhere under them take no part: geneD_te_like is a bare gene line
    fully containing the coding geneE, and neither that pair nor a mask on geneE comes of it."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/overlapping_loci_chain.fa',
                          'testdata/overlapping_loci_chain.gff3', clean_gff=True)

    names = dict(controller.session.query(SuperLocus.id, SuperLocus.given_name))
    pairs = {tuple(sorted((names[row.super_locus_id], names[row.partner_id]))): (row.start, row.end)
             for row in controller.session.execute(orm.super_locus_overlap.select())}

    assert pairs == {('geneA', 'geneB'): (199, 298),
                     ('geneA', 'geneC'): (899, 999),
                     ('geneF', 'geneG'): (3899, 5000),
                     ('geneF', 'geneH'): (3949, 3960),
                     ('geneG', 'geneH'): (3949, 3960)}
    # the transcript-less geneD_te_like is in no pair, and leaves geneE unmasked
    assert not any('geneD_te_like' in pair or 'geneE' in pair for pair in pairs)
    assert types.SL_OVERLAP_ERROR not in errors_of(controller.session, 'rnaE')
    # geneI overlaps nothing, so it is in no pair either
    assert not any('geneI' in pair for pair in pairs)


def test_missing_utr_errors_of_a_nested_gene_are_still_counted():
    """A gene nested inside another has no unclaimed sequence to extend a missing-UTR mask
    into, so its range comes out empty there. The finding is still counted in the import
    statistics, which are tallied from the error types detected rather than from what ended up
    masked (see clean_and_insert)."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/overlapping_loci_chain.fa',
                          'testdata/overlapping_loci_chain.gff3', clean_gff=True)

    # every transcript in the file is a single CDS with no UTR on either side
    n_coding = controller.stats.total_coding_transcripts
    assert controller.stats.errors[types.MISSING_UTR_5P] == n_coding
    assert controller.stats.errors[types.MISSING_UTR_3P] == n_coding
    # every overlapping locus here is in a chain, so none of the pairs can be resolved
    assert controller.stats.overlap_pairs_recorded == 5
    assert controller.stats.super_loci_in_overlap_pairs == 6
    assert controller.stats.overlap_pairs_resolved == 0
    assert controller.stats.overlap_pairs_both_masked_outright == 0
    assert controller.stats.overlap_pairs_nested == 0
    assert controller.stats.overlap_pairs_in_chains == 5
    # geneA, geneB and geneC in one chain, geneF, geneG and geneH in the other
    assert controller.stats.overlap_loci_in_chains == 6


def test_cds_starting_phase_is_reset_not_trusted():
    """The file's own CDS starting phase is not trusted: a complete CDS always starts a
    fresh codon (phase 0) by definition, so the persisted phase is always resets as 0
    regardless of what the file says. The file's value is still used to flag a
    WRONG_PHASE_5P diagnostic error when it disagrees, it just no longer corrupts the
    stored phase."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/wrong_phase.fa', 'testdata/wrong_phase.gff3', clean_gff=True)

    cds = controller.session.query(Feature).filter(
        Feature.type == types.GeenuffFeature.geenuff_cds).one()
    assert cds.phase == 0  # not the file's (wrong) phase of 1

    wrong_phase_errors = controller.session.query(TranscriptError).filter(
        TranscriptError.type == types.Errors.wrong_starting_phase).all()
    assert len(wrong_phase_errors) == 1


def test_codon_split_by_an_intron_is_not_reported_as_missing():
    """A start or stop codon can be cut in two by an intron, in which case it exists only in the
    spliced CDS and never as three contiguous genomic bases. The check reads the spliced CDS
    sequence for that reason; reading three bases from the CDS boundary instead reads into the
    intron and reports a codon that is there as missing."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/split_codon.fa', 'testdata/split_codon.gff3', clean_gff=True)

    # the three transcripts splice to ATG AAA, filler codons and TAA, two on the plus strand and one
    # on the minus
    assert controller.stats.total_coding_transcripts == 3
    assert controller.stats.longest_error_free_transcripts == 3
    assert not controller.stats.errors
    assert not controller.session.query(TranscriptError).all()
    assert not controller.session.query(Feature).filter(
        Feature.type == types.GeenuffFeature.geenuff_mask).all()


def test_codon_split_by_an_intron_is_still_reported_when_absent():
    """The mirror of the above: where the split codon is genuinely absent from the spliced CDS,
    the error is still raised, even though the genome does read a valid codon at the boundary of
    the CDS. Each gene here draws exactly one error, the one it is missing."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/split_codon_absent.fa',
                          'testdata/split_codon_absent.gff3', clean_gff=True)

    # the genome reads ATG over the first gene's CDS boundary and TGA over the second's, but
    # both of those codons are half intron, so neither survives into the transcript
    assert errors_of(controller.session, 'rnaNoStart') == {types.MISSING_START_CODON}
    assert errors_of(controller.session, 'rnaNoStop') == {types.MISSING_STOP_CODON}
    assert controller.stats.errors == {types.MISSING_START_CODON: 1, types.MISSING_STOP_CODON: 1}


def test_gff_gen():
    gff_organizer = OrganizedGFFEntries('testdata/testerSl.gff3')
    x = list(gff_organizer._gff_gen())
    assert len(x) == 103  # every feature type is known, NCBI's 'lnc_RNA' included
    assert x[0].type == 'region'
    assert x[-1].type == 'CDS'


def test_gff_useful_gen():
    gff_organizer = OrganizedGFFEntries('testdata/testerSl.gff3')
    x = list(gff_organizer._useful_gff_entries())
    assert len(x) == 100  # started at 103, should drop the 3 region entries
    assert x[0].type == 'gene'
    assert x[-1].type == 'CDS'


def test_gff_grouper():
    """Lines are grouped into genes by ID and Parent. gene0's two exons name rna0, an 'lnc_RNA'
    line, NCBI's spelling of lncRNA, so they are placed under it; the pseudogene gene3754's six
    exons name the gene itself and are left out."""
    gff_organizer = OrganizedGFFEntries('testdata/testerSl.gff3')
    gff_organizer.load_organized_entries()
    n_genes_seqid = {'NC_015438.2': 2, 'NC_015439.2': 2, 'NC_015440.2': 1}
    for seqid, count in n_genes_seqid.items():
        assert len(gff_organizer.organized_entries[seqid]) == count
        for group in gff_organizer.organized_entries[seqid]:
            assert group['super_locus'].type == 'gene'
    assert gff_organizer.stats.dropped_lines == {'gene_parented': {'exon': 6}}


def test_gff_grouper_leaves_out_genes_sharing_an_id():
    """gene1's ID is used by two gene lines, so which of them a transcript names is undecidable:
    both are left out, with their transcripts and those transcripts' exon, CDS and UTR lines, all
    counted under shared_id, the two coding transcripts among them as lost, from one gene ID, the
    two gene lines sharing it. Only gene2 remains."""
    gff_organizer = OrganizedGFFEntries('testdata/discontinuous_gene.gff3')
    gff_organizer.load_organized_entries()
    groups = gff_organizer.organized_entries['NC_TEST.1']
    assert [group['super_locus'].get_ID() for group in groups] == ['gene2']
    assert gff_organizer.stats.dropped_lines == {'shared_id': {'gene': 2}}
    assert gff_organizer.stats.dropped_lines_below == {'shared_id': {'mRNA': 2, 'exon': 6, 'CDS': 6,
                                                                     'five_prime_UTR': 2,
                                                                     'three_prime_UTR': 2}}
    assert gff_organizer.stats.coding_transcripts_dropped == {'shared_id': 2}
    assert gff_organizer.stats._dropped_lines_entry('shared_id')[1].endswith(
        '2 coding transcripts lost, from 1 gene ID')


@pytest.fixture(scope='module')
def grouping_organizer():
    gff_organizer = OrganizedGFFEntries('testdata/grouping.gff3')
    gff_organizer.load_organized_entries()
    return gff_organizer


def groups_by_gene_id(gff_organizer, seqid):
    return {group['super_locus'].get_ID(): group for group in gff_organizer.organized_entries[seqid]}


def test_transcripts_naming_a_missing_parent_share_a_gene_inferred_for_them(grouping_organizer):
    """Transcripts naming the same Parent ID that matches no line, all on one sequence and strand,
    become isoforms of one gene inferred for them, spanning them all. See testdata/grouping.gff3."""
    missing = groups_by_gene_id(grouping_organizer, 'A')['missingGene']
    gene = missing['super_locus']
    assert (gene.start, gene.end, gene.strand) == (101, 300, '+')
    assert {t.get_ID() for t in missing['transcripts']} == {'phantomA', 'phantomB'}
    assert grouping_organizer.stats.genes_inferred_for_missing_parents == 1


def test_an_exon_naming_several_transcripts_is_put_under_each(grouping_organizer):
    """GFF3 allows several parents, as for an exon shared by isoforms. See
    testdata/grouping.gff3."""
    transcripts = groups_by_gene_id(grouping_organizer, 'A')['g3']['transcripts']
    assert {t.get_ID(): [(e.start, e.end) for e in t_entries['exons']]
            for t, t_entries in transcripts.items()} == {'iso1': [(901, 1100)], 'iso2': [(901, 1100)]}


def test_lines_that_cannot_be_placed_are_left_out_and_counted(grouping_organizer):
    """A gene without an ID, transcripts naming several genes or another transcript, transcripts
    naming one missing parent on different sequences or strands, an exon on another sequence than
    its parent, and every exon below a transcript left out, counted under the reason its transcript
    was left out for. g1 and g2 remain without transcripts, and nothing remains on sequence B. See
    testdata/grouping.gff3."""
    groups = groups_by_gene_id(grouping_organizer, 'A')
    assert set(groups) == {'missingGene', 'g1', 'g2', 'g3'}
    assert not groups['g1']['transcripts'] and not groups['g2']['transcripts']
    assert 'B' not in grouping_organizer.organized_entries
    assert grouping_organizer.stats.dropped_lines == {
        'no_id': {'gene': 1},
        'several_genes': {'mRNA': 1},
        'parent_not_gene': {'mRNA': 1},
        'missing_parent_unplaceable': {'mRNA': 4},
        'other_sequence': {'exon': 1},
    }
    # twoGenes's, phantomC's and phantomD's exons
    assert grouping_organizer.stats.dropped_lines_below == {
        'several_genes': {'exon': 1},
        'missing_parent_unplaceable': {'exon': 2},
    }


def test_coding_transcripts_naming_no_parent_stop_the_import(tmp_path):
    """Coding transcripts without a Parent (and Parent= attribute), as gene predictors write them
    without gene lines, stop the import with an error naming them and pointing to add_gene_lines.py."""
    gff_path = tmp_path / 'no_gene.gff3'
    gff_path.write_text('##gff-version 3\n'
                        'CHR\tT\ttranscript\t21\t53\t.\t+\t.\tID=g1.t1;geneID=g1\n'
                        'CHR\tT\tCDS\t21\t53\t.\t+\t0\tParent=g1.t1\n')
    organizer = OrganizedGFFEntries(str(gff_path))
    with pytest.raises(MissingGeneLines, match=r'1 of 1 coding transcript lines \(100\.0%\) name no Parent gene.*'
                                               r'on 1 sequence \(\w+: 1\).*e\.g\. g1\.t1.*add_gene_lines\.py'):
        organizer.load_organized_entries()


def test_non_coding_transcripts_naming_no_parent_are_left_out(tmp_path):
    """A transcript without a Parent and without CDS lines gives no training data: it is left out
    with its exon line, counted under no_gene, while the gene beside it is grouped as usual."""
    gff_path = tmp_path / 'stray_ncrna.gff3'
    gff_path.write_text('##gff-version 3\n'
                        'CHR\ttest\tncRNA\t1\t24\t.\t+\t.\tID=rna1\n'
                        'CHR\ttest\texon\t1\t24\t.\t+\t.\tParent=rna1\n'
                        'CHR\ttest\tgene\t101\t400\t.\t+\t.\tID=g1\n'
                        'CHR\ttest\tmRNA\t101\t400\t.\t+\t.\tID=t1;Parent=g1\n'
                        'CHR\ttest\tCDS\t101\t400\t.\t+\t0\tParent=t1\n')
    organizer = OrganizedGFFEntries(str(gff_path))
    organizer.load_organized_entries()
    assert organizer.stats.dropped_lines == {'no_gene': {'ncRNA': 1}}
    assert organizer.stats.dropped_lines_below == {'no_gene': {'exon': 1}}
    assert [group['super_locus'].get_ID() for group in organizer.organized_entries['CHR']] == ['g1']


@pytest.fixture(scope='module')
def cds_without_transcript_controller():
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/cds_without_transcript.fa', 'testdata/cds_without_transcript.gff3',
                          clean_gff=True)
    return controller


def transcripts_by_name(controller):
    return {t.given_name: t for t in controller.session.query(Transcript)}


def feature_ranges(transcript, feature_type):
    return [(f.start, f.end) for f in transcript.transcript_pieces[0].features if f.type.value == feature_type]


def test_exons_are_built_from_the_cds_and_utr_lines_of_a_transcript_without_exon_lines(
        cds_without_transcript_controller):
    """Touching UTR and CDS lines form one exon, so tA's only intron lies between its two CDS
    lines and its UTRs are not missing. See testdata/cds_without_transcript.gff3."""
    transcripts = transcripts_by_name(cds_without_transcript_controller)
    assert feature_ranges(transcripts['tA'], types.GEENUFF_INTRON) == [(150, 200)]
    assert feature_ranges(transcripts['tA'], types.GEENUFF_TRANSCRIPT) == [(100, 300)]
    assert feature_ranges(transcripts['tA'], types.GEENUFF_CDS) == [(120, 230)]
    assert not error_types(transcripts['tA'])
    assert not feature_ranges(transcripts['tA'], types.GEENUFF_MASK)


def test_cds_lines_without_a_transcript_line_are_masked_whole(cds_without_transcript_controller):
    """CDS lines naming a gene without transcripts, a parent matching no line or no parent each get
    a transcript spanning them, selected for export and masked over at least that span, under the
    gene they name or one made for them. See testdata/cds_without_transcript.gff3."""
    transcripts = transcripts_by_name(cds_without_transcript_controller)
    spans = {'floating_geneC': (920, 1030), 'floating_ghostT': (1200, 1230),
             'floating_cds-lost': (1300, 1340), 'floating_CHR:1371-1380': (1370, 1380)}
    for name, (start, end) in spans.items():
        assert transcripts[name].longest
        assert feature_ranges(transcripts[name], types.GEENUFF_TRANSCRIPT) == [(start, end)]
        assert error_types(transcripts[name]) == {types.FLOATING_CDS}
        masks = feature_ranges(transcripts[name], types.GEENUFF_MASK)
        assert len(masks) == 1 and masks[0][0] <= start and masks[0][1] >= end
    assert transcripts['floating_geneC'].super_locus.given_name == 'geneC'
    genes = {sl.given_name for sl in cds_without_transcript_controller.session.query(SuperLocus)}
    assert genes == {'geneA', 'geneB', 'geneC', 'floating_gene_ghostT', 'floating_gene_cds-lost',
                     'floating_gene_CHR:1371-1380'}
    assert cds_without_transcript_controller.stats.floating_cds_transcripts == 4
    assert cds_without_transcript_controller.stats.floating_cds_genes == 3


def test_lines_without_a_transcript_line_that_are_not_masked_are_left_out(cds_without_transcript_controller):
    """Lines naming a gene that has a transcript, exon and UTR lines naming a parent matching no
    line, and CDS lines to be masked together that lie on opposite strands. geneB keeps rnaB1
    alone; its three CDS lines are counted per gene, one of them, cds-B1-echo, duplicating rnaB1's
    CDS. See testdata/cds_without_transcript.gff3."""
    stats = cds_without_transcript_controller.stats
    assert stats.dropped_lines == {
        'gene_parented_cds': {'CDS': 3},
        'gene_parented': {'five_prime_UTR': 1},
        'unknown_parent': {'three_prime_UTR': 1, 'exon': 1},
        'floating_cds_unplaceable': {'CDS': 2},
    }
    assert (len(stats.gene_parented_cds_genes), stats.gene_parented_cds_duplicates) == (1, 1)
    gene_b = cds_without_transcript_controller.session.query(SuperLocus).filter_by(given_name='geneB').one()
    assert [t.given_name for t in gene_b.transcripts] == ['rnaB1']


def test_add_gene_lines_gives_parentless_transcripts_their_gene(tmp_path):
    """g1.t1 and g1.t2, without a Parent, get a gene line g1 spanning both, written to the default
    name; x.t1, naming itself as Parent, gets the gene line x already in the file instead, no second
    one being written."""
    gff = tmp_path / 'genes.gff3'
    gff.write_text('##gff-version 3\n'
                   'CHR\ttest\tmRNA\t101\t400\t.\t+\t.\tID=g1.t1\n'
                   'CHR\ttest\tCDS\t101\t400\t.\t+\t0\tParent=g1.t1\n'
                   'CHR\ttest\tmRNA\t201\t500\t.\t+\t.\tID=g1.t2\n'
                   'CHR\ttest\tCDS\t201\t500\t.\t+\t0\tParent=g1.t2\n'
                   'CHR\ttest\tgene\t701\t900\t.\t+\t.\tID=x\n'
                   'CHR\ttest\tmRNA\t701\t900\t.\t+\t.\tID=x.t1;Parent=x.t1\n'
                   'CHR\ttest\tCDS\t701\t900\t.\t+\t0\tParent=x.t1\n')
    run_script('add_gene_lines.py', '--gff', str(gff))
    lines = feature_lines(tmp_path / 'genes_with_genes.gff3')
    assert [(c[2], c[3], c[4], c[8]) for c in lines if c[2] == 'gene'] == [('gene', '101', '500', 'ID=g1'),
                                                                           ('gene', '701', '900', 'ID=x')]
    assert {c[8] for c in lines if c[2] == 'mRNA'} == {'ID=g1.t1;Parent=g1', 'ID=g1.t2;Parent=g1',
                                                       'ID=x.t1;Parent=x'}


def test_add_gene_lines_takes_gene_ids_from_an_attribute_or_a_split_of_the_id(tmp_path):
    """a.t1, without a Parent, gives its gene ID in its geneID attribute and gets a new gene line
    geneA; G10-mRNA-1, naming itself as Parent, is split at '-mRNA-' into G10, the ID of a gene line
    already in the file, which becomes its Parent, no second one being written."""
    gff = tmp_path / 'genes.gff3'
    gff.write_text('##gff-version 3\n'
                   'CHR\ttest\tmRNA\t1\t9\t.\t+\t.\tID=a.t1;geneID=geneA\n'
                   'CHR\ttest\tCDS\t1\t9\t.\t+\t0\tParent=a.t1\n'
                   'CHR\ttest\tgene\t11\t19\t.\t+\t.\tID=G10\n'
                   'CHR\ttest\tmRNA\t11\t19\t.\t+\t.\tID=G10-mRNA-1;Parent=G10-mRNA-1\n'
                   'CHR\ttest\tCDS\t11\t19\t.\t+\t0\tParent=G10-mRNA-1\n')
    run_script('add_gene_lines.py', '--gff', str(gff), '--split-at=-mRNA-')
    lines = feature_lines(tmp_path / 'genes_with_genes.gff3')
    assert [(c[3], c[4], c[8]) for c in lines if c[2] == 'gene'] == [('1', '9', 'ID=geneA'), ('11', '19', 'ID=G10')]
    assert [c[8] for c in lines if c[2] == 'mRNA'] == ['ID=a.t1;geneID=geneA;Parent=geneA',
                                                       'ID=G10-mRNA-1;Parent=G10']


def test_add_gene_lines_reads_parents_as_the_import_does(tmp_path):
    """An empty Parent names no parent, and one naming the transcript twice names only itself, as the
    import reads them: both transcripts get a gene."""
    gff = tmp_path / 'genes.gff3'
    gff.write_text('##gff-version 3\n'
                   'CHR\ttest\tmRNA\t101\t400\t.\t+\t.\tID=a.t1;Parent=\n'
                   'CHR\ttest\tCDS\t101\t400\t.\t+\t0\tParent=a.t1\n'
                   'CHR\ttest\tmRNA\t701\t900\t.\t+\t.\tID=b.t1;Parent=b.t1,b.t1\n'
                   'CHR\ttest\tCDS\t701\t900\t.\t+\t0\tParent=b.t1\n')
    run_script('add_gene_lines.py', '--gff', str(gff))
    lines = feature_lines(tmp_path / 'genes_with_genes.gff3')
    # gene a & b created by the default split at "." from the mRNA name
    assert [c[8] for c in lines if c[2] == 'gene'] == ['ID=a', 'ID=b']
    assert {c[8] for c in lines if c[2] == 'mRNA'} == {'ID=a.t1;Parent=a', 'ID=b.t1;Parent=b'}


def test_rename_gff_seqids_renames_to_fasta_ids_and_leaves_out_unmatched_sequences(tmp_path):
    """The GFF3 sequence ID scaffold_7 is a word of the FASTA header 'Chr1 scaffold_7', so its lines
    are renamed to Chr1; lost matches no FASTA sequence and its lines are left out."""
    fasta, gff = tmp_path / 'genome.fa', tmp_path / 'genes.gff3'
    fasta.write_text('>Chr1 scaffold_7\nATGCCCTAACCCCCCCCCCC\n')
    gff.write_text('##gff-version 3\n'
                   'scaffold_7\ttest\tgene\t1\t9\t.\t+\t.\tID=g1\n'
                   'scaffold_7\ttest\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n'
                   'scaffold_7\ttest\tCDS\t1\t9\t.\t+\t0\tParent=t1\n'
                   'lost\ttest\tgene\t1\t9\t.\t+\t.\tID=g2\n'
                   'lost\ttest\tmRNA\t1\t9\t.\t+\t.\tID=t2;Parent=g2\n'
                   'lost\ttest\tCDS\t1\t9\t.\t+\t0\tParent=t2\n')
    run_script('rename_gff_seqids.py', '--fasta', str(fasta), '--gff', str(gff), '--max-missing-percent', '50')
    lines = feature_lines(tmp_path / 'genes_renamed.gff3')
    assert [(c[0], c[2]) for c in lines] == [('Chr1', 'gene'), ('Chr1', 'mRNA'), ('Chr1', 'CDS')]


# a 30 bp sequence holding a spliced gene, its CDS at 1-6 and 17-22: start codon ATG at 1-3, the
# intron's GT at 7-8 and AG at 15-16, stop codon TAA at 20-22, i.e. 4 checks for rename_gff_seqids.py
SPLICED = 'ATGCCC' + 'GT' + 'CCCCCC' + 'AG' + 'CCCTAA' + 'CCCCCCCC'


def spliced_gene(seqid: str, gene_id: str, offset: int = 0) -> str:
    """The GFF3 lines of the gene of SPLICED on the sequence seqid, offset bases into it."""
    return (f'{seqid}\ttest\tgene\t{1 + offset}\t{22 + offset}\t.\t+\t.\tID={gene_id}\n'
            f'{seqid}\ttest\tmRNA\t{1 + offset}\t{22 + offset}\t.\t+\t.\tID={gene_id}.t1;Parent={gene_id}\n'
            f'{seqid}\ttest\tCDS\t{1 + offset}\t{6 + offset}\t.\t+\t0\tParent={gene_id}.t1\n'
            f'{seqid}\ttest\tCDS\t{17 + offset}\t{22 + offset}\t.\t+\t0\tParent={gene_id}.t1\n')


def test_rename_gff_seqids_pairs_by_sequence_where_no_id_or_header_word_matches(tmp_path):
    """scaffoldA is neither a FASTA ID nor a word of a FASTA header; its gene's 4 checks all pass on
    Chr2 and fail on Chr1, so it is paired with Chr2 by sequence and renamed to it."""
    fasta, gff, pairs = tmp_path / 'genome.fa', tmp_path / 'genes.gff3', tmp_path / 'pairs.tsv'
    fasta.write_text('>Chr1 first\n' + 'C' * 30 + '\n>Chr2 second\n' + SPLICED + '\n')
    gff.write_text('##gff-version 3\n' + spliced_gene('scaffoldA', 'g1'))
    run_script('rename_gff_seqids.py', '--fasta', str(fasta), '--gff', str(gff), '--pairs-out', str(pairs))
    assert {c[0] for c in feature_lines(tmp_path / 'genes_renamed.gff3')} == {'Chr2'}
    assert [line.split('\t')[:3] for line in pairs.read_text().splitlines()[1:]] == [['scaffoldA', 'Chr2', 'sequence']]


def test_rename_gff_seqids_leaves_out_a_duplicate_of_a_paired_sequence(tmp_path):
    """Chr2 pairs by ID; old_Chr2 holds the same gene, which fits only Chr2, paired already: it is
    taken for a duplicate sequence removed from the FASTA file, and its lines are left out."""
    fasta, gff = tmp_path / 'genome.fa', tmp_path / 'genes.gff3'
    fasta.write_text('>Chr1\n' + 'C' * 30 + '\n>Chr2\n' + SPLICED + '\n')
    gff.write_text('##gff-version 3\n' + spliced_gene('Chr2', 'g1') + spliced_gene('old_Chr2', 'g2'))
    run_script('rename_gff_seqids.py', '--fasta', str(fasta), '--gff', str(gff))
    lines = feature_lines(tmp_path / 'genes_renamed.gff3')
    assert {(c[0], c[8].split(';')[0]) for c in lines if c[2] == 'gene'} == {('Chr2', 'ID=g1')}


def test_rename_gff_seqids_stops_where_changed_sequences_hold_too_many_genes(tmp_path):
    """The genes of Chr1 end at 70, more than its own 30 bp past the end of the FASTA sequence, so it
    changed after annotation; it holds half the coding transcripts, over the 5% that may be left out,
    so the script stops with an error and writes nothing."""
    fasta, gff, log = tmp_path / 'genome.fa', tmp_path / 'genes.gff3', tmp_path / 'rename.log'
    fasta.write_text('>Chr1\n' + 'C' * 30 + '\n>Chr2\n' + SPLICED + '\n')
    gff.write_text('##gff-version 3\n' + spliced_gene('Chr2', 'g1')
                   + 'Chr1\ttest\tgene\t41\t70\t.\t+\t.\tID=g2\n'
                   + 'Chr1\ttest\tmRNA\t41\t70\t.\t+\t.\tID=g2.t1;Parent=g2\n'
                   + 'Chr1\ttest\tCDS\t41\t70\t.\t+\t0\tParent=g2.t1\n')
    result = subprocess.run([sys.executable, str(Path(__file__).parents[2] / 'scripts' / 'rename_gff_seqids.py'),
                             '--fasta', str(fasta), '--gff', str(gff), '--log-file', str(log)])
    assert result.returncode == 1
    assert not (tmp_path / 'genes_renamed.gff3').exists()
    assert 'too many to leave out' in log.read_text()


def test_rename_gff_seqids_pairs_by_a_word_of_a_tab_separated_header(tmp_path):
    """The FASTA header 'GWHAAAA00000001<TAB>OriSeqID=TD01<TAB>Len=30', as GWH writes it, holds the GFF3
    sequence ID TD01 as a word, split at the tab and the '=', so TD01 is renamed to its FASTA ID."""
    fasta, gff, pairs = tmp_path / 'genome.fa', tmp_path / 'genes.gff3', tmp_path / 'pairs.tsv'
    fasta.write_text('>GWHAAAA00000001\tOriSeqID=TD01\tLen=30\n' + SPLICED + '\n')
    gff.write_text('##gff-version 3\n' + spliced_gene('TD01', 'g1'))
    run_script('rename_gff_seqids.py', '--fasta', str(fasta), '--gff', str(gff), '--pairs-out', str(pairs))
    assert {c[0] for c in feature_lines(tmp_path / 'genes_renamed.gff3')} == {'GWHAAAA00000001'}
    assert pairs.read_text().splitlines()[1].split('\t')[:3] == ['TD01', 'GWHAAAA00000001', 'full_name']


def test_rename_gff_seqids_counts_identical_checks_of_isoforms_once(tmp_path):
    """21 isoforms of one gene share their CDS, whose start and stop codons are not in the sequence:
    counted once, they are 2 checks, too few to judge the pair by ID, which stands, so nothing is
    written."""
    fasta, gff = tmp_path / 'genome.fa', tmp_path / 'genes.gff3'
    fasta.write_text('>Chr1\n' + 'C' * 20 + '\n')
    gff.write_text('##gff-version 3\nChr1\ttest\tgene\t5\t13\t.\t+\t.\tID=g1\n'
                   + ''.join(f'Chr1\ttest\tmRNA\t5\t13\t.\t+\t.\tID=t{i};Parent=g1\n'
                             f'Chr1\ttest\tCDS\t5\t13\t.\t+\t0\tParent=t{i}\n' for i in range(21)))
    run_script('rename_gff_seqids.py', '--fasta', str(fasta), '--gff', str(gff))
    assert not (tmp_path / 'genes_renamed.gff3').exists()


def test_rename_gff_seqids_reads_genes_crossing_the_origin_of_a_circular_sequence(tmp_path):
    """A CDS from 16 to 24 on a 20 bp sequence crosses its origin, as GFF3 writes it: read across the
    origin, its start codon at 16-18 and stop codon at 2-4 are in place, so the pair by ID stands and
    nothing is written, with a warning to check that the sequence is circular."""
    fasta, gff, log = tmp_path / 'genome.fa', tmp_path / 'genes.gff3', tmp_path / 'rename.log'
    fasta.write_text('>Plastid\nCTAACCCCCCCCCCCATGCC\n')
    gff.write_text('##gff-version 3\n'
                   'Plastid\ttest\tgene\t16\t24\t.\t+\t.\tID=g1\n'
                   'Plastid\ttest\tmRNA\t16\t24\t.\t+\t.\tID=t1;Parent=g1\n'
                   'Plastid\ttest\tCDS\t16\t24\t.\t+\t0\tParent=t1\n')
    run_script('rename_gff_seqids.py', '--fasta', str(fasta), '--gff', str(gff), '--log-file', str(log))
    assert not (tmp_path / 'genes_renamed.gff3').exists()
    assert 'read as crossing the origin of a circular sequence' in log.read_text()


def test_rename_gff_seqids_warns_of_a_sequence_region_one_bp_off(tmp_path):
    """A ##sequence-region written 0-based, 0 to 150 for a 150 bp sequence, is 1 bp off, so it does not
    judge the pair; the 5 genes give 20 checks, enough to judge the pair by them, and all pass, so the
    pair by ID stands and nothing is to be renamed, and a warning asks to check the file."""
    fasta, gff, log = tmp_path / 'genome.fa', tmp_path / 'genes.gff3', tmp_path / 'rename.log'
    fasta.write_text('>Chr1\n' + SPLICED * 5 + '\n')
    gff.write_text('##gff-version 3\n##sequence-region Chr1 0 150\n'
                   + ''.join(spliced_gene('Chr1', f'g{i}', 30 * i) for i in range(5)))
    run_script('rename_gff_seqids.py', '--fasta', str(fasta), '--gff', str(gff), '--log-file', str(log))
    assert not (tmp_path / 'genes_renamed.gff3').exists()
    text = log.read_text()
    assert '1 GFF3 sequences have a ##sequence-region 1 bp off' in text
    assert 'their genes fit their FASTA sequences: nothing to rename' in text


def test_rename_gff_seqids_leaves_out_lines_whose_parents_lie_on_unpaired_sequences(tmp_path):
    """The gene g2 and its mRNA t2 lie on lost, which matches no FASTA sequence, but a CDS line naming
    t2 lies on scaffold_7: it is left out with them, as the import leaves out a line on another
    sequence than its parent."""
    fasta, gff = tmp_path / 'genome.fa', tmp_path / 'genes.gff3'
    fasta.write_text('>Chr1 scaffold_7\nATGCCCTAACCCCCCCCCCC\n')
    gff.write_text('##gff-version 3\n'
                   'scaffold_7\ttest\tgene\t1\t9\t.\t+\t.\tID=g1\n'
                   'scaffold_7\ttest\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n'
                   'scaffold_7\ttest\tCDS\t1\t9\t.\t+\t0\tParent=t1\n'
                   'scaffold_7\ttest\tCDS\t12\t20\t.\t+\t0\tParent=t2\n'
                   'lost\ttest\tgene\t1\t9\t.\t+\t.\tID=g2\n'
                   'lost\ttest\tmRNA\t1\t9\t.\t+\t.\tID=t2;Parent=g2\n')
    run_script('rename_gff_seqids.py', '--fasta', str(fasta), '--gff', str(gff), '--max-missing-percent', '50')
    lines = feature_lines(tmp_path / 'genes_renamed.gff3')
    assert [(c[0], c[2], c[8]) for c in lines] == [('Chr1', 'gene', 'ID=g1'), ('Chr1', 'mRNA', 'ID=t1;Parent=g1'),
                                                   ('Chr1', 'CDS', 'Parent=t1')]


def test_filter_sequence_lengths_removes_a_short_sequence_and_its_lines(tmp_path):
    """With the hard floor alone, lowered to 10 bp, the 5 bp sequence tiny is removed from the FASTA file
    and its lines from the GFF3 file, written to the default names; the longer sequences and their lines
    stay."""
    fasta, gff = tmp_path / 'genome.fa', tmp_path / 'genes.gff3'
    fasta.write_text(''.join(f'>long{i}\n' + 'A' * 20 + '\n' for i in range(3)) + '>tiny\n' + 'A' * 5 + '\n')
    gff.write_text('##gff-version 3\n'
                   '##sequence-region tiny 1 5\n'
                   'long0\ttest\tgene\t1\t5\t.\t+\t.\tID=g1\n'
                   'tiny\ttest\tgene\t1\t5\t.\t+\t.\tID=g2\n')
    run_script('filter_sequence_lengths.py', '--fasta', str(fasta), '--gff', str(gff), '--no-adaptive',
               '--min-seq-len', '10')
    out_fasta = (tmp_path / 'genome_length_filtered.fa').read_text()
    assert [line for line in out_fasta.splitlines() if line.startswith('>')] == ['>long0', '>long1', '>long2']
    out_gff = (tmp_path / 'genes_length_filtered.gff3').read_text()
    assert 'tiny' not in out_gff and [c[8] for c in feature_lines(tmp_path / 'genes_length_filtered.gff3')] == ['ID=g1']


def test_filter_sequence_lengths_removes_a_cluster_of_small_sequences_and_those_below_the_floor(tmp_path):
    """Of six 200-250 bp sequences, two of 15 and 18 bp and one of 5 bp, the two small ones form a
    cluster below the largest gap between the log10 lengths, standing out from the other gaps and
    holding under 5% of the bases, and are removed although above the floor, lowered to 10 bp; the
    5 bp one is removed by the floor. Their GFF3 lines go with them."""
    lengths = {f'chr{i}': 200 + 10 * i for i in range(6)} | {'small1': 15, 'small2': 18, 'tiny': 5}
    fasta, gff = tmp_path / 'genome.fa', tmp_path / 'genes.gff3'
    fasta.write_text(''.join(f'>{seq_id}\n' + 'A' * length + '\n' for seq_id, length in lengths.items()))
    gff.write_text('##gff-version 3\n' + ''.join(f'{seq_id}\ttest\tgene\t1\t5\t.\t+\t.\tID=g_{seq_id}\n'
                                                 for seq_id in lengths))
    run_script('filter_sequence_lengths.py', '--fasta', str(fasta), '--gff', str(gff), '--min-seq-len', '10')
    out_fasta = (tmp_path / 'genome_length_filtered.fa').read_text()
    assert [line[1:] for line in out_fasta.splitlines() if line.startswith('>')] == [f'chr{i}' for i in range(6)]
    assert [c[0] for c in feature_lines(tmp_path / 'genes_length_filtered.gff3')] == [f'chr{i}' for i in range(6)]


def test_filter_sequence_lengths_leaves_out_lines_whose_parents_lie_on_removed_sequences(tmp_path):
    """The gene g2 and its mRNA t2 lie on the removed sequence tiny, but a CDS line naming t2 lies on
    long0: it is left out as well, as the import leaves out a line on another sequence than its
    parent. An exon naming both t1 on long0 and t2 stays, with t1 alone as its Parent."""
    fasta, gff = tmp_path / 'genome.fa', tmp_path / 'genes.gff3'
    fasta.write_text(''.join(f'>long{i}\n' + 'A' * 20 + '\n' for i in range(3)) + '>tiny\n' + 'A' * 5 + '\n')
    gff.write_text('##gff-version 3\n'
                   'long0\ttest\tgene\t1\t9\t.\t+\t.\tID=g1\n'
                   'long0\ttest\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n'
                   'long0\ttest\texon\t1\t9\t.\t+\t.\tParent=t1,t2\n'
                   'long0\ttest\tCDS\t10\t15\t.\t+\t0\tParent=t2\n'
                   'tiny\ttest\tgene\t1\t5\t.\t+\t.\tID=g2\n'
                   'tiny\ttest\tmRNA\t1\t5\t.\t+\t.\tID=t2;Parent=g2\n')
    run_script('filter_sequence_lengths.py', '--fasta', str(fasta), '--gff', str(gff), '--no-adaptive',
               '--min-seq-len', '10')
    lines = feature_lines(tmp_path / 'genes_length_filtered.gff3')
    assert [(c[0], c[2], c[8]) for c in lines] == [('long0', 'gene', 'ID=g1'), ('long0', 'mRNA', 'ID=t1;Parent=g1'),
                                                   ('long0', 'exon', 'Parent=t1')]


def test_add_transcript_lines_gives_genes_without_one_an_mrna(tmp_path):
    """add_transcript_lines.py gives g1, whose exon and CDS lines name the gene itself, an mRNA line
    spanning it as their Parent, which then imports as a gene giving labels instead of masked
    floating CDS lines; g2, having a transcript line, stays as it is."""
    gff, out = tmp_path / 'genes_only.gff3', tmp_path / 'with_transcripts.gff3'
    gff.write_text('##gff-version 3\n'
                   'CHR\ttest\tgene\t101\t400\t.\t+\t.\tID=g1\n'
                   'CHR\ttest\texon\t101\t400\t.\t+\t.\tParent=g1\n'
                   'CHR\ttest\tCDS\t101\t400\t.\t+\t0\tParent=g1\n'
                   'CHR\ttest\tgene\t601\t900\t.\t+\t.\tID=g2\n'
                   'CHR\ttest\tmRNA\t601\t900\t.\t+\t.\tID=t2;Parent=g2\n'
                   'CHR\ttest\texon\t601\t900\t.\t+\t.\tParent=t2\n'
                   'CHR\ttest\tCDS\t601\t900\t.\t+\t0\tParent=t2\n')
    run_script('add_transcript_lines.py', '--gff', str(gff), '-o', str(out))
    lines = feature_lines(out)
    assert [(columns[2], columns[3], columns[4], columns[8]) for columns in lines[:4]] == [
        ('gene', '101', '400', 'ID=g1'), ('mRNA', '101', '400', 'ID=g1_mRNA;Parent=g1'),
        ('exon', '101', '400', 'Parent=g1_mRNA'), ('CDS', '101', '400', 'Parent=g1_mRNA')]
    assert out.read_text().splitlines()[5:] == gff.read_text().splitlines()[4:]

    fasta = tmp_path / 'genes_only.fa'
    fasta.write_text('>CHR\n' + 'C' * 100 + 'ATG' + 'C' * 294 + 'TAA' + 'C' * 200 + 'ATG' + 'C' * 294 + 'TAA'
                     + 'C' * 100 + '\n')
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome(str(fasta), str(out), clean_gff=True)
    assert controller.stats.exported_outcomes['floating'] == 0
    assert controller.stats.labeled_genes() == 2


def test_coding_transcripts_past_a_sequence_edge_are_clipped_and_masked_whole():
    """A coding transcript reaching past its sequence keeps the part on it, selected for export
    and masked over at least that part; one with no CDS line on the sequence is left out with its
    gene. See testdata/sequence_edge.gff3."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/sequence_edge.fa', 'testdata/sequence_edge.gff3', clean_gff=True)
    transcripts = transcripts_by_name(controller)

    for name, span, cds in [('tStart', (0, 60), (10, 50)), ('tEdge', (100, 200), (120, 200))]:
        assert transcripts[name].longest
        assert feature_ranges(transcripts[name], types.GEENUFF_TRANSCRIPT) == [span]
        assert feature_ranges(transcripts[name], types.GEENUFF_CDS) == [cds]
        assert error_types(transcripts[name]) == {types.BEYOND_SEQUENCE_EDGE}
        masks = feature_ranges(transcripts[name], types.GEENUFF_MASK)
        assert len(masks) == 1 and masks[0][0] <= span[0] and masks[0][1] >= span[1]
    assert 'tGone' not in transcripts
    gene_gone = controller.session.query(SuperLocus).filter_by(given_name='geneGone').one()
    assert gene_gone.excluded_from_export == types.OUTSIDE_SEQUENCE
    assert controller.stats.transcripts_outside_sequence_dropped == 1
    assert controller.stats.errors[types.BEYOND_SEQUENCE_EDGE] == 2


def test_a_cds_missing_only_its_stop_codon_gets_the_one_after_it():
    """GTF leaves the stop codon out of the CDS, and GFF3 converted from GTF can keep it that way
    (GTF itself is not read). A CDS of whole codons without any stop codon is extended by the next
    3 bases of its exons where they are a stop codon, also across an intron and on the minus
    strand, and by the 3 bases right after its last exon where the exons end with the CDS, the
    exon and the transcript then extended with it; it is left as it is where they are none or its
    frame is off. See testdata/stop_codon.gff3."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/stop_codon.fa', 'testdata/stop_codon.gff3', clean_gff=True)
    transcripts = transcripts_by_name(controller)

    expected = {'tPlus': (20, 80), 'tSplit': (171, 261), 'tMinus': (396, 336),
                'tNoStop': (470, 530), 'tFrame': (620, 682), 'tCdsOnly': (760, 820), 'tExonEnd': (936, 876),
                'tShort': (1020, 1077), 'tExact': (1170, 1230)}
    assert {name: feature_ranges(t, types.GEENUFF_CDS)[0] for name, t in transcripts.items()} == expected
    # 57 bp of CDS in the file and 60 once extended, so none of them is too short either
    for name in ['tPlus', 'tSplit', 'tMinus']:
        assert not error_types(transcripts[name])
    # the stop codon recovered, only the UTRs their exons lack remain
    assert error_types(transcripts['tCdsOnly']) == {types.MISSING_UTR_5P, types.MISSING_UTR_3P}
    assert error_types(transcripts['tExonEnd']) == {types.MISSING_UTR_3P}
    assert feature_ranges(transcripts['tCdsOnly'], types.GEENUFF_TRANSCRIPT) == [(760, 820)]
    assert feature_ranges(transcripts['tExonEnd'], types.GEENUFF_TRANSCRIPT) == [(949, 876)]
    assert types.MISSING_STOP_CODON in error_types(transcripts['tNoStop'])
    assert types.TRUNCATED_CDS in error_types(transcripts['tFrame'])
    assert controller.stats.stop_codons_recovered == 5
    assert controller.stats.stop_codons_recovered_past_exons == 2


def test_a_cds_shorter_than_60_bp_masks_its_gene_whole():
    """A spliced CDS below min_cds_length (60 bp) is too_short_cds, masking the whole gene with both
    flanks like any error masking the coding sequence, while a CDS of exactly 60 bp is fine. See
    testdata/stop_codon.gff3: gShort (bases 1010-1099, 57 bp of CDS) gets the flank of
    min(gap // 2, int(sqrt(gap)) * 10) on each side, 95 bp of the 190 bp gap back to tCdsOnly's
    extended end at 820 and 30 bp of the 60 bp gap to gExact at 1160; gExonEnd lies in the first gap
    but on the minus strand, so it bounds nothing."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/stop_codon.fa', 'testdata/stop_codon.gff3', clean_gff=True)
    transcripts = transcripts_by_name(controller)

    short = transcripts['tShort']
    assert error_types(short) == {types.TOO_SHORT_CDS}
    assert feature_ranges(short, types.GEENUFF_MASK) == [(915, 1130)]
    # exported, masked outright, so it would never be kept over an overlapping gene
    assert short.longest and short.super_locus.excluded_from_export is None
    assert error_severity(error_types(short)) == 'masked_outright'
    assert controller.stats.errors[types.TOO_SHORT_CDS] == 1

    exact = transcripts['tExact']
    assert not error_types(exact)
    assert not feature_ranges(exact, types.GEENUFF_MASK)


def test_the_minimum_cds_length_comes_from_the_config():
    """min_cds_length can be set like min_intron_length: at 57 bp, gShort's CDS of 57 bp is not
    shorter than it."""
    controller = ImportController(database_path='sqlite:///:memory:', config={'min_cds_length': 57})
    controller.add_genome('testdata/stop_codon.fa', 'testdata/stop_codon.gff3', clean_gff=True)
    short = transcripts_by_name(controller)['tShort']
    assert not error_types(short)
    assert not feature_ranges(short, types.GEENUFF_MASK)


# section: import summary
# from the GFF gene lines on, each step's count follows from the one before
FUNNEL = SUMMARY_STEPS[1:]


@pytest.mark.parametrize('name', ['cds_outside_exons', 'cds_without_transcript', 'dummyloci', 'exonexonCDS',
                                  'exporting', 'gene_parented_duplicate_exons', 'hanging_intron', 'helixer_e2e',
                                  'intron_at_end', 'overlap_grading', 'overlapping_loci_chain',
                                  'overlapping_loci_pairs', 'sequence_edge', 'split_codon_absent', 'split_codon',
                                  'stop_codon', 'trans_spliced', 'unplaceable_strand', 'wrong_phase'])
def test_import_summary_funnel_adds_up(name):
    """Each step of the import summary, its count plus what is added and less what is subtracted
    below it, is the next step's count, from the GFF gene lines down to the genes giving labels,
    whose breakdown adds up to them."""
    gff = f'testdata/{name}.gff3' if os.path.exists(f'testdata/{name}.gff3') else f'testdata/{name}.gff'
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome(f'testdata/{name}.fa', gff, clean_gff=True)
    steps = {header: (count, entries) for header, count, _, entries, _ in controller.stats.summary_steps()}
    for step, following in zip(FUNNEL, FUNNEL[1:]):
        count, entries = steps[step]
        change = sum({'+': n, '-': -n}.get(sign, 0) for sign, n, _, _ in entries)
        assert count + change == steps[following][0], f'{step} -> {following}'
    labeled, breakdown = steps[GENES_GIVING_LABELS]
    assert sum(n for _, n, _, _ in breakdown) == labeled


def test_a_transcript_with_nested_exons_is_selected_and_masked_in_full(tmp_path):
    """Exon lines each running from the transcript start to the end of one CDS line overlap one
    another; the gene's only coding transcript is still its longest, exported and masked in full
    for overlapping_exons."""
    fasta, gff = tmp_path / 'nested.fa', tmp_path / 'nested.gff3'
    fasta.write_text('>CHR\n' + 'C' * 1000 + '\n')
    gff.write_text('##gff-version 3\n'
                   'CHR\ttest\tgene\t101\t700\t.\t+\t.\tID=g1\n'
                   'CHR\ttest\tmRNA\t101\t700\t.\t+\t.\tID=t1;Parent=g1\n'
                   + ''.join(f'CHR\ttest\texon\t101\t{end}\t.\t+\t.\tParent=t1\n'
                             for end in (700, 600, 500, 400, 300))
                   + ''.join(f'CHR\ttest\tCDS\t{start}\t{end}\t.\t+\t0\tParent=t1\n'
                             for start, end in ((101, 300), (351, 400), (451, 500), (551, 600), (651, 700))))
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome(str(fasta), str(gff), clean_gff=True)
    assert controller.session.query(Transcript).filter(Transcript.longest.is_(True)).count() == 1
    assert controller.stats.longest_transcripts == 1
    assert controller.stats.errors_masking_in_full == {types.OVERLAPPING_EXONS: 1}


def test_import_summary_counts_the_masked_bases_within_exported_genes_only():
    """Of the three exported 600 bp genes of cds_outside_exons, rnaUtrOnly is masked in full with
    both flanks: its 600 bases count as masked gene bases, the flanks beside it not."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/cds_outside_exons.fa', 'testdata/cds_outside_exons.gff3', clean_gff=True)
    extra = {header: extra for header, _, _, _, extra in controller.stats.summary_steps()}
    assert any(line.startswith('600 of the 1800 bases of the exported genes masked (33.33%)')
               for line in extra[GENES_GIVING_LABELS])


# section: types
def test_gff_types_are_imported_or_skipped():
    """Every known GFF3 type is either imported or skipped: an SO term GeenuFF has no use for, such
    as region, and Augustus' tss are skipped; exon and the non-SO lnc_RNA and five_prime_utr are
    imported."""
    assert {'region', 'tss'} <= types.IGNORABLE_GFF_TYPES
    assert {'exon', 'lnc_RNA', 'five_prime_utr'} <= types.USEFUL_GFF_TYPES
    assert not types.IGNORABLE_GFF_TYPES & types.USEFUL_GFF_TYPES
    assert types.KNOWN_GFF_TYPES == types.IGNORABLE_GFF_TYPES | types.USEFUL_GFF_TYPES


def test_enums_name_val_match():
    for level in (types.SuperLocusAll, types.TranscriptLevel, types.ExonLevel, types.CDSLevel, types.UTRLevel):
        for x in level:
            assert x.name == x.value


# section: helpers
def test_seqid_matching_pairs_gff_ids_with_fasta_ids():
    # the GFF covering only some of the FASTA sequences
    assert helpers.match_seqids({'a', 'b', 'c'}, {'a', 'c'}) == {'a': 'a', 'c': 'c'}
    # GFF IDs abbreviated from the FASTA IDs
    assert helpers.match_seqids({'a.seq', 'b.seq', 'c.seq'}, {'a', 'b'}) == {'a': 'a.seq', 'b': 'b.seq'}
    # an identical FASTA ID is taken over others containing the GFF ID
    assert helpers.match_seqids({'1', '10', '2.seq'}, {'1', '2'}) == {'1': '1', '2': '2.seq'}


def test_seqid_matching_raises_for_gff_ids_without_one_fasta_id():
    # a GFF ID with no FASTA sequence
    with pytest.raises(helpers.NonMatchableIDs, match=r'matching no FASTA header \(1\): d'):
        helpers.match_seqids({'a', 'b', 'c'}, {'a', 'b', 'c', 'd'})
    # GFF IDs contained in several FASTA IDs
    with pytest.raises(helpers.NonMatchableIDs, match=r'contained in several FASTA headers \(2\)'):
        helpers.match_seqids({'ab.seq', 'ba.seq', 'c.seq', 'a.seq', 'b.seq'}, {'a', 'b', 'c'})
    # two GFF IDs contained in the same FASTA ID
    with pytest.raises(helpers.NonMatchableIDs, match=r'containing several GFF sequence IDs \(1\): ab \(a, b\)'):
        helpers.match_seqids({'ab'}, {'a', 'b'})


@pytest.mark.parametrize('names, fasta_text, gff_text, message', [
    # each file named like the other kind
    (('genes.gff3', 'genome.fa'), '>Chr1\nACGT\n', 'Chr1\tt\tgene\t1\t4\t.\t+\t.\tID=g1\n',
     'look swapped by their names'),
    # a FASTA file not starting with '>'
    (('genome.fa', 'genes.gff3'), 'ACGT\n', 'Chr1\tt\tgene\t1\t4\t.\t+\t.\tID=g1\n', 'does not start with ">"'),
    # a sequence ID holding a tab, giving the first feature line 10 columns
    (('genome.fa', 'genes.gff3'), '>Chr1\nACGT\n', '##gff-version 3\nChr\t1\tt\tgene\t1\t4\t.\t+\t.\tID=g1\n',
     'has 10 tab-separated columns instead of 9'),
])
def test_require_fasta_and_gff_stops_on_files_not_of_their_kind(tmp_path, capsys, names, fasta_text, gff_text, message):
    fasta, gff = tmp_path / names[0], tmp_path / names[1]
    fasta.write_text(fasta_text)
    gff.write_text(gff_text)
    with pytest.raises(SystemExit):
        helpers.require_fasta_and_gff(argparse.ArgumentParser(), str(fasta), str(gff))
    assert message in capsys.readouterr().err


def test_gff_to_seqids():
    x = gffhelper.GFFParser.seqids('testdata/testerSl.gff3')
    assert x == {'NC_015438.2', 'NC_015439.2', 'NC_015440.2'}
