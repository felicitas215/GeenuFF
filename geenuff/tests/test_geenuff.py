import os
import logging
from collections import defaultdict

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
from geenuff.applications.importer import ImportController, InsertCounterHolder, OrganizedGFFEntries
from geenuff.applications.exporter import GeenuffExportController
from geenuff.applications.exporters.gff3 import FilteredGff3ExportController


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


def test_dummyloci_errors():
    """The errors recorded for each transcript of dummyloci{.gff|.fa}, and the ranges they mask,
    merged into one geenuff_mask feature where they meet. Transcripts not listed have neither."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/dummyloci.fa', 'testdata/dummyloci.gff', clean_gff=True)

    expected = {
        # test case 1 - see gff file for more documentation. Every mask reaches int(sqrt(1199)) * 10
        # = 340 of the 1199bp gap to gene_no_ATG: gene_empty and gene_non_coding lie between the two
        # but have no CDS, so they do not bound a mask (see _buffered_span)
        # x1's CDS ends where the transcript ends: the 3' flank from the CDS end
        'x1': ({types.MISSING_UTR_3P}, [(120, 740)]),
        # y1's spliced CDS (11-21, 111-120, 201-301) is 122bp (not a multiple of 3) and contains
        # a premature stop codon; neither was designed on purpose, this is just what falls out of
        # the arbitrary dummy CDS boundaries chosen to test the other errors above. A wrong reading
        # frame masks the whole gene (0-400) and the flank on both sides: nothing lies before it
        'y1': ({types.TRUNCATED_CDS, types.INFRAME_STOP_CODON}, [(0, 740)]),
        # z1's CDS is the whole transcript, has no start codon, and its single 10bp piece (111-120)
        # is not a multiple of 3, masking the whole gene, which takes in both missing UTR flanks
        'z1': ({types.MISSING_UTR_5P, types.MISSING_UTR_3P, types.MISSING_START_CODON,
                types.TRUNCATED_CDS}, [(0, 740)]),
        # test case 4 (test cases 2 and 3 are without errors): the whole gene (1599-1800) with 340
        # of the same 1199bp gap as test case 1 before it and none of the 1bp left to the sequence
        # end after it, 1 // 2 being 0
        'y4': ({types.MISSING_START_CODON}, [(1259, 1800)]),
        # test case 5: x5's spliced CDS (40-151, 152-182) is 143bp, not a multiple of 3, masking the
        # gene 0-300 and 249 // 2 = 124 of the 249bp gap to test case 6
        'x5': ({types.TRUNCATED_CDS}, [(0, 424)]),
        # test case 6: x6's starting phase is only recorded, the importer setting the phase itself
        'x6': ({types.WRONG_PHASE_5P}, []),
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
    assert orm_object_in_list(feature, sl_objects), '{} not found in \n{}'.format(feature, '\n'.join([str(o) for o in sl_objects if isinstance(o, Feature)]))
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
                      end=1699,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=False,
                      phase=0,
                      coordinate=coords[1])
    assert orm_object_in_list(feature, sl_objects)
    feature = Feature(given_name=None,
                      type=types.GeenuffFeature.geenuff_intron,
                      start=1678,
                      end=1599,
                      start_is_biological_start=True,
                      end_is_biological_end=True,
                      is_plus_strand=False,
                      phase=0,
                      coordinate=coords[1])
    assert orm_object_in_list(feature, sl_objects)

    # test if we have no extra objects
    assert not sl_objects


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
    assert (controller.stats.gene_parented_cds_genes, controller.stats.gene_parented_cds_duplicates) == (1, 0)
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

    # the three transcripts splice to ATG AAA TAA, two on the plus strand and one on the minus
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
    assert len(x) == 102  # started at 103, one 'lnc_RNA' line is not a known feature type
    assert x[0].type == 'region'
    assert x[-1].type == 'CDS'


def test_gff_useful_gen():
    gff_organizer = OrganizedGFFEntries('testdata/testerSl.gff3')
    x = list(gff_organizer._useful_gff_entries())
    assert len(x) == 99  # started at 102, should drop the 3 region entries
    assert x[0].type == 'gene'
    assert x[-1].type == 'CDS'


def test_gff_grouper():
    """Lines are grouped into genes by ID and Parent. gene0's two exons name rna0, whose
    'lnc_RNA' line is skipped as a feature type GeenuFF does not know, so their parent matches no
    line; the pseudogene gene3754's six exons name the gene itself."""
    gff_organizer = OrganizedGFFEntries('testdata/testerSl.gff3')
    gff_organizer.load_organized_entries()
    n_genes_seqid = {'NC_015438.2': 2, 'NC_015439.2': 2, 'NC_015440.2': 1}
    for seqid, count in n_genes_seqid.items():
        assert len(gff_organizer.organized_entries[seqid]) == count
        for group in gff_organizer.organized_entries[seqid]:
            assert group['super_locus'].type == 'gene'
    assert gff_organizer.stats.dropped_lines == {'unknown_parent': {'exon': 2},
                                                 'gene_parented': {'exon': 6}}


def test_gff_grouper_leaves_out_genes_sharing_an_id():
    """gene1's ID is used by two gene lines, so which of them a transcript names is undecidable:
    both are left out, with their transcripts and those transcripts' exon, CDS and UTR lines, all
    counted under shared_id, the two coding transcripts among them as lost. Only gene2 remains."""
    gff_organizer = OrganizedGFFEntries('testdata/discontinuous_gene.gff3')
    gff_organizer.load_organized_entries()
    groups = gff_organizer.organized_entries['NC_TEST.1']
    assert [group['super_locus'].get_ID() for group in groups] == ['gene2']
    assert gff_organizer.stats.dropped_lines == {'shared_id': {'gene': 2}}
    assert gff_organizer.stats.dropped_lines_below == {'shared_id': {'mRNA': 2, 'exon': 6, 'CDS': 6,
                                                                     'five_prime_UTR': 2,
                                                                     'three_prime_UTR': 2}}
    assert gff_organizer.stats.coding_transcripts_dropped == {'shared_id': 2}


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
    """A gene without an ID, transcripts naming no parent, several genes or another transcript,
    transcripts naming one missing parent on different sequences or strands, an exon on another
    sequence than its parent, and every exon below a transcript left out, counted under the reason
    its transcript was left out for. g1 and g2 remain without transcripts, and nothing remains on
    sequence B. See testdata/grouping.gff3."""
    groups = groups_by_gene_id(grouping_organizer, 'A')
    assert set(groups) == {'missingGene', 'g1', 'g2', 'g3'}
    assert not groups['g1']['transcripts'] and not groups['g2']['transcripts']
    assert 'B' not in grouping_organizer.organized_entries
    assert grouping_organizer.stats.dropped_lines == {
        'no_id': {'gene': 1},
        'no_parent': {'mRNA': 1},
        'several_genes': {'mRNA': 1},
        'parent_not_gene': {'mRNA': 1},
        'missing_parent_unplaceable': {'mRNA': 4},
        'other_sequence': {'exon': 1},
    }
    # parentless's, twoGenes's, phantomC's and phantomD's exons
    assert grouping_organizer.stats.dropped_lines_below == {
        'no_parent': {'exon': 1},
        'several_genes': {'exon': 1},
        'missing_parent_unplaceable': {'exon': 2},
    }


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
    assert (stats.gene_parented_cds_genes, stats.gene_parented_cds_duplicates) == (1, 1)
    gene_b = cds_without_transcript_controller.session.query(SuperLocus).filter_by(given_name='geneB').one()
    assert [t.given_name for t in gene_b.transcripts] == ['rnaB1']


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
    strand; it is left as it is where they are none or its frame is off. See
    testdata/stop_codon.gff3."""
    controller = ImportController(database_path='sqlite:///:memory:')
    controller.add_genome('testdata/stop_codon.fa', 'testdata/stop_codon.gff3', clean_gff=True)
    transcripts = transcripts_by_name(controller)

    expected = {'tPlus': (20, 53), 'tSplit': (121, 181), 'tMinus': (269, 236),
                'tNoStop': (320, 350), 'tFrame': (420, 449)}
    assert {name: feature_ranges(t, types.GEENUFF_CDS)[0] for name, t in transcripts.items()} == expected
    for name in ['tPlus', 'tSplit', 'tMinus']:
        assert not error_types(transcripts[name])
    assert types.MISSING_STOP_CODON in error_types(transcripts['tNoStop'])
    assert types.TRUNCATED_CDS in error_types(transcripts['tFrame'])
    assert controller.stats.stop_codons_recovered == 3


# section: types
def test_enum_non_inheritance():
    allknowngff = [x.name for x in list(types.AllKnownGFFFeatures)]
    # check that some random bits made it in to all
    assert 'region' in allknowngff
    assert 'exon' in allknowngff
    # check nothing is there twice
    assert len(set(allknowngff)) == len(allknowngff)


def test_enums_name_val_match():
    for x in types.AllKnownGFFFeatures:
        assert x.name == x.value


# section: helpers
def test_key_matching():
    # identical
    known = {'a', 'b', 'c'}
    mapper, is_forward = helpers.two_way_key_match(known, known)
    assert is_forward
    assert mapper('a') == 'a'
    assert isinstance(mapper, helpers.CheckMapper)
    with pytest.raises(KeyError):
        mapper('d')

    # subset, should behave as before
    mapper, is_forward = helpers.two_way_key_match(known, {'c'})
    assert is_forward
    assert mapper('a') == 'a'
    assert isinstance(mapper, helpers.CheckMapper)
    with pytest.raises(KeyError):
        mapper('d')

    # superset, should flip ordering
    mapper, is_forward = helpers.two_way_key_match({'c'}, known)
    assert not is_forward
    assert mapper('a') == 'a'
    assert isinstance(mapper, helpers.CheckMapper)
    with pytest.raises(KeyError):
        mapper('d')

    # other is abbreviated from known
    set1 = {'a.seq', 'b.seq', 'c.seq'}
    set2 = {'a', 'b', 'c'}
    mapper, is_forward = helpers.two_way_key_match(set1, set2)
    assert is_forward
    assert mapper('a') == 'a.seq'
    assert isinstance(mapper, helpers.DictMapper)
    with pytest.raises(KeyError):
        mapper('d')

    # cannot be safely differentiated
    set1 = {'ab.seq', 'ba.seq', 'c.seq', 'a.seq', 'b.seq'}
    set2 = {'a', 'b', 'c'}
    with pytest.raises(helpers.NonMatchableIDs):
        mapper, is_forward = helpers.two_way_key_match(set1, set2)
        print(mapper.key_vals)


def test_gff_to_seqids():
    x = helpers.get_seqids_from_gff('testdata/testerSl.gff3')
    assert x == {'NC_015438.2', 'NC_015439.2', 'NC_015440.2'}
