import gzip
import zipfile

import pytest

from geenuff.parsing import fastahelper
from geenuff.parsing.fastahelper import FastaParser
from geenuff.parsing.gffhelper import GFFObject, GFFParser

FASTA = b'>s1 first sequence\nAC GT\nnn\n>s2\nA\tC\n'


def write_as(path, data: bytes, file_format: str) -> str:
    if file_format == 'gzip':
        path = f'{path}.gz'
        with gzip.open(path, 'wb') as f:
            f.write(data)
    elif file_format == 'zip':
        path = f'{path}.zip'
        with zipfile.ZipFile(path, 'w') as archive:
            archive.writestr('inner', data)
    else:
        with open(path, 'wb') as f:
            f.write(data)
    return path


@pytest.mark.parametrize('file_format', ['text', 'gzip', 'zip'])
def test_read_fasta_plain_gzipped_and_zipped(tmp_path, file_format):
    path = write_as(tmp_path / 'genome.fa', FASTA, file_format)
    assert list(FastaParser.read_fasta(path)) == [('s1 first sequence', 'ACGTnn'), ('s2', 'AC')]


def test_read_fasta_rejects_text_before_the_first_header(tmp_path):
    path = write_as(tmp_path / 'genome.fa', b'text\n' + FASTA, 'text')
    with pytest.raises(ValueError, match='does not start with ">"'):
        list(FastaParser.read_fasta(path))


def test_read_fasta_rejects_an_id_given_twice(tmp_path):
    path = write_as(tmp_path / 'genome.fa', FASTA + b'>s1 again\nACGT\n', 'text')
    with pytest.raises(ValueError, match='sequence ID "s1" more than once'):
        list(FastaParser.read_fasta(path))


@pytest.mark.parametrize('file_format', ['text', 'gzip', 'zip'])
def test_gff_feature_lines_plain_gzipped_and_zipped(tmp_path, file_format):
    gff = b'##gff-version 3\nchr1\t.\tgene\t1\t9\t.\t+\t.\tID=g1\n\n##FASTA\n>chr1\nACGT\n'
    path = write_as(tmp_path / 'genes.gff3', gff, file_format)
    assert list(GFFParser.feature_lines(path)) ==['chr1\t.\tgene\t1\t9\t.\t+\t.\tID=g1\n']


def test_fasta_pieces_rebuild_file_across_chunk_borders(tmp_path, monkeypatch):
    data = b'>s1 a long header line\nACGTACGTAC\nGTA\n>s2\nAC'
    path = write_as(tmp_path / 'genome.fa', data, 'text')
    for chunk_size in range(1, len(data) + 1):
        # monkeypatch: here the test needs to invoke functionality which depends on global settings
        # i.e. test multiple chunk sizes (usually one fixed size)
        monkeypatch.setattr(fastahelper, 'FASTA_CHUNK_SIZE', chunk_size)
        pieces = list(FastaParser.pieces(path))
        assert b''.join(piece for _, piece in pieces) == data
        assert [piece for is_header, piece in pieces if is_header] == [b'>s1 a long header line\n', b'>s2\n']


@pytest.mark.parametrize('column, expected', [
    ('ID=t1;Parent=g1', {'ID': 't1', 'Parent': 'g1'}),
    ('ID=t1; Parent=g1 ', {'ID': 't1', 'Parent': 'g1'}),  # spaces around the pairs
    ('ID=t1;', {'ID': 't1'}),  # trailing ';'
    # unescaped '=' in a value, technically illegal in gff but unlikely a real ID -> some will crash somewhere else
    ('ID=a=b', {'ID': 'a=b'}),
    ('Parent=t1;Parent=t2', {'Parent': 't1,t2'}),  # a tag given twice
])
def test_parse_attributes(column, expected):
    assert GFFObject.parse_attributes(column) == expected


def test_with_parents_replaces_every_parent_attribute():
    assert GFFObject.with_parents('ID=e1;Parent=t1; Parent=t2;Note=a;', ['g1', 'g2']) == 'ID=e1;Note=a;Parent=g1,g2'


def test_gff_object_parents_and_id():
    """Deduplication of parents only happens in get_Parent() not directly after parsing the GFF line."""
    line = GFFObject('chr1\t.\texon\t1\t9\t.\t+\t.\tParent=t1,,t1,t2;Parent=t2,t3\n')
    assert line.get_Parent() == ['t1', 't2', 't3']
    assert line.get_ID() is None
    assert GFFObject('chr1\t.\tgene\t1\t9\t.\t+\t.\tID=g1').get_Parent() is None
