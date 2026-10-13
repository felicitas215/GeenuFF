"""End-to-end test of GFF3 -> GeenuFF database -> Helixer h5 through the two command line
scripts, import2geenuff.py and Helixer's geenuff2h5.py, run with their default settings.

The expectations come from docs/overlap_masking.md and docs/spec_vs_gff.md, not from the current
behavior of either side; an assertion that fails marks a place where code and concept have
diverged. testdata/helixer_e2e.gff3 describes each case beside the lines producing it.

Every coordinate below is GeenuFF's: 0-based, start inclusive, end exclusive, and runs along the
plus strand on both strands.

Helixer is not a dependency of GeenuFF, only the other way round, so the module is skipped where
geenuff2h5.py is absent and is marked 'helixer' so that a general test run can deselect it with
`pytest -m 'not helixer'`.
"""
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

h5py = pytest.importorskip('h5py', reason='reading the Helixer h5 needs h5py')
GEENUFF2H5 = shutil.which('geenuff2h5.py')  # not working on Windows, but GeenuFF/Helixer are mostly for Linux
if GEENUFF2H5 is None:
    pytest.skip('Helixer is not installed; it requires GeenuFF, not the other way round, so this '
                'test is an optional extra', allow_module_level=True)

pytestmark = pytest.mark.helixer

IMPORT2GEENUFF = Path(__file__).parents[2] / 'import2geenuff.py'
TESTDATA = Path(__file__).parents[1] / 'testdata'
LENGTH = 6000  # of both sequences

INTERGENIC, UTR, CDS, INTRON = range(4)
UNDEFINED = -1  # the stretch two exported genes both claim, where whichever is written second wins

# (sequence, strand): masked ranges, from overlap_masking.md and spec_vs_gff.md
EXPECTED_MASKED = {
    ('CHR', 'plus'): [
        (2000, 2400),  # case 2 'Crossing, one gene kept': geneCrossDropped beyond geneCrossKept
        (2600, 4000),  # case 3 'Nested, both masked whole': geneNestOuter and its flanks
        (4200, 4600),  # case 4 'Which gene is kept': geneBothCodeDropped beyond geneBothCodeKept
        (5799, 5899),  # case 5 no 3' UTR: toward the sequence end, 201 // 2 = 100 bp of the gap
    ],
    ('CHR', 'minus'): [],
    ('CHR2', 'plus'): [
        (200, 2000),  # case 6 'Refused, both masked whole': both genes and their outer flanks
        (2800, 3900),  # case 7 geneUtrDropped beyond geneUtrKept and its flank, 500 bp of 2600
    ],
    ('CHR2', 'minus'): [
        (3560, 4200),  # case 8 no 3' UTR: toward the minus strand's 3' end, 640 bp of the 4200 bp gap
    ],
}
# (sequence, strand): (start, end, label), intergenic everywhere else
EXPECTED_LABELS = {
    ('CHR', 'plus'): [
        (200, 300, UTR), (300, 699, CDS), (699, 800, UTR),  # case 1 geneClean
        (1200, 1300, UTR), (1300, 1900, CDS), (1900, 2000, UTR),  # case 2 geneCrossKept only
        (2800, 2900, UTR), (2900, 3050, CDS), (3050, 3100, INTRON),  # case 3 both genes
        (3100, 3400, UNDEFINED), (3400, 3450, INTRON), (3450, 3699, CDS), (3699, 3800, UTR),
        (4600, 4650, UTR), (4650, 5148, CDS), (5148, 5200, UTR),  # case 4 geneBothCodeKept only
        (5500, 5550, UTR), (5550, 5799, CDS),  # case 5 geneNoThreeUtr
    ],
    ('CHR', 'minus'): [],
    ('CHR2', 'plus'): [
        (400, 500, UTR), (500, 1000, CDS), (1000, 1200, UNDEFINED),  # case 6 both genes
        (1200, 1700, CDS), (1700, 1800, UTR),
        (2200, 2300, UTR), (2300, 2699, CDS), (2699, 2800, UTR),  # case 7 geneUtrKept only
    ],
    ('CHR2', 'minus'): [
        (4200, 4599, CDS), (4599, 4800, UTR),  # case 8 geneMinusNoThreeUtr
    ],
}


def run_script(*args):
    result = subprocess.run([sys.executable, *map(str, args)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def read_strands(h5_path):
    """Gives per (sequence, strand) the label per base (argmax of y) and the sample weights, both
    indexed by plus strand coordinate. The test sequences are shorter than the default chunk size,
    so each strand is one chunk: its first LENGTH bases, padding following, stored reversed on the
    minus strand."""
    with h5py.File(h5_path, 'r') as h5:
        rows = zip(h5['data/seqids'][:], h5['data/start_ends'][:], h5['data/y'][:],
                   h5['data/sample_weights'][:])
    strands = {}
    for seqid, (start, end), y, weights in rows:
        assert {start, end} == {0, LENGTH}, 'each strand is expected in one chunk'
        step = 1 if start < end else -1
        strands[(seqid.decode(), 'plus' if step == 1 else 'minus')] = (
            np.argmax(y[:LENGTH], axis=1)[::step], weights[:LENGTH][::step])
    return strands


def masked_ranges(weights):
    """The half open ranges where the sample weights are 0."""
    edges = np.diff(np.concatenate([[0], (weights == 0).astype(int), [0]]))
    return list(zip(np.nonzero(edges == 1)[0].tolist(), np.nonzero(edges == -1)[0].tolist()))


@pytest.fixture(scope='module')
def strands(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp('helixer_e2e')
    db_path, h5_path = tmp_path / 'e2e.sqlite3', tmp_path / 'e2e.h5'
    run_script(IMPORT2GEENUFF, '--fasta', TESTDATA / 'helixer_e2e.fa',
               '--gff3', TESTDATA / 'helixer_e2e.gff3', '--db-path', db_path,
               '--log-file', tmp_path / 'import.log', '--species', 'dummy')
    run_script(GEENUFF2H5, '--input-db-path', db_path, '--h5-output-path', h5_path)
    return read_strands(h5_path)


@pytest.mark.parametrize('key', EXPECTED_MASKED, ids='-'.join)
def test_each_strand_is_masked_where_documented(strands, key):
    """A mask covers only its own strand, so a strand without genes stays unmasked."""
    _, weights = strands[key]
    assert masked_ranges(weights) == EXPECTED_MASKED[key]


@pytest.mark.parametrize('key', EXPECTED_LABELS, ids='-'.join)
def test_each_strand_carries_the_labels_of_the_exported_genes(strands, key):
    """Dropped genes (geneCrossDropped, geneBothCodeDropped, geneUtrDropped) leave no label
    behind."""
    labels, _ = strands[key]
    expected = np.full(LENGTH, INTERGENIC)
    for start, end, label in EXPECTED_LABELS[key]:
        expected[start:end] = label
    defined = expected != UNDEFINED
    assert np.array_equal(labels[defined], expected[defined])
