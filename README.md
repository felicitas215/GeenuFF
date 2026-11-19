# GeenuFF

Schema and API for a relational db that encodes gene models in an explicit, structured,
and robust fashion. GeenuFF is currently mostly intended for
[Helixer](https://github.com/usadellab/Helixer). As such, since Version 0.4.0 it shifted
away from the previous aim of providing a way of unambiguously encoding gene models
(the way the DNA sequence is interpreted to produce proteins) that is both
robust to partial information and biological complexity. An example of partial information
still kept are trans-spliced genes (see the [trans-splicing documentation](docs/trans_splicing.md)).
An example of left out data is genes sharing an ID.

We're still keen to get feedback or ideas from the community (even if it's just whether
you think this could be useful to you if developed further), but if you build on GeenuFF
as it is now for anything other than getting training data for Helixer, you're doing so
at your own risk.

A more extensive description can be found [here](https://usadellab.github.io/GeenuFF/).

## What's new in 0.4.0

GeenuFF 0.4.0 is a major rework of the import, focused on reliable training data for Helixer. It
is released together with Helixer 0.3.8, which requires it.

> [!WARNING]
> Databases created before 0.4.0 fail to import into Helixer 0.3.8, and many GeenuFF scripts do
> not work with them either. **Please reimport** your genomes.

- **Masking:** a gene with an error in its structure (e.g. missing start/stop codon, in-frame stop
  codon) is masked in full instead of only its flank. A missing UTR masks only a flank of a
  reasonable size, so far less intergenic sequence is masked. A bug placing the 5' flank mask from
  the upstream gene's end was fixed. New checks: too short or truncated CDS, in-frame stop codons,
  overlapping CDS, CDS outside exons, CDS without transcript, transcripts past their sequence.
- **Overlapping genes:** of two overlapping coding genes, the better one is kept and its partner
  left out of the export, or both are masked (see [docs/overlap_masking.md](docs/overlap_masking.md)).
- **Stricter reading:** lines are grouped into genes by `ID`/`Parent` instead of line order; GFF3
  sequence IDs must be identical to the FASTA IDs; duplicate FASTA IDs or invalid sequence
  characters stop the import, as do GTF files and damaged gzipped or zipped files; unknown feature
  types are skipped. A stop codon left out of the CDS, as in GFF3 converted from GTF, is added.
- **No chromosome size limit:** chromosomes longer than 2,147,483,647 bp can be imported.
- **Faster:** imports are 1.7 to 3.7 times faster and need up to 4 times less memory (see
  [docs/performance.md](docs/performance.md)).
- **Import summary:** the import ends with a summary of what was left out or masked and why;
  `summarize_geenuff_db.py` recreates it from a database.
- **New scripts** to fix common GFF3/FASTA problems before import (`check_seqids.py`,
  `rename_gff_seqids.py`, `add_gene_lines.py`, `add_transcript_lines.py`, `rebuild_exons.py`,
  `filter_sequence_lengths.py`) and to read a database or compare against a reference
  (`dump_filtered_gff3.py`, `filter_longest_transcripts.py`, `list_transcripts.py`), all installed
  onto the PATH (see [docs/scripts.md](docs/scripts.md)).
- `import2geenuff.py`: `--db-path` must end in `.sqlite3`, and the log is written beside the
  database by default. Python 3.10.12 or newer is required; `dustdas` is no longer a dependency.

## Details on the relation to Helixer

GeenuFF is the annotation store [Helixer](https://github.com/usadellab/Helixer) trains from.
Importing a GFF3 checks the gene models, repairs what can be repaired unambiguously, and records
what stays ambiguous as errors of the transcript, merged into its geenuff_mask features. Helixer
reads one transcript per gene through `GeenuffExportController.genome_query`, the masks keeping
unreliable annotation out of training.

Only the standard genetic code is checked against: ATG = start codon and TAA, TAG and
TGA = stop codons. Genes with an alternative start codon, and species using one or two of
the stop codons as regular codons, are masked as erroneous for now; organelle genomes are
out of scope. The exported transcript of a gene is its longest coding one, whatever its
errors, so in this release a gene whose longest isoform is masked gives no labels even where
a shorter one is error-free.

The database holds the annotation as given; what a consumer is handed, and what is masked, is
decided on export. Worth reading before relying on the output:

* [docs/spec_vs_gff.md](docs/spec_vs_gff.md): how GeenuFF's features differ from a GFF3's, the
  error types, what each of them masks, and errors that mask nothing
* [docs/overlap_masking.md](docs/overlap_masking.md): what happens where two genes claim the same
  sequence, which of them is kept, and why
* [docs/trans_splicing.md](docs/trans_splicing.md): why a gene whose pieces sit on different strands
  is stored but never exported, and why nothing is masked in its place
* [docs/scripts.md](docs/scripts.md): what each script in `scripts/` reads out of a database
* [docs/performance.md](docs/performance.md): import speed, memory and masking of GeenuFF 0.4.0
  against 0.3.2

## Install

GeenuFF needs Python 3.10.12 or newer. Install it preferably in a
[virtual environment](https://docs.python-guide.org/dev/virtualenvs/).

```bash
git clone https://github.com/usadellab/GeenuFF.git
cd GeenuFF
pip install .  # pip install -e . for editable installs
cd ..
```

Alternatively, install directly from GitHub via pip:

```bash
pip install git+https://github.com/usadellab/GeenuFF.git
```

To run the tests:

```bash
cd GeenuFF/geenuff
py.test
cd ../..
```

`tests/test_helixer_h5.py` (inside `GeenuFF/geenuff`) runs the import and Helixer's `geenuff2h5.py`
end to end. It is skipped where Helixer is not installed, and `py.test -m 'not helixer'` leaves it
out where it is.

## Usage

You can run `bash example.sh` for a quick start with public data.
This will set up the folder 'three_algae', download public data in
the expected format, and import it into a GeenuFF sqlite database for each
species. For more information please see
[the api docs](https://usadellab.github.io/GeenuFF/api.html).

Each import ends with a summary that counts down from the GFF lines to the genes giving Helixer
labels, step by step: what was left out or masked on the way, and why, including the errors of the
exported transcripts. It is the quickest way to see whether an annotation is in the state you
expected. Where no gene, or fewer than one in ten of the coding transcripts read, gives labels, a
warning after it names the largest loss and the scripts in [docs/scripts.md](docs/scripts.md)
that prepare an incorrectly formatted GFF3 file for import.

## Future plans

* Add a validation module to check structure of gene models.
* Visualization.

## Thanks

To @janinamass for discussion and advice.