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

## Details on the relation to Helixer

GeenuFF is the annotation store [Helixer](https://github.com/usadellab/Helixer) trains from.
Importing a GFF3 checks the gene models, repairs what can be repaired unambiguously, and records
what stays ambiguous as errors of the transcript, merged into its geenuff_mask features. Helixer
reads one transcript per gene through `GeenuffExportController.genome_query`, the masks keeping
unreliable annotation out of training.

The database holds the annotation as given; what a consumer is handed, and what is masked, is
decided on export. Worth reading before relying on the output:

* [docs/spec_vs_gff.md](docs/spec_vs_gff.md): how GeenuFF's features differ from a GFF3's, the
  error types, what each of them masks, and errors that mask nothing
* [docs/overlap_masking.md](docs/overlap_masking.md): what happens where two genes claim the same
  sequence, which of them is kept, and why
* [docs/trans_splicing.md](docs/trans_splicing.md): why a gene whose pieces sit on different strands
  is stored but never exported, and why nothing is masked in its place
* [docs/scripts.md](docs/scripts.md): what each script in `scripts/` reads out of a database

## Install

GeenuFF needs Python 3.10.12 or newer. Install it preferably in a
[virtual environment](https://docs.python-guide.org/dev/virtualenvs/), using `pip install -e .`
for an editable installation:

```bash
git clone https://github.com/usadellab/GeenuFF.git
cd GeenuFF
pip install .
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

`tests/test_helixer_h5.py` runs the import and Helixer's `geenuff2h5.py` end to end. It is
skipped where Helixer is not installed, and `py.test -m 'not helixer'` leaves it out where it is.

## Usage

You can run `bash example.sh` for a quick start with public data.
This will set up the folder 'three_algae', download public data in
the expected format, and import it into a GeenuFF sqlite database for each
species. For more information please see
[the api docs](https://usadellab.github.io/GeenuFF/api.html).

Each import ends with a summary that counts down from the GFF lines to the genes giving Helixer
labels, step by step: what was left out or masked on the way, and why, including the errors of the
exported transcripts. It is the quickest way to see whether an annotation is in the state you
expected.

## Future plans

* Add a validation module to check structure of gene models.
* Visualization.

## Thanks

To @janinamass for discussion and advice.