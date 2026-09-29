# GeenuFF

Schema and API for a relational db that encodes gene models in an explicit, structured, and robust fashion.

## beta disclaimer

GeenuFF is currently *extremely* beta and very unstable.
We're keen to get feed back or ideas from the community
(even if it's just whether you think this could be useful to you
if developed further), but if you build on GeenuFF as it is now,
you're doing so at your own risk.

## Motivation

We developed this to provide a way of unambiguously encoding gene models,
(the way the DNA sequence is interpreted to produce proteins) that is both
robust to partial information and biological complexity.

A more extensive description can be found [here](https://weberlab-hhu.github.io/GeenuFF/).

## Relation to Helixer

GeenuFF is the annotation store [Helixer](https://github.com/usadellab/Helixer) trains from.
Importing a GFF3 file here does the interpreting: gene models are checked, repaired where that
can be done unambiguously, and whatever stays ambiguous is recorded as an error feature marking
the range it covers. Helixer then reads one transcript per gene through
`GeenuffExportController.genome_query` and turns it into its per-base-pair matrices, where those
error features become the mask that keeps unreliable annotation out of training.

That split is deliberate: the database holds the annotation as given, and nothing is deleted for
being wrong. Which parts of it a consumer is handed, and which are masked, is decided on the way
out. Two pieces of that are worth reading before relying on the output:

* [docs/spec_vs_gff.md](docs/spec_vs_gff.md): how GeenuFF's features differ from a GFF3's, the error types, and what a zero length error feature means
* [docs/overlap_masking.md](docs/overlap_masking.md): what happens where two genes claim the same sequence, which of them is kept, and why
* [docs/trans_splicing.md](docs/trans_splicing.md): why a gene whose pieces sit on different strands is stored but never exported, and why nothing is masked in its place
* [docs/scripts.md](docs/scripts.md): what each script in `scripts/` reads out of a database

## Install

GeenuFF needs python3.10.12 or newer, which `pyproject.toml` enforces via `requires-python`.

I would recommend installation in a virtual environment.
https://docs.python-guide.org/dev/virtualenvs/

From a directory of your choice (and preferably in a virtualenv):

Clone and install GeenuFF:

```bash
git clone https://github.com/weberlab-hhu/GeenuFF.git
cd GeenuFF
pip install .
cd ..
```

Add `-e` to install it editable, i.e. so that changes in the working tree take effect without
reinstalling.

And you might want to run the tests:

```bash
cd GeenuFF/geenuff
py.test
cd ../..
```

Alternatively, install directly from GitHub via pip:

```bash
pip install git+https://github.com/weberlab-hhu/GeenuFF.git
```

## usage

You can run `bash example.sh` for a quick start with public data.
This will set up the folder 'three_algae', download public data in
the expected format, and import it into a geenuff spec db for each
species. For more information please see
[the api docs](https://weberlab-hhu.github.io/GeenuFF/api.html).

Each import ends with a summary of what was found: how many genes and transcripts the file held,
how many survived to be exported, which genes overlap another and how that was settled, and every
error type with the number of transcripts it was found in. It is worth reading, as it is the
quickest way to see whether an annotation is in the state you expected.

## Major plans

* Add a validation module to check structure of gene models.
* Add extraction of raw & mature transcript, CDS, and protein sequence as a demo application.
* Visualization.

## Thanks

To @janinamass for discussion and advice.