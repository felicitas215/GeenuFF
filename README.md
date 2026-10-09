# GeenuFF

Schema and API for a relational db that encodes gene models in an explicit, structured,
and robust fashion.

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
[virtual environment](https://docs.python-guide.org/dev/virtualenvs/), adding `-e` to `pip install`
for an editable installation:

```bash
git clone https://github.com/weberlab-hhu/GeenuFF.git
cd GeenuFF
pip install .
cd ..
```

To run the tests:

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

Each import ends with a summary: the GFF lines left out and what they cost, what Helixer gets from
each coding gene (labeled in full, masked, or not exported, and why), and the errors of the
exported transcripts. It is the quickest way to see whether an annotation is in the state you
expected.

## Major plans

* Add a validation module to check structure of gene models.
* Add extraction of raw & mature transcript, CDS, and protein sequence as a demo application.
* Visualization.

## Thanks

To @janinamass for discussion and advice.