# Command line tools

__Warning! little about this is stable or tested yet__

## import a species into the database

Imports a eukaryotic genome and annotation into a GeenuFF database, logging the import; `<>` marks
what you specify:

```
import2geenuff.py --fasta <PATH_TO_GENOME_FASTA_FILE> --gff3 <PATH_TO_GFF3_FILE> \
    --db-path <YOUR_DATA_NAME>.sqlite3 --log-file <YOUR_DATA_NAME>.import.log \
    --species <SPECIES_NAME>
```

Or with some of the testdata filled in:
```
geenuff_path=<PATH/TO/GeenuFF/>

import2geenuff.py --fasta $geenuff_path/geenuff/testdata/exporter.fa \
    --gff3 $geenuff_path/geenuff/testdata/exporter.gff3 \
    --db-path dummy.sqlite3 --log-file dummy.import.log --species dummy
```


With the files arranged as below, `--basedir` sets input and output in one go:

```
# target structure
<BASEDIR>/input/<YOUR_FILE>.fa
<BASEDIR>/input/<YOUR_FILE>.gff3

# simplified import
import2geenuff.py --basedir <BASEDIR> --species <SPECIES_NAME>

# the output files "<SPECIES_NAME>.sqlite3" and "import.log" will be written in
# the directory <BASEDIR>/output/
```

Explicit paths override those from `--basedir`. Import one species into one database once; to
redo it, delete the database or pass `--replace-db`.

## extract fasta sequences from a database

`scripts/dump_to_fasta.py` writes transcripts, CDS and other sequence breakdowns as FASTA (see
`docs/scripts.md` for all modes):

```
# to obtain the CDS (ignoring phase) you can run
python $geenuff_path/scripts/dump_to_fasta.py --db-path-in <GENUFF_DB> --mode CDS

# to obtain the final mRNA sequence you can run
python $geenuff_path/scripts/dump_to_fasta.py --db-path-in <GENUFF_DB> --mode mRNA

# to obtain the unspliced transcript you can run
python $geenuff_path/scripts/dump_to_fasta.py --db-path-in <GENUFF_DB> --mode pre-mRNA
```

## extract sequence lengths from a database

As above with `scripts/dump_lengthinfo.py`, or `--stats-only` for summary statistics.

# GeenuFF API

For more flexibility, use the python module directly. Nothing about it is stable yet.

## import
To change how and where a genome is imported (as the tests do), call the importer directly, e.g.
from `$geenuff_path/geenuff`:

```{python}
from geenuff.applications.importer import ImportController

controller = ImportController(database_path='sqlite:///' + EXPORTING_DB)
controller.add_genome('testdata/exporting.fa', 'testdata/exporting.gff3', clean_gff=True,
                      genome_args={'species': 'dummy'})
```

An in-memory database, `'sqlite:///:memory:'`, is handy for testing or exploring. A prokaryotic
importer and hooks for non-conforming GFFs are planned, but not there yet.

## sequence output
See `geenuff.applications.exporter` and `geenuff.applications.exporters.sequence`, e.g. use or
extend `FastaExportController` to do something with the sequence breakdowns other than writing
FASTA:

```{python}
from geenuff.applications.exporter import MODES
from geenuff.applications.exporters.sequence import FastaExportController

controller = FastaExportController(PATH_TO_GEENUFF_DB)
controller.prep_ranges(MODES['pre-mRNA'])
for export_group in controller.export_ranges:
    pre_mrna_seq = controller.get_seq(export_group)
    # your code here
    # check for your motif of interest, count kmers, send to custom output, etc...
```

For other breakdowns, see `RangeMaker` in `geenuff.applications.exporter`.

## gff3 output
look at `geenuff.applications.exporters.gff3` (class `FilteredGff3ExportController`)

Writes back a plain GFF3 file with just the longest CDS-containing (i.e. protein-coding)
transcript per super locus that is not masked whole, i.e. the same set old Helixer's h5
export gives labels from. Super loci with no coding transcript at all are skipped entirely.

```
python $geenuff_path/scripts/dump_filtered_gff3.py --db-path-in <GENUFF_DB> -o filtered.gff3
```

```{python}
from geenuff.applications.exporters.gff3 import FilteredGff3ExportController

controller = FilteredGff3ExportController(PATH_TO_GEENUFF_DB)
controller.write_filtered_gff3('filtered.gff3')
```

## json output

A simple query of one region, returning flattened JSON, e.g. for a database imported from
`geenuff/testdata/exporting.*` as above:

```
from geenuff.applications.exporters.json import JsonExportController

controller = JsonExportController(PATH_TO_GEENUFF_DB)
json_out = controller.coordinate_range_to_json(species='dummy',
                                               seqid='Chr1:195000-199000',
                                               start=1, end=3900, 
                                               is_plus_strand=True)
```
The range is ascending, `[start, end)`, on either strand. It returns every gene on that strand
whose gene model overlaps the range, whether or not it is exported, with all of its isoforms, or
with `longest` set on the controller only the transcript selected for export (the one Helixer
uses), in this format:

```
[{"coordinate_piece":
    {"id": int, "seqid": str, "sequence": str, "start": int, "end": int},
 "super_loci":
    [{"id": int,
      "given_name": str,
      "type": str,
      "is_fully_contained": bool,
      "overlaps": bool,
      "excluded_from_export": str or null (why the gene is left out of exports),
      "exported": bool (not excluded and with a transcript selected for export),
      "transcripts": [{"id": int,
                       "given_name": str,
                       "type": str,
                       "is_fully_contained": bool,
                       "overlaps": bool,
                       "selected_for_export": bool (the one transcript per gene an export uses),
                       "exported": bool (selected and its gene not excluded),
                       "errors": [str, ...] (the error types found, see spec_vs_gff.md),
                       "error_severity": "none", "flank_masked" or "masked_outright",
                       "masks": [[start, end], ...] (its geenuff_mask features),
                       "masked_in_full": bool (a mask covering the whole transcript),
                       "features": [{"id": int,
                                     "given_name": str,
                                     "protein_id": str,
                                     "type": str,
                                     "start": int,
                                     "start_is_biological_start": bool,
                                     "end": int,
                                     "end_is_biological_end": bool,
                                     "score": float,
                                     "source": str,
                                     "phase": int (in {0, 1, 2}),
                                     "is_plus_strand": bool,
                                     "is_fully_contained": bool,
                                     "overlaps": bool
                                    }, ...]
                      }, ...]
     }, ...]

}, ...]
```

"is_fully_contained" and "overlaps" are relative to the queried range: whether the gene model
lies wholly within it, or touches it at all, the geenuff_mask features beside it not counting.
With the default range, the whole sequence, every gene is fully contained. "error_severity" grades
the transcript's own errors as an overlap decision does (see overlap_masking.md), overlaps not
counting. A gene lies on one sequence, the importer leaving out lines on another sequence than
their parent, so one query covers it. Features in a transcript are always in 5'-3' order.
