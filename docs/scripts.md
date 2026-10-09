# Command line scripts

Five scripts in `scripts/` read a GeenuFF database and write something out of it. None of them
modifies the database. All take `--db-path-in` and write to stdout; all but
`summarize_geenuff_db.py` write to a file instead when given `-o`.

Importing a genome is a separate script, `import2geenuff.py` in the repository root, which is the
only one installed onto the PATH (see `pyproject.toml`). The rest are run from the checkout.

## dump_filtered_gff3.py

Writes a plain GFF3 of one transcript per gene, the longest coding one, for comparing a Helixer
prediction against the reference it was trained on.

| written                     | by default | `--include-erroneous` |
|-----------------------------|------------|-----------------------|
| labeled in full             | yes        | yes                   |
| masked in full              | no         | yes                   |
| dropped for an overlap      | no         | yes                   |
| features not placeable      | no         | no                    |

The default is what the h5 export gives labels from; a gene with only a missing UTR's flank or an
overlap partner's overhang masked beside it counts as labeled in full. `--include-erroneous` is
for comparing a prediction gene by gene rather than on sound genes only. The log says how many
genes were left out and why. An mRNA line names its transcript's errors, e.g.
`geenuff_errors=missing_utr_5p,truncated_cds`, and that of a gene dropped for an overlap the
reason, `geenuff_excluded=overlap_dropped`.

## dump_to_fasta.py

Writes sequence to FASTA, one record per range. `--mode` is required and picks which ranges:

| mode         | what it writes                                  |
|--------------|-------------------------------------------------|
| `mRNA`       | the spliced transcript, exons joined            |
| `pre-mRNA`   | the transcript as transcribed, introns included |
| `CDS`        | the spliced coding sequence                     |
| `exons`      | each exon separately                            |
| `introns`    | each intron separately                          |
| `UTR`        | the spliced untranslated regions                |
| `pre-UTR`    | the untranslated regions before splicing        |
| `intergenic` | what lies between the genes                     |

`--longest` restricts it to one transcript per gene, the same one the h5 export uses, rather than
every isoform.

## dump_lengthinfo.py

Takes the same `--mode` and `--longest` as `dump_to_fasta.py` but writes the lengths of those
ranges rather than their sequence, one per line. `--stats-only` writes summary statistics over
them instead of the lengths themselves. Useful for asking how long the introns or UTRs of an
annotation actually are before deciding what a pipeline should tolerate.

## dump_json.py

Writes one strand of one sequence as JSON, for looking at a gene model in full rather than in
aggregate. It takes `--species`, `--seqid` and `--strand`, optionally narrowed by `--start` and
`--end` (0-based, end exclusive, ascending on either strand), `--longest` as above, and
`--pretty` for indented output. By default, every isoform of a gene is written.

Unlike the others this is a reading tool rather than an export: it shows every gene, exported or
not, with every feature GeenuFF built for it, geenuff_mask features included. Genes and
transcripts carry flags to filter or sort by: whether they are exported and why not, the
transcript selected for export, its errors and their severity, its mask ranges and whether it is
masked in full (see `api.md` for the format).

## summarize_geenuff_db.py

Rebuilds the import summary from the database, in the same layout and wording, for when the import
log is lost. It needs a database with the `transcript_error` table.

Its counts are identical to the log's for the gene and transcript totals, the training data
section and the errors; the overlap outcomes are recovered by running the overlap sweep again over
the exported spans. What was only counted while reading the GFF3 cannot be recovered and is listed
at the end: lines left out while grouping and the coding transcripts lost through them, created
genes, exons built, stop codons added and dropped transcripts.
