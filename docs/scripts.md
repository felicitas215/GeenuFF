# Command line scripts

Five scripts in `scripts/` read a GeenuFF database and write something out of it. None of them
modifies the database. All take `--db-path-in` and write to stdout unless given `-o`.

Importing a genome is a separate script, `import2geenuff.py` in the repository root, which is the
only one installed onto the PATH (see `pyproject.toml`). The rest are run from the checkout.

## dump_filtered_gff3.py

Writes a plain GFF3 of one transcript per gene, the longest coding one, for comparing a Helixer
prediction against the reference it was trained on.

By default, it writes exactly what the h5 export trains on: only genes that reach that export, and
of those only the ones with no error recorded at all.

`--include-erroneous` writes every gene that can be written at all, whatever is wrong with it,
which is the set to compare against when the question is what Helixer predicted per gene rather
than how it did on sound ones. Genes dropped from the h5 export for overlapping another are
included: nothing is wrong with them beyond sharing sequence, which a GFF3 holds without trouble.
Only genes whose features cannot be placed are left out either way, and the log says how many and
why (see `overlap_masking.md` and `trans_splicing.md`). The mRNA line of an erroneous transcript
names its error types, e.g. `geenuff_errors=missing_utr_5p,truncated_cds`, and that of a gene
dropped from the h5 export for an overlap the reason, `geenuff_excluded=overlap_dropped`.

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

Writes one region of one sequence as JSON, for looking at a gene model in full rather than in
aggregate. It takes `--species`, `--seqid` and `--strand`, optionally narrowed by `--start` and
`--end` (0-based, end exclusive, as everywhere in GeenuFF), `--longest` as above, and `--pretty`
for indented output.

Unlike the others this is a reading tool rather than an export: it is the quickest way to see
every feature GeenuFF built for a locus, geenuff_mask features included, and every transcript's
error types, when a summary count says something is wrong but not what.

## summarize_geenuff_db.py

Rebuilds the import summary from the database, for when the import log has been lost or the
database arrived without one. It reads through plain sqlite rather than the ORM, so a database
written by an older GeenuFF opens as well; sections resting on tables or columns that version did
not have are simply left out.

Most of the summary is recoverable, and the counts it prints are identical to the log's. Gene and
transcript totals, the empty and the reused-ID genes, the transcripts selected for export and how
many of those are error free, the errors per type and the contents of the `super_locus_overlap`
table are all read straight out. The overlap resolution is not stored as counts, but the outcome
is: the sweep is run again over the spans in the database, and the pairs split into resolved,
refused and chained by which genes carry `excluded_from_export = 'overlap_dropped'`.

What the summary counted while reading the GFF3 is gone for good, because the lines it counted
were never written: genes on no definite strand, transcripts dropped for having no storable range,
the exon and CDS lines parented to a gene instead of a transcript, and the lines skipped for an
unused feature type. The script lists these at the end.
