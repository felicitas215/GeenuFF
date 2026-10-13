# Command line scripts

Six scripts in `scripts/` read a GeenuFF database and write something out of it. None of them
modifies the database. All take `--db-path-in` and write to stdout; all but
`summarize_geenuff_db.py` write to a file instead when given `-o`.

Importing a genome is a separate script, `import2geenuff.py` in the repository root, which is the
only one installed onto the PATH (see `pyproject.toml`). The rest are run from the checkout, as are
`check_seqids.py`, `rename_gff_seqids.py`, `add_gene_lines.py` and `filter_sequence_lengths.py`,
which check and prepare the files for import.

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

## list_transcripts.py

Lists transcripts by what GeenuFF found and decided, as a table or as gene or transcript IDs to
filter other files with. By default every coding transcript is listed; filters combine with AND.

| filter                                             | lists the transcripts                                                          |
|----------------------------------------------------|--------------------------------------------------------------------------------|
| `--error TYPE ...`                                 | with any of these error types                                                  |
| `--with-errors` / `--error-free`                   | with at least one error / checked and without any                              |
| `--severity none\|flank_masked\|masked_outright`   | by what their own errors take away                                             |
| `--masked-in-full` / `--labeled`                   | covered entirely by a mask / giving labels                                     |
| `--exported` / `--not-exported`                    | the one per gene an export takes / all others                                  |
| `--excluded-reason REASON ...`                     | of genes left out of exports for these reasons                                 |
| `--selected`, `--seqid ...`, `--include-noncoding` | the selected transcript per gene only, on these sequences, non-coding ones too |

`--output` chooses a table (default: IDs, sequence, GFF3 coordinates, export state, errors,
severity, masks), `gene_ids` or `transcript_ids` (one per line, e.g. for `grep -f`), or `counts`
(transcripts and genes per error type and severity). Non-coding transcripts and those of genes left
out before the error checks show "not checked" and match no error filter. E.g. the genes masked in
full although exported: `--masked-in-full --exported --output gene_ids`.

## check_seqids.py

Checks whether every sequence ID of a GFF3 file is a FASTA ID (the header up to the first space).
Where some are not, GeenuFF's import and many other tools fail or leave out the genes on them; the
warning names them and points to `rename_gff_seqids.py` for sequences renamed after annotation. It
exits with 1 then, 0 where the IDs match, so it can gate a pipeline.

## rename_gff_seqids.py

Renames the sequence IDs of a GFF3 file to those of its FASTA file. Take the FASTA and GFF3 files
from the same source wherever possible; this script is for when it is not, e.g. a FASTA file
uploaded to NCBI, which renames and may trim sequences, beside the GFF3 file made before.

Each GFF3 sequence ID is paired with the first FASTA sequence found by:

| by        | FASTA sequence                                                               |
|-----------|------------------------------------------------------------------------------|
| ID        | of the same ID (its header line up to the first space)                       |
| full name | whose whole header line holds the GFF3 ID as a word, e.g. `scaffold_12`      |
| sequence  | the only one the genes fit, where they have enough checks to rule out chance |

A check is one position where the GFF3 file places a start codon (ATG), a stop codon (TAA, TAG,
TGA, also right after the CDS as in GFF3 converted from GTF) or an intron end (GT or GC, AG), so a
transcript with n introns gives 2 + 2n. Up to 250 per sequence ID are used: those of its first and
last genes (20 each) and the rest spread evenly between them. A right sequence passes nearly all,
a wrong or shifted one about 9% by chance (a pair needs at least 50% in every window, 80% in all
when found by sequence alone).

A pair by ID or full name does not fit, the FASTA sequence having been changed after annotation,
if the genes reach past its end, a `##sequence-region` length differs from its own, or, with at
least 20 checks, fewer than half pass in any window of them (the first 20, the last 20, and about
25 each in between). A cut or insertion shifts every check after it, so the windows show it
wherever it lies before the last gene. With 20 checks a wrong sequence passes half in about 2 of a
million cases, and a right one with a poor annotation passing 70% fails in about 2%.

| case                                                                      | result                                                                                                                                               |
|---------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------|
| all IDs match and fit                                                     | checked only, no file written                                                                                                                        |
| sequences changed after annotation, up to 5% of the coding transcripts    | these are left out of the GFF3 and FASTA files, listed in a warning; the FASTA copy without them is written beside the output GFF3 file (`.fa(.gz)`) |
| GFF3 IDs with no FASTA sequence, up to 15% of the coding transcripts lost | these lines are left out, counted in a warning                                                                                                       |
| GFF3 IDs whose genes fit only a FASTA sequence paired already             | left out as duplicates, counted apart                                                                                                                |
| most pairs by ID or full name not fitting                                 | error: the files belong to different assemblies                                                                                                      |
| coding transcripts changed above 5% or missing above 15%                  | error                                                                                                                                                |
| stop codons right after the CDS                                           | counted (don't cause a check to not pass); GeenuFF adds them on import (see `spec_vs_gff.md`), those past the exons listed in a warning              |

Every warning gives the share of coding transcripts lost. `--pairs-out` writes every pair with its
method, checks passed and FASTA full name.

A possible extension: where most pairs by ID or full name fail, search a sample of them by sequence
first, to pair all by sequence where the names were merely reassigned (e.g. scaffolds renumbered
after annotation) rather than stopping as for different assemblies.

## add_gene_lines.py

Adds gene lines for transcripts without a `Parent`, as gene predictors such as AUGUSTUS, BRAKER
and Tiberius write them; GeenuFF imports nothing of a transcript without a gene (unless it
names a parent `Parent=gene1`, then the missing line is inferred). Each such transcript's gene
ID is taken from the first source that gives one:

| source                                                                                  | e.g. transcript line            | gene ID   |
|-----------------------------------------------------------------------------------------|---------------------------------|-----------|
| an attribute, by default `geneID`, `gene_id` or `gene` (`--gene-id-attributes`)         | `ID=g1.t1;geneID=g1` (Tiberius) | `g1`      |
| the ID up to the last occurrence of a split text, by default `.` (`--split-at`)         | `ID=g653.t1` (AUGUSTUS)         | `g653`    |
| neither: a gene of its own                                                              | `ID=X1`                         | `X1_gene` |

The split text is cut at its last occurrence, so `.` gives `Glyma.01G000100` for
`Glyma.01G000100.1`; for MAKER's `gene7-mRNA-2` give `--split-at -mRNA-`. The log shows a few
example splits and warns where transcripts put into one gene do not overlap each other, as isoforms
do, a sign of the wrong split text or attribute.

Transcripts of one gene ID share a gene line spanning them, written before the first of them, and
each gets it as `Parent`; a gene ID naming a gene line already is only set as `Parent`. A gene ID
on several sequences or strands, or already the ID of another line, gets the sequence and strand
appended (e.g. `g2_CHR_minus`), with a warning, as do transcripts given a gene of their own.
Without transcripts lacking a `Parent`, no file is written.

## filter_sequence_lengths.py

Removes short junk sequences, e.g. organelle genomes or small unplaced scaffolds, from a FASTA
file and their lines from its GFF3 file before import:

| removed                                                 | how                                                                                                                                                                        |
|---------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| a cluster of small sequences (off with `--no-adaptive`) | below the largest gap between the sorted log10 lengths of all sequences, where it stands out from the other gaps (robust z-score of at least 3); needs 5 sequences or more |
| sequences below the hard floor (`--min-seq-len`)        | always, by default below 1000 bp                                                                                                                                           |

Nothing is removed where no FASTA sequence or no sequence with coding transcripts would remain.
The output files default to the input files' names with `_length_filtered` before the extension,
e.g. `genome.fna.gz` gives `genome_length_filtered.fna.gz`; `--fasta-out` and `--gff-out` set
others. The log gives the removed sequences, their length range and the coding genes and
transcripts on them; genes are the better measure of the training data lost, as GeenuFF exporting
one transcript per gene. `--removed-out` lists every removed ID with its length. The FASTA file
is read in fixed-size chunks, so a sequence written on one line is never held in memory whole.

Side note on the clustering itself: The cluster is the part below the largest gap between the
sorted log10 lengths of all sequences, where that gap stands out from the others (a robust z-score
of at least MIN_GAP_RANK (3.0), from their median and median absolute deviation); with fewer than
MIN_SEQUENCES_FOR_CLUSTERING (5) sequences, or no such gap, only the floor applies. Nothing is
removed where no FASTA sequence or no sequence with coding transcripts would remain.
