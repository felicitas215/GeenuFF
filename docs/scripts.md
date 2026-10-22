# Command line scripts

Six scripts in `scripts/` read a GeenuFF database and write something out of it. None of them
modifies the database. All take `--db-path-in`. The four `dump_*` scripts write to the file given
with `-o`, which they require; `list_transcripts.py` writes to stdout or, given `-o`, to a file,
and `summarize_geenuff_db.py` to stdout.

Importing a genome is a separate script, `import2geenuff.py` in the repository root, and
`check_seqids.py`, `rename_gff_seqids.py`, `add_gene_lines.py`, `add_transcript_lines.py`,
`rebuild_exons.py` and `filter_sequence_lengths.py` (i.e. the preparation scripts), check and
prepare the files for import.  
Installing GeenuFF, directly or as a dependency of Helixer, puts all of them onto the PATH
(see `pyproject.toml`), so they are called by name, e.g.
`check_seqids.py --fasta genome.fa --gff genome.gff3`.

`dump_filtered_gff3.py`, `rename_gff_seqids.py`, `add_gene_lines.py`, `add_transcript_lines.py`,
`rebuild_exons.py` and `filter_sequence_lengths.py` log what they leave out or change; `--log-file`
writes that log to a file as well as the terminal, to keep it beside the output. No script writes
over its input files or database: an output path naming one stops it before anything is written,
as do an input file or database that does not exist and an output path in a directory that does
not exist. The scripts reading a FASTA or GFF3 file, like `import2geenuff.py`, also stop where the
two look swapped by their names, the FASTA file does not start with `>`, or the GFF3 file's first
feature line does not have 9 tab-separated columns (e.g. 10 from a sequence ID holding a tab).

The preparation scripts read gzipped or plain files alike. Each file they write is gzipped where its
input file is, its name given `.gz` or not to match; `--gzip` or `--no-gzip` can be used to override
this for all of them. Gzipped output is written at `--compression-level` 4 by default (1 fastest,
9 smallest), compressing a genome nearly as well as gzip's default 6, several times faster. It is
compressed in `--threads` 4 threads by default (1 for a single one, -1 for all the system has).
Zipped files holding a single file are read as well; the output of a zipped input file is gzipped
and named `.gz`.

Without an output path (`-o`, or `--fasta-out` and `--gff-out` for `filter_sequence_lengths.py`),
each preparation script writes beside its input, the input's name with a suffix before the
extension, so a chain of them names every step, e.g.
`genes_renamed_length_filtered_with_genes.gff3`:

| script                       | suffix              |
|------------------------------|---------------------|
| `rename_gff_seqids.py`       | `_renamed`          |
| `filter_sequence_lengths.py` | `_length_filtered`  |
| `add_gene_lines.py`          | `_with_genes`       |
| `add_transcript_lines.py`    | `_with_transcripts` |
| `rebuild_exons.py`           | `_exons_rebuilt`    |

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

| mode         | what it writes                                                                          |
|--------------|-----------------------------------------------------------------------------------------|
| `mRNA`       | the spliced transcript, exons joined                                                    |
| `pre-mRNA`   | the transcript as transcribed, introns included                                         |
| `CDS`        | the spliced coding sequence, needs `--coding` or `--longest`                            |
| `exons`      | each exon separately                                                                    |
| `introns`    | each intron separately                                                                  |
| `UTR`        | the spliced untranslated regions, needs `--coding` or `--longest`                       |
| `pre-UTR`    | the untranslated regions before splicing, needs `--coding` or `--longest`               |
| `intergenic` | what lies between the genes, named by GFF3 coordinates and strand, e.g. `Chr1:1-100(+)` |

| option      | transcripts written                                                                                                                          |
|-------------|----------------------------------------------------------------------------------------------------------------------------------------------|
| neither     | every isoform of genes typed `gene` with transcripts typed `mRNA` or `transcript`                                                            |
| `--longest` | one per gene, the coding one the h5 export uses, whatever the GFF3 types (e.g. `protein_coding_gene`); genes dropped for an overlap left out |
| `--coding`  | every transcript with a CDS, whatever its GFF3 types                                                                                         |

Genes whose strand or coordinates cannot be placed are always left out, as in the h5 export.
`intergenic` is what lies outside the transcripts of the chosen option, genes left out included.

## dump_lengthinfo.py

Takes the same `--mode`, `--longest` and `--coding` as `dump_to_fasta.py` but writes the lengths
of those ranges rather than their sequence, one per line. `--stats-only` writes summary
statistics over them instead of the lengths themselves. Useful for asking how long the introns or
UTRs of an annotation actually are before deciding what a pipeline should tolerate.

## dump_json.py

Writes one strand of one sequence as JSON, for looking at a gene model in full rather than in
aggregate. By default, every isoform of genes typed `gene` with transcripts typed `mRNA` or
`transcript` is written.

| option                                         | what it does                                                                                                                                                 |
|------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--species`, `--seqid`, `--strand` (`+` / `-`) | required: the strand of one sequence to write                                                                                                                |
| `--start`, `--end`                             | narrow it to a range, 0-based, end exclusive, ascending on either strand (default: the whole sequence; GFF3 example: `CDS  22  89  .  -` -> start=21 end=89) |
| `--longest`                                    | the selected coding transcript per gene only, as for `dump_to_fasta.py`                                                                                      |
| `--show-all`                                   | genes and transcripts of every GFF3 type, e.g. tRNA or lncRNA genes; for curiosity only, as Helixer's training data is coding only                           |
| `--pretty`                                     | indented output                                                                                                                                              |

Unlike the others this is a reading tool rather than an export: it shows every gene (of every type
with `--show-all`), exported or not, with every feature GeenuFF built for it, geenuff_mask
features included. Genes and transcripts carry flags to filter or sort by: whether they are
exported and why not, the transcript selected for export, its errors and their severity, its
mask ranges and whether it is masked in full (see `api.md` for the format).

## summarize_geenuff_db.py

Rebuilds the import summary from the database, in the same layout and wording, for when the import
log is lost. It needs a database with the `transcript_error` table.

From "genes imported" on, its counts are identical to the log's (each gene stores how many of its
transcripts were left out); the overlap outcomes are recovered by running the overlap sweep again
over the exported spans. What was only counted while reading the GFF3 cannot be recovered and is
listed at the end: the GFF lines and GFF gene lines steps, with the lines left out while grouping,
the coding transcripts lost through them and the genes created, and the exons built and stop
codons added.

## list_transcripts.py

Lists transcripts by what GeenuFF found and decided, as a table or as gene or transcript IDs to
filter other files with. By default, every coding transcript is listed; filters combine with AND.

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

Checks whether every sequence ID of a GFF3 file is a FASTA ID (the header up to the first
whitespace). A GFF3 ID contained in exactly one FASTA ID, e.g. `Chr1` in `Chr1.v2`, is paired with
it by GeenuFF's import; the warning lists these pairs to look at, as they are no exact match. Where
IDs cannot be paired so, GeenuFF's import and many other tools fail or leave out the genes on them;
the warning names them and points to `rename_gff_seqids.py` for sequences renamed after
annotation. It exits with 1 then, else 0, so it can gate a pipeline.

## rename_gff_seqids.py

Renames the sequence IDs of a GFF3 file to those of its FASTA file. Take the FASTA and GFF3 files
from the same source wherever possible; this script is for when it is not, e.g. a FASTA file
uploaded to NCBI, which renames and may trim sequences, beside the GFF3 file made before.

It also serves as a filter where the IDs match already: the lines of GFF3 sequence IDs that are in
no FASTA sequence, e.g. sequences removed from the FASTA file after annotation, are left out, up
to 15% of the coding transcripts lost (`--max-missing-percent`); above that, nothing is written.

Each GFF3 sequence ID is paired with the first FASTA sequence found by:

| by        | FASTA sequence                                                                                                         |
|-----------|------------------------------------------------------------------------------------------------------------------------|
| ID        | of the same ID (its header line up to the first whitespace)                                                            |
| full name | whose whole header line holds the GFF3 ID as a word, e.g. `scaffold_12`                                                |
| sequence  | the one the genes fit best, where they have enough checks to rule out chance; a tie is logged and goes to the first    |

A check is one position where the GFF3 file places a start codon (ATG), a stop codon (TAA, TAG,
TGA, also right after the CDS as in GFF3 converted from GTF) or an intron end (GT or GC, AG), so a
transcript with n introns gives 2 + 2n. Codons are read as GeenuFF's import reads them: along the
CDS lines and, for a stop codon right after the CDS, along the exons, across an intron where one
splits the codon. Up to 250 per sequence ID are used: those of its first and
last genes (20 each) and the rest spread evenly between them. A right sequence passes nearly all,
a wrong or shifted one about 9% by chance (a pair needs at least 50% in every window, 80% in all
when found by sequence alone). Identical checks, as from isoforms sharing a CDS start or an intron,
count once, so that one gene with many isoforms cannot decide a sequence alone.

A pair by ID or full name does not fit, its FASTA sequence changed after annotation, if:
- the genes reach past its end by more than its length; by less, they are read as crossing the
  origin of a circular sequence, as GFF3 writes such features (their end plus the sequence's
  length), and checked so, with a warning to check that the sequence is circular or was cut at its
  end after annotation (e.g. trimmed on upload to NCBI), the import masking genes past the end
- a `##sequence-region` length differs from its own by >= 2 bp
- or, with 20 checks or more, fewer than half pass in any window (the first 20,
  the last 20, about 25 each in between)

A cut or insertion shifts every check after it, so a window shows it wherever it lies
before the last gene. With 20 checks a wrong sequence passes half in about 3 of a million cases; a
right one whose annotation passes only 70% fails in about 2%. A `##sequence-region` 1 bp off, as
when the regions are written 0-based (starting at 0), is not used to judge the pair; a warning
asks to check the GFF3 file's coordinates.

| case                                                                                                | result                                                                                                                                                                                              |
|-----------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| all IDs match and fit                                                                               | checked only, no file written, not even `--pairs-out`                                                                                                                                               |
| sequences changed after annotation, up to 5% of the coding transcripts                              | these are left out of the GFF3 and FASTA files, listed in a warning; the FASTA copy without them is written beside the output GFF3 file, e.g. `genome_sequences_dropped.fna.gz` for `genome.fna.gz` |
| GFF3 IDs with no FASTA sequence, up to 15% of the coding transcripts lost (`--max-missing-percent`) | these lines are left out, counted in a warning                                                                                                                                                      |
| GFF3 IDs whose genes fit a FASTA sequence paired with another GFF3 ID                               | left out as duplicates, counted apart                                                                                                                                                               |
| most pairs by ID or full name not fitting                                                           | error: the files belong to different assemblies                                                                                                                                                     |
| coding transcripts changed above 5% or missing above 15%                                            | error, nothing written                                                                                                                                                                              |
| stop codons right after the CDS                                                                     | pass their check and are counted, those past the exons listed in a warning; GeenuFF adds them on import (see `spec_vs_gff.md`)                                                                      |

Each warning about left-out lines gives the share of coding transcripts lost. `--pairs-out` writes
every pair with its method, checks passed and FASTA full name. Lines on a kept sequence whose
parents all lie on a left-out one go with them and are counted in the log, as in
`filter_sequence_lengths.py` (see there for why).

A possible extension: where most pairs by ID or full name fail, search a sample of them by sequence
first, to pair all by sequence where the names were merely reassigned (e.g. scaffolds renumbered
after annotation) rather than stopping as for different assemblies.

## add_gene_lines.py

Adds gene lines for transcripts without a `Parent`, as gene predictors such as AUGUSTUS, BRAKER
and Tiberius write them, which stop GeenuFF's import where they are coding (a transcript naming a
`Parent` that matches no line gets a gene inferred instead). A transcript naming itself as
`Parent`, e.g. `ID=g7-mRNA-1;Parent=g7-mRNA-1`, is treated as having none:
GeenuFF leaves it out as naming a transcript as parent. Transcripts no CDS line names, e.g. the
tRNAs and rRNAs of NCBI organelle annotations, are left as they are: GeenuFF leaves them out either
way, as they give no training data. Each coding transcript's gene ID is taken from the first of
these sources that gives one:

| source                                                                                                                                | e.g. transcript line            | gene ID   |
|---------------------------------------------------------------------------------------------------------------------------------------|---------------------------------|-----------|
| an attribute, by default `geneID`, `gene_id` or `gene` (`--gene-id-attributes`)                                                       | `ID=g1.t1;geneID=g1` (Tiberius) | `g1`      |
| the ID up to the last occurrence of a split text, by default `.` (`--split-at`)                                                       | `ID=g653.t1` (AUGUSTUS)         | `g653`    |
| neither: a gene of its own (set `--one-gene-per-transcript` to apply this method to every transcript instead of trying the first two) | `ID=X1`                         | `X1_gene` |

Cutting at the last occurrence, `.` gives `Glyma.01G000100` for `Glyma.01G000100.1`; for
`gene7-mRNA-2` give `--split-at=-mRNA-` (with `=`, as the text starts with a dash). The log shows
a few example splits and warns where transcripts put into one gene do not overlap each other as
isoforms do, a sign of the wrong split text or attribute.

Transcripts of one gene ID get it as `Parent` and share a gene line spanning them, written before
the first of them. Where the file has a gene line of that ID already, whichever source gave the ID,
the transcripts only get it as `Parent` and no second gene line is written; a transcript naming
itself as `Parent` gets that gene in place of itself. A gene ID on several sequences or strands, or
already the ID of a line other than a gene line, gets the sequence and strand appended (e.g.
`g2_CHR_minus`) with a warning. Transcripts falling back to a gene of their own are warned about
too, unless `--one-gene-per-transcript` is set. Where every coding transcript names a `Parent`
other than itself, nothing needs adding and no file is written.

## add_transcript_lines.py

Puts exon, CDS and UTR lines that name the gene itself as `Parent` under a transcript of that gene.
GeenuFF masks CDS lines without a transcript line (`floating_cds`) rather than guess a gene model
from them, and leaves out those of a gene with transcripts, which then count as non-coding: such a
file gives little or no training data. Two cases are handled, the gene's lines getting the
transcript as `Parent` in place of the gene:

| case                                                         | e.g.                                                   | transcript                                                    |
|--------------------------------------------------------------|--------------------------------------------------------|---------------------------------------------------------------|
| the gene has no transcript line (some MAKER files)           | `CDS ID=g1:CDS:1;Parent=g1`                            | a new `mRNA ID=g1_mRNA;Parent=g1` spanning the gene, after it |
| the gene's only transcript has no exon, CDS or UTR lines yet | `mRNA ID=g1.1;Parent=g1`, `CDS ID=g1.1.CDS1;Parent=g1` | that transcript, `g1.1`                                       |

Genes with several transcripts, or one that has exon, CDS or UTR lines already, stay as they are;
so do genes with no CDS line naming them, e.g. NCBI pseudogenes with exon lines only, which are
non-coding and give no labels either way; the import summary counts the non-coding genes with CDS
lines naming them. The log warns where a gene has overlapping exon or CDS lines, a sign of several
isoforms written as one gene, which GeenuFF masks. It also counts the lines naming a `Parent` that
matches no line at all, which stay as they are, there being neither a gene nor a transcript line
to put them under. Where no line needs moving, no file is written.  
This should in general only be used on files where each gene parented feature belongs to one isoform
only, as multi-isoform features all parented to one gene can not be sorted into the correct isoform
reliably (automatically); which is why genes with several transcripts are left as they are.

## rebuild_exons.py

Rebuilds the exon lines of transcripts whose CDS lines reach outside them, for files listing exon
lines for the UTRs only (happens rarely). The new exons are the union of a transcript's
exon, CDS and UTR lines, touching lines forming one exon. Use it only where this is systematic,
i.e. GeenuFF's import warns about `cds_outside_exons` for more than 5% of the coding transcripts:
elsewhere the CDS or the exons are wrong, which GeenuFF masks and rebuilding would hide. The log
gives the share rebuilt and warns at 5% or below. Without such transcripts, no file is written.

## filter_sequence_lengths.py

Removes short junk sequences, e.g. organelle genomes or small unplaced scaffolds, from a FASTA
file and their lines from its GFF3 file before import:

| removed                                                 | how                                                                                                                                                                                                                                                                                                                                                                                                                     |
|---------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| a cluster of small sequences (off with `--no-adaptive`) | below the largest gap between the sorted log10 lengths of all sequences, where it stands out from the other gaps (robust z-score from their median and median absolute deviation, at least `MIN_GAP_RANK` = 3) and is at least a 2.5-fold step in length (`MIN_CLUSTER_STEP`); needs 5 sequences or more (`MIN_SEQUENCES_FOR_CLUSTERING`) and may hold at most 5% of the total length (`--max-cluster-length-fraction`) |
| sequences below the hard floor (`--min-seq-len`)        | always, by default below 1000 bp                                                                                                                                                                                                                                                                                                                                                                                        |

Nothing is removed where no FASTA sequence or no sequence with coding transcripts would remain.
Where the sequences below the floor hold over 5% of the bases or coding genes, they are still
removed, with a warning that a lower `--min-seq-len` keeps them.
In a fragmented scaffold-level assembly, the gaps in the dense middle of the length distribution
are tiny, so the largest gap, among the few longest scaffolds, stands out without separating junk
from the rest; the length share keeps such a cluster (e.g. 9458 of 9463 scaffolds, 97.6% of the
bp, in one moss assembly, *Ulota hutchinsiae*, NCBI GCA_050002435.1).  
In a chromosome-level assembly all gaps are small, so one barely larger stands out as well; the
2.5-fold step keeps the smallest chromosome there (e.g. 12.0 Mb beside 14.2 Mb in *Physcomitrium
patens* v7.1 from Phytozome, or 8.8 Mb beside 21.3 Mb in *Racopilum cuspidigerum*, NCBI
GCA_050083975.1), while organelle genomes and junk scaffolds lie 2.68-fold or more below the
shortest sequence kept. The log says where a gap stands out but is too small a step, and lists
what was removed; `--no-adaptive`, or not filtering at all, keeps it.  
Only the largest gap is cut at, so a short sequence among others of similar length stays, e.g. the
16.4 kb mitochondrion of the blue whale (*Balaenoptera musculus*, NCBI RefSeq GCF_009873245.2)
beside unplaced scaffolds of 18.5 kb and more. Sequences known by ID are easily removed with
[seqkit](https://bioinf.shenwei.me/seqkit/); `rename_gff_seqids.py` then leaves out their GFF3
lines, which the import would otherwise stop on:

```bash
# exclude.txt: the IDs of the sequences to remove, one per line
seqkit grep -v -f exclude.txt genome.fa.gz -o genome_excluded.fa.gz
rename_gff_seqids.py --fasta genome_excluded.fa.gz --gff genes.gff3.gz  # or use awk
```

| output                         | option          | default                                                                                            |
|--------------------------------|-----------------|----------------------------------------------------------------------------------------------------|
| filtered FASTA file            | `--fasta-out`   | the input name with `_length_filtered` before the extension, e.g. `genome_length_filtered.fna.gz`  |
| filtered GFF3 file             | `--gff-out`     | the input name with `_length_filtered` before the extension                                        |
| removed IDs with their lengths | `--removed-out` | not written                                                                                        |

The log gives the removed sequences, their length range and the coding genes and transcripts on
them, genes being the better measure of the training data lost as GeenuFF exports one transcript
per gene. The FASTA file is read in fixed-size chunks, so a sequence written on one line is never
held in memory whole.

Lines on a kept sequence whose parents all lie on removed ones, through any number of Parent steps,
are left out as well and counted in the log; a line naming only some of them keeps the others as
`Parent`. Such genes, spread over two sequences, are likely wrong, and a mask for their lines on
the kept sequence would sit where the file places part of a gene GeenuFF no longer sees. Leaving
them out matches the import, which leaves out a line on another sequence than its parent unmasked,
so a kept sequence imports the same with or without filtering. Left in, the CDS lines among them
would name a parent matching no line and be masked as `floating_cds`. `rename_gff_seqids.py` does
the same for the sequences it leaves out.

Without `--gff` only the FASTA file is filtered, e.g. before predicting with Helixer, which saves
compute otherwise spent on junk fragments.
