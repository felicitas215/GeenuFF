## Specific changes vs gff
### features

#### types
Feature types differ from a GFF's; the "geenuff_" prefix marks where the meaning is not identical.

##### core types:
* geenuff_transcript: range of pre-mRNA / unspliced transcript
* geenuff_cds: range between the start and stop codon (ignoring introns)
* geenuff_intron: range between a donor and acceptor splice site
* geenuff_mask: range not to be trained on, everything the errors of its transcript mask merged
  into as few ranges as possible

##### error types:
Errors mark what is ambiguous or unknown about a gene model; only coding transcripts are checked.
Each error type found is recorded once per transcript in the `transcript_error` table, without a
range; what the errors mask is worked out once per transcript and stored as its geenuff_mask
features (see the table below). Annotation stored directly in GeenuFF rather than via GFF could
give errors more precise ranges.

An error leaving a gene's end unknown is extended into the unclaimed sequence beside the gene, the
*flank*, reaching `min(gap // 2, int(sqrt(gap)) * 10)` bp into the gap toward the closest edge of
another exported coding gene on that side and strand, or toward the end of the sequence. A gene
dropped for an overlap counts as a neighbor, one overlapping the gene does not. Genes without a
CDS never bound a flank, so that a transposon or lncRNA overlapping a coding gene cannot leave it
without one. Like every mask, the flank covers only the gene's own strand.

**To watch: masking in dense genomes.** Up to a gap of about 400 bp the flanks of two genes nearly
meet: an erroneous gene masks half of every gap next to it, and a gap between two erroneous genes
is masked entirely, save one base. Beyond it the middle of the gap stays usable as intergenic:

| gap    | flank per erroneous side | gap masked, one erroneous neighbor  | gap masked, both erroneous |
|--------|--------------------------|-------------------------------------|----------------------------|
| 100 bp | 50 bp                    | 50 %                                | 100 %                      |
| 400 bp | 200 bp                   | 50 %                                | 100 %                      |
| 1 kb   | 310 bp                   | 31 %                                | 62 %                       |
| 2.5 kb | 500 bp                   | 20 %                                | 40 %                       |
| 10 kb  | 1000 bp                  | 10 %                                | 20 %                       |

Within one gap this rarely masks clearly intergenic sequence, 50-200 bp being a typical UTR length.
The risk is a bias: in compact genomes (fungi, many algae, gene-dense plant regions) short
intergenic stretches are underrepresented in training, unmasked intergenic sequence coming mostly
from long gaps or gaps beside clean genes. To check on a dense and a sparse genome: the share of
intergenic bp masked by gap size; a possible remedy is a share below half the gap.

| type                     | cause                                                | masked                                  |
|--------------------------|------------------------------------------------------|-----------------------------------------|
| missing_utr_5p           | CDS starts where the transcript starts               | 5' flank, up to the CDS start           |
| missing_utr_3p           | CDS ends where the transcript ends                   | 3' flank, from the CDS end              |
| missing_start_codon      | spliced CDS does not begin with ATG                  | whole gene and both flanks              |
| missing_stop_codon       | spliced CDS does not end with a stop codon           | whole gene and both flanks              |
| truncated_cds            | spliced CDS length is not a multiple of 3            | whole gene and both flanks              |
| inframe_stop_codon       | a stop codon in frame before the end of the CDS      | whole gene and both flanks              |
| too_short_cds            | spliced CDS shorter than `min_cds_length` (60 bp)    | whole gene and both flanks              |
| truncated_intron         | the transcript line reaches past its outermost exon  | whole gene and both flanks              |
| too_short_intron         | full intron shorter than `min_intron_length` (20 bp) | whole gene and both flanks              |
| overlapping_exons        | two exon lines of one transcript overlap             | whole gene and both flanks              |
| overlapping_cds          | two CDS lines of one transcript overlap              | whole gene and both flanks              |
| cds_outside_exons        | a CDS line reaches outside the exon lines            | whole gene and both flanks              |
| floating_cds             | CDS lines without a transcript line (see grouping)   | whole gene and both flanks              |
| beyond_sequence_edge     | coding transcript reaching past its sequence         | whole gene, as clipped, and both flanks |
| wrong_starting_phase     | phase of the first CDS piece in the file is not 0    | nothing, only recorded                  |
| super_loci_overlap_error | two coding genes share sequence                      | see below and `overlap_masking.md`      |

A missing UTR leaves only the transcript's end unknown; the CDS stays labeled. Every other error
that masks means the CDS boundaries cannot be trusted, so the gene is masked whole with both
flanks: masking only part would leave a hole whose edges read as transitions that are not there
(see `overlap_masking.md`). A too short intron typically stands in for an assembly frameshift an
annotation pipeline bridged. An intron the transcript starts or ends in is a truncated intron,
its length not checked, and the transcript's start or end counts as not biological. A CDS below
60 bp, the stop codon included and counted after one left out of the CDS is added, is more likely a
fragment or an artifact than a gene; Helixer also drops predicted genes with less coding sequence
by default (HelixerPost's `--min-coding-length` of 60). Masked rather than left out, a real
micropeptide gene is not taught as intergenic either. Both minimums can be changed with a config
file given to `import2geenuff.py --config-file` (see `config/import.yml`).

Overlapping exon or CDS lines cannot exist, so both ends of the transcript and CDS count as not
biological and nothing else is checked, the spliced CDS repeating the shared bases. The same holds
for a CDS outside the exons, which lies in an intron. Above 5% of the coding transcripts with
either error, the import warns: for a CDS outside the exons, the file may list exon lines for the
UTRs only (see `rebuild_exons.py`); for overlapping exon lines, these may be written wrongly
throughout, e.g. nested, and where the CDS lines are intact, dropping the exon lines lets the
import build the exons from the CDS and UTR lines.  
A wrong starting phase changes no label, the importer setting every phase itself, so it masks nothing; the
phases of later CDS pieces are not checked, a CDS leaving its frame being caught as truncated_cds.

###### super_loci_overlap_error is positioned differently:
It marks sequence two exported coding genes both claim, every shared base carrying two
incompatible labels, rather than anything wrong with one gene, and is recorded for the exported
transcript only. Of two crossing genes the better one is kept whole and masked over what its
dropped partner covered beyond it, without a flank of its own; a nested or refused pair is masked
whole with both flanks. See `overlap_masking.md` for the rules and the `super_locus_overlap`
table, which records every overlapping pair of coding genes without masking.

###### errors masking nothing:
An error is recorded whether or not it masks anything: a wrong starting phase never does, and an
error extended into a flank masks nothing where there is no unclaimed sequence, e.g. the missing
5' UTR of a gene nested inside another. It still counts in the import statistics and is named in
the filtered GFF3 export's geenuff_errors attribute. A gene dropped for an overlap has its errors
recorded too, without a mask and outside the import statistics, so that why it was dropped in
favor of its partner can be looked up.

##### start_is_biological_start and end_is_biological_end:
When `True`, these attributes mean the start and end attributes
of a feature correspond to a meaningful biological transition.
* start: start of a region, inclusive
  * geenuff_transcript --> transcription start site (1st transcribed bp)
  * geenuff_cds --> start codon, the A of the ATG
  * geenuf_intron --> donor splice site, first bp of intron
* end: end of a region, exclusive (i.e. start of one there after)
  * geenuff_transcript --> 1 after transcription termination site (1st non-transcribed bp)
  * geenuf_cds --> 1 after stop codon, first non-coding bp, e.g. the N in TGAN
  * geenuff_intron --> 1 after acceptor splice site (1st bp that is part of final transcript)

When `False`, the start or end is not, or not known to be, a biological transition, though the
region is confidently of its type. E.g. where the first exon and the first CDS start at the same
position (+ strand), the 5' UTR is missing: the transcript's start_is_biological_start is False
and the range upstream of the CDS is masked, as it is unclear which part is UTR and which
intergenic.

#### feature start/end/at numbering

Features delineate a range counted from 0, start inclusive, end exclusive: "geenuff_cds, start" is
the A of the ATG, the first coding bp, and "geenuff_cds, end" the first non-coding bp after the
stop codon.

##### matching sequence IDs

Each GFF sequence ID is paired with the identical FASTA ID (the header up to its first
whitespace), else with the only FASTA ID containing it (`Chr1` with `Chr1.v2`), which the import
logs as a warning. A GFF ID matching no FASTA ID or several stops the import with an error listing
them. Two GFF IDs contained in the same FASTA ID are an error as well, and reading the GFF stops at
a `##FASTA` line or a line starting with `>`.

##### lines outside their sequence

Lines starting before position 1 or ending past the length of their sequence come from an
annotation of another assembly version, or from a gene crossing the origin of a circular molecule
written with an end past its length. What lies beyond cannot be read to check it or label it, and
a line starting before the sequence could not even be stored.

A coding transcript with such an mRNA, exon or CDS line is clipped to the part on the sequence, its
exon and CDS lines wholly beyond being left out, and masked whole with both flanks
(`beyond_sequence_edge`): where it ends is unknown, and its CDS cannot be checked. A transcript with
no CDS line on the sequence is left out entirely, and its gene is kept out of exports
(`excluded_from_export = 'outside_sequence'`), counted in the import summary. Trans-spliced genes,
written with a start past their end or on no definite strand, stay out of exports unmasked (see
docs/trans_splicing.md): whether and where the intron between their parts lies is unknown.

##### stop codons left out of the CDS

GTF writes the stop codon as a line of its own and leaves it out of the CDS, and GFF3 converted from
GTF can keep it that way. Where a CDS is whole (no inframe stop, truncation, ...) but for its stop
codon, and the next 3 bases along the transcript's exons are one, the CDS is extended by them,
across an intron if need be: translation ends there whatever the file says. Where the exons end
first (CDS lines only, or exons ending with the CDS), the 3 bases after the last exon are read and
the last exon and the transcript are extended to hold them. A stop codon further on is not looked
for, as it would make a longer protein than the file annotates. The start codon is not recovered,
GTF counting it as part of the CDS, so a missing one is likely a real error. The import summary
counts the CDS extended, so a file converted from GTF shows as such.

##### grouping lines into genes

Lines are grouped into genes by their `ID` and `Parent` attributes, whatever their order in the
file. The import summary counts every line left out per reason and GFF type, the lines below it
under the same reason (e.g. the transcripts of a gene line with a shared ID under `shared_id`),
and the coding transcripts, i.e. potential training data, each reason cost.

- Gene and transcript lines need an `ID`, being the lines others name as their parent; one without
  an `ID`, or sharing it with another gene or transcript line, is left out. Exon, CDS and UTR lines
  need no `ID`, and CDS lines sharing one `ID` are normal.
- Feature types are the Sequence Ontology's, plus three common non-SO spellings: NCBI's `lnc_RNA`
  as a transcript, and AUGUSTUS' lowercase `five_prime_utr` and `three_prime_utr` as UTRs, so that
  exons built from CDS and UTR lines keep the UTRs.
- A transcript names exactly one gene; one naming several genes or a transcript, even itself, is
  left out (`add_gene_lines.py` gives one naming itself its gene instead). Coding transcripts
  naming no parent at all, as gene predictors such as AUGUSTUS and Tiberius write them, stop the
  import with an error pointing to `add_gene_lines.py`; non-coding ones, e.g. organelle ncRNAs in
  NCBI annotations, are left out, giving no training data. Transcripts naming a `Parent` that
  matches no line share a gene inferred for them, spanning them all, unless they lie on different
  sequences or strands: then which belong together cannot be told, and all are left out.
- An exon, CDS or UTR line is put under every transcript it names. One naming a gene that has
  transcripts is left out unmasked, whether or not it duplicates a line of that gene's transcripts:
  it could equally belong to an isoform without a transcript line (the summary counts such CDS
  lines per gene, and the genes left non-coding by them; where the gene's only transcript has no
  lines of its own, `add_transcript_lines.py` puts them under it). An exon or UTR line naming no
  transcript otherwise is left out as well.
- CDS lines naming no transcript otherwise are masked (`floating_cds`) rather than turned into gene
  models, whose CDS, exons and UTRs would have to be guessed from IDs every source writes
  differently; masked, their sequence is not taught as intergenic. They are grouped by a shared
  `ID` (else one by one) where they name no parent, by the `Parent` where it matches no line, and
  by the gene where it has no transcripts (`add_transcript_lines.py` gives such genes an mRNA
  before import; only use for single isoform gff files!). Each group gets a transcript spanning
  it, under its gene or one created for it, masked whole with both flanks, and overlaps like any
  gene masked outright. A group is left out if its lines lie on different sequences or strands,
  or on another one than the gene they name.
- A transcript with CDS lines but no exon lines gets exons built from its CDS and UTR lines, lines
  that touch forming one exon. Overlapping lines stay apart, so that the overlap is found as an
  error.
- A line on another sequence than its parent is left out.

##### import summary

The summary is a funnel from the GFF lines down to Helixer's training labels:
1. GFF lines
2. GFF gene lines
3. Genes imported
4. Genes with transcripts
5. Coding genes
6. Exported coding genes (one transcript each)
7. Genes giving labels

Below each step stands what was subtracted on the way to the next one, and why, so the counts
add up; counts of 0 are left out, transcript counts are given in brackets. The **Exported coding
genes** list the errors masking genes in full, once per type and transcript, and their UTRs
each gene in one line: both annotated, only missing 5' or 3', or both missing. The **Genes giving
labels** are split by the UTR flanks masked beside them: none, one or both. The length of the UTR
mask can be zero for lack of room, e.g. gene starts/ends exactly at a sequence edge or lies
directly beside another exported gene. As the masks are stored merged, a 0-length flank mask can
be missed only in the very rare (probably impossible) case where a 1 bp gap beside the gene is
covered by the overhang of a gene dropped in its favor. These overhangs in general are counted
apart, being no part of the gene itself.

The summary ends with:
- the total number of overlapping masked (overlaps only counted on the same strand +/-)
- the bases masked as a share of both strands
- the bases masked as a share of the bases of the exported genes (both strands combined, for
  overlapping genes the union of their spans is used); missing UTRs have no known length and are
  not counted

The complete lists of the GFF types skipped are logged with `--verbose`. Where no gene, or fewer
than one in ten of the coding transcripts read, gives labels, a warning after the summary names
the largest loss and the scripts preparing a GFF3 file for import (see `docs/scripts.md`).

##### reverse complement

The coding start always points to the A of the ATG, whatever the strand. On the + strand the
sequence \[1, 4) can be sliced with 1 and 4 directly as python coordinates; on the minus strand it
cannot:

```
 0  1  2  3  4  5
.N [A .T .G )N .N
 |  |  |  |  |  |
 N. T. A. C. N. N.
```

For the reverse complement on the minus strand, the inclusive start is 3 and the exclusive end 0,
one off from the python coordinates:


```
 0  1  2  3  4  5
.N [A .T .G )N .N
 |  |  |  |  |  |
 N( T. A. C] N. N.
```

##### differences vs gff
How the features' coordinates compare to the closest GFF coordinate. First and last are as in a
coordinate-sorted GFF, so reversed on the - strand.

Plus strand (+)

| Common Name                         | GFF               | GFF start | GFF end | geenuff type       | bearing | position |
|-------------------------------------|:------------------|----------:|--------:|:-------------------|:--------|---------:|
| TSS, Transcription start site       | start 1st exon    |         x |         | geenuff_transcript | start   |    x - 1 |
| TTS, Transcription termination site | end last exon     |           |       x | geenuff_transcript | end     |        x |
| 1st bp of start codon               | start 1st CDS     |         x |         | geenuff_cds        | start   |    x - 1 |
| coding end                          | end last CDS      |           |       x | geenuff_cds        | end     |        x |
| donor splice site (5' of intron)    | end non-last exon |           |       x | geenuff_intron     | start   |        x |
| acceptor splice site (3' of intron) | start 2nd+ exon   |         x |         | geenuff_intron     | end     |    x - 1 |


Minus strand (-)

| Common Name                         | GFF               | GFF start | GFF end | genuff type        | bearing | position |
|-------------------------------------|:------------------|----------:|--------:|:-------------------|:--------|---------:|
| TSS, Transcription start site       | end last exon     |           |       x | geenuff_transcript | start   |    x - 1 |
| TTS, Transcription termination site | start 1st exon    |         x |         | geenuff_transcript | end     |    x - 2 |
| 1st bp of start codon               | end last CDS      |           |       x | geenuff_cds        | start   |    x - 1 |
| coding end                          | start 1st CDS     |         x |         | geenuff_cds        | end     |    x - 2 |
| donor splice site (5' of intron)    | start 2nd+ exon   |         x |         | geenuff_intron     | start   |    x - 2 |
| acceptor splice site (3' of intron) | end non-last exon |           |       x | geenuff_intron     | end     |    x - 1 |

