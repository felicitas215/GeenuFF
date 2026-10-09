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
*flank*. It is measured from the gene's selected transcript toward the closest edge of another
coding gene on that side, or toward the end of the sequence: a gene dropped for an overlap counts
as a neighbor, one overlapping the gene does not, and genes without a CDS never bound a flank, so
that a transposon or lncRNA overlapping a coding gene cannot leave it without one. The flank
reaches `min(gap // 2, int(sqrt(gap)) * 10)` bp into the gap: up to a gap of about 400 bp the
flanks of two genes nearly meet, beyond it the middle of the gap stays usable as intergenic.

**To watch: masking in dense genomes.** Up to a gap of about 400 bp an erroneous gene masks half of
every gap next to it, and a gap between two erroneous genes is masked entirely, save one base:

| gap    | flank per erroneous side | gap masked, one erroneous neighbor  | gap masked, both erroneous |
|--------|--------------------------|-------------------------------------|----------------------------|
| 100 bp | 50 bp                    | 50 %                                | 100 %                      |
| 400 bp | 200 bp                   | 50 %                                | 100 %                      |
| 1 kb   | 310 bp                   | 31 %                                | 62 %                       |
| 2.5 kb | 500 bp                   | 20 %                                | 40 %                       |
| 10 kb  | 1000 bp                  | 10 %                                | 20 %                       |

Within one gap this rarely masks clearly intergenic sequence, 50-200 bp being a typical UTR length.
The risk is a bias: in compact genomes (fungi, many algae, gene-dense plant regions) short
intergenic stretches are underrepresented in training, the unmasked intergenic sequence coming
mostly from long gaps or gaps beside clean genes. To check on real data: the share of intergenic bp
masked, by gap size, on a dense and a sparse genome; a possible remedy is a share below half the
gap.

| type                     | cause                                                | masked                                  |
|--------------------------|------------------------------------------------------|-----------------------------------------|
| missing_utr_5p           | CDS starts where the transcript starts               | 5' flank, up to the CDS start           |
| missing_utr_3p           | CDS ends where the transcript ends                   | 3' flank, from the CDS end              |
| missing_start_codon      | spliced CDS does not begin with ATG                  | whole gene and both flanks              |
| missing_stop_codon       | spliced CDS does not end with a stop codon           | whole gene and both flanks              |
| truncated_cds            | spliced CDS length is not a multiple of 3            | whole gene and both flanks              |
| inframe_stop_codon       | a stop codon in frame before the end of the CDS      | whole gene and both flanks              |
| truncated_intron         | the transcript line reaches past its outermost exon  | whole gene and both flanks              |
| too_short_intron         | full intron shorter than `min_intron_length` (20 bp) | whole gene and both flanks              |
| overlapping_exons        | two exon lines of one transcript overlap             | whole gene and both flanks              |
| overlapping_cds          | two CDS lines of one transcript overlap              | whole gene and both flanks              |
| floating_cds             | CDS lines without a transcript line (see grouping)   | whole gene and both flanks              |
| beyond_sequence_edge     | coding transcript reaching past its sequence         | whole gene, as clipped, and both flanks |
| wrong_starting_phase     | phase of the first CDS piece in the file is not 0    | nothing, only recorded                  |
| super_loci_overlap_error | two coding genes share sequence                      | see below and `overlap_masking.md`      |

A missing UTR leaves only the transcript's end unknown; the CDS stays labeled. Every other error
that masks means the CDS boundaries cannot be trusted, so the gene is masked whole with both
flanks: masking only part would leave a hole whose edges read as transitions that are not there
(see `overlap_masking.md`). A too short intron typically stands in for an assembly frameshift an
annotation pipeline bridged. An intron the transcript starts or ends in is a truncated intron,
its length not checked, and the transcript's start or end counts as not biological.

Overlapping exon or CDS lines cannot exist, so both ends of the transcript and CDS count as not
biological and nothing else is checked, the spliced CDS repeating the shared bases. A wrong
starting phase changes no label, the importer setting every phase itself, so it masks nothing; the
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
GTF can keep it that way. Where a CDS is whole (no inframe stop, truncation, ...), only the stop
codon is missing and the next 3 bases along the transcript's exons are a stop codon, the CDS is
extended by them, across an intron if need be: translation ends there whatever the file says. Only
these 3 bases are read; a stop codon further on would make a longer protein than the file
annotates. The start codon is not recovered, GTF counting it as part of the CDS, so a missing one
is likely a legit error. The import summary counts the CDS extended, so a file converted from GTF
shows as such.

##### grouping lines into genes

Lines are grouped into genes by their `ID` and `Parent` attributes, whatever their order in the
file. The import summary counts every line left out per reason and GFF type, the lines below it
under the same reason (e.g. the transcripts of a gene line with a shared ID under `shared_id`),
and states the coding transcripts, i.e. potential training data, each reason cost. CDS lines
naming a gene that has transcripts are counted per gene, as perhaps an isoform without a
transcript line.

The training data section accounts for every coding gene, of which Helixer gets one transcript:
exported and labeled in full, with or without sequence beside it masked, or masked in full, each
by cause; or not exported, by reason, with how many genes dropped for an overlap would otherwise
have been labeled. Errors are counted over the exported transcripts only, and counts of 0 are
shown only in that section.

- Gene and transcript lines need an `ID`, being the lines others name as their parent; one without
  an `ID`, or sharing it with another gene or transcript line, is left out. Exon, CDS and UTR lines
  need no `ID`, and CDS lines sharing one `ID` are normal. They need a `Parent` naming a transcript
  to be part of a gene model; an exon or UTR line without one is left out, a CDS line without one
  is masked (see below).
- UTR lines are the Sequence Ontology's UTR types, e.g. `five_prime_UTR` and `three_prime_UTR`;
  the lowercase `five_prime_utr` and `three_prime_utr` some tools write are not SO types and are
  ignored.
- A transcript names exactly one gene. One naming no parent, several genes or another transcript
  is left out. Transcripts naming a `Parent` that matches no line share a gene inferred for them,
  spanning them all. If they lie on different sequences or strands, they are all left out: which of
  them belong together cannot be told.
- An exon, CDS or UTR line is put under every transcript it names. One naming a gene that has
  transcripts is left out, whether or not it duplicates a line of one of that gene's transcripts:
  it could equally belong to an isoform the file gives no transcript line of its own. An exon or
  UTR line naming no parent, a gene without transcripts or a `Parent` that matches no line is left
  out as well.
- CDS lines without a transcript line are masked (`floating_cds`) rather than turned into gene
  models, whose CDS, exons and UTRs would have to be guessed from IDs every source writes
  differently; masked, their sequence is not taught as intergenic. This covers CDS lines naming no
  parent (grouped by a shared `ID`, else one by one), a `Parent` matching no line (grouped by it)
  or a gene without transcripts (grouped by the gene). Each group gets a transcript spanning it,
  under its gene or one created for it, masked whole with both flanks, and overlaps like any gene
  masked outright. A group is left out if its lines lie on different sequences or strands, or on
  another one than the gene they name. CDS lines naming a gene that has a transcript are left out
  unmasked, the gene's transcripts deciding the region.
- A transcript with CDS lines but no exon lines gets exons built from its CDS and
  UTR lines, lines that touch forming one exon. Overlapping lines stay apart, so that the overlap
  is found as an error.
- A line on another sequence than its parent is left out.

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

