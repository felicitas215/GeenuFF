## Specific changes vs gff
### features

#### types
Feature types differ somewhat from a gff, "geenuff_" has been 
appended to names to reduce confusion where the meaning is not
identical.

##### core types:
* geenuff_transcript: range of pre-mRNA / unspliced transcript
* geenuff_cds: range between the start and stop codon (ignoring introns)
* geenuff_intron / trans intron: range between a donor and acceptor splice site

##### error types:
These do not exist in a gff, but are used in geenuff to denote things
that might be ambiguous or unknown about a gene model. 

Currently, errors are assigned when any obvious gene model inconsistency
is encountered during gff parsing. Only coding transcripts are checked. If gff format were not
used as an intermediary and gene annotation was performed and stored directly in a geenuff
structured database, all the errors could be assigned more precise ranges for any ambiguity.

An error that leaves a gene's end unknown is extended from that end into the gap of unclaimed
sequence beside the gene, the *flank*. The flank is measured from the gene with the error:
`min(gap // 2, int(sqrt(gap)) * 10)` bp into the gap toward the next exported coding gene that
way, or toward the end of the sequence where there is none. Two genes both extending into the
gap between them therefore never meet in the middle of a large gap, and sequence far enough from
any gene stays usable as intergenic. Genes without a CDS do not bound a flank, as a transposon or
lncRNA overlapping a coding gene would otherwise leave that gene no flank at all.

| type                     | cause                                                   | masked                               |
|--------------------------|---------------------------------------------------------|--------------------------------------|
| missing_utr_5p           | CDS starts where the transcript starts                  | 5' flank, up to the CDS start        |
| missing_utr_3p           | CDS ends where the transcript ends                      | 3' flank, from the CDS end           |
| missing_start_codon      | spliced CDS does not begin with ATG                     | whole gene and both flanks           |
| missing_stop_codon       | spliced CDS does not end with a stop codon              | whole gene and both flanks           |
| truncated_cds            | spliced CDS length is not a multiple of 3               | whole gene and both flanks           |
| inframe_stop_codon       | a stop codon in frame before the end of the CDS         | whole gene and both flanks           |
| truncated_intron         | the transcript line reaches past its outermost exon     | whole gene and both flanks           |
| too_short_intron         | intron shorter than `min_intron_length` (default 20 bp) | whole gene and both flanks           |
| overlapping_exons        | two exons of one transcript overlap                     | the transcript                       |
| wrong_starting_phase     | phase of the first CDS piece in the file is not 0       | nothing, only recorded               |
| super_loci_overlap_error | two coding genes share sequence                         | see below and `overlap_masking.md`   |

A missing UTR leaves only where the transcript ends unknown; the CDS itself is sound and stays
labelled. Every other error masking a flank means the CDS boundaries cannot be trusted: a missing
start or stop codon, a reading frame that is wrong (truncated CDS, in-frame stop codon), a gene
that is partial (truncated intron) or an intron too short to be spliced, which typically stands in
for a frameshift in the assembly that an annotation pipeline bridged and may or may not be right.
Such a gene is masked whole together with the flank on both sides. Masking only part of it would
leave a hole, whose edges read as transitions that are not there (see `overlap_masking.md`).

Overlapping exons are a structure that cannot exist, so none of the labels inside the transcript
can be trusted; for now only the transcript is masked. A wrong starting phase changes no label,
as the importer sets every CDS phase itself whatever the file says, so it is recorded as a zero
length error feature (see below). The file's phases of the other CDS pieces are not checked at
all: they are never used, and a CDS that really leaves its frame is caught as truncated_cds.

###### super_loci_overlap_error is positioned differently:
Unlike the types above, it is not extended into a flank of its own, and it marks
sequence two genes both claim rather than anything wrong with one gene. Every shared base
carries two mutually exclusive true labels (e.g. CDS for one gene and intron or UTR for the
other), which a one-class-per-base consumer cannot represent.

Of two overlapping genes the better one is exported whole and its partner is left out of exports
altogether; the error then covers what that partner occupied beyond the kept gene, or where it
sat if it lay wholly inside it. Where both are masked outright for their own errors, both genes
are exported and each is covered over its whole length. Only genes with a CDS
take part, so a coding gene annotated inside a transposable-element or other non-coding record
is not masked for it. Every overlapping pair of coding genes is recorded without masking in the
`super_locus_overlap` table as well. See `overlap_masking.md` for both records, the rules
deciding which gene is kept, and how this has changed.

###### zero length error features:
An error feature whose start equals its end masks nothing (consumers apply an error over
`[start, end)`), but still records that the error was found. This happens where an error type
that would normally be extended into a flank has no unclaimed sequence to
extend into, e.g. the missing 5' UTR of a gene nested inside another gene: the neighbouring
locus already reaches past the boundary the mask would start from, so the range available for
it is empty. The finding is real and the correct amount to mask is zero, so the feature is
kept rather than dropped, which also keeps it visible to consumers that treat any error
feature as disqualifying (e.g. the filtered GFF3 export). The import statistics count error
types as they are detected, so they agree with these features rather than with the subset
that masks a non-empty range.

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

When `False`, these attributes mean the start and end attributes
of a feature either do not, or it is not known if they correspond
to a biological transition, yet the region they delineate is
confidently of the given type. 

For instance, if the parser finds a gene model in a gff where the
start of the first exon and the start of the first CDS (+ strand)
have the same position (the A in ATG), then it is apparent that
we are missing the 5' UTR, so for the geenuff_transcript feature
the start_is_biological_start will be set to False, and an
error mask will be added upstream of the CDS. We are still confident
that all of the CDS must occur within the transcript, we know
the start codon is part of the transcript region, but we mark that the start point
itself is probably wrong, and mask the upstream range as it's 
unclear what part of this is intergenic and which part UTR.
 
#### feature start/end/at numbering

Features have start and end coordinates that
delineate a range. 

The positioning of these features is in keeping with the common
coordinate system: count from 0, start inclusive, end exclusive. 
So, the "geenuff_cds, start", is at the A, of the ATG, AKA the first
coding base pair; while in contrast, the "geenuff_cds, end" is
after the stop-codon, AKA, the first non-coding bp.

##### reverse complement

Importantly, the coding-start should always point to the first
A, of ATG, regardless of strand. This means the numeric coordinates
have to change and unfortunately while one could take the
sequence \[1, 4) on the + strand, and directly use 1 and 4 as python coordinates
and get the sequence; the same is not going to work on the minus strand.
Instead: 

```
 0  1  2  3  4  5
.N [A .T .G )N .N
 |  |  |  |  |  |
 N. T. A. C. N. N.
```

To get the reverse complement of this on the minus strand, we set the
inclusive start to 3, and exclusive end to 0. Note this is now off by 
one from the python coordinates


```
 0  1  2  3  4  5
.N [A .T .G )N .N
 |  |  |  |  |  |
 N( T. A. C] N. N.
```

##### differences vs gff
Cheat sheet for how the Features compare to the gff (in particular any discrepancy
between the closest coordinate in the gff, and the now standardized, consistent coordinate).

First and last for gff are reported as they are typically in gff (coordinate sorted),
so reverse to the interpretation when on the - strand. 

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

