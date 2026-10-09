# GeenuFF

## Source & Install
See [github repository](https://github.com/weberlab-hhu/GeenuFF)

## What
Relational db Schema & Api to store and interpret gene structure

### Conceptual summary
GeenuFF defines ranges on the genomic sequence of a gene-structure "type" (transcribed, coding,
...), each delineated by a feature with a type, start and end. Which features make up a processed
molecule (mRNA, protein) is recorded in links to the 'outer' tables ("transcript", "protein").

GeenuFF encodes even partial information unambiguously: "start_is_biological_start" and
"end_is_biological_end" say whether a feature covers its full biological range or only part of
it, so a gene model states whether it is complete and exactly _where_ our knowledge of it ends.

It encodes both biological complexity and technical limitations: transcript pieces group the
features of a transcript originating from several genomic loci, whether split across two
scaffolds of a fragmented assembly or truly derived from two or more loci by trans-splicing. The
schema allows this; the importer, built for eukaryotic GFF3 and Helixer, writes one transcript
piece per transcript on one sequence, and keeps trans-spliced genes flattened and out of exports
(see [trans_splicing.md](trans_splicing.md)).

## Why
### general goal
A gene annotation (which parts of a genome are transcribed into RNA, how the RNA is processed into
mature mRNA, and which part of it is translated into protein) needs a tailored data structure to
recreate at will the coordinates of a gene, its parts, and the sequences made from them. It has
to handle cases like:
* several ways to 'put together' the transcripts or proteins of one locus
* several transcripts producing the same protein, differing only in mRNA
* one transcript (in prokaryotes) translated into several proteins (Note: currently GeenuFF is
  mostly tailored towards eukaryotes)
* rarely, one mRNA derived from several genomic loci by trans-splicing
* all of the above with partial information:
  * incomplete genomic sequence
  * fragmented genomic sequence
  * incomplete information about the gene annotation itself

### specific advantages over alternatives
GeenuFF is not the first attempt to encode gene structures; if you know e.g. GFF and wonder
whether this is one more standard on the [pile](https://xkcd.com/927/), this section explains
_why_ we think it is worth it.

Today GeenuFF is mainly the store of eukaryotic annotation that
[Helixer](https://github.com/usadellab/Helixer) trains from: it imports GFF3, checks and repairs
the gene models, and records what stays ambiguous as errors and masks. The schema below is more
general than what the importer and the exports build from it; where they fall short of it, this
section says so.

Gene annotations are usually encoded in a variant of the GFF format (GTF, GFF, GFF3), see
[http://gmod.org/wiki/GFF3](http://gmod.org/wiki/GFF3). Superficial drawbacks, such as the custom
key-value encoding of relationships and meta info in the last column, are solved by a good parser
or tools such as [gffutils](https://daler.github.io/gffutils/). This section is about the more
fundamental issues of the underlying _structure_.

#### Gff-like implicit encodings
In GFF-like formats genes and their pieces are ranges, e.g. `[ exon ]`, which leaves many
components implicit. An intron is the gap between two exons:
```
# gff features
 [      transcript        ]
 [ exon ]          [ exon ]
# interpretation
 [ exon ]( intron )[ exon ]
```
The transcription start site (TSS) is implied as well:
```
# at the start of the 1st exon (+ strand)
 [    transcript    ]
 [ exon ]    [ exon ]
 ^
 TSS

# or, at the end of the last exon (- strand)
 [    transcript    ]
 [ exon ]    [ exon ]
                    ^
                    TSS
```
Beyond the extra parsing, the real problem is partial information. Say homology and truncated
RNAseq mappings show that at least the first exon of a gene is missing. A GFF-like file can either
include the known exons, wrongly _implying_ the TSS lies where there is actually an acceptor
splice site, or skip the gene, wrongly _implying_ the region is intergenic. There is no right way
to record what is known and what is not.

##### GeenuFF more explicit encodings
GeenuFF also encodes ranges, with some _implicit_ encodings left to avoid redundancy, but the
elements encoded reflect biology, and a feature's start and end mean the same for every type. The
two-exon, + strand transcript from above becomes:
```
# gff-like
 [      transcript        ]
 [ exon ]          [ exon ]

# geenuff
 [       transcript        )
         [ intron  )
 ^
 start, TSS
```
Or on the - strand:
```
# geenuff
(       transcript        ]
        ( intron  ]
                          ^
                          start, TSS
```
No separate rule is needed to find the TSS on the minus strand. More importantly, start and end
carry a boolean (start_is_biological_start, end_is_biological_end) saying whether the transition
is a biological one (`True`) or the edge of a partially known gene model (`False`). The gene model
missing its first exon becomes:
```
# geenuff
                       [     transcript      )
                       ^                     ^
                       start                 end

start_is_biological_start=False
end_is_biological_end=True

# to mark the region before as unknown too (an intron, or an assembly error?), the error is
# recorded for the transcript (e.g. missing_utr_5p) and a geenuff_mask feature covers the region
[    geenuff_mask      )
^                      ^
start                  end
```

#### Gff-like misassignment of attributes and relations
A GFF-like gene model maps biology onto an ill-fitting model.

In GFF a transcript is always the child of a gene. That works for most eukaryotic genes, but a
prokaryotic transcript yielding several proteins needs an awkward patch: the genes derived from
it become children of a new feature, the operon, the transcript has several parent genes, and the
CDS pieces need a prokaryote-specific attribute (Derives_from=<a gene ID>) to name their gene.
This changes the meaning of the _implicit_ features above and needs different code to interpret.
The importer does not read this encoding.

Corner cases like trans-splicing, where one mature mRNA is ligated from two distant transcripts,
have no standard and no good encoding in GFF at all; the importer keeps such genes flattened and
never exports them (see [trans_splicing.md](trans_splicing.md)).

Attributes sit at the wrong level, too. Trans-splicing makes clear that neither the gene nor the
protein has on-genome coordinates; their pieces do, such as transcription or coding start and end
sites. Likewise, the protein ID is commonly (but not by any standard) derived from the gene or
transcript ID rather than assigned to the protein.

Some of this confusion comes from the loosely defined biological concept of a gene...
* a 'gene' is often a genomic locus, one ID over the related transcripts alternative splicing
  makes from it, i.e. an umbrella over all mRNA _transcribed_ from the locus and all proteins
  _translated_ from it.
* in prokaryotes, a 'gene' rather denotes a protein, a sub-section of what is _transcribed_.
* whether a protein from trans-splicing belongs to one or two 'genes' is unclear.

... but that is no excuse not to encode what is unambiguous clearly and consistently.

##### GeenuFF restructuring to bring the map closer to the territory
GeenuFF is a relational database schema (with an API to interpret it). It drops the artificial
connection rules of GFF-like formats, and its many-to-many fields directly allow several
transcripts per protein or several proteins per transcript. Its key tables:

* SuperLocus: an abstract holder for related transcripts, proteins and all the pieces that might
  be combined to make them, delineating the largest graph one might walk to put any component in
  context. It links (directly or indirectly) to:
  * Feature: "geenuff_transcript", "geenuff_cds", "geenuff_intron" and the like, with coordinates
    on the genome; the "geenuff_" prefix sets them apart from similarly named GFF features.
  * Transcript: (AKA pre-mRNA) has an ID and links (via [transcript_piece](#transcript_piece)) to
    all features making the transcript and any protein translated from it.
  * Protein: has an ID and links to the geenuff_cds features making it.

Eukaryotic, prokaryotic and even trans-spliced gene structures are then encoded the same way and
parsed by the same code:
* eukaryotic: a SuperLocus points to several Transcripts, each with one "geenuff_transcript" and
  one "geenuff_cds" feature and connected to one Protein.
* prokaryotic: a SuperLocus points to one Transcript with one "geenuff_transcript" feature,
  connected to several Proteins and with several "geenuff_cds" features.
* trans-splicing: a Transcript has one "geenuff_transcript" feature at each of its loci, ordered
  by the "position" attribute of its TranscriptPieces (see [transcript_piece](#transcript_piece)
  and [feature](#feature)).
* any of these across artificial breaks in the sequence (e.g. a fragmented assembly), using the
  "<>\_is_biological\_<>" attributes, masks where needed and TranscriptPieces.

("geenuff_cds", start) is always interpreted the same way, unlike a GFF CDS edge whose meaning
depends on eukaryote or prokaryote, on other CDS lines and on the overlapping exon. Likewise, the
same code reads a split-locus gene model whether the split is biological (trans-splicing) or
artificial (fragmented assembly): the features differ, the logic stays the same.

The importer builds only the first of these cases: one transcript piece and one protein per
transcript, from eukaryotic GFF3. It leaves out lines on another sequence than their parent, so no
transcript spans two scaffolds, and keeps trans-spliced genes flattened and out of exports. The
exports, Helixer's above all, rely on that: one transcript per gene, on one strand of one sequence.

#### gff-like coordinate troubles.
GFF's `[inclusive start, inclusive end]` coordinates work for the sequences of the transcripts and
proteins, but the _implicit_ components come out inverted and tedious: an intron is
`(exclude end exon_i, exclude start exon_i+1)`, the 3' UTR `[inclusive start exon, exclude start
CDS)`.

##### GeenuFF increased coordinate consistency
In GeenuFF a feature's start is always inclusive and its end always exclusive, as in python. A
component left implicit to avoid redundancy (e.g. the UTR as "geenuff_transcript" -
"geenuff_cds") is therefore always `[inclusive start, exclusive end <or start next>)`, e.g. the
first coding exon is `[geenuff_cds start, geenuff_intron start)`.

__Caveat:__
Ranges on the minus strand are off by one from the pythonic coordinates.

#### Extensible
Being a relational database, GeenuFF can be _extended_ with tables of its own, linked by foreign
keys to e.g. Coordinate or Genome, without changing the shared core.

## What (with details)

### Comparison of spec to gff
Quick start / reference for 1:1 comparison with gff:
[spec_vs_gff.html](spec_vs_gff.html)

### Coordinate system
GeenuFF coordinates count from 0, with the _start_ attribute marking the _inclusive_ start of a
range and the _end_ attribute its _exclusive_ end. On the plus strand they match python
coordinates; the same logic holds on the _minus_ strand.

To select the elements `{2, 3}` of `[0, 1, 2, 3, 4, ...]`: on the plus strand this is `[2, 4)`,
exactly the pythonic coordinates; on the minus strand (`{3, 2}`) it includes the start and not the
end, `[3, 1)`. With the genomic sequence as a python list, a minus strand range is selected with:
```
genomic_sequence[end + 1:start + 1]
```
It usually also needs reversing and complementing, which `applications.exporters.sequence` does.

### Schema summary
For the details see `base/orm.py`; this describes the overall structure and what it means.

#### Tables & relations
Indentation inside a piece indicates a one-to-? relation

* genome
  * coordinate
* super_locus
  * transcript, < many2many to protein >
    * transcript_piece, < many2many to feature >
    * transcript_error
  * protein, < many2many with transcript, feature >
  * feature, < many2many with protein, transcript_piece; many2one to coordinate >

(and linkage-only association tables for the many2many fields)

##### genome
This _mostly_ holds meta information

##### coordinate
Sequence meta info: seqid, length, sha1 hash of the sequence the annotation is for, optionally
the full sequence.

#### super_loci and children
###### super_loci
Delineates the graph of things that might possibly be combined, with given_name and type for the
~gene.

###### feature
Geenuff-specific types and ranges, otherwise resembling GFF features. Each has:

* a sequence (foreign key to coordinates)
* a position (start and end)
* a type, e.g. {geenuff_transcript, geenuff_cds, geenuff_intron, geenuff_mask}; the error types of
  a transcript are recorded per transcript instead (transcript_error)
* start_is_biological_start and end_is_biological_end, saying whether start and end mark the
  biologically meaningful transition or just the edge of what we know
* a boolean is_plus_strand
* a score (confidence)
* a phase, 0 for every feature: a CDS starts a fresh codon, the importer setting its phase itself
  whatever the file says (a different starting phase in the file is recorded as an error)
* a given_name (often "None", nothing being available)
* a source
* a many to many relationship with protein (mostly to assign the 'protein_id')
* a many to many relationship with transcript_piece

__Feature self consistency:__
All features on a transcript piece must be interpretable when sorted 5' to 3', and share the same
"is_plus_strand". An intronic (_cis_ or _trans_) range can only "start" or "end" in a transcribed
region, and a coding range only inside a transcript's non-intronic region. A feature whose
<>\_is_biological\_<> attribute is False should lie at the edge of the piece, or come with a
geenuff_mask feature masking the ambiguous area.

###### transcript_error
One row per error type found for a transcript, without a range; what the errors mask is merged
into the transcript's geenuff_mask features.

###### protein
Points to one protein's worth of "geenuff_cds" features (and associated transcript & super_locus),
and has a given_name.

###### transcript_piece
TranscriptPieces delineate the features that can be interpreted (5'-3') together. In the standard
case (_cis_- or no splicing, a full gene model) a Transcript has a single TranscriptPiece pointing
to all its Features.

In cases where the Transcript (final mRNA) is split for either biological or technical reasons,
each involved locus should have its own transcript_piece.

Each TranscriptPiece should have a single "geenuff_transcript" feature covering its whole range,
its _start_ at the most 5' and its _end_ at the most 3' part of the piece; only geenuff_mask
features may exceed it.

transcript_piece has a many2one relationship with transcript. The 5' to 3' ordering of
transcript_pieces within a transcript can be accomplished using the `position` attribute.

###### transcript
Transcripts define what is interpreted together to produce the final molecule (pre-mRNA, mRNA,
protein, ...). With their transcript_pieces sorted by ascending 'position' and the features within
each piece by coordinate, a transcript reads 5'-3', and whatever is of interest (the transcript
range, the splice sites, the start codon, ...) can be extracted. Example logic is in
`geenuff.applications.exporter.RangeMaker`.
