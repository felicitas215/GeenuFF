# Overlapping genes

Two genes cannot share a base in an export: one label is written per base per strand, so
whichever is written second silently wins. Of two overlapping genes the better one is kept whole
and its partner is left out of exports, whichever of them has coding sequence where they
overlap. Only where both are masked outright for their own errors are both masked over their
whole length.

One coherent gene beats two with a hole through them, and a gene with a hole is worth no more
than no gene at all. A hole costs more than its base pairs: the gene either side of it has no
learnable terminus, and a consumer deriving feature boundaries by diffing the label array reads
the edges of the hole as transitions that are not there. That is why a pair that cannot be
settled loses both genes entirely rather than only the sequence they share.

Only genes with a CDS take part in any of this. A bare `gene` line with no transcript under it,
a purely non-coding gene, and a non-coding isoform of a coding gene are all invisible here: none
of them writes a CDS, intron or UTR label anywhere, so none of them can collide with anything.

## The two records

|               | `super_locus_overlap` table        | `super_loci_overlap_error` feature                    |
|---------------|------------------------------------|-------------------------------------------------------|
| covers        | every pair of coding genes         | only pairs of exported genes                          |
| measured over | every coding transcript a gene has | the one transcript it exports                         |
| masks         | nothing                            | the range it covers                                   |
| written for   | consumers other than Helixer       | consumers labelling one class per base (i.e. Helixer) |

The table is the wider record and carries no judgement about whether either gene is usable. All
overlaps of one gene are `WHERE super_locus_id = X OR partner_id = X`.

## How a pair is settled

Only an **isolated pair** is decided, i.e. one where neither gene overlaps a third. In a chain a
dropped gene can sit between two kept ones, and the single masked region recorded per gene
cannot express a hole split into pieces. Chains are a few percent of overlaps in a nuclear
genome and nearly all of them in an organellar one (organelles = something that should be excluded
from Helixer training).

`=` the exported transcript's span, `|<->|` the sequence the two share, `#` masked, `x` dropped
from the export entirely.

### Crossing, one gene kept

The better gene is kept (see "Which gene is kept"). The dropped gene's overhang is masked; the
shared stretch keeps the kept gene's labels, which are true for it, even where the dropped gene
had coding sequence there. Identical spans need no mask at all.

```
  geneA   ==================            the better gene
  geneB             ==================
  shared            |<---->|

  result  ==================xxxxxxxxxx  geneA exported whole, geneB not exported
                            ##########  masked: what geneB covered beyond geneA
```

Where the dropped gene's far end is not where the gene really ended, the mask runs on past it,
as far as any other error mask would: its flank into the gap toward the next exported gene,
measured from the dropped gene itself (see `spec_vs_gff.md`). An unannotated UTR on that side,
or any error masking the gene whole (a missing start or stop codon, a truncated CDS or intron, an
in-frame stop or a too short intron), leaves it unknown how much further the gene ran, so that
sequence cannot be taught as intergenic either.

```
  geneA   ==================            geneB has no annotated 3' UTR, so
  geneB             ==================  where it stopped is unknown
  geneC                                            ==========  next exported gene

  result  ==================xxxxxxxxxx             ==========  kept: gene A (and C)
                            ###############
                                           ^ the mask reaches part way into the gap,
                                             not only to geneB's annotated end
```

Where no exported gene lies that way at all, the end of the sequence stands in for one and the
mask reaches part way toward it by the same rule. It does not run to the end. The offset is
`min(gap // 2, int(sqrt(gap)) * 10)`, whose square root term grows slowly enough to be a ceiling in
itself: a gene at the edge of a 10 Mb chromosome arm masks some 32 kb of it, not all 10 Mb. One
unannotated UTR on the last gene of a chromosome should not cost a whole telomere, and sequence
far enough from a gene is intergenic whatever that gene did.

```
  geneA   ==========                          | end of the sequence
                                              |
  result  ==========####                      |
                        ^ part way toward the end, the same offset as toward a next gene
                        (gene A is missing a 3' UTR -> flank mask only)
```

### Nested, one gene kept

Where the outer gene is the better one, the inner gene is masked where it sat, without touching
the outer gene's start or end, the hole being strictly interior. This holds even where the outer
gene has CDS inside that span, which the hole then covers.

```
  geneA   ================================  outer, the better gene
  geneB             xxxxxxxxxx              inner, not exported
  shared            |<------>|

  result  ==========##########============  geneA exported whole, masked only where geneB sat
```

Where the inner gene is the better one, the outer gene is dropped and what it covers beyond the
inner gene is masked on both sides. Each side runs on past the outer gene's end where that end is
unknown, as for a crossing pair.

```
  geneA   ================================  outer, not exported
  geneB             ==========              inner, the better gene

  result  ##########==========##########    geneB exported whole, geneA's overhangs masked
```

### Refused, both masked whole

Where both genes are masked outright for their own errors, keeping either recovers nothing.

Both genes are then masked over their **whole length**, not only over what they share. Masking
just the shared stretch would leave each of them with a hole through it, and a gene with a hole
teaches neither terminus nor continuity across it. Neither gene is usable, so there is nothing
left for a smaller mask to save.

```
  geneA   ==================
  geneB             ==================
  shared            |<---->|

  result  ##################                both exported, but each masked end to end
                    ##################
```

The same rule about unknown ends applies here as to a dropped gene: where a gene's outer end is
not where it really ended, its mask runs on past its span toward the next exported gene. The
ends the two face each other with need no such treatment, the partner leaving no room.

```
  geneA          ==================            neither has an annotated outer UTR or other severe errors
  geneB                    ==================

  result    #######################                 geneA masked, reaching out to its left
                           #######################  geneB masked, reaching out to its right
```

A pair left alone as part of a chain is masked the same way, for the same reason.

## Which gene is kept

One ranking, whatever the shape of the pair and wherever either gene has coding sequence:

1. **The less damaged gene.** A gene masked outright is never kept; where both are, the pair is
   refused.
2. **On a tie, the longer spliced CDS**, the same measure that picks the longest isoform.
3. **On a tie of both, the 5'-most gene.**

Keeping a gene means its labels cover the shared stretch, so coding sequence of the dropped gene
lying there reads as the kept gene's UTR, intron or CDS. Of two annotations that cannot both be
written, the better one is taken whole.

Damage is graded, not a yes/no:

| level           | what its own errors mask                       | may be kept?            |
|-----------------|------------------------------------------------|-------------------------|
| none            | nothing                                        | yes                     |
| flank masked    | the flank only, coding sequence still labelled | yes, below a clean gene |
| masked outright | the coding sequence itself                     | no, nothing to recover  |

*flank masked* is a missing UTR; *masked outright* (+ flank if applicable) is a missing start or stop codon, a truncated
CDS, an in-frame stop codon, a truncated intron, a too short intron or overlapping exons. A wrong
starting phase masks nothing and counts as none. The same errors decide whether a gene's end is
known, i.e. whether its overlap mask runs on past that end (see the table in `spec_vs_gff.md`).

**A missing UTR never disqualifies a gene from being kept.** Its unknown boundary is already
covered by its own missing-UTR mask, which on the side facing the partner simply runs into the
dropped gene's masked region. Keeping such a gene gives exactly what GeenuFF gives any gene with
an unannotated UTR: the coding sequence labelled, the uncertain flank masked.

## What a dropped gene leaves behind

Nothing is deleted. The gene, its transcripts and all their features stay in the database
exactly as annotated; only exports change. `SuperLocus.excluded_from_export` is set to
`overlap_dropped`, and both `GeenuffExportController.genome_query` and the filtered GFF3 export
skip such a locus. Dropping is an export decision, not an annotation one.

That also means masks measured toward "the next gene" skip it. A mask running part of the way to
its neighbour is measuring toward whatever a consumer will be handed, so an unexported gene does
not bound one; the border falls between the two genes actually written.

