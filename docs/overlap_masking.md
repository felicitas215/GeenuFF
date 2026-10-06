# Overlapping genes

Two genes cannot share a base in an export: one label is written per base per strand, so
whichever is written second silently wins. Of two crossing genes the better one is kept whole
and its partner is left out of exports, whichever of them has coding sequence where they
overlap. Where both are masked outright for their own errors, or one lies inside the other, both
are masked over their whole length.

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

`=` the exported transcript's span, `|<->|` the sequence the two share, `#` masked, `~` masked
flank (see `spec_vs_gff.md`), `x` dropped from the export entirely.

### Crossing, one gene kept

The better gene is kept (see "Which gene is kept"). The dropped gene's overhang is masked; the
shared stretch keeps the kept gene's labels, which are true for it, even where the dropped gene
had coding sequence there. Identical spans need no mask at all.

```
  geneA   ==================            the better gene
  geneB             ==================  error-free
  shared            |<---->|

  result  ==================xxxxxxxxxx  geneA exported whole, geneB not exported
                            ##########  masked: what geneB covered beyond geneA, no flank
```

Where the dropped gene has any error that masks something (a missing UTR or any error masking it
whole), its outer ends are no more trustworthy than those of any other erroneous gene, so the
mask runs on past it into the flank, as far as any other error mask would (see `spec_vs_gff.md`),
on every side where it sticks out past the kept gene. Never toward the kept gene, whose labels
are trusted.

```
  geneA   ==================                                    the better gene
  geneB             ==================                          no annotated 3' UTR, so where it stopped is unknown
  geneC                                               ========  next gene

  result  ==================xxxxxxxxxx                ========  kept: geneA (and geneC)
                            ##########~~~~~                     flank part way into the gap, never toward geneA
```

Where no gene lies that way at all, the end of the sequence stands in for one and the mask
reaches part way toward it by the same rule. It does not run to the end. The flank reaches
`min(int(10 * sqrt(gap)), midpoint of the gap)`, whose square root term grows slowly enough to be
a ceiling in itself: a gene at the edge of a 10 Mb chromosome arm masks some 32 kb of it, not all
10 Mb. One
unannotated UTR on the last gene of a chromosome should not cost a whole telomere, and sequence
far enough from a gene is intergenic whatever that gene did.

```
  geneA   ==========                   |  no annotated 3' UTR; | end of the sequence
                                       |
  result  ==========~~~~               |  flank only, part way toward the end, the same reach
                                       |  as toward a next gene
```

### Nested, both masked whole

A gene lying wholly inside another on the same strand is more often an annotation mistake than
two real genes, so a nested pair is never decided: both are masked whole, as in the next
section. Identical spans are not nested, there being no inner gene; they are settled like a
crossing pair, which needs no mask at all.

```
  geneP   ===                                                previous gene
  geneA             ==========================               outer
  geneB                     ==========                       inner
  geneN                                                 ===  next gene

  result  ===   ~~~~##########################~~~~      ===  geneA masked end to end plus both flanks
                        ~~~~##########~~~~                   geneB likewise, inside geneA's mask
```

Both get the flank on both sides whether they have errors of their own or not, like any pair
that is not decided (see below).

### Refused, both masked whole

Where both genes are masked outright for their own errors, keeping either recovers nothing.

Both genes are then masked over their **whole length**, not only over what they share. Masking
just the shared stretch would leave each of them with a hole through it, and a gene with a hole
teaches neither terminus nor continuity across it. Neither gene is usable, so there is nothing
left for a smaller mask to save.

```
  geneA        ==================
  geneB                  ==================
  shared                 |<---->|

  result   ~~~~##################~~~~            geneA masked end to end plus both flanks
                     ~~~~##################~~~~  geneB likewise
```

Like any gene masked whole, each also gets the flank on both sides. On the side facing its
partner the flank is measured past the partner, which overlaps it and so does not bound it. That
flank ends inside the partner's own mask, which starts from further out and so reaches further,
so it adds nothing to what is masked.

A nested pair and a pair left alone as part of a chain are masked the same way, flanks included,
whether the genes have errors of their own or not.

## Which gene is kept

One ranking for a crossing pair, wherever either gene has coding sequence:

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
starting phase masks nothing and counts as none. A dropped gene graded anything but none gets
the flank on the sides it sticks out past the kept gene (see the table in `spec_vs_gff.md`).

**A missing UTR never disqualifies a gene from being kept.** Its unknown boundary is already
covered by its own missing-UTR mask, which on the side facing the partner simply runs into the
dropped gene's masked region. Keeping such a gene gives exactly what GeenuFF gives any gene with
an unannotated UTR: the coding sequence labelled, the uncertain flank masked.

## What a dropped gene leaves behind

Nothing is deleted. The gene, its transcripts and all their features stay in the database
exactly as annotated; only exports change. `SuperLocus.excluded_from_export` is set to
`overlap_dropped`, and both `GeenuffExportController.genome_query` and the filtered GFF3 export
skip such a locus. Dropping is an export decision, not an annotation one.

A dropped gene still bounds the flank of any gene it does not overlap. It is an annotated gene
all the same, which another gene's true end is not assumed to run through, and its own region is
masked already.

