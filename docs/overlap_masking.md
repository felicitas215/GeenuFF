# Overlapping genes

Two genes cannot share a base in an export: one label is written per base per strand, so
whichever is written second silently wins. Of two crossing genes the better one is kept whole and
its partner is left out of exports. Where both are masked outright for their own errors, or one
lies inside the other, both are masked over their whole length.

One coherent gene beats two with a hole through them: the genes either side of a hole have no
learnable terminus, and a consumer deriving boundaries from the label array reads its edges as
transitions that are not there. So a pair that cannot be settled loses both genes entirely rather
than only the sequence they share.

Only genes with a CDS take part. A bare `gene` line, a non-coding gene and a non-coding isoform
write no CDS, intron or UTR label, so they cannot collide with anything.

## The two records

|               | `super_locus_overlap` table        | `super_loci_overlap_error`                                         |
|---------------|------------------------------------|--------------------------------------------------------------------|
| covers        | every pair of coding genes         | only pairs of exported genes                                       |
| measured over | every coding transcript a gene has | the one transcript it exports                                      |
| masks         | nothing                            | its range, merged into that transcript's geenuff_mask features     |
| written for   | consumers other than Helixer       | consumers labelling one class per base (i.e. Helixer)              |

The table carries no judgement about whether either gene is usable. All overlaps of one gene are
`WHERE super_locus_id = X OR partner_id = X`.

## How a pair is settled

Only an **isolated pair** is decided, one where neither gene overlaps a third: in a chain a
dropped gene can sit between two kept ones, which one masked region per gene cannot express.
Chains are a few percent of overlaps in a nuclear genome and nearly all of them in an organellar
one, which should be excluded from Helixer training anyway.

`=` the exported transcript's span, `|<->|` the sequence the two share, `#` masked, `~` masked
flank (see `spec_vs_gff.md`), `x` dropped from the export entirely.

### Crossing, one gene kept

The better gene is kept (see "Which gene is kept") and the dropped gene's overhang is masked. The
shared stretch keeps the kept gene's labels, true for it, even where the dropped gene had coding
sequence there. Identical spans need no mask at all.

```
  geneA   ==================            the better gene
  geneB             ==================  error-free
  shared            |<---->|

  result  ==================xxxxxxxxxx  geneA exported whole, geneB not exported
                            ##########  masked: what geneB covered beyond geneA, no flank
```

Where the dropped gene has an error that masks something (a missing UTR or any error masking it
whole), its outer ends are untrustworthy, so the mask runs on into the flank on every side where
it sticks out past the kept gene, never toward the kept gene.

```
  geneA   ==================                                    the better gene
  geneB             ==================                          3' UTR missing, end unknown
  geneC                                               ========  next gene

  result  ==================xxxxxxxxxx                ========  kept: geneA (and geneC)
                            ##########~~~~~                     flank part way into the gap
```

Where no gene lies that way, the flank reaches toward the end of the sequence by the same rule,
`min(gap // 2, int(sqrt(gap)) * 10)`: a gene at the edge of a 10 Mb chromosome arm masks some
32 kb of it. One unannotated UTR should not cost a whole telomere.

```
  geneA   ==========                   |  no annotated 3' UTR; | end of the sequence
                                       |
  result  ==========~~~~               |  flank only, part way toward the end, the same reach
                                       |  as toward a next gene
```

### Nested, both masked whole

A gene lying wholly inside another on the same strand is more often an annotation mistake than
two real genes, so a nested pair is never decided: both are masked whole. Identical spans are not
nested and are settled like a crossing pair.

```
  geneP   ===                                                previous gene
  geneA             ==========================               outer
  geneB                     ==========                       inner
  geneN                                                 ===  next gene

  result  ===   ~~~~##########################~~~~      ===  geneA masked plus both flanks
                        ~~~~##########~~~~                   geneB likewise, inside geneA's mask
```

### Refused, both masked whole

Where both genes are masked outright for their own errors, keeping either recovers nothing, so
both are masked over their **whole length**: masking only the shared stretch would leave each with
a hole and save nothing.

```
  geneA        ==================
  geneB                  ==================
  shared                 |<---->|

  result   ~~~~##################~~~~            geneA masked end to end plus both flanks
                     ~~~~##################~~~~  geneB likewise
```

Nested, refused and chained genes all get the flank on both sides, errors of their own or not.
On the side facing the partner the flank is measured past it, so it ends inside the partner's
own mask and adds nothing.

## Which gene is kept

One ranking for a crossing pair, wherever either gene has coding sequence:

1. **The less damaged gene.** A gene masked outright is never kept; where both are, the pair is
   refused.
2. **On a tie, the longer spliced CDS**, the same measure that picks the longest isoform.
3. **On a tie of both, the 5'-most gene.**

The kept gene's labels cover the shared stretch, so coding sequence of the dropped gene there
reads as the kept gene's UTR, intron or CDS.

| level           | what its own errors mask                       | may be kept?            | errors                                                                                                                                                                                       |
|-----------------|------------------------------------------------|-------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| none            | nothing                                        | yes                     | none, or a wrong starting phase                                                                                                                                                              |
| flank masked    | the flank only, coding sequence still labelled | yes, below a clean gene | a missing UTR                                                                                                                                                                                |
| masked outright | the coding sequence itself                     | no, nothing to recover  | a missing start or stop codon, a truncated CDS, an in-frame stop codon, a truncated or too short intron, overlapping exon or CDS lines, floating CDS lines, a transcript beyond its sequence |

A missing UTR never disqualifies a gene from being kept: its own mask already covers the unknown
boundary, and the coding sequence stays labelled, as for any gene with an unannotated UTR. A
dropped gene graded anything but none gets the flank on the sides it sticks out past the kept gene.

## What a dropped gene leaves behind

Nothing is deleted: the gene and its features stay in the database as annotated, with
`SuperLocus.excluded_from_export = 'overlap_dropped'`, which `GeenuffExportController.genome_query`
skips, and the filtered GFF3 export too unless given `--include-erroneous`. Its errors are
recorded, but it gets no mask. It still bounds the flank of any gene it does not overlap, being an
annotated gene all the same.
