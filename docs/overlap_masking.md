# Overlapping genes

Two genes cannot share a base in an export: one label is written per base per strand, so
whichever is written second silently wins. Where one of two overlapping genes can be kept
without mislabelling the other's coding sequence, it is kept whole and its partner is left out
of exports. Where neither can be, both are masked over their whole length.

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

The gene owning coding sequence in the shared stretch is the one kept, so those bases stay
labelled as the coding sequence they are. The dropped gene's overhang is masked; the shared
stretch keeps the kept gene's labels, which are true for it. Identical spans need no mask at all.

```
  geneA   ==================            has CDS in the shared stretch
  geneB             ==================  its CDS starts further right
  shared            |<---->|

  result  ==================xxxxxxxxxx  geneA exported whole, geneB not exported
                            ##########  masked: what geneB covered beyond geneA
```

Where the dropped gene's far end is not where the gene really ended, the mask runs on past it,
as far as any other error mask would: the border with the next exported gene. A truncated CDS,
an in-frame stop, a missing start codon (or start codon in A if B would be kept) or an
unannotated UTR on that side all leave it unknown how much further the gene ran, so that sequence
cannot be taught as intergenic either.

```
  geneA   ==================            geneB has no annotated 3' UTR, so
  geneB             ==================  where it stopped is unknown
  geneC                                            ==========  next exported gene

  result  ==================xxxxxxxxxx             ==========
                            ###############
                                           ^ the mask reaches part way into the gap,
                                             not only to geneB's annotated end
```

### Nested, the outer kept

A nested gene can be masked away without touching the outer gene's start or end, the hole being
strictly interior. Not attempted when the outer gene also has CDS inside that span, which the
mask would cut into.

```
  geneA   ================================  outer, its intron spanning geneB
  geneB             xxxxxxxxxx              inner, not exported
  shared            |<------>|

  result  ==========##########============  geneA exported whole, masked only where geneB sat
```

### Refused, both masked whole

With coding sequence from both genes in the shared stretch there is nothing to choose: whichever
were kept, the other's CDS would read as UTR or intron. The same applies when the only gene that
could be kept is masked outright for its own errors, since keeping it recovers nothing.

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
  geneA          ==================            neither has an annotated outer UTR
  geneB                    ==================

  result    #######################                 geneA masked, reaching out to its left
                           #######################  geneB masked, reaching out to its right
```

A pair left alone as part of a chain is masked the same way, for the same reason.

## Which gene is kept

Two rules, in this order:

1. **Coding sequence wins.** If exactly one of the pair has CDS in the shared stretch, it is the
   keeper. If that gene cannot be kept, the pair is refused rather than the other one kept.
2. **Then the less damaged gene**, and on a tie the longer spliced CDS, which is the same measure
   that picks the longest isoform.

Damage is graded, not a yes/no:

| level           | what its own errors mask                       | may be kept?            |
|-----------------|------------------------------------------------|-------------------------|
| none            | nothing                                        | yes                     |
| flank masked    | the flank only, coding sequence still labelled | yes, below a clean gene |
| masked outright | the coding sequence itself                     | no, nothing to recover  |

*flank masked* is a missing UTR or a missing start or stop codon; *masked outright* is a
truncated CDS or an in-frame stop codon.

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

## What changed

- **GeenuFF v0.3.2** masked only when a missing-UTR error happened to coincide with the
  overlap, and placed that mask *beside* the shared stretch rather than on it, covering none of
  it. A nested gene got no mask at all, its missing-UTR errors deleted instead.
- **Now** an isolated pair keeps one gene whole where that can be done honestly, and where it
  cannot, both genes go entirely rather than being left with a hole each.

See `spec_vs_gff.md` for the error type list and for zero length error features.
