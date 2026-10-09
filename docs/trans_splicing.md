# Trans-spliced genes

A trans-spliced mRNA is assembled from pieces transcribed separately, so a gene annotation can
give pieces on opposite strands, out of order, or wrapping the origin of a circular molecule.
NCBI writes `exception=trans-splicing` on such a gene and often a strand of `?` on its mRNA line,
while the exon and CDS lines keep definite strands of their own.

GeenuFF keeps whatever of these genes it can and never exports any of it
(see [What survives when a strand is unusable](#what-survives-when-a-strand-is-unusable)).

## Why not exported

The exported arrays hold one label per base pair **per strand**, and HelixerPost's HMM walks one
strand in one direction. A transcript that leaves its strand halfway through has no representation
in either, so there is nothing to export it into.

## Why the introns are the part that breaks

Each exon and CDS piece has a definite strand, and its bases really are coding or untranslated on
it. The intron is the problem: an intron is sequence removed from *one* pre-mRNA, while
trans-spliced pieces come from separate pre-mRNAs, so the sequence between them is spliced out of
nothing and has no strand.

GeenuFF builds introns by chopping the exons out of one transcript span on one strand (see
`_parse_gff_entries`):

```python
itree[inv_t_start:inv_t_end] = t_is_plus_strand   # one span, one strand
for exon in exons:
    itree.chop(inv_e_start, inv_e_end)            # whatever is left over is an intron
```

For a trans-spliced gene that leftover includes the sequence between pieces, which is not intronic
and often holds other genes.

## What is actually stored today

The schema could represent these genes: strand lives on `Feature`, not on `SuperLocus`, and
`docs/index.md` describes `TranscriptPiece` for exactly this case:

> In cases where the Transcript (final mRNA) is split for either biological or technical reasons,
> each involved locus should have its own transcript_piece. […] The 5' to 3' ordering of
> transcript_pieces within a transcript can be accomplished using the `position` attribute.

The importer does not use it: one piece per transcript, every feature on the transcript's strand
(the gene's where the transcript has none), and no introns at all once an exon's strand disagrees.
Enough to see that a gene was annotated there, not enough to reconstruct it.

#### What survives when a strand is unusable

The parser reads a strand of `?` or `.` as no strand, and the importer meets it at four points:

| where `?` or `.` appears | what is stored                                                        | exported?                                 |
|--------------------------|-----------------------------------------------------------------------|-------------------------------------------|
| the **gene** line        | the gene record only, no transcripts and no features                  | no, nothing to select                     |
| the **mRNA** line        | the transcript and its features, on the gene's strand                 | no, `excluded_from_export` where it codes |
| an **exon** line         | the transcript and its features, without introns, none being invented | no, `excluded_from_export` where it codes |
| a **CDS** line           | the transcript and its features, introns included                     | no, `excluded_from_export` is set         |

A gene line on no strand is counted as `unstranded_super_loci`. A non-coding transcript on no
strand is stored on the gene's strand unmarked, as non-coding transcripts are never exported.

Storing them faithfully would take one `TranscriptPiece` per same-strand block, each with its own
`geenuff_transcript` feature and introns chopped only within the block. The database is ready for
that; the importer, the per-strand grouping in `clean_and_insert` and the one-transcript-per-gene
export are not.

## Does masking make sense

Not for the gene, and not for the gaps between its pieces:

* the gene's **envelope** spans from the first piece to the last and routinely contains unrelated,
  correctly annotated genes; masking it would destroy those.
* the **gaps between pieces** are neither introns nor part of the gene, so there is nothing to mask.

The pieces themselves would be a well-defined mask, but an excluded gene exports none of its
features, geenuff_mask included; masking the pieces would need a third state between exported and
excluded. So the pieces are taught as intergenic although they code. That mislabelling is accepted:
it affects a handful of genes per genome, a few kb each, the same kind of error any incompletely
annotated genome presents in far greater quantity.

See `spec_vs_gff.md` for the error types and `overlap_masking.md` for the other case where two
annotations claim the same sequence.
