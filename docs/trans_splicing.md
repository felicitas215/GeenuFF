# Trans-spliced genes

A trans-spliced mRNA is assembled from pieces transcribed separately, so a gene annotation can
give pieces on opposite strands, out of order, or wrapping the origin of a circular molecule.
NCBI writes `exception=trans-splicing` on such a gene and often a strand of `?` on its mRNA line,
while the exon and CDS lines keep definite strands of their own.

GeenuFF keeps whatever of these genes it can and never exports any of it
(see [What survives when a strand is unusable](#what-survives-when-a-strand-is-unusable)).

## Why not exported

Neither Helixer nor HelixerPost's HMM can represent a transcript assembled from pieces on
different strands. The exported arrays hold one label per base pair **per strand**, and the HMM
walks one strand in one direction; a transcript that leaves the strand halfway through is not
something either can be handed in a form it could learn or decode. Exclusion is not a shortcut
around a hard case, it is the absence of any representation to export it into.

## Why the introns are the part that breaks

The exon and CDS pieces are not the problem: each carries a definite strand, and its bases
really are coding or untranslated on that strand.

The intron is the problem. An intron is sequence removed from *one* pre-mRNA. Trans-spliced pieces come
from separate pre-mRNAs, so the sequence lying between them was never transcribed as part of
either and is not spliced out of anything. Asking which strand that intron belongs on has no
answer because the intron does not exist.

GeenuFF only has to ask because of how it builds introns, by taking one transcript span on one
strand and chopping the exons out of it (see `_parse_gff_entries`):

```python
itree[inv_t_start:inv_t_end] = t_is_plus_strand   # one span, one strand
for exon in exons:
    itree.chop(inv_e_start, inv_e_end)            # whatever is left over is an intron
```

Everything left over becomes an intron. For a trans-spliced gene that leftover includes the
sequence between pieces, which is not intronic, and often other genes entirely.

## What is actually stored today

The schema can represent these genes properly. `SuperLocus` has no strand column at all, strand
living on `Feature`, and `docs/index.md` describes `TranscriptPiece` as the construct for exactly
this case:

> In cases where the Transcript (final mRNA) is split for either biological or technical reasons,
> each involved locus should have its own transcript_piece. […] The 5' to 3' ordering of
> transcript_pieces within a transcript can be accomplished using the `position` attribute.

The importer does not use it. It creates one piece per transcript, puts every feature on the
gene's strand, and gives up on introns entirely as soon as an exon's strand disagrees. So what is
written for a mixed-strand gene is a flattened approximation: the transcript and its CDS, some of
it on a strand it does not belong on, and no introns at all. It is enough to see that a gene was
annotated there, not enough to reconstruct it.


#### What survives when a strand is unusable
How much survives depends on where the unusable strand sits. GFF3 writes a strand that is
neither `+` nor `-` as `?` or `.`, both of which the parser reads as no strand at all, and the
importer meets them at three different points:

| where `?` or `.` appears    | what is stored                                                                                                                                                            | exported?                                                            |
|-----------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------|
| the **gene** line           | nothing but the gene record: the importer returns before building any transcript, so the locus has no transcripts and no features whatsoever                              | no, there is nothing to export and no `longest` transcript to select |
| the **mRNA** line           | the transcript and its features, on the gene's strand                                                                                                                     | no, `excluded_from_export` is set where the transcript codes         |
| an **exon** or **CDS** line | the same, minus the introns: without every exon placed there is nothing to subtract them from, so none are built rather than one being invented over the whole transcript | no, `excluded_from_export` is set                                    |

The gene line therefore breaks earliest and hardest, and it is the only one of the three that
needs no `excluded_from_export` at all: a locus with no transcripts has nothing an export could
select. It is counted separately, as `unstranded_super_loci` rather than under the unexported
reasons, which is why the two numbers in the import summary do not overlap.

A non-coding transcript with an unusable strand is stored on the gene's strand and *not* marked,
since it is never exported anyway, no non-coding transcript being selected in the first place.

Storing them faithfully would mean one `TranscriptPiece` per maximal same-strand block, each with
its own `geenuff_transcript` feature on that block's strand and introns chopped only within the
block. Every intron then has an unambiguous strand, its block's, and nothing is invented between
blocks. The database is ready for that; the importer, the per-strand grouping in
`clean_and_insert` and the one-transcript-per-gene assumption in the export path are not.

## Does masking make sense

Not for the gene, and not for the gaps between its pieces:

* the gene's **envelope** spans everything from the first piece to the last, which for a
  trans-spliced gene routinely contains unrelated, correctly annotated genes. Masking it would
  destroy those to hide one.
* the **gaps between pieces** are not introns and not part of the gene, as above. There is
  nothing there to mask.

The one range that *would* be well-defined is the pieces themselves: each has a definite strand
and definite bounds, so "something is annotated here that cannot be interpreted" is a true and
bounded statement about them. That is the honest mask, if one were wanted.

It cannot be written today, though, and the reason is worth knowing: excluding a locus removes
**all** of its features from the export, geenuff_mask features included, so an excluded gene cannot
carry a mask. Masking the pieces would need a third state between exported and excluded, one
where a locus contributes its masks but not its annotation. That does not exist.

So the cost of exclusion is that the pieces of such a gene are taught as intergenic, when they
are really coding. That is a genuine mislabelling, and it is accepted on grounds of scale: it
affects a handful of genes per genome, a few kb each, and it is the same kind of error as any
incompletely annotated genome already presents in far greater quantity. It is not worth a new
export state to fix, but it should be known rather than assumed away.

See `spec_vs_gff.md` for the error types and `overlap_masking.md` for the other case where two
annotations claim the same sequence.
