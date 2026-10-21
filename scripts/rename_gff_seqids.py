#! /usr/bin/env python3
"""Rewrites the sequence IDs of a GFF3 file to those of the FASTA file it annotates, for when the
sequences were renamed after annotation, e.g. on upload to NCBI. A GFF3 sequence ID is paired with
the FASTA sequence of the same ID (its header line up to the first whitespace), else with the one
whose full name (the whole header line) contains it as a word, else with the one its genes fit; every
pair is checked against the sequence, at the start and stop codons and the splice sites the GFF3
places on it. Lines on sequence IDs without a FASTA sequence are left out. Sequences changed after
annotation are left out of both files, the FASTA file then being copied without them.
"""
import argparse
import bisect
import itertools
import logging
import math
import os
import re
from collections import defaultdict
from dataclasses import dataclass, field

from geenuff.base.helpers import (DEFAULT_COMPRESSION_LEVEL, DEFAULT_THREADS, Compression, add_output_compression_args,
                                  check_output_paths, codon_after_cds, compression_for, insert_before_extension,
                                  open_out, output_compression, output_path, require_fasta_and_gff,
                                  require_input_files, setup_logging, write_fasta_without)
from geenuff.base.types import PIECE_KEYS
from geenuff.parsing.fastahelper import FastaParser, text_or_gzip_open
from geenuff.parsing.gffhelper import GFFObject, GFFParser

logger = logging.getLogger(__name__)

MOTIFS = {'start': {'ATG'}, 'stop': {'TAA', 'TAG', 'TGA'}, 'donor': {'GT', 'GC'}, 'acceptor': {'AG'}}
# a stop codon at the end of the CDS or right after it, where GFF3 converted from GTF leaves it out
STOP_OR_NEXT = 'stop_or_next'
MOTIF_LENGTHS = {'donor': 2, 'acceptor': 2}
# the chance of a random position passing each check
RANDOM_PASS = {'start': 1 / 64, 'stop': 3 / 64, STOP_OR_NEXT: 6 / 64, 'donor': 2 / 16, 'acceptor': 1 / 16}
COMPLEMENT = str.maketrans('ACGT', 'TGCA')
# (kind, 0-based position on the plus strand, is_plus), a codon check with the 0-based positions of
# its bases in reading order added, read across an intron where it is split by one: the start codon,
# or for STOP_OR_NEXT the stop codon and the codon after the CDS (see cds_end_codons, codon_after)
Check = tuple[str, int, bool] | tuple[str, int, bool, tuple[int, ...]]
# A check is one position where the GFF3 file places a start codon, a stop codon or one end of an
# intron, tested for the bases expected there; a transcript with n introns gives 2 + 2n checks.
# At most MAX_CHECKS per sequence ID, spread evenly along it, far more than needed to tell its sequence
# from others
MAX_CHECKS = 250
MIN_NAMED_SCORE = 0.5  # pairs by ID or full name; a random sequence passes about a tenth of the checks
MIN_SEARCHED_SCORE = 0.8
MAX_RANDOM_HITS = 1e-3  # expected random sequences passing the checks of a pair found by search
# the number of checks (not a share) a pair by ID or full name, or a window of its checks, needs to be
# judged by them. A wrong sequence passes about 9% of the checks, so with 20 it passes half by chance
# in about 2e-6 of cases, and a right one with a poor annotation passing 70% fails half in about 2%;
# with fewer, partial genes at a sequence end could fail most checks of a right pair
MIN_DECISIVE_CHECKS = 20
# a pair by ID or full name is judged in windows of about this many consecutive checks along the
# sequence, each of which must pass; after a piece cut out of or put into the FASTA sequence all
# checks fail, however many before it pass
WINDOW_CHECKS = 25
MAX_CHANGED_SHARE = 0.05  # of the coding transcripts, on sequences changed after annotation
MAX_MISSING_SHARE = 0.15  # of the coding transcripts, on sequences matching no FASTA sequence
OUT_SUFFIX = '_renamed'  # inserted before the extension of the GFF3 file's name for the default output
ADVICE = ('Renaming cannot fix this; use the FASTA file the GFF3 file was made on, or transfer the annotation '
          'to this assembly (e.g. with Liftoff).')


class MismatchedFiles(Exception):
    pass


@dataclass
class StopSite:
    transcript_id: str
    codon: tuple[int, int, int]  # the last codon of the CDS (see cds_end_codons)
    is_plus: bool
    exon_beyond: bool  # whether its exons reach 3 bases or more past the CDS
    after: tuple[int, int, int]  # the codon following the CDS (see codon_after)


@dataclass
class GffSequence:
    max_end: int = 0
    region_length: int | None = None  # from a ##sequence-region line
    n_lines: int = 0
    coding_transcripts: int = 0
    checks: list[Check] = field(default_factory=list)
    stops: list[StopSite] = field(default_factory=list)


@dataclass
class Pair:
    fasta_id: str
    method: str  # 'id', 'full_name' or 'sequence'
    passed: int
    checks: int


@dataclass
class Pairing:
    pairs: dict[str, Pair]  # by GFF3 sequence ID
    changed: dict[str, tuple[str, str]]  # GFF3 sequence ID: (FASTA ID, why their genes do not fit it)
    duplicates: dict[str, str] = field(default_factory=dict)  # GFF3 sequence ID: FASTA ID paired already


def cds_end_codons(cds: list[tuple[int, int]], is_plus: bool) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """The 0-based plus-strand positions of a CDS's first and last codon, each in reading order and read
    along the CDS lines, so across an intron where it is split by one."""
    # the 3 lowest and the 3 highest CDS bases, taken from the CDS lines without listing them all
    ascending = (p for start, end in sorted(cds) for p in range(start, end + 1))
    descending = (p for start, end in sorted(cds, reverse=True) for p in range(end, start - 1, -1))
    low, high = list(itertools.islice(ascending, 3)), list(itertools.islice(descending, 3))
    # reading runs up the plus strand and down the minus strand, so the start codon is the low end on
    # the plus strand and the high end on the minus strand
    first, last = (low, high[::-1]) if is_plus else (high, low[::-1])
    return tuple(p - 1 for p in first), tuple(p - 1 for p in last)


def codon_after(cds: list[tuple[int, int]], exons: list[tuple[int, int]], is_plus: bool) -> tuple[int, int, int]:
    """The 0-based plus-strand positions of the 3 bases following a transcript's CDS, in reading order,
    as GeenuFF's import reads them (see geenuff.base.helpers.codon_after_cds)."""
    return tuple(p - 1 for p in codon_after_cds(cds, exons, is_plus)[0])


def transcript_checks(cds: list[tuple[int, int]], exons: list[tuple[int, int]], is_plus: bool) -> list[Check]:
    """The start and stop codon of a transcript's CDS and the splice sites of the introns between its
    exons, or its CDS lines without exon lines, as checks (see Check).

    Args:
        cds: the GFF3 (start, end) of its CDS lines
        exons: the GFF3 (start, end) of its exon lines
        is_plus: whether it lies on the plus strand
    """
    checks = []
    if cds:
        start_codon, stop_codon = cds_end_codons(cds, is_plus)
        after = codon_after(cds, exons, is_plus)
        # each sorted along the sequence by its lowest base
        checks += [('start', min(start_codon), is_plus, start_codon),
                   (STOP_OR_NEXT, min(stop_codon), is_plus, stop_codon + after)]
    pieces = sorted(exons or cds)
    # the intron's low end is the donor site (GT) on the plus strand, the acceptor site (AG) on the
    # minus strand, read there as the reverse complement
    left, right = ('donor', 'acceptor') if is_plus else ('acceptor', 'donor')
    for (_, end), (start, _) in zip(pieces, pieces[1:]):
        # a gap of under 4 bases is no intron; the checks lie on its first 2 bases (0-based, the base
        # after the piece's end) and on its last 2, right before the next piece, i.e. GTAG (both splice
        # sites lumped together, technically too short, but good enough for this test)
        if start - end > 4:
            checks += [(left, end, is_plus), (right, start - 3, is_plus)]
    return checks


def stop_site(transcript_id: str, cds: list[tuple[int, int]], exons: list[tuple[int, int]],
              is_plus: bool) -> StopSite:
    """Where a transcript's CDS ends, and whether its exons go on 3 bases or more past it, which GeenuFF
    needs to extend the CDS by a stop codon left out of it."""
    codon, after = cds_end_codons(cds, is_plus)[1], codon_after(cds, exons, is_plus)
    # the exon bases lying past the 3' end of the CDS, summed over the exons
    if is_plus:
        cds_end = max(end for _, end in cds)
        beyond = sum(max(0, end - max(start - 1, cds_end)) for start, end in exons)
    else:
        cds_start = min(start for start, _ in cds)
        beyond = sum(max(0, min(end + 1, cds_start) - start) for start, end in exons)
    return StopSite(transcript_id, codon, is_plus, beyond >= 3, after)


def sequence_regions(gff_path: str) -> dict[str, int]:
    """{sequence ID: length} from the ##sequence-region lines anywhere in a GFF3 file."""
    regions = {}
    with text_or_gzip_open(gff_path) as f:
        for line in f:
            if line.startswith('##sequence-region'):
                _, seqid, start, end = line.split()[:4]
                regions[seqid] = int(end) - int(start) + 1
    return regions


def read_gff(gff_path: str) -> dict[str, GffSequence]:
    """The sequence IDs of a GFF3 file's feature lines, with what is needed to pair and check them."""
    sequences = defaultdict(GffSequence)
    # by (seqid, Parent); the exon and CDS types GeenuFF imports, a line under each transcript it names
    pieces = defaultdict(lambda: {'cds': [], 'exons': [], 'strands': set()})
    for line in GFFParser.feature_lines(gff_path):
        seqid, _, feature_type, start, end, _, strand, _, column = line.rstrip('\n').split('\t')[:9]
        sequences[seqid].max_end = max(sequences[seqid].max_end, int(end))
        sequences[seqid].n_lines += 1
        key = PIECE_KEYS.get(feature_type)
        if key not in ('cds', 'exons'):
            continue
        for parent in GFFObject.parents(GFFObject.parse_attributes(column)):
            transcript = pieces[(seqid, parent)]
            transcript[key].append((int(start), int(end)))
            transcript['strands'].add(strand)
    for (seqid, transcript_id), transcript in pieces.items():
        sequence = sequences[seqid]
        sequence.coding_transcripts += bool(transcript['cds'])
        if transcript['strands'] not in ({'+'}, {'-'}):
            continue
        is_plus = transcript['strands'] == {'+'}
        sequence.checks += transcript_checks(transcript['cds'], transcript['exons'], is_plus)
        if transcript['cds']:
            sequence.stops.append(stop_site(transcript_id, transcript['cds'], transcript['exons'], is_plus))
    for sequence in sequences.values():
        # MAX_CHECKS spread evenly along the sequence, so a change anywhere on it shows, between its first
        # and last MIN_DECISIVE_CHECKS checks, those of its first and last genes (see failing_window)
        checks = sorted(sequence.checks, key=lambda check: check[1])
        edge = MIN_DECISIVE_CHECKS
        if len(checks) > MAX_CHECKS:
            # the two edge checks are kept whole, from the middle every len(middle) / n_spread-th one
            # (i.e. sampled at even steps)
            middle, n_spread = checks[edge:-edge], MAX_CHECKS - 2 * edge
            checks = (checks[:edge] + [middle[i * len(middle) // n_spread] for i in range(n_spread)]
                      + checks[-edge:])
        sequence.checks = checks
    for seqid, length in sequence_regions(gff_path).items():
        if seqid in sequences:
            sequences[seqid].region_length = length
    return sequences


def codon_at(sequence: str, positions: tuple[int, ...], is_plus: bool) -> str:
    """The bases at 0-based plus-strand positions in reading order, complemented on the minus strand,
    or '' where one lies outside the sequence."""
    if not all(0 <= p < len(sequence) for p in positions):
        return ''
    bases = ''.join(sequence[p] for p in positions)
    return bases if is_plus else bases.translate(COMPLEMENT)


def passes(sequence: str, kind: str, position: int, is_plus: bool, codons: tuple[int, ...] | None = None) -> bool:
    """Whether the bases a check tests are the motif expected there (see Check): a codon check reads its
    codon's positions, a STOP_OR_NEXT check passing where the stop codon or the codon after it is one."""
    if kind == 'start':
        return codon_at(sequence, codons, is_plus) in MOTIFS['start']
    if kind == STOP_OR_NEXT:
        # codons holds the stop codon's 3 positions, then those of the codon after the CDS
        # (in case the stop codon lies after the CDS, i.e. GTF -> GFF3 conversion)
        return any(codon_at(sequence, codons[i:i + 3], is_plus) in MOTIFS['stop'] for i in (0, 3))
    if position < 0:
        return False
    # a splice site's 2 contiguous bases, read as the reverse complement on the minus strand
    motif = sequence[position:position + MOTIF_LENGTHS[kind]]
    if not is_plus:
        motif = motif.translate(COMPLEMENT)[::-1]
    return motif in MOTIFS[kind]


def n_passed(sequence: str, checks: list[Check]) -> int:
    return sum(passes(sequence, *check) for check in checks)


def random_pass_log10(checks: list[Check], n: int) -> float:
    """log10 of an upper bound on the chance that a random sequence passes n of the checks."""
    # (number of ways to pick the n passing checks) times (the chance of the n easiest to pass all
    # passing), which bounds every pick from above
    log_chances = sorted(math.log10(RANDOM_PASS[kind]) for kind, *_ in checks)
    return math.log10(math.comb(len(checks), n)) + sum(log_chances[len(checks) - n:])


def failing_window(sequence: str, checks: list[Check]) -> str | None:
    """The first window of consecutive checks along the sequence (see WINDOW_CHECKS) of which fewer
    than MIN_NAMED_SCORE pass, described, or None where all pass; checks are sorted by position. The
    first and the last MIN_DECISIVE_CHECKS checks, those of the first and last genes, are windows of
    their own, a change near an end of the sequence failing too few checks to sink a whole window."""
    edge = MIN_DECISIVE_CHECKS
    # windows of nearly equal size, about WINDOW_CHECKS each, covering all checks, plus the two edges
    n_windows = max(1, len(checks) // WINDOW_CHECKS)
    bounds = [len(checks) * i // n_windows for i in range(n_windows + 1)]
    windows = [(0, edge)] + list(zip(bounds, bounds[1:])) + [(len(checks) - edge, len(checks))]
    for lo, hi in windows:
        window = checks[lo:hi]
        passed = n_passed(sequence, window)
        if passed / len(window) < MIN_NAMED_SCORE:
            return (f'{passed} of {len(window)} start and stop codons and splice sites in place from '
                    f'{window[0][1] + 1} to {window[-1][1] + 3} bp')
    return None


def full_name_words(full_names: dict[str, str]) -> dict[str, list[str]]:
    """{word of a FASTA full name: [IDs of the FASTA sequences whose full name holds it]}"""
    words = defaultdict(list)
    for fasta_id, full_name in full_names.items():
        for word in set(re.split(r'[\s,;:|=]+', full_name)):
            words[word].append(fasta_id)
    return words


def pair_by_id_or_full_name(gff_sequences: dict[str, GffSequence], sequences: dict[str, str],
                            full_names: dict[str, str]) -> tuple[dict[str, Pair], dict[str, tuple[str, str]], int]:
    """Pairs each GFF3 sequence ID with the FASTA sequence of the same ID, else the only one whose
    full name holds it as a word and whose ID is no GFF3 sequence ID, which would pair with it by
    ID. A pair fails where the GFF3 genes reach past the FASTA sequence,
    its ##sequence-region length differs from the FASTA sequence's, or, where the GFF3 sequence has
    at least MIN_DECISIVE_CHECKS checks (a count), fewer than MIN_NAMED_SCORE of them (a share) pass
    in any window of them (see failing_window).

    Returns:
        the pairs, {GFF3 sequence ID: (FASTA ID, why it failed)} for the failed ones, and how many
        pairs were tried
    """
    words = full_name_words(full_names)
    pairs, claimed, failed = {}, set(), {}
    n_tried = 0
    for gff_id, gff_sequence in gff_sequences.items():
        if gff_id in sequences:
            fasta_id, method = gff_id, 'id'
        # one full name only holding it, of a FASTA sequence no other GFF3 ID pairs with by ID
        elif len(words.get(gff_id, ())) == 1 and words[gff_id][0] not in gff_sequences:
            fasta_id, method = words[gff_id][0], 'full_name'
        else:
            continue
        sequence, checks = sequences[fasta_id], gff_sequence.checks
        n_tried += 1
        window = failing_window(sequence, checks) if len(checks) >= MIN_DECISIVE_CHECKS else None
        if gff_sequence.max_end > len(sequence):
            failed[gff_id] = fasta_id, f'genes end at {gff_sequence.max_end}, past its {len(sequence)} bp'
        elif gff_sequence.region_length not in (None, len(sequence)):
            failed[gff_id] = fasta_id, (f'{gff_sequence.region_length} bp by ##sequence-region, {len(sequence)} bp '
                                        f'in the FASTA file')
        elif window is not None:
            failed[gff_id] = fasta_id, window
        # a FASTA sequence pairs with one GFF3 ID; a second is left to the search by sequence, which
        # takes it for a duplicate
        elif fasta_id not in claimed:
            passed = n_passed(sequence, checks)
            pairs[gff_id] = Pair(fasta_id, method, passed, len(checks))
            claimed.add(fasta_id)
    return pairs, failed, n_tried


def share_passed(pair: Pair) -> float:
    return pair.passed / pair.checks


def pair_by_sequence(gff_sequences: dict[str, GffSequence], sequences: dict[str, str],
                     claimed: set[str]) -> tuple[dict[str, Pair], dict[str, str]]:
    """Pairs each GFF3 sequence ID with the FASTA sequence its genes fit in and pass the largest share
    of their checks on, at least MIN_SEARCHED_SCORE, where a random sequence doing so is unlikely
    enough for the number of sequences searched. Where several GFF3 sequence IDs fit one FASTA
    sequence, the best fitting one pairs with it. A tie is logged and goes to the first.

    Returns:
        the pairs, and {GFF3 sequence ID: FASTA sequence ID} for those whose fitting sequence was
        claimed already or pairs with a better fitting GFF3 sequence, i.e. duplicates of another
    """
    by_length = sorted((len(seq), fasta_id) for fasta_id, seq in sequences.items())
    lengths = [length for length, _ in by_length]
    hits = {}
    for gff_id, gff_sequence in gff_sequences.items():
        checks = gff_sequence.checks
        # the FASTA sequences long enough to hold its genes, found by bisection in those sorted by length
        candidates = [fasta_id for _, fasta_id in by_length[bisect.bisect_left(lengths, gff_sequence.max_end):]]
        min_passed = math.ceil(MIN_SEARCHED_SCORE * len(checks))
        if not checks or not candidates:
            continue
        if random_pass_log10(checks, min_passed) + math.log10(len(candidates)) > math.log10(MAX_RANDOM_HITS):
            continue  # too few checks to tell its sequence from a random one among this many
        # a candidate passing under half of 8 checks spread along the sequence, so from different genes
        # and not sunk by one poor or partial gene, is skipped without testing the rest
        quick = checks[::max(1, len(checks) // 8)][:8]
        found = []
        for fasta_id in candidates:
            sequence = sequences[fasta_id]
            if n_passed(sequence, quick) * 2 >= len(quick):
                passed = n_passed(sequence, checks)
                if passed >= min_passed:
                    found.append(Pair(fasta_id, 'sequence', passed, len(checks)))
        if found:
            best = max(found, key=share_passed)
            tied = [pair.fasta_id for pair in found if share_passed(pair) == share_passed(best)]
            if len(tied) > 1:
                logger.warning(f'The genes of GFF3 sequence ID {gff_id} fit the FASTA sequences {", ".join(tied)} '
                               f'equally well; paired with the first')
            hits[gff_id] = best
    # a hit on a FASTA sequence paired by ID or full name already is a duplicate; of several GFF3 IDs
    # hitting one FASTA sequence, the best fitting pairs with it and the others are duplicates
    duplicates = {gff_id: pair.fasta_id for gff_id, pair in hits.items() if pair.fasta_id in claimed}
    by_fasta_id = defaultdict(list)
    for gff_id, pair in hits.items():
        if gff_id not in duplicates:
            by_fasta_id[pair.fasta_id].append(gff_id)
    pairs = {}
    for fasta_id, gff_ids in by_fasta_id.items():
        best = max(gff_ids, key=lambda gff_id: share_passed(hits[gff_id]))
        tied = [gff_id for gff_id in gff_ids if share_passed(hits[gff_id]) == share_passed(hits[best])]
        if len(tied) > 1:
            logger.warning(f'The genes of GFF3 sequence IDs {", ".join(tied)} fit the FASTA sequence {fasta_id} '
                           f'equally well; paired the first')
        pairs[best] = hits[best]
        duplicates.update({gff_id: fasta_id for gff_id in gff_ids if gff_id != best})
    return pairs, duplicates


def pair_seqids(gff_sequences: dict[str, GffSequence], sequences: dict[str, str], full_names: dict[str, str],
                gff_name: str, fasta_name: str, max_missing_share: float = MAX_MISSING_SHARE) -> Pairing:
    """Pairs the GFF3 sequence IDs by ID or full name, then the rest by sequence. GFF3 sequences
    whose pair by ID or full name fails are taken for changed after annotation and left unpaired.

    Args:
        gff_sequences: see read_gff
        sequences: {FASTA ID: sequence}
        full_names: {FASTA ID: full name}
        gff_name: the GFF3 file name, for messages
        fasta_name: the FASTA file name, for messages
        max_missing_share: the largest share of the coding transcripts the sequences matching no
            FASTA sequence may hold

    Raises:
        MismatchedFiles: if most pairs by ID or full name fail, the changed sequences hold more than
            MAX_CHANGED_SHARE of the coding transcripts, or the sequences matching no FASTA sequence
            more than max_missing_share
    """
    pairs, changed, n_tried = pair_by_id_or_full_name(gff_sequences, sequences, full_names)
    examples = [f'{gff_id} ({fasta_id}): {why}' for gff_id, (fasta_id, why) in changed.items()]
    if len(changed) * 2 > n_tried:
        raise MismatchedFiles(
            f'The GFF3 file "{gff_name}" was made on another assembly than the FASTA file "{fasta_name}": '
            f'{len(changed)} of the {n_tried} GFF3 sequences found in it by ID or full name do not fit their '
            f'FASTA sequence, e.g. {examples[0]}. {ADVICE}')
    n_coding = sum(seq.coding_transcripts for seq in gff_sequences.values())
    n_changed_coding = sum(gff_sequences[gff_id].coding_transcripts for gff_id in changed)
    if n_changed_coding > MAX_CHANGED_SHARE * n_coding:
        raise MismatchedFiles(
            f'{len(changed)} GFF3 sequences of "{gff_name}" found in "{fasta_name}" by ID or full name do not '
            f'fit their FASTA sequence, which was likely changed after annotation (e.g. trimmed on upload to '
            f'NCBI), and hold {coding_share(n_changed_coding, n_coding)}, too many to leave out: '
            f'{"; ".join(examples[:5])}{"; ..." if len(examples) > 5 else ""}. {ADVICE}')

    unpaired = {gff_id: seq for gff_id, seq in gff_sequences.items() if gff_id not in pairs and gff_id not in changed}
    if not unpaired:
        return Pairing(pairs, changed)
    logger.info(f'Searching the FASTA sequences for {len(unpaired)} GFF3 sequence IDs not paired by ID or '
                f'full name')
    # the FASTA sequence of a changed GFF3 sequence may pair with another, whose genes fit it
    claimed = {pair.fasta_id for pair in pairs.values()}
    searched_pairs, duplicates = pair_by_sequence(unpaired, sequences, claimed)
    missing = sorted(unpaired.keys() - searched_pairs.keys() - duplicates.keys(),
                     key=lambda gff_id: -gff_sequences[gff_id].coding_transcripts)
    n_missing_coding = sum(gff_sequences[gff_id].coding_transcripts for gff_id in missing)
    if n_missing_coding > max_missing_share * n_coding:
        raise MismatchedFiles(
            f'{len(missing)} GFF3 sequences of "{gff_name}" match no sequence of "{fasta_name}" and hold '
            f'{coding_share(n_missing_coding, n_coding)}, more than the {100 * max_missing_share:g}% allowed '
            f'to leave out (--max-missing-percent), the most in '
            f'{", ".join(missing[:5])}. The FASTA file lacks sequences the GFF3 file was made on, or they were '
            f'renamed and changed after annotation (e.g. trimmed on upload to NCBI). {ADVICE}')
    return Pairing(pairs | searched_pairs, changed, duplicates)


def write_gff(gff_path: str, out_path: str, pairs: dict[str, Pair], compression: Compression | None) -> None:
    """Writes the GFF3 file with its sequence IDs renamed, leaving out the lines of unpaired ones."""
    with open_out(out_path, compression) as f_out:
        for line in GFFParser.lines(gff_path):
            # the sequence ID of a ##sequence-region or feature line renamed, the line dropped (None)
            # where the ID is unpaired; other comment lines are kept as they are
            if line.startswith('##sequence-region'):
                _, seqid, rest = line.split(maxsplit=2)
                line = f'##sequence-region {pairs[seqid].fasta_id} {rest}' if seqid in pairs else None
            elif line.strip() and not line.startswith('#'):
                seqid, rest = line.split('\t', 1)
                line = f'{pairs[seqid].fasta_id}\t{rest}' if seqid in pairs else None
            if line is not None:
                f_out.write(line)


def fasta_out_path(fasta_path: str, gff_out_path: str) -> str:
    """The FASTA copy's path beside the output GFF3 file: the FASTA file's name with
    _sequences_dropped before its extension, e.g. genome_sequences_dropped.fna.gz."""
    return os.path.join(os.path.dirname(gff_out_path),
                        insert_before_extension(os.path.basename(fasta_path), '_sequences_dropped'))


def coding_share(n: int, total: int) -> str:
    return f'{n} of {total} coding transcripts ({100 * n / total if total else 0:.2f}%)'


def log_summary(gff_sequences: dict[str, GffSequence], pairing: Pairing, gff_name: str, fasta_name: str) -> None:
    pairs, changed, duplicates = pairing.pairs, pairing.changed, pairing.duplicates
    methods = defaultdict(int)
    for pair in pairs.values():
        methods[pair.method] += 1
    logger.info(f'Paired {len(pairs)} of {len(gff_sequences)} GFF3 sequence IDs of "{gff_name}" with a FASTA '
                f'sequence of "{fasta_name}": {methods["id"]} by identical ID, {methods["full_name"]} by a word of '
                f'the FASTA full name (header line), {methods["sequence"]} by sequence alone')
    unchecked = [gff_id for gff_id, pair in pairs.items() if not pair.checks]
    if unchecked:
        logger.info(f'  {len(unchecked)} of them hold no CDS or intron to check, paired on ID or full name and '
                    f'length alone')

    n_coding = sum(seq.coding_transcripts for seq in gff_sequences.values())

    def log_left_out(gff_ids: list[str], why: str, examples: list[str]) -> int:
        n_left_out = sum(gff_sequences[gff_id].coding_transcripts for gff_id in gff_ids)
        if gff_ids:
            n_lines = sum(gff_sequences[gff_id].n_lines for gff_id in gff_ids)
            logger.warning(f'Left out {len(gff_ids)} GFF3 sequence IDs {why}: {n_lines} lines, '
                           f'{coding_share(n_left_out, n_coding)}, {examples}')
        return n_left_out

    changed_ids = sorted(changed)
    n_lost = log_left_out(changed_ids, 'found in the FASTA file by ID or full name whose FASTA sequence does not fit '
                                       'their genes, likely changed after annotation (e.g. trimmed on upload to '
                                       'NCBI); their FASTA sequences are left out as well, unless paired with '
                                       'another GFF3 sequence ID',
                          'each: ' + '; '.join(f'{gff_id} ({changed[gff_id][0]}): {changed[gff_id][1]}'
                                               for gff_id in changed_ids))
    gff_id_by_fasta_id = {pair.fasta_id: gff_id for gff_id, pair in pairs.items()}
    duplicate_ids = sorted(duplicates)
    n_lost += log_left_out(duplicate_ids, 'whose genes fit a FASTA sequence paired with another GFF3 sequence ID, '
                                          'likely duplicate sequences removed from the FASTA file',
                           'e.g. ' + ', '.join(f'{gff_id} (like {gff_id_by_fasta_id[duplicates[gff_id]]})'
                                               for gff_id in duplicate_ids[:5]))
    unpaired = sorted(set(gff_sequences) - pairs.keys() - changed.keys() - duplicates.keys())
    n_lost += log_left_out(unpaired, 'matching no FASTA sequence', 'e.g. ' + ', '.join(unpaired[:5]))
    logger.info(f'Kept {coding_share(n_coding - n_lost, n_coding)}')


def log_stop_codons(gff_sequences: dict[str, GffSequence], pairs: dict[str, Pair], sequences: dict[str, str]) -> None:
    """Logs the paired transcripts whose stop codon lies right after their CDS instead of in it."""
    outside, unextendable = 0, []
    for gff_id, pair in pairs.items():
        sequence = sequences[pair.fasta_id]
        for stop in gff_sequences[gff_id].stops:
            if (codon_at(sequence, stop.codon, stop.is_plus) not in MOTIFS['stop']
                    and codon_at(sequence, stop.after, stop.is_plus) in MOTIFS['stop']):
                outside += 1
                if not stop.exon_beyond:
                    unextendable.append(stop.transcript_id)
    if outside:
        logger.info(f'{outside} transcripts have their stop codon right after their CDS instead of as its last '
                    f'codon, as GFF3 converted from GTF does (GTF writes the stop codon as a line of its own, '
                    f'outside the CDS); on import, GeenuFF extends the CDS of such a transcript by the 3 bases '
                    f'following it along its exons')
    if unextendable:
        logger.warning(f'{len(unextendable)} of them have no exon reaching past their CDS (CDS lines only, or exon '
                       f'lines ending with the CDS), so GeenuFF takes the 3 bases right after their last exon and '
                       f'extends that exon and the transcript to hold the stop codon; check that this suits the '
                       f'GFF3 file: {", ".join(unextendable[:5])}{", ..." if len(unextendable) > 5 else ""}')


def write_pairs(path: str, pairs: dict[str, Pair], full_names: dict[str, str]) -> None:
    with open(path, 'w') as f:
        f.write('gff_seqid\tfasta_id\tpaired_by\tchecks_passed\tchecks\tfasta_full_name\n')
        for gff_id, pair in sorted(pairs.items()):
            f.write(f'{gff_id}\t{pair.fasta_id}\t{pair.method}\t{pair.passed}\t{pair.checks}\t'
                    f'{full_names[pair.fasta_id]}\n')


def rename(fasta_path: str, gff_path: str, out_path: str, pairs_out_path: str | None = None,
           max_missing_share: float = MAX_MISSING_SHARE, gzipped: bool | None = None,
           compression_level: int = DEFAULT_COMPRESSION_LEVEL,
           threads: int = DEFAULT_THREADS) -> tuple[dict[str, GffSequence], Pairing, bool]:
    """Pairs the sequence IDs of a GFF3 file with those of its FASTA file and writes the renamed GFF3
    file, with a FASTA copy where changed sequences are left out (see the module docstring).

    Args:
        fasta_path: the FASTA file
        gff_path: the GFF3 file
        out_path: the renamed GFF3 file to write
        pairs_out_path: a TSV file to write the pairs to, if given
        max_missing_share: see pair_seqids
        gzipped: whether to gzip the renamed GFF3 file and the FASTA copy, None for each as its input
        compression_level: their gzip level
        threads: the threads compressing them (see compression_for)

    Returns:
        the GFF3 sequences (see read_gff), their pairing, and whether a file was written, none being
        where all IDs match already

    Raises:
        MismatchedFiles: see pair_seqids
    """
    gff_name, fasta_name = os.path.basename(gff_path), os.path.basename(fasta_path)
    logger.info(f'Pairing the sequence IDs of the GFF3 file "{gff_name}" with those of the FASTA file '
                f'"{fasta_name}"')
    logger.info(f'Reading "{fasta_path}"')
    sequences, full_names = {}, {}
    for full_name, sequence in FastaParser.read_fasta(fasta_path):
        fasta_id = FastaParser.header_id(full_name)
        sequences[fasta_id], full_names[fasta_id] = sequence.upper(), full_name
    logger.info(f'Reading "{gff_path}"')
    gff_sequences = read_gff(gff_path)

    pairing = pair_seqids(gff_sequences, sequences, full_names, gff_name, fasta_name, max_missing_share)
    log_stop_codons(gff_sequences, pairing.pairs, sequences)
    pairs = pairing.pairs
    if len(pairs) == len(gff_sequences) and all(pair.method == 'id' for pair in pairs.values()):
        logger.info(f'All {len(pairs)} sequence IDs of "{gff_name}" match "{fasta_name}" and their genes fit '
                    f'their FASTA sequences: nothing to rename, no file written')
        return gff_sequences, pairing, False
    log_summary(gff_sequences, pairing, gff_name, fasta_name)

    compression = compression_for(gff_path, gzipped, compression_level, threads)
    write_gff(gff_path, output_path(out_path, compression), pairs, compression)
    # a FASTA sequence paired with one GFF3 sequence ID stays, whichever other one it did not fit
    dropped = {fasta_id for fasta_id, _ in pairing.changed.values()} - {pair.fasta_id for pair in pairs.values()}
    if dropped:
        compression = compression_for(fasta_path, gzipped, compression_level, threads)
        fasta_out = output_path(fasta_out_path(fasta_path, out_path), compression)
        write_fasta_without(fasta_path, fasta_out, dropped, compression)
        logger.warning(f'Wrote the FASTA file without the {len(dropped)} changed sequences to "{fasta_out}"; '
                       f'import it with the renamed GFF3 file')
    if pairs_out_path:
        write_pairs(pairs_out_path, pairs, full_names)
    return gff_sequences, pairing, True


def main(args: argparse.Namespace) -> None:
    setup_logging(args.log_file)
    try:
        rename(args.fasta, args.gff, args.out, args.pairs_out, args.max_missing_percent / 100, args.gzip,
               args.compression_level, args.threads)
    except MismatchedFiles as e:
        logger.error(e)
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--fasta', type=str, required=True, help='FASTA file the GFF3 file annotates (may be gzipped)')
    parser.add_argument('--gff', type=str, required=True, help='GFF3 file to rename (may be gzipped)')
    parser.add_argument('-o', '--out', type=str,
                        help=(f'output GFF3 file, gzipped where the GFF3 input is (see --gzip) (default: the GFF3 '
                              f'file\'s name with {OUT_SUFFIX} before the extension, e.g. genes{OUT_SUFFIX}.gff3); '
                              f'not written when all sequence IDs match already. Where sequences changed after '
                              f'annotation are left out, the FASTA file is copied without them beside it, named like '
                              f'the FASTA file with _sequences_dropped before its extension and gzipped where the '
                              f'FASTA input is'))
    parser.add_argument('--pairs-out', type=str,
                        help=('TSV file listing each pair made, how, how many of its checks passed and the '
                              'full name of the FASTA sequence'))
    parser.add_argument('--max-missing-percent', type=float, default=100 * MAX_MISSING_SHARE,
                        help=('the most coding transcripts, in percent, the GFF3 sequences matching no FASTA '
                              'sequence may hold to be left out; above it, nothing is written '
                              f'(default: {100 * MAX_MISSING_SHARE:g})'))
    add_output_compression_args(parser)
    parser.add_argument('--log-file', type=str, help='file to write the log to, in addition to the terminal')
    args = parser.parse_args()
    require_input_files(parser, [args.gff, args.fasta])
    require_fasta_and_gff(parser, args.fasta, args.gff)
    default_out = insert_before_extension(args.gff, OUT_SUFFIX)
    args.out = output_path(args.out or default_out, output_compression(args, args.gff))
    check_output_paths(parser, [args.out, args.pairs_out, args.log_file], [args.gff, args.fasta])
    main(args)
