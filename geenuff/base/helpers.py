import argparse
import copy
import hashlib
import itertools
import logging
import enum
import os
import sys
from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from types import GeneratorType
from typing import TYPE_CHECKING, BinaryIO, TextIO

from zlib_ng import gzip_ng_threaded

from geenuff.parsing.fastahelper import FastaParser, detect_format, text_or_gzip_open
from geenuff.parsing.gffhelper import GFFObject, GFFParser

if TYPE_CHECKING:
    from geenuff.base.orm import Coordinate

logger = logging.getLogger(__name__)


# Types related
# ----------------------------------


# Enum makers
def join_to_enum(name, *args):
    """joins enums from args into returned enum"""
    enum_bits = []
    for cls in args:
        enum_bits += list(cls)
    out = enum.Enum(name, [(x.name, x.value) for x in enum_bits])
    return out


def make_enum(name, *args):
    """makes enum from list of strings"""
    return enum.Enum(name, [(x, x) for x in args])


# General
# ----------------------------------
# bases encoded and hashed at a time by sequence_hash; slices up to 64 KiB hashed 500 Mb in 0.30 s,
# from 256 KiB on in 0.66-0.79 s, as larger slices no longer stay in the CPU cache between encoding
# and hashing; 64 KiB is the largest of the fast sizes, keeping the loop short (one machine)
HASH_SLICE = 1 << 16


def sequence_hash(sequence: str) -> str:
    """The SHA-1 hex digest of a sequence, hashed in slices so no encoded copy of it all is held."""
    sha1 = hashlib.sha1()
    for start in range(0, len(sequence), HASH_SLICE):
        sha1.update(sequence[start:start + HASH_SLICE].encode())
    return sha1.hexdigest()


def in_enum_values(x, enum):
    return x in [item.value for item in enum]


def none_to_list(x):
    if x is None:
        return []
    else:
        assert isinstance(x, list)
        return x


def convert2list(obj):
    if isinstance(obj, list):
        out = obj
    elif isinstance(obj, set) or isinstance(obj, GeneratorType) or isinstance(obj, tuple):
        out = list(obj)
    else:
        out = [obj]
    return out


def full_db_path(path):
    sqlite_prefix = 'sqlite:///'
    if not path.startswith(sqlite_prefix):
        return sqlite_prefix + path
    return path


def db_attr_as_dict(orm_obj):
    """removes known, shared non-db entry attributes from a copy of orm_obj.__dict__
    can be used, e.g. to set up new db entry matching the here-by filtered one
    """
    exclude = ['_sa_instance_state', 'handler']
    out = copy.copy(orm_obj.__dict__)
    for item in exclude:
        if item in out:
            del out[item]
    return out


def chunk_str(string, length):
    for i in range(0, len(string), length):
        yield string[i:(i + length)]


def get_repr(class_name, params, addition=''):
    param_str = ', '.join(f'{k}:{v}' for k, v in params.items())
    if addition:
        return class_name + '[' + param_str + ', ' + addition + ']'
    else:
        return class_name + '[' + param_str + ']'


def insert_before_extension(path: str, text: str) -> str:
    """The path with text inserted before its extension, a final .gz or .zip kept last, e.g.
    genome.fna.gz and _length_filtered give genome_length_filtered.fna.gz."""
    compressed = next((suffix for suffix in ('.gz', '.zip') if path.endswith(suffix)), '')
    stem, extension = os.path.splitext(path[:len(path) - len(compressed)])
    return f'{stem}{text}{extension}{compressed}'


def setup_logging(log_file: str | None = None) -> None:
    """Logs INFO and above to stdout and, given a path, to that file as well, overwriting it.

    Args:
        log_file: the file to write the log to besides stdout, if given
    """
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file, mode='w'))
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s', handlers=handlers)


# Mapper
# ----------------------------------
class NonMatchableIDs(Exception):
    pass


def _listed(ids: list[str], n_shown: int = 5) -> str:
    return ', '.join(ids[:n_shown]) + (', ...' if len(ids) > n_shown else '')


def match_seqids(fasta_ids: Iterable[str], gff_ids: Iterable[str]) -> dict[str, str]:
    """Pairs each GFF sequence ID with a FASTA sequence ID: the identical one if there is one, else
    the only one of the remaining FASTA IDs containing it (e.g. 'Chr1.v2' for 'Chr1').

    Args:
        fasta_ids: the FASTA sequence IDs, i.e. the headers up to their first whitespace
        gff_ids: the sequence IDs of the GFF feature lines

    Returns:
        {GFF ID: FASTA ID}

    Raises:
        NonMatchableIDs: if a GFF ID matches no FASTA ID or several, or two GFF IDs the same one.
    """
    gff_ids = set(gff_ids)
    matched = {gff_id: gff_id for gff_id in gff_ids & set(fasta_ids)}
    unclaimed = set(fasta_ids) - matched.keys()
    substring_matched, unmatched, ambiguous = {}, [], []
    for gff_id in sorted(gff_ids - matched.keys()):
        candidates = sorted(fa_id for fa_id in unclaimed if gff_id in fa_id)
        if len(candidates) == 1:
            substring_matched[gff_id] = candidates[0]
        elif candidates:
            ambiguous.append(f'{gff_id} (in {_listed(candidates)})')
        else:
            unmatched.append(gff_id)
    gff_ids_by_fasta_id = defaultdict(list)
    for gff_id, fa_id in substring_matched.items():
        gff_ids_by_fasta_id[fa_id].append(gff_id)
    claimed_twice = [f'{fa_id} ({", ".join(ids)})' for fa_id, ids in sorted(gff_ids_by_fasta_id.items())
                     if len(ids) > 1]

    problems = [f'{label} ({len(ids)}): {_listed(ids)}'
                for label, ids in (('GFF sequence IDs matching no FASTA header', unmatched),
                                   ('GFF sequence IDs contained in several FASTA headers', ambiguous),
                                   ('FASTA headers containing several GFF sequence IDs', claimed_twice))
                if ids]
    if problems:
        raise NonMatchableIDs('; '.join(problems))
    if substring_matched:
        pairs = [f'{gff_id} -> {fa_id}' for gff_id, fa_id in substring_matched.items()]
        logger.warning(f'GFF sequence IDs paired with the one FASTA header containing them, for lack of an '
                       f'identical one ({len(pairs)}); check that these pairs are right: {_listed(pairs)}')
    return matched | substring_matched


# Script options and file checks
# ----------------------------------
# gzip level of the scripts' output files; 4 compresses a genome nearly as well as gzip's default 6,
# several times faster
DEFAULT_COMPRESSION_LEVEL = 4
# threads compressing the scripts' gzipped output files, each about as fast as plain gzip
DEFAULT_THREADS = 4


@dataclass(frozen=True)
class Compression:
    """How a gzipped output file is written: its gzip level and the threads compressing it."""
    level: int
    threads: int


def add_output_compression_args(parser: argparse.ArgumentParser) -> None:
    """Adds --gzip/--no-gzip, --compression-level and --threads, read by output_compression, to a
    script's arguments."""
    parser.add_argument('--gzip', action=argparse.BooleanOptionalAction,
                        help=('gzip the output files, or write them as plain text; by default each output file is '
                              'gzipped where its input file is, its name given .gz or not to match'))
    parser.add_argument('--compression-level', type=int, default=DEFAULT_COMPRESSION_LEVEL, choices=range(1, 10),
                        metavar='{1-9}', help=(f'gzip level of the output files, 1 fastest, 9 smallest (default: '
                                               f'{DEFAULT_COMPRESSION_LEVEL})'))
    parser.add_argument('--threads', type=int, default=DEFAULT_THREADS,
                        help=(f'threads compressing the gzipped output files; 1 for a single one, -1 for all cores '
                              f'the system has (default: {DEFAULT_THREADS})'))


def compression_for(input_path: str, gzipped: bool | None, compression_level: int,
                    threads: int = DEFAULT_THREADS) -> Compression | None:
    """How an output file is gzipped, None for plain text: gzipped as given, else where its input file
    is gzipped or zipped."""
    if gzipped is None:
        gzipped = detect_format(input_path) != 'text'
    return Compression(compression_level, threads) if gzipped else None


def output_compression(args: argparse.Namespace, input_path: str) -> Compression | None:
    """compression_for with a script's --gzip, --compression-level and --threads (see
    add_output_compression_args)."""
    return compression_for(input_path, args.gzip, args.compression_level, args.threads)


def output_path(path: str, compression: Compression | None) -> str:
    """The path ending in .gz for gzipped output, without it for plain text, a final .zip dropped
    either way."""
    stem = path.removesuffix('.gz').removesuffix('.zip')
    return stem if compression is None else f'{stem}.gz'


def check_output_paths(parser: argparse.ArgumentParser, outputs: Iterable[str | None],
                       inputs: Iterable[str | None]) -> None:
    """Stops with a usage error, before anything is written, where an output's directory does not
    exist, as for a mistyped path, or where an output names an input file, which writing it would
    empty before it is read. Paths given as None, i.e. options left out, are skipped."""
    input_paths = {os.path.realpath(path.removeprefix('sqlite:///')) for path in inputs if path}
    for path in outputs:
        if not path:
            continue
        directory = os.path.dirname(path.removeprefix('sqlite:///')) or '.'
        if not os.path.isdir(directory):
            parser.error(f'the directory "{directory}" of the output "{path}" does not exist')
        if os.path.realpath(path.removeprefix('sqlite:///')) in input_paths:
            parser.error(f'the output "{path}" is an input file, which writing it would destroy; give another path')


def require_input_files(parser: argparse.ArgumentParser, paths: Iterable[str | None]) -> None:
    """Stops with a usage error where an input file does not exist, before any output is created (a
    database opened at a mistyped path would be created empty). Paths given as None are skipped."""
    for path in paths:
        if path and not os.path.isfile(path.removeprefix('sqlite:///')):
            parser.error(f'the input file "{path}" does not exist')


# the endings of FASTA and GFF3 file names, each possibly followed by .gz or .zip
FASTA_EXTENSIONS = ('.fa', '.fasta', '.fna', '.fas')
GFF_EXTENSIONS = ('.gff', '.gff3')


def is_fasta_name(path: str) -> bool:
    return path.lower().removesuffix('.gz').removesuffix('.zip').endswith(FASTA_EXTENSIONS)


def is_gff_name(path: str) -> bool:
    return path.lower().removesuffix('.gz').removesuffix('.zip').endswith(GFF_EXTENSIONS)


def require_fasta_and_gff(parser: argparse.ArgumentParser, fasta: str | None, gff: str | None) -> None:
    """Stops with a usage error where the FASTA or GFF3 file name ends like the other kind's, as when
    the two are swapped, where the FASTA file does not start with '>', or where the GFF3 file's first
    line that is no comment does not have the 9 tab-separated columns of a feature line, e.g. 10 from
    a sequence ID holding a tab. Either file may be None, for a script reading only the other."""
    fasta_named_gff = fasta and is_gff_name(fasta)
    gff_named_fasta = gff and is_fasta_name(gff)
    if fasta_named_gff and gff_named_fasta:
        parser.error(f'the FASTA file "{fasta}" and the GFF3 file "{gff}" look swapped by their names')
    if fasta_named_gff:
        parser.error(f'"{fasta}", given as the FASTA file, is named like a GFF3 file')
    if gff_named_fasta:
        parser.error(f'"{gff}", given as the GFF3 file, is named like a FASTA file')
    if fasta:
        with text_or_gzip_open(fasta, 'rb') as f:
            if f.read(1) != b'>':
                parser.error(f'"{fasta}" is no FASTA file: it does not start with ">"')
    if gff:
        with text_or_gzip_open(gff, 'rb') as f:
            line = next((line for line in f if line.strip() and not line.startswith(b'#')), b'')
        n_columns = len(line.rstrip(b'\r\n').split(b'\t'))
        if n_columns != 9:
            parser.error(f'"{gff}" is no valid GFF3 file: its first feature line has {n_columns} tab-separated '
                         f'columns instead of 9: {line.decode(errors="replace").strip()[:200]}')


def open_out(path: str, compression: Compression | None, binary: bool = False) -> TextIO | BinaryIO:
    """Opens a file for writing, gzipped as compression says (with zlib-ng, in threads), or as plain
    text where it is None."""
    if compression is None:
        return open(path, 'wb' if binary else 'w')
    return gzip_ng_threaded.open(path, 'wb' if binary else 'wt', compresslevel=compression.level,
                                 threads=compression.threads)


def write_fasta_without(fasta_path: str, out_path: str, left_out: set[str], compression: Compression | None) -> None:
    """Copies a FASTA file byte for byte, leaving out the sequences of the given IDs."""
    keep = True
    with open_out(out_path, compression, binary=True) as f_out:
        for is_header, piece in FastaParser.pieces(fasta_path):
            if is_header:
                keep = FastaParser.header_id(piece[1:].decode()) not in left_out
            if keep:
                f_out.write(piece)


def _ids_on_removed(gff_path: str, removed: set[str]) -> set[str]:
    """Returns the IDs of the GFF3 lines on the removed sequences, and those whose lines are on kept
    sequences but their parents all lie on removed ones, i.e. of genes spread over a removed and a
    kept sequence."""
    on_removed, on_kept, parents_of = set(), set(), {}
    for line in GFFParser.feature_lines(gff_path):
        columns = line.split('\t')
        attributes = GFFObject.parse_attributes(columns[8])
        line_id = attributes.get('ID')
        if line_id is None:
            continue  # parent line without ID, children can have no ID and still be removed
        if columns[0] in removed:
            on_removed.add(line_id)
        else:
            on_kept.add(line_id)
            if parents := GFFObject.parents(attributes):
                parents_of[line_id] = parents
    dropped = on_removed - on_kept
    # one round per level of the Parent chain, e.g. mRNA lines below a gene, then exons below them
    while below := {line_id for line_id, parents in parents_of.items()
                    if line_id not in dropped and all(parent in dropped for parent in parents)}:
        dropped |= below
    return dropped


def write_gff_without(gff_path: str, out_path: str, removed: set[str], compression: Compression | None,
                      renamed: dict[str, str] | None = None) -> int:
    """Copies a GFF3 file without the feature and ##sequence-region lines of the removed sequences,
    and without the lines on kept sequences all of whose parents lie on removed ones, through any
    number of Parent steps, which the import would leave out as on another sequence than their parent
    (see _ids_on_removed); a line naming only some of those as Parent keeps the others. The kept
    sequence IDs in renamed get its new names. Logs and returns the number of lines left out on kept
    sequences."""
    renamed = renamed or {}
    dropped_ids = _ids_on_removed(gff_path, removed)
    n_below = 0
    with open_out(out_path, compression) as f_out:
        for line in GFFParser.lines(gff_path):
            if line.startswith('##sequence-region'):
                _, seqid, rest = line.split(maxsplit=2)
                if seqid in removed:
                    continue
                if seqid in renamed:
                    line = f'##sequence-region {renamed[seqid]} {rest}'
            elif line.strip() and not line.startswith('#'):
                columns = line.rstrip('\n').split('\t')
                if columns[0] in removed:
                    continue
                parents = GFFObject.parents(GFFObject.parse_attributes(columns[8]))
                kept = [parent for parent in parents if parent not in dropped_ids]
                if parents and not kept:
                    n_below += 1
                    continue
                if len(kept) < len(parents):
                    columns[8] = GFFObject.with_parents(columns[8], kept)
                columns[0] = renamed.get(columns[0], columns[0])
                line = '\t'.join(columns) + '\n'
            f_out.write(line)
    if n_below:
        logger.info(f'Left out {n_below} lines on kept sequences whose parents all lie on removed ones, i.e. of '
                    f'genes spread over several sequences, which the import leaves out unmasked as well '
                    f'(i.e. the GeenuFF reason "on another sequence than their parent")')
    return n_below


# GFF start/end to GeenuFF
# ----------------------------------
def strand_as_bool(strand):
    # anything but '+' or '-' raises; use strand_or_none below where that is not fatal
    if strand == '+':
        return True
    elif strand == '-':
        return False
    else:
        raise ValueError(f'cannot interpret strand "{strand}"')


def get_strand_direction(gffentry):
    return strand_as_bool(gffentry.strand)


def strand_or_none(gffentry):
    """The entry's strand, or None where the GFF gives neither '+' nor '-', e.g. NCBI's '?' for
    a trans-spliced feature."""
    try:
        return get_strand_direction(gffentry)
    except ValueError:
        return None


# todo: could this be simpler?
def get_geenuff_start_end(gff_start, gff_end, is_plus_strand):
    gff_start, gff_end = to_count_from_0(gff_start), to_count_from_0(gff_end)

    if is_plus_strand:
        start = gff_start
        end = to_exclusive_end(gff_end, is_plus_strand)
    else:
        start = gff_end
        end = to_exclusive_end(gff_start, is_plus_strand)
    return start, end


def to_count_from_0(coord):
    return coord - 1


def to_exclusive_end(end, is_plus_strand):
    if is_plus_strand:
        return end + 1
    else:
        return end - 1


def geenuff_to_gff_start_end(start: int, end: int, is_plus_strand: bool) -> tuple[int, int]:
    """inverse of get_geenuff_start_end: turns a GeenuFF-internal (0-based, direction-aware,
    half-open) start/end back into a standard GFF3 1-based, inclusive, ascending start/end"""
    if is_plus_strand:
        return start + 1, end
    else:
        return end + 2, start + 1


def is_masked_whole(transcript_feature, masks) -> bool:
    """Whether one of the masks covers the whole transcript feature, leaving no labels from it.
    The masks of a transcript are merged, so one covering it suffices, and any other lies beside
    it. Takes anything with a GeenuFF start and end, on either strand."""
    t_lo, t_hi = sorted((transcript_feature.start, transcript_feature.end))
    return any(min(m.start, m.end) <= t_lo and max(m.start, m.end) >= t_hi for m in masks)


# maps "normal" (0, 1, 2; cumulative CDS length so far, mod 3) phase counting to GFF3's
# biological (0, 2, 1) convention. Self-inverse (0<->0, 1<->2, 2<->1), so the same mapping
# converts a phase in either direction
GFF_PHASE_FROM_NORMAL: dict[int, int] = {0: 0, 1: 2, 2: 1}


# Reverse complement
# ----------------------------------

def mk_rc_key():
    fw = 'ACGTMRWSYKVHDBN'
    rv = 'TGCAKYWSRMBDHVN'
    key = {}
    for f, r in zip(fw, rv):
        key[f] = r
    return key


# so one doesn't recalculate it for every call of revers_complement
REV_COMPLEMENT_KEY = mk_rc_key()


def reverse_complement(seq):
    key = REV_COMPLEMENT_KEY
    rc_seq = []
    for base in reversed(seq):
        try:
            rc_seq.append(key[base])
        except KeyError as e:
            raise KeyError(f'{e} caused by non DNA character {base}')
    return rc_seq


# Start/Stop codon detection
# ----------------------------------
START_CODON = 'ATG'
STOP_CODONS = ['TAG', 'TGA', 'TAA']


def spliced_cds_sequence(coordinate: 'Coordinate', cds_pieces: Sequence[GFFObject], is_plus_strand: bool) -> str:
    """joins raw GFF CDS pieces (each with a 1-based, inclusive .start/.end, sorted by
    ascending genomic start) into one continuous, 5' to 3' oriented sequence"""
    ordered = cds_pieces if is_plus_strand else reversed(cds_pieces)
    parts = []
    for piece in ordered:
        start, end = get_geenuff_start_end(piece.start, piece.end, is_plus_strand)
        part = (coordinate.sequence_range(start, end) if is_plus_strand
                else coordinate.sequence_range(end + 1, start + 1))
        parts.append(part if is_plus_strand else ''.join(reverse_complement(part)))
    return ''.join(parts)


def has_inframe_stop_codon(cds_seq: str) -> bool:
    """checks a fully spliced, 5' to 3' oriented CDS sequence for a stop codon before the
    final codon; a truncated (not a multiple of 3) trailing partial codon is ignored"""
    last_full_codon_start = len(cds_seq) - len(cds_seq) % 3 - 3
    codons = (cds_seq[i:i + 3] for i in range(0, last_full_codon_start, 3))
    return any(codon in STOP_CODONS for codon in codons)


def codon_after_cds(cds: Sequence[tuple[int, int]], exons: Sequence[tuple[int, int]],
                    is_plus_strand: bool) -> tuple[list[int], int]:
    """The 1-based positions of the 3 bases following a CDS in reading order, along its exons, across
    an intron if need be, and straight on past the last exon where fewer are left in them, with how
    many of the 3 lie past the exons; where a stop codon GTF leaves out of the CDS lies. The positions
    can repeat, for overlapping exons, or lie past the sequence end."""
    if is_plus_strand:
        cds_end = max(end for _, end in cds)
        in_exons = (p for start, end in sorted(exons) for p in range(max(start, cds_end + 1), end + 1))
        past = itertools.count(max(end for _, end in [*exons, *cds]) + 1)
    else:
        cds_start = min(start for start, _ in cds)
        in_exons = (p for start, end in sorted(exons, reverse=True)
                    for p in range(min(end, cds_start - 1), start - 1, -1))
        past = itertools.count(min(start for start, _ in [*exons, *cds]) - 1, -1)
    following = list(itertools.islice(in_exons, 3))
    n_past = 3 - len(following)
    return following + list(itertools.islice(past, n_past)), n_past


# Span/range helpers
# ----------------------------------
def merge_ranges(ranges: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merges ascending, half-open (lo, hi) ranges into a minimal, sorted list of ranges."""
    merged = []
    for lo, hi in sorted(ranges):
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return [tuple(r) for r in merged]


def merge_gff_spans(spans: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """GFF spans (1-based, inclusive) sorted and merged where they touch or overlap."""
    merged = []
    for start, end in sorted(spans):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


# CDS position check
# ----------------------------------
def cds_outside_exons(cds: Iterable[tuple[int, int]], exons: Iterable[tuple[int, int]]) -> bool:
    """Whether a CDS span reaches outside the exon spans (GFF, 1-based, inclusive), touching or
    overlapping exons counting as one."""
    covered = merge_gff_spans(exons)
    return any(not any(lo <= start and end <= hi for lo, hi in covered) for start, end in cds)


# String/import summary formatting
# ----------------------------------
def plural(n: int, noun: str) -> str:
    """'1 transcript', '2 transcripts'."""
    return f'{n} {noun}{"" if n == 1 else "s"}'


# SQL alchemy core queue control
# ----------------------------------
class Counter(object):
    def __init__(self, cl=None, at=0):
        self.cl = cl
        self._at = at

    def __call__(self, *args, **kwargs):
        self._at += 1
        return self._at


class QueueController(object):
    def __init__(self, session, engine):
        self.session = session
        self.engine = engine
        self.ordered_queues = []

    def total_size(self):
        return sum([len(q.queue) for q in self.ordered_queues])

    def execute_so_far(self):
        conn = self.engine.connect()
        n_elements = self.total_size()
        for queue in self.ordered_queues:
            if queue.queue:
                conn.execute(queue.action, queue.queue)
                del queue.queue[:]
        self.session.commit()
        logger.debug(f'All core queues executed ({n_elements} elements)')


class CoreQueue(object):
    def __init__(self, action):
        self.queue = []
        self.action = action
