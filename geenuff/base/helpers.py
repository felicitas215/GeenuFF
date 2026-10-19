import argparse
import copy
import gzip
import hashlib
import logging
import enum
import os
import sys
from collections import defaultdict
from collections.abc import Iterable, Iterator, Sequence
from types import GeneratorType
from typing import BinaryIO, TextIO

import dustdas.fastahelper
from dustdas.gffhelper import GFFObject

logger = logging.getLogger(__name__)


##### types related #######


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


##### General #####
def sequence_hash(sequence):
    sha1 = hashlib.sha1(sequence.encode())
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
    exclude = ["_sa_instance_state", "handler"]
    out = copy.copy(orm_obj.__dict__)
    for item in exclude:
        if item in out:
            del out[item]
    return out


def chunk_str(string, length):
    for i in range(0, len(string), length):
        yield string[i:(i + length)]


def get_repr(class_name, params, addition=''):
    param_str = ', '.join('{}:{}'.format(k, v) for k, v in params.items())
    if addition:
        return class_name + '[' + param_str + ', ' + addition + ']'
    else:
        return class_name + '[' + param_str + ']'


def insert_before_extension(path: str, text: str) -> str:
    """The path with text inserted before its extension, a final .gz kept last, e.g.
    genome.fna.gz and _length_filtered give genome_length_filtered.fna.gz."""
    gz = '.gz' if path.endswith('.gz') else ''
    stem, extension = os.path.splitext(path[:len(path) - len(gz)])
    return f'{stem}{text}{extension}{gz}'


def setup_logging(log_file: str | None = None) -> None:
    """Logs INFO and above to stdout and, given a path, to that file as well, overwriting it.

    Args:
        log_file: the file to write the log to besides stdout, if given
    """
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file, mode='w'))
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s', handlers=handlers)


##### Mapper #####

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


def gff_feature_lines(gff_file: str) -> Iterator[str]:
    """Yields the feature lines of a (gzipped) GFF3 file, skipping comments and blank lines, up to an
    embedded FASTA section (a ##FASTA line or a line starting with '>')."""
    with dustdas.fastahelper.text_or_gzip_open(gff_file) as f:
        for line in f:
            if line.startswith(('##FASTA', '>')):
                return
            if line.strip() and not line.startswith('#'):
                yield line


def gff_attributes(column: str) -> dict[str, str]:
    """{tag: value} of a GFF3 attribute column, values left unsplit (e.g. 'a,b' for several parents)."""
    return dict(pair.split('=', 1) for pair in column.strip().strip(';').split(';') if '=' in pair)


# gzip level of the scripts' output files; 4 compresses a genome nearly as well as gzip's default 6,
# several times faster
DEFAULT_COMPRESSION_LEVEL = 4


def add_output_compression_args(parser: argparse.ArgumentParser) -> None:
    """Adds --gzip/--no-gzip and --compression-level, read by output_compression, to a script's arguments."""
    parser.add_argument('--gzip', action=argparse.BooleanOptionalAction,
                        help=('gzip the output files, or write them as plain text; by default each output file is '
                              'gzipped where its input file is, its name given .gz or not to match'))
    parser.add_argument('--compression-level', type=int, default=DEFAULT_COMPRESSION_LEVEL, choices=range(1, 10),
                        metavar='{1-9}', help=(f'gzip level of the output files, 1 fastest, 9 smallest (default: '
                                               f'{DEFAULT_COMPRESSION_LEVEL})'))


def compression_for(input_path: str, gzipped: bool | None, compression_level: int) -> int | None:
    """The gzip level of an output file, None for plain text: gzipped as given, else as its input file."""
    if gzipped is None:
        gzipped = is_gzipped(input_path)
    return compression_level if gzipped else None


def output_compression(args: argparse.Namespace, input_path: str) -> int | None:
    """compression_for with a script's --gzip and --compression-level (see add_output_compression_args)."""
    return compression_for(input_path, args.gzip, args.compression_level)


def output_path(path: str, compression_level: int | None) -> str:
    """The path ending in .gz for gzipped output, without it for plain text."""
    stem = path.removesuffix('.gz')
    return stem if compression_level is None else f'{stem}.gz'


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


def open_out(path: str, compression_level: int | None, binary: bool = False) -> TextIO | BinaryIO:
    """Opens a file for writing, gzipped at compression_level, or as plain text where it is None."""
    if compression_level is None:
        return open(path, 'wb' if binary else 'w')
    return gzip.open(path, 'wb' if binary else 'wt', compresslevel=compression_level)


FASTA_CHUNK_SIZE = 1 << 20


def is_gzipped(path: str) -> bool:
    """Whether a file is gzipped, by its first two bytes rather than its name."""
    with open(path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'


def open_binary(path: str) -> BinaryIO:
    """Opens a file, gzipped or not, for reading bytes."""
    return gzip.open(path, 'rb') if is_gzipped(path) else open(path, 'rb')


def fasta_lines(path: str) -> Iterator[tuple[bool, bytes]]:
    """Yields (whether it begins a line, piece) for the pieces of each line of a (gzipped) FASTA file,
    without line breaks. The file is read in chunks of FASTA_CHUNK_SIZE, a longer sequence line coming
    in several pieces, so a sequence written on a single line is never held in memory whole; header
    lines always come whole."""
    at_line_start, rest = True, b''
    with open_binary(path) as f:
        while chunk := f.read(FASTA_CHUNK_SIZE):
            lines = (rest + chunk).split(b'\n')
            rest = lines.pop()
            for line in lines:
                # a line ending right after the piece of it yielded from the last chunk leaves nothing
                if line or at_line_start:
                    yield at_line_start, line
                at_line_start = True
            if rest and not (at_line_start and rest.startswith(b'>')) and len(rest) > FASTA_CHUNK_SIZE:
                yield at_line_start, rest
                at_line_start, rest = False, b''
    if rest:
        yield at_line_start, rest


def fasta_header_id(header: str) -> str:
    """The sequence ID of a FASTA header line without its '>': the text up to the first whitespace.

    Raises:
        ValueError: for a header holding no ID
    """
    words = header.split(maxsplit=1)
    if not words:
        raise ValueError(f'FASTA header without a sequence ID: ">{header.rstrip()}"')
    return words[0]


def fasta_ids(path: str) -> list[str]:
    """The sequence IDs (the header up to the first whitespace) of a (gzipped) FASTA file, in file order."""
    return [fasta_header_id(piece[1:].decode()) for at_line_start, piece in fasta_lines(path)
            if at_line_start and piece.startswith(b'>')]


def get_seqids_from_gff(gff_file: str) -> set[str]:
    return {line.split('\t')[0] for line in gff_feature_lines(gff_file)}


##### GFF start/end to GeenuFF #####
def strand_as_bool(strand):
    # anything but '+' or '-' raises; use strand_or_none below where that is not fatal
    if strand == '+':
        return True
    elif strand == '-':
        return False
    else:
        raise ValueError('cannot interpret strand "{}"'.format(strand))


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


##### Reverse complement #####

def mk_rc_key():
    fw = "ACGTMRWSYKVHDBN"
    rv = "TGCAKYWSRMBDHVN"
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
            raise KeyError('{} caused by non DNA character {}'.format(e, base))
    return rc_seq


##### Start/Stop codon detection #####

START_CODON = 'ATG'
STOP_CODONS = ['TAG', 'TGA', 'TAA']


def substr_seq(seq, start, end, is_plus_strand):
    """returns a substring of sequence according to geenuff coordinates and strand type"""
    if is_plus_strand:
        return seq[start:end]
    else:
        return seq[end + 1:start + 1]


def spliced_cds_sequence(seq: str, cds_pieces: Sequence[GFFObject], is_plus_strand: bool) -> str:
    """joins raw GFF CDS pieces (each with a 1-based, inclusive .start/.end, sorted by
    ascending genomic start) into one continuous, 5' to 3' oriented sequence"""
    ordered = cds_pieces if is_plus_strand else reversed(cds_pieces)
    parts = []
    for piece in ordered:
        start, end = get_geenuff_start_end(piece.start, piece.end, is_plus_strand)
        part = substr_seq(seq, start, end, is_plus_strand)
        parts.append(part if is_plus_strand else ''.join(reverse_complement(part)))
    return ''.join(parts)


def has_inframe_stop_codon(cds_seq: str) -> bool:
    """checks a fully spliced, 5' to 3' oriented CDS sequence for a stop codon before the
    final codon; a truncated (not a multiple of 3) trailing partial codon is ignored"""
    last_full_codon_start = len(cds_seq) - len(cds_seq) % 3 - 3
    codons = (cds_seq[i:i + 3] for i in range(0, last_full_codon_start, 3))
    return any(codon in STOP_CODONS for codon in codons)


##### SQL alchemy core queue control #####

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
        logger.info(f'All core queues executed ({n_elements} elements)')


class CoreQueue(object):
    def __init__(self, action):
        self.queue = []
        self.action = action
