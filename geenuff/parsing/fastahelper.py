"""Opens plain, gzipped or zipped FASTA and GFF3 files and reads FASTA files in chunks.

text_or_gzip_open and FastaParser are adapted from https://github.com/janinamass/dustdas; credit for
the original code: Janina Mass.
"""
import io
import zipfile
from collections.abc import Iterator
from typing import BinaryIO, TextIO

from zlib_ng import gzip_ng_threaded

# bytes read at a time; from 64 KiB to 1 MiB reading is about equally fast, smaller chunks adding
# Python overhead per byte and larger ones getting slower again; 256 KiB was fastest for plain and
# gzipped files alike (small test on one machine)
FASTA_CHUNK_SIZE = 1 << 18
# whitespace in a FASTA sequence line, never part of the sequence
SEQUENCE_WHITESPACE = b' \t\n\r\v\f'


def detect_format(path: str) -> str:
    """'gzip', 'zip' or 'text', by a file's first bytes rather than its name."""
    with open(path, 'rb') as f:
        start = f.read(4)
    if start.startswith(b'\x1f\x8b'):
        return 'gzip'
    if start == b'PK\x03\x04':
        return 'zip'
    return 'text'


def text_or_gzip_open(path: str, mode: str = 'r') -> TextIO | BinaryIO:
    """Opens a plain, gzipped or zipped file for reading, as text ('r') or bytes ('rb'); a gzipped
    one is decompressed by zlib-ng in a thread of its own.

    Raises:
        ValueError: for a mode other than 'r' or 'rb', or a zip file not holding exactly one file
    """
    if mode not in ('r', 'rb'):
        raise ValueError(f'mode must be "r" or "rb", not "{mode}"')
    binary = mode == 'rb'
    file_format = detect_format(path)
    if file_format == 'gzip':
        return gzip_ng_threaded.open(path, 'rb' if binary else 'rt', threads=1)
    if file_format == 'zip':
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            if len(names) != 1:
                raise ValueError(f'"{path}" holds {len(names)} files; a zip file must hold exactly one')
            f = archive.open(names[0])
        return f if binary else io.TextIOWrapper(f)
    return open(path, mode)


class FastaParser:
    @staticmethod
    def pieces(path: str) -> Iterator[tuple[bool, bytes]]:
        """Yields (whether it is a header line, piece) for a (gzipped) FASTA file, the pieces joined
        giving the file byte for byte: each header line whole, with its '>' and line break, and the
        sequence lines between them in pieces of up to FASTA_CHUNK_SIZE, line breaks included, so a
        sequence is never held in memory whole.

        Raises:
            ValueError: for a file not starting with '>', a header holding no ID (see header_id) or an
                ID given twice
        """
        seen = set()

        def checked(header: bytes) -> bytes:
            seq_id = FastaParser.header_id(header[1:].decode())
            if seq_id in seen:
                raise ValueError(f'"{path}" holds the sequence ID "{seq_id}" more than once; FASTA sequence IDs '
                                 f'must be unique')
            seen.add(seq_id)
            return header

        at_line_start, rest = True, b''
        with text_or_gzip_open(path, 'rb') as f:
            chunk = f.read(FASTA_CHUNK_SIZE)
            if chunk and not chunk.startswith(b'>'):
                raise ValueError(f'"{path}" is no FASTA file: it does not start with ">"')
            while chunk:
                block, position, rest = rest + chunk, 0, b''
                # past the block's start, each piece starts a line
                while position < len(block):
                    is_header = (position > 0 or at_line_start) and block.startswith(b'>', position)
                    if is_header:
                        end = block.find(b'\n', position) + 1
                        if not end:  # a header line reaching into the next chunk waits for the rest of it
                            rest = block[position:]
                            break
                    else:
                        end = block.find(b'\n>', position) + 1 or len(block)  # up to the next header line
                    piece = block[position:end]
                    yield is_header, checked(piece) if is_header else piece
                    position = end
                at_line_start = bool(rest) or block.endswith(b'\n')
                chunk = f.read(FASTA_CHUNK_SIZE)
        if rest:
            yield True, checked(rest)

    @staticmethod
    def header_id(header: str) -> str:
        """The sequence ID of a FASTA header line without its '>': the text up to the first whitespace.

        Raises:
            ValueError: for a header holding no ID
        """
        words = header.split(maxsplit=1)
        if not words:
            raise ValueError(f'FASTA header without a sequence ID: ">{header.rstrip()}"')
        return words[0]

    @staticmethod
    def ids(path: str) -> list[str]:
        """The sequence IDs (see header_id) of a (gzipped) FASTA file, in file order."""
        return [FastaParser.header_id(piece[1:].decode()) for is_header, piece in FastaParser.pieces(path)
                if is_header]

    @staticmethod
    def read_fasta(path: str) -> Iterator[tuple[str, str]]:
        """Yields (header line without '>' and trailing whitespace, sequence without whitespace) for
        each sequence of a (gzipped) FASTA file.

        Raises:
            ValueError: see pieces
        """
        header, parts = None, []
        for is_header, piece in FastaParser.pieces(path):
            if not is_header:
                parts.append(piece.translate(None, SEQUENCE_WHITESPACE))
                continue
            if header is not None:
                yield header, FastaParser._joined(parts)
            header, parts = piece[1:].decode().rstrip(), []
        if header is not None:
            yield header, FastaParser._joined(parts)

    @staticmethod
    def _joined(parts: list[bytes]) -> str:
        """The pieces of a sequence as one string, emptying the list, so no more than the bytes and
        the string of one sequence are held at once."""
        sequence = b''.join(parts)
        parts.clear()
        return sequence.decode()
