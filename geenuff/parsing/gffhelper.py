"""Reads the feature lines and attributes of GFF3 files.

GFFObject is adapted from https://github.com/janinamass/dustdas; credit for the original code: Janina
Mass.
"""
from collections.abc import Iterator, Sequence

from geenuff.parsing.fastahelper import text_or_gzip_open


class GFFParser:
    @staticmethod
    def lines(gff_file: str) -> Iterator[str]:
        """Yields every line of a (gzipped) GFF3 file, comments included, up to an embedded FASTA
        section (a ##FASTA line or a line starting with '>')."""
        with text_or_gzip_open(gff_file) as f:
            for line in f:
                if line.startswith(('##FASTA', '>')):
                    return
                yield line

    @staticmethod
    def feature_lines(gff_file: str) -> Iterator[str]:
        """lines without comments and blank lines."""
        return (line for line in GFFParser.lines(gff_file) if line.strip() and not line.startswith('#'))

    @staticmethod
    def seqids(gff_file: str) -> set[str]:
        """The sequence IDs of a GFF3 file's feature lines."""
        return {line.split('\t')[0] for line in GFFParser.feature_lines(gff_file)}


class GFFObject:
    """One feature line of a GFF3 file: its first 9 columns as stripped strings. Of the 9th, only the ID
    and the parents are kept parsed (get_ID, get_Parent); attributes parses the whole column on each
    access. A genome's lines are held by the hundreds of thousands, and the other attributes, often
    most of a line, are rarely read.

    Sources: http://www.ensembl.org/info/website/upload/gff.html (caution, unsecure link) and
            https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

        seqID - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix.
                Rules from the specification are: IDs may contain any characters, but must escape any characters
                not in the set [a-zA-Z0-9.:^*$@!+_?-|]. In particular, IDs may not contain unescaped whitespace
                and must not begin with an unescaped ">".
        source - name of the program that generated this feature, or the data source (database or project name)
                 Examples: "BRAKER", "Helixer", "." (for unknown/unspecified)
        feature - feature type name, e.g. Gene, Variation, Similarity
                  Examples: "gene", "mRNA", "exon"
        start - Start position of the feature, with sequence numbering starting at 1.
        end - End position of the feature, with sequence numbering starting at 1.
        score - A floating point value. It is strongly recommended that E-values be used for sequence similarity
                features, and that P-values be used for ab initio gene prediction features. In reality many GFF
                files don't actually use scores and the column is filled with "." (unknown).
        strand - defined as + (forward) or - (reverse).
        frame - Only applies to CDS features. One of '0', '1' or '2'. '0' indicates that the first base of the
                feature is the first base of a codon, '1' that the second base is the first base of a codon,
                and so on... (biological phase counting: 0, 2, 1)
        attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
            ID, Name, Alias, Parent, Target, Gap, Derives_from, Note,
            Dbxref, Ontology_term, Is_circular
            (All attributes that begin with an uppercase letter are reserved for later use. Attributes that begin with
            a lowercase letter can be used freely by applications.)
            ID - Indicates the ID of the feature. The ID attribute is required for features that have children
                 (e.g. gene and mRNAs), or for those that span multiple lines, but are optional for other features.
                 IDs for each feature must be unique within the scope of the GFF file. In the case of discontinuous
                 features (i.e. a single feature that exists over multiple genomic locations) the same ID may appear
                 on multiple lines. All lines that share an ID must collectively represent a single feature.
                 Example: ID=mrna0001
            Name - Display name for the feature. This is the name to be displayed to the user. Unlike IDs, there is
                   no requirement that the Name be unique within the file.
            Alias - A secondary name for the feature. It is suggested that this tag be used whenever a secondary
                    identifier for the feature is needed, such as locus names and accession numbers. Unlike ID,
                    there is no requirement that Alias be unique within the file.
            Parent - Indicates the parent of the feature. A parent ID can be used to group exons into transcripts,
                     transcripts into genes, and so forth. A feature may have multiple parents. Parent can only be
                     used to indicate a "partof" relationship.
                     Example: Parent=mrna0001
            Target - Indicates the target of a nucleotide-to-nucleotide or protein-to-nucleotide alignment.
                     The format of the value is "target_id start end [strand]", where strand is optional and may be
                     "+" or "-". If the target_id contains spaces, they must be escaped as hex escape %20.
            Gap - The alignment of the feature to the target if the two are not collinear (e.g. contain gaps).
                  The alignment format is taken from the CIGAR format described in the Exonerate documentation.
                  (http://cvsweb.sanger.ac.uk/cgi-bin/cvsweb.cgi/exonerate?cvsroot=Ensembl).
            Derives_from - Used to disambiguate the relationship between one feature and another when the relationship
                           is a temporal one rather than a purely structural "part of" one. This is needed for
                           polycistronic genes. See "PATHOLOGICAL CASES" for further discussion.
            Note - A free text note.
            Dbxref - A database cross-reference. See the section "Ontology Associations and Db Cross-References"
                     here: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md for details.
                     Example: Dbxref=taxon:3702 (NCBI taxon ID for Arabidopsis thaliana.)
            Ontology_term - A cross-reference to an ontology term. See the section "Ontology Associations and Db
                            Cross-References" for details.
            Is_circular - A flag to indicate whether a feature is circular. See extended discussion online.

            Parent, the Alias, Note, DBxref and Ontology_term attributes can have multiple values.
            There can be a variety of custom lowercase attributes."""

    __slots__ = ('seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attribute', '_id',
                 '_parents')

    def __init__(self, gffline: str):
        columns = [column.strip() for column in gffline.split('\t')[:9]]
        (self.seqid, self.source, self.type, self.start, self.end, self.score, self.strand, self.phase,
         self.attribute) = columns
        attributes = GFFObject.parse_attributes(self.attribute)
        self._id = attributes.get('ID')
        self._parents = GFFObject.parents(attributes) or None

    @property
    def attributes(self) -> dict[str, str]:
        """The attribute column parsed (see parse_attributes), anew on each access."""
        return GFFObject.parse_attributes(self.attribute)

    @staticmethod
    def parse_attributes(column: str) -> dict[str, str]:
        """{tag: value} of a GFF3 attribute column, values left unsplit (e.g. 'a,b' for several
        parents). Each tag=value pair is stripped of surrounding whitespace and split at its first '=';
        pairs without '=', as after a trailing ';', are skipped, and the values of a tag given twice are
        joined by ','."""
        attributes = {}
        for pair in column.split(';'):
            tag, separator, value = pair.strip().partition('=')
            if separator:
                attributes[tag] = f'{attributes[tag]},{value}' if tag in attributes else value
        return attributes

    @staticmethod
    def parents(attributes: dict[str, str]) -> list[str]:
        """The distinct Parent IDs of a line's attributes (see parse_attributes) in file order, blank
        ones (as in 'Parent=,tx1') skipped."""
        return [parent for parent in dict.fromkeys(attributes.get('Parent', '').split(',')) if parent]

    @staticmethod
    def with_parents(column: str, parents: Sequence[str]) -> str:
        """A GFF3 attribute column with its Parent attribute replaced by one naming the given parents,
        put last; the other attributes are kept as written."""
        pairs = [pair for pair in column.strip().strip(';').split(';')
                 if pair.strip().partition('=')[0] != 'Parent']
        return ';'.join(pairs + [f'Parent={",".join(parents)}'])

    def get_ID(self) -> str | None:
        return self._id

    def get_Parent(self) -> list[str] | None:
        """The line's parents (see parents), None where there are none."""
        return self._parents

    def __repr__(self) -> str:
        return (f'{self.seqid},{self.source},{self.type},{self.start},{self.end},{self.score},{self.strand},'
                f'{self.phase},{self.attribute}')
