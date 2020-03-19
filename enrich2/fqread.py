from __future__ import print_function
from sys import stderr
import os.path
import re
import string
import itertools
import bz2
import gzip
from array import array

# The following regex is referenced by line number in the class documentation.
# Matches FASTQ headers based on the following pattern (modify as needed):
# @<MachineName>:<Lane>:<Tile>:<X>:<Y>:<Chastity>#<IndexRead>/<ReadNumber>
header_pattern = re.compile(
    "@(?P<MachineName>.+)"
    ":(?P<Lane>\d+)"
    ":(?P<Tile>\d+)"
    ":(?P<X>\d+)"
    ":(?P<Y>\d+)"
    ":(?P<Chastity>[01])"
    "#(?P<IndexRead>\d)"
    "/(?P<ReadNumber>\d)"
)


BUFFER_SIZE = 100000  # empirically optimized for reading FASTQ files


dna_trans = string.maketrans("actgACTG", "tgacTGAC")


class FQRead(object):
    """
    Stores a single record from a FASTQ_ file. Quality values are stored 
    internally as a list of integer `Phred quality scores \
    <http://www.phrap.com/phred/#qualityscores>`_. The *qbase* parameter is 
    the ASCII value that correponds to Phred score of 0. The *sequence* and 
    *quality* strings must be the same length. 
    """

    # use slots for memory efficiency
    __slots__ = ("header", "sequence", "header2", "quality", "qbase")

    def __init__(self, header, sequence, header2, quality, qbase=33):
        if len(sequence) != len(quality):
            raise ValueError("different lengths for sequence and quality")
        elif header[0] != "@" or header2[0] != "+":
            raise ValueError("improperly formatted FASTQ record")
        else:
            self.header = header
            self.sequence = sequence
            self.header2 = header2
            # quality is a list of integers
            self.quality = [x - qbase for x in array("b", quality).tolist()]
            self.qbase = qbase

    def __str__(self):
        """
        Reformat as a four-line FASTQ_ record. This method converts the 
        integer quality values back into a string.
        """
        return "\n".join(
            [
                self.header,
                self.sequence,
                self.header2,
                array("b", [x + self.qbase for x in self.quality]).tostring(),
            ]
        )

    def __len__(self):
        """
        Object length is the length of the sequence.
        """
        return len(self.sequence)

    def trim(self, start=1, end=None):
        """
        Trims this :py:class:`~fqread.FQRead` to contain bases between 
        *start* and *end* (inclusive). Bases are numbered starting at 1.
        """
        self.sequence = self.sequence[start - 1 : end]
        self.quality = self.quality[start - 1 : end]

    def trim_length(self, length, start=1):
        """
        Trims this :py:class:`~fqread.FQRead` to contain *length* bases, 
        beginning with *start*. Bases are numbered starting at 1.
        """
        self.trim(start=start, end=start + length - 1)

    def revcomp(self):
        """
        Reverse-complement the sequence in place. Also reverses the array of 
        quality values.
        """
        self.sequence = self.sequence.translate(dna_trans)[::-1]
        self.quality = self.quality[::-1]

    def header_information(self, pattern=header_pattern):
        """header_information(pattern=header_pattern)

        Parses the first FASTQ_ header (@ header) and returns a dictionary. 
        Dictionary keys are the named groups in the regular expression 
        *pattern*. Unnamed matches are ignored. Integer values are converted 
        from strings to integers.

        The default pattern matches a header in the format::

            @<MachineName>:<Lane>:<Tile>:<X>:<Y>:<Chastity>#<IndexRead>/<ReadNumber>

        """
        match = pattern.match(self.header)
        if match is None:
            return None
        else:
            header_dict = match.groupdict()
            for key in header_dict:
                if header_dict[key].isdigit():
                    header_dict[key] = int(header_dict[key])
            return header_dict

    def min_quality(self):
        """
        Return the minimum Phred-like quality score.
        """
        return min(self.quality)

    def mean_quality(self):
        """
        Return the average Phred-like quality score.
        """
        return float(sum(self.quality)) / len(self)

    def is_chaste(self, raises=True):
        """
        Returns ``True`` if the chastity bit is set in the header. The 
        regular experession used by :py:meth:`header_information` must  
        include a ``'Chastity'`` match that equals ``1`` if the read is 
        chaste.

        If ``raises`` is ``True``, raises an informative error if the 
        chastity information in the header is not found. Otherwise, a 
        read without chastity information is treated as unchaste.
        """
        try:
            if self.header_information()["Chastity"] == 1:
                return True
            else:
                return False
        except KeyError:  # no 'Chastity' in pattern
            if raises:
                raise KeyError("No chastity bit in FASTQ header pattern")
            else:
                return False
        except TypeError:  # no header match (unexpected format)
            if raises:
                raise ValueError("Unexpected FASTQ header format")
            else:
                return False


def split_fastq_path(fname):
    """
    Check that *fname* exists and has a valid FASTQ_ file extension. Valid 
    file extensions are ``.fastq`` or ``.fq``, optionally followed by ``.gz`` 
    or ``.bz2`` if the file is compressed. 

    Returns a tuple containing the directory, the file base name with no 
    extension, the FASTQ_ file extension used, and the compression format 
    (``"gz"``, ``"bz2"``, or ``None``).

    Raises an ``IOError`` if the file doesn't exist. Returns ``None`` if the 
    file extension is not recognized.
    """
    if os.path.isfile(fname):
        compression = None
        head, tail = os.path.split(fname)
        base, ext = os.path.splitext(tail)
        if ext.lower() == ".bz2":
            compression = "bz2"
            base, ext = os.path.splitext(base)
        elif ext.lower() == ".gz":
            compression = "gz"
            base, ext = os.path.splitext(base)
        if ext.lower() in (".fq", ".fastq"):
            return (head, base, ext, compression)
        else:
            print(
                "Warning: unexpected file extension for '{fname}'".format(fname=fname),
                file=stderr,
            )
            return None
    else:
        raise IOError("file '{fname}' doesn't exist".format(fname=fname))


def create_compressed_outfile(fname, compression):
    """
    Utility function for opening compressed output files. Accepted values for 
    *compression* are ``"gz"``, ``"bz2"``, or ``None``. Returns a file handle 
    of the appropriate type opened for writing. Existing files with the same 
    name are overwritten.
    """
    if compression == "bz2":
        handle = bz2.BZ2File(fname + ".bz2", "w")
    elif compression == "gz":
        handle = gzip.GzipFile(fname + ".gz", "w")
    elif compression is None:
        handle = open(fname, "w")
    else:
        raise IOError("unrecognized compression mode '{mode}'".format(mode=compression))
    return handle


def read_fastq(fname, filter_function=None, buffer_size=BUFFER_SIZE, qbase=33):
    """
    Generator function for reading from FASTQ_ file *fname*. Yields an 
    :py:class:`~fqread.FQRead` object for each FASTQ_ record in the file. The 
    *filter_function* must operate on an :py:class:`~fqread.FQRead` object 
    and return ``True`` or ``False``. If the result is ``False``, the record 
    will be skipped silently.

    .. note:: To read multiple files in parallel (such as index or \
        forward/reverse reads), use :py:func:`read_fastq_multi` instead.
    """
    _, _, _, compression = split_fastq_path(fname)
    if compression is None:  # raw FASTQ
        handle = open(fname, "rU")
    elif compression == "bz2":
        handle = bz2.BZ2File(fname, "rU")
    elif compression == "gz":
        handle = gzip.GzipFile(fname, "rU")
    else:
        raise IOError("unrecognized compression mode '{mode}'".format(mode=compression))

    eof = False
    leftover = ""

    while not eof:
        buf = handle.read(buffer_size)
        if len(buf) < buffer_size:
            eof = True

        buf = leftover + buf  # prepend partial record from previous buffer
        lines = buf.split("\n")
        fastq_count = len(lines) / 4

        if not eof:  # handle lines from the trailing partial FASTQ record
            dangling = len(lines) % 4
            if dangling == 0:  # quality line (probably) incomplete
                dangling = 4
                fastq_count = fastq_count - 1
            # join the leftover lines back into a string
            leftover = "\n".join(lines[len(lines) - dangling :])

        # index into the list of lines to pull out the FASTQ records
        for i in xrange(fastq_count):
            # (header, sequence, header2, quality)
            fq = FQRead(*lines[i * 4 : (i + 1) * 4], qbase=qbase)
            if filter_function is None:  # no filtering
                yield fq
            elif filter_function(fq):  # passes filtering
                yield fq
            else:  # fails filtering
                continue

    handle.close()


def read_fastq_multi(
    fnames, filter_function=None, buffer_size=BUFFER_SIZE, match_lengths=True, qbase=33
):
    """
    Generator function for reading from multiple FASTQ_ files in parallel. 
    The argument *fnames* is an iterable of FASTQ_ file names. Yields a 
    tuple of :py:class:`~fqread.FQRead` objects, one for each file in 
    *fnames*. The *filter_function* must operate on an :py:class:`FQRead` 
    object and return ``True`` or ``False``. If the result is ``False`` for 
    any :py:class:`FQRead` in the tuple, the entire tuple will be skipped.

    If *match_lengths* is ``True``, the generator will yield ``None`` if the 
    files do not contain the same number of FASTQ_ records. Otherwise, it 
    will silently ignore partial records.
    """
    fq_generators = list()
    for f in fnames:
        fq_generators.append(
            read_fastq(f, filter_function=None, buffer_size=BUFFER_SIZE, qbase=qbase)
        )

    for records in itertools.izip_longest(*fq_generators, fillvalue=None):
        if None in records:  # mismatched file lengths
            if match_lengths:
                yield None
            else:
                break  # shortest FASTQ file is empty, so we're done
        if filter_function is None:  # no filtering
            yield records
        elif all(filter_function(x) for x in records):  # pass filtering
            yield records
        else:  # fail filtering
            continue


def fastq_filter_chastity(fq):
    """
    Filtering function for :py:func:`read_fastq` and 
    :py:func:`read_fastq_multi`. Returns ``True`` if the 
    :py:class:`~fqread.FQRead` object *fq* is chaste.
    """
    return fq.is_chaste()
