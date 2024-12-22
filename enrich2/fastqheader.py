import re

# Matches FASTQ headers based on the following pattern (modify as needed):
# @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number>

# Example: @M02564:876:000000000-L3775:1:1101:16862:1800 1:N:0:TCACTCGA+TAACGGTT
# Sample number contains indexes if they are present.

# See: https://help.basespace.illumina.com/files-used-by-basespace/fastq-files

# Note: this regex is currently unused since it is not needed to support the legacy chastity filtering feature
new_header_pattern = re.compile(
    r"""
    @(?P<Instrument>[^:]+):
    (?P<RunNumber>\d+):
    (?P<FlowcellID>[^:]+):
    (?P<Lane>\d+):
    (?P<Tile>\d+):
    (?P<XPos>\d+):
    (?P<YPos>\d+)
    \s
    (?P<Read>\d+):
    (?P<IsFiltered>[YN]):
    (?P<ControlNumber>[^:]+):
    (?P<SampleNumber>[^:]+)
    """,
    re.VERBOSE,
)

# Matches FASTQ headers based on the following pattern (modify as needed):
# @<MachineName>:<Lane>:<Tile>:<X>:<Y>:<Chastity>#<IndexRead>/<ReadNumber>
old_header_pattern = re.compile(
    r"""
    @(?P<MachineName>.+):
    (?P<Lane>\d+):
    (?P<Tile>\d+):
    (?P<X>\d+):
    (?P<Y>\d+):
    (?P<Chastity>[01])#
    (?P<IndexRead>\d)/
    (?P<ReadNumber>\d)
    """,
    re.VERBOSE,
)

def parse_fastq_header(fq, pattern=old_header_pattern):
    """Parse the read's FASTQ_ header and return key-value pairs.

    Parses the first FASTQ_ header (@ header) and returns a dictionary.
    Dictionary keys are the named groups in the regular expression
    *pattern*. Unnamed matches are ignored. Integer values are converted
    from strings to integers.

    The default pattern matches a header in the format::

        @<MachineName>:<Lane>:<Tile>:<X>:<Y>:<Chastity>#<IndexRead>/<ReadNumber>

    """
    match = pattern.match(fq.header)
    if match is None:
        return None
    else:
        header_dict = match.groupdict()
        for key in header_dict:
            if header_dict[key].isdigit():
                header_dict[key] = int(header_dict[key])
        return header_dict


def fastq_read_is_chaste(self, raises=True):
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
