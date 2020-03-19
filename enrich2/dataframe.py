import pandas as pd
import numpy as np
import logging
from .constants import WILD_TYPE_VARIANT
import collections
from .variant import mutation_count, re_protein, re_coding, re_noncoding
from .barcodemap import re_barcode, re_identifier
from .constants import AA_CODES
from .storemanager import ELEMENT_LABELS
from .sfmap import AA_LIST, NT_LIST

logger = logging.getLogger(__name__)

SingleMut = collections.namedtuple("SingleMut", ["pre", "post", "pos", "key"])


def validate_index(index, element):
    """
    Return a boolean list for which index values are valid for the given
    element type.
    """
    if element not in ELEMENT_LABELS:
        raise ValueError("Invalid element label '{}'".format(element))

    if element == "barcodes":
        retval = [re_barcode.match(x) is not None for x in index]
    elif element == "identifiers":
        retval = [re_identifier.match(x) is not None for x in index]
    elif element == "variants":
        pass
    elif element == "synonymous":
        pass
    else:
        raise NotImplementedError("Unimplemented element type '{}'" "".format(element))


def single_mutation_index(index):
    """
    Return a filtered pandas Index containing only single mutations. Filtering
    also removes unrecognized amino acids (denoted by ``"???"``) caused by
    some indels.

    *index* the index to be filtered for single mutations.
    """
    return pd.Index(x for x in index if mutation_count(x) == 1)


def filter_coding_index(index):
    """
    Return a filtered pandas Index with any unrecognized amino acids (denoted
    by ``"???"``) removed. These are caused by some frame shift mutations.

    *index* the index to be filtered.
    """
    return pd.Index(x for x in index if "???" not in x)


def single_mutations_to_tuples(index):
    """
    Return a list of SingleMut namedtuples for each single mutation in the
    *index*. The type of index (noncoding DNA, coding DNA, or protein) is
    automatically detected.

    Position value in the tuple is stored as an integer.

    If the *index* is a protein index, the amino acids are referred to by
    single-letter codes not three-letter codes.

    *index* is the index to convert to SingleMut tuples.

    Raises a ValueError if non-single mutations are included in *index*.

    Raises a ValueError if one of the *index* entries cannot be parsed.

    Raises an IndexError if the *index* is empty.
    """
    if any(mutation_count(x) != 1 for x in index):
        raise ValueError(
            "Non-single mutations cannot be converted into " "SingleMut tuples."
        )

    # identify the type of index
    try:
        if re_noncoding.match(index[0]):
            is_protein = False
            expression = re_noncoding
        elif re_coding.match(index[0]):
            is_protein = False
            expression = re_coding
        elif re_protein.match(index[0]):
            is_protein = True
            expression = re_protein
        else:
            raise ValueError("Unrecognized HGVS string.")
    except IndexError:
        raise IndexError("Cannot convert empty index to tuples.")

    # perform the regular expression matches and create the SingleMut tuples
    tuples = list()
    for x in index:
        m = expression.match(x)
        if m is None:
            raise ValueError("Unrecognized HGVS string.")
        else:
            if is_protein:  # convert to single-letter amino acid code
                tuples.append(
                    SingleMut(
                        AA_CODES[m.group("pre")],
                        AA_CODES[m.group("post")],
                        int(m.group("pos")),
                        m.group("match"),
                    )
                )
            else:
                tuples.append(
                    SingleMut(
                        m.group("pre"),
                        m.group("post"),
                        int(m.group("pos")),
                        m.group("match"),
                    )
                )

    return tuples


def fill_position_gaps(positions, gap_size):
    """
    Create a list of integer positions with gaps filled in. Used by
    :py:func:`singleton_dataframe`.

    Args:
        positions (list): integer positions
        gap_size (int): maximum length of gap that will be filled

    Returns:
        list: sorted list of unique integer positions with gaps filled
    """
    if len(positions) == 0:
        raise ValueError("Empty positions list.")

    # uniqify and sort
    positions = sorted(list(set(positions)))

    # fill in short gaps
    fill = set()
    for i in xrange(len(positions) - 1):
        delta = positions[i + 1] - positions[i]
        if delta > 1 and delta <= gap_size:
            fill.update(positions[i] + n + 1 for n in xrange(delta))
    fill.update(positions)

    return sorted(list(fill))


def singleton_dataframe(
    values, wt, gap_size=5, coding=True, plot_wt_score=True, aa_list=AA_LIST
):
    """
    Prepare data for plotting as a sequence-function map. Returns a data frame
    suitable for plotting as heat map data and a wild type sequence extracted
    from the variant information.

    The type of variants stored is automatically detected, and the index will
    be filtered for single mutations.

    The data frame has amino acids or nucleotides as columns and positions with
    rows. If there are no mutations at a given position, it will not appear in
    the data frame unless this gap is filled with rows containing no data. The
    wild type sequence entry for these rows will be blank.

    Args:
        values (|pd_Series|): data values (typically scores or counts)

        wt (WildTypeSequence): wild type for the data

        gap_size (int): maximum length of missing data gap that will be filled

        coding (bool): True for amino acid data, False for nucleotide

        plot_wt_score (bool): True if the wild type positions should have the
            wild type score, False if they should be missing

    Returns:
        tuple: two-element tuple containing a |pd_DataFrame| filled with the
        data values and a list of single-character wild type values
    """
    if len(values.index) == 0:
        raise ValueError(
            "Cannot process an empty data frame [{}]".format(wt.parent_name)
        )

    # save the wild type score for later
    if plot_wt_score:
        try:
            wt_score = values[WILD_TYPE_VARIANT]
        except KeyError:
            logger.warning("Wild type score not measured, will be missing in " "plots")
            wt_score = np.nan

    # select only rows with singleton mutations
    values = values[filter_coding_index(single_mutation_index(values.index))]

    # parse out the information from the index
    index_tuples = single_mutations_to_tuples(values.index)

    # create and populate the DataFrame
    # get sorted, unique list of positions that have a mutation
    positions = fill_position_gaps([x.pos for x in index_tuples], gap_size=gap_size)
    # initialize the DataFrame
    if coding:
        columns = aa_list
    else:
        columns = NT_LIST
    frame = pd.DataFrame(np.nan, columns=columns, index=positions)
    # populate the DataFrame
    for x in index_tuples:
        frame.loc[x.pos, x.post] = values.loc[x.key]

    # create a dictionary of position->nucleotide/amino acid
    wt_dict = dict(wt.position_tuples(protein=coding))

    # convert subset of the wild type dictionary into sequence
    try:
        wt_sequence = "".join(wt_dict[x] for x in positions)
    except KeyError:
        raise ValueError("Inconsistent wild type positions [{}]".format(wt.parent_name))

    # double-check that the wild type is consistent with the data frame
    for x in index_tuples:
        if x.pos in wt_dict:
            if x.pre != wt_dict[x.pos]:
                raise ValueError(
                    "Inconsistent wild type sequence [{}]".format(wt.parent_name)
                )

    # add wild type scores if desired
    if plot_wt_score:
        for p in positions:
            frame.loc[p, wt_dict[p]] = wt_score

    return (frame, wt_sequence)
