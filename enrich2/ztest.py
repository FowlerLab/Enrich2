"""z-score helper functions.

This module contains functions for calculating z-scores and corresponding
p-values for pairs of scores.
"""

import numpy as np
import scipy.stats as stats


def ztest_pair(df1, df2):
    """z-test for elements in two data frames.

    Takes two data frames with ``'score'`` and ``'SE'`` columns and returns
    a new data frame containing the scores and standard errors and result for
    pairwise comparison of elements in both input data frames.

    Parameters
    ----------
    df1 : pandas.DataFrame
        The first data frame. Must contain ``'score'`` and ``'SE'`` columns.

    df2 : pandas.DataFrame
        The second data frame. Must contain ``'score'`` and ``'SE'`` columns.

    Returns
    -------
    pandas.DataFrame
        Result data frame containing scores (``'score1'`` and ``'score2'``),
        standard errors (``'SE1'``, ``'SE2'``), z-scores (``'z'``), and
        p-values (``'pvalue_raw'``) for each element found in both data frames.

    Raises
    ------
    To be added.
    """
    shared = df1.loc[:, ("score", "SE")].merge(
        df2.loc[:, ("score", "SE")],
        how="inner",
        left_index=True,
        right_index=True,
        suffixes=("1", "2"),
    )
    shared["z"] = np.abs(shared["score1"] - shared["score2"]) / np.sqrt(
        shared["SE1"] ** 2 + shared["SE2"] ** 2
    )
    shared["pvalue_raw"] = 2 * stats.norm.sf(shared["z"])
    return shared


def ztest_single(df, score, se):
    """z-test for comparing elements in a data frame to a single score and SE.

    Takes a data frames with ``'score'`` and ``'SE'`` columns and a score and
    standard error to compare them to and returns a new data frame containing
    the scores, standard errors, and result for the pairwise comparisons.

    Parameters
    ----------
    df : pandas.DataFrame
        The data frame. Must contain ``'score'`` and ``'SE'`` columns.

    score : float
        Score used for comparison to elements in the data frame.

    se : float
        Standard error used for comparison to elements in the data frame.

    Returns
    -------
    pandas.DataFrame
        Result data frame containing score (``'score'``), standard error
        (``'SE'``), z-score (``'z'``), and p-value (``'pvalue_raw'``) for each
        element in the input data frame.

    Raises
    ------
    To be added.
    """
    result = df.loc[:, ("score", "SE")]
    result["z"] = np.abs(result["score"] - score) / np.sqrt(result["SE"] ** 2 + se ** 2)
    result["pvalue_raw"] = 2 * stats.norm.sf(result["z"])
    return result
