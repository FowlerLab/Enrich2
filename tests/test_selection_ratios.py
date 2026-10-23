"""
Tests for the ratio calculations in :py:class:`enrich2.selection.Selection`.

These read the first and last timepoint columns out of an HDF5 store with
``HDFStore.select``. pandas 3.1.0 rejects the where-clause
``"columns in ['c_0', c_last]"`` (a list mixing a literal and a local variable)
with ``NotImplementedError: unable to collapse Joint Filters``, so these tests
exercise every ``calc_ratios`` log-ratio method and ``calc_simple_ratios``.
"""

import numpy as np
import pandas as pd
import pytest

from enrich2.constants import WILD_TYPE_VARIANT
from enrich2.selection import Selection

LABEL = "variants"

# Three timepoints so that the middle column, c_1, is a distractor: the
# calculations must use only the first (c_0) and last (c_2) columns.
COUNTS = pd.DataFrame(
    {
        "c_0": [100, 20, 30, 40],
        "c_1": [7, 7, 7, 7],
        "c_2": [90, 50, 10, 60],
    },
    index=[WILD_TYPE_VARIANT, "A1T", "G2C", "T3A"],
)
# Reads that were filtered out of COUNTS; only used by the "full" method.
COUNTS_UNFILTERED = pd.DataFrame(
    {
        "c_0": [100, 20, 30, 40, 5],
        "c_1": [7, 7, 7, 7, 7],
        "c_2": [90, 50, 10, 60, 15],
    },
    index=[WILD_TYPE_VARIANT, "A1T", "G2C", "T3A", "C4G"],
)


@pytest.fixture
def make_selection(tmp_path):
    """Build a Selection whose HDF5 store holds hand-written count tables."""
    stores = []

    def _make(logr_method):
        selection = Selection()
        selection.name = "test_selection"
        # Selection.timepoints is derived from the keys of libraries.
        selection.libraries = {0: [], 1: [], 2: []}
        selection._labels = [LABEL]
        selection.logr_method = logr_method
        selection.store = pd.HDFStore(str(tmp_path / "selection.h5"))
        stores.append(selection.store)
        for suffix, counts in (
            ("counts", COUNTS),
            ("counts_unfiltered", COUNTS_UNFILTERED),
        ):
            selection.store.put(
                f"/main/{LABEL}/{suffix}",
                counts,
                format="table",
                data_columns=counts.columns,
            )
        return selection

    yield _make
    for store in stores:
        store.close()


def expected_log_ratios(shared_counts):
    """Reference implementation of the ratio calculation, independent of Enrich2."""
    observed = COUNTS[["c_0", "c_2"]].to_numpy(dtype=float) + 0.5
    log_ratios = np.log(observed) - np.log(shared_counts)
    logratio = log_ratios[:, 1] - log_ratios[:, 0]
    variance = np.sum(1.0 / observed, axis=1) + np.sum(1.0 / shared_counts)
    return pd.DataFrame(
        {
            "score": logratio,
            "SE": np.sqrt(variance),
            "logratio": logratio,
            "variance": variance,
        },
        index=COUNTS.index,
    )


@pytest.mark.parametrize(
    "logr_method, shared_counts",
    [
        # Counts of the wild type variant in the first and last timepoints.
        (
            "wt",
            COUNTS.loc[[WILD_TYPE_VARIANT], ["c_0", "c_2"]].to_numpy(dtype=float) + 0.5,
        ),
        # Total counts over the filtered table.
        ("complete", COUNTS[["c_0", "c_2"]].sum().to_numpy(dtype=float) + 0.5),
        # Total counts over the unfiltered table.
        ("full", COUNTS_UNFILTERED[["c_0", "c_2"]].sum().to_numpy(dtype=float) + 0.5),
    ],
)
def test_calc_ratios(make_selection, logr_method, shared_counts):
    selection = make_selection(logr_method)

    selection.calc_ratios(LABEL)

    scores = selection.store.select(f"/main/{LABEL}/scores")
    pd.testing.assert_frame_equal(
        scores.sort_index(),
        expected_log_ratios(shared_counts).sort_index(),
        check_dtype=False,
    )


def test_calc_simple_ratios(make_selection):
    selection = make_selection("wt")

    selection.calc_simple_ratios(LABEL)

    scores = selection.store.select(f"/main/{LABEL}/scores")
    c_0 = COUNTS["c_0"].to_numpy(dtype=float)
    c_2 = COUNTS["c_2"].to_numpy(dtype=float)
    ratio = (c_2 / c_2.sum()) / (c_0 / c_0.sum())
    expected = pd.DataFrame(
        {"ratio": ratio, "score": np.log2(ratio), "SE": np.nan},
        index=COUNTS.index,
    )
    pd.testing.assert_frame_equal(
        scores.sort_index()[["ratio", "score", "SE"]],
        expected.sort_index(),
        check_dtype=False,
    )
