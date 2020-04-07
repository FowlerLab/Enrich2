from __future__ import print_function
from .barcode import BarcodeSeqLib
from .barcodevariant import BcvSeqLib
from .barcodeid import BcidSeqLib
from .basic import BasicSeqLib
from .idonly import IdOnlySeqLib
from .overlap import OverlapSeqLib
from .config_check import seqlib_type
from .storemanager import StoreManager
import os
import pandas as pd
import numpy as np
import logging
import statsmodels.api as sm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scipy.stats as stats
from .sfmap import sfmap_plot
from .plots import (
    fit_axes,
    fit_axes_text,
    volcano_plot,
    configure_axes,
    plot_colors,
    weights_plot,
)
from .constants import WILD_TYPE_VARIANT, SYNONYMOUS_VARIANT
from .variant import protein_variant
from .dataframe import singleton_dataframe


#: map class names to class definitions to avoid use of globals()
SEQLIB_CLASSES = {
    "BarcodeSeqLib": BarcodeSeqLib,
    "BcvSeqLib": BcvSeqLib,
    "BcidSeqLib": BcidSeqLib,
    "BasicSeqLib": BasicSeqLib,
    "IdOnlySeqLib": IdOnlySeqLib,
    "OverlapSeqLib": OverlapSeqLib,
}


def regression_apply(row, timepoints, weighted):
    """
    :py:meth:`pandas.DataFrame.apply` apply function for calculating 
    enrichment using linear regression. If *weighted* is ``True`` perform
    weighted least squares; else perform ordinary least squares.

    Weights for weighted least squares are included in *row*.

    Returns a :py:class:`pandas.Series` containing regression coefficients,
    residuals, and statistics.
    """
    # retrieve log ratios from the row
    y = row[["L_{}".format(t) for t in timepoints]]

    # re-scale the x's to fall within [0, 1]
    xvalues = [x / float(max(timepoints)) for x in timepoints]

    # perform the fit
    X = sm.add_constant(xvalues)  # fit intercept
    if weighted:
        W = row[["W_{}".format(t) for t in timepoints]]
        fit = sm.WLS(y, X, weights=W).fit()
    else:
        fit = sm.OLS(y, X).fit()

    # re-format as a data frame row
    values = np.concatenate(
        [fit.params, [fit.bse["x1"], fit.tvalues["x1"], fit.pvalues["x1"]], fit.resid]
    )
    index = ["intercept", "slope", "SE_slope", "t", "pvalue_raw"] + [
        "e_{}".format(t) for t in timepoints
    ]
    return pd.Series(data=values, index=index)


class Selection(StoreManager):
    """
    Class for a single selection replicate, consisting of multiple 
    timepoints. This class coordinates :py:class:`~seqlib.seqlib.SeqLib` 
    objects.
    """

    store_suffix = "sel"
    treeview_class_name = "Selection"

    def __init__(self):
        StoreManager.__init__(self)
        self.libraries = dict()
        self.barcode_maps = dict()
        self._wt = None
        self.logger = logging.getLogger("{}.{}".format(__name__, self.__class__))

    def _children(self):
        """
        Return the :py:class:`~seqlib.seqlib.SeqLib` objects as a list, 
        sorted by timepoint and then by name.
        """
        libs = list()
        for tp in self.timepoints:
            libs.extend(sorted(self.libraries[tp], key=lambda x: x.name))
        return libs

    def remove_child_id(self, tree_id):
        """
        Remove the reference to a :py:class:`~seqlib.seqlib.SeqLib` with 
        Treeview id *tree_id*. Deletes empty time points.
        """
        empty = None
        for tp in self.libraries:
            tp_ids = [lib.treeview_id for lib in self.libraries[tp]]
            if tree_id in tp_ids:
                del self.libraries[tp][tp_ids.index(tree_id)]
                if len(self.libraries[tp]) == 0:
                    empty = tp
                break  # found the id, stop looking
        if empty is not None:
            del self.libraries[empty]

    @property
    def timepoints(self):
        return sorted(self.libraries.keys())

    @property
    def wt(self):
        if self.has_wt_sequence():
            if self._wt is None:
                self._wt = self.children[0].wt.duplicate(self.name)
            return self._wt
        else:
            if self._wt is not None:
                raise ValueError(
                    "Selection should not contain wild type sequence [{}]".format(
                        self.name
                    )
                )
            else:
                return None

    def configure(self, cfg, configure_children=True):
        """
        Set up the :py:class:`~selection.Selection` using the *cfg* object, 
        usually from a ``.json`` configuration file.

        If *configure_children* is false, do not configure the children in 
        *cfg*.
        """
        StoreManager.configure(self, cfg)
        self.logger = logging.getLogger(
            "{}.{} - {}".format(__name__, self.__class__.__name__, self.name)
        )
        if configure_children:
            if "libraries" not in cfg:
                raise KeyError(
                    "Missing required config value {} [{}]".format(
                        "libraries", self.name
                    )
                )

            for lib_cfg in cfg["libraries"]:
                libtype = seqlib_type(lib_cfg)
                if libtype is None:
                    raise ValueError("Unrecognized SeqLib config")
                elif libtype in ("BcvSeqLib", "BcidSeqLib"):
                    lib = SEQLIB_CLASSES[libtype]()
                    # don't re-parse the barcode maps if possible
                    mapfile = lib_cfg["barcodes"]["map file"]
                    if mapfile in self.barcode_maps.keys():
                        lib.configure(lib_cfg, barcode_map=self.barcode_maps[mapfile])
                    else:
                        lib.configure(lib_cfg)
                        self.barcode_maps[mapfile] = lib.barcode_map
                    self.add_child(lib)
                else:
                    # requires that the SeqLib derived classes be imported into the
                    # module namespace using "from x import y" style
                    lib = SEQLIB_CLASSES[libtype]()
                    lib.configure(lib_cfg)
                    self.add_child(lib)

    def validate(self):
        """
        Raises an informative ``ValueError`` if the time points in the analysis are not suitable.

        Calls validate method on all child SeqLibs.
        """
        # check the time points
        if 0 not in self.timepoints:
            raise ValueError("Missing timepoint 0 [{}]".format(self.name))
        if self.timepoints[0] != 0:
            raise ValueError("Invalid negative timepoint [{}]".format(self.name))
        if len(self.timepoints) < 2:
            raise ValueError("Multiple timepoints required [{}]".format(self.name))
        elif len(self.timepoints) < 3 and self.scoring_method in ("WLS", "OLS"):
            raise ValueError(
                "Insufficient number of timepoints for regression scoring [{}]".format(
                    self.name
                )
            )

        # check the wild type sequences
        if self.has_wt_sequence():
            for child in self.children[1:]:
                if self.children[0].wt != child.wt:
                    self.logger.warning("Inconsistent wild type sequences")
                    break

        # check that we're not doing wild type normalization on something with no wild type
        # if not self.has_wt_sequence() and self.logr_method == "wt":
        #    raise ValueError("No wild type sequence for wild type normalization [{}]".format(self.name))

        # validate children
        for child in self.children:
            child.validate()

    def serialize(self):
        """
        Format this object (and its children) as a config object suitable for dumping to a config file.
        """
        cfg = StoreManager.serialize(self)
        cfg["libraries"] = [child.serialize() for child in self.children]
        return cfg

    def add_child(self, child):
        if child.name in self.child_names():
            raise ValueError(
                "Non-unique sequencing library name '{}' [{}]".format(
                    child.name, self.name
                )
            )

        child.parent = self

        # add it to the libraries dictionary
        try:
            self.libraries[child.timepoint].append(child)
        except KeyError:
            self.libraries[child.timepoint] = [child]

    def is_barcodevariant(self):
        """
        Return ``True`` if all :py:class:`~seqlib.seqlib.SeqLib` in the 
        :py:class:`~selection.Selection` are 
        :py:class:`~barcodevariant.BcvSeqLib` objects with 
        the same barcode map, else ``False``.
        """
        return (
            all(isinstance(lib, BcvSeqLib) for lib in self.children)
            and len(self.barcode_maps.keys()) == 1
        )

    def is_barcodeid(self):
        """
        Return ``True`` if all :py:class:`~seqlib.SeqLib` in the 
        :py:class:`~selection.Selection` are 
        :py:class:`~barcodeid.BcidSeqLib` objects with 
        the same barcode map, else ``False``.
        """
        return (
            all(isinstance(lib, BcidSeqLib) for lib in self.children)
            and len(self.barcode_maps.keys()) == 1
        )

    def is_coding(self):
        """
        Return ``True`` if the all :py:class:`~seqlib.seqlib.SeqLib` in the 
        :py:class:`~selection.Selection` count protein-coding variants, else 
        ``False``.
        """
        return all(x.is_coding() for x in self.children)

    def has_wt_sequence(self):
        """
        Return ``True`` if the all :py:class:`~seqlib.seqlib.SeqLib` in the 
        :py:class:`~selection.Selection` have a wild type sequence, else 
        ``False``.
        """
        return all(x.has_wt_sequence() for x in self.children)

    def merge_counts_unfiltered(self, label):
        """
        Counts :py:class:`~seqlib.seqlib.SeqLib` objects and tabulates counts 
        for each timepoint. :py:class:`~seqlib.seqlib.SeqLib` objects from 
        the same timepoint are combined by summing the counts.

        Stores the unfiltered counts under ``/main/label/counts_unfiltered``.
        """
        if self.check_store("/main/{}/counts_unfiltered".format(label)):
            return

        # calculate counts for each SeqLib
        self.logger.info("Counting for each time point ({})".format(label))
        for lib in self.children:
            lib.calculate()

        # combine all libraries for a given timepoint
        self.logger.info("Aggregating SeqLib data")

        destination = "/main/{}/counts_unfiltered".format(label)
        if destination in self.store.keys():
            # need to remove the current destination table because we are using append
            # append is required because it takes the "min_itemsize" argument, and put doesn't
            self.logger.info("Replacing existing '{}'".format(destination))
            self.store.remove(destination)

        # seqlib count table name for this element type
        lib_table = "/main/{}/counts".format(label)

        # create an index of all elements in the analysis
        complete_index = pd.Index([])
        for tp in self.timepoints:
            for lib in self.libraries[tp]:
                complete_index = complete_index.union(
                    pd.Index(lib.store.select_column(lib_table, "index"))
                )
        self.logger.info(
            "Created shared index for count data ({} {})".format(
                len(complete_index), label
            )
        )

        # min_itemsize value
        max_index_length = complete_index.map(len).max()

        # perform operation in chunks
        tp_frame = None
        for i in xrange(0, len(complete_index), self.chunksize):
            # don't duplicate the index if the chunksize is large
            if self.chunksize < len(complete_index):
                index_chunk = complete_index[i : i + self.chunksize]
            else:
                index_chunk = complete_index
            self.logger.info(
                "Merging counts for chunk {} ({} rows)".format(
                    i / self.chunksize + 1, len(index_chunk)
                )
            )

            for tp in self.timepoints:
                c = self.libraries[tp][0].store.select(lib_table, "index = index_chunk")
                for lib in self.libraries[tp][1:]:
                    c = c.add(
                        lib.store.select(lib_table, "index = index_chunk"), fill_value=0
                    )
                c.columns = ["c_{}".format(tp)]
                if tp == 0:
                    tp_frame = c
                else:
                    tp_frame = tp_frame.join(c, how="outer")

            # save the unfiltered counts
            if "/main/{}/counts_unfiltered".format(label) not in self.store:
                self.store.append(
                    "/main/{}/counts_unfiltered".format(label),
                    tp_frame.astype(float),
                    min_itemsize={"index": max_index_length},
                    data_columns=list(tp_frame.columns),
                )
            else:
                self.store.append(
                    "/main/{}/counts_unfiltered".format(label), tp_frame.astype(float)
                )

    def filter_counts(self, label):
        """
        Converts unfiltered counts stored in ``/main/label/counts_unfiltered`` 
        into filtered counts calculated from complete cases (elements with a 
        non-zero count in each time point).

        For the most basic element type (variant or barcode, depending on the 
        experimental design), the result of this operation simply drops any 
        rows that have missing counts. For other element types, such as 
        synonymous variants, the counts are re-aggregated using only the 
        complete cases in the underlying element type.
        """
        if (self.is_barcodeid() or self.is_barcodevariant()) and label != "barcodes":
            # calculate proper combined counts
            # df = self.store.select("/main/barcodes/counts") # this should exist because of the order of label calculations
            # redo the barcode->variant/id mapping with the filtered counts
            # NOT YET IMPLEMENTED
            df = self.store.select(
                "/main/{}/counts_unfiltered".format(label)
            )  # just do this for now
        else:
            df = self.store.select("/main/{}/counts_unfiltered".format(label))
        df.dropna(axis="index", how="any", inplace=True)
        self.store.put(
            "/main/{}/counts".format(label),
            df.astype(float),
            format="table",
            data_columns=df.columns,
        )

    def combine_barcode_maps(self):
        if self.check_store("/main/barcodemap"):
            return

        bcm = None
        for lib in self.children:
            if bcm is None:
                bcm = lib.store["/raw/barcodemap"]
            else:
                bcm = bcm.join(
                    lib.store["/raw/barcodemap"], rsuffix=".drop", how="outer"
                )
                new = bcm.loc[pd.isnull(bcm)["value"]]
                bcm.loc[new.index, "value"] = new["value.drop"]
                bcm.drop("value.drop", axis="columns", inplace=True)
        bcm.sort_values("value", inplace=True)
        self.store.put(
            "/main/barcodemap", bcm, format="table", data_columns=bcm.columns
        )

    def calculate(self):
        """
        Wrapper method to calculate counts and enrichment scores 
        for all data in the :py:class:`~selection.Selection`.
        """
        if len(self.labels) == 0:
            raise ValueError(
                "No data present across all sequencing libraries [{}]".format(self.name)
            )
        for label in self.labels:
            self.merge_counts_unfiltered(label)
            self.filter_counts(label)
        if self.is_barcodevariant() or self.is_barcodeid():
            self.combine_barcode_maps()
        if self.scoring_method == "counts":
            pass
        elif self.scoring_method == "ratios":
            for label in self.labels:
                self.calc_ratios(label)
        elif self.scoring_method == "simple":
            for label in self.labels:
                self.calc_simple_ratios(label)
        elif self.scoring_method in ("WLS", "OLS"):
            if len(self.timepoints) <= 2:
                raise ValueError(
                    "Regression-based scoring requires three or more time points."
                )
            for label in self.labels:
                self.calc_log_ratios(label)
                if self.scoring_method == "WLS":
                    self.calc_weights(label)
                self.calc_regression(label)
        else:
            raise ValueError(
                'Invalid scoring method "{}" [{}]'.format(
                    self.scoring_method, self.name
                )
            )

        if self.scoring_method in ("ratios", "WLS", "OLS") and self.component_outliers:
            if self.is_barcodevariant() or self.is_barcodeid():
                self.calc_outliers("barcodes")
            if self.is_coding():
                self.calc_outliers("variants")

    def calc_simple_ratios(self, label):
        """
        Calculate simplified (original Enrich) ratios scores. This method does not produce standard errors.
        """
        if self.check_store("/main/{}/scores".format(label)):
            return

        self.logger.info("Calculating simple ratios ({})".format(label))
        c_last = "c_{}".format(self.timepoints[-1])
        df = self.store.select(
            "/main/{}/counts".format(label), "columns in ['c_0', c_last]"
        )

        # perform operations on the numpy values of the data frame for easier broadcasting
        ratios = (df[c_last].values.astype("float") / df[c_last].sum(axis="index")) / (
            df["c_0"].values.astype("float") / df["c_0"].sum(axis="index")
        )
        # make it a data frame again
        ratios = pd.DataFrame(data=ratios, index=df.index, columns=["ratio"])
        ratios["score"] = np.log2(ratios["ratio"])
        ratios["SE"] = np.nan
        ratios = ratios[["score", "SE", "ratio"]]  # re-order columns

        self.store.put(
            "/main/{}/scores".format(label),
            ratios,
            format="table",
            data_columns=ratios.columns,
        )

    def calc_ratios(self, label):
        """
        Calculate frequency ratios and standard errors between the last timepoint and the input.

        Ratios can be calculated using one of three methods:

        * wt
        * complete
        * full

        """
        if self.check_store("/main/{}/scores".format(label)):
            return

        self.logger.info("Calculating ratios ({})".format(label))
        c_last = "c_{}".format(self.timepoints[-1])
        df = self.store.select(
            "/main/{}/counts".format(label), "columns in ['c_0', c_last]"
        )

        if self.logr_method == "wt":
            if "variants" in self.labels:
                wt_label = "variants"
            elif "identifiers" in self.labels:
                wt_label = "identifiers"
            else:
                raise ValueError(
                    "Failed to use wild type log ratio method, suitable data table not present [{}]".format(
                        self.name
                    )
                )
            shared_counts = self.store.select(
                "/main/{}/counts".format(wt_label),
                "columns in ['c_0', c_last] & index='{}'".format(WILD_TYPE_VARIANT),
            )
            if len(shared_counts) == 0:  # wild type not found
                raise ValueError(
                    "Failed to use wild type log ratio method, wild type sequence not present [{}]".format(
                        self.name
                    )
                )
            shared_counts = shared_counts.values + 0.5
        elif self.logr_method == "complete":
            shared_counts = (
                self.store.select(
                    "/main/{}/counts".format(label), "columns in ['c_0', c_last]"
                )
                .sum(axis="index")
                .values
                + 0.5
            )
        elif self.logr_method == "full":
            shared_counts = (
                self.store.select(
                    "/main/{}/counts_unfiltered".format(label),
                    "columns in ['c_0', c_last]",
                )
                .sum(axis="index", skipna=True)
                .values
                + 0.5
            )
        else:
            raise ValueError(
                'Invalid log ratio method "{}" [{}]'.format(self.logr_method, self.name)
            )

        ratios = np.log(df[["c_0", c_last]].values + 0.5) - np.log(shared_counts)
        ratios = ratios[:, 1] - ratios[:, 0]  # selected - input

        ratios = pd.DataFrame(data=ratios, index=df.index, columns=["logratio"])

        shared_variance = np.sum(1.0 / shared_counts)
        ratios["variance"] = (
            np.sum(1.0 / (df[["c_0", c_last]].values + 0.5), axis=1) + shared_variance
        )

        ratios["score"] = ratios["logratio"]
        ratios["SE"] = np.sqrt(ratios["variance"])
        ratios = ratios[["score", "SE", "logratio", "variance"]]  # re-order columns

        self.store.put(
            "/main/{}/scores".format(label),
            ratios,
            format="table",
            data_columns=ratios.columns,
        )

    def calc_log_ratios(self, label):
        """
        Calculate the log ratios that will be fit using the linear model.
        """
        if self.check_store("/main/{}/log_ratios".format(label)):
            return

        self.logger.info("Calculating log ratios ({})".format(label))
        ratios = self.store.select("/main/{}/counts".format(label))
        index = ratios.index
        c_n = ["c_{}".format(x) for x in self.timepoints]
        ratios = np.log(ratios + 0.5)

        # perform operations on the numpy values of the data frame for easier broadcasting
        ratios = ratios[c_n].values
        if self.logr_method == "wt":
            if "variants" in self.labels:
                wt_label = "variants"
            elif "identifiers" in self.labels:
                wt_label = "identifiers"
            else:
                raise ValueError(
                    "Failed to use wild type log ratio method, suitable data table not present [{}]".format(
                        self.name
                    )
                )
            wt_counts = self.store.select(
                "/main/{}/counts".format(wt_label),
                "columns=c_n & index='{}'".format(WILD_TYPE_VARIANT),
            )
            if len(wt_counts) == 0:  # wild type not found
                raise ValueError(
                    "Failed to use wild type log ratio method, wild type sequence not present [{}]".format(
                        self.name
                    )
                )
            ratios = ratios - np.log(wt_counts.values + 0.5)
        elif self.logr_method == "complete":
            ratios = ratios - np.log(
                self.store.select("/main/{}/counts".format(label), "columns=c_n")
                .sum(axis="index")
                .values
                + 0.5
            )
        elif self.logr_method == "full":
            ratios = ratios - np.log(
                self.store.select(
                    "/main/{}/counts_unfiltered".format(label), "columns=c_n"
                )
                .sum(axis="index", skipna=True)
                .values
                + 0.5
            )
        else:
            raise ValueError(
                'Invalid log ratio method "{}" [{}]'.format(self.logr_method, self.name)
            )

        # make it a data frame again
        ratios = pd.DataFrame(
            data=ratios,
            index=index,
            columns=["L_{}".format(x) for x in self.timepoints],
        )
        self.store.put(
            "/main/{}/log_ratios".format(label),
            ratios,
            format="table",
            data_columns=ratios.columns,
        )

    def calc_weights(self, label):
        """
        Calculate the regression weights (1 / variance).
        """
        if self.check_store("/main/{}/weights".format(label)):
            return

        self.logger.info("Calculating regression weights ({})".format(label))
        variances = self.store.select("/main/{}/counts".format(label))
        c_n = ["c_{}".format(x) for x in self.timepoints]
        index = variances.index
        # perform operations on the numpy values of the data frame for easier broadcasting
        # variances = 1.0 / (variances[c_n].values + 0.5) + 1.0 / (variances[['c_0']].values + 0.5)
        variances = 1.0 / (variances[c_n].values + 0.5)
        if self.logr_method == "wt":
            if "variants" in self.labels:
                wt_label = "variants"
            elif "identifiers" in self.labels:
                wt_label = "identifiers"
            else:
                raise ValueError(
                    "Failed to use wild type log ratio method, suitable data table not present [{}]".format(
                        self.name
                    )
                )
            wt_counts = self.store.select(
                "/main/{}/counts".format(wt_label),
                "columns=c_n & index='{}'".format(WILD_TYPE_VARIANT),
            )
            if len(wt_counts) == 0:  # wild type not found
                raise ValueError(
                    "Failed to use wild type log ratio method, wild type sequence not present [{}]".format(
                        self.name
                    )
                )
            variances = variances + 1.0 / (wt_counts.values + 0.5)
        elif self.logr_method == "complete":
            variances = variances + 1.0 / (
                self.store.select("/main/{}/counts".format(label), "columns=c_n")
                .sum(axis="index")
                .values
                + 0.5
            )
        elif self.logr_method == "full":
            variances = variances + 1.0 / (
                self.store.select(
                    "/main/{}/counts_unfiltered".format(label), "columns=c_n"
                )
                .sum(axis="index", skipna=True)
                .values
                + 0.5
            )
        else:
            raise ValueError(
                'Invalid log ratio method "{}" [{}]'.format(self.logr_method, self.name)
            )
        variances = 1.0 / variances  # weights are reciprocal of variances
        # make it a data frame again
        variances = pd.DataFrame(
            data=variances,
            index=index,
            columns=["W_{}".format(x) for x in self.timepoints],
        )
        self.store.put(
            "/main/{}/weights".format(label),
            variances,
            format="table",
            data_columns=variances.columns,
        )

    def calc_regression(self, label):
        """
        Calculate least squares regression for *label*. If *weighted* is ``True``, calculates weighted least squares; else ordinary least squares.

        Regression results are stored in ``'/main/label/scores'``

        """
        if self.check_store("/main/{}/scores".format(label)):
            return
        elif "/main/{}/scores".format(label) in self.store.keys():
            # need to remove the current keys because we are using append
            self.store.remove("/main/{}/scores".format(label))

        self.logger.info(
            "Calculating {} regression coefficients ({})".format(
                self.scoring_method, label
            )
        )
        # append is required because it takes the "min_itemsize" argument, and put doesn't
        longest = (
            self.store.select("/main/{}/log_ratios".format(label), "columns='index'")
            .index.map(len)
            .max()
        )
        chunk = 1
        if self.scoring_method == "WLS":
            for data in self.store.select_as_multiple(
                ["/main/{}/log_ratios".format(label), "/main/{}/weights".format(label)],
                chunksize=self.chunksize,
            ):
                self.logger.info(
                    "Calculating weighted least squares for chunk {} ({} rows)".format(
                        chunk, len(data.index)
                    )
                )
                result = data.apply(
                    regression_apply, args=[self.timepoints, True], axis="columns"
                )
                self.store.append(
                    "/main/{}/scores".format(label),
                    result,
                    min_itemsize={"index": longest},
                )
                chunk += 1
        elif self.scoring_method == "OLS":
            for data in self.store.select(
                "/main/{}/log_ratios".format(label), chunksize=self.chunksize
            ):
                self.logger.info(
                    "Calculating ordinary least squares for chunk {} ({} rows)".format(
                        chunk, len(data.index)
                    )
                )
                result = data.apply(
                    regression_apply, args=[self.timepoints, False], axis="columns"
                )
                self.store.append(
                    "/main/{}/scores".format(label),
                    result,
                    min_itemsize={"index": longest},
                )
                chunk += 1
        else:
            raise ValueError(
                'Invalid regression scoring method "{}" [{}]'.format(
                    self.scoring_method, self.name
                )
            )

        # need to read from the file, calculate percentiles, and rewrite it
        self.logger.info(
            "Calculating slope standard error percentiles ({})".format(label)
        )
        data = self.store["/main/{}/scores".format(label)]
        data["score"] = data["slope"]
        data["SE"] = data["SE_slope"]
        data["SE_pctile"] = [
            stats.percentileofscore(data["SE"], x, "weak") for x in data["SE"]
        ]
        data = data[
            [
                "score",
                "SE",
                "SE_pctile",
                "slope",
                "intercept",
                "SE_slope",
                "t",
                "pvalue_raw",
            ]
        ]  # reorder columns
        self.store.put(
            "/main/{}/scores".format(label),
            data,
            format="table",
            data_columns=data.columns,
        )

    def wt_plot(self, pdf):
        """
        Create a plot of the linear fit of the wild type variant.

        *pdf* is an open PdfPages instance.

        Only created for selections that use WLS or OLS scoring and have a wild type specified. 
        Uses :py:func:`~plots.fit_axes` for the plotting.
        """
        self.logger.info("Creating wild type fit plot")

        # get the data and calculate log ratios
        if "variants" in self.labels:
            wt_label = "variants"
        elif "identifiers" in self.labels:
            wt_label = "identifiers"
        else:
            raise ValueError(
                "Invalid selection type for wild type fit plot [{}]".format(self.name)
            )
        data = self.store.select(
            "/main/{}/counts".format(wt_label),
            where='index = "{}"'.format(WILD_TYPE_VARIANT),
        ).ix[0]
        sums = self.store["/main/{}/counts".format(wt_label)].sum(
            axis="index"
        )  # sum of complete cases (N')
        yvalues = np.log(data + 0.5) - np.log(sums + 0.5)
        xvalues = [tp / float(max(self.timepoints)) for tp in self.timepoints]

        # fit the line
        X = sm.add_constant(xvalues)  # fit intercept
        if self.scoring_method == "WLS":
            W = 1 / (1 / (data + 0.5) + 1 / (sums + 0.5))
            fit = sm.WLS(yvalues, X, weights=W).fit()
        elif self.scoring_method == "OLS":
            fit = sm.OLS(yvalues, X).fit()
        else:
            raise ValueError(
                'Invalid regression scoring method "{}" [{}]'.format(
                    self.scoring_method, self.name
                )
            )
        intercept, slope = fit.params
        slope_se = fit.bse["x1"]

        # make the plot
        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        fit_axes(ax, xvalues, yvalues, slope, intercept, xlabels=self.timepoints)
        fit_axes_text(ax, cornertext="Slope {:3.2f}\nSE {:.1f}".format(slope, slope_se))
        ax.set_title("Wild Type Shape\n{}".format(self.name))
        ax.set_ylabel("Log Ratio (Complete Cases)")

        pdf.savefig(fig)
        plt.close(fig)

    def se_pctile_plot(self, label, pdf):
        """
        Create plots of the linear fit of 21 selected variants, evenly distributed based on their standard error percentiles (0, 5, 10, ..., 100).

        *label* is the data label (barcode, variant, etc.)

        *pdf* is an open PdfPages instance.

        Uses :py:func:`~plots.fit_axes` for the plotting.
        """
        self.logger.info("Creating representative fit plots ({})".format(label))

        se_data = self.store.select(
            "/main/{}/scores".format(label),
            where="columns in ['slope', 'intercept', 'SE_pctile'] & index!='{}' & index!='{}'".format(
                WILD_TYPE_VARIANT, SYNONYMOUS_VARIANT
            ),
        )
        se_data.sort_values("SE_pctile", inplace=True)
        se_data.dropna(axis="index", how="any", inplace=True)

        indices = np.linspace(0, len(se_data.index) - 1, 21).astype("int")
        se_data = se_data.ix[indices]

        # retrieves the whole DF because one case was hanging when trying to use select
        # totally unexplained, should fix later
        ratio_data = self.store.select("/main/{}/log_ratios".format(label)).loc[
            se_data.index
        ]

        fig, axarr = plt.subplots(7, 3, sharex="all", sharey="all")
        fig.subplots_adjust(
            hspace=0, wspace=0
        )  # eliminate white space between the subplots
        fig.set_size_inches((10, 17))
        fig.suptitle(
            "Representative {} Fits\n{} ({})".format(
                self.scoring_method, self.name, label.title()
            )
        )
        fig.subplots_adjust(top=0.958)  # eliminate white space after the title

        # turn off all tick labels
        plt.setp([ax.get_xticklabels() for ax in axarr.reshape(-1)], visible=False)
        plt.setp([ax.get_yticklabels() for ax in axarr.reshape(-1)], visible=False)

        # tick labels back on for some plots
        plt.setp(
            [ax.get_xticklabels() for ax in (axarr[-1, 0], axarr[-1, 2])], visible=True
        )
        plt.setp(
            [
                ax.get_yticklabels()
                for ax in (axarr[0, 0], axarr[2, 0], axarr[4, 0], axarr[6, 0])
            ],
            visible=True,
        )
        plt.setp(
            [ax.get_yticklabels() for ax in (axarr[1, -1], axarr[3, -1], axarr[5, -1])],
            visible=True,
        )
        plt.setp(
            [ax.yaxis for ax in (axarr[1, -1], axarr[3, -1], axarr[5, -1])],
            ticks_position="right",
        )

        # create the fit plots and add text to the individual plots
        for i, ax in enumerate(axarr.reshape(-1)):
            index = se_data.index[i]
            fit_axes(
                ax,
                xvalues=[x / float(max(self.timepoints)) for x in self.timepoints],
                yvalues=ratio_data.loc[index],
                slope=se_data.loc[index, "slope"],
                intercept=se_data.loc[index, "intercept"],
                xlabels=self.timepoints,
            )
            fit_axes_text(
                ax,
                cornertext="Slope {:3.2f}\nSE Pctile {:.1f}".format(
                    se_data.loc[index, "slope"], se_data.loc[index, "SE_pctile"]
                ),
                centertext=index,
            )

        # turn off the subplot axis labels
        [ax.set_xlabel("") for ax in axarr.reshape(-1)]
        [ax.set_ylabel("") for ax in axarr.reshape(-1)]

        # add x and y labels in the middle
        axarr[-1, 1].set_xlabel("Time Point")
        axarr[3, 0].set_ylabel("Log Ratio")

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    def timepoint_counts_plot(self, label, pdf):
        """
        Create barplot of the number of items counted for each time point.

        *label* is the data label (barcode, variant, etc.)

        *pdf* is an open PdfPages instance.
        """
        self.logger.info("Creating time point count plots ({})".format(label))

        counts = self.store["/main/{}/counts".format(label)].sum(axis="index")

        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        configure_axes(ax)

        xpos = np.arange(len(self.timepoints))
        width = 0.8
        ax.bar(xpos, counts, width, color=plot_colors["bright1"])
        ax.set_title("Total {}\n{}".format(label.title(), self.name))
        ax.set_ylabel("Count")
        ax.set_xlabel("Timepoint")
        ax.set_xticks(xpos + width / 2.0)
        ax.set_xticklabels(self.timepoints)

        pdf.savefig(fig)
        plt.close(fig)

    def volcano_plot(self, label, pdf, colors="YlGnBu_r", log_bins=True):
        """
        Create a volcano plot (p-value vs. functional score).

        *label* is the data label (barcode, variant, etc.)

        *pdf* is an open PdfPages instance.

        The p-values used are the regression p-values (p-value of non-zero slope). Due to the large number of points, we use a hexbin plot showing the density instead of a scatter plot.
        """
        self.logger.info("Creating volcano plot ({})".format(label))

        # get the data
        data = self.store.select(
            "/main/{}/scores".format(label), "columns=['score', 'pvalue_raw']"
        )
        volcano_plot(
            data,
            pdf,
            title="{} ({})".format(self.name, label.title()),
            colors=colors,
            log_bins=log_bins,
        )

    def make_plots(self):
        """
        Create plots for this entity.

        This function handles opening and closing the various PDF files for multi-page plots, as well as plotting similar plots for different data labels.
        """
        if self.plots_requested:
            self.logger.info("Creating plots")

            # counts per time point
            pdf = PdfPages(os.path.join(self.plot_dir, "timepoint_counts.pdf"))
            for label in self.labels:
                self.timepoint_counts_plot(label, pdf)
            pdf.close()

            # wild type shape
            if self.logr_method == "wt" and self.scoring_method in ("WLS", "OLS"):
                pdf = PdfPages(os.path.join(self.plot_dir, "wt_shape.pdf"))
                self.wt_plot(pdf)
                pdf.close()

            # regression weights
            if self.scoring_method == "WLS":
                pdf = PdfPages(os.path.join(self.plot_dir, "regression_weights.pdf"))
                for label in self.labels:
                    weights_plot(self, label, pdf)
                pdf.close()

            # linear fits by standard error percentile
            if self.scoring_method in ("WLS", "OLS"):
                pdf = PdfPages(os.path.join(self.plot_dir, "se_pctile.pdf"))
                for label in self.labels:
                    self.se_pctile_plot(label, pdf)
                pdf.close()

            # volcano plots
            # if self.scoring_method in ("WLS", "OLS", "ratios"):
            if self.scoring_method in ("WLS", "OLS") and "variants" in self.labels:
                pdf = PdfPages(os.path.join(self.plot_dir, "volcano.pdf"))
                for label in self.labels:
                    self.volcano_plot(label, pdf, log_bins=True)
                pdf.close()

            # library diversity for each time point (amino acid)
            if "synonymous" in self.labels:
                pdf = PdfPages(os.path.join(self.plot_dir, "diversity_aa.pdf"))
                for tp in self.timepoints:
                    self.sfmap_wrapper(
                        cname="c_{}".format(tp), pdf=pdf, coding=True, log10=True
                    )
                pdf.close()

            # library diversity for each time point (nucleotide)
            if "variants" in self.labels:
                pdf = PdfPages(os.path.join(self.plot_dir, "diversity_nt.pdf"))
                for tp in self.timepoints:
                    self.sfmap_wrapper(
                        cname="c_{}".format(tp), pdf=pdf, coding=False, log10=True
                    )
                pdf.close()

            # sequence-function maps
            if self.scoring_method != "counts":
                if "synonymous" in self.labels:
                    pdf = PdfPages(
                        os.path.join(self.plot_dir, "sequence_function_map_aa.pdf")
                    )
                    self.sfmap_wrapper(cname="score", pdf=pdf, coding=True)
                    pdf.close()
                if "variants" in self.labels:
                    pdf = PdfPages(
                        os.path.join(self.plot_dir, "sequence_function_map_nt.pdf")
                    )
                    self.sfmap_wrapper(cname="score", pdf=pdf, coding=False)
                    pdf.close()

        # SeqLib plots
        for lib in self.children:
            lib.make_plots()

    def write_tsv(self):
        """
        Write each table from the store to its own tab-separated file.

        Files are written to a ``tsv`` directory in the default output location. 
        File names are the HDF5 key with ``'_'`` substituted for ``'/'``.
        """
        if self.tsv_requested:
            self.logger.info("Generating tab-separated output files")
            for k in self.store.keys():
                self.write_table_tsv(k)
        for lib in self.children:
            lib.write_tsv()

    def synonymous_variants(self):
        """
        Populate and return a dictionary mapping synonymous variants to the 
        list of associated variants in ``/main/variants/counts``.
        """
        mapping = dict()
        try:
            variants = self.store.select_column("/main/variants/counts", "index")
        except KeyError:
            raise KeyError("No variant counts found [{}]".format(self.name))
        for v in variants:
            pv = protein_variant(v)
            try:
                mapping[pv].append(v)
            except KeyError:
                mapping[pv] = [v]
        return mapping

    def sfmap_wrapper(self, cname, pdf, coding, log10=False):
        """
        Create a sequence function map for either scores or library diversity.

        Uses :py:func:`~sfmap.sfmap_plot` for the plotting.
        """
        plot_options = self.get_root().plot_options

        if cname.startswith("c_"):
            counts = True
        elif cname == "score":
            counts = False
        else:
            raise ValueError("Invalid sequence-function map data column.")

        if coding:
            label = "amino acid"
        else:
            label = "nucleotide"

        if counts:
            self.logger.info("Creating diversity map ({})".format(label))
        else:
            self.logger.info("Creating sequence-function map ({})".format(label))

        # build the data frame name and get the data
        df_name = "/main/"
        if coding:
            df_name += "synonymous/"
        else:
            df_name += "variants/"
        if counts:
            df_name += "counts_unfiltered"
        else:
            df_name += "scores"
        if plot_options is not None:
            data, wtseq = singleton_dataframe(
                self.store[df_name][cname],
                self.wt,
                coding=coding,
                aa_list=plot_options["aa_list"],
            )
        else:
            data, wtseq = singleton_dataframe(
                self.store[df_name][cname], self.wt, coding=coding
            )
        if counts:
            data_se = None
        else:
            if plot_options is not None:
                data_se, _ = singleton_dataframe(
                    self.store[df_name]["SE"],
                    self.wt,
                    coding=coding,
                    aa_list=plot_options["aa_list"],
                )
            else:
                data_se, _ = singleton_dataframe(
                    self.store[df_name]["SE"], self.wt, coding=coding
                )

        # format the title
        if coding:
            title = "Amino Acid"
        else:
            title = "Nucleotide"
        if counts:
            title += " Diversity Map\n{} (Time {})".format(
                self.name, cname[2:]
            )  # trim the "c_"
        else:
            if self.scoring_method in ("WLS", "OLS"):
                title += " Sequence-Function Map\n{} ({} Slope)".format(
                    self.name, self.scoring_method
                )
            elif self.scoring_method == "ratios":
                title += " Sequence-Function Map\n{} ({})".format(
                    self.name, "Enrich2 Ratio"
                )
            elif self.scoring_method == "simple":
                title += " Sequence-Function Map\n{} ({})".format(
                    self.name, "Simplified Ratio"
                )
            else:
                raise ValueError("Invalid scoring method", self.name)

        if counts and log10:
            style = "logcounts"
        elif counts:
            style = "counts"
        else:
            style = "scores"

        if plot_options is not None:
            sfmap_plot(
                df=data,
                pdf=pdf,
                style=style,
                df_se=data_se,
                dimensions="tall",
                wt=wtseq,
                title=title,
                aa_list=plot_options["aa_list"],
                aa_label_groups=plot_options["aa_label_groups"],
            )
        else:
            sfmap_plot(
                df=data,
                pdf=pdf,
                style=style,
                df_se=data_se,
                dimensions="tall",
                wt=wtseq,
                title=title,
            )

    def barcodemap_mapping(self):
        mapping = dict()
        for bc, v in self.store["/main/barcodemap"].iterrows():
            v = v["value"]
            try:
                mapping[v].update([bc])
            except KeyError:
                mapping[v] = {bc}
        return mapping

    def calc_outliers(self, label, minimum_components=4, log_chunksize=20000):
        """
        Test whether an element's individual components have significantly different 
        scores from the element. Results are stored in ``'/main/<label>/outliers'``.

        Args:
            label (str): label for the component (``'variants'`` or ``'barcodes'``)

            minimum_components (int): minimum number of componenents required for any statistics to be calculated

            log_chunksize (int): will output a log message every *n* rows

        """
        if self.check_store("/main/{}/outliers".format(label)):
            return

        if label == "variants":
            label2 = "synonymous"
        elif label == "barcodes":
            if self.is_barcodevariant():
                label2 = "variants"
            elif self.is_barcodeid():
                label2 = "identifiers"
            else:
                # this should never happen
                raise AttributeError(
                    "Failed to identify parent label when calculating barcode outliers [{}]".format(
                        self.name
                    )
                )
        else:
            raise KeyError(
                "Invalid label '{}' for calc_outliers [{}]".format(label, self.name)
            )

        self.logger.info("Identifying outliers ({}-{})".format(label, label2))

        self.logger.info("Mapping {} to {}".format(label, label2))
        if label == "variants":
            mapping = self.synonymous_variants()
        elif label == "barcodes":
            mapping = self.barcodemap_mapping()
        else:
            raise KeyError(
                "Invalid label '{}' for calc_outliers [{}]".format(label, self.name)
            )

        # get the scores
        df1 = self.store.select(
            "/main/{}/scores".format(label), "columns=['score', 'SE']"
        )
        df2 = self.store.select(
            "/main/{}/scores".format(label2), "columns=['score', 'SE']"
        )

        # pre-calculate variances
        df1["var"] = df1["SE"] ** 2
        df1.drop("SE", axis=1, inplace=True)
        df2["var"] = df2["SE"] ** 2
        df2.drop("SE", axis=1, inplace=True)

        # set up empty results DF
        result_df = pd.DataFrame(
            np.nan, index=df1.index, columns=["z", "pvalue_raw", "parent"]
        )

        # because this step can be slow, output chunk-style logging messages
        # pre-calculate the lengths for the log messages
        log_chunksize_list = [log_chunksize] * (len(df2) / log_chunksize) + [
            len(df2) % log_chunksize
        ]
        log_chunk = 1

        for i, x in enumerate(df2.index):
            if i % log_chunksize == 0:
                self.logger.info(
                    "Calculating outlier p-values for chunk {} ({} rows) ({}-{})".format(
                        log_chunk, log_chunksize_list[log_chunk - 1], label, label2
                    )
                )
                log_chunk += 1
            try:
                components = df1.loc[mapping[x]].dropna(axis="index", how="all")
            except KeyError:
                # none of the components were in the index
                continue
            if len(components.index) >= minimum_components:
                for c in components.index:
                    zvalue = np.absolute(
                        df2.loc[x, "score"] - df1.loc[c, "score"]
                    ) / np.sqrt(df2.loc[x, "var"] + df1.loc[c, "var"])
                    result_df.loc[c, "z"] = zvalue
                    result_df.loc[c, "pvalue_raw"] = 2 * stats.norm.sf(zvalue)
                    result_df.loc[c, "parent"] = x
        if WILD_TYPE_VARIANT in result_df.index:
            result_df.loc[WILD_TYPE_VARIANT, "z"] = np.nan
            result_df.loc[WILD_TYPE_VARIANT, "pvalue_raw"] = np.nan
        result_df["z"] = result_df["z"].astype(float)
        result_df["pvalue_raw"] = result_df["pvalue_raw"].astype(float)

        self.store.put(
            "/main/{}/outliers".format(label),
            result_df,
            format="table",
            data_columns=result_df.columns,
        )
