#  Copyright 2016-2017 Alan F Rubin
#
#  This file is part of Enrich2.
#
#  Enrich2 is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Enrich2 is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Enrich2.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function
import logging
import pandas as pd
import numpy as np
import scipy.stats as stats
import re
from matplotlib.backends.backend_pdf import PdfPages
import os.path
from .storemanager import StoreManager
from .condition import Condition
from .constants import WILD_TYPE_VARIANT, SYNONYMOUS_VARIANT
from .sfmap import sfmap_plot
from .dataframe import singleton_dataframe
from .random_effects import rml_estimator


class Experiment(StoreManager):
    """
    Class for a coordinating multiple :py:class:`~.selection.Selection` 
    objects. Creating an 
    :py:class:`~experiment.Experiment` requires a valid *config* object, 
    usually from a ``.json`` configuration file.
    """

    store_suffix = "exp"
    treeview_class_name = "Experiment"

    def __init__(self):
        StoreManager.__init__(self)
        self.conditions = list()
        self._wt = None


    @property
    def wt(self):
        if self.has_wt_sequence():
            if self._wt is None:
                self._wt = self.selection_list()[0].wt.duplicate(self.name)
            return self._wt
        else:
            if self._wt is not None:
                raise ValueError("Experiment should not contain wild type sequence [{}]".format(self.name))
            else:
                return None


    def configure(self, cfg, configure_children=True):
        """
        Set up the :py:class:`~experiment.Experiment` using the *cfg* object, 
        usually from a ``.json`` configuration file.
        """
        StoreManager.configure(self, cfg)        
        if configure_children:
            if 'conditions' not in cfg:
                raise KeyError("Missing required config value {} [{}]".format('conditions', self.name))
            for cnd_cfg in cfg['conditions']:
                cnd = Condition()
                cnd.configure(cnd_cfg)
                self.add_child(cnd)

        selection_names = [x.name for x in self.selection_list()]
        if len(set(selection_names)) != len(selection_names):
            raise ValueError("Non-unique selection names [{}]".format(self.name))


    def serialize(self):
        """
        Format this object (and its children) as a config object suitable for dumping to a config file.
        """
        cfg = StoreManager.serialize(self)
        cfg['conditions'] = [child.serialize() for child in self.children]
        return cfg


    def _children(self):
        """
        Method bound to the ``children`` property. Returns a list of all 
        :py:class:`~condition.Condition` objects belonging to this object, 
        sorted by name.
        """
        return sorted(self.conditions, key=lambda x: x.name)


    def add_child(self, child):
        """
        Add a selection.
        """
        if child.name in self.child_names():
            raise ValueError("Non-unique condition name '{}' [{}]".format(child.name, self.name))
        child.parent = self
        self.conditions.append(child)


    def remove_child_id(self, tree_id):
        """
        Remove the reference to a :py:class:`~condition.Condition` with 
        Treeview id *tree_id*.
        """
        self.conditions = [x for x in self.conditions if x.treeview_id != tree_id]


    def selection_list(self):
        """
        Return the :py:class:`~selection.Selection` objects as a list.
        """
        selections = list()
        for cnd in self.children:
            selections.extend(cnd.children)
        return selections


    def validate(self):
        """
        Calls validate on all child Conditions. Also checks the wild type sequence status.
        """
        # check the wild type sequences
        if self.has_wt_sequence():
            for child in self.selection_list()[1:]:
                if self.selection_list()[0].wt != child.wt:
                    logging.warning("Inconsistent wild type sequences", extra={'oname' : self.name})
                    break

        for child in self.children:
            child.validate()


    def is_coding(self):
        """
        Return ``True`` if the all :py:class:`~selection.Selection` in the 
        :py:class:`~experiment.Experiment` count protein-coding variants, else 
        ``False``.
        """
        return all(x.is_coding() for x in self.selection_list())


    def has_wt_sequence(self):
        """
        Return ``True`` if the all :py:class:`~selection.Selection` in the 
        :py:class:`~experiment.Experiment` have a wild type sequence, else 
        ``False``.
        """
        return all(x.has_wt_sequence() for x in self.selection_list())


    def calculate(self):
        """
        Calculate scores for all :py:class:`~selection.Selection` objects.
        """
        if len(self.labels) == 0:
            raise ValueError("No data present across all conditions [{}]".format(self.name))
        for s in self.selection_list():
            s.calculate()
        self.combine_barcode_maps()
        for label in self.labels:
            self.calc_counts(label)
            if self.scoring_method != "counts":
                self.calc_shared_full(label)
                self.calc_shared(label)
                self.calc_scores(label)
                if label != "barcodes":
                    self.calc_pvalues_wt(label)


    def combine_barcode_maps(self):
        """
        Combine all barcode maps for :py:class:`~selection.Selection` objects into a single 
        data frame and store it in ``'/main/barcodemap'``.

        If multiple variants or IDs map to the same barcode, only the first one will be present
        in the barcode map table.

        The ``'/main/barcodemap'`` table is not created if no :py:class:`~selection.Selection` 
        has barcode map information.
        """
        if self.check_store("/main/barcodemap"):
            return

        bcm = None
        for sel in self.selection_list():
            if '/main/barcodemap' in sel.store.keys():
                if bcm is None:
                    bcm = sel.store['/main/barcodemap']
                else:
                    bcm = bcm.join(sel.store['/main/barcodemap'], rsuffix=".drop", how="outer")
                    new = bcm.loc[pd.isnull(bcm)['value']]
                    bcm.loc[new.index, 'value'] = new['value.drop']
                    bcm.drop("value.drop", axis="columns", inplace=True)
        if bcm is not None:
            bcm.sort_values("value", inplace=True)
            self.store.put("/main/barcodemap", bcm, format="table", data_columns=bcm.columns)


    def calc_counts(self, label):
        """
        Create a data frame of all counts in this Experiment. This data frame is not used for any calculations, but is provided to facilitate exploration of the data set.
        """
        if self.check_store("/main/{}/counts".format(label)):
            return

        # create columns multi-index
        # has to be lex-sorted for multi-slicing to work
        logging.info("Creating column multi-index for counts ({})".format(label), extra={'oname' : self.name})
        conditions_index = list()
        selections_index = list()
        values_index = list()

        for cnd in self.children:
            for sel in cnd.children:
                conditions_index.extend([cnd.name] * len(sel.timepoints))
                selections_index.extend([sel.name] * len(sel.timepoints))
                values_index.extend(["c_{}".format(x) for x in sorted(sel.timepoints)])
        columns = pd.MultiIndex.from_tuples(zip(conditions_index, selections_index, values_index), names=["condition", "selection", "timepoint"])

        # create union index
        logging.info("Creating row index for counts ({})".format(label), extra={'oname' : self.name})
        combined = None
        first = True
        for s in self.selection_list():
            if first:
                combined = s.store.select("/main/{}/counts_unfiltered".format(label), "columns='index'").index
                first = False
            else:
                combined = combined.join(s.store.select("/main/{}/counts_unfiltered".format(label), "columns='index'").index, how="outer")

        # create and fill the data frames
        logging.info("Populating Experiment data frame with counts ({})".format(label), extra={'oname' : self.name})
        data = pd.DataFrame(index=combined, columns=columns)
        for cnd in self.children:
            for sel in cnd.children:
                sel_data = sel.store.select("/main/{}/counts_unfiltered".format(label))
                for tp in sel.timepoints:
                    data.loc[:][cnd.name, sel.name, "c_{}".format(tp)] = sel_data["c_{}".format(tp)]
        self.store.put("/main/{}/counts".format(label), data, format="table")


    def calc_shared_full(self, label):
        """
        Use joins to create a data frame containing all scores across all Selections in the Experiment.
        """
        if self.check_store("/main/{}/scores_shared_full".format(label)):
            return

        # create columns multi-index
        # has to be lex-sorted for multi-slicing to work
        logging.info("Creating column multi-index for scores ({})".format(label), extra={'oname' : self.name})
        conditions_index = list()
        selections_index = list()
        values_index = list()

        if self.scoring_method == "simple":
            values_list = ["score"]
        else:
            values_list = ["score", "SE"]

        for cnd in self.children:
            for sel in cnd.children:
                conditions_index.extend([cnd.name] * len(values_list))
                selections_index.extend([sel.name] * len(values_list))
                values_index.extend(sorted(values_list))
        columns = pd.MultiIndex.from_tuples(zip(conditions_index, selections_index, values_index), names=["condition", "selection", "value"])

        # create union index
        logging.info("Creating row index for scores ({})".format(label), extra={'oname' : self.name})
        combined = None
        first = True
        for s in self.selection_list():
            if first:
                combined = s.store.select("/main/{}/scores".format(label), "columns='index'").index
                first = False
            else:
                combined = combined.join(s.store.select("/main/{}/scores".format(label), "columns='index'").index, how="outer")

        # create and fill the data frames
        logging.info("Populating Experiment data frame with scores ({})".format(label), extra={'oname' : self.name})
        data = pd.DataFrame(index=combined, columns=columns)
        for cnd in self.children:
            for sel in cnd.children:
                sel_data = sel.store.select("/main/{}/scores".format(label))
                for v in values_list:
                    data.loc[:][cnd.name, sel.name, v] = sel_data[v]
        self.store.put("/main/{}/scores_shared_full".format(label), data, format="table")


    def calc_shared(self, label):
        """
        Get the subset of scores that are shared across all Selections in each Condition.
        """
        if self.check_store("/main/{}/scores_shared".format(label)):
            return

        idx = pd.IndexSlice

        logging.info("Identifying subset shared across all Selections ({})".format(label), extra={'oname' : self.name})
        data = self.store.select("/main/{}/scores_shared_full".format(label))

        # identify variants that are found in all selections in at least one condition
        complete = np.full(len(data.index), False, dtype=bool)
        for cnd in data.columns.levels[0]:
            complete = np.logical_or(complete, data.loc[:, idx[cnd, :, :]].notnull().all(axis='columns'))

        data = data.loc[complete]

        self.store.put("/main/{}/scores_shared".format(label), data, format="table")


    def calc_scores(self, label):
        """
        Combine the scores and standard errors within each condition.
        """
        if self.check_store("/main/{}/scores".format(label)):
            return

        logging.info("Calculating per-condition scores ({})".format(label), extra={'oname' : self.name})

        # set up new data frame
        shared_index = self.store.select("/main/{}/scores_shared".format(label), "columns='index'").index
        columns = pd.MultiIndex.from_product([sorted(self.child_names()), sorted(["score", "SE", "epsilon"])], names=["condition", "value"])
        data = pd.DataFrame(np.nan, index=shared_index, columns=columns)
        del shared_index
        del columns

        # set up local variables
        idx = pd.IndexSlice

        score_df = self.store.select("/main/{}/scores_shared".format(label))
        if self.scoring_method == "simple":
            # special case for simple ratios that have no SE
            # calculates the average score
            for cnd in score_df.columns.levels[0]:
                data.loc[:, idx[cnd, 'score']] = score_df.loc[:, idx[cnd, :, 'score']].mean(axis=1)
        else:
            for cnd in score_df.columns.levels[0]:
                y = np.array(score_df.loc[:, idx[cnd, :, 'score']].values).T
                sigma2i = np.array(score_df.loc[:, idx[cnd, :, 'SE']].values **2).T

                # single replicate of the condition
                if y.shape[0] == 1:
                    data.loc[:, idx[cnd, 'score']] = y.ravel()
                    data.loc[:, idx[cnd, 'SE']] = np.sqrt(sigma2i).ravel()
                    data.loc[:, idx[cnd, 'epsilon']] = 0.
                # multiple replicates
                else:
                    betaML, sigma2ML, eps = rml_estimator(y, sigma2i)
                    data.loc[:, idx[cnd, 'score']] = betaML
                    data.loc[:, idx[cnd, 'SE']] = np.sqrt(sigma2ML)
                    data.loc[:, idx[cnd, 'epsilon']] = eps

                # special case for wild type variant and wild type normalization
                if self.logr_method == "wt" and WILD_TYPE_VARIANT in data.index:
                    data.loc[WILD_TYPE_VARIANT, idx[:, 'SE']] = 0.
                    data.loc[WILD_TYPE_VARIANT, idx[:, 'score']] = 0.
                    data.loc[WILD_TYPE_VARIANT, idx[:, 'epsilon']] = 0.

        # store the data
        self.store.put("/main/{}/scores".format(label), data, format="table")


    def calc_pvalues_wt(self, label):
        """
        Calculate uncorrected pvalue for each variant compared to wild type.
        """
        if self.check_store("/main/{}/scores_pvalues_wt".format(label)):
            return

        idx = pd.IndexSlice

        wt = self.store.select("/main/{}/scores".format(label), "index=WILD_TYPE_VARIANT")
        if len(wt) == 0:    # no wild type score
            logging.info("Failed to find wild type score, skipping wild type p-value calculations", extra={'oname' : self.name})
            return
        data = self.store.select("/main/{}/scores".format(label), "index!=WILD_TYPE_VARIANT")

        columns = pd.MultiIndex.from_product([sorted(self.child_names()), sorted(["z", "pvalue_raw"])], names=["condition", "value"])
        result_df = pd.DataFrame(index=data.index, columns=columns)

        condition_labels = data.columns.levels[0]
        for cnd in condition_labels:
            result_df.loc[:, idx[cnd, 'z']] = \
                np.absolute(wt.loc[WILD_TYPE_VARIANT, idx[cnd, 'score']] - data.loc[:, idx[cnd, 'score']]) / np.sqrt(wt.loc[WILD_TYPE_VARIANT, idx[cnd, 'SE']] ** 2 + data.loc[:, idx[cnd, 'SE']] ** 2)
            result_df.loc[:, idx[cnd, 'pvalue_raw']] = 2 * stats.norm.sf(result_df.loc[:, idx[cnd, 'z']])

        self.store.put("/main/{}/scores_pvalues_wt".format(label), result_df, format="table")


    def calc_pvalues_pairwise(self, label):
        """
        Calculate pvalues for each variant in each pair of Conditions.
        """
        if self.check_store("/main/{}/scores_pvalues".format(label)):
            return

        data = self.store['/main/{}/scores'.format(label)]

        cnd1_index = list()
        cnd2_index = list()
        values_index = list()

        values_list = ["z", "pvalue_raw"]

        condition_labels = data.columns.levels[0]
        for i, cnd1 in enumerate(condition_labels):
            for cnd2 in condition_labels[i + 1:]:
                cnd1_index.extend([cnd1] * len(values_list))
                cnd2_index.extend([cnd2] * len(values_list))
                values_index.extend(sorted(values_list))
        columns = pd.MultiIndex.from_tuples(zip(cnd1_index, cnd2_index, values_index), names=["condition1", "condition2", "value"])

        idx = pd.IndexSlice
        result_df = pd.DataFrame(np.nan, index=data.index, columns=columns)
        for i, cnd1 in enumerate(condition_labels):
            for cnd2 in condition_labels[i + 1:]:
                result_df.loc[:, idx[cnd1, cnd2, 'z']] = np.absolute(data.loc[:, idx[cnd1, 'score']] - data.loc[:, idx[cnd2, 'score']]) / np.sqrt(data.loc[:, idx[cnd1, 'SE']] ** 2 + data.loc[:, idx[cnd2, 'SE']] ** 2)
                result_df.loc[:, idx[cnd1, cnd2, 'pvalue_raw']] = 2 * stats.norm.sf(result_df.loc[:, idx[cnd1, cnd2, 'z']])

        self.store.put("/main/{}/scores_pvalues".format(label), result_df, format="table")


    def make_plots(self):
        if self.plots_requested:
            logging.info("Creating plots", extra={'oname' : self.name})

            # sequence-function maps
            if self.scoring_method != "counts":
                if "synonymous" in self.labels:
                    pdf = PdfPages(os.path.join(self.plot_dir, "sequence_function_map_aa.pdf"))
                    for condition in self.children: 
                        self.sfmap_wrapper(condition=condition.name, pdf=pdf, coding=True)
                    pdf.close()
                if "variants" in self.labels:
                    pdf = PdfPages(os.path.join(self.plot_dir, "sequence_function_map_nt.pdf"))
                    for condition in self.children: 
                        self.sfmap_wrapper(condition=condition.name, pdf=pdf, coding=False)
                    pdf.close()

        for s in self.selection_list():
            s.make_plots()


    def write_tsv(self):
        """
        Write each table from the store to its own tab-separated file.

        Files are written to a ``tsv`` directory in the default output location. 
        File names are the HDF5 key with ``'_'`` substituted for ``'/'``.
        """
        if self.tsv_requested:
            logging.info("Generating tab-separated output files", extra={'oname' : self.name})
            for k in self.store.keys():
                self.write_table_tsv(k)
        for s in self.selection_list():
            s.write_tsv()


    def is_coding(self):
        """
        Return ``True`` if the all :py:class:`~selection.Selection` in the 
        :py:class:`~experiment.Experiment` score protein-coding variants, else 
        ``False``.
        """
        return all(x.is_coding() for x in self.selection_list())


    def sfmap_wrapper(self, condition, pdf, coding):
        """
        Create a sequence function map for scores in *condition*.

        Uses :py:func:`~sfmap.sfmap_plot` for the plotting.
        """
        if coding:
            label = "amino acid"
        else:
            label = "nucleotide"

        logging.info("Creating sequence-function map ({}, {})".format(condition, label), extra={'oname' : self.name})

        idx = pd.IndexSlice

        if coding:
            df_name = '/main/synonymous/scores'
        else:
            df_name = '/main/variants/scores'
        data, wtseq = singleton_dataframe(self.store[df_name][idx[condition, 'score']], self.wt, coding=coding)
        data_se, _ = singleton_dataframe(self.store[df_name][idx[condition, 'SE']], self.wt, coding=coding)

        # format the title
        if coding:
            title = "Amino Acid"
        else:
            title = "Nucleotide"

        if self.scoring_method in ("WLS", "OLS"):
            title += " Sequence-Function Map\n{} ({} Slope)".format(condition, self.scoring_method)
        elif self.scoring_method == "ratios":
            title += " Sequence-Function Map\n{} ({})".format(condition, "Enrich2 Ratio")
        elif self.scoring_method == "simple":
            title += " Sequence-Function Map\n{} ({})".format(condition, "Simplified Ratio")
        else:
            raise ValueError("Invalid scoring method", self.name)

        sfmap_plot(df=data, pdf=pdf, style="scores", df_se=data_se, dimensions="tall", wt=wtseq, title=title)


    def correlation_plot(self, pdf, label):
        """
        Create a triangular heatmap showing the Pearson correlation coefficient
        for each pairwise comparison of replicate scores.
        """
        pass

