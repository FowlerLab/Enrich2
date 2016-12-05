#  Copyright 2016 Alan F Rubin
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
import time
import logging
import os.path
import pandas as pd
import numpy as np
from collections import OrderedDict
from matplotlib.backends.backend_pdf import PdfPages
import sys
from .plots import counts_plot
from .storemanager import StoreManager, fix_filename

class SeqLib(StoreManager):
    """
    Abstract class for handling count data from a single sequencing library.
    
    """
    
    # Note: the following block is referenced by line number above
    # When adding new messages, update the documentation line numbers also!
    filter_messages = OrderedDict([
                ('min quality', "single-base quality"),
                ('avg quality', "average quality"),
                ('max N', "excess N bases"),
                ('chastity', "not chaste"),
                ('max mutations', "excess mutations"),
                ('remove unresolvable', "unresolvable mismatch"),
                ('merge failure', "unable to merge reads"),
                ('total', "total")
            ])

    store_suffix = "lib"
    
    def __init__(self):
        StoreManager.__init__(self)
        self.timepoint = None
        self.counts_file = None
        self.report_filtered = None
        self._filters = dict()
        self.filter_stats = dict()
        self.default_filters = dict()
        self.default_filters.update({'min quality' : 0})
        self.default_filters.update({'max N' : sys.maxsize})
        self.default_filters.update({'avg quality' : 0})
        self.default_filters.update({'chastity' : False})


    @property
    def filters(self):
        return self._filters

    @filters.setter
    def filters(self, config_filters):
        """
        Set up the filters dictionary using the options selected in *config_filters*, filling in missing entries with defaults.
        """
        self._filters.clear()
        self._filters.update(self.default_filters)

        unused = list()
        for key in config_filters:
            if key in self._filters:
                if config_filters[key] is not None:
                    self._filters[key] = config_filters[key]
            else:
                unused.append(key)
        if len(unused) > 0:
            logging.warning("Unused filter parameters ({})".format(', '.join(unused)), extra={'oname' : self.name})

        self.filter_stats.clear()
        for key in self._filters:
            self.filter_stats[key] = 0
        self.filter_stats['total'] = 0


    def serialize_filters(self):
        """
        Return a dictionary of filtering options that have non-default values.
        """
        cfg = dict()
        for key in self.filters.keys():
            if self.filters[key] != self.default_filters[key]:
                cfg[key] = self.filters[key]
        return cfg


    def _children(self):
        """
        These objects have no children. Returns ``None``.
        """
        return None


    def add_child(self, child):
        """
        No children, raises an AttributeError.
        """
        raise AttributeError("SeqLib objects do not support adding children")


    def remove_child_id(self, tree_id):
        """
        No children, raises an AttributeError.
        """
        raise AttributeError("SeqLib objects do not support removing children")


    def validate(self):
        """
        Validates paramaters for a configured SeqLib. Currently does nothing.
        """
        pass


    def has_wt(self):
        """
        Returns whether or not the object has a wild type sequence. Returns ``False``
        unless overloaded by a derived class (such as :py:class:`~seqlib.seqlib.VariantSeqLib`).
        """
        return False


    def configure(self, cfg):
        """
        Set up the object using the config object *cfg*, usually derived from 
        a ``.json`` file.
        """
        StoreManager.configure(self, cfg)
        try:
            self.timepoint = int(cfg['timepoint'])
            if 'report filtered reads' in cfg:
                self.report_filtered = cfg['report filtered reads']
            else:
                self.report_filtered = False
            if 'counts file' in cfg:
                self.counts_file = cfg['counts file']
            else:
                self.counts_file = None
        except KeyError as key:
            raise KeyError("Missing required config value {key}".format(key=key), 
                              self.name)
        except ValueError as value:
            raise ValueError("Invalid parameter value {value}".format(value=value), self.name)


    def serialize(self):
        """
        Format this object (and its children) as a config object suitable for dumping to a config file.
        """
        cfg = StoreManager.serialize(self)
        cfg['timepoint'] = self.timepoint
        cfg['report filtered reads'] = self.report_filtered
        if self.counts_file is not None:
            cfg['counts file'] = self.counts_file

        return cfg


    def calculate(self):
        """
        Pure virtual method that defines how the data are counted.
        """
        raise NotImplementedError("must be implemented by subclass")


    def report_filtered_read(self, fq, filter_flags):
        """
        Write the :py:class:`~fqread.FQRead` object *fq* to the ``DEBUG``
        logging . The dictionary *filter_flags* contains ``True`` 
        values for each filtering option that applies to *fq*. Keys in 
        *filter_flags* are converted to messages using the 
        ``StoreManager.filter_messages`` dictionary.
        """
        logging.debug("Filtered read ({messages})\n{read!s}".format(
                      messages=', '.join(StoreManager.filter_messages[x] 
                                for x in filter_flags if filter_flags[x]), 
                      name=self.name, read=fq), extra={'oname' : self.name})


    def save_counts(self, label, df_dict, raw):
        """
        Convert the counts in the dictionary *df_dict* into a DataFrame object and save it to the data store.

        If *raw* is ``True``, the counts are stored under ``"/raw/label/counts"``; else ``"/main/label/counts"``.
        """
        if len(df_dict.keys()) == 0:
            raise ValueError("Failed to count {} [{}]".format(label, self.name))
        df = pd.DataFrame.from_dict(df_dict, orient="index", dtype="int32")
        df.columns = ['count']
        df.sort_values('count', ascending=False, inplace=True)
        logging.info("Counted {n} {label} ({u} unique)".format(
                n=df['count'].sum(), u=len(df.index), label=label), extra={'oname' : self.name})
        if raw:
            key = "/raw/{}/counts".format(label)
        else:
            key = "/main/{}/counts".format(label)
        self.store.put(key, df, format="table", data_columns=df.columns)
        del df


    def save_filtered_counts(self, label, query):
        """
        Filter the counts in ``"/raw/label/counts"`` using the *query* string and store the result in ``"/main/label/counts"``

        For more information on building query strings, see http://pandas.pydata.org/pandas-docs/stable/io.html#querying-a-table
        """
        logging.info("Converting raw {} counts to main counts".format(label), extra={'oname' : self.name})
        raw_table = "/raw/{}/counts".format(label)
        main_table = "/main/{}/counts".format(label)
        self.map_table(source=raw_table, destination=main_table, source_query=query)
        logging.info("Counted {n} {label} ({u} unique) after query".format(
                n=self.store[main_table]['count'].sum(),
                u=len(self.store[main_table].index),
                label=label), extra={'oname' : self.name})


    def report_filter_stats(self):
        """
        Create report file for the number of filtered reads.

        The report file is located in the output directory, named ``SeqLibName.filter.txt``. 
        It contains the number of reads filtered for each category, plus the total number filtered.

        .. note:: Reads are checked for all quality-based criteria before filtering.
        """
        with open(os.path.join(self.output_dir, fix_filename(self.name) + ".filter.txt"), "w") as handle:
            elements = list()
            for key in sorted(self.filter_stats, key=self.filter_stats.__getitem__, reverse=True):
                if key != 'total' and self.filter_stats[key] > 0:
                    print(SeqLib.filter_messages[key], self.filter_stats[key], sep="\t", file=handle)
            print("total", self.filter_stats['total'], sep="\t", file=handle)
        logging.info("Wrote filtering statistics", extra={'oname' : self.name})


    def save_filter_stats(self):
        """
        Save a DataFrame containing the number of filtered reads under ``'/raw/filter'``.

        This DataFrame contains the same information as ``report_filter_stats``
        """
        df = pd.DataFrame(index=SeqLib.filter_messages.values(), columns=['count'])
        for key in self.filter_stats.keys():
            if self.filter_stats[key] > 0 or key == 'total':
                df.loc[SeqLib.filter_messages[key], 'count'] = self.filter_stats[key]
        df.dropna(inplace=True)
        self.store.put('/raw/filter', df.astype(int), format="table", data_columns=df.columns)


    def read_quality_filter(self, fq):
        """
        Check the quality of the FQRead object *fq*.

        Checks ``'chastity'``, ``'min quality'``, ``'avg quality'``, ``'max N'``, and ``'remove unresolvable'``.
        Counts failed reads for later output and reports the filtered read if desired. 
        Returns ``True`` if the read passes all filters, else ``False``.
        """
        filter_flags = dict()
        for key in self.filters:
            filter_flags[key] = False

        if self.filters['chastity']:
            if not fq.is_chaste():
                self.filter_stats['chastity'] += 1
                filter_flags['chastity'] = True

        if self.filters['min quality'] > 0:
            if fq.min_quality() < self.filters['min quality']:
                self.filter_stats['min quality'] += 1
                filter_flags['min quality'] = True

        if self.filters['avg quality'] > 0:
            if fq.mean_quality() < self.filters['avg quality']:
                self.filter_stats['avg quality'] += 1
                filter_flags['avg quality'] = True

        if self.filters['max N'] > 0:
            if fq.sequence.upper().count('N') > self.filters['max N']:
                self.filter_stats['max N'] += 1
                filter_flags['max N'] = True

        if 'remove unresolvable' in self.filters:   # OverlapSeqLib only
            if self.filters['remove unresolvable']:
                if 'X' in fq.sequence:
                    self.filter_stats['remove unresolvable'] += 1
                    filter_flags['remove unresolvable'] = True

        # update total and report if failed 
        if any(filter_flags.values()):
            self.filter_stats['total'] += 1
            if self.report_filtered:
                self.report_filtered_read(fq, filter_flags)
            return False
        else:
            return True


    def make_plots(self):
        """
        Make plots that are shared by all :py:class:`~seqlib.seqlib.SeqLib` objects.

        Creates counts histograms for all labels.
        """
        if self.plots_requested:
            logging.info("Creating plots", extra={'oname' : self.name})

            pdf = PdfPages(os.path.join(self.plot_dir, "counts.pdf"))
            for label in self.labels:
                counts_plot(self, label, pdf, log=True)
                counts_plot(self, label, pdf, log=False)
            pdf.close()


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


    def copy_raw(self):
        """
        If an HDF store containing raw counts has been specified, open the store,
        copy those counts into this store, and close the counts store.

        Copies all tables in the ``'/raw'`` group along with their metadata.
        """
        if not os.path.exists(self.counts_file):
            raise IOError("Counts file '{}' not found [{}]".format(self.counts_file, self.name))
        elif os.path.splitext(self.counts_file)[-1].lower() not in  (".h5", ".hdf5", "hdf"):
            raise ValueError("Unrecognized counts file extension for '{}' [{}]".format(self.counts_file, self.name))
        else:
            store = pd.HDFStore(self.counts_file)
            logging.info("Using existing HDF5 data store '{}' for raw data".format(self.counts_file), extra={'oname' : self.name})
            # this could probably be much more efficient, but the PyTables docs don't explain copying subsets of files adequately
            raw_keys = [key for key in store.keys() if key.startswith("/raw/")]
            if len(raw_keys) == 0:
                raise ValueError("No raw counts found in '{}' [{}]".format(self.counts_file, self.name))
            else:
                for k in raw_keys:
                    # copy the data table
                    raw = store[k]
                    self.store.put(k, raw, format="table", data_columns=raw.columns)
                    # copy the metadata
                    self.set_metadata(k, self.get_metadata(k, store=store), update=False)
                    logging.info("Copied raw data '{}'".format(k), extra={'oname' : self.name})
            store.close()



