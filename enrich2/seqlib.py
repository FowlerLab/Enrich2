from __future__ import print_function
import logging
import os.path
import pandas as pd
import numpy as np
from collections import OrderedDict
from matplotlib.backends.backend_pdf import PdfPages
import sys
from .plots import counts_plot
from .storemanager import StoreManager, fix_filename, ELEMENT_LABELS


class SeqLib(StoreManager):
    """
    Abstract class for handling count data from a single sequencing library.
    """

    # Note: the following block is referenced by line number above
    # When adding new messages, update the documentation line numbers also!
    filter_messages = OrderedDict(
        [
            ("min quality", "single-base quality"),
            ("avg quality", "average quality"),
            ("max N", "excess N bases"),
            ("chastity", "not chaste"),
            ("remove unresolvable", "unresolvable mismatch"),
            ("merge failure", "unable to merge reads"),
            ("total", "total"),
        ]
    )

    store_suffix = "lib"

    def __init__(self):
        StoreManager.__init__(self)
        self.logger = logging.getLogger("{}.{}".format(__name__, self.__class__))
        self.timepoint = None
        self.counts_file = None
        self.report_filtered = None
        self._filters = dict()
        self.filter_stats = dict()
        self.default_filters = dict()
        self.default_filters.update({"min quality": 0})
        self.default_filters.update({"max N": sys.maxsize})
        self.default_filters.update({"avg quality": 0})
        self.default_filters.update({"chastity": False})

    @property
    def filters(self):
        return self._filters

    @filters.setter
    def filters(self, config_filters):
        """
        Set up the filters dictionary using the options selected in
        *config_filters*, filling in missing entries with defaults.
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
            self.logger.warning(
                "Unused filter parameters ({})" "".format(", ".join(unused))
            )

        self.filter_stats.clear()
        for key in self._filters:
            self.filter_stats[key] = 0
        self.filter_stats["total"] = 0

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

    def has_wt_sequence(self):
        """
        Returns whether or not the object has a wild type sequence. Returns
        ``False`` unless overloaded by a derived class (such as
        :py:class:`~seqlib.seqlib.VariantSeqLib`).
        """
        return False

    def configure(self, cfg):
        """
        Set up the object using the config object *cfg*, usually derived from
        a ``.json`` file.
        """
        StoreManager.configure(self, cfg)
        self.logger = logging.getLogger(
            "{}.{} - {}".format(__name__, self.__class__.__name__, self.name)
        )
        try:
            self.timepoint = int(cfg["timepoint"])
            if "report filtered reads" in cfg:
                self.report_filtered = cfg["report filtered reads"]
            else:
                self.report_filtered = False
            if "counts file" in cfg:
                self.counts_file = cfg["counts file"]
            else:
                self.counts_file = None
        except KeyError as key:
            raise KeyError(
                "Missing required config value {key}" "".format(key=key), self.name
            )
        except ValueError as value:
            raise ValueError(
                "Invalid parameter value {value}" "".format(value=value), self.name
            )

    def serialize(self):
        """
        Format this object (and its children) as a config object suitable for
        dumping to a config file.
        """
        cfg = StoreManager.serialize(self)
        cfg["timepoint"] = self.timepoint
        cfg["report filtered reads"] = self.report_filtered
        if self.counts_file is not None:
            cfg["counts file"] = self.counts_file

        return cfg

    def calculate(self):
        """
        Pure virtual method that defines how the data are counted.
        """
        raise NotImplementedError("must be implemented by subclass")

    def report_filtered_read(self, fq, filter_flags):
        """
        Write the :py:class:`~fqread.FQRead` object *fq* to the ``DEBUG``
        log. The dictionary *filter_flags* contains ``True``
        values for each filtering option that applies to *fq*. Keys in
        *filter_flags* are converted to messages using the
        ``SeqLib.filter_messages`` dictionary.
        """
        self.logger.debug(
            "Filtered read ({messages})\n{read!s}".format(
                messages=", ".join(
                    SeqLib.filter_messages[x] for x in filter_flags if filter_flags[x]
                ),
                name=self.name,
                read=fq,
            )
        )

    def save_counts(self, label, df_dict, raw):
        """
        Convert the counts in the dictionary *df_dict* into a DataFrame object
        and save it to the data store.

        If *raw* is ``True``, the counts are stored under
        ``"/raw/label/counts"``; else ``"/main/label/counts"``.
        """
        if len(df_dict.keys()) == 0:
            raise ValueError("Failed to count {} [{}]".format(label, self.name))
        df = pd.DataFrame.from_dict(df_dict, orient="index", dtype=np.int32)
        df.columns = ["count"]
        df.sort_values("count", ascending=False, inplace=True)
        self.logger.info(
            "Counted {n} {label} ({u} unique)".format(
                n=df["count"].sum(), u=len(df.index), label=label
            )
        )
        if raw:
            key = "/raw/{}/counts".format(label)
        else:
            key = "/main/{}/counts".format(label)
        self.store.put(key, df, format="table", data_columns=df.columns)
        del df

    def save_filtered_counts(self, label, query):
        """
        Filter the counts in ``"/raw/label/counts"`` using the *query* string
        and store the result in ``"/main/label/counts"``

        For more information on building query strings, see
        http://pandas.pydata.org/pandas-docs/stable/io.html#querying-a-table
        """
        self.logger.info("Converting raw {} counts to main counts".format(label))
        raw_table = "/raw/{}/counts".format(label)
        main_table = "/main/{}/counts".format(label)
        self.map_table(source=raw_table, destination=main_table, source_query=query)
        self.logger.info(
            "Counted {n} {label} ({u} unique) after query".format(
                n=self.store[main_table]["count"].sum(),
                u=len(self.store[main_table].index),
                label=label,
            )
        )

    def report_filter_stats(self):
        """
        Create report file for the number of filtered reads.

        The report file is located in the output directory, named
        ``SeqLibName.filter.txt``.
        It contains the number of reads filtered for each category, plus the
        total number filtered.

        .. note:: Reads are checked for all quality-based criteria before \
        filtering.
        """
        with open(
            os.path.join(self.output_dir, fix_filename(self.name) + ".filter.txt"), "w"
        ) as handle:
            for key in sorted(
                self.filter_stats, key=self.filter_stats.__getitem__, reverse=True
            ):
                if key != "total" and self.filter_stats[key] > 0:
                    print(
                        SeqLib.filter_messages[key],
                        self.filter_stats[key],
                        sep="\t",
                        file=handle,
                    )
            print("total", self.filter_stats["total"], sep="\t", file=handle)
        self.logger.info("Wrote filtering statistics")

    def save_filter_stats(self):
        """
        Save a DataFrame containing the number of filtered reads under
        ``'/raw/filter'``.

        This DataFrame contains the same information as ``report_filter_stats``
        """
        df = pd.DataFrame(index=SeqLib.filter_messages.values(), columns=["count"])
        for key in self.filter_stats.keys():
            if self.filter_stats[key] > 0 or key == "total":
                df.loc[SeqLib.filter_messages[key], "count"] = self.filter_stats[key]
        df.dropna(inplace=True)
        self.store.put(
            "/raw/filter", df.astype(int), format="table", data_columns=df.columns
        )

    def read_quality_filter(self, fq):
        """
        Check the quality of the FQRead object *fq*.

        Checks ``'chastity'``, ``'min quality'``, ``'avg quality'``,
        ``'max N'``, and ``'remove unresolvable'``.
        Counts failed reads for later output and reports the filtered read if
        desired. 
        Returns ``True`` if the read passes all filters, else ``False``.
        """
        filter_flags = dict()
        for key in self.filters:
            filter_flags[key] = False

        if self.filters["chastity"]:
            if not fq.is_chaste():
                self.filter_stats["chastity"] += 1
                filter_flags["chastity"] = True

        if self.filters["min quality"] > 0:
            if fq.min_quality() < self.filters["min quality"]:
                self.filter_stats["min quality"] += 1
                filter_flags["min quality"] = True

        if self.filters["avg quality"] > 0:
            if fq.mean_quality() < self.filters["avg quality"]:
                self.filter_stats["avg quality"] += 1
                filter_flags["avg quality"] = True

        if self.filters["max N"] >= 0:
            if fq.sequence.upper().count("N") > self.filters["max N"]:
                self.filter_stats["max N"] += 1
                filter_flags["max N"] = True

        if "remove unresolvable" in self.filters:  # OverlapSeqLib only
            if self.filters["remove unresolvable"]:
                if "X" in fq.sequence:
                    self.filter_stats["remove unresolvable"] += 1
                    filter_flags["remove unresolvable"] = True

        # update total and report if failed
        if any(filter_flags.values()):
            self.filter_stats["total"] += 1
            if self.report_filtered:
                self.report_filtered_read(fq, filter_flags)
            return False
        else:
            return True

    def make_plots(self):
        """
        Make plots that are shared by all :py:class:`~seqlib.seqlib.SeqLib`
        objects.

        Creates counts histograms for all labels.
        """
        if self.plots_requested:
            self.logger.info("Creating plots")

            pdf = PdfPages(os.path.join(self.plot_dir, "counts.pdf"))
            for label in self.labels:
                counts_plot(self, label, pdf, log=True)
                counts_plot(self, label, pdf, log=False)
            pdf.close()

    def write_tsv(self):
        """
        Write each table from the store to its own tab-separated file.

        Files are written to a ``tsv`` directory in the default output
        location.
        File names are the HDF5 key with ``'_'`` substituted for ``'/'``.
        """
        if self.tsv_requested:
            self.logger.info("Generating tab-separated output files")
            for k in self.store.keys():
                self.write_table_tsv(k)

    def counts_from_file_h5(self, fname):
        """
        If an HDF store containing raw counts has been specified, open the
        store, copy those counts into this store, and close the counts store.

        Copies all tables in the ``'/raw'`` group along with their metadata.
        """
        store = pd.HDFStore(fname)
        self.logger.info(
            "Using existing HDF5 data store '{}' for raw data" "".format(fname)
        )
        # this could probably be much more efficient, but the PyTables docs
        # don't explain copying subsets of files adequately
        raw_keys = [key for key in store.keys() if key.startswith("/raw/")]
        if len(raw_keys) == 0:
            raise ValueError(
                "No raw counts found in '{}' [{}]" "".format(fname, self.name)
            )
        else:
            for k in raw_keys:
                # copy the data table
                raw = store[k]
                self.store.put(k, raw, format="table", data_columns=raw.columns)
                # copy the metadata
                self.set_metadata(k, self.get_metadata(k, store=store), update=False)
                self.logger.info("Copied raw data '{}'".format(k))
        store.close()

    def counts_from_file_tsv(self, fname):
        """
        If a counts file in tsv format has been specified, read the counts into
        a new dataframe and save as raw counts.
        """
        df = pd.read_table(fname, sep="\t", header=0, index_col=0)
        if df.columns != ["count"]:
            raise ValueError(
                "Invalid column names for counts file [{}]" "".format(self.name)
            )
        if len(df) == 0:
            raise ValueError("Empty counts file [{}]".format(self.name))
        label = None
        for elem in ELEMENT_LABELS:
            if elem in self.labels:
                label = elem
                break
        if label is None:
            raise ValueError("No valid element labels [{}]".format(self.name))
        key = "/raw/{}/counts".format(label)
        self.store.put(key, df, format="table", data_columns=df.columns, dtype=np.int32)

    def counts_from_file(self, fname):
        """Get raw counts from a counts file instead of FASTQ_ file.

        The ``'/raw/<element>/counts'`` table will be populated using the given
        input file. The input file should be a two-column file readable by
        ``pandas`` as a series or two-column dataframe or an Enrich2 HDF5 file.

        If the input file is a two-column file, the index will be checked using
        the SeqLib's ``validate_index()`` method.

        If the input file is an HDF5 file, the entire set of ``'/raw'`` tables
        will be copied over, with the metadata intact.
        """
        if not os.path.exists(fname):
            raise IOError("Counts file '{}' not found [{}]" "".format(fname, self.name))
        elif os.path.splitext(fname)[-1].lower() in (".h5"):
            self.counts_from_file_h5(self.counts_file)
        elif os.path.splitext(fname)[-1].lower() in (".txt", ".tsv", ".csv"):
            self.counts_from_file_tsv(self.counts_file)
        else:
            raise ValueError(
                "Unrecognized counts file extension for '{}' "
                "[{}]".format(fname, self.name)
            )
