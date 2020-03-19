from __future__ import print_function
import os
import logging
import pandas as pd
import collections
import getpass
import time


#: Dictionary specifying available scoring methods for the analysis
#: Key is the internal name of the method, value is the GUI label
#: For command line options, internal name is used for the option string itself
#: and the value is the help string
SCORING_METHODS = collections.OrderedDict(
    [
        ("WLS", "weighted least squares"),
        ("ratios", "log ratios (Enrich2)"),
        ("counts", "counts only"),
        ("OLS", "ordinary least squares"),
        ("simple", "log ratios (old Enrich)"),
    ]
)


#: Dictionary specifying available scoring methods for the analysis
#: Key is the internal name of the method, value is the GUI label
#: For command line options, internal name is used for the option string itself
#: and the value is the help string
LOGR_METHODS = collections.OrderedDict(
    [
        ("wt", "wild type"),
        ("complete", "library size (complete cases)"),
        ("full", "library size (all reads)"),
    ]
)


#: List specifying valid labels in their sorted order
#: Sorted order is the order in which they should be calculated in
ELEMENT_LABELS = ["barcodes", "identifiers", "variants", "synonymous"]


def fix_filename(s):
    """
    Clean up a file name by removing invalid characters and converting 
    spaces to underscores.

    :param str s: file name
    :return: cleaned file name
    :rtype: str
    """
    fname = "".join(c for c in s if c.isalnum() or c in (" ._~"))
    fname = fname.replace(" ", "_")
    return fname


class StoreManager(object):
    """
    Abstract class for all data-containing classes
    (:py:class:`~enrich2.seqlib.seqlib.SeqLib`,
    :py:class:`~enrich2.selection.selection.Selection`, and
    :py:class:`~enrich2.experiment.experiment.Experiment`).

    Contains common operations for all data-containing classes, such as HDF5
    store and directory management.
    """

    store_suffix = None
    has_store = True
    treeview_class_name = None

    def __init__(self):
        self.logger = logging.getLogger("{}.{}".format(__name__, self.__class__))

        # general data members
        self._name = None
        self.name = "Unnamed" + self.__class__.__name__
        self.parent = None
        self._labels = list()

        # metadata members
        self.username = getpass.getuser()
        self.creationtime = time.asctime()

        # HDF5 data
        self.store_cfg = None
        self.store_path = None
        self.store = None
        self.chunksize = 100000

        # output locations
        self._output_dir = None
        self._output_dir_override = None
        self._plot_dir = None
        self._tsv_dir = None

        # analysis parameters
        self._scoring_method = None
        self._logr_method = None

        self._force_recalculate = None
        self._component_outliers = None
        self._plots_requested = None
        self._tsv_requested = None

        # GUI variables
        self.treeview_id = None
        self.treeview_info = None

        # Plotting options
        self.plot_options = None

    def child_labels(self):
        """
        Returns a list of labels shared by every child.
        """
        shared = list()
        for x in self.children:
            shared.extend(x.labels)
        shared = collections.Counter(shared)
        shared = [x for x in shared.keys() if shared[x] == len(self.children)]
        return sorted(shared, key=lambda a: ELEMENT_LABELS.index(a))

    @property
    def labels(self):
        """
        Property for returning labels shared by all children.

        If the object has no children (such as a SeqLib), the object's _labels
        are returned.
        """
        if self.children is None:
            return self._labels
        else:
            if len(self.children) > 0:
                return self.child_labels()
            else:
                return self._labels

    @property
    def force_recalculate(self):
        """
        This property should only be set for the root element. All other
        elements in the analysis should have ``None``.

        Recursively traverses up the config tree to find the root element.
        """
        if self._force_recalculate is None:
            if self.parent is not None:
                return self.parent.force_recalculate
            else:
                raise ValueError(
                    "Forced recalculation option not specified "
                    "at root [{}]".format(self.name)
                )
        else:
            return self._force_recalculate

    @force_recalculate.setter
    def force_recalculate(self, value):
        """
        Make sure the *value* is valid and set it.
        """
        if value in (True, False):
            self._force_recalculate = value
        else:
            raise ValueError(
                "Invalid setting '{}' for force_recalculate "
                "[{}]".format(value, self.name)
            )

    @property
    def component_outliers(self):
        """
        This property should only be set for the root element. All other
        elements in the analysis should have ``None``.

        Recursively traverses up the config tree to find the root element.
        """
        if self._component_outliers is None:
            if self.parent is not None:
                return self.parent.component_outliers
            else:
                raise ValueError(
                    "Calculate component outliers option not "
                    "specified at root [{}]".format(self.name)
                )
        else:
            return self._component_outliers

    @component_outliers.setter
    def component_outliers(self, value):
        """
        Make sure the *value* is valid and set it.
        """
        if value in (True, False):
            self._component_outliers = value
        else:
            raise ValueError(
                "Invalid setting '{}' for component_outliers "
                "[{}]".format(value, self.name)
            )

    @property
    def plots_requested(self):
        """
        This property should only be set for the root element. All other
        elements in the analysis should have ``None``.

        Recursively traverses up the config tree to find the root element.
        """
        if self._plots_requested is None:
            if self.parent is not None:
                return self.parent.plots_requested
            else:
                raise ValueError(
                    "Make plots option not specified at root " "[{}]".format(self.name)
                )
        else:
            return self._plots_requested

    @plots_requested.setter
    def plots_requested(self, value):
        """
        Make sure the *value* is valid and set it.
        """
        if value in (True, False):
            self._plots_requested = value
        else:
            raise ValueError(
                "Invalid setting '{}' for plots_requested "
                "[{}]".format(value, self.name)
            )

    @property
    def tsv_requested(self):
        """
        This property should only be set for the root element. All other
        elements in the analysis should have ``None``.

        Recursively traverses up the config tree to find the root element.
        """
        if self._tsv_requested is None:
            if self.parent is not None:
                return self.parent.tsv_requested
            else:
                raise ValueError(
                    "Write tsv option not specified at root " "[{}]".format(self.name)
                )
        else:
            return self._tsv_requested

    @tsv_requested.setter
    def tsv_requested(self, value):
        """
        Make sure the *value* is valid and set it.
        """
        if value in (True, False):
            self._tsv_requested = value
        else:
            raise ValueError(
                "Invalid setting '{}' for tsv_requested [{}]".format(value, self.name)
            )

    @property
    def scoring_method(self):
        """
        This property should only be set for the root element. All other
        elements in the analysis should have ``None``.

        Recursively traverses up the config tree to find the root element.
        """
        if self._scoring_method is None:
            if self.parent is not None:
                return self.parent.scoring_method
            else:
                raise ValueError(
                    "Scoring method not specified at root [{}]".format(self.name)
                )
        else:
            return self._scoring_method

    @scoring_method.setter
    def scoring_method(self, value):
        """
        Make sure the *value* is valid and set it.
        """
        if value in SCORING_METHODS.keys():
            self._scoring_method = value
        else:
            raise ValueError(
                "Invalid setting for scoring_method [{}]".format(self.name)
            )

    @property
    def output_dir(self):
        """
        Elements that have ``None`` as their output directory value (*i.e.*
        no value set in the confic object) will automatically use the output
        directory from the parent.
        """
        if self._output_dir is None:
            if self.parent is not None:
                return self.parent.output_dir
            else:
                raise ValueError(
                    "No output directory specified at top level "
                    "[{}]".format(self.name)
                )
        else:
            return self._output_dir

    @output_dir.setter
    def output_dir(self, dirname):
        """
        Sets the object's base output directory to *dirname* and creates the
        directory if it doesn't exist.
        """
        try:
            dirname = os.path.expanduser(dirname)  # handle leading '~'
        except AttributeError as e:
            raise AttributeError("Invalid input for output directory: {}".format(e))
        try:
            if not os.path.exists(dirname):
                os.makedirs(dirname)
        except OSError as e:
            raise OSError("Failed to create output directory: {}".format(e))
        self._output_dir = dirname

    @property
    def output_dir_override(self):
        """
        This property should only be set for the root element. All other
        elements in the analysis should have ``None``.

        Recursively traverses up the config tree to find the root element.
        """
        if self._output_dir_override is None:
            if self.parent is not None:
                return self.parent.output_dir_override
            else:
                raise ValueError(
                    "Output directory override not specified at root [{}]".format(
                        self.name
                    )
                )
        else:
            return self._output_dir_override

    @output_dir_override.setter
    def output_dir_override(self, value):
        """
        Make sure the *value* is valid and set it.
        """
        if value in (True, False):
            self._output_dir_override = value
        else:
            raise ValueError(
                "Invalid setting '{}' for output_dir_override [{}]".format(
                    value, self.name
                )
            )

    @property
    def plot_dir(self):
        """
        Creates the plot directory when first requested if it doesn't exist.

        The plot directory is ``<output directory>/plots/<object name>/``.
        """
        if self._plot_dir is None:
            dirname = os.path.join(
                self.output_dir,
                "plots",
                "{}_{}".format(fix_filename(self.name), self.store_suffix),
            )
            try:
                if not os.path.exists(dirname):
                    os.makedirs(dirname)
            except OSError as e:
                raise OSError("Failed to create plotting directory: {}".format(e))
            self._plot_dir = dirname
        return self._plot_dir

    @property
    def tsv_dir(self):
        """
        Creates the tsv directory when first requested if it doesn't exist.

        The plot directory is ``<output directory>/tsv/<object name>/``.
        """
        if self._tsv_dir is None:
            dirname = os.path.join(
                self.output_dir,
                "tsv",
                "{}_{}".format(fix_filename(self.name), self.store_suffix),
            )
            try:
                if not os.path.exists(dirname):
                    os.makedirs(dirname)
            except OSError as e:
                raise OSError("Failed to create tsv directory: {}".format(e))
            self._tsv_dir = dirname
        return self._tsv_dir

    @property
    def logr_method(self):
        """
        This property should only be set for the root element. All other
        elements in the analysis should have ``None``.

        Recursively traverses up the config tree to find the root element.
        """
        if self._logr_method is None:
            if self.parent is not None:
                return self.parent.logr_method
            else:
                raise AttributeError(
                    "Log-ratio method not specified for root "
                    "element [{}]".format(self.name)
                )
        else:
            return self._logr_method

    @logr_method.setter
    def logr_method(self, value):
        """
        Make sure the *value* is valid and set it.
        """
        if value in LOGR_METHODS.keys():
            self._logr_method = value
        else:
            raise ValueError(
                "Invalid setting '{}' for log-ratio method "
                "[{}]".format(value, self.name)
            )

    @property
    def children(self):
        """
        Property for returning all children of this element.

        Calls the objects :py:meth:`_children` method. This allows the
        property to be overloaded without creating a new property in the
        derived classes.
        """
        return self._children()

    def _children(self):
        """
        Pure virtual method for returning children.
        """
        raise NotImplementedError("must be implemented by subclass")

    def remove_child_id(self, tree_id):
        """
        Pure virtual method for removing a child with Treeview id *tree_id*
        from this element's children.
        """
        raise NotImplementedError("must be implemented by subclass")

    def add_child(self, child):
        """
        Pure virtual method for adding a child.
        """
        raise NotImplementedError("must be implemented by subclass")

    def child_names(self):
        """
        Return a list of the ``name`` attributes for all children.
        """
        try:
            names = [x.name for x in self.children]
        except AttributeError:
            raise AttributeError("No name set for child [{}]".format(self.name))
        else:
            return names

    def add_label(self, x):
        """
        Add element label to this object.
        """
        labels = set(self._labels)
        if isinstance(x, str):
            if x in ELEMENT_LABELS:
                labels.update([x])
            else:
                raise ValueError("Invalid element label '{}' [{}]".format(x, self.name))
        else:
            raise AttributeError("Failed to add labels [{}]".format(self.name))
        # sort based on specified order
        self._labels = sorted(labels, key=lambda a: ELEMENT_LABELS.index(a))

    def configure(self, cfg):
        """
        Set up the object using the config object *cfg*, usually derived from
        a ``.json`` file.
        """
        try:
            self.name = cfg["name"]
            if "output directory" in cfg:
                if self.output_dir_override:
                    self.logger.warning(
                        "Using command line supplied output "
                        "directory instead of config file output "
                        "directory"
                    )
                else:
                    self.output_dir = cfg["output directory"]
            if "store" in cfg:
                self.store_cfg = True
                if not os.path.exists(cfg["store"]):
                    raise IOError(
                        'Specified store file "{}" not found'.format(cfg["store"]),
                        self.name,
                    )
                elif os.path.splitext(cfg["store"])[-1].lower() != ".h5":
                    raise ValueError(
                        "Unrecognized store file extension for "
                        '"{}"'.format(cfg["store"]),
                        self.name,
                    )
                else:
                    self.store_path = cfg["store"]
                    self.logger.info(
                        'Using specified HDF5 data store "{}"'.format(self.store_path)
                    )
            else:
                self.store_cfg = False
                self.store_path = None
        except KeyError as key:
            raise KeyError("Missing required config value {}".format(key), self.name)

    def serialize(self):
        """
        Format this object (and its children) as a config object suitable for
        dumping to a config file.
        """
        cfg = {"name": self.name}
        if self.store_cfg:
            cfg["store"] = self.store_path
        if self.output_dir is not None:
            if self.parent is not None:
                if self.output_dir != self.parent.output_dir:
                    cfg["output directory"] = self.output_dir
            else:
                cfg["output directory"] = self.output_dir
        return cfg

    def validate(self):
        """
        Pure virtual method for making sure configured object is valid.
        """
        raise NotImplementedError("must be implemented by subclass")

    def store_open(self, children=False, force_delete=True):
        """
        Open the HDF5 file associated with this object. If the
        ``force_recalculate`` option is selected and ``force_delete`` is
        ``True``, the existing tables under ``'/main'`` will be deleted upon
        opening.

        This method needs a lot more error checking.
        """
        if self.has_store:
            if not self.store_cfg:
                self.store_path = os.path.join(
                    self.output_dir,
                    "{}_{}.h5".format(fix_filename(self.name), self.store_suffix),
                )
                if os.path.exists(self.store_path):
                    self.logger.info(
                        'Found existing HDF5 data store "{}"'.format(self.store_path)
                    )
                else:
                    self.logger.info(
                        'Creating new HDF5 data store "{}"'.format(self.store_path)
                    )
            self.store = pd.HDFStore(self.store_path)
            if self.force_recalculate and force_delete:
                if "/main" in self.store:
                    self.logger.info("Deleting existing calculated values")
                    self.store.remove("/main")
                else:
                    self.logger.warning("No existing calculated values in file")

        if children and self.children is not None:
            for child in self.children:
                child.store_open(children=True)

    def store_close(self, children=False):
        # needs more error checking
        if self.has_store:
            self.store.close()
            self.store = None
        if children and self.children is not None:
            for child in self.children:
                child.store_close(children=True)

    def get_metadata(self, key, store=None):
        """
        Retrieve the Enrich2 metadata dictionary from the HDF5 store.

        *key* is the name of the group or node in the HDF5 data store.

        Returns the metadata dictionary for *key*. If no metadata has been set
        for *key*, returns ``None``.

        *store* can be an external open HDFStore (used when copying metadata
        from raw counts). If it is ``None``, use this object's store.

        """
        if store is None:
            store = self.store
        try:
            metadata = store.get_storer(key).attrs["enrich2"]
        except AttributeError:
            if store is self.store:  # store parameter was None
                raise AttributeError(
                    "Invalid HDF store node '{}' [{}]".format(key, self.name)
                )
            else:
                raise AttributeError(
                    "Invalid external HDF store node '{}' in "
                    "'{}' [{}]".format(key, store.filename, self.name)
                )
        except KeyError:  # no enrich2 metadata
            return None
        else:
            return metadata

    def set_metadata(self, key, d, update=True):
        """
        Replace or update the metadata dictionary from the HDF5 store.

        *key* is the name of the group or node in the HDF5 data store.

        *d* is the dictionary containing the new metadata.

        If *update* is ``False``, *d* completely replaces the existing metadata
        for *key*. Otherwise, *d* updates the existing metadata using standard
        dictionary update.
        """
        # get existing for update
        # also performs check for the node, so we don't check here
        metadata = self.get_metadata(key)
        if metadata is None or not update:
            metadata = d
        else:
            metadata.update(d)
        self.store.get_storer(key).attrs["enrich2"] = metadata

    def calculate(self):
        """
        Pure virtual method that defines how the data are calculated.
        """
        raise NotImplementedError("must be implemented by subclass")

    def make_plots(self):
        """
        Pure virtual method for creating plots.
        """
        raise NotImplementedError("must be implemented by subclass")

    def write_tsv(self):
        """
        Pure virtual method for writing tsv files.
        """
        raise NotImplementedError("must be implemented by subclass")

    def write_table_tsv(self, key):
        """
        Write the table under *key* as a tsv file.

        Files are written to a ``tsv`` directory in the default output
        location.

        File names are the HDF5 key with ``'_'`` substituted for ``'/'``.
        """
        fname = key.strip("/")  # remove leading slash
        fname = fname.replace("/", "_") + ".tsv"
        self.store[key].to_csv(os.path.join(self.tsv_dir, fname), sep="\t", na_rep="NA")

    def check_store(self, key):
        """
        Checks to see if a particular data frame in the HDF5 store already
        exists.

        Args:
            key (str): key for the requested data frame

        Returns:
            bool: True if the key exists in the HDF5 store, else False.
        """
        if key in self.store.keys():
            self.logger.info("Found existing '{}'".format(key))
            return True
        else:
            return False

    def map_table(
        self,
        source,
        destination,
        source_query=None,
        row_callback=None,
        row_callback_args=None,
        destination_data_columns=None,
    ):
        """
        Converts source table into destination table.

        This method really needs a better name.
        """
        if destination in self.store.keys():
            # remove the current destination table because we are using append
            # append takes the "min_itemsize" argument, and put doesn't
            self.logger.info("Overwriting existing '{}'".format(destination))
            self.store.remove(destination)

        # turn the single table name into a list to use select_as_multiple
        if isinstance(source, str):
            source = [source]

        # assumes the source tables all have the same index
        # find the min_itemsize
        max_index_length = self.store.select_column(source[0], "index").map(len).max()
        for df in self.store.select_as_multiple(
            source, source_query, chunksize=self.chunksize, selector=source[0]
        ):
            if row_callback is not None:
                df = df.apply(row_callback, args=row_callback_args, axis="columns")
            if destination not in self.store:
                if destination_data_columns is None:
                    # if not specified, index all columns
                    destination_data_columns = list(df.columns)
                self.store.append(
                    destination,
                    df,
                    min_itemsize={"index": max_index_length},
                    data_columns=destination_data_columns,
                )
            else:
                self.store.append(destination, df)

    def combined_index(self, tables):
        """
        Return an index containing all elements in *tables*
        """
        shared = pd.Index()
        for t in tables:
            shared = shared.union(pd.Index(self.store.select_column(t, "index")))
        return shared

    def get_root(self):
        if self.parent is not None:
            return self.parent.get_root()
        else:
            return self
