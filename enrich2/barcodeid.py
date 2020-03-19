import logging
from .seqlib import SeqLib
from .barcode import BarcodeSeqLib
from .barcodemap import BarcodeMap
import pandas as pd
from .plots import barcodemap_plot
from matplotlib.backends.backend_pdf import PdfPages
import os.path


class BcidSeqLib(BarcodeSeqLib):
    """
    Class for counting data from barcoded sequencing libraries with non-variant
    identifiers. 
    Creating a :py:class:`BcidSeqLib` requires a valid *config* 
    object with an ``'barcodes'`` entry and information.

    The ``barcode_map`` keyword argument can be used to pass an existing 
    :py:class:`~seqlib.barcodemap.BarcodeMap`. Ensuring this is the 
    right :py:class:`~seqlib.barcodemap.BarcodeMap` is the responsibility 
    of the caller.
    """

    treeview_class_name = "Barcoded ID SeqLib"

    def __init__(self):
        BarcodeSeqLib.__init__(self)
        self.barcode_map = None
        self.identifier_min_count = None
        self.add_label("identifiers")
        self.logger = logging.getLogger("{}.{}".format(__name__, self.__class__))

    def configure(self, cfg, barcode_map=None):
        """
        Set up the object using the config object *cfg*, usually derived from 
        a ``.json`` file.
        """
        BarcodeSeqLib.configure(self, cfg)
        self.logger = logging.getLogger(
            "{}.{} - {}".format(__name__, self.__class__.__name__, self.name)
        )
        try:
            if "min count" in cfg["identifiers"]:
                self.identifier_min_count = int(cfg["identifiers"]["min count"])
            else:
                self.identifier_min_count = 0

            if barcode_map is not None:
                if barcode_map.filename == cfg["barcodes"]["map file"]:
                    self.barcode_map = barcode_map
                else:
                    raise ValueError(
                        "Attempted to assign non-matching barcode map [{}]".format(
                            self.name
                        )
                    )
            else:
                self.barcode_map = BarcodeMap(
                    cfg["barcodes"]["map file"], is_variant=False
                )
        except KeyError as key:
            raise KeyError(
                "Missing required config value {key} [{name}]".format(
                    key=key, name=self.name
                )
            )

    def serialize(self):
        """
        Format this object (and its children) as a config object suitable for dumping to a config file.
        """
        cfg = BarcodeSeqLib.serialize(self)

        cfg["identifiers"] = dict()
        if self.identifier_min_count > 0:
            cfg["identifiers"]["min count"] = self.identifier_min_count

        if self.barcode_map is not None:  # required for creating new objects in GUI
            cfg["barcodes"]["map file"] = self.barcode_map.filename

        return cfg

    def calculate(self):
        """
        Counts the barcodes using :py:meth:`BarcodeSeqLib.count` and combines them into 
        identifier counts using the :py:class:`BarcodeMap`.
        """
        if not self.check_store("/main/identifiers/counts"):
            BarcodeSeqLib.calculate(self)  # count the barcodes
            df_dict = dict()
            barcode_identifiers = dict()

            self.logger.info("Converting barcodes to identifiers")
            # store mapped barcodes
            self.save_filtered_counts(
                "barcodes",
                "index in self.barcode_map.keys() & count >= self.barcode_min_count",
            )

            # count identifiers associated with the barcodes
            for bc, count in self.store["/main/barcodes/counts"].iterrows():
                count = count["count"]
                identifier = self.barcode_map[bc]
                try:
                    df_dict[identifier] += count
                except KeyError:
                    df_dict[identifier] = count
                barcode_identifiers[bc] = identifier

            # save counts, filtering based on the min count
            self.save_counts(
                "identifiers",
                {
                    k: v
                    for k, v in df_dict.iteritems()
                    if v >= self.identifier_min_count
                },
                raw=False,
            )
            del df_dict

            # write the active subset of the BarcodeMap to the store
            barcodes = barcode_identifiers.keys()
            barcode_identifiers = pd.DataFrame(
                {"value": [barcode_identifiers[bc] for bc in barcodes]}, index=barcodes
            )
            del barcodes
            barcode_identifiers.sort_values("value", inplace=True)
            self.store.put(
                "/raw/barcodemap",
                barcode_identifiers,
                data_columns=barcode_identifiers.columns,
                format="table",
            )
            del barcode_identifiers

            # self.report_filter_stats()
            self.save_filter_stats()

    def make_plots(self):
        """
        Make plots for :py:class:`~seqlib.seqlib.BcidSeqLib` objects.

        Creates plot of the number of barcodes mapping to each identifier.
        """
        if self.plots_requested:
            SeqLib.make_plots(self)
            # open the PDF file
            pdf = PdfPages(os.path.join(self.plot_dir, "barcodes_per_identifier.pdf"))
            barcodemap_plot(self, pdf)
            pdf.close()
