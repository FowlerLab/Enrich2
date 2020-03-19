import logging
from .seqlib import SeqLib
from .variant import VariantSeqLib
from .barcode import BarcodeSeqLib
from .barcodemap import BarcodeMap
import pandas as pd
from .plots import barcodemap_plot
from matplotlib.backends.backend_pdf import PdfPages
import os.path


class BcvSeqLib(VariantSeqLib, BarcodeSeqLib):
    """
    Class for counting variant data from barcoded sequencing libraries. 
    Creating a :py:class:`BcvSeqLib` requires a valid *config* 
    object with an ``'barcodes'`` entry and information about the wild type 
    sequence.

    The ``barcode_map`` keyword argument can be used to pass an existing 
    :py:class:`~seqlib.barcodemap.BarcodeMap`. Ensuring this is the 
    right :py:class:`~seqlib.barcodemap.BarcodeMap` is the responsibility 
    of the caller.
    """

    treeview_class_name = "Barcoded Variant SeqLib"

    def __init__(self):
        VariantSeqLib.__init__(self)
        BarcodeSeqLib.__init__(self)
        self.barcode_map = None
        self.logger = logging.getLogger("{}.{}".format(__name__, self.__class__))

    def configure(self, cfg, barcode_map=None):
        """
        Set up the object using the config object *cfg*, usually derived from 
        a ``.json`` file.
        """
        VariantSeqLib.configure(self, cfg)
        BarcodeSeqLib.configure(self, cfg)
        self.logger = logging.getLogger(
            "{}.{} - {}".format(__name__, self.__class__.__name__, self.name)
        )
        try:
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
                    cfg["barcodes"]["map file"], is_variant=True
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
        cfg = VariantSeqLib.serialize(self)
        cfg.update(BarcodeSeqLib.serialize(self))

        if self.barcode_map is not None:  # required for creating new objects in GUI
            cfg["barcodes"]["map file"] = self.barcode_map.filename

        return cfg

    def calculate(self):
        """
        Counts the barcodes using :py:meth:`BarcodeSeqLib.count` and combines them into 
        variant counts using the :py:class:`BarcodeMap`.
        """
        if not self.check_store("/main/variants/counts"):
            BarcodeSeqLib.calculate(self)  # count the barcodes
            df_dict = dict()
            barcode_variants = dict()

            self.logger.info("Converting barcodes to variants")
            # store mapped barcodes
            self.save_filtered_counts(
                "barcodes",
                "index in self.barcode_map.keys() & count >= self.barcode_min_count",
            )

            # count variants associated with the barcodes
            max_mut_barcodes = 0
            max_mut_variants = 0
            for bc, count in self.store["/main/barcodes/counts"].iterrows():
                count = count["count"]
                variant = self.barcode_map[bc]
                mutations = self.count_variant(variant)
                if mutations is None:  # variant has too many mutations
                    max_mut_barcodes += 1
                    max_mut_variants += count
                    if self.report_filtered:
                        self.report_filtered_variant(variant, count)
                else:
                    try:
                        df_dict[mutations] += count
                    except KeyError:
                        df_dict[mutations] = count
                    barcode_variants[bc] = mutations

            # save counts, filtering based on the min count
            self.save_counts(
                "variants",
                {k: v for k, v in df_dict.iteritems() if v >= self.variant_min_count},
                raw=False,
            )
            del df_dict

            # write the active subset of the BarcodeMap to the store
            barcodes = barcode_variants.keys()
            barcode_variants = pd.DataFrame(
                {"value": [barcode_variants[bc] for bc in barcodes]}, index=barcodes
            )
            del barcodes
            barcode_variants.sort_values("value", inplace=True)
            self.store.put(
                "/raw/barcodemap",
                barcode_variants,
                data_columns=barcode_variants.columns,
                format="table",
            )
            del barcode_variants

            if self.aligner is not None:
                self.logger.info("Aligned {} variants".format(self.aligner.calls))
                self.aligner_cache = None
            # self.report_filter_stats()
            self.logger.info(
                "Removed {} unique barcodes ({} total variants) "
                "with excess mutations".format(max_mut_barcodes, max_mut_variants)
            )
            self.save_filter_stats()

        self.count_synonymous()

    def make_plots(self):
        """
        Make plots for :py:class:`~seqlib.seqlib.BcvSeqLib` objects.

        Creates plot of the number of barcodes mapping to each variant.
        """
        if self.plots_requested:
            SeqLib.make_plots(self)
            # open the PDF file
            pdf = PdfPages(os.path.join(self.plot_dir, "barcodes_per_variant.pdf"))
            barcodemap_plot(self, pdf)
            pdf.close()
