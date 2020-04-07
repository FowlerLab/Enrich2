#!/usr/bin/env python
#
from __future__ import print_function
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import logging
import json
import sys
import platform
import os.path

if platform.system() == "Darwin":
    # Explicitly set the backend to avoid the NSInvalidArgumentException when
    # running in GUI mode. Advanced users who want to use another matplotlib
    # backend when running in MacOS on the command line can modify this section
    # accordingly.
    import matplotlib

    matplotlib.use("TkAgg")
elif os.path.exists("/.dockerenv"):
    # Explicitly set the backend for running inside Docker. This may fail for
    # older versions of docker or alternative containerization tools such as
    # Singularity.
    import matplotlib

    matplotlib.use("Agg")
import enrich2.config_check as config_check
from enrich2.experiment import Experiment
from enrich2.selection import Selection
from enrich2.barcode import BarcodeSeqLib
from enrich2.barcodeid import BcidSeqLib
from enrich2.barcodevariant import BcvSeqLib
from enrich2.basic import BasicSeqLib
from enrich2.overlap import OverlapSeqLib
from enrich2.idonly import IdOnlySeqLib
from enrich2.storemanager import SCORING_METHODS, LOGR_METHODS
from enrich2.gui.configurator import Configurator
from enrich2.sfmap import parse_aa_list


__author__ = "Alan F Rubin"
__copyright__ = "Copyright 2016-2020, Alan F Rubin"
__license__ = "BSD-3-Clause"
__version__ = "1.3.1"
__maintainer__ = "Alan F Rubin"
__email__ = "alan.rubin@wehi.edu.au"


#: Name of the driver script. Used for logging output.
DRIVER_NAME = os.path.basename(sys.argv[0])


#: Format string for log entries (console or file).
LOG_FORMAT = "%(asctime)-15s [%(name)s] %(message)s"

#: Default log level
LOG_LEVEL = logging.INFO


def main_gui():
    """
    Entry point for GUI.

    """
    logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT)
    app = Configurator(__version__)
    app.mainloop()


def main_cmd():
    """
    Entry point for command line.

    """
    # build description string based on available methods
    desc_string = (
        "Command-line driver for Enrich2 v{}".format(__version__)
        + "\n\nscoring methods:\n"
        + "\n".join(["  {:22}{}".format(k, v) for k, v in SCORING_METHODS.items()])
        + "\n\nlog ratio methods:\n"
        + "\n".join(["  {:22}{}".format(k, v) for k, v in LOGR_METHODS.items()])
    )

    # create parser and add description
    parser = ArgumentParser(
        prog="Enrich2",
        description=desc_string,
        formatter_class=RawDescriptionHelpFormatter,
    )

    # add command line arguments
    parser.add_argument("config", help="JSON configuration file")
    parser.add_argument(
        "scoring_method", help="scoring method", choices=SCORING_METHODS.keys()
    )
    parser.add_argument(
        "logr_method", help="log ratio method", choices=LOGR_METHODS.keys()
    )

    # add support for semantic version checking
    parser.add_argument(
        "--version", action="version", version="{}".format(__version__)
    )

    # add analysis options
    parser.add_argument(
        "--log", metavar="FILE", dest="log_file", help="path to log file"
    )
    parser.add_argument(
        "--no-plots",
        dest="plots_requested",
        action="store_false",
        default=True,
        help="don't make plots",
    )
    parser.add_argument(
        "--no-tsv",
        dest="tsv_requested",
        action="store_false",
        default=True,
        help="don't generate tsv files",
    )
    parser.add_argument(
        "--recalculate",
        dest="force_recalculate",
        action="store_true",
        default=False,
        help="force recalculation",
    )
    parser.add_argument(
        "--component-outliers",
        dest="component_outliers",
        action="store_true",
        default=False,
        help="calculate component outlier stats",
    )
    parser.add_argument(
        "--output-dir",
        metavar="DIR",
        dest="output_dir_override",
        help="override the config file's output directory",
    )
    parser.add_argument(
        "--sfmap-aa-file",
        metavar="FILE",
        dest="sfmap_aa_file",
        help="amino acid groups for sequence-function maps",
    )

    args = parser.parse_args()

    # start the logs
    if args.log_file is not None:
        logging.basicConfig(filename=args.log_file, level=LOG_LEVEL, format=LOG_FORMAT)
    else:
        logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT)
    logger = logging.getLogger(__name__)

    # read the JSON file
    try:
        cfg = json.load(open(args.config, "U"))
    except IOError:
        raise IOError("Failed to open '{}' [{}]".format(args.config, DRIVER_NAME))
    except ValueError:
        raise ValueError("Improperly formatted .json file [{}]".format(DRIVER_NAME))

    # identify config file type and create the object
    if config_check.is_experiment(cfg):
        logger.info("Detected an Experiment config file")
        obj = Experiment()
    elif config_check.is_selection(cfg):
        logger.info("Detected a Selection config file")
        obj = Selection()
    elif config_check.is_seqlib(cfg):
        seqlib_type = config_check.seqlib_type(cfg)
        logger.info("Detected a %s config file", seqlib_type)
        if seqlib_type == "BarcodeSeqLib":
            obj = BarcodeSeqLib()
        elif seqlib_type == "BcidSeqLib":
            obj = BcidSeqLib()
        elif seqlib_type == "BcvSeqLib":
            obj = BcvSeqLib()
        elif seqlib_type == "BasicSeqLib":
            obj = BasicSeqLib()
        elif seqlib_type == "OverlapSeqLib":
            obj = OverlapSeqLib()
        elif seqlib_type == "IdOnlySeqLib":
            obj = IdOnlySeqLib()
        else:
            raise ValueError(
                "Unrecognized SeqLib type '{}' [{}]".format(seqlib_type, DRIVER_NAME)
            )
    else:
        raise ValueError("Unrecognized .json config [{}]".format(DRIVER_NAME))

    # set analysis options
    obj.force_recalculate = args.force_recalculate
    obj.component_outliers = args.component_outliers
    obj.scoring_method = args.scoring_method
    obj.logr_method = args.logr_method
    obj.plots_requested = args.plots_requested
    obj.tsv_requested = args.tsv_requested

    if args.output_dir_override is not None:
        obj.output_dir_override = True
        obj.output_dir = args.output_dir_override
    else:
        obj.output_dir_override = False

    if args.sfmap_aa_file is not None:
        obj.plot_options = dict()
        obj.plot_options["aa_list"], obj.plot_options[
            "aa_label_groups"
        ] = parse_aa_list(args.sfmap_aa_file)

    # configure the object
    obj.configure(cfg)

    # make sure objects are valid
    try:
        obj.validate()
    except ValueError:
        logger.exception("Invalid configuration")
    else:
        # open HDF5 files for the object and all child objects
        obj.store_open(children=True)

        # perform the analysis
        obj.calculate()

        # generate desired output
        obj.make_plots()
        try:
            obj.make_plots()
        except Exception:
            logger.exception("Calculations completed, but plotting failed.")
        try:
            obj.write_tsv()
        except Exception:
            logger.exception("Calculations completed, but TSV output failed.")

        # clean up
        obj.store_close(children=True)


if __name__ == "__main__":
    gui_mode = False

    try:
        if sys.argv[1] == "gui":
            gui_mode = True
    except IndexError:
        pass

    if gui_mode:
        main_gui()
    else:
        main_cmd()
