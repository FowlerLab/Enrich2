#!/usr/bin/env python
#
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
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import logging
import json
import sys
import os.path
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
__copyright__ = "Copyright 2016-2017, Alan F Rubin"
__license__ = "GPLv3"
__version__ = "1.1.0"
__maintainer__ = "Alan F Rubin"
__email__ = "alan.rubin@wehi.edu.au"


#: Name of the driver script. Used for logging output.
DRIVER_NAME = os.path.basename(sys.argv[0])


#: Format string for log entries (console or file).
LOG_FORMAT = "%(asctime)-15s [%(oname)s] %(message)s"


def start_logging(log_file, log_level):
    """
    Begin logging. This function should be called by the driver at the start of program execution. 
    Message format is defined by :py:const:`LOG_FORMAT`.

    Args:
        log_file (str): Name of the log output file, or ``None`` to output to console.

        log_level: Requested logging level. 
            See :py:class:`~logging.Logger` for a detailed description of the options. Most program 
            status messages are output at the ``INFO`` level.

    """
    if log_file is not None:
        logging.basicConfig(filename=log_file, level=log_level, format=LOG_FORMAT)
    else:
        logging.basicConfig(level=log_level, format=LOG_FORMAT)


def main_gui():
    """
    Entry point for GUI.

    """
    start_logging(None, logging.DEBUG)
    app = Configurator()
    app.mainloop()


def main_cmd():
    """
    Entry point for command line.

    """
    # build description string based on available methods
    desc_string = "Command-line driver for Enrich2 v{}".format(__version__) + \
        "\n\nscoring methods:\n" + \
        "\n".join(["  {:22}{}".format(k, v) for k, v in
                   SCORING_METHODS.items()]) + \
        "\n\nlog ratio methods:\n" + \
        "\n".join(["  {:22}{}".format(k, v) for k, v in
                   LOGR_METHODS.items()])

    # create parser and add description
    parser = ArgumentParser(prog="Enrich2", description=desc_string,
                            formatter_class=RawDescriptionHelpFormatter)

    # add command line arguments
    parser.add_argument("config", help="JSON configuration file")
    parser.add_argument("scoring_method", help="scoring method",
                        choices=SCORING_METHODS.keys())
    parser.add_argument("logr_method", help="log ratio method",
                        choices=LOGR_METHODS.keys())

    # add support for semantic version checking
    parser.add_argument("--version", action="version",
                        version="%(prog)s {}".format(__version__))

    # add analysis options
    parser.add_argument("--log", metavar="FILE", dest="log_file",
                        help="path to log file")
    parser.add_argument("--no-plots", dest="plots_requested",
                        action="store_false", default=True,
                        help="don't make plots")
    parser.add_argument("--no-tsv", dest="tsv_requested",
                        action="store_false", default=True,
                        help="don't generate tsv files")
    parser.add_argument("--recalculate", dest="force_recalculate",
                        action="store_true", default=False,
                        help="force recalculation")
    parser.add_argument("--component-outliers", dest="component_outliers",
                        action="store_true", default=False,
                        help="calculate component outlier stats")
    parser.add_argument("--output-dir", metavar="DIR",
                        dest="output_dir_override",
                        help="override the config file's output directory")
    parser.add_argument("--sfmap-aa-file", metavar="FILE",
                        dest="sfmap_aa_file",
                        help="amino acid groups for sequence-function maps")

    args = parser.parse_args()

    # start the logs
    start_logging(args.log_file, logging.DEBUG)

    # read the JSON file
    try:
        cfg = json.load(open(args.config, "U"))
    except IOError:
        raise IOError("Failed to open '{}' [{}]".format(
            args.config, DRIVER_NAME))
    except ValueError:
        raise ValueError("Improperly formatted .json file [{}]".format(
            DRIVER_NAME))

    # identify config file type and create the object
    if config_check.is_experiment(cfg):
        logging.info("Detected an Experiment config file",
                     extra={'oname': DRIVER_NAME})
        obj = Experiment()
    elif config_check.is_selection(cfg):
        logging.info("Detected a Selection config file",
                     extra={'oname': DRIVER_NAME})
        obj = Selection()
    elif config_check.is_seqlib(cfg):
        seqlib_type = config_check.seqlib_type(cfg)
        logging.info("Detected a %s config file", seqlib_type,
                     extra={'oname': DRIVER_NAME})
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
            raise ValueError("Unrecognized SeqLib type '{}' [{}]".format(
                seqlib_type, DRIVER_NAME))
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
        obj.plot_options['aa_list'], obj.plot_options['aa_label_groups'] = \
            parse_aa_list(args.sfmap_aa_file)

    # configure the object
    obj.configure(cfg)

    # make sure objects are valid
    try:
        obj.validate()
    except ValueError:
        logging.exception("Invalid configuration",
                          extra={'oname': DRIVER_NAME})
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
            logging.exception("Calculations completed, but plotting failed.",
                              extra={'oname': DRIVER_NAME})
        try:
            obj.write_tsv()
        except Exception:
            logging.exception("Calculations completed, but TSV ouput failed.",
                              extra={'oname': DRIVER_NAME})

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
