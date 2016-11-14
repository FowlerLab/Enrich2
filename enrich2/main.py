#!/usr/bin/env python
#
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
import argparse
import logging
import json
import sys
import os.path
import enrich2.config_check as config_check
from enrich2.experiment import Experiment
from enrich2.selection import Selection
from enrich2.basic import BasicSeqLib
from enrich2.barcodevariant import BcvSeqLib
from enrich2.barcode import BarcodeSeqLib
from enrich2.overlap import OverlapSeqLib
from enrich2.storemanager import available_scoring_methods, available_logr_methods
from enrich2.gui.configurator import Configurator

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
    desc_string = "Command-line driver for Enrich2." + \
        "\n\nscoring methods:\n" + \
        "\n".join(["  {:22}{}".format(k, v) for k, v in available_scoring_methods.items()]) + \
        "\n\nlog ratio methods:\n" + \
        "\n".join(["  {:22}{}".format(k, v) for k, v in available_logr_methods.items()]) 

    # create parser and add description
    parser = argparse.ArgumentParser(description=desc_string, formatter_class=argparse.RawDescriptionHelpFormatter)

    # add command line arguments
    parser.add_argument("config", help="JSON configuration file")
    parser.add_argument("scoring_method", help="scoring method", choices=available_scoring_methods.keys())
    parser.add_argument("logr_method", help="log ratio method", choices=available_logr_methods.keys())

    # add analysis options
    parser.add_argument("--log", dest="log_file", metavar="FILE", help="path to log file")
    parser.add_argument("--no-plots", help="don't make plots", dest="plots_requested", action="store_false", default=True)
    parser.add_argument("--no-tsv", help="don't generate tsv files", dest="tsv_requested", action="store_false", default=True)
    parser.add_argument("--recalculate", help="force recalculation", dest="force_recalculate", action="store_true", default=False)
    parser.add_argument("--component-outliers", help="calculate component outlier stats", dest="component_outliers", action="store_true", default=False)
    parser.add_argument("--output-dir", help="override the config file's output directory", dest="output_dir_override", metavar="DIR")
    
    args = parser.parse_args()

    # start the logs
    start_logging(args.log_file, logging.DEBUG)

    # read the JSON file
    try:
        cfg = json.load(open(args.config, "U"))
    except IOError:
        raise IOError("Failed to open '{}' [{}]".format(args.config, DRIVER_NAME))
    except ValueError:
        raise ValueError("Improperly formatted .json file [{}]".format(DRIVER_NAME))

    # identify config file type and create the object
    if config_check.is_experiment(cfg):
        logging.info("Detected an Experiment config file", extra={'oname' : DRIVER_NAME})
        obj = Experiment()
    elif config_check.is_selection(cfg):
        logging.info("Detected a Selection config file", extra={'oname' : DRIVER_NAME})
        obj = Selection()
    elif config_check.is_seqlib(cfg):
        seqlib_type = config_check.seqlib_type(cfg)
        logging.info("Detected a {} config file".format(seqlib_type), extra={'oname' : DRIVER_NAME})
        obj = globals()[seqlib_type]()
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

    # configure the object
    obj.configure(cfg)

    # make sure objects are valid
    try:
        obj.validate()
    except ValueError, e:
        logging.error("Invalid settings: {}".format(e), extra={'oname' : DRIVER_NAME})
    else:
        # open HDF5 files for the object and all child objects
        obj.store_open(children=True)

        # perform the analysis
        obj.calculate()
        
        # generate desired output
        obj.make_plots()
        try:
            obj.make_plots()
        except Exception, e:
            logging.warning("Calculations completed, but plotting failed: {}".format(e), extra={'oname' : DRIVER_NAME})
        try:
            obj.write_tsv()
        except Exception, e:
            logging.warning("Calculations completed, but tsv output failed: {}".format(e), extra={'oname' : DRIVER_NAME})

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
