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
import logging
import sys
from .seqlib import SeqLib
from .fqread import read_fastq, split_fastq_path


class BarcodeSeqLib(SeqLib):
    """
    Class for count data from barcoded sequencing libraries. Designed for
    barcode-only scoring or as a parent class for
    :py:class:`~seqlib.barcodevariant.BcvSeqLib` and
    :py:class:`~seqlib.barcodeid.BcidSeqLib`.
    """

    treeview_class_name = "Barcode SeqLib"

    def __init__(self):
        # Init step handled by VariantSeqLib's init for Barcode-variant
        if type(self).__name__ != "BcvSeqLib":
            SeqLib.__init__(self)
        self.reads = None
        self.revcomp_reads = None
        self.trim_start = None
        self.trim_length = None
        self.barcode_min_count = None
        self.add_label('barcodes')

    def configure(self, cfg):
        """
        Set up the object using the config object *cfg*, usually derived from
        a ``.json`` file.
        """
        SeqLib.configure(self, cfg)

        # handle non-FASTQ config options
        try:
            if 'min count' in cfg['barcodes']:
                self.barcode_min_count = int(cfg['barcodes']['min count'])
            else:
                self.barcode_min_count = 0
        except KeyError as key:
            raise KeyError("Missing required config value {}".format(key),
                           self.name)

        # if counts are specified, copy them later
        # else handle the FASTQ config options and check the files
        if self.counts_file is None:
            self.configure_fastq(cfg)
            try:
                if split_fastq_path(self.reads) is None:
                    raise ValueError("FASTQ file error: unrecognized file "
                                     "extension", self.name)
            except IOError as fqerr:
                raise IOError("FASTQ file error: {}".format(fqerr), self.name)

    def serialize(self):
        """
        Format this object (and its children) as a config object suitable for
        dumping to a config file.
        """
        cfg = SeqLib.serialize(self)

        cfg['barcodes'] = dict()
        if self.barcode_min_count > 0:
            cfg['barcodes']['min count'] = self.barcode_min_count

        cfg['fastq'] = self.serialize_fastq()

        return cfg

    def configure_fastq(self, cfg):
        """
        Set up the object's FASTQ_ file handling and filtering options.
        """
        try:
            self.reads = cfg['fastq']['reads']
            self.revcomp_reads = cfg['fastq']['reverse']

            if 'start' in cfg['fastq']:
                self.trim_start = cfg['fastq']['start']
            else:
                self.trim_start = 1

            if 'length' in cfg['fastq']:
                self.trim_length = cfg['fastq']['length']
            else:
                self.trim_length = sys.maxsize

            self.filters = cfg['fastq']['filters']
        except KeyError as key:
            raise KeyError("Missing required config value {}".format(key),
                           self.name)

    def serialize_fastq(self):
        """
        Serialize this object's FASTQ_ file handling and filtering options.
        """
        fastq = {
            'reads': self.reads,
            'reverse': self.revcomp_reads,
            'filters': self.serialize_filters()
        }
        if self.trim_start > 1:
            fastq['start'] = self.trim_start

        if self.trim_length < sys.maxsize:
            fastq['length'] = self.trim_length

        return fastq

    def counts_from_reads(self):
        """
        Reads the forward or reverse FASTQ_ file (reverse reads are
        reverse-complemented), performs quality-based filtering, and counts
        the barcodes.

        Barcode counts after read-level filtering are stored under
        ``"/raw/barcodes/counts"``.
        """
        df_dict = dict()

        filter_flags = dict()
        for key in self.filters:
            filter_flags[key] = False

        # count all the barcodes
        logging.info("Counting barcodes", extra={'oname': self.name})
        for fqr in read_fastq(self.reads):
            fqr.trim_length(self.trim_length, start=self.trim_start)
            if self.revcomp_reads:
                fqr.revcomp()

            if self.read_quality_filter(fqr):  # passed filtering
                try:
                    df_dict[fqr.sequence.upper()] += 1
                except KeyError:
                    df_dict[fqr.sequence.upper()] = 1

        self.save_counts('barcodes', df_dict, raw=True)
        del df_dict

    def calculate(self):
        """
        Counts the barcodes from the FASTQ file or from the provided counts
        file depending on the config.

        Barcodes that pass the minimum count
        filtering are stored under ``"/main/barcodes/counts"``.

        If ``"/main/barcodes/counts"`` already exists, those will be used
        instead of re-counting.
        """
        if self.check_store('/main/barcodes/counts'):
            return

        # no raw counts present
        if not self.check_store('/raw/barcodes/counts'):
            if self.counts_file is not None:
                self.counts_from_file(self.counts_file)
            else:
                self.counts_from_reads()

        if len(self.labels) == 1:  # only barcodes
            self.save_filtered_counts('barcodes',
                                      "count >= self.barcode_min_count")
            self.save_filter_stats()

