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

from __future__ import absolute_import
from .variant import VariantSeqLib
from .fqread import read_fastq, split_fastq_path
import pandas as pd
import logging
import sys


class BasicSeqLib(VariantSeqLib):
    """
    Class for count data from sequencing libraries with a single read for 
    each variant. Creating a :py:class:`BasicSeqLib` requires a valid 
    *config* object, usually from a ``.json`` configuration file.
    """

    treeview_class_name = "Basic SeqLib"

    def __init__(self):
        VariantSeqLib.__init__(self)
        self.reads = None
        self.revcomp_reads = None
        self.trim_start = None
        self.trim_length = None


    def configure(self, cfg):
        """
        Set up the object using the config object *cfg*, usually derived from 
        a ``.json`` file.
        """
        VariantSeqLib.configure(self, cfg)

        # if counts are specified, copy them later
        # else handle the FASTQ config options and check the files
        if self.counts_file is None:
            self.configure_fastq(cfg)
            try:
                if split_fastq_path(self.reads) is None:
                    raise IOError("FASTQ file error: unrecognized extension [{}]".format(self.name)) 
            except IOError as fqerr:
                raise IOError("FASTQ file error [{}]: {}".format(self.name, fqerr))


    def serialize(self):
        """
        Format this object (and its children) as a config object suitable for dumping to a config file.
        """
        cfg = VariantSeqLib.serialize(self)

        cfg['fastq'] = self.serialize_fastq()

        return cfg


    def configure_fastq(self, cfg):
        """
        Set up the object's FASTQ_ file handling and filtering options.
        """
        try:
            self.reads = cfg['fastq']['reads']

            if 'reverse' in cfg['fastq']:
                self.revcomp_reads = cfg['fastq']['reverse']
            else:
                self.revcomp_reads = False

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
            raise KeyError("Missing required config value {key} [{name}]".format(key=key, name=self.name))


    def serialize_fastq(self):
        """
        Serialize this object's FASTQ_ file handling and filtering options.
        """
        fastq = {
            'filters' : self.serialize_filters()
        }
        fastq['reads'] = self.reads

        if self.revcomp_reads:
            fastq['reverse'] = True
        else:
            fastq['reverse'] = False

        if self.trim_start > 1:
            fastq['start'] = self.trim_start

        if self.trim_length < sys.maxsize:
            fastq['length'] = self.trim_length

        return fastq


    def calculate(self):
        """
        Reads the forward or reverse FASTQ file (reverse reads are reverse-complemented),
        performs quality-based filtering, and counts the variants.
        """
        if not self.check_store("/main/variants/counts"):
            if not self.check_store("/raw/variants/counts"):
                if self.counts_file is not None:
                    self.copy_raw()
                else:   # count everything
                    df_dict = dict()

                    logging.info("Counting variants", extra={'oname' : self.name})
                    for fq in read_fastq(self.reads):
                        fq.trim_length(self.trim_length, start=self.trim_start)
                        if self.revcomp_reads:
                            fq.revcomp()

                        if self.read_quality_filter(fq): # passed quality filtering
                            mutations = self.count_variant(fq.sequence)
                            if mutations is None: # read has too many mutations
                                self.filter_stats['max mutations'] += 1
                                self.filter_stats['total'] += 1
                                if self.report_filtered:
                                    self.report_filtered_read(fq, filter_flags)
                            else:
                                try:
                                    df_dict[mutations] += 1
                                except KeyError:
                                    df_dict[mutations] = 1

                    self.save_counts('variants', df_dict, raw=True)
                    del df_dict

                    if self.aligner is not None:
                        logging.info("Aligned {} variants".format(self.aligner.calls), extra={'oname' : self.name})
                        self.aligner_cache = None
                    #self.report_filter_stats()
                    self.save_filter_stats()

            self.save_filtered_counts('variants', "count >= self.variant_min_count")
        self.count_synonymous()


