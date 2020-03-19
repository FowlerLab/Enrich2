from .variant import VariantSeqLib
from .fqread import read_fastq, split_fastq_path
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
        self.logger = logging.getLogger("{}.{}".format(__name__, self.__class__))

    def configure(self, cfg):
        """
        Set up the object using the config object *cfg*, usually derived from
        a ``.json`` file.
        """
        VariantSeqLib.configure(self, cfg)
        self.logger = logging.getLogger(
            "{}.{} - {}".format(__name__, self.__class__.__name__, self.name)
        )

        # if counts are specified, copy them later
        # else handle the FASTQ config options and check the files
        if self.counts_file is None:
            self.configure_fastq(cfg)
            try:
                if split_fastq_path(self.reads) is None:
                    raise IOError(
                        "FASTQ file error: unrecognized extension "
                        "[{}]".format(self.name)
                    )
            except IOError as fqerr:
                raise IOError("FASTQ file error [{}]: {}".format(self.name, fqerr))

    def serialize(self):
        """
        Format this object (and its children) as a config object suitable for
        dumping to a config file.
        """
        cfg = VariantSeqLib.serialize(self)

        cfg["fastq"] = self.serialize_fastq()

        return cfg

    def configure_fastq(self, cfg):
        """
        Set up the object's FASTQ_ file handling and filtering options.
        """
        try:
            self.reads = cfg["fastq"]["reads"]

            if "reverse" in cfg["fastq"]:
                self.revcomp_reads = cfg["fastq"]["reverse"]
            else:
                self.revcomp_reads = False

            if "start" in cfg["fastq"]:
                self.trim_start = cfg["fastq"]["start"]
            else:
                self.trim_start = 1

            if "length" in cfg["fastq"]:
                self.trim_length = cfg["fastq"]["length"]
            else:
                self.trim_length = sys.maxsize

            self.filters = cfg["fastq"]["filters"]
        except KeyError as key:
            raise KeyError(
                "Missing required config value {key} [{name}]"
                "".format(key=key, name=self.name)
            )

    def serialize_fastq(self):
        """
        Serialize this object's FASTQ_ file handling and filtering options.
        """
        fastq = {"filters": self.serialize_filters()}
        fastq["reads"] = self.reads

        if self.revcomp_reads:
            fastq["reverse"] = True
        else:
            fastq["reverse"] = False

        if self.trim_start > 1:
            fastq["start"] = self.trim_start

        if self.trim_length < sys.maxsize:
            fastq["length"] = self.trim_length

        return fastq

    def counts_from_reads(self):
        """
        Reads the forward or reverse FASTQ_ file (reverse reads are
        reverse-complemented), performs quality-based filtering, and counts
        the variants.
        """
        df_dict = dict()

        self.logger.info("Counting variants")
        max_mut_variants = 0
        for fq in read_fastq(self.reads):
            fq.trim_length(self.trim_length, start=self.trim_start)
            if self.revcomp_reads:
                fq.revcomp()

            if self.read_quality_filter(fq):
                mutations = self.count_variant(fq.sequence)
                if mutations is None:  # too many mutations
                    max_mut_variants += 1
                    if self.report_filtered:
                        self.report_filtered_variant(fq.sequence, 1)
                else:
                    try:
                        df_dict[mutations] += 1
                    except KeyError:
                        df_dict[mutations] = 1

        self.save_counts("variants", df_dict, raw=True)
        del df_dict

        if self.aligner is not None:
            self.logger.info("Aligned {} variants".format(self.aligner.calls))
            self.aligner_cache = None
        self.logger.info(
            "Removed {} total variants with excess mutations"
            "".format(max_mut_variants)
        )
        self.save_filter_stats()

    def calculate(self):
        """
        Counts variants from counts file or FASTQ.
        """
        if not self.check_store("/main/variants/counts"):
            if not self.check_store("/raw/variants/counts"):
                if self.counts_file is not None:
                    self.counts_from_file(self.counts_file)
                else:
                    self.counts_from_reads()
            self.save_filtered_counts("variants", "count >= self.variant_min_count")
        self.count_synonymous()
