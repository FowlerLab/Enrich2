from __future__ import print_function
import pandas as pd
import logging
from matplotlib.backends.backend_pdf import PdfPages
import os.path
from .plots import overlap_merge_plot
from .seqlib import SeqLib
from .variant import VariantSeqLib
from .fqread import read_fastq_multi, split_fastq_path, FQRead


class OverlapSeqLib(VariantSeqLib):
    """
    Class for count data from sequencing libraries with overlapping paired-end 
    reads for each variant. Creating a 
    :py:class:`~seqlib.overlap.OverlapSeqLib` requires a valid *config* object 
    with an ``'overlap'`` entry.

    The ``"fastq"`` config entry must contain two read files, with the keys 
    ``"forward"`` and ``"reverse"``. Information about how to combine these 
    reads is in the ``"overlap"`` config entry.

    The ``"overlap"`` config entry contains the following keys:

    * ``"forward start"`` --- position in the forward read where the \
        overlapping region begins
    * ``"reverse start"`` --- position in the reverse read where the \
        overlapping region begins (before being reverse-complemented)
    * ``"length"`` --- number of bases in the overlapping region
    * ``"max mismatches"`` --- maximum number of mismatches tolerated in the \
        overlapping region before discarding the read
    * ``"overlap only"`` --- whether to trim the merged read to contain only \
        the overlapping region (optional, defaults to ``False``)

    Here is a schematic of the case in the above JSON example::

        forward ---> 1   
                     CGACGCAAGGA
                       |||||||||
                       ACTCCTTGCGTCG
                                   1 <--- reverse

    Note that the merged sequence is identical to the wild type sequence given 
    in the JSON file.
    """

    treeview_class_name = "Overlap SeqLib"

    def __init__(self):
        VariantSeqLib.__init__(self)
        self.forward = None
        self.reverse = None
        self.fwd_start = None
        self.rev_start = None
        self.overlap_length = None
        self.trim = None
        self.max_overlap_mismatches = None
        self.merge_mismatches = None
        self.default_filters.update({"merge failure": True})
        self.default_filters.update({"remove unresolvable": False})
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
                self.fwd_start = int(cfg["overlap"]["forward start"])
                self.rev_start = int(cfg["overlap"]["reverse start"])
                self.overlap_length = int(cfg["overlap"]["length"])
                self.trim = cfg["overlap"]["trim"]
                self.max_overlap_mismatches = int(cfg["overlap"]["max mismatches"])

                forward_error = False
                reverse_error = False
                if split_fastq_path(self.forward) is None:
                    forward_error = True
                if split_fastq_path(self.reverse) is None:
                    reverse_error = True
                if forward_error and reverse_error:
                    raise IOError(
                        "FASTQ file error: unrecognized extension (forward and reverse) [{}]".format(
                            self.name
                        )
                    )
                elif forward_error:
                    raise IOError(
                        "FASTQ file error: unrecognized extension (forward) [{}]".format(
                            self.name
                        )
                    )
                elif reverse_error:
                    raise IOError(
                        "FASTQ file error: unrecognized extension (reverse) [{}]".format(
                            self.name
                        )
                    )
            except IOError as fqerr:
                raise IOError("FASTQ file error [{}]: {}".format(self.name, fqerr))
            except KeyError as key:
                raise KeyError(
                    "Missing required config value {key} [{name}]".format(
                        key=key, name=self.name
                    )
                )
            except ValueError as value:
                raise ValueError(
                    "Invalid parameter value {value} [{name}]".format(
                        value=value, name=self.name
                    )
                )

    def serialize(self):
        """
        Format this object (and its children) as a config object suitable for dumping to a config file.
        """
        cfg = VariantSeqLib.serialize(self)

        cfg["fastq"] = self.serialize_fastq()
        cfg["overlap"] = {
            "forward start": self.fwd_start,
            "reverse start": self.rev_start,
            "length": self.overlap_length,
            "trim": self.trim,
            "max mismatches": self.max_overlap_mismatches,
        }

        return cfg

    def configure_fastq(self, cfg):
        """
        Set up the object's FASTQ_ file handling and filtering options.
        """
        try:
            self.forward = cfg["fastq"]["forward reads"]
            self.reverse = cfg["fastq"]["reverse reads"]

            if "merge failure" in cfg["fastq"]["filters"]:
                raise ValueError(
                    "'merge failure' is not user-configurable [{}]".format(self.name)
                )
            self.filters = cfg["fastq"]["filters"]
        except KeyError as key:
            raise KeyError(
                "Missing required config value {key} [{name}]".format(
                    key=key, name=self.name
                )
            )

    def serialize_fastq(self):
        """
        Serialize this object's FASTQ_ file handling and filtering options.
        """
        fastq = {
            "forward reads": self.forward,
            "reverse reads": self.reverse,
            "filters": self.serialize_filters(),
        }

        return fastq

    def merge_reads(self, fwd, rev):
        """
        Combines the *fwd* and *rev* :py:class:`~fqread.FQRead` objects into a 
        single :py:class:`~fqread.FQRead` with the same header information as 
        *fwd*. Mismatches are resolved by taking the highest quality base. If 
        discrepant bases have the same quality value, this position is 
        unresolvable and an ``'X'`` is inserted. Quality values in the 
        resulting :py:class:`~fqread.FQRead` are the maximum quality for the 
        given base at that position. Returns ``None`` if the maximum number of 
        mismatches in the overlap region is exceded.
        """
        rev.revcomp()

        # print(fwd.sequence, "-" * (self.rev_start - 1), sep="")
        # print("-" * (self.fwd_start - 1), rev.sequence, sep="")
        rev_extra_start = len(rev) - self.rev_start + 1
        fwd_end = self.fwd_start + self.overlap_length - 1
        merge = FQRead(
            header=fwd.header + "|" + rev.header,
            sequence="A",
            header2=fwd.header2 + "|" + rev.header2,
            quality="#",
            qbase=fwd.qbase,
        )
        merge.sequence = list(fwd.sequence[:fwd_end] + rev.sequence[rev_extra_start:])
        merge.quality = fwd.quality[:fwd_end] + rev.quality[rev_extra_start:]

        mismatches = 0
        first = True
        for i in xrange(self.overlap_length):
            a = self.fwd_start - 1 + i
            b = len(rev) - self.rev_start - self.overlap_length + i + 1
            try:
                if fwd.sequence[a] == rev.sequence[b]:
                    # take the highest quality value
                    if rev.quality[b] > fwd.quality[a]:
                        merge.quality[a] = rev.quality[b]
                else:
                    if fwd.quality[a] == rev.quality[b]:
                        merge.sequence[a] = "X"  # unresolvable
                        self.merge_mismatches.iloc[i]["unresolved"] += 1
                    elif rev.quality[b] > fwd.quality[a]:
                        merge.sequence[a] = rev.sequence[b]
                        merge.quality[a] = rev.quality[b]
                        self.merge_mismatches.iloc[i]["resolved"] += 1
                    else:
                        # overlap region already same as fwd
                        self.merge_mismatches.iloc[i]["resolved"] += 1
                    mismatches += 1
                    if first:
                        self.merge_mismatches.iloc[i]["first"] += 1
                        first = False
            except IndexError:
                raise IndexError(
                    "Failed to calculate overlap (a={a}, len(a)={lena}, b={b}, len(b)={lenb}) [{name}]".format(
                        a=a,
                        b=b,
                        lena=len(fwd.sequence),
                        lenb=len(rev.sequence),
                        name=self.name,
                    )
                )

        if mismatches > self.max_overlap_mismatches:
            return None  # merge failed

        merge.sequence = "".join(merge.sequence)
        if self.trim:
            merge.trim_length(self.overlap_length, self.fwd_start)
        return merge

    def counts_from_reads(self):
        df_dict = dict()

        self.merge_mismatches = pd.DataFrame(
            data=0,
            index=[
                x + self.fwd_start + self.wt.dna_offset
                for x in xrange(0, self.overlap_length)
            ],
            columns=["resolved", "unresolved", "first"],
        )

        self.logger.info("Counting variants")
        max_mut_variants = 0
        for fwd, rev in read_fastq_multi([self.forward, self.reverse]):
            # filter on chastity before merge
            chaste = True
            if self.filters["chastity"]:
                if not fwd.is_chaste():
                    chaste = False
                    if self.report_filtered:
                        self.report_filtered_read(fwd, "chastity")
                if not rev.is_chaste():
                    chaste = False
                    if self.report_filtered:
                        self.report_filtered_read(rev, "chastity")
                if not chaste:
                    self.filter_stats["chastity"] += 1
                    self.filter_stats["total"] += 1
                    continue

            merge = self.merge_reads(fwd, rev)
            if merge is None:  # merge failed
                self.filter_stats["merge failure"] += 1
                self.filter_stats["total"] += 1
                if self.report_filtered:
                    self.report_filtered_read(fwd, filter_flags)
                    self.report_filtered_read(rev, filter_flags)
            else:
                if self.read_quality_filter(merge):
                    mutations = self.count_variant(merge.sequence)
                    if mutations is None:  # merge read has too many mutations
                        max_mut_variants += 1
                        if self.report_filtered:
                            self.report_filtered_variant(merge.sequence, 1)
                    else:
                        try:
                            df_dict[mutations] += 1
                        except KeyError:
                            df_dict[mutations] = 1

        self.store.put(
            "/raw/overlap_mismatches",
            self.merge_mismatches,
            format="table",
            data_columns=self.merge_mismatches.columns,
        )
        self.merge_mismatches = None
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
        Reads the forward and reverse reads, merges them, performs 
        quality-based filtering, and counts the variants.
        """
        if not self.check_store("/main/variants/counts"):
            if not self.check_store("/raw/variants/counts"):
                if self.counts_file is not None:
                    self.counts_from_file(self.counts_file)
                else:  # count everything
                    self.counts_from_reads()
            self.save_filtered_counts("variants", "count >= self.variant_min_count")

        self.count_synonymous()

    def make_plots(self):
        """
        Make plots for :py:class:`~seqlib.seqlib.OverlapSeqLib` objects.

        Creates plots of the location of merged read mismatches.
        """
        if self.plots_requested:
            SeqLib.make_plots(self)
            pdf = PdfPages(os.path.join(self.plot_dir, "overlap_mismatches.pdf"))
            overlap_merge_plot(self, pdf)
            pdf.close()
