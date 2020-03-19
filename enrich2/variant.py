from __future__ import print_function
import re
from .aligner import Aligner
from .seqlib import SeqLib
from .wildtype import WildTypeSequence
import logging
from .constants import WILD_TYPE_VARIANT, SYNONYMOUS_VARIANT
from .constants import CODON_TABLE, AA_CODES

#: Default number of maximum mutation.
#: Must be set to avoid data frame performance errors.
DEFAULT_MAX_MUTATIONS = 10

#: Matches a single amino acid substitution in HGVS_ format.
re_protein = re.compile(
    "(?P<match>p\.(?P<pre>[A-Z][a-z][a-z])(?P<pos>-?\d+)" "(?P<post>[A-Z][a-z][a-z]))"
)

#: Matches a single nucleotide substitution (coding or noncoding)
#: in HGVS_ format.
re_nucleotide = re.compile(
    "(?P<match>[nc]\.(?P<pos>-?\d+)(?P<pre>[ACGT])>(?P<post>[ACGT]))"
)

#: Matches a single coding nucleotide substitution in HGVS_ format.
re_coding = re.compile(
    "(?P<match>c\.(?P<pos>-?\d+)(?P<pre>[ACGT])>(?P<post>[ACGT]) "
    "\(p\.(?:=|[A-Z][a-z][a-z]-?\d+[A-Z][a-z][a-z])\))"
)

#: Matches a single noncoding nucleotide substitution in HGVS_ format.
re_noncoding = re.compile(
    "(?P<match>n\.(?P<pos>-?\d+)(?P<pre>[ACGT])>(?P<post>[ACGT]))"
)


def valid_variant(s, is_coding=True):
    """
    Returns True if s is a valid coding or noncoding variant, else False.
    """
    if s == WILD_TYPE_VARIANT:
        return True
    else:
        if is_coding:
            for mut in s.split(", "):
                match = re_coding.match(mut)
                if match is None:
                    return False
            return True
        else:
            for mut in s.split(", "):
                match = re_noncoding.match(mut)
                if match is None:
                    return False
            return True


def hgvs2single(s):
    """
    Convert HGVS string from Enrich2 to <pre><pos><post>
    single-letter amino acid variant string.

    Returns a list of single-letter variants.
    """
    t = re_protein.findall(s)
    return ["{}{}{}".format(AA_CODES[m[1]], m[2], AA_CODES[m[3]]) for m in t]


def single2hgvs(s):
    """
    Convert single-letter amino acid changes in the form
    <pre><pos><post> into HGVS strings that match Enrich2 
    output.

    Searches the string s for all instances of the above
    pattern and returns a list of Enrich2 variants.
    """
    t = re.findall("[A-Z*]\d+[A-Z*]", s)
    return ["p.{}{}{}".format(AA_CODES[x[0]], x[1:-1], AA_CODES[x[-1]]) for x in t]


def get_variant_type(variant):
    """
    Use regular expressions to determine whether the variant is protein,
    coding, or noncoding.

    :param str variant: variant string
    :return: ``'protein'``, ``'coding'``, or ``'noncoding'`` depending on \
    which regular expression matches, else ``None``. Note that both wild type \
    and synonymous special variants will return ``None``.
    :rtype: str
    """
    v = variant.split(", ")[0]  # test first token of multi-mutant
    if re_protein.match(v) is not None:
        return "protein"
    elif re_coding.match(v) is not None:
        return "coding"
    elif re_noncoding.match(v) is not None:
        return "noncoding"
    else:
        return None


def mutation_count(variant):
    """
    Counts the number of mutations in the HGVS_ *variant*.

    :param str variant: variant string
    :return: number of variants (0 if wild type or synonymous)
    :rtype: int
    """
    if len(variant) == 0:
        raise ValueError("Empty variant string.")
    elif variant == WILD_TYPE_VARIANT:
        return 0
    elif variant == SYNONYMOUS_VARIANT:
        return 0
    else:
        return len(variant.split(","))


def has_indel(variant):
    """
    Tests if the HGVS_ *variant* contains an indel mutation.

    :param str variant: variant string
    :return: ``True`` if there is an indel, else ``False``
    :rtype: bool
    """
    if len(variant) == 0:
        raise ValueError("Empty variant string.")
    elif variant == WILD_TYPE_VARIANT:
        return False
    elif variant == SYNONYMOUS_VARIANT:
        return False
    else:
        return any(x in variant for x in ("ins", "del", "dup"))


def has_unresolvable(variant):
    """
    Tests if the HGVS_ *variant* has an unresolvable amino acid change.

    Unresolvable amino acid changes are most commonly caused by the presence
    of N or X nucleotides, resulting in a non-translatable codon. They are
    also found when a frameshift causes the last part of the coding sequence
    to not have three nucleotides.

    :param str variant: variant string
    :return: ``True`` if there is an unresolvable change, else ``False``
    :rtype: bool
    """
    if AA_CODES["?"] in variant:
        return True
    else:
        return False


def protein_variant(variant):
    """
    Return an HGVS_ variant string containing only the protein changes in a
    coding HGVS_ variant string. If all variants are synonymous, returns the
    synonymous variant code. If the variant is wild type, returns the wild
    type variant.

    :param str variant: coding variant string
    :return: protein variant string (or synonymous or wild type)
    :rtype: str
    """
    if len(variant) == 0:
        raise ValueError("Empty variant string.")
    elif variant == WILD_TYPE_VARIANT:
        return WILD_TYPE_VARIANT
    elif variant == SYNONYMOUS_VARIANT:
        return SYNONYMOUS_VARIANT
    else:
        matches = re.findall("\((p\.\S*)\)", variant)
        if len(matches) == 0:
            raise ValueError("Invalid coding variant string.")
        # uniqify and remove synonymous
        seen = {"p.=": True}
        unique_matches = list()
        for v in matches:
            if v in seen:
                continue
            else:
                seen[v] = True
                unique_matches.append(v)
        if len(unique_matches) == 0:
            return SYNONYMOUS_VARIANT
        else:
            return ", ".join(unique_matches)


class VariantSeqLib(SeqLib):
    """
    Abstract :py:class:`~seqlib.seqlib.SeqLib` class for for Enrich libraries
    containing variants. Implements core functionality for assessing variants,
    either coding or noncoding. Subclasess must evaluate the variant DNA
    sequences that are being counted.
    """

    def __init__(self):
        SeqLib.__init__(self)
        self.wt = WildTypeSequence(self.name)
        self.aligner = None
        self.aligner_cache = None
        self.variant_min_count = None
        self.max_mutations = None
        # 'synonymous' label may be added in configure() if wt is coding
        self.add_label("variants")
        self.logger = logging.getLogger("{}.{}".format(__name__, self.__class__))

    def configure(self, cfg):
        """
        Set up the object using the config object *cfg*, usually derived from
        a ``.json`` file.
        """
        SeqLib.configure(self, cfg)
        self.logger = logging.getLogger(
            "{}.{} - {}".format(__name__, self.__class__.__name__, self.name)
        )

        self.wt.configure(cfg["variants"]["wild type"])
        if self.is_coding():
            self.add_label("synonymous")
        try:
            if "use aligner" in cfg["variants"]:
                if cfg["variants"]["use aligner"]:
                    self.aligner = Aligner()
                    self.aligner_cache = dict()
                else:
                    self.aligner = None
                    self.aligner_cache = None

            if "min count" in cfg["variants"]:
                self.variant_min_count = int(cfg["variants"]["min count"])
            else:
                self.variant_min_count = 0

            if "max mutations" in cfg["variants"]:
                self.max_mutations = int(cfg["variants"]["max mutations"])
            else:
                self.max_mutations = DEFAULT_MAX_MUTATIONS

        except KeyError as key:
            raise KeyError(
                "Missing required config value {key} [{name}]"
                "".format(key=key, name=self.name)
            )

    def serialize(self):
        """
        Format this object (and its children) as a config object suitable for
        dumping to a config file.
        """
        cfg = SeqLib.serialize(self)

        cfg["variants"] = dict()
        cfg["variants"]["wild type"] = self.wt.serialize()
        cfg["variants"]["use aligner"] = self.aligner is not None
        if self.max_mutations != DEFAULT_MAX_MUTATIONS:
            cfg["variants"]["max mutations"] = self.max_mutations
        if self.variant_min_count > 0:
            cfg["variants"]["min count"] = self.variant_min_count

        return cfg

    def is_coding(self):
        """
        Return ``True`` if the variants are protein-coding, else ``False``.
        """
        return self.wt.is_coding()

    def has_wt_sequence(self):
        """
        Returns ``True``, because :py:class:`~seqlib.seqlib.VariantSeqLib`
        objects have a wild type sequence. Raises a ValueError if the wild type
        sequence is not set properly.
        """
        if self.wt is not None:
            return True
        else:
            raise ValueError("Wild type not set properly [{}]".format(self.name))

    def align_variant(self, variant_dna):
        """
        Use the local :py:class:`~seqlib.aligner.Aligner` instance to align the
        *variant_dna* to the wild type sequence. Returns a list of HGVS_
        variant strings.

        Aligned variants are stored in a local dictionary to avoid recomputing
        alignments. This dictionary should be cleared after all variants are
        counted, to save memory.

        .. warning:: Using the :py:class:`~seqlib.aligner.Aligner` \
        dramatically increases runtime.
        """
        if variant_dna in self.aligner_cache.keys():
            return self.aligner_cache[variant_dna]

        mutations = list()
        traceback = self.aligner.align(self.wt.dna_seq, variant_dna)
        for x, y, cat, length in traceback:
            if cat == "match":
                continue
            elif cat == "mismatch":
                mut = "{pre}>{post}".format(pre=self.wt.dna_seq[x], post=variant_dna[y])
            elif cat == "insertion":
                if y > length:
                    dup = variant_dna[y : y + length]
                    if dup == variant_dna[y - length : y]:
                        mut = "dup{seq}".format(seq=dup)
                    else:
                        mut = "_{pos}ins{seq}".format(pos=x + 2, seq=dup)
                else:
                    mut = "_{pos}ins{seq}".format(
                        pos=x + 2, seq=variant_dna[y : y + length]
                    )
            elif cat == "deletion":
                mut = "_{pos}del".format(pos=x + length)
            else:
                raise ValueError("Unable to resolve mutation [{}]".format(self.name))
            mutations.append((x, mut))

        self.aligner_cache[variant_dna] = mutations
        return mutations

    def count_variant(self, variant_dna, include_indels=True):
        """
        Identifies mutations and counts the *variant_dna* sequence.
        The algorithm attempts to call variants by comparing base-by-base.
        If the *variant_dna* and wild type DNA are different lengths, or if
        there are an excess of mismatches (indicating a possible indel), local
        alignment is performed using :py:meth:`align_variant` if this option
        has been selected in the configuration.

        Each variant is stored as a tab-delimited string of mutations in HGVS_
        format. Returns a list of HGVS_ variant strings. Returns an empty list
        if the variant is wild type. Returns None if the variant was discarded
        due to excess mismatches.
        """
        if not re.match("^[ACGTNXacgtnx]+$", variant_dna):
            raise ValueError(
                "Variant DNA sequence contains unexpected "
                "characters [{}]".format(self.name)
            )

        variant_dna = variant_dna.upper()

        if len(variant_dna) != len(self.wt.dna_seq):
            if self.aligner is not None:
                mutations = self.align_variant(variant_dna)
            else:
                return None
        else:
            mutations = list()
            for i in xrange(len(variant_dna)):
                if variant_dna[i] != self.wt.dna_seq[i]:
                    mutations.append(
                        (
                            i,
                            "{pre}>{post}".format(
                                pre=self.wt.dna_seq[i], post=variant_dna[i]
                            ),
                        )
                    )
                    if len(mutations) > self.max_mutations:
                        if self.aligner is not None:
                            mutations = self.align_variant(variant_dna)
                            if len(mutations) > self.max_mutations:
                                # too many mutations post-alignment
                                return None
                            else:
                                # stop looping over this variant
                                break
                        else:
                            # too many mutations and not using aligner
                            return None

        mutation_strings = list()
        if self.is_coding():
            variant_protein = ""
            for i in xrange(0, len(variant_dna), 3):
                try:
                    variant_protein += CODON_TABLE[variant_dna[i : i + 3]]
                except KeyError:  # garbage codon due to indel, X, or N
                    variant_protein += "?"

            for pos, change in mutations:
                ref_dna_pos = pos + self.wt.dna_offset + 1
                ref_pro_pos = pos / 3 + self.wt.protein_offset + 1
                mut = "c.{pos}{change}".format(pos=ref_dna_pos, change=change)
                if has_indel(change):
                    mut += " (p.{pre}{pos}fs)".format(
                        pre=AA_CODES[self.wt.protein_seq[pos / 3]], pos=ref_pro_pos
                    )
                elif variant_protein[pos / 3] == self.wt.protein_seq[pos / 3]:
                    mut += " (p.=)"
                else:
                    mut += " (p.{pre}{pos}{post})".format(
                        pre=AA_CODES[self.wt.protein_seq[pos / 3]],
                        pos=ref_pro_pos,
                        post=AA_CODES[variant_protein[pos / 3]],
                    )
                mutation_strings.append(mut)
        else:
            for pos, change in mutations:
                ref_dna_pos = pos + self.wt.dna_offset + 1
                mut = "n.{pos}{change}".format(pos=ref_dna_pos, change=change)
                mutation_strings.append(mut)

        if len(mutation_strings) > 0:
            variant_string = ", ".join(mutation_strings)
        else:
            variant_string = WILD_TYPE_VARIANT
        return variant_string

    def count_synonymous(self):
        """
        Combine counts for synonymous variants (defined as variants that differ
        at the nucleotide level but not at the amino acid level) and store them 
        under the ``synonymous`` label.

        This method should be called only after ``variants`` have been counted.

        .. note:: The total number of ``synonymous`` variants may be greater \
        than the total number of ``variants`` after filtering. This is \
        because ``variants`` are combined into ``synonymous`` entries at the \
        :py:class:`~seqlib.seqlib.SeqLib` level before count-based filtering, \
        allowing filtered ``variants`` to contribute counts to their \
        ``synonymous`` entry.
        """
        if not self.is_coding():
            self.logger.warning("Cannot count synonymous mutations in noncoding data")
            return

        if self.check_store("/main/synonymous/counts"):
            return

        self.logger.info("Counting synonymous variants")

        df_dict = dict()

        for variant, count in self.store["/main/variants/counts"].iterrows():
            if variant == WILD_TYPE_VARIANT:
                df_dict[variant] = count["count"]
            else:
                variant = protein_variant(variant)
                if len(variant) == 0:
                    variant = SYNONYMOUS_VARIANT
                try:
                    df_dict[variant] += count["count"]
                except KeyError:
                    df_dict[variant] = count["count"]

        self.save_counts("synonymous", df_dict, raw=False)
        del df_dict

    def report_filtered_variant(self, variant, count):
        """
        Outputs a summary of the filtered variant to *handle*. Internal filter
        names are converted to messages using the ``SeqLib.filter_messages``
        dictionary. Related to :py:meth:`SeqLib.report_filtered`.
        """
        self.logger.debug(
            "Filtered variant (quantity={n}) (excess mutations)"
            "\n{read!s}".format(n=count, read=variant)
        )
