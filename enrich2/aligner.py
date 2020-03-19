"""
Module for alignment of variants to the wild type sequence.

This module is optional, and using it will dramatically increase runtime when
counting variants. It is only recommended for users who need to count
insertion and deletion variants (i.e. not coding sequences).
"""

import numpy as np

#: Default similarity matrix used by the aligner.
#: User-defined matrices must have this format.
_simple_similarity = {
    "A": {"A": 1, "C": -1, "G": -1, "T": -1, "N": 0, "X": 0},
    "C": {"A": -1, "C": 1, "G": -1, "T": -1, "N": 0, "X": 0},
    "G": {"A": -1, "C": -1, "G": 1, "T": -1, "N": 0, "X": 0},
    "T": {"A": -1, "C": -1, "G": -1, "T": 1, "N": 0, "X": 0},
    "N": {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "X": 0},
    "X": {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "X": 0},
    "gap": -1,
}


class Aligner(object):
    """
    Class for performing local alignment of two DNA sequences.

    This class implements `Needleman-Wunsch <http://en.wikipedia.org/wiki/
    Needleman%E2%80%93Wunsch_algorithm>`_ local alignment.

    The :py:class:`~aligner.Aligner` requires a scoring matrix when
    created. The format is a nested dictionary, with a special ``'gap'`` entry
    for the gap penalty (this value is used for both gap opening and gap
    extension).

    The ``'X'`` nucleotide is a special case for unresolvable mismatches in
    :py:class:`~overlap.OverlapSeqLib` variant data.
    """

    _MAT = 1  # match
    _INS = 2  # insertion (with respect to wild type)
    _DEL = 3  # deletion (with respect to wild type)
    _END = 4  # end of traceback

    def __init__(self, similarity=_simple_similarity):
        similarity_keys = similarity.keys()
        if "gap" in similarity_keys:
            similarity_keys.remove("gap")
        for key in similarity_keys:
            if not all(x in similarity[key] for x in similarity_keys) or len(
                similarity[key]
            ) != len(similarity_keys):
                raise ValueError("Asymmetrical alignment scoring matrix")

        self.similarity = similarity
        if "gap" not in self.similarity:
            raise ValueError("No gap penalty in alignment scoring matrix.")

        self.matrix = None
        self.seq1 = None
        self.seq2 = None
        self.calls = 0

    def align(self, seq1, seq2):
        """
        Aligns the two sequences, *seq1* and *seq2* and returns a list of
        tuples describing the differences between the sequences.

        The tuple format is ``(i, j, type, length)``, where ``i`` and ``j``
        are the positions in *seq1* and *seq2*, respectively, and type is one
        of ``"match"``, ``"mismatch"``, ``"insertion"``, or ``"deletion"``.
        For indels, the ``length`` value is the number of bases inserted or
        deleted with respect to *seq1* starting at ``i``.
        """
        self.matrix = np.ndarray(
            shape=(len(seq1) + 1, len(seq2) + 1),
            dtype=np.dtype([("score", np.int), ("trace", np.byte)]),
        )
        seq1 = seq1.upper()
        seq2 = seq2.upper()

        # build matrix of scores/traceback information
        for i in xrange(len(seq1) + 1):
            self.matrix[i, 0] = (self.similarity["gap"] * i, Aligner._DEL)
        for j in xrange(len(seq2) + 1):
            self.matrix[0, j] = (self.similarity["gap"] * j, Aligner._INS)
        for i in xrange(1, len(seq1) + 1):
            for j in xrange(1, len(seq2) + 1):
                match = (
                    self.matrix[i - 1, j - 1]["score"]
                    + self.similarity[seq1[i - 1]][seq2[j - 1]],
                    Aligner._MAT,
                )
                delete = (
                    self.matrix[i - 1, j]["score"] + self.similarity["gap"],
                    Aligner._DEL,
                )
                insert = (
                    self.matrix[i, j - 1]["score"] + self.similarity["gap"],
                    Aligner._INS,
                )
                self.matrix[i, j] = max(delete, insert, match, key=lambda x: x[0])
        self.matrix[0, 0] = (0, Aligner._END)

        # calculate alignment from the traceback
        i = len(seq1)
        j = len(seq2)
        traceback = list()
        while i > 0 or j > 0:
            if self.matrix[i, j]["trace"] == Aligner._MAT:
                if seq1[i - 1] == seq2[j - 1]:
                    traceback.append((i - 1, j - 1, "match", None))
                else:
                    traceback.append((i - 1, j - 1, "mismatch", None))
                i -= 1
                j -= 1
            elif self.matrix[i, j]["trace"] == Aligner._INS:
                traceback.append((i - 1, j - 1, "insertion", 1))
                j -= 1
            elif self.matrix[i, j]["trace"] == Aligner._DEL:
                traceback.append((i - 1, j - 1, "deletion", 1))
                i -= 1
            elif self.matrix[i, j]["trace"] == Aligner._END:
                pass
            else:
                raise RuntimeError("Invalid value in alignment traceback.")
        traceback.reverse()

        # combine indels
        indel = None
        traceback_combined = list()
        for t in traceback:
            if t[2] == "insertion" or t[2] == "deletion":
                if indel is not None:
                    if t[2] == indel[2]:
                        indel[3] += t[3]
                    else:
                        raise RuntimeError(
                            "Aligner failed to combine indels. " "Check gap penalty."
                        )
                else:
                    indel = list(t)
            else:
                if indel is not None:
                    traceback_combined.append(tuple(indel))
                    indel = None
                traceback_combined.append(t)
        if indel is not None:
            traceback_combined.append(tuple(indel))

        self.calls += 1
        return traceback_combined
