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
from __future__ import absolute_import
import logging
import re
from .constants import CODON_TABLE
from six.moves import range


class WildTypeSequence(object):
    """
    Container class for wild type sequence information. Used by :py:class:`~seqlib.seqlib.VariantSeqLib` objects and 
    :py:class:`~enrich2.selection.Selection` or :py:class:`~enrich2.experiment.Experiment` objects that contain 
    variant information.

    Requires a *parent_name* that associates this object with a StoreManager object for the 
    purposes of error reporting and logging.
    """
    def __init__(self, parent_name):
        self.parent_name = parent_name
        self.dna_seq = None
        self.protein_seq = None
        self.dna_offset = None
        self.protein_offset = None


    def __eq__(self, other):
        # note we don't need to check protein_offset, since it depends on dna_offset and protein_seq
        return self.dna_seq == other.dna_seq and self.protein_seq == other.protein_seq and self.dna_offset == other.dna_offset


    def __ne__(self, other):
        return not self == other


    def configure(self, cfg):
        try:
            # remove whitespace from WT DNA sequence and capitalize
            self.dna_seq = "".join(cfg['sequence'].split()).upper()

            # check that only valid characters are included (ACGT)
            if not re.match("^[ACGT]+$", self.dna_seq):
                raise ValueError("WT DNA sequence contains unexpected "
                                  "characters [{}]".format(self.parent_name))

            # set the reference offset
            if 'reference offset' in cfg:
                try:
                    self.dna_offset = int(cfg['reference offset'])
                except ValueError:
                    raise ValueError("Invalid reference offset value [{}]".format(self.parent_name))
            else:
                self.dna_offset = 0

            # handle coding sequences
            if cfg['coding']:
                # require coding sequences are in-frame
                if len(self.dna_seq) % 3 != 0:
                    raise ValueError("WT DNA sequence contains incomplete codons [{}]".format(self.parent_name))

                # perform translation
                self.protein_seq = ""
                for i in range(0, len(self.dna_seq), 3):
                    self.protein_seq += CODON_TABLE[self.dna_seq[i:i + 3]]

                # set the reference offset if it's a multiple of three
                if self.dna_offset % 3 == 0:
                    self.protein_offset = self.dna_offset // 3
                else:
                    logging.warning("Ignoring reference offset for protein changes (not a multiple of three)", extra={'oname' : self.parent_name})
                    self.protein_offset = 0
            else:
                self.protein_seq = None
                self.protein_offset = None

        except KeyError as key:
            raise KeyError("Missing required config value {key} [{name}]".format(key=key, name=self.parent_name))


    def serialize(self):
        """
        Format this object as a config object suitable for dumping to a config file.
        """
        cfg = {
                'sequence' : self.dna_seq, 
                'coding' : self.is_coding(), 
                'reference offset' : self.dna_offset
              }
        return cfg


    def is_coding(self):
        return self.protein_seq is not None


    def duplicate(self, new_parent_name):
        """
        Create a copy of this object with the *new_parent_name*.

        Uses the configure and serialize methods to perform the copy.
        """
        new = WildTypeSequence(new_parent_name)
        new.configure(self.serialize())

        if new != self:
            raise ValueError("Failed to duplicate wild type sequence [{}]".format(self.parent_name))
        else:
            return new
    

    def position_tuples(self, protein=False):
        """
        Return a list of tuples containing the position number (after offset adjustment) and 
        single-letter symbol (nucleotide or amino acid) for each position the wild type sequence.
        """
        if protein:
            if not self.is_coding():
                raise AttributeError("Cannot return wild type protein position tuples for non-coding wild type [{}]".format(self.parent_name))
            else:
                seq = self.protein_seq
                offset = self.protein_offset
        else:
            seq = self.dna_seq
            offset = self.dna_offset

        return [(i + offset + 1, seq[i]) for i in range(len(seq))]

