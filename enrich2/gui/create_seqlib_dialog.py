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
import Tkinter as tk
import ttk
import tkSimpleDialog
import tkMessageBox
import tkFileDialog
import json
from copy import deepcopy
from collections import OrderedDict
from ..basic import BasicSeqLib
from ..barcodevariant import BcvSeqLib
from ..barcodeid import BcidSeqLib
from ..barcode import BarcodeSeqLib
from ..overlap import OverlapSeqLib
from ..seqlib import SeqLib
from ..variant import VariantSeqLib


seqlib_label_text = OrderedDict([("BcvSeqLib", "Barcoded Variant"),
                                 ("BcidSeqLib", "Barcoded Identifier"),
                                 ("OverlapSeqLib", "Overlap"),
                                 ("BasicSeqLib", "Basic"), 
                                 ("BarcodeSeqLib", "Barcodes Only")])



class CreateSeqLibDialog(tkSimpleDialog.Dialog):
    """
    Dialog box for creating a new SeqLib.
    """
    def __init__(self, parent_window, title="New SeqLib"):
        self.element_tkstring = tk.StringVar()
        self.element_type = None
        tkSimpleDialog.Dialog.__init__(self, parent_window, title)


    def body(self, master):
        message = ttk.Label(master, text="SeqLib type:")
        message.grid(column=0, row=0)

        for i, k in enumerate(seqlib_label_text.keys()):
            rb = ttk.Radiobutton(master, text=seqlib_label_text[k], variable=self.element_tkstring, value=k)
            rb.grid(column=0, row=(i + 1), sticky="w")
            if i == 0:
                rb.invoke()


    def buttonbox(self):
        """
        Display only one button.
        """
        box = tk.Frame(self)

        w = tk.Button(box, text="OK", width=10, command=self.ok, default="active")
        w.pack(side="left", padx=5, pady=5)

        self.bind("<Return>", self.ok)

        box.pack()


    def apply(self):
        try:
            self.element_type = globals()[self.element_tkstring.get()]
        except KeyError:
            raise KeyError("Unrecognized element type.")
