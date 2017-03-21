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

import tkinter as tk
import tkinter.ttk
import tkinter.simpledialog
import tkinter.messagebox
import tkinter.filedialog

import json
from copy import deepcopy
from collections import OrderedDict
from .dialog_elements import FileEntry, StringEntry, DEFAULT_COLUMNS
from ..experiment import Experiment
from ..condition import Condition
from ..selection import Selection
from ..basic import BasicSeqLib
from ..barcodevariant import BcvSeqLib
from ..barcodeid import BcidSeqLib
from ..barcode import BarcodeSeqLib
from ..overlap import OverlapSeqLib
from ..seqlib import SeqLib
from ..variant import VariantSeqLib
from .create_seqlib_dialog import seqlib_label_text


class CreateRootDialog(tkinter.simpledialog.Dialog):
    """
    Dialog box for creating a new root element.
    """
    def __init__(self, parent_window, title="Create Root Object"):
        self.element_tkstring = tk.StringVar()
        self.cfg_dict = dict()
        self.output_directory_tk = FileEntry("Output Directory", self.cfg_dict, 'output directory', optional=False, directory=True)
        self.name_tk = StringEntry("Name", self.cfg_dict, 'name', optional=False)
        self.element = None
        tkinter.simpledialog.Dialog.__init__(self, parent_window, title)


    def body(self, master):
        row_no = self.name_tk.body(master, 0)
        row_no += self.output_directory_tk.body(master, row_no)

        element_types = tkinter.ttk.Frame(master, padding=(3, 3, 12, 12))
        element_types.grid(column=0, row=row_no, sticky="nsew", columnspan=DEFAULT_COLUMNS)

        message = tkinter.ttk.Label(element_types, text="Root object type:")
        message.grid(column=0, row=0)

        label = tkinter.ttk.Label(element_types, text="Experiment")
        label.grid(column=0, row=1, sticky="w")
        rb = tkinter.ttk.Radiobutton(element_types, text="Experiment", variable=self.element_tkstring, value="Experiment")
        rb.grid(column=0, row=2, sticky="w")
        rb.invoke()

        label = tkinter.ttk.Label(element_types, text="Selection")
        label.grid(column=0, row=3, sticky="w")
        rb = tkinter.ttk.Radiobutton(element_types, text="Selection", variable=self.element_tkstring, value="Selection")
        rb.grid(column=0, row=4, sticky="w")

        label = tkinter.ttk.Label(element_types, text="SeqLib")
        label.grid(column=0, row=5, sticky="w")
        for i, k in enumerate(seqlib_label_text.keys()):
            rb = tkinter.ttk.Radiobutton(element_types, text=seqlib_label_text[k], variable=self.element_tkstring, value=k)
            rb.grid(column=0, row=(i + 6), sticky="w")


    def buttonbox(self):
        """
        Display only one button.
        """
        box = tk.Frame(self)

        w = tk.Button(box, text="OK", width=10, command=self.ok, default="active")
        w.pack(side="left", padx=5, pady=5)

        self.bind("<Return>", self.ok)

        box.pack()


    def validate(self):
        # check the fields
        return self.output_directory_tk.validate() and self.name_tk.validate()


    def apply(self):
        # apply the fields
        self.output_directory_tk.apply()
        self.name_tk.apply()

        # create the object
        try:
            self.element = globals()[self.element_tkstring.get()]()
        except KeyError:
            raise KeyError("Unrecognized element type '{}'".format(self.element_tkstring.get()))

        # set the properties from this dialog
        self.element.output_dir_override = False
        self.element.output_dir = self.cfg_dict['output directory']
        self.element.name = self.cfg_dict['name']

