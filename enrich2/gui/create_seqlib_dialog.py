from __future__ import print_function
import Tkinter as tk
import ttk
import tkSimpleDialog
from collections import OrderedDict


seqlib_label_text = OrderedDict([("BcvSeqLib", "Barcoded Variant"),
                                 ("BcidSeqLib", "Barcoded Identifier"),
                                 ("OverlapSeqLib", "Overlap"),
                                 ("BasicSeqLib", "Basic"),
                                 ("BarcodeSeqLib", "Barcodes Only"),
                                 ("IdOnlySeqLib", "Identifiers Only"),
                                 ])


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
            rb = ttk.Radiobutton(master, text=seqlib_label_text[k],
                                 variable=self.element_tkstring, value=k)
            rb.grid(column=0, row=(i + 1), sticky="w")
            if i == 0:
                rb.invoke()

    def buttonbox(self):
        """
        Display only one button.
        """
        box = tk.Frame(self)

        w = tk.Button(box, text="OK", width=10, command=self.ok, 
                      default="active")
        w.pack(side="left", padx=5, pady=5)

        self.bind("<Return>", self.ok)

        box.pack()

    def apply(self):
        try:
            self.element_type = globals()[self.element_tkstring.get()]
        except KeyError:
            raise KeyError("Unrecognized element type.")
