from __future__ import print_function
import Tkinter as tk
import ttk
import tkSimpleDialog
from .dialog_elements import FileEntry, StringEntry, DEFAULT_COLUMNS
from .create_seqlib_dialog import SEQLIB_LABEL_TEXT
from ..barcode import BarcodeSeqLib
from ..barcodevariant import BcvSeqLib
from ..barcodeid import BcidSeqLib
from ..basic import BasicSeqLib
from ..idonly import IdOnlySeqLib
from ..overlap import OverlapSeqLib
from ..selection import Selection
from ..experiment import Experiment


#: map class names to class definitions to avoid use of globals()
ELEMENT_CLASSES = {
    "BarcodeSeqLib": BarcodeSeqLib,
    "BcvSeqLib": BcvSeqLib,
    "BcidSeqLib": BcidSeqLib,
    "BasicSeqLib": BasicSeqLib,
    "IdOnlySeqLib": IdOnlySeqLib,
    "OverlapSeqLib": OverlapSeqLib,
    "Selection": Selection,
    "Experiment": Experiment,
}


class CreateRootDialog(tkSimpleDialog.Dialog):
    """
    Dialog box for creating a new root element.
    """

    def __init__(self, parent_window, title="Create Root Object"):
        self.element_tkstring = tk.StringVar()
        self.cfg_dict = dict()
        self.output_directory_tk = FileEntry(
            "Output Directory",
            self.cfg_dict,
            "output directory",
            optional=False,
            directory=True,
        )
        self.name_tk = StringEntry("Name", self.cfg_dict, "name", optional=False)
        self.element = None
        tkSimpleDialog.Dialog.__init__(self, parent_window, title)

    def body(self, master):
        row_no = self.name_tk.body(master, 0)
        row_no += self.output_directory_tk.body(master, row_no)

        element_types = ttk.Frame(master, padding=(3, 3, 12, 12))
        element_types.grid(
            column=0, row=row_no, sticky="nsew", columnspan=DEFAULT_COLUMNS
        )

        message = ttk.Label(element_types, text="Root object type:")
        message.grid(column=0, row=0)

        label = ttk.Label(element_types, text="Experiment")
        label.grid(column=0, row=1, sticky="w")
        rb = ttk.Radiobutton(
            element_types,
            text="Experiment",
            variable=self.element_tkstring,
            value="Experiment",
        )
        rb.grid(column=0, row=2, sticky="w")
        rb.invoke()

        label = ttk.Label(element_types, text="Selection")
        label.grid(column=0, row=3, sticky="w")
        rb = ttk.Radiobutton(
            element_types,
            text="Selection",
            variable=self.element_tkstring,
            value="Selection",
        )
        rb.grid(column=0, row=4, sticky="w")

        label = ttk.Label(element_types, text="SeqLib")
        label.grid(column=0, row=5, sticky="w")
        for i, k in enumerate(SEQLIB_LABEL_TEXT.keys()):
            rb = ttk.Radiobutton(
                element_types,
                text=SEQLIB_LABEL_TEXT[k],
                variable=self.element_tkstring,
                value=k,
            )
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
            self.element = ELEMENT_CLASSES[self.element_tkstring.get()]()
        except KeyError:
            raise KeyError(
                "Unrecognized element type '{}'".format(self.element_tkstring.get())
            )

        # set the properties from this dialog
        self.element.output_dir_override = False
        self.element.output_dir = self.cfg_dict["output directory"]
        self.element.name = self.cfg_dict["name"]
