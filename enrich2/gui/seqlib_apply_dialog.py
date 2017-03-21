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

class SeqLibApplyDialog(tkinter.simpledialog.Dialog):
    """
    Confirmation dialog box for applying FASTQ filtering options to selected SeqLibs from the Treeview.
    """
    def __init__(self, parent_window, tree, source_id, title="Confirm Filtering Changes"):
        self.tree = tree
        self.source_id = source_id
        self.target_ids = [x for x in self.tree.treeview.selection() if x != source_id and type(self.tree.get_element(self.source_id)) == type(self.tree.get_element(x))]
        tkinter.simpledialog.Dialog.__init__(self, parent_window, title)


    def body(self, master):
        """
        Generates the required text listing all SeqLibs that will have their FASTQ options updated.

        Displays the "OK" and "Cancel" buttons.
        """
        if len(self.target_ids) == 0:
            message_string = "No elegible SeqLibs selected."
        elif len(self.target_ids) == 1:
            message_string = 'Apply FASTQ filtering options from "{}" to "{}"?'.format(self.tree.get_element(self.source_id).name, self.tree.get_element(self.target_ids[0]).name)
        else:
            bullet = "    " + u"\u25C6"
            message_string = 'Apply FASTQ filtering options from "{}"" to the following?\n'.format(self.tree.get_element(self.source_id).name)
            for x in self.target_ids:
                message_string += u"{bullet} {name}\n".format(bullet=bullet, name=self.tree.get_element(x).name)
        message = tkinter.ttk.Label(master, text=message_string, justify="left")
        message.grid(row=0, sticky="w")


    def buttonbox(self):
        """
        Display only one button if there's no selection. Otherwise, use the default method to display two buttons.
        """
        if len(self.target_ids) == 0:
            box = tk.Frame(self)

            w = tk.Button(box, text="OK", width=10, command=self.cancel, default="active")
            w.pack(side="left", padx=5, pady=5)

            self.bind("<Return>", self.cancel)

            box.pack()
        else:
            tkinter.simpledialog.Dialog.buttonbox(self)


    def apply(self):
        """
        Called when the user chooses "OK". Performs the FASTQ filtering update.
        """
        filter_cfg = self.tree.get_element(self.source_id).serialize_filters()
        for x in self.target_ids:
            self.tree.get_element(x).filters = filter_cfg
        self.tree.refresh_treeview()
