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
import re
import os.path

DEFAULT_COLUMNS = 3

class SectionLabel(object):
    def __init__(self, text):
        self.text = text


    def body(self, master, row, columns=DEFAULT_COLUMNS, **kwargs):
        label = tkinter.ttk.Label(master, text=self.text)
        label.grid(row=row, column=0, columnspan=columns, sticky="w")
        return 1


    def validate(self):
        return True


    def apply(self):
        return None



class Checkbox(object):
    def __init__(self, text, cfg, key):
        self.value = tk.BooleanVar()
        self.text = text
        self.cfg = cfg
        self.key = key
        try:
            if self.cfg[self.key] not in (True, False):
                self.value.set(False)
            else:
                self.value.set(self.cfg[self.key])
        except KeyError:
            self.value.set(False)   # default to False


    def body(self, master, row, columns=DEFAULT_COLUMNS, **kwargs):
        """
        Place the required elements using the grid layout method.

        Returns the number of rows taken by this element.
        """
        checkbox = tkinter.ttk.Checkbutton(master, text=self.text, variable=self.value)
        checkbox.grid(row=row, column=0, columnspan=columns, sticky="w")
        return 1


    def validate(self):
        return True


    def apply(self):
        self.cfg[self.key] = self.value.get()



class MyEntry(object):
    """
    Base class for labeled Entry fields.

    *text* is the Label/error box text.
    """
    def __init__(self, text, cfg, key, optional=False):
        self.value = tk.StringVar()
        self.text = text
        self.cfg = cfg
        self.key = key
        self.optional = optional
        try:
            if self.cfg[self.key] is None:
                self.value.set("")
            else:
                self.value.set(self.cfg[self.key])
        except KeyError:
            self.value.set("")


    def body(self, master, row, columns=DEFAULT_COLUMNS, **kwargs):
        """
        Place the required elements using the grid layout method.

        Returns the number of rows taken by this element.
        """
        label = tkinter.ttk.Label(master, text=self.text)
        label.grid(row=row, column=0, columnspan=1, sticky="e")
        entry = tkinter.ttk.Entry(master, textvariable=self.value)
        entry.grid(row=row, column=1, columnspan=columns - 1, sticky="ew")
        return 1


    def validate(self):
        """
        Validates the input. Returns ``True`` unless the field is blank and *optional* is ``False``.
        """
        if not self.optional and len(self.value.get()) == 0:
            tkinter.messagebox.showwarning("", "{} not specified.".format(self.text))
            return False
        else:
            return True


    def apply(self):
        if len(self.value.get()) > 0:
            self.cfg[self.key] = self.value.get()
        else:
            self.cfg[self.key] = None



class FileEntry(MyEntry):
    """
    Creates a labeled Entry field for a file or directory.

    *text* is the Label/error box text.
    *directory* is ``True`` if the Entry is for selecting a directory (instead of a file).
    *extensions* is a list of valid file endings

    """
    def __init__(self, text, cfg, key, optional=False, directory=False, extensions=None):
        MyEntry.__init__(self, text, cfg, key, optional)
        self.directory = directory
        if extensions is not None:
            self.extensions = [x.lower() for x in extensions]
        else:
            self.extensions = None


    def body(self, master, row, columns=DEFAULT_COLUMNS, **kwargs):
        """
        Place the required elements using the grid layout method.

        Returns the number of rows taken by this element.
        """
        label = tkinter.ttk.Label(master, text=self.text)
        label.grid(row=row, column=0, columnspan=1, sticky="e")
        entry = tkinter.ttk.Entry(master, textvariable=self.value)
        entry.grid(row=row, column=1, columnspan=columns - 1, sticky="ew")
        if self.directory:
            choose = tkinter.ttk.Button(master, text="Choose...", command=lambda: self.value.set(tkinter.filedialog.askdirectory()))
        else:
            choose = tkinter.ttk.Button(master, text="Choose...", command=lambda: self.value.set(tkinter.filedialog.askopenfilename()))
        choose.grid(row=row + 1, column=1, sticky="w")
        if self.optional:
            clear = tkinter.ttk.Button(master, text="Clear", command=lambda: self.value.set(""))
            clear.grid(row=row + 1, column=2, sticky="e")
        return 2


    def validate(self):
        if len(self.value.get()) == 0:
            if not self.optional:
                tkinter.messagebox.showwarning("", "{} not specified.".format(self.text))
                return False
            else:
                return True
        else:
            if os.path.exists(self.value.get()):
                if self.extensions is not None:
                    if any(self.value.get().lower().endswith(x) for x in self.extensions):
                        return True
                    else:
                        tkinter.messagebox.showwarning("", "Invalid file extension for {}.".format(self.text))
                        return False
                else: # no extension restriction
                    return True
            else:
                tkinter.messagebox.showwarning("", "{} file does not exist.".format(self.text))
                return False



class StringEntry(MyEntry):
    """
    Creates a labeled Entry field for a string.

    *text* is the Label/error box text.
    """
    def __init__(self, text, cfg, key, optional=False):
        MyEntry.__init__(self, text, cfg, key, optional)


    def body(self, master, row, columns=DEFAULT_COLUMNS, **kwargs):
        """
        Place the required elements using the grid layout method.

        Returns the number of rows taken by this element.
        """
        label = tkinter.ttk.Label(master, text=self.text)
        label.grid(row=row, column=0, columnspan=1, sticky="e")
        entry = tkinter.ttk.Entry(master, textvariable=self.value)
        entry.grid(row=row, column=1, columnspan=columns - 1, sticky="ew")
        return 1



class IntegerEntry(MyEntry):
    """
    Creates a labeled Entry field for an integer.

    *text* is the Label/error box text.
    """
    def __init__(self, text, cfg, key, optional=False, minvalue=0):
        MyEntry.__init__(self, text, cfg, key, optional)
        self.minvalue = minvalue


    def body(self, master, row, columns=DEFAULT_COLUMNS, width=4, left=False, **kwargs):
        """
        Add the labeled entry to the Frame *master* using grid at *row*.

        *width* controls the width of the Entry.
        *left* is ``True`` if the Entry is to the left of the Label.
        *columns* is the number of columns in *master*.

        Returns the number of rows taken by this element.
        """
        if left:
            entry_column = 0
            entry_sticky = "e"
            entry_width = 1
            label_column = 1
            label_sticky = "w"
            label_width = columns - 1
        else:
            entry_column = 1
            entry_sticky = "w"
            entry_width = columns - 1
            label_column = 0
            label_sticky = "e"
            label_width = 1

        label = tkinter.ttk.Label(master, text=self.text)
        label.grid(row=row, column=label_column, columnspan=label_width, sticky=label_sticky)
        entry = tkinter.ttk.Entry(master, textvariable=self.value, width=width)
        entry.grid(row=row, column=entry_column, columnspan=entry_width, sticky=entry_sticky)
        return 1


    def validate(self):
        """
        Returns ``True`` if the value entered validates; else ``False``.

        If *self.optional* is ``True``, the field can be empty.
        Checks the *self.minvalue* that was passed on creation.
        """
        try:
            intvalue = int(self.value.get())
        except ValueError:
            if len(self.value.get()) == 0:
                if not self.optional:
                    tkinter.messagebox.showwarning("", "{} not specified.".format(self.text))
                    return False
                else:
                    return True
            else:
                tkinter.messagebox.showwarning("", "{} is not an integer.".format(self.text))
                return False
        else:
            if intvalue < self.minvalue:
                tkinter.messagebox.showwarning("", "{} lower than minimum value ({}).".format(self.text, self.minvalue))
                return False
            else:
                return True


    def apply(self):
        if len(self.value.get()) > 0:
            self.cfg[self.key] = int(self.value.get())
        else:
            self.cfg[self.key] = None

