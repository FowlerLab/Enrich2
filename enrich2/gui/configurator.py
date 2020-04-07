from __future__ import print_function
import Tkinter as tk
import ttk
import tkFileDialog
import tkMessageBox
import platform
import json
from ..config_check import is_seqlib, is_experiment, is_selection, seqlib_type
from .delete_dialog import DeleteDialog
from .seqlib_apply_dialog import SeqLibApplyDialog
from .edit_dialog import EditDialog, clear_nones
from .create_root_dialog import CreateRootDialog
from .create_seqlib_dialog import CreateSeqLibDialog
from .runner_window import RunnerWindow, RunnerSavePrompt
from ..experiment import Experiment
from ..condition import Condition
from ..selection import Selection
from ..barcode import BarcodeSeqLib
from ..barcodevariant import BcvSeqLib
from ..barcodeid import BcidSeqLib
from ..basic import BasicSeqLib
from ..idonly import IdOnlySeqLib
from ..overlap import OverlapSeqLib
from ..seqlib import SeqLib
from ..storemanager import SCORING_METHODS, LOGR_METHODS


#: map class names to class definitions to avoid use of globals()
SEQLIB_CLASSES = {
    "BarcodeSeqLib": BarcodeSeqLib,
    "BcvSeqLib": BcvSeqLib,
    "BcidSeqLib": BcidSeqLib,
    "BasicSeqLib": BasicSeqLib,
    "IdOnlySeqLib": IdOnlySeqLib,
    "OverlapSeqLib": OverlapSeqLib,
}


def write_json(d, handle):
    """
    Write the contents of dictionary *d* to an open file *handle* in json format.

    *d* is passed to ``clear_nones`` before output.
    """
    json.dump(clear_nones(d), handle, indent=2, sort_keys=True)


class Configurator(tk.Tk):
    def __init__(self, version):
        tk.Tk.__init__(self)
        self.title("Enrich 2 Configurator v{}".format(version))

        self.treeview = None
        self.treeview_popup_target_id = None

        self.cfg_file_name = tk.StringVar()

        self.element_dict = dict()
        self.root_element = None

        # analysis options
        self.scoring_method = tk.StringVar()
        self.logr_method = tk.StringVar()
        self.force_recalculate = tk.BooleanVar()
        self.component_outliers = tk.BooleanVar()
        self.plots_requested = tk.BooleanVar()
        self.tsv_requested = tk.BooleanVar()

        # allow resizing
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)

        # create UI elements
        self.create_main_frame()
        self.create_menubar()
        self.create_treeview_context_menu()

    def create_treeview_context_menu(self):
        self.treeview_popup = tk.Menu(self, tearoff=0)
        self.treeview_popup.add_command(
            label="Apply FASTQ...", command=self.apply_seqlib_fastq
        )

    def apply_seqlib_fastq(self):
        SeqLibApplyDialog(self, self, self.treeview_popup_target_id)

    def treeview_context_menu(self, click):
        target = self.treeview.identify_row(click.y)
        if target != "":
            self.treeview_popup_target_id = target
            self.treeview_popup.post(click.x_root, click.y_root)
        self.treeview_popup_target_id = None

    def create_main_frame(self):
        # Frame for the Treeview and New/Edit/Delete buttons
        main = ttk.Frame(self, padding=(3, 3, 12, 12))
        main.rowconfigure(0, weight=1)
        main.columnconfigure(0, weight=1)
        main.grid(row=0, column=0, sticky="nsew")

        # Frame for the Treeview and its scrollbars
        tree_frame = ttk.Frame(main, padding=(3, 3, 12, 12))
        tree_frame.rowconfigure(0, weight=1)
        tree_frame.rowconfigure(1, weight=0)
        tree_frame.columnconfigure(0, weight=1)
        tree_frame.columnconfigure(1, weight=0)
        tree_frame.grid(row=0, column=0, sticky="nsew")

        # Treeview with column headings
        self.treeview = ttk.Treeview(tree_frame)
        self.treeview["columns"] = ("class", "barcodes", "variants")
        self.treeview.column("class", width=120)
        self.treeview.heading("class", text="Type")
        self.treeview.column("barcodes", width=25, stretch=tk.NO, anchor=tk.CENTER)
        self.treeview.heading("barcodes", text="BC")
        self.treeview.column("variants", width=25, stretch=tk.NO, anchor=tk.CENTER)
        self.treeview.heading("variants", text="V")
        self.treeview.grid(row=0, column=0, sticky="nsew")

        # Treeview context menu bindings
        self.treeview.bind("<Button-2>", self.treeview_context_menu)

        # Treeview scrollbars
        tree_ysb = tk.Scrollbar(
            tree_frame, orient="vertical", command=self.treeview.yview
        )
        tree_xsb = tk.Scrollbar(
            tree_frame, orient="horizontal", command=self.treeview.xview
        )
        tree_ysb.grid(row=0, column=1, sticky="nsw")
        tree_xsb.grid(row=1, column=0, sticky="ewn")
        self.treeview.config(yscroll=tree_ysb.set, xscroll=tree_xsb.set)

        # Frame for New/Edit/Delete buttons
        button_frame = ttk.Frame(main, padding=(3, 3, 12, 12))
        button_frame.grid(row=1, column=0)
        new_button = ttk.Button(
            button_frame, text="New...", command=self.new_button_press
        )
        new_button.grid(row=0, column=0)
        edit_button = ttk.Button(
            button_frame, text="Edit...", command=self.edit_button_press
        )
        edit_button.grid(row=0, column=1)
        delete_button = ttk.Button(
            button_frame, text="Delete", command=self.delete_button_press
        )
        delete_button.grid(row=0, column=2)

        # Frame for Analysis Options
        options_frame = ttk.Frame(main, padding=(3, 3, 12, 12))
        options_frame.grid(row=0, column=1, rowspan=2, sticky="nsew")

        row = 0
        heading = ttk.Label(options_frame, text="Analysis Options")
        heading.grid(column=0, row=row)
        row += 1

        scoring_heading = ttk.Label(options_frame, text="Scoring Method")
        scoring_heading.grid(column=0, row=row)
        row += 1
        for i, k in enumerate(SCORING_METHODS.keys()):
            rb = ttk.Radiobutton(
                options_frame,
                text=SCORING_METHODS[k].title(),
                variable=self.scoring_method,
                value=k,
            )
            rb.grid(column=0, row=row, sticky="w")
            row += 1
            if i == 0:
                rb.invoke()

        logr_heading = ttk.Label(options_frame, text="Normalization Method")
        logr_heading.grid(column=0, row=row)
        row += 1
        for i, k in enumerate(LOGR_METHODS.keys()):
            rb = ttk.Radiobutton(
                options_frame,
                text=LOGR_METHODS[k].title(),
                variable=self.logr_method,
                value=k,
            )
            rb.grid(column=0, row=row, sticky="w")
            row += 1
            if i == 0:
                rb.invoke()

        other_heading = ttk.Label(options_frame, text="Other Options")
        other_heading.grid(column=0, row=row)
        row += 1

        # force recalculate
        force_recalculate = ttk.Checkbutton(
            options_frame, text="Force Recalculation", variable=self.force_recalculate
        )
        force_recalculate.grid(column=0, row=row, sticky="w")
        row += 1

        # component outliers
        component_outliers = ttk.Checkbutton(
            options_frame,
            text="Component Outlier Statistics",
            variable=self.component_outliers,
        )
        component_outliers.grid(column=0, row=row, sticky="w")
        row += 1

        # make plots
        plots_requested = ttk.Checkbutton(
            options_frame, text="Make Plots", variable=self.plots_requested
        )
        plots_requested.grid(column=0, row=row, sticky="w")
        plots_requested.invoke()
        row += 1

        # write tsv
        tsv_requested = ttk.Checkbutton(
            options_frame, text="Write TSV Files", variable=self.tsv_requested
        )
        tsv_requested.grid(column=0, row=row, sticky="w")
        tsv_requested.invoke()
        row += 1

        go_button = ttk.Button(
            options_frame, text="Run Analysis", command=self.go_button_press
        )
        go_button.grid(column=0, row=row, sticky="sew")

    def go_button_press(self):
        if self.root_element is None:
            tkMessageBox.showwarning("", "No experimental design specified.")
        else:
            RunnerSavePrompt(self)
            RunnerWindow(self)

    def create_new_element(self):
        """
        Create and return a new element based on the current selection.

        This element is not added to the treeview. 
        """
        element = None
        parent_element = self.get_focused_element()
        if isinstance(parent_element, Experiment):
            element = Condition()
            element.parent = parent_element
        elif isinstance(parent_element, Condition):
            element = Selection()
            element.parent = parent_element
        elif isinstance(parent_element, Selection):
            element = CreateSeqLibDialog(self).element_type()
            element.parent = parent_element
        elif isinstance(parent_element, SeqLib):
            # special case: creates a copy of the selected SeqLib as a sibling
            element = type(parent_element)()
            element.configure(parent_element.serialize())
            element.parent = parent_element.parent
            # clear out the seqlib-specific values
            element.name = None
            element.timepoint = None
            element.counts_file = None
            if isinstance(element, OverlapSeqLib):
                element.forward = None
                element.reverse = None
            else:
                element.reads = None
        else:
            raise ValueError(
                "Unrecognized parent object type '{}'".format(type(parent_element))
            )
        return element

    def new_button_press(self):
        if self.treeview.focus() == "" and self.root_element is not None:
            tkMessageBox.showwarning(None, "No parent element selected.")
        else:
            if self.treeview.focus() == "" and self.root_element is None:
                element = CreateRootDialog(self).element
                if isinstance(element, SeqLib):
                    EditDialog(self, self, element)
                self.root_element = element
            else:
                element = self.create_new_element()
                EditDialog(self, self, element)

            # refresh the treeview and re-assign treeview id's
            self.refresh_treeview()

            # select the newly added element if it was successfully added
            if element.treeview_id in self.element_dict.keys():
                self.treeview.focus(element.treeview_id)
                self.treeview.selection_set(element.treeview_id)
            else:
                if element.parent is not None:
                    self.treeview.focus(element.parent.treeview_id)
                    self.treeview.selection_set(element.parent.treeview_id)
                del element

    def edit_button_press(self):
        if self.treeview.focus() == "":
            tkMessageBox.showwarning(None, "No element selected.")
        else:
            EditDialog(self, self, self.get_focused_element())

    def delete_button_press(self):
        if self.treeview.focus() == "":
            tkMessageBox.showwarning(None, "No element selected.")
        else:
            DeleteDialog(self, self)

    def create_menubar(self):
        # make platform-specific keybinds
        if platform.system() == "Darwin":
            accel_string = "Command+"
            accel_bind = "Command-"
        else:
            accel_string = "Ctrl+"
            accel_bind = "Control-"

        # create the menubar
        menubar = tk.Menu(self)

        # file menu
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(
            label="Open...",
            accelerator="{}O".format(accel_string),
            command=self.menu_open,
        )
        filemenu.add_command(
            label="Save", accelerator="{}S".format(accel_string), command=self.menu_save
        )
        filemenu.add_command(
            label="Save As...",
            accelerator="{}Shift+S".format(accel_string),
            command=self.menu_saveas,
        )
        menubar.add_cascade(label="File", menu=filemenu)

        # edit menu
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(
            label="Select All",
            accelerator="{}A".format(accel_string),
            command=self.menu_selectall,
        )
        menubar.add_cascade(label="Edit", menu=filemenu)

        # add the menubar
        self.config(menu=menubar)

        # add file menu keybinds
        self.bind("<{}o>".format(accel_bind), lambda event: self.menu_open())
        self.bind("<{}s>".format(accel_bind), lambda event: self.menu_save())
        self.bind("<{}Shift-s>".format(accel_bind), lambda event: self.menu_saveas())

        # add edit menu keybinds
        self.bind("<{}a>".format(accel_bind), lambda event: self.menu_selectall())

    def menu_open(self):
        fname = tkFileDialog.askopenfilename()
        if len(fname) > 0:  # file was selected
            try:
                with open(fname, "rU") as handle:
                    cfg = json.load(handle)
            except ValueError:
                tkMessageBox.showerror(None, "Failed to parse config file.")
            except IOError:
                tkMessageBox.showerror(None, "Could not read config file.")
            else:
                if is_experiment(cfg):
                    obj = Experiment()
                elif is_selection(cfg):
                    obj = Selection()
                elif is_seqlib(cfg):
                    obj = SEQLIB_CLASSES[seqlib_type(cfg)]()
                else:
                    tkMessageBox.showerror(None, "Unrecognized config format.")
                    return
                obj.output_dir_override = False
                try:
                    obj.configure(cfg)
                except Exception as e:
                    tkMessageBox.showerror(
                        None, "Failed to process config file:\n{}".format(e)
                    )
                else:
                    self.root_element = obj
                    self.cfg_file_name.set(fname)
                    self.refresh_treeview()

    def menu_save(self):
        if len(self.cfg_file_name.get()) == 0:
            self.menu_saveas()
        elif self.root_element is None:
            tkMessageBox.showwarning(None, "Cannot save empty configuration.")
        else:
            try:
                with open(self.cfg_file_name.get(), "w") as handle:
                    write_json(self.root_element.serialize(), handle)
            except IOError:
                tkMessageBox.showerror(None, "Failed to save config file.")
            else:
                tkMessageBox.showinfo(
                    None, "Save successful:\n{}".format(self.cfg_file_name.get())
                )

    def menu_saveas(self):
        if self.root_element is None:
            tkMessageBox.showwarning(None, "Cannot save empty configuration.")
        else:
            fname = tkFileDialog.asksaveasfilename()
            if len(fname) > 0:  # file was selected
                try:
                    with open(fname, "w") as handle:
                        write_json(self.root_element.serialize(), handle)
                except IOError:
                    tkMessageBox.showerror(None, "Failed to save config file.")
                else:
                    self.cfg_file_name.set(fname)
                    tkMessageBox.showinfo(
                        None, "Save successful:\n{}".format(self.cfg_file_name.get())
                    )

    def menu_selectall(self):
        """
        Add all elements in the Treeview to the selection.
        """
        for k in self.element_dict.keys():
            self.treeview.selection_add(k)

    def delete_element(self, tree_id):
        """
        Delete element with Treeview id *tree_id* from the tree, from the element 
        dictionary, and from the associated data structure. Recursively 
        deletes all children of *tree_id*.

        The tree should be refreshed using :py:meth:`refresh_tree` after 
        each deletion. This is the responsibility of the caller.

        """
        # if self.treeview.exists(tree_id):
        if tree_id in self.element_dict:
            # recursively delete children
            if self.element_dict[tree_id].children is not None:
                for child in self.element_dict[tree_id].children:
                    self.delete_element(child.treeview_id)

            # check if deleting the root element
            if self.root_element.treeview_id == tree_id:
                # clear the root element
                self.root_element = None
            else:
                try:
                    # remove the element from its parent's list of children
                    self.element_dict[tree_id].parent.remove_child_id(tree_id)
                except AttributeError:
                    raise AttributeError("Non-root element lacks proper parent")

            # delete the element from the dictionary
            del self.element_dict[tree_id]

    def refresh_treeview(self):
        """
        Clears the Treeview and repopulates it with the current contents of the tree.
        """
        # clear the entries in the Treeview
        for x in self.treeview.get_children():
            self.treeview.delete(x)

        # clear the id-element dictionary
        # elements may be given new id's after repopulation
        self.element_dict.clear()

        # repopulate
        if self.root_element is not None:
            self.populate_tree(self.root_element)

    def set_treeview_properties(self, element):
        """
        Set the information text for the Treeview *element*.
        """
        # set class property
        self.treeview.set(element.treeview_id, "class", element.treeview_class_name)

        # add the check marks for barcodes/variants
        if "variants" in element.labels:
            self.treeview.set(element.treeview_id, "variants", u"\u2713")
        else:
            self.treeview.set(element.treeview_id, "variants", "")
        if "barcodes" in element.labels:
            self.treeview.set(element.treeview_id, "barcodes", u"\u2713")
        else:
            self.treeview.set(element.treeview_id, "barcodes", "")

        self.treeview.set(element.treeview_id, "class", element.treeview_class_name)

    def populate_tree(self, element, parent_id=""):
        """
        Recursively populate the Treeview.

        Also populates the *id_cfgstrings*.
        """
        # insert into the Treeview
        element.treeview_id = self.treeview.insert(
            parent_id, "end", text=element.name, open=True
        )
        # add id-element pair to dictionary
        self.element_dict[element.treeview_id] = element
        # set information fields
        self.set_treeview_properties(element)

        # populate for children
        if element.children is not None:
            for child in element.children:
                self.populate_tree(child, parent_id=element.treeview_id)

    def get_element(self, treeview_id):
        """
        Returns the element with *treeview_id*.
        """
        return self.element_dict[treeview_id]

    def get_focused_element(self):
        """
        Returns the element that is currently being focused in the Treeview.

        Returns ``None`` if nothing is focused.
        """
        if self.treeview.focus() != "":
            return self.get_element(self.treeview.focus())
        else:
            return None

    def get_selected_elements(self):
        """
        Returns a list of elements that are currently selected in the Treeview.

        If no elements are selected, it returns an empty list.
        """
        return [self.get_element(x) for x in self.treeview.selection()]
