from __future__ import print_function
import Tkinter as tk
import ttk
import tkSimpleDialog
import tkMessageBox
import logging

logger = logging.getLogger(__name__)


class RunnerSavePrompt(tkSimpleDialog.Dialog):
    """
    Dialog box for prompting the user to save before running.
    """

    def __init__(self, parent_window, title="Enrich2"):
        self.pw = parent_window

        self.dialog_text = tk.StringVar()
        self.dialog_text.set("Would you like to save your config changes?")

        tkSimpleDialog.Dialog.__init__(self, parent_window, title)

    def body(self, master):
        frame = ttk.Frame(master, padding=(12, 6, 12, 6))
        frame.pack()

        dialog_text_label = ttk.Label(frame, textvariable=self.dialog_text)
        dialog_text_label.grid(column=0, row=0, sticky="nsew")

    def apply(self):
        self.pw.menu_save()


class RunnerWindow(tkSimpleDialog.Dialog):
    """
    Dialog box for blocking input while running the analysis.
    """

    def __init__(self, parent_window, title="Enrich2"):
        self.pw = parent_window
        self.run_button = None

        self.dialog_text = tk.StringVar()
        self.dialog_text.set("Ready to start analysis...")

        tkSimpleDialog.Dialog.__init__(self, parent_window, title)

    def body(self, master):
        frame = ttk.Frame(master, padding=(12, 6, 12, 6))
        frame.pack()

        dialog_text_label = ttk.Label(frame, textvariable=self.dialog_text)
        dialog_text_label.grid(column=0, row=0, sticky="nsew")

        self.run_button = tk.Button(
            frame, text="Begin", width=10, command=self.runner, default="active"
        )
        self.run_button.grid(column=0, row=1, sticky="nsew")

    def buttonbox(self):
        """
        Display no buttons.
        """
        pass

    def runner(self):
        # gray out the "Run" button
        self.run_button.config(state="disabled")
        self.update_idletasks()

        # set the analysis options
        self.pw.root_element.force_recalculate = self.pw.force_recalculate.get()
        self.pw.root_element.component_outliers = self.pw.component_outliers.get()
        self.pw.root_element.scoring_method = self.pw.scoring_method.get()
        self.pw.root_element.logr_method = self.pw.logr_method.get()
        self.pw.root_element.plots_requested = self.pw.plots_requested.get()
        self.pw.root_element.tsv_requested = self.pw.tsv_requested.get()

        # run the analysis, catching any errors to display in a dialog box
        try:
            # ensure that all objects are valid
            self.pw.root_element.validate()

            # open HDF5 files for the root and all child objects
            self.pw.root_element.store_open(children=True)

            # perform the analysis
            self.pw.root_element.calculate()

        except Exception as e:
            # display error
            logger.error(e)
            tkMessageBox.showerror(
                "Enrich2 Error", "Enrich2 encountered an error:\n{}".format(e)
            )

        else:
            # no exception occurred during calculation and setup
            # generate desired output
            if self.pw.plots_requested.get():
                try:
                    self.pw.root_element.make_plots()
                except Exception as e:
                    tkMessageBox.showwarning(
                        None,
                        "Calculations completed, but plotting failed:\n{}".format(e),
                    )
            if self.pw.tsv_requested.get():
                try:
                    self.pw.root_element.write_tsv()
                except Exception as e:
                    tkMessageBox.showwarning(
                        None,
                        "Calculations completed, but tsv output failed:\n{}".format(e),
                    )

            # show the dialog box
            tkMessageBox.showinfo("", "Analysis completed.")

        finally:
            # close the HDF5 files
            self.pw.root_element.store_close(children=True)

            # close this window
            self.destroy()
