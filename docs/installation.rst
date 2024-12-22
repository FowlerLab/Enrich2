Getting started
=======================================================

.. _required packages:

Required packages
-------------------------------------------------------

Enrich2 runs on Python 3 and has the following dependencies:

* `NumPy <http://www.numpy.org/>`_
* `SciPy <http://www.scipy.org/>`_
* `pandas <http://pandas.pydata.org/>`_
* `PyTables <http://www.pytables.org/>`_
* `Statsmodels <http://statsmodels.sourceforge.net/>`_
* `matplotlib <http://matplotlib.org/>`_
* `fqfa <https://fqfa.readthedocs.io/>`_

The configuration GUI requires `Tkinter <https://docs.python.org/2/library/tkinter.html>`_. Building a local copy of the documentation requires `Sphinx <http://sphinx-doc.org/>`_.

.. note:: PyTables may not be installed by your distribution. If you encounter errors, check that the ``tables`` module is present.

.. note:: Tkinter may not be installed by your distribution. If you encounter errors, try installing ``python3-tk`` or similar using your system package manager.

Installation and example dataset
-------------------------------------------------------

You can install Enrich2 in a new `virtual environment <https://docs.python.org/3/library/venv.html>`_ using `pip <https://docs.python.org/3/installing/index.html>`_::

    python3 -m venv e2env
    source e2env/bin/activate
    pip install enrich2

To download the example dataset, visit the `Enrich2-Example GitHub repository <https://github.com/FowlerLab/Enrich2-Example/>`_.
Running this preconfigured analysis will create several :ref:`plots`. The :ref:`example-notebooks` demonstrate how to explore the :ref:`hdf5-files`.

Enrich2 executables
-------------------------------------------------------

The Enrich2 installer places two executable scripts into the user's path. Both executables run the same analysis, but through different interfaces.

* ``enrich_gui`` launches the Enrich2 graphical user interface. This is the recommended way to create a configuration file for Enrich2. See :ref:`gui-documentation` for a step-by-step guide.

* ``enrich_cmd`` launches the program from the command line. This is recommended for users performing analyses on a remote server who have already created configuration files. For a detailed list of command line options, type ``enrich_cmd --help``


