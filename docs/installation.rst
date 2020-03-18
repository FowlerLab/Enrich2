Getting started
=======================================================

.. _required packages:

Required packages
-------------------------------------------------------

Enrich2 runs on Python 2.7 and has the following dependencies:

* `NumPy <http://www.numpy.org/>`_ version 1.10.4 or higher
* `SciPy <http://www.scipy.org/>`_ version 0.16.0 or higher
* `pandas <http://pandas.pydata.org/>`_ version 0.18 or 0.19
* `PyTables <http://www.pytables.org/>`_ version 3.2.0 or higher
* `Statsmodels <http://statsmodels.sourceforge.net/>`_ version 0.6.1 or higher
* `matplotlib <http://matplotlib.org/>`_ version 1.4.3 or higher

The configuration GUI requires `Tkinter <https://docs.python.org/2/library/tkinter.html>`_. Building a local copy of the documentation requires `Sphinx <http://sphinx-doc.org/>`_.

.. note:: We recommend using a scientific Python distribution such as `Anaconda <https://store.continuum.io/cshop/anaconda/>`_ to install and manage dependencies.

.. note:: PyTables may not be installed when using the default settings for your distribution. If you encounter errors, check that the ``tables`` module is present. 

Installation and example dataset
-------------------------------------------------------

1. Make sure the `required packages`_ are installed.

    a. To set up a new `Anaconda <https://store.continuum.io/cshop/anaconda/>`_ environment for Enrich2, use the following environment file (:download:`click to download <_static/enrich2_env.yml>`).
Detailed instructions for setting up conda environments can be found in the `conda documentation <https://conda.io/docs/using/envs.html>`_.
    
        .. literalinclude:: _static/enrich2_env.yml

#. `Download Enrich2 <https://github.com/FowlerLab/Enrich2/archive/master.zip>`_ from the `GitHub repository <https://github.com/FowlerLab/Enrich2/>`_ and unzip it.

#. Using the terminal, navigate to the Enrich2 directory and run the setup script by typing ``python setup.py install``

To download the example dataset, visit the `Enrich2-Example GitHub repository <https://github.com/FowlerLab/Enrich2-Example/>`_.
Running this preconfigured analysis will create several :ref:`plots`. The :ref:`example-notebooks` demonstrate how to explore the :ref:`hdf5-files`.

Enrich2 executables
-------------------------------------------------------

The Enrich2 installer places two executable scripts into the user's path. Both executables run the same analysis, but through different interfaces.

* ``enrich_gui`` launches the Enrich2 graphical user interface. This is the recommended way to create a configuration file for Enrich2. See :ref:`gui-documentation` for a step-by-step guide.

* ``enrich_cmd`` launches the program from the command line. This is recommended for users performing analyses on a remote server who have already created configuration files. For a detailed list of command line options, type ``enrich_cmd --help``


