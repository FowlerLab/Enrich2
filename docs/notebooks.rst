.. _example-notebooks:

Example notebooks
====================================

Begin exploring Enrich2 datasets with the following notebooks. They rely on the `Enrich2 example dataset <https://github.com/FowlerLab/Enrich2-Example/>`_, so please perform that analysis before running any of these notebooks locally.

The notebooks can be run interactively by using the command line to navigate to the "Enrich2/docs/notebooks" directory and enter ``jupyter notebook <notebook.ipynb>`` where ``<notebook.ipynb>`` is the notebook file name.

The first two notebooks demonstrate using pandas to open an HDF5 file, extract its contents into a data frame, and perform queries on tables in the HDF5 file. For more information, see the `pandas HDF5 documentation <http://pandas.pydata.org/pandas-docs/stable/io.html#hdf5-pytables>`_.

.. include:: exported_notebooks/min_count.rst

.. include:: exported_notebooks/unique_barcodes.rst

For more information on Enrich2 data tables, see :ref:`hdf5-files`.
