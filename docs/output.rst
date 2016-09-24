.. _hdf5-files:

Output HDF5 files
=======================================

Enrich2 stores data in an HDF5 file for each Experiment, Selection, and SeqLib analysis object. The name of the HDF5 file is the object's name plus the suffix "_<obj>.h5", where <obj> is the object type ("exp", "sel", or "lib").
Each file has multiple tables that can be queried and retrieved as pandas data frames (see :ref:`example-notebooks`). 

Each Experiment, Selection, and SeqLib has its own directory inside "Results/tsv/" containing tab-separated value files for users who want to work with other tools, such as R or Excel.

.. _output-table-organization:

Table organization
---------------------------------------------------

HDF5 files organize tables into groups like directories in a file system. Enrich2 has two top-level groups, "/main" (used for most tables) and "/raw" (used exclusively in SeqLibs to store raw counts). The first subgroup is typically the element type (variant, barcode, etc.), followed by the kind of data (counts, scores, etc.).

.. note:: When the "Force recalculation" analysis option is chosen, the "/main" tables are deleted from all HDF5 files in this analysis, and regenerated based on the "/raw" count data.

Enrich2 uses NaN (Not a Number) values to represent missing data, such as zero counts or scores that could not be calculated.

List of tables by object type
-------------------------------------------------------

Experiment
+++++++++++++++++++++++

Most experiment tables use a pandas MultiIndex for column names. The MultiIndex levels are: condition, selection (if applicable), and data value. See the `pandas advanced indexing documentation <http://pandas.pydata.org/pandas-docs/stable/advanced.html>`_ for more information on how to work with these objects.

* "/main/<element>/counts"

    Counts of elements that appear in at least one time point in the experiment.

* "/main/<element>/scores"

    Condition-level scores, standard errors, and epsilon (change in the standard error after the last iteration of the random-effects model) for all elements scored in all selections of at least one condition.

* "/main/<element>/scores_shared"

    Selection-level scores and standard errors for each element with at least one condition-level score.

* "/main/<element>/scores_shared_full"

    Selection-level scores and standard errors for each element scored in at least one selection.

* "/main/<element>/scores_pvalues_wt"

    z-scores and p-values for each variant or synonymous element with a condition-level score. The null hypothesis is that the element's score is equal to wild type.

* "/main/barcodemap"

    Barcode-variant or barcode-identifier map for all barcodes that appear in the Experiment. Only present for Barcoded Variant or Barcoded Identifier SeqLibs. 

Selection
+++++++++++++++++++++++

* "/main/<element>/counts"

    Counts of elements that appear in all time points in the selection.

* "/main/<element>/counts_unfiltered"

    Counts of elements that appear in at least one time point in the selection.

* "/main/<element>/scores"

    Scores, standard errors, standard error percentiles, and method-specific values (e.g. regression slope and intercept) for all elements counted in all time points in the selection.

* "/main/<element>/weights"

    Regression weights for each element at each time point in weighted least squares regression.

* "/main/<element>/log_ratios"

    Y-values for each element at each time point in weighted and ordinary least squares regression. 

* "/main/barcodemap"

    Barcode-variant or barcode-identifier map for all barcodes that appear in the Selection. Only present for Barcoded Variant or Barcoded Identifier SeqLibs. 

SeqLib
+++++++++++++++++++++++

* "/main/<element>/counts"

    Counts of elements after minimum count filtering and barcode mapping.

* "/raw/<element>/counts"

    Counts of elements taken directly from the FASTQ_ data.

* "/raw/filter"

    Number of reads removed for each FASTQ_ filtering option.

* "/raw/barcodemap"

    Barcode-variant or barcode-identifier map for barcodes that appear in this SeqLib. Only present for Barcoded Variant or Barcoded Identifier SeqLibs. 
