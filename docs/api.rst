Appendix: API documentation
#############################################

This page contains automatically generated documentation from the Enrich2 codebase. It is intended for developers and advanced users.

:py:mod:`~enrich2.storemanager` --- Abstract class for Enrich2 data
===================================================================

.. py:module:: storemanager
	:synopsis: Abstract class for Enrich2 data.

This module contains the class definition for the :py:class:`~enrich2.storemanager.storemanager.StoreManager` abstract class, the shared base class for most classes in the `Enrich2 <index.html>`_ project. This class provides general behavior for the GUI and for handling HDF5 data files.

:py:class:`~enrich2.storemanager.StoreManager` class
-----------------------------------------------------------------
.. autoclass:: enrich2.storemanager.StoreManager
	:members:


:py:mod:`~enrich2.seqlib` --- Sequencing library file handling and element counting
===================================================================================

.. py:module:: seqlib
	:synopsis: Sequencing library file handling and element counting.

This module provides class definitions for the various types of sequencing library designs usable by `Enrich2 <index.html>`_. Data for each FASTQ_ file (or pair of overlapping FASTQ_ files for overlapping paired-end data) is read into its own :py:class:`~enrich2.seqlib.SeqLib` object. If necessary, FASTQ_ files should be split by index read before being read by a :py:class:`~enrich2.seqlib.SeqLib` object. :py:class:`~enrich2.seqlib.SeqLib` objects are coordinated by :py:mod:`~enrich2.selection.Selection` objects.

:py:class:`~enrich2.seqlib.SeqLib` and :py:class:`~enrich2.variant.VariantSeqLib` are abstract classes. 

:py:class:`~enrich2.seqlib.SeqLib` class
-------------------------------------------------------
.. autoclass:: enrich2.seqlib.SeqLib
	:members:

:py:class:`~enrich2.variant.VariantSeqLib` class
-------------------------------------------------------
.. autoclass:: enrich2.variant.VariantSeqLib
	:members:

:py:class:`~enrich2.barcode.BarcodeSeqLib` class
-------------------------------------------------------
.. autoclass:: enrich2.barcode.BarcodeSeqLib
	:members:

:py:class:`~enrich2.barcodevariant.BcvSeqLib` class
-------------------------------------------------------------
.. autoclass:: enrich2.barcodevariant.BcvSeqLib
	:members:

:py:class:`~enrich2.barcodeid.BcidSeqLib` class
-------------------------------------------------------------
.. autoclass:: enrich2.barcodeid.BcidSeqLib
	:members:

:py:class:`~enrich2.basic.BasicSeqLib` class
-----------------------------------------------------
.. autoclass:: enrich2.basic.BasicSeqLib
	:members:

:py:class:`~enrich2.overlap.OverlapSeqLib` class
--------------------------------------------------------
.. autoclass:: enrich2.overlap.OverlapSeqLib
	:members:

:py:class:`~enrich2.idonly.IdOnlySeqLib` class
-------------------------------------------------------------
.. autoclass:: enrich2.idonly.IdOnlySeqLib
	:members:

:py:class:`~enrich2.seqlib.SeqLib` helper classes
-------------------------------------------------------

:py:class:`~enrich2.aligner.Aligner` class
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: enrich2.aligner.Aligner
	:members:

:py:class:`~enrich2.wildtype.WildTypeSequence` class
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: enrich2.wildtype.WildTypeSequence
	:members:

:py:class:`~enrich2.barcodemap.BarcodeMap` class
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: enrich2.barcodemap.BarcodeMap
	:members:

:py:mod:`~enrich2.selection` --- Functional score calculation using SeqLib count data
=====================================================================================

.. py:module:: selection
	:synopsis: Functional score calculation using :py:class:`~enrich2.seqlib.SeqLib` count data.

This module provides class definitions for the :py:class:`~enrich2.selection.Selection` class. This is where functional scores are calculated from the :py:class:`~enrich2.seqlib.SeqLib` count data. For time series data, each time point in the selection can have multiple :py:class:`~enrich2.seqlib.SeqLib` assigned to it, in which case the counts for each element will be added together. Each time series selection must have a time point 0 (the "input library").

:py:class:`~enrich2.selection.Selection` class
----------------------------------------------------------
.. autoclass:: enrich2.selection.Selection
	:members:

:py:class:`~enrich2.selection.Selection` helpers
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autofunction:: enrich2.selection.regression_apply

:py:mod:`~enrich2.condition` --- Dummy class for GUI
=======================================================================

.. py:module:: condition
	:synopsis: Dummy class for GUI.

This module provides class definitions for the :py:class:`~enrich2.experiment.condition.Condition` classes. This class is required for proper GUI operation. All condition-related behaviors are in the :py:class:`~enrich2.experiment.Experiment` class.

:py:class:`~enrich2.condition.Condition` class
-----------------------------------------------------------
.. autoclass:: enrich2.condition.Condition
	:members:

:py:mod:`~enrich2.experiment` --- Aggregation of replicate selections
=======================================================================

.. py:module:: experiment
	:synopsis: Aggregation of replicate selections.

This module provides class definitions for the :py:class:`~enrich2.experiment.Experiment`. Functional scores for selections within the same condition are combined to generate a single functional score (and associated error) for each element in each experimental condition.

:py:class:`~enrich2.experiment.Experiment` class
--------------------------------------------------------------
.. autoclass:: enrich2.experiment.Experiment
	:members:

Enrich2 plotting
===================================================================

.. py:module:: plots
	:synopsis: Library for general Enrich2 plotting.

Text goes here.

.. automodule:: enrich2.plots
	:members:

Sequence-function map plotting
--------------------------------------------------------------------

.. py:module:: sfmap
	:synopsis: Library for sequence-function map plotting.

Text goes here.

.. automodule:: enrich2.sfmap
	:members:

Utility functions
====================================================================

Configuration object type detection
---------------------------------------------------------------------------

.. automodule:: enrich2.config_check
	:members:

Dataframe and index helper functions
----------------------------------------------------------------------------

.. automodule:: enrich2.dataframe
	:members:

.. _api-variant-helper:

Variant helper functions
----------------------------------------------------------------------------

.. automodule:: enrich2.variant
	:members: mutation_count, has_indel, has_unresolvable, protein_variant

HGVS_ variant regular expressions
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. autodata:: enrich2.variant.re_protein

.. autodata:: enrich2.variant.re_coding

.. autodata:: enrich2.variant.re_noncoding

Enrich2 entry points
====================================================================

.. automodule:: enrich2.main
	:members:


