Enrich2: deep mutational scanning data analysis
============================================================================

Enrich2 is a general software tool for processing, analyzing, and visualizing data from deep mutational scanning experiments.

The software is freely available from https://github.com/FowlerLab/Enrich2/ under the BSD 3-clause license.

For an example dataset, visit https://github.com/FowlerLab/Enrich2-Example/.

To cite Enrich2, please reference `A statistical framework for analyzing deep mutational scanning data <https://doi.org/10.1186/s13059-017-1272-5>`_.

Enrich2 was written by `Alan F Rubin <mailto:alan.rubin@wehi.edu.au>`_ |ORCID_icon| http://orcid.org/0000-0003-1474-605X

.. error:: Important notice for users of Enrich2 v1.0 or v1.1

Enrich2 v1.2.0 corrected an error in the software that, for most datasets, resulted in the standard errors for combined scores being over-estimated.
The counts, scores, and replicate-wise standard errors are unaffected.

If you have analyzed datasets that contain replicates with a previous version of Enrich2, the easiest way to get the correct standard error values is to delete the experiment HDF5_ file (the file name ends with ``'_exp.h5'``) and re-run the program.
This will recalculate combined scores and standard errors without redoing other parts of the analysis.

.. |ORCID_icon| image:: _static/iD_icon.png
    :target: http://orcid.org


.. toctree::
    :hidden:
    :maxdepth: 0
    
    installation
    introduction
    gui
    seqlib_config
    output
    plots
    notebooks
    api

