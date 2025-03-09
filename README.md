[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14681278.svg)](https://doi.org/10.5281/zenodo.14681278)
[![PyPI version](https://badge.fury.io/py/Enrich2.svg)](https://badge.fury.io/py/Enrich2)

Enrich2
=======

Enrich2 is a general software tool for processing, analyzing, and visualizing data from deep mutational scanning experiments.
For more information or to cite Enrich2, please refer to [A statistical framework for analyzing deep mutational scanning data](https://doi.org/10.1186/s13059-017-1272-5).

[Enrich2 documentation](https://enrich2.readthedocs.io) is available on [Read the Docs](https://readthedocs.org/).

An example dataset is available at the [Enrich2-Example GitHub repository](https://github.com/FowlerLab/Enrich2-Example/).

Thanks to the efforts of [Chris Macdonald](https://github.com/odcambc), Enrich2 is now able to run under modern versions of Python as of v2.0.0!

Installation and dependencies
-----------------------------

Enrich2 runs on Python 3 (v2.0.0 and higher) and requires the following packages:

* [NumPy](http://www.numpy.org/)
* [SciPy](http://www.scipy.org/)
* [pandas](http://pandas.pydata.org/)
* [PyTables](http://www.pytables.org/)
* [Statsmodels](http://statsmodels.sourceforge.net/)
* [matplotlib](http://matplotlib.org/)
* [fqfa](https://fqfa.readthedocs.io/)

The configuration GUI requires [Tkinter](https://docs.python.org/2/library/tkinter.html).
Building a local copy of the documentation requires [Sphinx](http://sphinx-doc.org/).

Enrich2 can be installed in a new virtual environment using pip:

    python3 -m venv e2env
    source e2env/bin/activate
    pip install enrich2

You should now be able to launch the Enrich2 graphical user interface by typing `enrich_gui` or the command line interface by typing `enrich_cmd`.

For additional information consult the [online documentation](https://enrich2.readthedocs.io/).

Questions?
----------

Please use the [GitHub Issue Tracker](https://github.com/FowlerLab/Enrich2/issues) to file bug reports or request features. 

Enrich2 was written by [Alan F Rubin](mailto:alan.rubin@wehi.edu.au).
