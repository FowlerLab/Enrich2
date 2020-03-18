[![DOI](https://zenodo.org/badge/69113902.svg)](https://zenodo.org/badge/latestdoi/69113902)

Enrich2
=======

Enrich2 is a general software tool for processing, analyzing, and visualizing data from deep mutational scanning experiments.
For more information or to cite Enrich2, please refer to [A statistical framework for analyzing deep mutational scanning data](https://doi.org/10.1186/s13059-017-1272-5).

[Enrich2 documentation](https://enrich2.readthedocs.io) is available on [Read the Docs](https://readthedocs.org/).

An example dataset is available at the [Enrich2-Example GitHub repository](https://github.com/FowlerLab/Enrich2-Example/).

Installation and dependencies
-----------------------------

Enrich2 runs on Python 2.7 and requires the following packages:

* [NumPy](http://www.numpy.org/) version 1.10.4 or higher
* [SciPy](http://www.scipy.org/) version 0.16.0 or higher
* [pandas](http://pandas.pydata.org/) version 0.18 or higher
* [PyTables](http://www.pytables.org/) version 3.2.0 or higher
* [Statsmodels](http://statsmodels.sourceforge.net/) version 0.6.1 or higher
* [matplotlib](http://matplotlib.org/) version 1.4.3 or higher

The configuration GUI requires [Tkinter](https://docs.python.org/2/library/tkinter.html).
Building a local copy of the documentation requires [Sphinx](http://sphinx-doc.org/).

We recommend using a scientific Python distribution such as [Anaconda](https://store.continuum.io/cshop/anaconda/) to install and manage dependencies.
A [Conda environment file](https://conda.io/docs/using/envs.html#use-environment-from-file) is included with the documentation.

To install Enrich2 using Anaconda (or Miniconda), git clone or download and unzip the Enrich2 code.
From the Enrich2 directory, create the environment, activate the environment, and install Enrich2:

    conda env create -f docs/_static/enrich2_env.yml
    conda activate enrich2
    python setup.py install

You should now be able to launch the Enrich2 graphical user interface by typing `enrich_gui` or the command line interface by typing `enrich_cmd`.

Questions?
----------

Please use the [GitHub Issue Tracker](https://github.com/FowlerLab/Enrich2/issues) to file bug reports or request features. 

Enrich2 was written by [Alan F Rubin](mailto:alan.rubin@wehi.edu.au).
