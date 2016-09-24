
Selecting variants by number of unique barcodes
-----------------------------------------------

This notebook gets scores for the variants in an Experiment that are
linked to multiple barcodes, and plots the relationship between each
variant's score and number of unique barcodes.

.. code:: python

    % matplotlib inline

.. code:: python

    from __future__ import print_function
    import os.path
    from collections import Counter
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from enrich2.variant import WILD_TYPE_VARIANT
    import enrich2.plots as enrich_plot
    pd.set_option("display.max_rows", 10) # rows shown when pretty-printing

Modify the ``results_path`` variable in the next cell to match the
output directory of your Enrich2-Example dataset.

.. code:: python

    results_path = "/path/to/Enrich2-Example/Results/"

Open the Experiment HDF5 file.

.. code:: python

    my_store = pd.HDFStore(os.path.join(results_path, "BRCA1_Example_exp.h5"))

The ``pd.HDFStore.keys()`` method returns a list of all the tables in
this HDF5 file.

.. code:: python

    my_store.keys()




.. parsed-literal::

    ['/main/barcodemap',
     '/main/barcodes/counts',
     '/main/barcodes/scores',
     '/main/barcodes/scores_shared',
     '/main/barcodes/scores_shared_full',
     '/main/synonymous/counts',
     '/main/synonymous/scores',
     '/main/synonymous/scores_pvalues_wt',
     '/main/synonymous/scores_shared',
     '/main/synonymous/scores_shared_full',
     '/main/variants/counts',
     '/main/variants/scores',
     '/main/variants/scores_pvalues_wt',
     '/main/variants/scores_shared',
     '/main/variants/scores_shared_full']



First we will work with the barcode-variant map for this analysis,
stored in the "/main/barcodemap" table. The index is the barcode and it
has a single column for the variant HGVS string.

.. code:: python

    bcm = my_store['/main/barcodemap']
    bcm




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>value</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>TTTTTTGTGTCTGTGA</th>
          <td>_wt</td>
        </tr>
        <tr>
          <th>GGGCACGTCTTTATAG</th>
          <td>_wt</td>
        </tr>
        <tr>
          <th>GTTACTGGTTAGTATT</th>
          <td>_wt</td>
        </tr>
        <tr>
          <th>GTTACTTGATCCGACC</th>
          <td>_wt</td>
        </tr>
        <tr>
          <th>GTTAGATGGATGTACG</th>
          <td>_wt</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
        </tr>
        <tr>
          <th>GAGTACTTTTTTGATT</th>
          <td>c.9T&gt;C (p.=), c.62T&gt;C (p.Leu21Ser), c.63A&gt;T (p...</td>
        </tr>
        <tr>
          <th>ATGATGACGTGTCTTG</th>
          <td>c.9T&gt;G (p.=)</td>
        </tr>
        <tr>
          <th>TCACCGGAACGTTGGT</th>
          <td>c.9T&gt;G (p.=)</td>
        </tr>
        <tr>
          <th>TGACGATGTTGCATTT</th>
          <td>c.9T&gt;G (p.=), c.17G&gt;C (p.Arg6Pro), c.18C&gt;T (p....</td>
        </tr>
        <tr>
          <th>GTTATCAGCGCCCCTT</th>
          <td>c.9T&gt;G (p.=), c.85T&gt;G (p.Leu29Gly), c.86T&gt;G (p...</td>
        </tr>
      </tbody>
    </table>
    <p>20325 rows × 1 columns</p>
    </div>



To find out how many unique barcodes are linked to each variant, we'll
count the number of times each variant appears in the barcode-variant
map using a `Counter data
structure <https://docs.python.org/2/library/collections.html#counter-objects>`__.
We'll then output the top ten variants by number of unique barcodes.

.. code:: python

    variant_bcs = Counter(bcm['value'])
    variant_bcs.most_common(10)




.. parsed-literal::

    [('_wt', 5844),
     ('c.63A>T (p.Leu21Phe)', 109),
     ('c.39C>A (p.=)', 91),
     ('c.61T>A (p.Leu21Ile), c.63A>T (p.Leu21Ile)', 77),
     ('c.62T>A (p.Leu21Tyr), c.63A>T (p.Leu21Tyr)', 77),
     ('c.63A>G (p.=)', 73),
     ('c.72C>A (p.=)', 72),
     ('c.62T>G (p.Leu21Cys), c.63A>T (p.Leu21Cys)', 71),
     ('c.13C>A (p.Leu5Ile)', 70),
     ('c.62T>A (p.Leu21Ter)', 63)]



Next we'll turn the Counter into a data frame.

.. code:: python

    bc_counts = pd.DataFrame(variant_bcs.most_common(), columns=['variant', 'barcodes'])
    bc_counts




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>variant</th>
          <th>barcodes</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>_wt</td>
          <td>5844</td>
        </tr>
        <tr>
          <th>1</th>
          <td>c.63A&gt;T (p.Leu21Phe)</td>
          <td>109</td>
        </tr>
        <tr>
          <th>2</th>
          <td>c.39C&gt;A (p.=)</td>
          <td>91</td>
        </tr>
        <tr>
          <th>3</th>
          <td>c.61T&gt;A (p.Leu21Ile), c.63A&gt;T (p.Leu21Ile)</td>
          <td>77</td>
        </tr>
        <tr>
          <th>4</th>
          <td>c.62T&gt;A (p.Leu21Tyr), c.63A&gt;T (p.Leu21Tyr)</td>
          <td>77</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>1958</th>
          <td>c.77G&gt;T (p.Cys26Leu), c.78C&gt;A (p.Cys26Leu), c....</td>
          <td>1</td>
        </tr>
        <tr>
          <th>1959</th>
          <td>c.67T&gt;A (p.Cys23Ile), c.68G&gt;T (p.Cys23Ile), c....</td>
          <td>1</td>
        </tr>
        <tr>
          <th>1960</th>
          <td>c.41T&gt;C (p.Ile14Thr), c.48T&gt;G (p.=)</td>
          <td>1</td>
        </tr>
        <tr>
          <th>1961</th>
          <td>c.55A&gt;C (p.Lys19Leu), c.56A&gt;T (p.Lys19Leu), c....</td>
          <td>1</td>
        </tr>
        <tr>
          <th>1962</th>
          <td>c.50T&gt;C (p.Met17Thr), c.78C&gt;T (p.=)</td>
          <td>1</td>
        </tr>
      </tbody>
    </table>
    <p>1963 rows × 2 columns</p>
    </div>



The data frame has the information we want, but it will be easier to use
later if it's indexed by variant rather than row number.

.. code:: python

    bc_counts.index = bc_counts['variant']
    bc_counts.index.name = None
    del bc_counts['variant']
    bc_counts




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>barcodes</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>_wt</th>
          <td>5844</td>
        </tr>
        <tr>
          <th>c.63A&gt;T (p.Leu21Phe)</th>
          <td>109</td>
        </tr>
        <tr>
          <th>c.39C&gt;A (p.=)</th>
          <td>91</td>
        </tr>
        <tr>
          <th>c.61T&gt;A (p.Leu21Ile), c.63A&gt;T (p.Leu21Ile)</th>
          <td>77</td>
        </tr>
        <tr>
          <th>c.62T&gt;A (p.Leu21Tyr), c.63A&gt;T (p.Leu21Tyr)</th>
          <td>77</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
        </tr>
        <tr>
          <th>c.77G&gt;T (p.Cys26Leu), c.78C&gt;A (p.Cys26Leu), c.81G&gt;T (p.=)</th>
          <td>1</td>
        </tr>
        <tr>
          <th>c.67T&gt;A (p.Cys23Ile), c.68G&gt;T (p.Cys23Ile), c.69T&gt;A (p.Cys23Ile)</th>
          <td>1</td>
        </tr>
        <tr>
          <th>c.41T&gt;C (p.Ile14Thr), c.48T&gt;G (p.=)</th>
          <td>1</td>
        </tr>
        <tr>
          <th>c.55A&gt;C (p.Lys19Leu), c.56A&gt;T (p.Lys19Leu), c.57A&gt;G (p.Lys19Leu), c.81G&gt;T (p.=)</th>
          <td>1</td>
        </tr>
        <tr>
          <th>c.50T&gt;C (p.Met17Thr), c.78C&gt;T (p.=)</th>
          <td>1</td>
        </tr>
      </tbody>
    </table>
    <p>1963 rows × 1 columns</p>
    </div>



We'll use a cutoff to choose variants with a minimum number of unique
barcodes, and store this subset in a new index. We'll also exclude the
wild type by dropping the first entry of the index.

.. code:: python

    bc_cutoff = 10

.. code:: python

    multi_bc_variants = bc_counts.loc[bc_counts['barcodes']  >= bc_cutoff].index[1:]
    multi_bc_variants




.. parsed-literal::

    Index([u'c.63A>T (p.Leu21Phe)', u'c.39C>A (p.=)',
           u'c.61T>A (p.Leu21Ile), c.63A>T (p.Leu21Ile)',
           u'c.62T>A (p.Leu21Tyr), c.63A>T (p.Leu21Tyr)', u'c.63A>G (p.=)',
           u'c.72C>A (p.=)', u'c.62T>G (p.Leu21Cys), c.63A>T (p.Leu21Cys)',
           u'c.13C>A (p.Leu5Ile)', u'c.62T>A (p.Leu21Ter)',
           u'c.63A>C (p.Leu21Phe)',
           ...
           u'c.88A>C (p.Ile30Arg), c.89T>G (p.Ile30Arg), c.90C>T (p.Ile30Arg)',
           u'c.76T>A (p.Cys26Lys), c.77G>A (p.Cys26Lys), c.78C>G (p.Cys26Lys)',
           u'c.22G>A (p.Glu8Ile), c.23A>T (p.Glu8Ile), c.24A>T (p.Glu8Ile)',
           u'c.49A>T (p.Met17Ser), c.50T>C (p.Met17Ser), c.51G>A (p.Met17Ser)',
           u'c.64G>A (p.Glu22Arg), c.65A>G (p.Glu22Arg)',
           u'c.77G>C (p.Cys26Ser), c.78C>G (p.Cys26Ser)',
           u'c.29T>A (p.Val10Glu), c.30A>G (p.Val10Glu)',
           u'c.50T>A (p.Met17Asn), c.51G>T (p.Met17Asn)',
           u'c.61T>A (p.Leu21Thr), c.62T>C (p.Leu21Thr), c.63A>G (p.Leu21Thr)',
           u'c.49A>G (p.Met17Ala), c.50T>C (p.Met17Ala)'],
          dtype='object', length=504)



We can use this index to get condition-level scores for these variants
by querying the "/main/variants/scores" table. Since we are working with
an Experiment HDF5 file, the data frame column names are a MultiIndex
with two levels, one for experimental conditions and one for data values
(see the `pandas
documentation <http://pandas.pydata.org/pandas-docs/stable/advanced.html>`__
for more information).

.. code:: python

    multi_bc_scores = my_store.select('/main/variants/scores', where='index in multi_bc_variants')
    multi_bc_scores




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr>
          <th>condition</th>
          <th colspan="3" halign="left">E3</th>
        </tr>
        <tr>
          <th>value</th>
          <th>SE</th>
          <th>epsilon</th>
          <th>score</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>c.10G&gt;A (p.Ala4Thr)</th>
          <td>1.435686e-01</td>
          <td>3.469447e-18</td>
          <td>-0.238174</td>
        </tr>
        <tr>
          <th>c.10G&gt;T (p.Ala4Ser)</th>
          <td>1.456404e-29</td>
          <td>2.087042e-57</td>
          <td>-0.177983</td>
        </tr>
        <tr>
          <th>c.11C&gt;A (p.Ala4Asp)</th>
          <td>5.309592e-01</td>
          <td>1.110223e-16</td>
          <td>0.027898</td>
        </tr>
        <tr>
          <th>c.13C&gt;A (p.Leu5Ile)</th>
          <td>1.333666e-01</td>
          <td>0.000000e+00</td>
          <td>-0.623652</td>
        </tr>
        <tr>
          <th>c.13C&gt;A (p.Leu5Ser), c.14T&gt;G (p.Leu5Ser)</th>
          <td>3.612046e-01</td>
          <td>2.775558e-17</td>
          <td>0.657916</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>c.89T&gt;G (p.Ile30Ser), c.90C&gt;T (p.Ile30Ser)</th>
          <td>6.069463e-01</td>
          <td>0.000000e+00</td>
          <td>-0.826140</td>
        </tr>
        <tr>
          <th>c.8C&gt;A (p.Ser3Tyr)</th>
          <td>3.785724e-01</td>
          <td>2.775558e-17</td>
          <td>-1.440477</td>
        </tr>
        <tr>
          <th>c.8C&gt;T (p.Ser3Phe)</th>
          <td>8.669053e-02</td>
          <td>2.602085e-18</td>
          <td>-0.091250</td>
        </tr>
        <tr>
          <th>c.90C&gt;A (p.=)</th>
          <td>9.681631e-02</td>
          <td>5.204170e-18</td>
          <td>-0.217977</td>
        </tr>
        <tr>
          <th>c.90C&gt;T (p.=)</th>
          <td>5.450037e-117</td>
          <td>8.793373e-229</td>
          <td>0.805631</td>
        </tr>
      </tbody>
    </table>
    <p>486 rows × 3 columns</p>
    </div>



There are fewer rows in ``multi_bc_scores`` than in
``multi_bc_variants`` because some of the variants were not scored in
all replicate selections, and therefore do not have a condition-level
score.

Now that we're finished getting data out of the HDF5 file, we'll close
it.

.. code:: python

    my_store.close()

We'll add a column to the ``bc_counts`` data frame that contains scores
from the ``multi_bc_scores`` data frame. To reference a column in a data
frame with a MultiIndex, we need to specify all column levels.

.. code:: python

    bc_counts['score'] = multi_bc_scores['E3', 'score']
    bc_counts




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>barcodes</th>
          <th>score</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>_wt</th>
          <td>5844</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>c.63A&gt;T (p.Leu21Phe)</th>
          <td>109</td>
          <td>1.387659</td>
        </tr>
        <tr>
          <th>c.39C&gt;A (p.=)</th>
          <td>91</td>
          <td>-0.189253</td>
        </tr>
        <tr>
          <th>c.61T&gt;A (p.Leu21Ile), c.63A&gt;T (p.Leu21Ile)</th>
          <td>77</td>
          <td>-1.031977</td>
        </tr>
        <tr>
          <th>c.62T&gt;A (p.Leu21Tyr), c.63A&gt;T (p.Leu21Tyr)</th>
          <td>77</td>
          <td>0.310854</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>c.77G&gt;T (p.Cys26Leu), c.78C&gt;A (p.Cys26Leu), c.81G&gt;T (p.=)</th>
          <td>1</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>c.67T&gt;A (p.Cys23Ile), c.68G&gt;T (p.Cys23Ile), c.69T&gt;A (p.Cys23Ile)</th>
          <td>1</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>c.41T&gt;C (p.Ile14Thr), c.48T&gt;G (p.=)</th>
          <td>1</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>c.55A&gt;C (p.Lys19Leu), c.56A&gt;T (p.Lys19Leu), c.57A&gt;G (p.Lys19Leu), c.81G&gt;T (p.=)</th>
          <td>1</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>c.50T&gt;C (p.Met17Thr), c.78C&gt;T (p.=)</th>
          <td>1</td>
          <td>NaN</td>
        </tr>
      </tbody>
    </table>
    <p>1963 rows × 2 columns</p>
    </div>



Many rows in ``bc_counts`` are missing scores (displayed as NaN) because
those variants were not in ``multi_bc_scores``. We'll drop them before
continuing.

.. code:: python

    bc_counts.dropna(inplace=True)
    bc_counts




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>barcodes</th>
          <th>score</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>c.63A&gt;T (p.Leu21Phe)</th>
          <td>109</td>
          <td>1.387659</td>
        </tr>
        <tr>
          <th>c.39C&gt;A (p.=)</th>
          <td>91</td>
          <td>-0.189253</td>
        </tr>
        <tr>
          <th>c.61T&gt;A (p.Leu21Ile), c.63A&gt;T (p.Leu21Ile)</th>
          <td>77</td>
          <td>-1.031977</td>
        </tr>
        <tr>
          <th>c.62T&gt;A (p.Leu21Tyr), c.63A&gt;T (p.Leu21Tyr)</th>
          <td>77</td>
          <td>0.310854</td>
        </tr>
        <tr>
          <th>c.63A&gt;G (p.=)</th>
          <td>73</td>
          <td>-0.406277</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>c.64G&gt;A (p.Glu22Arg), c.65A&gt;G (p.Glu22Arg)</th>
          <td>10</td>
          <td>-2.577200</td>
        </tr>
        <tr>
          <th>c.77G&gt;C (p.Cys26Ser), c.78C&gt;G (p.Cys26Ser)</th>
          <td>10</td>
          <td>-3.497939</td>
        </tr>
        <tr>
          <th>c.50T&gt;A (p.Met17Asn), c.51G&gt;T (p.Met17Asn)</th>
          <td>10</td>
          <td>-1.737378</td>
        </tr>
        <tr>
          <th>c.61T&gt;A (p.Leu21Thr), c.62T&gt;C (p.Leu21Thr), c.63A&gt;G (p.Leu21Thr)</th>
          <td>10</td>
          <td>-1.307432</td>
        </tr>
        <tr>
          <th>c.49A&gt;G (p.Met17Ala), c.50T&gt;C (p.Met17Ala)</th>
          <td>10</td>
          <td>-1.958962</td>
        </tr>
      </tbody>
    </table>
    <p>486 rows × 2 columns</p>
    </div>



Now that we have a data frame containing the subset of variants we're
interested in, we can make a plot of score vs. number of unique
barcodes. This example uses functions and colors from the Enrich2
plotting library.

.. code:: python

    fig, ax = plt.subplots()
    enrich_plot.configure_axes(ax, xgrid=True)
    ax.plot(bc_counts['barcodes'], 
            bc_counts['score'], 
            linestyle='none', marker='.', alpha=0.6,
            color=enrich_plot.plot_colors['bright5'])
    ax.set_xlabel("Unique Barcodes")
    ax.set_ylabel("Variant Score")




.. parsed-literal::

    <matplotlib.text.Text at 0xd91fe80>




.. image:: _static/notebook_plots/unique_barcodes_plot.png


