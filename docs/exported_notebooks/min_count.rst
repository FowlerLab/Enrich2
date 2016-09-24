
Selecting variants by input library count
-----------------------------------------

This notebook gets scores and standard errors for the variants in a
Selection that exceed a minimum count cutoff in the input time point,
and plots the relationship between each variant's score and input count.

.. code:: python

    % matplotlib inline

.. code:: python

    from __future__ import print_function
    import os.path
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

Open the Selection HDF5 file with the variants we are interested in.

.. code:: python

    my_store = pd.HDFStore(os.path.join(results_path, "Rep1_sel.h5"))

The ``pd.HDFStore.keys()`` method returns a list of all the tables in
this HDF5 file.

.. code:: python

    my_store.keys()




.. parsed-literal::

    ['/main/barcodemap',
     '/main/barcodes/counts',
     '/main/barcodes/counts_unfiltered',
     '/main/barcodes/log_ratios',
     '/main/barcodes/scores',
     '/main/barcodes/weights',
     '/main/synonymous/counts',
     '/main/synonymous/counts_unfiltered',
     '/main/synonymous/log_ratios',
     '/main/synonymous/scores',
     '/main/synonymous/weights',
     '/main/variants/counts',
     '/main/variants/counts_unfiltered',
     '/main/variants/log_ratios',
     '/main/variants/scores',
     '/main/variants/weights']



We will work with the "/main/variants/counts" table first. Enrich2
names the columns for counts ``c_n`` where ``n`` is the time point,
beginning with ``0`` for the input library.

We can use a query to extract the subset of variants in the table that
exceed the specified cutoff. Since we're only interested in variants,
we'll explicitly exclude the wild type. We will store the data we
extract in the ``variant_count`` data frame.

.. code:: python

    read_cutoff = 10

.. code:: python

    variant_counts = my_store.select('/main/variants/counts', where='c_0 > read_cutoff and index != WILD_TYPE_VARIANT')
    variant_counts




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>c_0</th>
          <th>c_2</th>
          <th>c_5</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>c.10G&gt;A (p.Ala4Arg), c.11C&gt;G (p.Ala4Arg), c.12T&gt;A (p.Ala4Arg)</th>
          <td>787.0</td>
          <td>106.0</td>
          <td>124.0</td>
        </tr>
        <tr>
          <th>c.10G&gt;A (p.Ala4Asn), c.11C&gt;A (p.Ala4Asn)</th>
          <td>699.0</td>
          <td>80.0</td>
          <td>114.0</td>
        </tr>
        <tr>
          <th>c.10G&gt;A (p.Ala4Asn), c.11C&gt;A (p.Ala4Asn), c.12T&gt;C (p.Ala4Asn)</th>
          <td>94.0</td>
          <td>8.0</td>
          <td>13.0</td>
        </tr>
        <tr>
          <th>c.10G&gt;A (p.Ala4Ile), c.11C&gt;T (p.Ala4Ile)</th>
          <td>1280.0</td>
          <td>137.0</td>
          <td>80.0</td>
        </tr>
        <tr>
          <th>c.10G&gt;A (p.Ala4Ile), c.11C&gt;T (p.Ala4Ile), c.12T&gt;A (p.Ala4Ile)</th>
          <td>717.0</td>
          <td>42.0</td>
          <td>27.0</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>c.9T&gt;A (p.=)</th>
          <td>327.0</td>
          <td>217.0</td>
          <td>284.0</td>
        </tr>
        <tr>
          <th>c.9T&gt;C (p.=)</th>
          <td>1947.0</td>
          <td>523.0</td>
          <td>1230.0</td>
        </tr>
        <tr>
          <th>c.9T&gt;C (p.=), c.49A&gt;T (p.Met17Ser), c.50T&gt;C (p.Met17Ser), c.51G&gt;A (p.Met17Ser)</th>
          <td>277.0</td>
          <td>43.0</td>
          <td>5.0</td>
        </tr>
        <tr>
          <th>c.9T&gt;C (p.=), c.62T&gt;C (p.Leu21Ser), c.63A&gt;T (p.Leu21Ser)</th>
          <td>495.0</td>
          <td>138.0</td>
          <td>55.0</td>
        </tr>
        <tr>
          <th>c.9T&gt;G (p.=)</th>
          <td>406.0</td>
          <td>18.0</td>
          <td>20.0</td>
        </tr>
      </tbody>
    </table>
    <p>1440 rows × 3 columns</p>
    </div>



The index of the data frame is the list of variants that exceeded the
cutoff.

.. code:: python

    variant_counts.index




.. parsed-literal::

    Index([u'c.10G>A (p.Ala4Arg), c.11C>G (p.Ala4Arg), c.12T>A (p.Ala4Arg)',
           u'c.10G>A (p.Ala4Asn), c.11C>A (p.Ala4Asn)',
           u'c.10G>A (p.Ala4Asn), c.11C>A (p.Ala4Asn), c.12T>C (p.Ala4Asn)',
           u'c.10G>A (p.Ala4Ile), c.11C>T (p.Ala4Ile)',
           u'c.10G>A (p.Ala4Ile), c.11C>T (p.Ala4Ile), c.12T>A (p.Ala4Ile)',
           u'c.10G>A (p.Ala4Ile), c.11C>T (p.Ala4Ile), c.12T>C (p.Ala4Ile)',
           u'c.10G>A (p.Ala4Lys), c.11C>A (p.Ala4Lys), c.12T>A (p.Ala4Lys)',
           u'c.10G>A (p.Ala4Met), c.11C>T (p.Ala4Met), c.12T>G (p.Ala4Met)',
           u'c.10G>A (p.Ala4Ser), c.11C>G (p.Ala4Ser)',
           u'c.10G>A (p.Ala4Ser), c.11C>G (p.Ala4Ser), c.12T>C (p.Ala4Ser)',
           ...
           u'c.8C>T (p.Ser3Phe), c.60C>T (p.=)',
           u'c.8C>T (p.Ser3Phe), c.9T>C (p.Ser3Phe)', u'c.90C>A (p.=)',
           u'c.90C>G (p.Ile30Met)', u'c.90C>T (p.=)', u'c.9T>A (p.=)',
           u'c.9T>C (p.=)',
           u'c.9T>C (p.=), c.49A>T (p.Met17Ser), c.50T>C (p.Met17Ser), c.51G>A (p.Met17Ser)',
           u'c.9T>C (p.=), c.62T>C (p.Leu21Ser), c.63A>T (p.Leu21Ser)',
           u'c.9T>G (p.=)'],
          dtype='object', length=1440)



We can use this index to get the scores for these variants by querying
the "/main/variants/scores" table. We'll store the result of the query
in a new data frame named ``variant_scores``, and keep only the score
and standard error (SE) columns.

.. code:: python

    variant_scores = my_store.select('/main/variants/scores', where='index in variant_counts.index')
    variant_scores = variant_scores[['score', 'SE']]
    variant_scores




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>score</th>
          <th>SE</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>c.10G&gt;A (p.Ala4Arg), c.11C&gt;G (p.Ala4Arg), c.12T&gt;A (p.Ala4Arg)</th>
          <td>-0.980091</td>
          <td>0.134873</td>
        </tr>
        <tr>
          <th>c.10G&gt;A (p.Ala4Asn), c.11C&gt;A (p.Ala4Asn)</th>
          <td>-0.972035</td>
          <td>0.268962</td>
        </tr>
        <tr>
          <th>c.10G&gt;A (p.Ala4Asn), c.11C&gt;A (p.Ala4Asn), c.12T&gt;C (p.Ala4Asn)</th>
          <td>-1.138667</td>
          <td>0.403767</td>
        </tr>
        <tr>
          <th>c.10G&gt;A (p.Ala4Ile), c.11C&gt;T (p.Ala4Ile)</th>
          <td>-1.875331</td>
          <td>0.014883</td>
        </tr>
        <tr>
          <th>c.10G&gt;A (p.Ala4Ile), c.11C&gt;T (p.Ala4Ile), c.12T&gt;A (p.Ala4Ile)</th>
          <td>-2.552289</td>
          <td>0.421699</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>c.9T&gt;A (p.=)</th>
          <td>0.705661</td>
          <td>0.774559</td>
        </tr>
        <tr>
          <th>c.9T&gt;C (p.=)</th>
          <td>0.438654</td>
          <td>0.014857</td>
        </tr>
        <tr>
          <th>c.9T&gt;C (p.=), c.49A&gt;T (p.Met17Ser), c.50T&gt;C (p.Met17Ser), c.51G&gt;A (p.Met17Ser)</th>
          <td>-1.930922</td>
          <td>1.085535</td>
        </tr>
        <tr>
          <th>c.9T&gt;C (p.=), c.62T&gt;C (p.Leu21Ser), c.63A&gt;T (p.Leu21Ser)</th>
          <td>-0.897249</td>
          <td>0.884321</td>
        </tr>
        <tr>
          <th>c.9T&gt;G (p.=)</th>
          <td>-2.314604</td>
          <td>0.671760</td>
        </tr>
      </tbody>
    </table>
    <p>1440 rows × 2 columns</p>
    </div>



Now that we're finished getting data out of the HDF5 file, we'll close
it.

.. code:: python

    my_store.close()

To more easily explore the relationship between input count and score,
we'll add a column to the ``variant_scores`` data frame that contains
input counts from the ``variant_counts`` data frame.

.. code:: python

    variant_scores['input_count'] = variant_counts['c_0']
    variant_scores




.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>score</th>
          <th>SE</th>
          <th>input_count</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>c.10G&gt;A (p.Ala4Arg), c.11C&gt;G (p.Ala4Arg), c.12T&gt;A (p.Ala4Arg)</th>
          <td>-0.980091</td>
          <td>0.134873</td>
          <td>787.0</td>
        </tr>
        <tr>
          <th>c.10G&gt;A (p.Ala4Asn), c.11C&gt;A (p.Ala4Asn)</th>
          <td>-0.972035</td>
          <td>0.268962</td>
          <td>699.0</td>
        </tr>
        <tr>
          <th>c.10G&gt;A (p.Ala4Asn), c.11C&gt;A (p.Ala4Asn), c.12T&gt;C (p.Ala4Asn)</th>
          <td>-1.138667</td>
          <td>0.403767</td>
          <td>94.0</td>
        </tr>
        <tr>
          <th>c.10G&gt;A (p.Ala4Ile), c.11C&gt;T (p.Ala4Ile)</th>
          <td>-1.875331</td>
          <td>0.014883</td>
          <td>1280.0</td>
        </tr>
        <tr>
          <th>c.10G&gt;A (p.Ala4Ile), c.11C&gt;T (p.Ala4Ile), c.12T&gt;A (p.Ala4Ile)</th>
          <td>-2.552289</td>
          <td>0.421699</td>
          <td>717.0</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>c.9T&gt;A (p.=)</th>
          <td>0.705661</td>
          <td>0.774559</td>
          <td>327.0</td>
        </tr>
        <tr>
          <th>c.9T&gt;C (p.=)</th>
          <td>0.438654</td>
          <td>0.014857</td>
          <td>1947.0</td>
        </tr>
        <tr>
          <th>c.9T&gt;C (p.=), c.49A&gt;T (p.Met17Ser), c.50T&gt;C (p.Met17Ser), c.51G&gt;A (p.Met17Ser)</th>
          <td>-1.930922</td>
          <td>1.085535</td>
          <td>277.0</td>
        </tr>
        <tr>
          <th>c.9T&gt;C (p.=), c.62T&gt;C (p.Leu21Ser), c.63A&gt;T (p.Leu21Ser)</th>
          <td>-0.897249</td>
          <td>0.884321</td>
          <td>495.0</td>
        </tr>
        <tr>
          <th>c.9T&gt;G (p.=)</th>
          <td>-2.314604</td>
          <td>0.671760</td>
          <td>406.0</td>
        </tr>
      </tbody>
    </table>
    <p>1440 rows × 3 columns</p>
    </div>



Now that all the information is in a single data frame, we can make a
plot of score vs. input count. This example uses functions and colors
from the Enrich2 plotting library. Taking the log10 of the counts makes
the data easier to visualize.

.. code:: python

    fig, ax = plt.subplots()
    enrich_plot.configure_axes(ax, xgrid=True)
    ax.plot(np.log10(variant_scores['input_count']), 
            variant_scores['score'], 
            linestyle='none', marker='.', alpha=0.6,
            color=enrich_plot.plot_colors['bright4'])
    ax.set_xlabel("log10(Input Count)")
    ax.set_ylabel("Variant Score")




.. parsed-literal::

    <matplotlib.text.Text at 0x9e796a0>




.. image:: _static/notebook_plots/min_count_plot.png


