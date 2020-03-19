from __future__ import print_function
import collections
import logging
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from statsmodels.nonparametric.kde import KDEUnivariate

logger = logging.getLogger(__name__)

tick_params = {"labelsize": 6, "pad": 2}

label_params = {"fontsize": 8, "labelpad": 2}

legend_params = {"fontsize": 8}

letter_params = {
    "fontsize": 10,
    "weight": "bold",
    "horizontalalignment": "left",
    "verticalalignment": "center",
}

#: Default colors for Enrich2 plots.
#: Standard defaults are based on ColorBrewer Set1 and Pastel1.
plot_colors = {
    "bright1": "#377eb8",
    "bright2": "#4daf4a",
    "bright3": "#ff7f00",
    "bright4": "#984ea3",
    "bright5": "#e41a1c",
    "bright6": "#a65628",
    "pastel1": "#b3cde3",
    "pastel2": "#ccebc5",
    "pastel3": "#fed9a6",
    "pastel4": "#decbe4",
    "pastel5": "#fbb4ae",
    "pastel6": "#e5d8bd",
    "accent1": "#000000",
    "accent2": "#545454",
    "accent3": "#bebebe",
    "accent4": "#999999",
    "missing": "#808080",
}


#: Default matplotlib color maps.
#: ``'diverging'`` is used for sequence-function maps.
#: ``'sequential'`` is used for library diversity plots.
#: ``'hexbin'`` is used for volcano plots.
plot_cmaps = {"diverging": "RdBu_r", "sequential": "BuPu", "hexbin": "YlGnBu_r"}


#: Default number of histogram bins.
DEFAULT_BINS = 30

#: Default number of hexbins in each direction (x and y).
DEFAULT_GRIDSIZE = 50


def configure_axes(ax, xgrid=False):
    """
    Apply the standard Enrich2 style to the axes *ax*. This includes invisible 
    tick marks, horizontal (and optionally vertical) grey grid lines at major 
    tick intervals, and scientific notation on numerical axes.

    *ax* is the axes object that will have its properties set.

    *xgrid* sets the drawing of vertical grid lines that match the horizontal 
    grid lines.

    """
    # turn off all tick marks
    ax.tick_params(
        axis="both", which="both", bottom=False, top=False, left=False, right=False
    )
    # turn on y-axis grid (horizontal lines)
    ax.yaxis.grid(color=plot_colors["accent3"], linestyle="-")
    # optionally turn on x-axis grid (vertical lines)
    if xgrid:  # for scatterplots, hexbin
        ax.xaxis.grid(color=plot_colors["accent3"], linestyle="-")
    else:  # for barplots, histograms, etc.
        ax.xaxis.grid(False)
    # put the axis and grid behind any plot elements
    ax.set_axisbelow(True)
    # use scientific notation for large or small tick label values
    ax.ticklabel_format(style="sci", scilimits=(-3, 3))


def fit_axes(ax, xvalues, yvalues, slope, intercept, xlabels=None):
    """
    Create a scatterplot with fit line for a particular element's 
    regression-based score. This is used for both WLS and OLS scoring methods.

    *ax* is the axes object used for plotting.

    *xvalues* are the regression's x-values (normalized time points).

    *yvalues* are the regression's y-values (log ratios).

    *slope* is the slope of the fit line.

    *intercept* is the intercept of the fit line.

    *xlabels* are optional tick labels for the x-axis (usually the 
    pre-normalization time points). If *xlabels* is ``None``, the *xvalues* 
    will be used as tick labels.

    """
    # set appropriate x axis limits
    xmin = sorted(xvalues)[0]
    xmax = sorted(xvalues)[-1]
    xlim_margin = (xmax - xmin) * 0.1
    ax.set_xlim([xmin - xlim_margin, xmax + xlim_margin])

    # add axis labels
    ax.set_xlabel("Time Point")
    ax.set_ylabel("Log Ratio")

    # add x tick labels
    ax.set_xticks(xvalues)
    if xlabels is None:
        ax.set_xticklabels(xvalues)
    else:
        ax.set_xticklabels(xlabels)

    # plot the log ratios
    ax.scatter(xvalues, yvalues, color=plot_colors["bright1"], marker="D")

    # plot residual lines
    fit_yvalues = [x * slope + intercept for x in xvalues]
    for x, y1, y2 in zip(xvalues, yvalues, fit_yvalues):
        ax.plot(
            (x, x), (y1, y2), color=plot_colors["bright1"], linestyle="-", alpha=0.5
        )

    # plot the fit line
    ax.plot(
        ax.get_xlim(),
        [x * slope + intercept for x in ax.get_xlim()],
        color=plot_colors["bright3"],
        linestyle="--",
    )


def fit_axes_text(ax, cornertext=None, centertext=None, centerwrap=40):
    """
    Add text to a fit plot in either the bottom left corner and/or the top
    middle.

    *ax* is the axes object containing the fit plot.

    *cornertext* is displayed in the lower left corner of the plot. Typically
    this includes details about the fit (slope, standard error, etc.). If it
    is ``None``, no text is displayed.

    *centertext* is displayed in a smaller font in the upper center of the
    plot. Typically this is the HGVS_ variant string or the barcode for the
    element whose fit is plotted. The string is split on commas into a
    multi-line string as required given its length (see *centerwrap*).

    *centerwrap* is the maximum length of each line in the *centertext* string.

    """
    if cornertext is not None:
        ax.text(
            0.05,
            0.05,
            cornertext,
            horizontalalignment="left",
            verticalalignment="bottom",
            transform=ax.transAxes,
        )

    if centertext is not None:
        splits = centertext.split(",")
        splits = [x.strip() for x in splits]
        length = len(splits[0])
        text = splits[0]
        for x in splits[1:]:
            if length + len(x) < centerwrap:
                text += ", " + x
                length += 2 + len(x)
            else:
                text += ",\n" + x
                length = len(x)
        ax.text(
            0.5,
            0.95,
            text,
            horizontalalignment="center",
            verticalalignment="top",
            transform=ax.transAxes,
            fontsize=8,
        )


def volcano_plot(df, pdf, title=None, colors=None, log_bins=True, logp_max=10):
    """
    Create a volcano plot of p-values vs. functional scores.

    *df* is a DataFrame containing ``'score'`` and ``'pvalue_raw'`` columns.
        Each element (barcode, variant, etc.) has its own row. Other columns
        will be ignored.

    *pdf* is the open PdfPages file object.

    *title* is the plot title. If ``None``, no title will be displayed.

    *colors* is the name of the `matplotlib cmap`_ used.
        If *colors* is ``None``, the default hexbin cmap will be used.

    *log_bins* sets whether or not the counts in each hexbin should be log10
    transformed.

    *logp_max* sets the maximum minus log-transformed p-value that will be
        plotted. p-values smaller than this will be set to the given max.

    """
    df = df.loc[:, ("pvalue_raw", "score")]
    df.dropna(axis="index", how="any", inplace=True)
    df["pvalue_raw"] = -np.log10(df["pvalue_raw"])
    df.loc[df["pvalue_raw"] > logp_max, "pvalue_raw"] = logp_max  # set max

    # create the figure
    fig, ax = plt.subplots()
    fig.set_tight_layout(True)
    configure_axes(ax, xgrid=True)

    # make the plot
    if log_bins:
        bins = "log"
    else:
        bins = None
    hexbin = ax.hexbin(
        x=df["score"].values,
        y=df["pvalue_raw"].values,
        cmap=colors,
        bins=bins,
        mincnt=1,
        edgecolor="none",
        gridsize=DEFAULT_GRIDSIZE,
    )

    # set the labels
    if title is not None:
        ax.set_title(title)
    ax.set_ylabel("-log10(raw p-value)")
    ax.set_xlabel("Score")

    # space the boundaries out so the points aren't on the edges
    ax.set_xlim([1.05 * df["score"].min(), 1.05 * df["score"].max()])
    ax.set_ylim([0.0, 1.05 * df["pvalue_raw"].max()])

    # create color bar
    cbar = fig.colorbar(hexbin)
    if log_bins:
        cbar.set_label("log10(Count)")
    else:
        cbar.set_label("Count")
    cbar.ax.tick_params(bottom=False, top=False, left=False, right=False)

    # save and clean up
    pdf.savefig(fig)
    plt.close(fig)


def barcodemap_plot(
    obj, pdf, log=False, bins=DEFAULT_BINS, color=plot_colors["bright3"]
):
    """
    Plot the number of barcodes assigned to each variant.

    Args:
        obj (:py:class:`~enrich2.barcodevariant.BcvSeqLib`): object with counts
            to plot

        pdf (|mpl_PdfPages|): destination file handle

        log (bool): whether to log10-transform the counts

        bins (int): number of histogram bins

        color (str): histogram bar color
    """
    try:
        data = collections.Counter(obj.store["/raw/barcodemap"]["value"])
    except KeyError:
        logger.warning(
            "Failed to generate barcode-variant histogram "
            "(barcode map data not found)"
        )
        return

    if len(data.keys()) <= 1:
        logger.warning("Not enough elements to make barcodemap plot")
        return

    # create the plot and set up the axes
    fig, ax = plt.subplots()
    fig.set_tight_layout(True)
    configure_axes(ax)

    # plot the histogram
    ax.hist(data.values(), bins=bins, log=log, color=color)

    # set the title and axes labels
    ax.set_title("Barcodes per Variant\n{}".format(obj.name))
    if log:
        ax.set_ylabel("log(Variants)")
    else:
        ax.set_ylabel("Variants")
    ax.set_xlabel("Barcode Count")

    # save and clean up
    pdf.savefig()
    plt.close(fig)


def overlap_merge_plot(obj, pdf):
    """
    Plot the location of all mismatches in the overlap region (split into
    resolved and unresolved) and the location of the first mismatch in the
    overlap region.

    Args:
        obj (OverlapSeqLib): object with counts to plot

        pdf (|mpl_PdfPages|): destination file handle
    """
    # read the data
    try:
        data = obj.store["/raw/overlap_mismatches"]
    except KeyError:
        logger.warning(
            "Failed to generate overlap mismatch plot " "(mismatch data not found)"
        )
        return

    # create the plot and set up the axes
    fig, ax = plt.subplots()
    fig.set_tight_layout(True)
    configure_axes(ax)

    # plotting constants
    width = 0.8  # width of the bars
    xpos = np.arange(len(data.index))  # x values starting at 0

    # plot the stacked barplot for unresolved and resolved mismatches
    ax.bar(
        xpos,
        data["resolved"],
        width,
        color=plot_colors["bright1"],
        label="Resolved Mismatches",
    )
    ax.bar(
        xpos,
        data["unresolved"],
        width,
        bottom=data["resolved"],
        color=plot_colors["bright2"],
        label="Unresolved Mismatches",
    )

    # plot the barplot for first mismatches
    ax.bar(
        xpos,
        data["first"],
        width,
        color=plot_colors["bright3"],
        label="First Mismatches",
    )

    # set title, axes labels, and figure legend
    ax.set_title("Overlap Mismatches\n{}".format(obj.name))
    ax.set_ylabel("Count")
    ax.set_xlabel("Nucleotide Position")
    ax.set_xticks(xpos + width / 2.0)
    ax.set_xticklabels(data.index, rotation=90)
    ax.legend(loc=4)  # bottom-right corner

    # save this figure and set up for next figure
    pdf.savefig(fig)
    plt.close(fig)


def counts_plot(
    seqlib, label, pdf, log=False, bins=DEFAULT_BINS, color=plot_colors["bright1"]
):
    """
    Plot a histogram of processed counts. Excludes counts for wild type and
    synonymous variants.

    Args:
        seqlib (SeqLib): object with counts to plot

        label (str): data label (barcode, variant, etc.)

        pdf (|mpl_PdfPages|): destination file handle

        bins (int): number of histogram bins

        color (str): histogram bar color
    """
    data = seqlib.store.select(
        "/main/{}/counts".format(label),
        where="columns=count & "
        "index != WILD_TYPE_VARIANT & "
        "index != SYNONYMOUS_VARIANT",
    ).values

    if len(data) <= 1:
        logger.warning("Not enough elements to make '{}' counts plot" "".format(label))
        return

    # create the plot and set up the axes
    fig, ax = plt.subplots()
    fig.set_tight_layout(True)
    configure_axes(ax)

    # plot the histogram
    ax.hist(data, bins=bins, log=log, color=color)

    # set title and axes labels
    ax.set_title("{} Counts\n{}".format(label.title(), seqlib.name))
    ax.set_ylabel("{}".format(label.title()))
    if log:
        ax.set_xlabel("log10(Count)")
    else:
        ax.set_xlabel("Count")

    # save figure and clean up
    pdf.savefig()
    plt.close()


def weights_plot(selection, label, pdf):
    """
    Create a boxplot of weights for weighted least squares (WLS) regression 
    scoring. Weights are converted to proportional weights (0 to 1) before 
    plotting.

    *selection* is a :py:class:`~selection.Selection`

    *label* is the data label (barcode, variant, etc.)

    *pdf* is an open PdfPages instance

    """
    # get the data and drop NA values
    data = selection.store.select(
        "/main/{}/weights".format(label),
        where="index != WILD_TYPE_VARIANT & " "index != SYNONYMOUS_VARIANT",
    )
    data.dropna(axis="index", how="any", inplace=True)

    # normalize the weights to [0..1]
    data = data.values / data.values.sum(axis=1)[:, None]

    # create the figure and axes
    fig, ax = plt.subplots()
    fig.set_tight_layout(True)
    configure_axes(ax)

    # add a horizontal line to show the weights if all weights were equal
    ax.axhline(
        y=1.0 / len(selection.timepoints),
        lw=1,
        linestyle="--",
        color=plot_colors["accent1"],
    )

    # plot the boxplot, setting the colors for each part of the plot
    ax.boxplot(
        data,
        notch=False,
        showmeans=True,
        meanline=True,
        sym=None,  # workaround for flierprops being ignored otherwise
        boxprops=dict(color=plot_colors["bright1"]),
        flierprops=dict(markeredgecolor=plot_colors["bright1"], marker="+"),
        medianprops=dict(color=plot_colors["bright2"], linestyle="-"),
        meanprops=dict(color=plot_colors["bright2"], linestyle=":"),
        capprops=dict(color=plot_colors["bright1"], linestyle="-"),
        whiskerprops=dict(color=plot_colors["bright1"], linestyle="--"),
    )

    # add the title and axes labels
    ax.set_title("Regression Weights for {}\n{}".format(label.title(), selection.name))
    ax.set_xlabel("Time Point")
    ax.set_ylabel("Proportional Weight")
    ax.set_xticklabels(selection.timepoints)

    # save the figure and clean up
    pdf.savefig(fig)
    plt.close(fig)


def corr_axes(ax, x, y, score_min, score_max, **kwargs):
    """
    Plot the correlation between scores in two replicates.

    Args:
        ax (matplotlib.axes.Axes): axes object used for plotting

        x (vector): scores for the first replicate

        y (vector): scores for the second replicate

        score_min: minimum score in all plots

        score_max: maximum score in all plots

        **kwargs: arguments to pass to ax.plot for plotting the points (point style, color, etc.)
    """
    configure_axes(ax, xgrid=True)

    # set the axes limits
    ax.set_xlim(left=score_min * 1.05, right=score_max * 1.05)
    ax.set_ylim(bottom=score_min * 1.05, top=score_max * 1.05)

    # plot bolder zero lines
    ax.plot(
        ax.get_xlim(),
        [0.0, 0.0],
        linestyle="-",
        color=plot_colors["accent3"],
        linewidth=1,
    )
    ax.plot(
        [0.0, 0.0],
        ax.get_ylim(),
        linestyle="-",
        color=plot_colors["accent3"],
        linewidth=1,
    )

    # plot the points
    ax.plot(x, y, **kwargs)

    # plot the fit line
    fit = np.polyfit(x, y, 1)
    ax.plot(
        ax.get_xlim(),
        np.poly1d(fit)(ax.get_xlim()),
        linestyle="--",
        color=plot_colors["accent2"],
    )

    # plot the text with regression line info
    r, _ = scipy.stats.pearsonr(x, y)
    ax.text(
        0.05,
        0.95,
        "r^2 = {:.2g}".format(r ** 2),
        horizontalalignment="left",
        verticalalignment="top",
        transform=ax.transAxes,
        fontsize=8,
    )


def hexbin_corr_axes(
    ax, x, y, score_min, score_max, cbar_ax=None, gridsize=50, cmap="YlGnBu_r"
):
    """
    Plot the correlation between scores in two replicates using a hexbin. This 
    function is preferable to corr_axes when plotting a large number of data 
    points.

    Args:
        ax (matplotlib.axes.Axes): axes object used for plotting

        x (vector): scores for the first replicate

        y (vector): scores for the second replicate

        score_min: minimum score in all plots

        score_max: maximum score in all plots

        cbar_ax (matplotlib.axes.Axes): axes object used for plotting the colorbar (optional)

    """
    configure_axes(ax, xgrid=True)

    # set the axes limits
    ax.set_xlim([score_min * 1.05, score_max * 1.05])
    ax.set_ylim([score_min * 1.05, score_max * 1.05])

    # plot bolder zero lines
    ax.plot(
        ax.get_xlim(),
        [0.0, 0.0],
        linestyle="-",
        color=plot_colors["accent3"],
        linewidth=2,
    )
    ax.plot(
        [0.0, 0.0],
        ax.get_ylim(),
        linestyle="-",
        color=plot_colors["accent3"],
        linewidth=2,
    )

    # plot the points
    hb = ax.hexbin(x, y, mincnt=1, cmap=cmap, gridsize=gridsize)
    hb.set_zorder(3)  # on top of the zero lines

    # plot the fit line
    fit = np.polyfit(x, y, 1)
    ax.plot(
        ax.get_xlim(),
        np.poly1d(fit)(ax.get_xlim()),
        linestyle="--",
        linewidth=1,
        color=plot_colors["accent2"],
        zorder=4,
    )

    # plot the text with regression line info
    r, _ = scipy.stats.pearsonr(x, y)
    ax.text(
        0.05,
        0.95,
        "r^2 = {:.3g}".format(r ** 2),
        horizontalalignment="left",
        verticalalignment="top",
        transform=ax.transAxes,
    )

    # plot the colorbar if desired
    if cbar_ax is not None:
        cbar = plt.colorbar(hb, cax=cbar_ax, orientation="horizontal")
        cbar.ax.tick_params(bottom=False, top=False, left=False, right=False)
        cbar.set_label("Count")
        cbar.outline.set_visible(False)

    return hb


def density_ax(ax, ys, xmin, xmax, xlabel, line_params, legend_loc="best"):

    if len(ys) != len(line_params):
        raise ValueError("All y-value sets must have a linestyle")

    configure_axes(ax, xgrid=True)

    d_ys = [KDEUnivariate(y.values) for y in ys]
    [d_y.fit() for d_y in d_ys]

    xs = np.linspace(xmin, xmax, 1000)

    for i in xrange(len(ys)):
        ax.plot(xs, d_ys[i].evaluate(xs), label=ys[i].name, **line_params[i])

    ax.legend(loc=legend_loc, **legend_params)

    ax.set_xlabel(xlabel, **label_params)
    ax.set_ylabel("Density", **label_params)
    ax.tick_params(axis="both", which="major", **tick_params)
