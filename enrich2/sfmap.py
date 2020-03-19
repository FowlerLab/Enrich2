from __future__ import print_function
import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from .plots import plot_cmaps, plot_colors

logger = logging.getLogger(__name__)

#: List of amino acids in row order for sequence-function maps.
AA_LIST = [
    "H",
    "K",
    "R",  # (+)
    "D",
    "E",  # (-)
    "C",
    "M",
    "N",
    "Q",
    "S",
    "T",  # Polar-neutral
    "A",
    "I",
    "L",
    "V",  # Non-polar
    "F",
    "W",
    "Y",  # Aromatic
    "G",
    "P",  # Unique
    "*",
]

#: List of tuples for amino acid physiochemical property groups.
#: Each tuple contains the label string and the corresponding start and end
#: indices in :py:const:`aa_list` (inclusive).
AA_LABEL_GROUPS = [
    ("(+)", 0, 2),
    ("(-)", 3, 4),
    ("Polar-neutral", 5, 10),
    ("Non-polar", 11, 14),
    ("Aromatic", 15, 17),
    ("Unique", 18, 19),
]

#: List of nucleotides in row order for sequence-function maps.
NT_LIST = ["A", "C", "G", "T"]


def parse_aa_list(fname):
    group_labels = list()
    groups = list()
    found_aa = list()
    with open(fname) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            if len(line) == 0:
                continue
            try:
                label, aa = line.split("\t")
            except ValueError:
                raise ValueError(
                    "Unexpected AA list file line format [{}]" "".format("SFMAP.PY")
                )
            aa = aa.split(",")
            aa = [x.strip() for x in aa]
            if any(x not in AA_LIST for x in aa):
                raise ValueError(
                    "Invalid amino acid in AA list file [{}]" "".format("SFMAP.PY")
                )
            else:
                group_labels.append(label)
                groups.append(aa)
                found_aa.extend(aa)
    if any(x not in found_aa for x in AA_LIST):
        raise ValueError(
            "Not all amino acids assigned in AA list file [{}]" "".format("SFMAP.PY")
        )
    if len(AA_LIST) != len(found_aa):
        raise ValueError(
            "Duplicate assignments in AA list file [{}]" "".format("SFMAP.PY")
        )
    pos = -1
    aa_label_groups = list()
    for i, label in enumerate(group_labels):
        aa_label_groups.append((label, pos + 1, pos + len(groups[i])))
        pos = pos + len(groups[i])

    return found_aa, aa_label_groups


def recentered_cmap(cmap, vmin, vmax):
    """
    Rescale the diverging color map *cmap* such that the center color is at 0
    in the data. Returns the rescaled cmap.

    Based on http://stackoverflow.com/a/20528097

    *cmap* is the color map to re-center. Should be a diverging brewer
    palette.

    *vmin* is the minimum value in the data being plotted with the *cmap*
    (must be negative).

    *vmax* is the maximum value in the data being plotted with the *cmap*
    (must be positive).

    """
    # regular index to compute the colors
    reg_index = np.linspace(0.0, 1.0, 257)

    # shifted index to match the data
    centerpoint = 1 - vmax / (vmax + abs(vmin))
    shift_index = np.hstack(
        [
            np.linspace(0.0, centerpoint, 128, endpoint=False),
            np.linspace(centerpoint, 1.0, 129, endpoint=True),
        ]
    )

    # re-map the colors at each index value using a dictionary
    cdict = {"red": [], "green": [], "blue": [], "alpha": []}
    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict["red"].append((si, r, r))
        cdict["green"].append((si, g, g))
        cdict["blue"].append((si, b, b))
        cdict["alpha"].append((si, a, a))

    # save the dictionary as a color map
    newcmap = LinearSegmentedColormap("recentered", cdict)

    return newcmap


def sfmap_axes(
    df,
    ax,
    style,
    wt,
    df_se=None,
    tall=False,
    colors=None,
    missing_color=None,
    vmin=None,
    vmax=None,
    vmax_se=None,
    show_positions=True,
    show_wt=True,
    show_variants=True,
    aa_list=AA_LIST,
    aa_label_groups=AA_LABEL_GROUPS,
):
    """
    Create heatmap of scores or counts for each position in *df* using the axes
    instance *ax*. Returns the |mpl_pcolormesh| (requred to create the color
    bar).

    Amino acids are ordered by their physiochemical properties:

    .. csv-table::
        :header: "Charge" , "Amino Acid" , "Abbreviation" , "Letter"

        Charged (+) , Histidine , His , H
        Charged (+) , Lysine , Lys , K
        Charged (+) , Arginine , Arg , R
        Charged (-) , Aspartate , Asp , D
        Charged (-) , Glutamate , Glu , E
        Polar-neutral , Cysteine , Cys , C
        Polar-neutral , Methionine , Met , M
        Polar-neutral , Asparagine , Asn , N
        Polar-neutral , Glutamine , Gln , Q
        Polar-neutral , Serine , Ser , S
        Polar-neutral , Threonine , Thr , T
        Non-polar , Alanine , Ala , A
        Non-polar , Glycine , Gly , G
        Non-polar , Isoleucine , Ile , I
        Non-polar , Leucine , Leu , L
        Non-polar , Proline , Pro , P
        Non-polar , Valine , Val , V
        Aromatic , Phenylalanine , Phe , F
        Aromatic , Tryptophan , Trp , W
        Aromatic , Tyrosine , Tyr , Y

    *df* is a DataFrame containing scores or counts. 
        Rows are labeled by the integer positions in the wild type sequence,
        numbered in ascending order. Columns are labeled by nucleotide or
        amino acid change.

    *ax* is the axes object used for the heatmap.

    *style* is one of ``'counts'``, ``'logcounts'``, or ``'scores'``.
        ``'counts'`` is used for plotting raw count data for library diversity
        maps.

        ``'logcounts'`` is also used for library diversity maps, but the counts
        are log10 transformed before plotting.

        ``'scores'`` is used for sequence-function maps.

    *wt* is the wild type sequence for the region included in *df*.

    *tall* sets the orientation of the heatmap. 
        ``True`` when the rows of the heatmap correspond to positions, which
        typically results in a plot that is taller than it is wide. Otherwise,
        the columns of the heatmap correspond to positions.

        The actual size of the heatmap is determined by the function that
        handles the returned axes object.

    *colors* is the name of the `matplotlib cmap`_ used.
        If *colors* is ``None``, the ``'counts'`` and ``'logcounts'`` *style*
        use the default sequential cmap and the ``'scores'`` *style* uses the
        default diverging cmap.

        When using the ``'scores'`` style, the cmap is automatically
        re-centered such that the center of the cmap corresponds to 0 in the
        data.

    *missing_color* is the color to use for missing data (``NaN``) in the plot.
        If *missing_color* is ``None``, the default missing color is used.

    *vmin* and *vmax* override the maximum and minimum values in *df*.
        These are used for plotting multiple plots with different ranges using
        the same color bar.

    *show_positions* sets the drawing of the integer position axis tick labels.

    *show_wt* sets the drawing of the wild type sequence axis tick labels.

    *show_variants* sets the drawing of the amino acid or nucleotide change
    axis tick labels.



    *df_se* contains the standard error of each item in *df*

    Standard errors are plotted as diagonal bars on each cell, scaled as a
    percentage of the highest standard error in the plot (or *vmax_se* if
    specified). To scale standard error bars to the highest standard error
    the dataset, calculate the maximum and send it to *vmax_se*. Standard
    errors that are less than 2% of the maximum are not plotted.
    """
    if style not in ("counts", "logcounts", "scores"):
        logger.warning("Invalid style specified for sfmap_axes")
        return None

    # rotate if necessary
    if not tall:
        df = df.transpose()
        if df_se is not None:
            df_se = df_se.transpose()

    # re-order the rows for plotting
    df = df.reindex(index=df.index[::-1])
    if df_se is not None:
        df_se = df_se.reindex(index=df_se.index[::-1])

    # log transform counts if necessary
    if style == "logcounts":
        df = np.log10(df)

    # mask NaN so we can color them separately
    masked = np.ma.array(df, mask=np.isnan(df))
    if df_se is not None:
        masked_se = np.ma.array(df_se, mask=np.isnan(df_se))

    # set colors and min/max data values
    if style == "counts" or style == "logcounts":
        try:
            if colors is not None:
                cmap = plt.get_cmap(colors)
            else:
                cmap = plt.get_cmap(plot_cmaps["sequential"])
        except ValueError:
            logger.warning(
                "Invalid sequential color map choice. " "Falling back to 'BuPu'"
            )
            cmap = plt.get_cmap("BuPu")
        if vmin is None:
            vmin = 0.0
        if vmax is None:
            vmax = masked.max()
    else:  # style == "scores"
        try:
            if colors is not None:
                cmap = plt.get_cmap(colors)
            else:
                cmap = plt.get_cmap(plot_cmaps["diverging"])
        except ValueError:
            logger.warning(
                "Invalid diverging color map choice. " "Falling back to 'RdYlBu_r'"
            )
            cmap = plt.get_cmap("RdYlBu_r")
        if vmin is None:
            vmin = masked.min()
        if vmax is None:
            vmax = masked.max()
        if vmin < 0.0 and vmax > 0.0:
            cmap = recentered_cmap(cmap, vmin, vmax)
        else:
            logger.warning("Unexpected range. Not recentering color map.")
    if missing_color is not None:
        cmap.set_bad(missing_color, 1.0)
    else:
        cmap.set_bad(plot_colors["missing"], 1.0)

    # draw the heatmap and get rid of extra whitespace
    mesh = ax.pcolormesh(masked, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xlim(0, len(df.columns))
    ax.set_ylim(0, len(df.index))

    # add marks on wild type positions
    wt = list(wt)
    if tall:
        wt_xy = zip((list(df.columns).index(x) for x in wt), reversed(xrange(len(wt))))
    else:
        wt_xy = zip(xrange(len(wt)), (list(df.index).index(x) for x in wt))
    for x, y in wt_xy:
        ax.add_patch(
            Circle(
                (x + 0.5, y + 0.5),
                0.1666,
                fill=True,
                facecolor="black",
                edgecolor="none",
                alpha=0.5,
            )
        )

    # add diagonal SE marks
    if df_se is not None:
        if vmax_se is None:
            vmax_se = masked_se.max()
        # rescale the SE's onto 0 .. 0.98
        # rescaling onto 0 .. 1.0 causes the corners to look funny
        masked_se = masked_se / vmax_se * 0.98
        for x in xrange(len(df.index)):
            for y in xrange(len(df.columns)):
                value = masked_se[x, y]
                if value and value >= 0.02:  # not masked, above threshold
                    corner_dist = (1.0 - value) / 2.0
                    diag = Line2D(
                        [y + corner_dist, y + 1 - corner_dist],
                        [x + corner_dist, x + 1 - corner_dist],
                        transform=ax.transData,
                        color="grey",
                    )
                    ax.add_line(diag)

    # calculate tick label positions
    ypos = np.arange(len(df.index)) + 0.5
    xpos = np.arange(len(df.columns)) + 0.5
    ax.set_yticks(ypos)
    ax.set_xticks(xpos)

    # create the secondary axis for the wild type sequence
    if tall:
        wt_ax = ax.twinx()
        wt_ax.set_yticks(ypos)
    else:
        wt_ax = ax.twiny()
        wt_ax.set_xticks(xpos)

    # add the integer position axis tick labels
    if show_positions:
        # add position axis label
        if tall:
            ax.set_ylabel("Position")
        else:
            ax.set_xlabel("Position")
        if tall:
            ax.set_yticklabels(list(df.index))
        else:
            ax.set_xticklabels(list(df.columns), rotation=90)

    # add the variant information
    # note this is added as text within the plot not tick labels
    if show_variants:
        if tall:
            # turn off existing xticklabels
            plt.setp(ax.get_xticklabels(), visible=False)
            # add fake xlabels as plot text
            for i, x in enumerate(list(df.columns)):
                ax.text(
                    i + 0.5,
                    len(df.index) + 0.4,
                    x,
                    horizontalalignment="center",
                    verticalalignment="center",
                    transform=ax.transData,
                )
            # add additional amino acid property information
            if list(df.columns) == aa_list:
                ax.set_ylim(0, len(df.index) + 2.0)
                wt_ax.set_ylim(0, len(df.index) + 2.0)
                for label, start, end in aa_label_groups:
                    if len(label) == 0:
                        continue
                    ax.text(
                        (end - start + 1) / 2.0 + start,
                        len(df.index) + 1.4,
                        label,
                        horizontalalignment="center",
                        verticalalignment="center",
                        transform=ax.transData,
                    )
                    bar = Line2D(
                        [start + 0.125, end + 1 - 0.125],
                        [len(df.index) + 0.9, len(df.index) + 0.9],
                        transform=ax.transData,
                        color="grey",
                    )
                    bar.set_clip_on(False)
                    ax.add_line(bar)
            # no amino acid properties, just add padding
            else:
                ax.set_ylim(0, len(df.index) + 1.0)
                wt_ax.set_ylim(0, len(df.index) + 1.0)
        else:
            # turn off existing yticklabels
            plt.setp(ax.get_yticklabels(), visible=False)
            # add fake ylabels as plot text
            for i, x in enumerate(list(df.index)):
                ax.text(
                    -0.4,
                    i + 0.5,
                    x,
                    horizontalalignment="center",
                    verticalalignment="center",
                    transform=ax.transData,
                )
            # add additional amino acid property information
            if list(df.index)[::-1] == aa_list:
                ax.set_xlim(-2, len(df.columns))
                wt_ax.set_xlim(-2, len(df.columns))
                for label, start, end in aa_label_groups:
                    if len(label) == 0:
                        continue
                    ax.text(
                        -1,
                        len(df.index) - ((end - start + 1) / 2.0 + start),
                        label,
                        rotation="vertical",
                        verticalalignment="center",
                        horizontalalignment="right",
                        transform=ax.transData,
                    )
                    bar = Line2D(
                        [-0.9, -0.9],
                        [
                            len(df.index) - (start + 0.125),
                            len(df.index) - (end + 1 - 0.125),
                        ],
                        transform=ax.transData,
                        color="grey",
                    )
                    bar.set_clip_on(False)
                    ax.add_line(bar)
            # no amino acid properties, just add padding
            else:
                ax.set_xlim(-1, len(df.columns))
                wt_ax.set_xlim(-1, len(df.columns))

    # add the wild type sequence labels
    if show_wt:
        if tall:
            wt_ax.set_yticklabels(wt[::-1])
            wt_ax.set_ylabel("Wild Type Sequence")
        else:
            wt_ax.set_xticklabels(wt)
            wt_ax.set_xlabel("Wild Type Sequence")

    # hide the tick marks
    # needs to be done after the ax.twiny() or ax.twinx() call for some reason
    ax.tick_params(bottom=False, top=False, left=False, right=False)
    ax.set_frame_on(False)
    wt_ax.tick_params(bottom=False, top=False, left=False, right=False)
    wt_ax.set_frame_on(False)

    return mesh


def sfmap_plot(
    df, pdf, style, wt, dimensions, df_se=None, title=None, show_colorbar=True, **kwargs
):
    """
    Create heatmap of scores or counts for each position in *df* and save it to
    the file *pdf*.

    *df* is a DataFrame containing scores or counts.
        Rows are labeled by the integer positions in the wild type sequence,
        numbered in ascending order. Columns are labeled by nucleotide or
        amino acid change.

    *pdf* is the open PdfPages file object.

    *style* is one of ``'counts'``, ``'logcounts'``, or ``'scores'``.
        ``'counts'`` is used for plotting raw count data for library diversity
        maps.

        ``'logcounts'`` is also used for library diversity maps, but the counts
        are log10 transformed before plotting.

        ``'scores'`` is used for sequence-function maps.

    *wt* is the wild type sequence for the region included in *df*.

    *dimensions* is either ``'tall'`` or ``'wide'``, or a two-tuple of floats
    giving the figure size in inches.

        ``'tall'`` will create a figure 8.5 inches wide that is tall enough to
        make all cells in the heatmap roughly square.

        ``'wide'`` will create a figure 8.5 inches tall that is wide enough to
        make all cells in the heatmap roughly square.

        The two-tuple of floats specifies the exact figure size (including
        title and color bar, if requested) in inches. The long dimension is
        determined based on the ratio between width and height.

        .. note:: Specifying a figure size may cause the axis labels to not \
        display properly.

    *title* is the figure title. It can include ``'\\n'`` characters for a
    multi-line title.

    *show_colorbar* is ``True`` to display a color bar, or ``False`` to not
    display a color bar. If the figure is in tall mode, the color bar appears
    along the bottom of the plot, otherwise it appears in the right side.

    *kwargs* is the dictionary of additional optional or required parameters to
    :py:func:`sfmap_axes`.

    """
    if style not in ("counts", "logcounts", "scores"):
        logger.warning("Invalid style specified for sfmap_axes.")
        return None

    # dimensions set explicitly
    if len(dimensions) == 2:
        try:
            dimx, dimy = [float(n) for n in dimensions]
        except ValueError:
            logger.warning("Invalid sequence-function map dimensions.")
            return None
    else:
        # scaling based on mode
        if dimensions == "tall":
            # headings + pos_ax + wt_ax
            colx = len(df.columns)
            # rows + cbar_ax + aa_headings + title_ax
            coly = len(df.index) + 2 + 2 + 4
        elif dimensions == "wide":
            # rows + aa_headings + cbar_ax
            colx = len(df.index) + 2 + 2
            # headings + pos_ax + wt_ax + title_ax
            coly = len(df.columns) + 2
        else:
            logger.warning("Invalid sequence-function map scaling mode.")
            return None
        base_scale = 8.5 / 21
        dimx = colx * base_scale
        dimy = coly * base_scale

    tall = dimy > dimx

    # create the figure
    fig = plt.figure()
    fig.set_size_inches(dimx, dimy)

    # set up the grid for the figure
    if tall:
        if title is not None and show_colorbar:
            grid = GridSpec(
                3,
                1,
                height_ratios=[1, len(df.index) + 2, 1],
                hspace=1.0 / len(df.index),
            )
            title_ax = plt.subplot(grid[0])
            mesh_ax = plt.subplot(grid[1])
            cbar_ax = plt.subplot(grid[2])
        elif title is not None:  # no colorbar
            grid = GridSpec(2, 1, height_ratios=[1, len(df.index)])
            title_ax = plt.subplot(grid[0])
            mesh_ax = plt.subplot(grid[1])
            cbar_ax = None
        elif show_colorbar:  # no title
            grid = GridSpec(2, 1, height_ratios=[len(df.index), 1])
            title_ax = None
            mesh_ax = plt.subplot(grid[0])
            cbar_ax = plt.subplot(grid[1])
        else:  # only the heatmap
            grid = GridSpec(1, 1)
            title_ax = None
            mesh_ax = plt.subplot(grid[0])
            cbar_ax = None
    else:
        if title is not None and show_colorbar:
            grid = GridSpec(
                2,
                2,
                height_ratios=[1, len(df.columns)],
                width_ratios=[len(df.index), 1],
                wspace=1.0 / len(df.index),
            )
            title_ax = plt.subplot(grid[0, :])
            mesh_ax = plt.subplot(grid[1, 0])
            cbar_ax = plt.subplot(grid[1, 1])
        elif title is not None:  # no colorbar
            grid = GridSpec(2, 1, height_ratios=[1, len(df.index)])
            title_ax = plt.subplot(grid[0])
            mesh_ax = plt.subplot(grid[1])
            cbar_ax = None
        elif show_colorbar:  # no title
            grid = GridSpec(1, 2, width_ratios=[len(df.columns), 1])
            title_ax = None
            mesh_ax = plt.subplot(grid[0])
            cbar_ax = plt.subplot(grid[1])
        else:  # only the heatmap
            grid = GridSpec(1, 1)
            title_ax = None
            mesh_ax = plt.subplot(grid[0])
            cbar_ax = None

    # create the heatmap
    mesh = sfmap_axes(df, mesh_ax, style, wt, df_se=df_se, tall=tall, **kwargs)

    # add the colorbar axis
    if cbar_ax is not None:
        if tall:
            cbar = fig.colorbar(mesh, cax=cbar_ax, orientation="horizontal")
        else:
            cbar = fig.colorbar(mesh, cax=cbar_ax)
        if style == "logcounts":
            cbar.set_label("log10(Count)")
        elif style == "counts":
            cbar.set_label("Count")
        else:  # style == "scores"
            cbar.set_label("Score")
        cbar.outline.set_visible(False)
        cbar.ax.tick_params(bottom=False, top=False, left=False, right=False)

    # add the title
    if title_ax is not None:
        title_ax.text(s=title, x=0.5, y=1.0, ha="center", va="top", size="large")
        title_ax.axis("off")

    # save the heatmap
    pdf.savefig(fig, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
