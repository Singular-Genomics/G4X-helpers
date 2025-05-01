import math
from typing import List, Literal, Optional

import matplotlib.pyplot as plt
import matplotlib_inline.backend_inline as mpl_inline
import numpy as np
from cmcrameri import cm  # noqa: F401
from matplotlib import gridspec
from matplotlib.transforms import ScaledTranslation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

style_dict = {
    'axes.grid': False,
    'grid.alpha': 0.35,
    'grid.color': '#ABB2BF',
    'grid.linewidth': 0.5,
    'text.color': '#ABB2BF',
    #'font.family': 'Lato',
    'font.size': 12,
    'figure.facecolor': '#282C34',
    'figure.figsize': [6, 6],
    'figure.dpi': 100,
    'figure.titleweight': 'bold',
    'axes.edgecolor': '#ABB2BF',
    'axes.facecolor': '#21252B',
    'axes.labelcolor': '#ABB2BF',
    'axes.labelpad': 8,
    'axes.labelsize': 14,
    'axes.labelweight': 'bold',
    'axes.titlecolor': '#ABB2BF',
    'axes.titlepad': 12,
    'axes.titlesize': 20,
    'axes.titleweight': 'bold',
    'axes.linewidth': 0.5,
    'axes.xmargin': 0.1,
    'axes.ymargin': 0.1,
    'xtick.color': '#ABB2BF',
    'ytick.color': '#ABB2BF',
    'xtick.major.size': 3.5,
    'ytick.major.size': 3.5,
    'xtick.minor.size': 2,
    'ytick.minor.size': 2,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'ytick.major.width': 0.5,
    'xtick.major.width': 0.5,
    'ytick.minor.width': 0.35,
    'xtick.minor.width': 0.35,
    'ytick.minor.visible': True,
    'xtick.minor.visible': True,
}

set25 = ['#72FFAB', '#A16CFD', '#FF7043', '#008FFF', '#D32F2F',
         '#7CB342', '#7F34BE', '#FFCA28', '#0C8668', '#FB4695',
         '#005EE1', '#90A4AE', '#28EDED', '#A17B64', '#FFFF58',
         '#BC29AE', '#006D8F', '#FFBAFF', '#FFD091', '#5C6BC0',
         '#F490B2', '#C6E1A6']

palD3 = ['#1F77B4FF', '#FF7F0EFF', '#2CA02CFF', '#D62728FF', 
         '#9467BDFF', '#8C564BFF', '#E377C2FF', '#7F7F7FFF', 
         '#BCBD22FF', '#17BECFFF', '#AEC7E8FF', '#FFBB78FF', 
         '#98DF8AFF', '#FF9896FF', '#C5B0D5FF', '#C49C94FF', 
         '#F7B6D2FF', '#C7C7C7FF', '#DBDB8DFF', '#9EDAE5FF']

plt.rcParams.update(style_dict)

retina = True
dpi = 100
inline_dpi = None

if inline_dpi is None:
    inline_dpi = dpi

# if inline_transparent == True:
#     inline_facecolor = "none"
# elif inline_transparent == False:
#     inline_facecolor = fig_facecolor

if retina:
    inl_format = 'retina'
    inline_dpi = inline_dpi * 2
else:
    inl_format = 'png'

mpl_inline.set_matplotlib_formats(inl_format, facecolor='none', bbox_inches='tight', dpi=inline_dpi)


def montage_plot(
    n_plots: int = None,
    n_rows: int = None,
    n_cols: int = None,
    design: Optional[List[List[int]]] = None,
    w_ratios: list = None,
    h_ratios: list = None,
    wspace: float = None,
    hspace: float = None,
    show_labels: bool = False,
    label_loc: str = 'ul',
    figure=None,
    title: str = None,
    layout: str = 'constrained',
    figsize: list = [6, 6],
    panel_size: float = None,
    return_fig=False,
):
    """
    Complex subplot layout from a list-of-list design.

    Example:
    --------
    design = [
        [0, 0, 3],
        [1, 2, 3],
        [4, 4, 4]
    ]

    placing -1 into the design will create a blank space
    """

    def _build_design(n_plots, n_rows, n_cols):
        if n_plots == 1:
            design = [[0]]

        else:
            # Create a default design based on the calculated layout
            design = []
            j = 0
            for i in range(n_rows):
                row = []
                for k in range(n_cols):
                    row.append(j if j < n_plots else -1)
                    j += 1
                design.append(row)

        return design

    if not figure:
        figure = plt.figure(figsize=figsize, layout=layout)

    figure.suptitle(title)

    if design is None and n_plots is None and n_rows is None and n_cols is None:
        n_plots = 1
        # raise ValueError('Please specify either a design or the number of plots via `n_plots` and/or `n_rows|n_cols`.')

    # If design is provided, use it to determine rows and columns
    if design is not None:
        if not all(isinstance(row, list) for row in design):
            raise ValueError('Design should be a list of lists.')
        n_rows = len(design)
        n_cols = len(design[0])

    elif design is None and n_plots is not None:
        # Derive n_rows and n_cols from n_plots with a preference for balanced dimensions
        if n_cols is None and n_rows is None:
            n_rows = 1 if n_plots < 4 else 2 if n_plots <= 10 else int(math.floor(math.sqrt(n_plots)))
            n_cols = math.ceil(n_plots / n_rows)
        elif n_rows is not None and n_cols is None:
            n_cols = math.ceil(n_plots / n_rows)
        elif n_cols is not None and n_rows is None:
            n_rows = math.ceil(n_plots / n_cols)

        design = _build_design(n_plots=n_plots, n_rows=n_rows, n_cols=n_cols)

    elif n_plots is None and n_rows is not None and n_cols is not None:
        n_plots = n_cols * n_rows

        design = _build_design(n_plots=n_plots, n_rows=n_rows, n_cols=n_cols)

    if panel_size is not None:
        figure.set_size_inches(n_cols * panel_size, n_cols * panel_size)

    grid = gridspec.GridSpec(
        nrows=n_rows,
        ncols=n_cols,
        figure=figure,
        width_ratios=w_ratios,
        height_ratios=h_ratios,
        wspace=wspace,
        hspace=hspace,
    )

    arr = np.array(design)
    ax_labels = np.unique(arr)
    ax_labels = ax_labels[ax_labels != -1]

    plot_layout = []

    for plot in ax_labels:
        if plot != -1:
            ax_row = np.where(arr == plot)[0]
            ax_col = np.where(arr == plot)[1]

            plot_grids = grid[min(ax_row) : max(ax_row) + 1, min(ax_col) : max(ax_col) + 1]
            plot_layout.append(plt.subplot(plot_grids))

    if show_labels is not False:
        if show_labels == 'abc':
            ax_labels = [chr(65 + num) for num in ax_labels]

        for i, ax in enumerate(plot_layout):
            place_text(ax, ax_labels[i], offset=15, size=20, loc=label_loc)

    if len(plot_layout) == 1:
        plot_layout = plot_layout[0]

    if return_fig:
        return figure, plot_layout
    else:
        return plot_layout


def downsample_mask(mask, mask_length=None, ax=None, display_size=None, downsample='auto'):
    def get_ax_size(ax):
        fig = ax.get_figure()

        bbox = ax.get_position()
        fig_width, fig_height = fig.get_size_inches()
        width_in = bbox.width * fig_width
        height_in = bbox.height * fig_height

        ax_dim = (width_in, height_in)
        return ax_dim

    if downsample is None:
        return mask

    elif isinstance(downsample, int):
        downsample_fct = int(downsample)

    elif downsample.startswith('auto'):
        if downsample == 'auto':
            target_ppi = 200
        elif downsample.startswith('auto_'):
            target_ppi = int(downsample.split('_')[-1])

        if ax is None and display_size is None:
            raise ValueError("Either ax or display_size must be provided when downsample is 'auto'")
        elif ax is None and display_size:
            panel_size = display_size
        elif ax and display_size is None:
            ax_dim = get_ax_size(ax)
            panel_size = np.array(ax_dim).max()

        if mask_length is None:
            mask_length = np.array(mask.shape).max()
        elif mask_length is not None:
            mask_length = mask_length

        downsample_fct = int((mask_length / panel_size) / target_ppi)
        # print(f'Downsample factor: {downsample_fct}')

    else:
        raise ValueError("downsample must be None, int or 'auto'")

    if downsample_fct > 0:
        sampled_mask = mask[::downsample_fct, ::downsample_fct].copy()
    else:
        sampled_mask = mask.copy()

    return sampled_mask


def place_text(ax, label, offset=20, loc='tl', size=20, weight='bold', color=None, **kwargs):
    if color is None:
        color = plt.rcParams['text.color']

    dpi = plt.gcf().dpi

    if loc == 'tl':
        loc_x, loc_y = 0, 1
        dx, dy = offset / dpi, -offset / dpi
        ha, va = 'left', 'top'
    elif loc == 'tr':
        loc_x, loc_y = 1, 1
        dx, dy = -offset / dpi, -offset / dpi
        ha, va = 'right', 'top'
    elif loc == 'c':
        loc_x, loc_y = 0.5, 0.5
        dx, dy = 0 / dpi, 0 / dpi
        ha, va = 'center', 'center'
    elif loc == 'br':
        loc_x, loc_y = 1, 0
        dx, dy = -offset / dpi, offset / dpi
        ha, va = 'right', 'bottom'
    elif loc == 'bl':
        loc_x, loc_y = 0, 0
        dx, dy = offset / dpi, offset / dpi
        ha, va = 'left', 'bottom'
    else:
        raise ValueError("Invalid location. Choose from ['tl', 'tr', 'c', 'br', 'bl'].")

    offset = ScaledTranslation(dx, dy, ax.figure.dpi_scale_trans)
    transformed = ax.transAxes + offset

    ax.text(
        loc_x,
        loc_y,
        label,
        horizontalalignment=ha,
        verticalalignment=va,
        transform=transformed,
        fontsize=size,
        weight=weight,
        color=color,
        zorder=100,
        **kwargs,
    )


_unit = Literal['micron', 'px']
_loc = Literal[
    'upper left',
    'upper center',
    'upper right',
    'center left',
    'center',
    'center right',
    'lower left',
    'lower center',
    'lower right',
]

def _round_scale(scale_size):
    # list of (upper‐bound, rounding increment)
    thresholds = [
        (25,   5),
        (100, 10),
        (500, 50),
        (1000,100),
        (10000,500),
    ]
    # find the first threshold your size falls under
    for bound, inc in thresholds:
        if scale_size < bound:
            return round(scale_size / inc) * inc
    # fallback for anything ≥ 10000
    return round(scale_size / 1000) * 1000

def _scale_bar(
    ax,
    length: float,
    unit: _unit = 'micron',
    lw: float = 0.5,
    theme: str = 'dark',
    contrast: float = None,
    fontsize: float = 10,
    fontweight: str = 'bold',
    loc: _loc = 'lower left',
    borderpad: float = 1,
    background: bool = False,
    background_border: float = 0,
    background_color: str = None,
    background_alpha: float = 0.5,
    **kwargs,
):
    import matplotlib.font_manager as fm
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

    if theme == 'dark':
        if background_color is None:
            background_color = '#21252B'
        if contrast is None:
            contrast = 0.6
        color = str(contrast)

    elif theme == 'light':
        if background_color is None:
            background_color = '#ABB2BF'
        if contrast is None:
            contrast = 0.9
        color = str(1 - contrast)

    if unit == 'micron':
        len_microns = length
        len_px = len_microns / 0.3125
        unit_txt = 'µm'
    elif unit == 'px':
        len_px = length
        unit_txt = 'px'

    # Define the scale bar parameters
    scalebar_length = len_px  # Length of the scale bar in data units
    scalebar_text = f'{length} {unit_txt}'  # Text displayed below the scale bar

    pad = kwargs.pop('pad', 0.5)
    sep = kwargs.pop('sep', 3)
    label_top = kwargs.pop('label_top', True)
    size_vertical = kwargs.pop('size_vertical', 0)

    # Create the scale bar using AnchoredSizeBar
    scale_bar = AnchoredSizeBar(
        transform=ax.transData,
        size=scalebar_length,
        label=scalebar_text,
        loc=loc,
        borderpad=borderpad,
        sep=sep,
        pad=pad,
        label_top=label_top,
        size_vertical=size_vertical,
        fill_bar=True,
        fontproperties=fm.FontProperties(size=fontsize, weight=fontweight),
        color=color,
        frameon=background,
        **kwargs,
    )

    scale_bar.patch.set_facecolor(background_color)
    scale_bar.patch.set_edgecolor(color)
    scale_bar.patch.set_linewidth(background_border)
    scale_bar.patch.set_alpha(background_alpha)
    scale_bar.size_bar.get_children()[0].set_linewidth(lw)

    ax.add_artist(scale_bar)

def _show_image(image, roi, ax, vmin=None, vmax=None, cmap=None, scale_bar=True):
    plot_image = downsample_mask(image, mask_length=None, ax=ax, downsample='auto')
    im = ax.imshow(plot_image, vmin=vmin, vmax=vmax, cmap=cmap, extent=roi.extent_array, origin='lower')
    roi._limit_ax(ax)

    arr = plot_image / np.max(plot_image)

    y, x = arr.shape[0], arr.shape[1]
    ym, xm = int(y / 5), int(x / 5)

    scale_corner_value = np.mean(arr[0:ym, 0:xm])

    if scale_bar:
        um = roi.width * 0.3125
        rounded_scale = _round_scale(um / 5)
        theme = 'light' if scale_corner_value > 0.75 else 'dark'
        _scale_bar(
            ax, length=rounded_scale, unit='micron', theme=theme, background=True, lw=1, background_alpha=0.5
        )

    return im

def _add_colorbar(ax, cax, loc, title=''):
    if loc == 'right':
        loc = 'upper left'
        width = 0.1
        height = '100%'
        anchor = (1.025, 0, 1, 1)
        orientation = 'vertical'

    if loc == 'bottom':
        loc = 'upper center'
        width = '100%'
        height = 0.1
        anchor = (0, 0, 1, -0.025)
        orientation = 'horizontal'

    axins = inset_axes(
        ax,
        width=width,
        height=height,
        loc=loc,
        bbox_to_anchor=anchor,
        bbox_transform=ax.transAxes,
        borderpad=0,
    )
    cbar = ax.figure.colorbar(cax, cax=axins, label='Test', orientation=orientation)
    cbar.set_label(title, fontsize=12, fontweight='bold')


_legend_pos = Literal['right', 'bottom', 'inset_br', 'inset_bl', 'inset_tr', 'inset_tl']
def _add_legend(
    ax,
    categories,
    title: str = None,
    palette=list,
    loc: _legend_pos = 'right',
    ncol: str = None,
    text_size=12,
    labelspacing=0.75,
    handletextpad=0,
    borderpad=0,
    columnspacing=1,
    handlelength=2,
    handleheight=1,
):
    def calculate_ncol(in_size):
        if in_size <= 2:
            out = 8
        elif in_size <= 4:
            out = 6
        elif in_size <= 6:
            out = 5
        elif in_size <= 8:
            out = 4
        elif in_size <= 10:
            out = 3
        else:
            out = 2

        return out

    color_dict = {cat: palette[i] for i, cat in enumerate(categories)}

    if loc in ['right', 'inset_tr', 'inset_tl', 'inset_br', 'inset_bl']:
        if ncol is None:
            ncol = 1 if len(categories) <= 14 else 2 if len(categories) <= 30 else 3
        anchor = {'right': (1, 1), 'inset_tr': (1, 1), 'inset_tl': (0, 1), 'inset_br': (1, 0), 'inset_bl': (0, 0)}[loc]
        loc = {'right': 'upper left', 'inset_tr': 'upper right', 'inset_tl': 'upper left', 'inset_br': 'lower right', 'inset_bl': 'lower left'}[loc]

    elif loc == 'bottom':
        if ncol is None:
            longest_cat = np.max([len(category) for category in color_dict])
            ncol = calculate_ncol(longest_cat)
        loc = 'upper center'
        anchor = (0.5, 0)

    for label in categories:
        ax.scatter([], [], c=color_dict[label], label=label)

        legend = ax.legend(
            title=title,
            fontsize=text_size,
            frameon=False,
            loc=loc,
            bbox_to_anchor=anchor,
            bbox_transform=ax.transAxes,
            ncol=ncol,
            labelspacing=labelspacing,  # Space between legend entries
            handletextpad=handletextpad,  # Space between handle and text
            borderpad=borderpad,  # Space between legend border and content
            columnspacing=columnspacing,  # Space between columns if ncol > 1
            handlelength=handlelength,  # Length of legend handles
            handleheight=handleheight,  # Height of legend handles
        )

    plt.setp(legend.get_title(), color=plt.rcParams['axes.labelcolor'], fontweight='bold', size=14)
    legend._legend_box.align = 'left'

