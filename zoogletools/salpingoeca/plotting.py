from difflib import get_close_matches
from enum import Enum, StrEnum

import arcadia_pycolor as apc
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from zoogletools.plotting import create_save_fig_config

apc.plotly.setup()


class SalpingoecaStage(StrEnum):
    """Enum for Salpingoeca expression stages"""

    SLOW_SWIMMER = "Slow Swimmer"
    FAST_SWIMMER = "Fast Swimmer"
    ROSETTE = "Rosette"
    THECATE = "Thecate"


SALPINGOECA_EXPRESSION_STAGES = [
    SalpingoecaStage.SLOW_SWIMMER,
    SalpingoecaStage.FAST_SWIMMER,
    SalpingoecaStage.ROSETTE,
    SalpingoecaStage.THECATE,
]
SALPINGOECA_CELLTYPE_COLORMAP = {
    SalpingoecaStage.SLOW_SWIMMER: apc.seaweed.hex_code,
    SalpingoecaStage.FAST_SWIMMER: apc.crow.hex_code,
    SalpingoecaStage.ROSETTE: apc.aster.hex_code,
    SalpingoecaStage.THECATE: apc.dragon.hex_code,
}

DEFAULT_SIGNIFICANCE_COLORS = [
    "white",
    apc.lapis,
    apc.aegean,
    apc.vital,
    apc.dove,
]


# This colorscale is used to map the significance values
# within different ranges to different colors.
# In plotly colormaps, overlapping values result in discrete color boundaries.
def _create_significance_colorscale(colors: list[str]) -> list[list[str]]:
    """Create a colorscale for the significance heatmap.

    Args:
        colors: List of colors to use for the significance heatmap
    """
    return [
        [0, colors[0]],  # These two entries show the background color,
        [0.01, colors[0]],  # used for the diagonal of the heatmap.
        [0.01, colors[1]],  # p < 0.001 lower bound
        [0.25, colors[1]],  # p < 0.001 upper bound
        [0.25, colors[2]],  # p < 0.01 lower bound
        [0.5, colors[2]],  # p < 0.01 upper bound
        [0.5, colors[3]],  # p < 0.05 lower bound
        [0.75, colors[3]],  # p < 0.05 upper bound
        [0.75, colors[4]],  # Not significant
        [1, colors[4]],  # Not significant
    ]


SIGNIFICANCE_COLORSCALE = _create_significance_colorscale(DEFAULT_SIGNIFICANCE_COLORS)

# These values are used to map the significance values
# to the colorbar.
SIGNIFICANCE_COLOR_VALUE_MAP = {
    "diagonal": 0,
    "p < 0.001": 0.125,
    "p < 0.01": 0.375,
    "p < 0.05": 0.625,
    "n.s.": 0.875,
}

SIGNIFICANCE_COLORBAR_PARAMS = dict(
    tickmode="array",
    ticktext=["p < 0.001", "p < 0.01", "p < 0.05", "n.s."],
    tickvals=[0.125, 0.375, 0.625, 0.875],
    thickness=25,
    outlinewidth=0,
    len=1.11,
)


def create_salpingoeca_id_mapping(
    results_filepath: str,
    diffex_filepath: str,
):
    """Create a mapping between HGNC gene symbols and Salpingoeca UniProt IDs.

    Args:
        results_filepath: Path to Zoogle results TSV file containing HGNC gene symbols
            and Salpingoeca protein IDs
        diffex_filepath: Path to Leon et al differential expression TSV file containing
            Uniprot IDs for expressed proteins

    Returns:
        dict: Mapping from HGNC gene symbols to their corresponding
            UniProt IDs for the proteins in the Leon et al TSV file.
    """
    zoogle_results = pd.read_csv(
        results_filepath, sep="\t", usecols=["hgnc_gene_symbol", "nonref_protein"]
    )
    salpingoeca_expression = pd.read_csv(diffex_filepath, sep="\t", usecols=["Uniprot ID"])

    salpingoeca_map = pd.merge(
        zoogle_results,
        salpingoeca_expression[["Uniprot ID"]].drop_duplicates(),
        left_on="nonref_protein",
        right_on="Uniprot ID",
        how="inner",
    )[["hgnc_gene_symbol", "nonref_protein"]].drop_duplicates()

    return dict(
        zip(salpingoeca_map["hgnc_gene_symbol"], salpingoeca_map["nonref_protein"], strict=True)
    )


def _convert_expression_rows_to_boxplot_df(expression_rows: pd.DataFrame) -> pd.DataFrame:
    """Convert expression data rows into a format suitable for boxplot visualization.

    Args:
        expression_rows: DataFrame containing expression data with condition and replicate columns

    Returns:
        DataFrame with columns 'Cell type' and 'TPM' suitable for boxplot creation
    """
    cell_types = []
    tpm_values = []

    # Replicate columns for each condition, as found in the Leon et al differential expression file.
    condition1_replicates = ["Cond 1a", "Cond 1b", "Cond 1c"]
    condition2_replicates = ["Cond 2a", "Cond 2b", "Cond 2c"]

    for condition in SALPINGOECA_EXPRESSION_STAGES:
        condition1_df = expression_rows[expression_rows["Cond 1"] == condition]
        condition2_df = expression_rows[expression_rows["Cond 2"] == condition]

        if not condition1_df.empty:
            values = condition1_df.iloc[0][condition1_replicates].values
            cell_types.extend([condition] * len(values))
            tpm_values.extend(values)
        elif not condition2_df.empty:
            values = condition2_df.iloc[0][condition2_replicates].values
            cell_types.extend([condition] * len(values))
            tpm_values.extend(values)

    return pd.DataFrame({"Cell type": cell_types, "TPM": tpm_values})


def _create_boxplot_traces(
    boxplot_df: pd.DataFrame,
    showlegend: bool = False,
    colormap: dict[SalpingoecaStage, str] = SALPINGOECA_CELLTYPE_COLORMAP,
) -> list[go.Box]:
    """Create boxplot traces for each cell type.

    Args:
        boxplot_df: DataFrame containing cell type and TPM values
        showlegend: Whether to show the legend for the boxplot traces
    """
    traces = []
    for stage in SALPINGOECA_EXPRESSION_STAGES:
        traces.append(
            go.Box(
                x=boxplot_df[boxplot_df["Cell type"] == stage]["Cell type"],
                y=boxplot_df[boxplot_df["Cell type"] == stage]["TPM"],
                line=dict(color=colormap[stage]),
                boxpoints="all",
                showlegend=showlegend,
            )
        )
    return traces


def _load_expression_rows(
    expression_filepath: str, salpingoeca_map_ids: dict, symbol: str
) -> pd.DataFrame:
    """Load expression rows for a given symbol.
    Each row in the differential expression TSV file contains the expression data for a single
        pairwise comparison of a given gene's expression between two conditions.

    Args:
        expression_filepath: Path to the differential expression TSV file
        salpingoeca_map_ids: Dictionary mapping HGNC gene symbols to Salpingoeca protein IDs
        symbol: HGNC gene symbol to search for

    Returns:
        pd.DataFrame: Expression rows for the given symbol
    """
    expression_df = pd.read_csv(expression_filepath, sep="\t")
    try:
        expression_rows = expression_df[expression_df["Uniprot ID"] == salpingoeca_map_ids[symbol]]
    except KeyError as e:
        possible_matches = [i for i in salpingoeca_map_ids.keys() if pd.notna(i)]
        close_matches = get_close_matches(symbol, possible_matches, n=3, cutoff=0.6)
        error_msg = f"Symbol {symbol} not found in {expression_filepath}."
        if close_matches:
            error_msg += f"\nDid you mean one of these?: {', '.join(close_matches)}"
        raise ValueError(error_msg) from e
    return expression_rows


def plot_expression_boxplot(
    symbol: str,
    expression_filepath: pd.DataFrame,
    salpingoeca_map_ids: dict,
    width: int = 800,
    height: int = 600,
    output_image_filepath: str = None,
    output_html_filepath: str = None,
    colormap: dict[SalpingoecaStage, str] = SALPINGOECA_CELLTYPE_COLORMAP,
):
    """Plot the expression of a given symbol in the form of a boxplot.

    Args:
        symbol: HGNC gene symbol to plot
        expression_filepath: Path to the differential expression TSV file
        salpingoeca_map_ids: Dictionary mapping HGNC gene symbols to Salpingoeca protein IDs
        width: Width of the figure
        height: Height of the figure
        output_image_filepath: Path to save the output image file
        output_html_filepath: Path to save the output HTML file
        colormap: Dictionary mapping Salpingoeca stages to colors

    Returns:
        fig: Plotly figure
    """
    expression_rows = _load_expression_rows(expression_filepath, salpingoeca_map_ids, symbol)
    boxplot_df = _convert_expression_rows_to_boxplot_df(expression_rows)

    traces = _create_boxplot_traces(boxplot_df, showlegend=False, colormap=colormap)
    fig = go.Figure(data=traces)

    fig.update_layout(
        title=f"{symbol} ({salpingoeca_map_ids[symbol]}) Expression in S. rosetta",
        width=width,
        height=height,
        showlegend=False,
        margin=dict(
            t=15,
        ),
        yaxis=dict(
            title="TPM",
        ),
        xaxis=dict(
            title="Cell type",
        ),
    )

    apc.plotly.set_yticklabel_monospaced(fig)
    apc.plotly.set_xaxis_categorical(fig)

    if output_image_filepath:
        fig.write_image(output_image_filepath)
    if output_html_filepath:
        fig.write_html(output_html_filepath, config=create_save_fig_config(width, height))

    return fig


class StatisticalMeasure(Enum):
    """Enum for statistical significance measures"""

    P_VALUE = "p-value"
    Q_VALUE = "q-value"


def _get_significance_matrix(
    expression_rows: pd.DataFrame, measure: StatisticalMeasure = StatisticalMeasure.Q_VALUE
) -> pd.DataFrame:
    """Create a matrix of p- or q-values for the pairwise comparisons of expression
        between conditions.

    Args:
        expression_rows: DataFrame containing expression data with condition and replicate columns

    Returns:
        pd.DataFrame: Significance matrix with cell types as rows and columns
    """
    significance_matrix = {
        stage1: {
            stage2: -1 if stage1 == stage2 else None for stage2 in SALPINGOECA_EXPRESSION_STAGES
        }
        for stage1 in SALPINGOECA_EXPRESSION_STAGES
    }

    for cond1 in SALPINGOECA_EXPRESSION_STAGES:
        for cond2 in SALPINGOECA_EXPRESSION_STAGES:
            if cond1 != cond2:
                pair_data = expression_rows[
                    (expression_rows["Cond 1"] == cond1) & (expression_rows["Cond 2"] == cond2)
                ]

                if not pair_data.empty:
                    significance_matrix[cond1][cond2] = pair_data.iloc[0][measure.value]

    sig_df = pd.DataFrame(significance_matrix)

    # Make the matrix symmetric by filling in missing values from the other triangle.
    sig_df = sig_df.combine_first(sig_df.T)

    return sig_df


def _categorize_pvalue(p: float) -> float:
    if p == -1:
        return SIGNIFICANCE_COLOR_VALUE_MAP["diagonal"]
    elif p < 0.001:
        return SIGNIFICANCE_COLOR_VALUE_MAP["p < 0.001"]
    elif p < 0.01:
        return SIGNIFICANCE_COLOR_VALUE_MAP["p < 0.01"]
    elif p < 0.05:
        return SIGNIFICANCE_COLOR_VALUE_MAP["p < 0.05"]
    else:
        return SIGNIFICANCE_COLOR_VALUE_MAP["n.s."]


def _create_pvalue_label(p: float) -> str:
    significance = f"<b>q-value:</b> {p}<br>"

    if p == -1:
        return None
    elif p < 0.001:
        return significance + "***"
    elif p < 0.01:
        return significance + "**"
    elif p < 0.05:
        return significance + "*"
    else:
        return significance + " (n.s.)"


def _create_heatmap_trace(
    sig_df: pd.DataFrame, significance_colors: list[str] = DEFAULT_SIGNIFICANCE_COLORS
) -> go.Heatmap:
    """Create a heatmap trace for a significance matrix.

    Args:
        sig_df: DataFrame containing p- or q-values

    Returns:
        go.Heatmap: Heatmap trace
    """
    categorical_df = sig_df.map(lambda x: _categorize_pvalue(x) if not pd.isna(x) else x)
    labels = sig_df.map(lambda x: _create_pvalue_label(x))

    colorscale = _create_significance_colorscale(significance_colors)

    return go.Heatmap(
        z=categorical_df.values,
        x=categorical_df.columns,
        y=categorical_df.index,
        colorscale=colorscale,
        hoverongaps=False,
        hoverinfo="text",
        text=labels,
        colorbar=SIGNIFICANCE_COLORBAR_PARAMS,
    )


def plot_significance_heatmap(
    target: str,
    expression_filepath: pd.DataFrame,
    salpingoeca_map_ids: dict,
    width: int = 470,
    height: int = 400,
    output_image_filepath: str = None,
    output_html_filepath: str = None,
    significance_colors: list[str] = DEFAULT_SIGNIFICANCE_COLORS,
):
    """Plot the significance of a given symbol in the form of a heatmap.

    Args:
        target: HGNC gene symbol to plot
        expression_filepath: Path to the differential expression TSV file
        salpingoeca_map_ids: Dictionary mapping HGNC gene symbols to Salpingoeca protein IDs
        width: Width of the figure
        height: Height of the figure
        output_image_filepath: Path to save the output image file
        output_html_filepath: Path to save the output HTML file

    Returns:
        fig: Plotly figure
    """
    expression_rows = _load_expression_rows(expression_filepath, salpingoeca_map_ids, target)

    significance_df = _get_significance_matrix(expression_rows)
    trace = _create_heatmap_trace(significance_df, significance_colors)

    fig = go.Figure(trace)

    fig.update_layout(
        xaxis=dict(title="Condition 2", side="top", scaleanchor="y"),
        yaxis=dict(title="Condition 1", autorange="reversed", scaleanchor="x"),
        width=width,
        height=height,
        title="P-value Significance Matrix",
        margin=dict(
            t=20,
        ),
    )

    apc.plotly.set_xaxis_categorical(fig)
    apc.plotly.set_yaxis_categorical(fig)
    apc.plotly.hide_axis_lines(fig)

    if output_image_filepath:
        fig.write_image(output_image_filepath)
    if output_html_filepath:
        fig.write_html(output_html_filepath, config=create_save_fig_config(width, height))

    return fig


def plot_expression_boxplot_and_heatmap(
    symbol: str,
    expression_filepath: pd.DataFrame,
    salpingoeca_map_ids: dict,
    width: int = 700,
    height: int = 325,
    output_image_filepath: str = None,
    output_html_filepath: str = None,
    celltype_colormap: dict[SalpingoecaStage, str] = SALPINGOECA_CELLTYPE_COLORMAP,
    significance_colors: list[str] = DEFAULT_SIGNIFICANCE_COLORS,
):
    """Plot the expression of a given symbol in the form of a boxplot and a heatmap.

    Args:
        symbol: HGNC gene symbol to plot
        expression_filepath: Path to the differential expression TSV file
        salpingoeca_map_ids: Dictionary mapping HGNC gene symbols to Salpingoeca protein IDs
        width: Width of the figure
        height: Height of the figure
        output_image_filepath: Path to save the output image file
        output_html_filepath: Path to save the output HTML file

    Returns:
        fig: Plotly figure
    """
    expression_rows = _load_expression_rows(expression_filepath, salpingoeca_map_ids, symbol)

    boxplot_df = _convert_expression_rows_to_boxplot_df(expression_rows)
    heatmap_df = _get_significance_matrix(expression_rows)

    # Create figure with secondary y-axis
    fig = make_subplots(
        rows=1, cols=2, column_widths=[0.65, 0.35], horizontal_spacing=0.28, subplot_titles=("", "")
    )

    boxplot_traces = _create_boxplot_traces(boxplot_df, colormap=celltype_colormap)
    for trace in boxplot_traces:
        fig.add_trace(trace, row=1, col=1)

    heatmap_trace = _create_heatmap_trace(heatmap_df, significance_colors)
    fig.add_trace(heatmap_trace, row=1, col=2)

    fig.update_layout(
        xaxis=dict(
            title="Cell type",
        ),
        yaxis=dict(
            title="TPM",
        ),
        xaxis2=dict(
            side="top",
            tickangle=-45,
        ),
        yaxis2=dict(
            scaleanchor="x2",
            autorange="reversed",
        ),
        width=width,
        height=height,
    )

    fig.update_layout(
        yaxis2=dict(),
        margin=dict(
            t=0,
            l=0,
            r=0,
        ),
    )

    fig.add_annotation(
        text=f"{symbol} ({salpingoeca_map_ids[symbol]}) expression in <i>S. rosetta</i>",
        xref="paper",
        yref="paper",
        x=0,
        y=1.3,
        showarrow=False,
        font=dict(
            size=apc.style_defaults.TITLE_FONT_SIZE,
            family="SuisseIntl-SemiBold",
        ),
    )

    apc.plotly.set_xaxis_categorical(fig, row=1, col=1)

    apc.plotly.set_xaxis_categorical(fig, row=1, col=2)
    apc.plotly.set_yaxis_categorical(fig, row=1, col=2)

    apc.plotly.hide_axis_lines(fig, row=1, col=2)

    if output_image_filepath:
        fig.write_image(output_image_filepath)
    if output_html_filepath:
        fig.write_html(output_html_filepath, config=create_save_fig_config(width, height))

    return fig
