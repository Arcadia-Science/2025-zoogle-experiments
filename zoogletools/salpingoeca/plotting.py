from difflib import get_close_matches

import arcadia_pycolor as apc
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

apc.plotly.setup()

SALPINGOECA_EXPRESSION_STAGES = ["Slow Swimmer", "Fast Swimmer", "Rosette", "Thecate"]
SALPINGOECA_CELLTYPE_COLORMAP = {
    "Slow Swimmer": apc.seaweed.hex_code,
    "Fast Swimmer": apc.crow.hex_code,
    "Rosette": apc.aster.hex_code,
    "Thecate": apc.dragon.hex_code,
}

SIGNIFICANCE_COLORSCALE = [
    [0, "white"],
    [0.01, "white"],
    [0.01, apc.lapis],  # p < 0.001
    [0.25, apc.lapis],
    [0.25, apc.aegean],  # p < 0.01
    [0.5, apc.aegean],
    [0.5, apc.vital],  # p < 0.05
    [0.75, apc.vital],
    [0.75, apc.dove],  # Not significant
    [1, apc.dove],
]
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
    salpingoeca_results = pd.read_csv(
        results_filepath, sep="\t", usecols=["hgnc_gene_symbol", "nonref_protein"]
    )
    salpingoeca_expression = pd.read_csv(diffex_filepath, sep="\t", usecols=["Uniprot ID"])

    salpingoeca_map = pd.merge(
        salpingoeca_results,
        salpingoeca_expression[["Uniprot ID"]].drop_duplicates(),
        left_on="nonref_protein",
        right_on="Uniprot ID",
        how="inner",
    )[["hgnc_gene_symbol", "nonref_protein"]].drop_duplicates()

    return dict(
        zip(salpingoeca_map["hgnc_gene_symbol"], salpingoeca_map["nonref_protein"], strict=True)
    )


def _convert_expression_rows_to_boxplot_df(expression_rows):
    cell_types = []
    tpm_values = []

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


def _create_boxplot_traces(boxplot_df, showlegend=False):
    traces = []
    for stage in SALPINGOECA_EXPRESSION_STAGES:
        traces.append(
            go.Box(
                x=boxplot_df[boxplot_df["Cell type"] == stage]["Cell type"],
                y=boxplot_df[boxplot_df["Cell type"] == stage]["TPM"],
                line=dict(color=SALPINGOECA_CELLTYPE_COLORMAP[stage]),
                boxpoints="all",
                showlegend=showlegend,
            )
        )
    return traces


def _create_save_fig_config(width, height):
    config = {
        "toImageButtonOptions": {
            "format": "svg",
            "width": width,
            "height": height,
        }
    }
    return config


def _load_expression_rows(expression_filepath: str, salpingoeca_map_ids: dict, symbol: str):
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
    output_image_filepath: str = None,
    output_html_filepath: str = None,
    width=800,
    height=600,
):
    expression_rows = _load_expression_rows(expression_filepath, salpingoeca_map_ids, symbol)
    boxplot_df = _convert_expression_rows_to_boxplot_df(expression_rows)

    traces = _create_boxplot_traces(boxplot_df, showlegend=False)
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
        fig.write_html(output_html_filepath, config=_create_save_fig_config(width, height))

    return fig


def _get_significance_matrix(expression_rows):
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
                    significance_matrix[cond1][cond2] = pair_data.iloc[0]["q-value"]

    sig_df = pd.DataFrame(significance_matrix)

    # Make the matrix symmetric by filling in missing values from the other triangle.
    sig_df = sig_df.combine_first(sig_df.T)

    return sig_df


def categorize_pvalue(p):
    if p == -1:
        return 0  # Diagonal
    elif p < 0.001:
        return 0.24
    elif p < 0.01:
        return 0.49
    elif p < 0.05:
        return 0.74
    else:
        return 1


def create_pvalue_label(p):
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


def _create_heatmap_trace(sig_df):
    categorical_df = sig_df.map(lambda x: categorize_pvalue(x) if not pd.isna(x) else x)
    labels = sig_df.map(lambda x: create_pvalue_label(x))

    return go.Heatmap(
        z=categorical_df.values,
        x=categorical_df.columns,
        y=categorical_df.index,
        colorscale=SIGNIFICANCE_COLORSCALE,
        hoverongaps=False,
        hoverinfo="text",
        text=labels,
        colorbar=SIGNIFICANCE_COLORBAR_PARAMS,
    )


def plot_significance_heatmap(
    target, expression_filepath, salpingoeca_map_ids, width=470, height=400
):
    expression_rows = _load_expression_rows(expression_filepath, salpingoeca_map_ids, target)

    significance_df = _get_significance_matrix(expression_rows)
    trace = _create_heatmap_trace(significance_df)

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

    return fig


def plot_expression_boxplot_and_heatmap(
    symbol: str,
    expression_filepath: pd.DataFrame,
    salpingoeca_map_ids: dict,
    width=700,
    height=325,
    save_image_filepath: str = None,
    save_html_filepath: str = None,
):
    expression_rows = _load_expression_rows(expression_filepath, salpingoeca_map_ids, symbol)

    boxplot_df = _convert_expression_rows_to_boxplot_df(expression_rows)
    heatmap_df = _get_significance_matrix(expression_rows)

    # Create figure with secondary y-axis
    fig = make_subplots(
        rows=1, cols=2, column_widths=[0.65, 0.35], horizontal_spacing=0.28, subplot_titles=("", "")
    )

    boxplot_traces = _create_boxplot_traces(boxplot_df)
    for trace in boxplot_traces:
        fig.add_trace(trace, row=1, col=1)

    heatmap_trace = _create_heatmap_trace(heatmap_df)
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

    if save_image_filepath:
        fig.write_image(save_image_filepath)
    if save_html_filepath:
        fig.write_html(save_html_filepath, config=_create_save_fig_config(width, height))

    return fig
