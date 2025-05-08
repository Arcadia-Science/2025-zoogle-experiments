from pathlib import Path

import arcadia_pycolor as apc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scipy.stats

from zoogletools.constants import DEFAULT_FUNNEL_COLOR_LIST, DEFAULT_SPECIES_COLORS
from zoogletools.dendropy_utils import get_species_distances

apc.plotly.setup()


def create_funnel_plot(
    filter_counts: dict[str, int],
    title: str = "Protein filtering funnel chart",
    height: int = 400,
    width: int = 600,
    label_offset: int = 300,
    color_list: list[str] = DEFAULT_FUNNEL_COLOR_LIST,
    save_filepath: str | None = None,
) -> go.Figure:
    """
    Create a funnel plot showing the number of items remaining after each filtering step.

    Args:
        filter_counts (dict): Dictionary mapping filter names to counts of remaining items
        title (str, optional): Title for the plot. Defaults to "Protein filtering funnel chart".
        height (int, optional): Height of the plot in pixels. Defaults to 400.
        width (int, optional): Width of the plot in pixels. Defaults to 600.
        label_offset (int, optional): Offset for labels from the funnel bars. Defaults to 300.
        color_list (list, optional): List of hex color codes for the funnel sections.
            Defaults to DEFAULT_FUNNEL_COLOR_LIST.

    Returns:
        plotly.graph_objects.Figure: The funnel plot figure
    """

    # Create a list of steps and counts
    steps = list(filter_counts.keys())
    counts = list(filter_counts.values())

    # Calculate percentage of items remaining after each filter
    percentages = [count / counts[0] * 100 for count in counts]
    percentage_texts = [f"{p:.2f}%" for p in percentages]

    # Create text to display on hover
    hover_texts = []
    for i, (step, count) in enumerate(zip(steps, counts, strict=False)):
        if i == 0:
            hover_texts.append(f"Starting items: {count}<br>100% remaining")
        else:
            previous_count = counts[i - 1]
            reduction = previous_count - count
            percent_reduction = (reduction / previous_count) * 100
            hover_texts.append(
                f"Filter: {step}<br>"
                + f"Items remaining: {count}<br>"
                + f"Removed: {reduction} ({percent_reduction:.1f}%)<br>"
                + f"{percentage_texts[i]} remaining"
            )

    # Create the funnel chart
    fig = go.Figure(
        go.Funnel(
            y=steps,
            x=counts,
            textinfo="none",  # Turn off default text
            marker={"color": color_list},
            connector={"fillcolor": apc.denim.hex_code},
            hovertemplate="%{customdata}<extra></extra>",
            customdata=hover_texts,
        )
    )

    # Add direct annotations for each bar with fixed offset
    for i, (step, count, percentage) in enumerate(zip(steps, counts, percentages, strict=False)):
        if i == 0:
            x_offset = 0.5
            xanchor = "center"
            font_color = "white"
        else:
            x_offset = count / 2 + label_offset
            xanchor = "left"
            font_color = "black"

        fig.add_annotation(
            x=x_offset,  # Use the configurable offset
            y=step,
            text=f"{count:,}({percentage:.1f}%)",
            showarrow=False,
            xanchor=xanchor,
            font=dict(size=14, color=font_color, family=apc.style_defaults.MONOSPACE_FONT_PLOTLY),
        )

    # Update layout for better visualization
    fig.update_layout(
        title={"text": title, "y": 0.95, "x": 0.5, "xanchor": "center", "yanchor": "top"},
        height=height,
        width=width,
    )

    apc.plotly.set_yaxis_categorical(fig)
    apc.plotly.hide_yaxis_line(fig)

    if save_filepath:
        config = {
            "toImageButtonOptions": {
                "format": "svg",
                "width": width,
                "height": height,
            }
        }
        fig.write_html(save_filepath, config=config)

    return fig


def _map_identifiers(
    identifiers_filepath: str | Path, input_id: str, input_type: str, output_type: str
) -> str | list[str] | None:
    """
    Map identifiers between different identifier types

    Args:
        input_id (str): Input identifier to map
    """

    identifiers = pd.read_csv(identifiers_filepath, sep="\t")

    if input_type not in identifiers.columns:
        raise ValueError(f"Input type {input_type} not found in identifiers")
    if output_type not in identifiers.columns:
        raise ValueError(f"Output type {output_type} not found in identifiers")

    hits = set(identifiers[identifiers[input_type] == input_id][output_type].values)
    if not hits:
        return None
    elif len(hits) == 1:
        return hits.pop()
    print(f"More than one hit found for input identifier: {input_id}")
    return list(hits)


def _create_regression_line(data: pd.DataFrame) -> tuple[go.Scatter, float]:
    x = data["species_dist"]
    y = data["trait_dist"]

    slope, intercept, r_value, _, _ = scipy.stats.linregress(x, y)
    r_squared_value = r_value**2

    x_range = np.linspace(0, data["species_dist"].max() * 1, 100)
    y_range = slope * x_range + intercept

    return (
        go.Scatter(
            x=x_range,
            y=y_range,
            mode="lines",
            line=dict(color="black", width=2, dash="dash"),
            name=f"Regression Line (r² = {r_squared_value:.3f})",
            hovertext=f"(r² = {r_squared_value**2:.3f})",
        ),
        r_squared_value,
    )


def _style_scatter_figure(
    fig: go.Figure,
    data: pd.DataFrame,
    gene_symbol: str,
    uniprot_id: str,
    width: int,
    height: int,
    lock_lower_ylimit: bool,
    invert_yaxis: bool,
    invert_xaxis: bool,
) -> None:
    fig.update_layout(
        title=f"{gene_symbol} ({uniprot_id}) phylogenetic distance scatterplot",
        xaxis_title="Distance from humans",
        yaxis_title="Trait distance",
        width=width,
        height=height,
    )

    if lock_lower_ylimit:
        fig.update_layout(
            yaxis_range=[0 - data["trait_dist"].max() * 0.1, data["trait_dist"].max() * 1.1]
        )

    if invert_yaxis:
        fig.update_layout(yaxis=dict(autorange="reversed"))

    if invert_xaxis:
        fig.update_layout(xaxis=dict(autorange="reversed"))

    fig.update_layout(legend=dict(orientation="v", yanchor="top", y=1.02, xanchor="left", x=1))

    apc.plotly.set_ticklabel_monospaced(fig)


SYMBOL_MAP = {"single": "circle", "multiple": "circle-cross"}


def create_save_fig_config(width: int, height: int) -> dict:
    """Create configuration for saving figures as SVGs.

    Args:
        width: Width of the figure
        height: Height of the figure

    Returns:
        dict: Configuration for saving figures
    """
    config = {
        "toImageButtonOptions": {
            "format": "svg",
            "width": width,
            "height": height,
        }
    }
    return config


def phylogenetic_distance_scatter(
    tree_filepath: str | Path,
    identifiers_filepath: str | Path,
    input_id: str,
    input_type: str,
    data_dirpath: str | Path,
    width: int = 750,
    height: int = 500,
    lock_lower_ylimit: bool = False,
    invert_yaxis: bool = False,
    invert_xaxis: bool = False,
    annotate_multiple_copies: bool = False,
    image_filepath: str = None,
    html_filepath: str = None,
    color_dictionary: dict = DEFAULT_SPECIES_COLORS,
) -> None:
    uniprot_id = _map_identifiers(identifiers_filepath, input_id, input_type, "uniprot_id")
    gene_symbol = _map_identifiers(identifiers_filepath, input_id, input_type, "symbol")

    data_dirpath = Path(data_dirpath)
    filepath = data_dirpath / f"per-ref-protein/{uniprot_id}.tsv"
    data = pd.read_csv(filepath, sep="\t")

    data = (
        data.groupby("nonref_species")
        .agg(
            trait_dist=("trait_dist", "first"),
            nonref_protein=("nonref_protein", "first"),
            pvalue_rowwise=("pvalue_rowwise", "first"),
            pvalue_colwise=("pvalue_colwise", "first"),
            gene_count=("trait_dist", "count"),
            copy_number=("trait_dist", lambda x: "single" if len(x) == 1 else "multiple"),
        )
        .reset_index()
    )

    distances = get_species_distances(tree_filepath)

    data = data.merge(distances, on="nonref_species", how="inner")
    data.sort_values(by="species_dist", ascending=True, inplace=True)

    if annotate_multiple_copies:
        symbol = "copy_number"
    else:
        symbol = None

    hover_data = {
        "nonref_species": True,
        "gene_count": True,
        "copy_number": True,
        "species_dist": ":.2f",
        "pvalue_rowwise": ":.2f",
        "pvalue_colwise": ":.2f",
    }

    fig = px.scatter(
        data,
        x="species_dist",
        y="trait_dist",
        hover_data=hover_data,
        color="nonref_species",
        color_discrete_map=color_dictionary,
        symbol=symbol,
        symbol_map=SYMBOL_MAP,
    )
    fig.update_traces(marker=dict(size=12, line_width=1, line_color=apc.paper))

    regression_line, r_squared_value = _create_regression_line(data)
    fig.add_trace(regression_line)

    fig.add_annotation(
        x=0.03,
        y=1,
        text=f"r² = {r_squared_value:.3f}",
        showarrow=False,
        xref="paper",
        yref="paper",
        font=dict(family=apc.style_defaults.MONOSPACE_FONT_PLOTLY, size=12),
    )

    _style_scatter_figure(
        fig,
        data,
        gene_symbol,
        uniprot_id,
        width,
        height,
        lock_lower_ylimit,
        invert_yaxis,
        invert_xaxis,
    )

    if image_filepath:
        fig.write_image(image_filepath)

    config = create_save_fig_config(width, height)

    if html_filepath:
        fig.write_html(html_filepath, config=config)

    return fig


def phylogenetic_distance_scatter_zoogle(
    data_filepath: str | Path,
    tree_filepath: str | Path,
    width: int = 750,
    height: int = 500,
    lock_lower_ylimit: bool = False,
    invert_yaxis: bool = False,
    invert_xaxis: bool = False,
    annotate_multiple_copies: bool = False,
    color_dictionary: dict = DEFAULT_SPECIES_COLORS,
    image_filepath: str = None,
    html_filepath: str = None,
) -> go.Figure:
    data = pd.read_csv(data_filepath, sep="\t")
    gene_symbol = data["hgnc_gene_symbol"].iloc[0]
    uniprot_id = data["ref_protein"].iloc[0]

    data = (
        data.groupby("nonref_species")
        .agg(
            # Because values are already sorted by minimum trait distance,
            # we can take the first entry.
            trait_dist=("trait_dist", "first"),
            nonref_protein=("nonref_protein", "first"),
            gene_count=("trait_dist", "count"),
            copy_number=("trait_dist", lambda x: "single" if len(x) == 1 else "multiple"),
        )
        .reset_index()
    )

    distances = get_species_distances(tree_filepath)

    data = data.merge(distances, on="nonref_species", how="left")
    data.sort_values(by="species_dist", ascending=True, inplace=True)

    if annotate_multiple_copies:
        symbol = "copy_number"
    else:
        symbol = None

    hover_data = {
        "nonref_species": True,
        "gene_count": True,
        "copy_number": True,
        "species_dist": ":.2f",
    }

    fig = px.scatter(
        data,
        x="species_dist",
        y="trait_dist",
        hover_data=hover_data,
        color="nonref_species",
        color_discrete_map=color_dictionary,
        symbol=symbol,
        symbol_map=SYMBOL_MAP,
    )
    fig.update_traces(marker=dict(size=12, line_width=1, line_color=apc.paper))

    regression_line, r_squared_value = _create_regression_line(data)
    fig.add_trace(regression_line)

    fig.add_annotation(
        x=0.03,
        y=1,
        text=f"r² = {r_squared_value:.3f}",
        showarrow=False,
        xref="paper",
        yref="paper",
        font=dict(family=apc.style_defaults.MONOSPACE_FONT_PLOTLY, size=12),
    )

    _style_scatter_figure(
        fig,
        data,
        gene_symbol,
        uniprot_id,
        width,
        height,
        lock_lower_ylimit,
        invert_yaxis,
        invert_xaxis,
    )
    if image_filepath:
        fig.write_image(image_filepath)

    config = create_save_fig_config(width, height)

    if html_filepath:
        fig.write_html(html_filepath, config=config)

    return fig
