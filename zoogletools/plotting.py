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
) -> str:
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

    hits = identifiers[identifiers[input_type] == input_id][output_type].values
    if len(hits) == 0:
        return None
    elif len(hits) == 1:
        return hits[0]
    elif len(hits) > 1:
        if len(set(hits)) == 1:
            return hits[0]
        else:
            print(f"More than one hit found for input identifier: {input_id}")
            return hits


def phylogenetic_distance_scatterplot(
    tree_filepath: str | Path,
    complexities_filepath: str | Path,
    identifiers_filepath: str | Path,
    data_dirpath: str | Path,
    input_id: str,
    input_type: str,
    xaxis_data: str = "species_dist",
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

    if isinstance(data_dirpath, str):
        data_dirpath = Path(data_dirpath)

    filepath = data_dirpath / f"per-ref-protein/{uniprot_id}.tsv"
    data = pd.read_csv(filepath, sep="\t")

    data = (
        data.groupby("nonref_species")
        .agg(
            # Because values are already sorted by minimum trait distance,
            # we can take the first entry.
            trait_dist=("trait_dist", "first"),
            nonref_protein=("nonref_protein", "first"),
            pvalue_rowwise=("pvalue_rowwise", "first"),
            pvalue_colwise=("pvalue_colwise", "first"),
            gene_count=("trait_dist", "count"),
            copy_number=("trait_dist", lambda x: "single" if len(x) == 1 else "multiple"),
        )
        .reset_index()
    )

    distances = get_species_distances(tree_filepath, complexities_filepath)

    data = data.merge(distances, on="nonref_species", how="left")
    data.sort_values(by="species_dist", ascending=True, inplace=True)

    symbol_map = {"single": "circle", "multiple": "circle-cross"}

    if annotate_multiple_copies:
        symbol = "copy_number"
    else:
        symbol = None

    fig = px.scatter(
        data,
        x=xaxis_data,
        y="trait_dist",
        hover_data=["nonref_species"],
        color="nonref_species",
        color_discrete_map=color_dictionary,
        symbol=symbol,
        symbol_map=symbol_map,
    )

    fig.update_traces(marker=dict(size=12, line_width=1, line_color=apc.paper))

    if xaxis_data == "species_dist":
        xaxis_title = "Phylogenetic distance"
    elif xaxis_data == "n_cells":
        xaxis_title = "Organismal complexity"
        # make xaxis log
        fig.update_layout(xaxis_type="log")
    elif xaxis_data == "n_cell_types":
        xaxis_title = "Organismal complexity"
        # make xaxis log
        fig.update_layout(xaxis_type="log")

        # Add linear regression line
    x = data["species_dist"]
    y = data["trait_dist"]

    # Calculate the linear regression
    slope, intercept, r_value, _, _ = scipy.stats.linregress(x, y)

    r_squared_value = r_value**2

    # Create x values for the line
    x_range = np.linspace(0, data["species_dist"].max() * 1, 100)
    y_range = slope * x_range + intercept

    # Add regression line
    fig.add_trace(
        go.Scatter(
            x=x_range,
            y=y_range,
            mode="lines",
            line=dict(color="black", width=2, dash="dash"),
            name=f"Regression Line (r² = {r_squared_value:.3f})",
            hovertext=f"(r² = {r_squared_value**2:.3f})",
        )
    )

    # Add r-squared value to top left of plot canvas
    fig.add_annotation(
        x=0.03,
        y=1,
        text=f"r² = {r_squared_value:.3f}",
        showarrow=False,
        xref="paper",
        yref="paper",
        font=dict(family=apc.style_defaults.MONOSPACE_FONT_PLOTLY, size=12),
    )

    fig.update_layout(
        title=f"{input_id} phylogenetic distance scatterplot",
        xaxis_title=xaxis_title,
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

    if image_filepath:
        fig.write_image(image_filepath)

    config = {
        "toImageButtonOptions": {
            "format": "svg",
            "width": width,
            "height": height,
        }
    }

    if html_filepath:
        fig.write_html(html_filepath, config=config)

    return fig
