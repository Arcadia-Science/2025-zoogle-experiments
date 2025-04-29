import arcadia_pycolor as apc
import plotly.graph_objects as go

apc.plotly.setup()


DEFAULT_FUNNEL_COLOR_LIST = [
    apc.aster,
    apc.lapis,
    apc.aegean,
    apc.vital,
    apc.seaweed,
    apc.teal,
    apc.glass,
    apc.matcha,
    apc.lime,
]


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
