from pathlib import Path

import arcadia_pycolor as apc
import pandas as pd
import plotly.graph_objects as go
import scanpy as sc

from zoogletools.ciona.constants import (
    CAO_DATA_DIRPATH,
    PIEKARZ_DATA_DIRPATH,
    STAGE_COLORS,
    CionaStage,
)
from zoogletools.ciona.data_processing import _append_tech_replicate_to_barcode, load_cell_clusters
from zoogletools.plotting import create_save_fig_config, hex_to_plotly_rgba


def load_cao_cell_barcodes(
    data_dir: str | Path = CAO_DATA_DIRPATH,
) -> pd.DataFrame:
    cell_clusters = load_cell_clusters(data_dir)
    cao_cell_barcodes = cell_clusters[["stage_replicate", "barcode"]].copy()
    cao_cell_barcodes["stage_replicate"] = "Cao_" + cao_cell_barcodes["stage_replicate"].astype(str)

    return cao_cell_barcodes


def load_piekarz_scrnaseq_data(
    stage: CionaStage,
    data_dir: str | Path = PIEKARZ_DATA_DIRPATH,
) -> sc.AnnData:
    """
    Load the scRNA-seq data for a given developmental stage directly from the Piekarz dataset.

    There are a number of inconsistencies in the Piekarz data, including:
      - inconsistent column names
      - incorrect replicate numbers
      - incorrect stage names
      - different cell barcode formats compared to Cao et al. for technical replicates

    This function resolves the inconsistent column names and modifies the barcodes
    to be able to merge the Cao and Piekarz datasets.

    Unlike the similar function in the `data_processing.py` module, this function does
    *not* resolve incorrect replicate numbers and stage names. It is simply used to load
    the scRNA-seq data without modifications for the purposes of the disambiguation demonstration.

    Parameters:
    stage (CionaStage): The developmental stage to load the data for.
    data_dir (str | Path): The directory containing the scRNA-seq data.

    Returns:
    adata (AnnData): The scRNA-seq data for the given developmental stage.
    """
    index = CionaStage.ordered_stages().index(stage)

    data_dir = Path(data_dir)
    filepath = data_dir / f"{index + 1}_{stage}/final/{stage}.h5ad"

    # The midTII stage file is named differently.
    if stage == "midTII":
        filepath = data_dir / f"{index + 1}_{stage}/final/{stage}_KY21.h5ad"

    adata = sc.read_h5ad(filepath)

    # For some reason, the replicate column is named "larva" in the iniG stage.
    if stage == "iniG":
        adata.obs.rename(columns={"larva": "rep"}, inplace=True)

    # For some reason, the replicate column is named "replicate"
    # (rather than "rep") in the midTII stage file.
    if stage == "midTII":
        adata.obs.rename(columns={"replicate": "rep"}, inplace=True)
        adata.obs["rep"] = adata.obs["rep"].cat.rename_categories(
            {"midTII-1": "rep1", "midTII-2": "rep2"}
        )

    # Unlike the other stages, the latTI and latTII files also include technical replicates,
    # which are labeled in an additional "tech" column.
    # To be able to merge these with the Cao cell annotations,
    # we need to append "-1" and "-2" to the barcode.
    if stage in ["latTI", "latTII"]:
        adata.obs["barcode"] = adata.obs.apply(_append_tech_replicate_to_barcode, axis=1)

    adata.obs["stage"] = stage
    adata.obs["stage_replicate"] = adata.obs["stage"] + "_" + adata.obs["rep"].astype(str)

    return adata


def load_piekarz_cell_barcodes(
    data_dir: str | Path = PIEKARZ_DATA_DIRPATH,
) -> pd.DataFrame:
    data_dir = Path(data_dir)

    piekarz_cell_barcodes = pd.DataFrame()

    for stage in CionaStage.ordered_stages():
        result = load_piekarz_scrnaseq_data(stage, data_dir)
        barcodes = result.obs[["stage_replicate", "barcode"]].copy().reset_index(drop=True)
        barcodes["stage_replicate"] = "Piekarz_" + barcodes["stage_replicate"].astype(str)
        piekarz_cell_barcodes = pd.concat([piekarz_cell_barcodes, barcodes])

    return piekarz_cell_barcodes


def _sort_by_stage_and_replicate(stage_replicate_labels: list[str]) -> list[str]:
    """
    Sort a list of stage replicate labels by stage and replicate number,
    where the stage replicate labels are of the form "<stage>_<replicate>".

    For example, "iniG_rep1" should come before "iniG_rep2".

    Parameters:
        stage_replicate_labels (list[str]): List of stage replicate labels to sort.

    Returns:
        list[str]: Sorted list of stage replicate labels.
    """
    return sorted(
        stage_replicate_labels,
        key=lambda x: (
            CionaStage.ordered_stages().index(x.split("_")[1]),  # Sort by stage.
            int(x.split("_")[2].replace("rep", "")),  # Sort by replicate number.
        ),
    )


def plot_stage_replicate_sankey(
    merged_cell_barcodes: pd.DataFrame,
    width: int = 500,
    height: int = 800,
    image_filepath: str | Path | None = None,
    html_filepath: str | Path | None = None,
):
    """Create a Sankey diagram showing mapping between Cao and Piekarz stage replicates.

    Args:
        merged_cell_barcodes (pd.DataFrame): DataFrame with stage_replicate_cao
            and stage_replicate_piekarz columns
        width (int): Width of the plot in pixels
        height (int): Height of the plot in pixels
        image_filepath (str | Path | None): Path to save the image file
        html_filepath (str | Path | None): Path to save the HTML file

    Returns:
        plotly.graph_objects.Figure: The Sankey diagram
    """

    source_values = _sort_by_stage_and_replicate(
        merged_cell_barcodes["stage_replicate_cao"].unique()
    )
    target_values = _sort_by_stage_and_replicate(
        merged_cell_barcodes["stage_replicate_piekarz"].unique()
    )

    # Convert to categorical for sorting
    for col, categories in [
        ("stage_replicate_cao", source_values),
        ("stage_replicate_piekarz", target_values),
    ]:
        merged_cell_barcodes[col] = pd.Categorical(
            merged_cell_barcodes[col], categories=categories, ordered=True
        )

    merged_cell_barcodes = merged_cell_barcodes.sort_values(
        ["stage_replicate_cao", "stage_replicate_piekarz"]
    )

    # Calculate link values and total number of cells
    value_counts = merged_cell_barcodes.groupby(
        ["stage_replicate_cao", "stage_replicate_piekarz"]
    ).size()
    total_cells = value_counts.sum()

    # Create node information
    node_labels = source_values + target_values
    node_colors = [STAGE_COLORS[CionaStage(label.split("_")[1])] for label in node_labels]

    # Calculate y positions
    def get_cumulative_positions(totals):
        positions = [0.001]
        running_sum = 0
        for val in totals[1:]:
            running_sum += val
            positions.append(running_sum / total_cells)
        return positions

    source_y = get_cumulative_positions(value_counts.groupby(level=0).sum())
    target_y = get_cumulative_positions(value_counts.groupby(level=1).sum())

    # Create link colors
    link_colors = [hex_to_plotly_rgba(apc.stone, 0.2)] * len(value_counts)
    value_counts_idx = list(value_counts.index)

    # Highlight the largest links from each source (to any number of targets)
    # and the largest links to each target (from any number of sources).
    for groups in [value_counts.groupby(level=0), value_counts.groupby(level=1)]:
        for _, group in groups:
            max_idx = value_counts_idx.index(group.idxmax())
            link_colors[max_idx] = hex_to_plotly_rgba(apc.chateau, 0.7)

    # Create Sankey diagram
    fig = go.Figure(
        data=[
            go.Sankey(
                arrangement="snap",
                node=dict(
                    pad=15,
                    thickness=20,
                    line=dict(color="white", width=0),
                    label=node_labels,
                    color=node_colors,
                    x=[0.001] * len(source_values) + [0.999] * len(target_values),
                    y=source_y + target_y,
                ),
                link=dict(
                    source=[
                        list(source_values).index(s) for s in value_counts.index.get_level_values(0)
                    ],
                    target=[
                        len(source_values) + list(target_values).index(t)
                        for t in value_counts.index.get_level_values(1)
                    ],
                    value=value_counts.values,
                    color=link_colors,
                ),
            )
        ]
    )

    fig.update_layout(
        title_text="Stage barcode mapping: Cao vs Piekarz",
        font_size=10,
        width=width,
        height=height,
        margin=dict(
            t=30,
            l=30,
            r=30,
        ),
    )

    if image_filepath:
        fig.write_image(image_filepath)

    if html_filepath:
        fig.write_html(html_filepath, config=create_save_fig_config(width, height))

    return fig
