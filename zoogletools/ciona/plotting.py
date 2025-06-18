import os
import warnings
from enum import StrEnum
from pathlib import Path
from typing import Literal

import arcadia_pycolor as apc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
from plotly.subplots import make_subplots
from tqdm import tqdm

from zoogletools.ciona.constants import (
    BULK_RNA_SEQ_DATA_DIRPATH,
    CAO_DATA_DIRPATH,
    PIEKARZ_DATA_DIRPATH,
    TISSUE_TYPE_GRADIENTS,
    TISSUE_TYPE_PALETTE,
    CionaStage,
)
from zoogletools.ciona.data_processing import load_cell_clusters, load_ciona_scrnaseq_data
from zoogletools.ciona.identifier_mapping import CionaIDTypes, IdentifierMapper
from zoogletools.constants import (
    GRAY_GRADIENT,
    PLOTLY_TITLE_FONT,
)
from zoogletools.plotting import (
    create_save_fig_config,
)


def extract_umap_coordinates(adata: sc.AnnData) -> pd.DataFrame:
    """
    Extract the UMAP coordinates from the scRNAseq data.
    """
    umap_embeddings = adata.obsm["X_umap"]
    umap_df = pd.DataFrame(umap_embeddings, columns=["UMAP_1", "UMAP_2"], index=adata.obs_names)
    umap_df = pd.concat([umap_df, adata.obs], axis=1)

    return umap_df


def _merge_cluster_annotations(
    umap_df: pd.DataFrame,
    stage: CionaStage,
    cluster_annotations_filepath: str | Path,
) -> pd.DataFrame:
    """Merge the cluster annotations with the UMAP dataframe."""

    index = CionaStage.ordered_stages().index(stage) + 1
    index_stage = f"{index}_{stage}"

    cluster_annotations = pd.read_csv(
        cluster_annotations_filepath,
        sep="\t",
        dtype={"seurat_clusters": str},
    )

    cluster_annotations = cluster_annotations[cluster_annotations["index_stage"] == index_stage]
    cluster_annotations["top_cluster_tissue_type"] = pd.Categorical(
        cluster_annotations["top_cluster_tissue_type"],
        categories=TISSUE_TYPE_GRADIENTS.keys(),
        ordered=True,
    )

    return umap_df.merge(cluster_annotations, on="seurat_clusters", how="left")


def _merge_tissue_types(
    umap_df: pd.DataFrame, stage: CionaStage, data_dirpath: str | Path = CAO_DATA_DIRPATH
) -> pd.DataFrame:
    cell_clusters = load_cell_clusters(data_dirpath)

    cell_clusters = cell_clusters[cell_clusters["stage"] == stage]
    tissue_types = cell_clusters[["merging_barcode", "Tissue Type"]]

    umap_df = umap_df.merge(tissue_types, on="merging_barcode", how="left")
    umap_df.replace(
        {
            "Tissue Type": {
                np.nan: "unannotated",
                "muscle & heart": "muscle-heart",
                "nervous system": "nervous-system",
            }
        },
        inplace=True,
    )

    return umap_df


def _add_gene_expression(umap_df: pd.DataFrame, adata: sc.AnnData, ky_id: str) -> pd.DataFrame:
    gene_index = np.where(adata.var.index == ky_id)[0][0]
    gene_expression = adata.X[:, gene_index].toarray().flatten()
    umap_df["gene_expression"] = gene_expression

    return umap_df


def append_cell_count_barchart(
    fig: go.Figure,
    umap_df: pd.DataFrame,
    cluster_color_map: dict,
    row: int,
    col: int,
) -> go.Figure:
    cluster_stats = umap_df.groupby("seurat_clusters", observed=True).agg(
        total_cells=("gene_expression", "count"),
        expressing_cells=("gene_expression", lambda x: (x > 0).sum()),
    )

    clusters = umap_df["seurat_clusters"].unique()
    for cluster in clusters:
        if cluster_stats.loc[cluster, "expressing_cells"] > 0:
            fig.add_trace(
                go.Bar(
                    x=[f"C{cluster}"],
                    y=[cluster_stats.loc[cluster, "total_cells"]],
                    name=f"Total cells (cluster {cluster})",
                    marker_color=cluster_color_map[cluster],
                    showlegend=False,
                ),
                row=row,
                col=col,
            )


def append_expression_proportion_barchart(
    fig: go.Figure,
    umap_df: pd.DataFrame,
    cluster_color_map: dict,
    row: int,
    col: int,
) -> go.Figure:
    cluster_stats = umap_df.groupby("seurat_clusters", observed=True).agg(
        total_cells=("gene_expression", "count"),
        expressing_cells=("gene_expression", lambda x: (x > 0).sum()),
    )
    cluster_stats["non_expressing_cells"] = (
        cluster_stats["total_cells"] - cluster_stats["expressing_cells"]
    )
    cluster_stats["expressing_cells_percentage"] = (
        cluster_stats["expressing_cells"] / cluster_stats["total_cells"]
    ) * 100
    cluster_stats["non_expressing_cells_percentage"] = (
        cluster_stats["non_expressing_cells"] / cluster_stats["total_cells"]
    ) * 100

    # Add stacked bar chart to the first subplot, only for clusters with expressing cells
    clusters = umap_df["seurat_clusters"].unique()
    for cluster in clusters:
        if cluster_stats.loc[cluster, "expressing_cells"] > 0:
            fig.add_trace(
                go.Bar(
                    x=[f"C{cluster}"],
                    y=[cluster_stats.loc[cluster, "expressing_cells_percentage"]],
                    name=f"Expressing cells (cluster {cluster})",
                    marker_color=cluster_color_map[cluster],
                    showlegend=False,
                ),
                row=row,
                col=col,
            )
            fig.add_trace(
                go.Bar(
                    x=[f"C{cluster}"],
                    y=[cluster_stats.loc[cluster, "non_expressing_cells_percentage"]],
                    name="Non-expressing cells",
                    marker_color=apc.ice.hex_code,
                    showlegend=False,
                ),
                row=row,
                col=col,
            )

    fig.update_layout(
        barmode="stack",
    )


def append_expression_violin(
    fig: go.Figure,
    umap_df: pd.DataFrame,
    cluster_color_map: dict,
    row: int,
    col: int,
) -> go.Figure:
    umap_df_nonzero = umap_df[umap_df["gene_expression"] > 0]

    clusters = umap_df["seurat_clusters"].unique()
    for cluster in clusters:
        cluster_data = umap_df_nonzero[umap_df_nonzero["seurat_clusters"] == cluster].copy()

        cluster_data["cluster_label"] = cluster_data.apply(
            lambda row: f"C{row['seurat_clusters']}: {row['formatted_cluster_annotation']}", axis=1
        )

        if len(cluster_data) > 0:
            fig.add_trace(
                go.Violin(
                    x=cluster_data["cluster_label"],
                    y=cluster_data["gene_expression"],
                    name=cluster,
                    box_visible=True,
                    points="all",
                    pointpos=0,
                    jitter=0.3,
                    showlegend=False,
                    marker_color=cluster_color_map[cluster],
                ),
                row=row,
                col=col,
            )


def _add_yaxis_title(
    fig: go.Figure,
    y_domain: str,
    yaxis_title: str,
    yaxis_title_offset: float,
) -> go.Figure:
    fig.add_annotation(
        xref="paper",
        yref=y_domain,
        x=yaxis_title_offset,
        y=0.5,
        text=yaxis_title,
        font=dict(
            size=14,
            family=apc.style_defaults.DEFAULT_FONT_PLOTLY + "-SemiBold",
        ),
        showarrow=False,
        textangle=-90,
    )


def plot_expression_violin(
    stage: CionaStage,
    input_id: str,
    mapper: IdentifierMapper,
    input_id_type: str = CionaIDTypes.HGNC_GENE_SYMBOL,
    data_dirpath: str | Path = PIEKARZ_DATA_DIRPATH,
    cluster_annotations_filepath: str | Path | None = None,
    color_mode: Literal["tissue", "cluster"] = "tissue",
    width: int = 950,
    height: int = 700,
    yaxis_title_offset: float = -0.12,
    image_filepath: str | Path | None = None,
    html_filepath: str | Path | None = None,
) -> go.Figure:
    # Load cluster annotations
    if cluster_annotations_filepath is None:
        cluster_annotations_filepath = (
            data_dirpath / "cluster_annotations" / "Ciona_scRNAseq_cluster_annotations.tsv"
        )

    adata = load_ciona_scrnaseq_data(stage, data_dir=data_dirpath)
    umap_df = extract_umap_coordinates(adata)

    # Get all IDs using the provided mapper
    all_ciona_ids = mapper.map_to_all(input_id, input_id_type)

    uniprot_id = all_ciona_ids[CionaIDTypes.NONREF_PROTEIN]
    ky_id = all_ciona_ids[CionaIDTypes.KY_ID]
    hgnc_gene_symbol = all_ciona_ids[CionaIDTypes.HGNC_GENE_SYMBOL]

    umap_df = _add_gene_expression(umap_df, adata, ky_id)

    umap_df = _merge_cluster_annotations(
        umap_df=umap_df, stage=stage, cluster_annotations_filepath=cluster_annotations_filepath
    )

    if color_mode == "tissue":
        umap_df.sort_values(["top_cluster_tissue_type", "top_cluster_suffix"], inplace=True)
        unique_pairs = umap_df[["seurat_clusters", "cluster_tissue_color"]].drop_duplicates()
        cluster_color_map = dict(
            zip(unique_pairs["seurat_clusters"], unique_pairs["cluster_tissue_color"], strict=True)
        )
    elif color_mode == "cluster":
        umap_df.sort_values(["seurat_clusters"], inplace=True)
        clusters = umap_df["seurat_clusters"].unique()

        cluster_colors = apc.palettes.all_ordered.colors * (
            1 + len(clusters) // len(apc.palettes.all_ordered.colors)
        )
        cluster_colors = cluster_colors[: len(clusters)]

        cluster_color_map = dict(zip(clusters, cluster_colors, strict=True))
    else:
        raise ValueError(f"Invalid color_mode: {color_mode}")

    fig = make_subplots(
        rows=3,
        cols=1,
        row_heights=[0.15, 0.15, 0.7],
        vertical_spacing=0.05,
    )

    append_cell_count_barchart(
        fig=fig,
        umap_df=umap_df,
        cluster_color_map=cluster_color_map,
        row=1,
        col=1,
    )

    append_expression_proportion_barchart(
        fig=fig,
        umap_df=umap_df,
        cluster_color_map=cluster_color_map,
        row=2,
        col=1,
    )

    append_expression_violin(
        fig=fig,
        umap_df=umap_df,
        cluster_color_map=cluster_color_map,
        row=3,
        col=1,
    )

    fig.update_layout(
        width=width,
        height=height,
        margin=dict(t=40),
    )

    fig.update_xaxes(
        tickangle=45,
    )

    fig.add_annotation(
        xref="paper",
        yref="paper",
        x=0.5,
        y=1.05,
        text=f"Expression of <i>{hgnc_gene_symbol}</i> ({uniprot_id}) at {stage}",
        font=dict(
            size=16,
            family=apc.style_defaults.DEFAULT_FONT_PLOTLY + "-SemiBold",
        ),
        showarrow=False,
    )

    yaxis_titles = {
        "y": "# cells",
        "y2": "% expr.",
        "y3": f"{hgnc_gene_symbol} ({ky_id})",
    }

    for y in yaxis_titles.keys():
        _add_yaxis_title(
            fig=fig,
            y_domain=f"{y} domain",
            yaxis_title=yaxis_titles[y],
            yaxis_title_offset=yaxis_title_offset,
        )

    apc.plotly.set_yticklabel_monospaced(fig, row=1, col=1)
    apc.plotly.hide_xaxis_ticks(fig, row=1, col=1)
    apc.plotly.hide_xaxis_ticks(fig, row=2, col=1)
    apc.plotly.set_yticklabel_monospaced(fig, row=1, col=1)
    apc.plotly.set_yticklabel_monospaced(fig, row=2, col=1)
    apc.plotly.set_yticklabel_monospaced(fig, row=3, col=1)

    if image_filepath is not None:
        os.makedirs(os.path.dirname(image_filepath), exist_ok=True)
        fig.write_image(image_filepath)

    if html_filepath is not None:
        os.makedirs(os.path.dirname(html_filepath), exist_ok=True)
        fig.write_html(
            html_filepath,
            config=create_save_fig_config(
                width=width,
                height=height,
            ),
        )

    return fig


def plot_expression_violin_for_all_stages(
    input_id: str,
    mapper: IdentifierMapper,
    input_id_type: str = CionaIDTypes.HGNC_GENE_SYMBOL,
    data_dirpath: str | Path = PIEKARZ_DATA_DIRPATH,
    color_mode: Literal["tissue", "cluster"] = "tissue",
    output_dirpath: str | Path = None,
) -> go.Figure:
    # Get all IDs using the provided mapper
    all_ciona_ids = mapper.map_to_all(input_id, input_id_type)
    uniprot_id = all_ciona_ids[CionaIDTypes.NONREF_PROTEIN]

    output_dirpath = Path(output_dirpath)

    for stage in tqdm(CionaStage.ordered_stages()):
        index = CionaStage.ordered_stages().index(stage) + 1
        index_stage = f"{index}_{stage}"

        plot_expression_violin(
            stage=stage,
            input_id=input_id,
            mapper=mapper,
            input_id_type=input_id_type,
            data_dirpath=data_dirpath,
            color_mode=color_mode,
            image_filepath=output_dirpath
            / f"{input_id}_{uniprot_id}_expression/"
            / "violin"
            / f"{index_stage}_{input_id}_{color_mode}_violin.svg",
            html_filepath=output_dirpath
            / f"{input_id}_{uniprot_id}_expression/"
            / "violin"
            / f"{index_stage}_{input_id}_{color_mode}_violin.html",
        )


SCATTER_HOVERTEMPLATE = (
    "<b>Barcode</b>: %{customdata[0]}<br>"
    + "<b>Replicate</b>: %{customdata[1]}<br>"
    + "<b>Seurat cluster</b>: %{customdata[2]}<br>"
    + "<b>Annotation</b>: %{customdata[3]}<br>"
    + "<b>Expression</b>: %{customdata[4]}<br>"
    + "<b>Tissue type</b>: %{customdata[5]}"
)

SCATTER_CUSTOMDATA = [
    "barcode",
    "rep",
    "seurat_clusters",
    "formatted_cluster_annotation",
    "gene_expression",
    "Tissue Type",
]


def _map_position_to_row_col(position: str) -> tuple[int, int]:
    if position == "top-left":
        return 1, 1
    elif position == "top-right":
        return 1, 2
    elif position == "bottom-left":
        return 2, 1
    elif position == "bottom-right":
        return 2, 2
    else:
        raise ValueError(f"Invalid position: {position}")


def _map_position_to_domain(position: str) -> str:
    if position == "top-left":
        return ""
    elif position == "top-right":
        return "2"
    elif position == "bottom-left":
        return "3"
    elif position == "bottom-right":
        return "4"
    else:
        raise ValueError(f"Invalid position: {position}")


def _append_scatter_trace(
    fig: go.Figure,
    umap_df: pd.DataFrame,
    color_values: list,
    position: str | None = None,
) -> go.Scattergl:
    trace = go.Scattergl(
        x=umap_df["UMAP_1"],
        y=umap_df["UMAP_2"],
        mode="markers",
        marker=dict(
            color=color_values,
        ),
        customdata=umap_df[SCATTER_CUSTOMDATA],
        hovertemplate=SCATTER_HOVERTEMPLATE,
        showlegend=False,
    )

    if position is None:
        fig.add_trace(trace)
        return

    row, col = _map_position_to_row_col(position)
    fig.add_trace(trace, row=row, col=col)


def _annotate_scatter_trace(
    fig: go.Figure,
    text: str,
    position: str | None = None,
) -> go.Figure:
    if position is None:
        position = "top-left"

    xref = f"x{_map_position_to_domain(position)} domain"
    yref = f"y{_map_position_to_domain(position)} domain"

    fig.add_annotation(
        text=text,
        x=0.01,
        y=1.01,
        xref=xref,
        yref=yref,
        xanchor="left",
        yanchor="top",
        showarrow=False,
        font=dict(
            size=14,
            family=PLOTLY_TITLE_FONT,
        ),
    )


def append_seurat_clusters_trace(
    fig: go.Figure,
    umap_df: pd.DataFrame,
    position: str | None = None,
    title: str = "Seurat cluster",
) -> go.Figure:
    clusters = [str(i) for i in sorted(umap_df["seurat_clusters"].astype(int).unique())]

    cluster_colors = apc.palettes.all_ordered.colors * (
        1 + len(clusters) // len(apc.palettes.all_ordered.colors)
    )
    cluster_colors = cluster_colors[: len(clusters)]
    cluster_color_map = dict(zip(clusters, cluster_colors, strict=True))

    _append_scatter_trace(
        fig,
        umap_df,
        umap_df["seurat_clusters"].map(cluster_color_map),
        position,
    )
    _annotate_scatter_trace(
        fig,
        title,
        position,
    )


DEFAULT_COLORBAR_PARAMS = dict(
    title=dict(
        text="Expression",
        font=dict(family=apc.style_defaults.DEFAULT_FONT_PLOTLY + "-SemiBold", size=14),
    ),
    outlinewidth=0,
    thickness=15,
    ypad=0,
    len=0.45,
    y=0.78,
    x=1,
    xref="paper",
    yref="paper",
    yanchor="middle",
    tickfont=dict(
        size=12,
        family=apc.style_defaults.MONOSPACE_FONT_PLOTLY,
    ),
)


def append_expression_trace(
    fig: go.Figure,
    umap_df: pd.DataFrame,
    position: str | None = None,
    title: str = "Expression",
    expression_colormap: apc.Gradient = GRAY_GRADIENT,
    null_expression_color: str = apc.gray,
    expression_colorbar_params: dict = DEFAULT_COLORBAR_PARAMS,
) -> go.Figure:
    umap_zero_df = umap_df[umap_df["gene_expression"] == 0]
    umap_nonzero_df = umap_df[umap_df["gene_expression"] != 0]

    _append_scatter_trace(
        fig,
        umap_zero_df,
        null_expression_color,
        position,
    )

    trace = go.Scattergl(
        x=umap_nonzero_df["UMAP_1"],
        y=umap_nonzero_df["UMAP_2"],
        mode="markers",
        marker=dict(
            color=umap_nonzero_df["gene_expression"],
            colorscale=expression_colormap.to_plotly_colorscale(),
            colorbar=expression_colorbar_params,
        ),
        customdata=umap_nonzero_df[SCATTER_CUSTOMDATA],
        hovertemplate=SCATTER_HOVERTEMPLATE,
        showlegend=False,
    )

    if position is None:
        fig.add_trace(trace)
        return

    row, col = _map_position_to_row_col(position)
    fig.add_trace(trace, row=row, col=col)

    _annotate_scatter_trace(
        fig,
        title,
        position,
    )


def append_tissue_type_trace(
    fig: go.Figure,
    umap_df: pd.DataFrame,
    position: str | None = None,
    title: str = "Cluster tissue type",
) -> go.Figure:
    _append_scatter_trace(
        fig,
        umap_df,
        umap_df["cluster_tissue_color"],
        position,
    )
    _annotate_scatter_trace(
        fig,
        title,
        position,
    )


def append_cell_tissue_type_trace(
    fig: go.Figure,
    umap_df: pd.DataFrame,
    position: str | None = None,
    title: str = "Cell tissue type",
) -> go.Figure:
    _append_scatter_trace(
        fig,
        umap_df,
        umap_df["Tissue Type"].map(TISSUE_TYPE_PALETTE),
        position,
    )
    _annotate_scatter_trace(
        fig,
        title,
        position,
    )


def _append_tissue_type_legend(
    fig: go.Figure,
    legend_x: float = 1.015,
    legend_y: float = 0.47,
    y_offset: float = 0.024,
    title_offset: float = 0.03,
) -> go.Figure:
    fig.add_annotation(
        text="Tissue color",
        x=legend_x,
        y=legend_y,
        yanchor="middle",
        xanchor="left",
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(color="black", family=PLOTLY_TITLE_FONT, size=14),
    )

    for i, (tissue, color) in enumerate(TISSUE_TYPE_PALETTE.items()):
        if tissue == "muscle-heart":
            formatted_tissue_text = "Muscle / heart"
        else:
            formatted_tissue_text = tissue.capitalize().replace("-", " ")

        fig.add_annotation(
            text="■",
            x=legend_x,
            y=legend_y - title_offset - (i * y_offset),
            yanchor="middle",
            xanchor="left",
            xref="paper",
            yref="paper",
            showarrow=False,
            font=dict(color=color, size=20),
        )

        fig.add_annotation(
            text=formatted_tissue_text,
            x=legend_x + 0.03,
            y=legend_y - title_offset - (i * y_offset) - 0.003,
            yanchor="middle",
            xanchor="left",
            xref="paper",
            yref="paper",
            showarrow=False,
            font=dict(color=color, size=12),
        )


def plot_expression_umap(
    stage: CionaStage,
    input_id: str,
    input_id_type: CionaIDTypes,
    mapper: IdentifierMapper,
    data_dirpath: str | Path = PIEKARZ_DATA_DIRPATH,
    annotation_data_dirpath: str | Path = CAO_DATA_DIRPATH,
    plot_order: tuple[str, ...] = (
        "seurat_clusters",
        "expression",
        "tissue_type",
        "cell_tissue_type",
    ),
    expression_colormap: apc.Gradient = GRAY_GRADIENT,
    null_expression_color: str = apc.gray,
    expression_colorbar_params: dict = DEFAULT_COLORBAR_PARAMS,
    marker_size: int = 3,
    width: int = 850,
    height: int = 800,
    image_filepath: str | Path | None = None,
    html_filepath: str | Path | None = None,
) -> go.Figure:
    """Plot UMAP visualizations of gene expression and cell type annotations.

    Args:
        stage: The developmental stage to plot
        input_id: The gene identifier to plot
        input_id_type: The type of the input identifier
        data_dirpath: Path to the Piekarz scRNA-seq data directory
        annotation_data_dirpath: Path to the Cao annotation data directory
        mapper: IdentifierMapper instance for ID conversion

    Returns:
        A plotly Figure object with 4 UMAP plots showing:
        - Seurat clusters
        - Gene expression
        - Cluster tissue types
        - Cell tissue types
    """
    data_dirpath = Path(data_dirpath)

    adata = load_ciona_scrnaseq_data(stage, data_dir=data_dirpath)

    umap_df = extract_umap_coordinates(adata)
    umap_df = _merge_cluster_annotations(
        umap_df,
        stage,
        data_dirpath / "cluster_annotations" / "Ciona_scRNAseq_cluster_annotations.tsv",
    )
    umap_df = _merge_tissue_types(umap_df, stage, annotation_data_dirpath)
    umap_df.sort_values(["top_cluster_tissue_type", "top_cluster_suffix"], inplace=True)

    all_ciona_ids = mapper.map_to_all(input_id, input_id_type)

    uniprot_id = all_ciona_ids[CionaIDTypes.NONREF_PROTEIN]
    ky_id = all_ciona_ids[CionaIDTypes.KY_ID]
    hgnc_gene_symbol = all_ciona_ids[CionaIDTypes.HGNC_GENE_SYMBOL]

    umap_df = _add_gene_expression(umap_df, adata, ky_id)

    fig = make_subplots(
        rows=2,
        cols=2,
        shared_xaxes="all",
        shared_yaxes="all",
        horizontal_spacing=0.03,
        vertical_spacing=0.03,
    )

    # Plot positions based on the specified plot order.
    plot_positions = ["top-left", "top-right", "bottom-left", "bottom-right"]

    for position, plot in zip(plot_positions, plot_order, strict=True):
        if plot == "seurat_clusters":
            append_seurat_clusters_trace(
                fig,
                umap_df,
                position,
                title="Seurat cluster",
            )
        elif plot == "expression":
            append_expression_trace(
                fig,
                umap_df,
                position,
                title=f"<i>{hgnc_gene_symbol}</i> ({ky_id})",
                null_expression_color=null_expression_color,
                expression_colormap=expression_colormap,
                expression_colorbar_params=expression_colorbar_params,
            )
        elif plot == "tissue_type":
            append_tissue_type_trace(fig, umap_df, position, title="Cluster tissue type")
        elif plot == "cell_tissue_type":
            append_cell_tissue_type_trace(fig, umap_df, position, title="Cell tissue type")

    # Add UMAP axis labels.
    fig.add_annotation(
        text="UMAP1",
        x=0,
        y=0,
        xref="x3 domain",
        yref="y3 domain",
        xanchor="left",
        yanchor="top",
        showarrow=False,
        font=dict(
            size=10,
            color=apc.steel,
        ),
    )
    fig.add_annotation(
        text="UMAP2",
        x=0,
        y=0,
        xref="x3 domain",
        yref="y3 domain",
        xanchor="right",
        yanchor="bottom",
        textangle=-90,
        showarrow=False,
        font=dict(
            size=10,
            color=apc.steel,
        ),
    )

    # Add title annotation.
    fig.add_annotation(
        text=f"Single cell expression of <i>{hgnc_gene_symbol}</i> ({uniprot_id}) at {stage}",
        x=0,
        y=1.05,
        xref="paper",
        yref="paper",
        xanchor="left",
        yanchor="top",
        showarrow=False,
        font=dict(
            size=16,
            family=PLOTLY_TITLE_FONT,
        ),
    )

    _append_tissue_type_legend(
        fig,
        legend_x=1.015,
        legend_y=0.47,
        y_offset=0.024,
        title_offset=0.03,
    )

    fig.update_layout(width=width, height=height, margin=dict(l=15, r=150))
    fig.update_traces(marker=dict(size=marker_size))

    apc.plotly.hide_xaxis_ticks(fig)
    apc.plotly.hide_yaxis_ticks(fig)

    if image_filepath:
        os.makedirs(os.path.dirname(image_filepath), exist_ok=True)
        fig.write_image(image_filepath)

    if html_filepath:
        os.makedirs(os.path.dirname(html_filepath), exist_ok=True)
        fig.write_html(html_filepath, config=create_save_fig_config(width=width, height=height))

    return fig


def plot_expression_umap_for_all_stages(
    input_id: str,
    input_id_type: CionaIDTypes,
    mapper: IdentifierMapper,
    data_dirpath: str | Path = PIEKARZ_DATA_DIRPATH,
    annotation_data_dirpath: str | Path = CAO_DATA_DIRPATH,
    output_dirpath: str | Path = None,
) -> go.Figure:
    # Get all IDs using the provided mapper
    all_ciona_ids = mapper.map_to_all(input_id, input_id_type)
    uniprot_id = all_ciona_ids[CionaIDTypes.NONREF_PROTEIN]

    output_dirpath = Path(output_dirpath)

    for stage in tqdm(CionaStage.ordered_stages()):
        index = CionaStage.ordered_stages().index(stage) + 1
        index_stage = f"{index}_{stage}"

        plot_expression_umap(
            stage,
            input_id,
            input_id_type,
            mapper,
            data_dirpath,
            annotation_data_dirpath,
            image_filepath=output_dirpath
            / f"{input_id}_{uniprot_id}_expression/"
            / "umap"
            / f"{index_stage}_{input_id}_umap.svg",
            html_filepath=output_dirpath
            / f"{input_id}_{uniprot_id}_expression/"
            / "umap"
            / f"{index_stage}_{input_id}_umap.html",
        )


PUB_COLORBAR_PARAMS = dict(
    outlinewidth=0,
    thickness=10,
    ypad=0,
    len=0.8,
    y=1.05,
    x=1,
    xref="paper",
    yref="paper",
    yanchor="bottom",
    xanchor="right",
    tickfont=dict(
        size=12,
        family=apc.style_defaults.MONOSPACE_FONT_PLOTLY,
    ),
    orientation="h",
)


def plot_single_expression_umap(
    stage: CionaStage,
    input_id: str,
    input_id_type: CionaIDTypes,
    mapper: IdentifierMapper,
    data_dirpath: str | Path = PIEKARZ_DATA_DIRPATH,
    annotation_data_dirpath: str | Path = CAO_DATA_DIRPATH,
    plot_type: Literal[
        "seurat_clusters", "expression", "tissue_type", "cell_tissue_type"
    ] = "seurat_clusters",
    expression_colormap: apc.Gradient = GRAY_GRADIENT,
    null_expression_color: str = apc.gray,
    expression_colorbar_params: dict = DEFAULT_COLORBAR_PARAMS,
    marker_size: int = 3,
    width: int = 200,
    height: int = 300,
    image_filepath: str | Path | None = None,
    html_filepath: str | Path | None = None,
) -> go.Figure:
    """Plot UMAP visualizations of gene expression and cell type annotations.

    Args:
        stage: The developmental stage to plot
        input_id: The gene identifier to plot
        input_id_type: The type of the input identifier
        data_dirpath: Path to the Piekarz scRNA-seq data directory
        annotation_data_dirpath: Path to the Cao annotation data directory
        mapper: IdentifierMapper instance for ID conversion
        plot_type: The type of plot to generate.
        expression_colormap: The colormap to use for the expression plot.
        null_expression_color: The color to use for cells with no expression.
        expression_colorbar_params: The parameters to use for the expression colorbar.
        marker_size: The size of the markers in the plot.
        width: The width of the plot.

    Returns:
        A plotly Figure object containing a single UMAP plot.
    """
    data_dirpath = Path(data_dirpath)

    adata = load_ciona_scrnaseq_data(stage, data_dir=data_dirpath)

    umap_df = extract_umap_coordinates(adata)
    umap_df = _merge_cluster_annotations(
        umap_df,
        stage,
        data_dirpath / "cluster_annotations" / "Ciona_scRNAseq_cluster_annotations.tsv",
    )
    umap_df = _merge_tissue_types(umap_df, stage, annotation_data_dirpath)
    umap_df.sort_values(["top_cluster_tissue_type", "top_cluster_suffix"], inplace=True)

    all_ciona_ids = mapper.map_to_all(input_id, input_id_type)

    ky_id = all_ciona_ids[CionaIDTypes.KY_ID]
    hgnc_gene_symbol = all_ciona_ids[CionaIDTypes.HGNC_GENE_SYMBOL]

    umap_df = _add_gene_expression(umap_df, adata, ky_id)

    fig = go.Figure()

    position = None
    t_margin = 40
    r_margin = 10
    title_offset = 1.1

    if plot_type == "seurat_clusters":
        append_seurat_clusters_trace(
            fig,
            umap_df,
            position,
            title="Seurat cluster",
        )
    elif plot_type == "expression":
        append_expression_trace(
            fig,
            umap_df,
            position,
            title=f"<i>{hgnc_gene_symbol}</i> ({ky_id})",
            expression_colormap=expression_colormap,
            null_expression_color=null_expression_color,
            expression_colorbar_params=expression_colorbar_params,
        )
        t_margin += 100
        height += 100
        title_offset += 0.25
    elif plot_type == "tissue_type":
        append_tissue_type_trace(fig, umap_df, position, title="Cluster tissue type")
        _append_tissue_type_legend(
            fig,
            legend_x=1.015,
            legend_y=1,
            y_offset=0.09,
            title_offset=0.03,
        )
        r_margin += 100
        width += 100
    elif plot_type == "cell_tissue_type":
        append_cell_tissue_type_trace(fig, umap_df, position, title="Cell tissue type")
        _append_tissue_type_legend(
            fig,
            legend_x=1.015,
            legend_y=1,
            y_offset=0.09,
            title_offset=0.03,
        )
        r_margin += 100
        width += 100

    # Add UMAP axis labels.
    fig.add_annotation(
        text="UMAP1",
        x=0,
        y=0,
        xref="x domain",
        yref="y domain",
        xanchor="left",
        yanchor="top",
        showarrow=False,
        font=dict(
            size=10,
            color=apc.steel,
        ),
    )
    fig.add_annotation(
        text="UMAP2",
        x=0,
        y=0,
        xref="x domain",
        yref="y domain",
        xanchor="right",
        yanchor="bottom",
        textangle=-90,
        showarrow=False,
        font=dict(
            size=10,
            color=apc.steel,
        ),
    )

    fig.add_annotation(
        text=f"<i>{hgnc_gene_symbol}</i> at {stage}",
        x=0.2,
        y=title_offset,
        xref="paper",
        yref="paper",
        xanchor="left",
        yanchor="top",
        showarrow=False,
        font=dict(
            size=16,
            family=PLOTLY_TITLE_FONT,
        ),
    )

    fig.update_layout(
        width=width,
        height=height,
        margin=dict(
            l=15,
            r=r_margin,
            t=t_margin,
        ),
    )
    fig.update_traces(marker=dict(size=marker_size))

    apc.plotly.hide_xaxis_ticks(fig)
    apc.plotly.hide_yaxis_ticks(fig)

    if image_filepath:
        os.makedirs(os.path.dirname(image_filepath), exist_ok=True)
        fig.write_image(image_filepath)

    if html_filepath:
        os.makedirs(os.path.dirname(html_filepath), exist_ok=True)
        fig.write_html(html_filepath, config=create_save_fig_config(width=width, height=height))

    return fig


def get_tissue_expression_data(
    input_id,
    input_id_type,
    data_dirpath=PIEKARZ_DATA_DIRPATH,
    annotation_data_dirpath=CAO_DATA_DIRPATH,
):
    """Get tissue-specific expression data for a given gene across developmental stages.

    Args:
        input_id (str): Gene identifier
        input_id_type (CionaIDTypes): Type of input identifier
        data_dirpath (Path): Path to Piekarz data directory
        annotation_data_dirpath (Path): Path to Cao data directory

    Returns:
        pd.DataFrame
    """
    mapper = IdentifierMapper()
    data_collector = pd.DataFrame()

    for stage in CionaStage.ordered_stages():
        adata = load_ciona_scrnaseq_data(stage, data_dir=data_dirpath)

        umap_df = extract_umap_coordinates(adata)
        umap_df = _merge_cluster_annotations(
            umap_df,
            stage,
            data_dirpath / "cluster_annotations" / "Ciona_scRNAseq_cluster_annotations.tsv",
        )
        umap_df = _merge_tissue_types(umap_df, stage, annotation_data_dirpath)
        umap_df.sort_values(["top_cluster_tissue_type", "top_cluster_suffix"], inplace=True)

        all_ciona_ids = mapper.map_to_all(input_id, input_id_type)

        ky_id = all_ciona_ids[CionaIDTypes.KY_ID]

        umap_df = _add_gene_expression(umap_df, adata, ky_id)

        with warnings.catch_warnings():
            # Sometimes this will throw a warning if a row doesn't show expression,
            # but the analysis will run just fine.
            warnings.simplefilter("ignore", category=RuntimeWarning)
            expression_by_tissue = (
                umap_df.groupby("Tissue Type")
                .agg(
                    mean_expression=(
                        "gene_expression",
                        lambda x: x[x > 0].mean() if (x > 0).any() else 0,
                    ),
                    n_expressing_cells=(
                        "gene_expression",
                        lambda x: x[x > 0].sum(),
                    ),
                    n_cells=("gene_expression", "count"),
                )
                .reset_index()
            )
            expression_by_tissue["stage"] = stage
            expression_by_tissue["percent_expressing"] = (
                expression_by_tissue["n_expressing_cells"] / expression_by_tissue["n_cells"] * 100
            )

        data_collector = pd.concat([data_collector, expression_by_tissue])

    return data_collector


def append_total_cell_bubbles(
    fig: go.Figure,
    data_collector: pd.DataFrame,
    size_scaling_function=lambda x: np.round(np.sqrt(x), 2),
):
    for tissue_type in TISSUE_TYPE_PALETTE.keys():
        tissue_rows = data_collector[data_collector["Tissue Type"] == tissue_type]
        trace = go.Scatter(
            x=tissue_rows["Tissue Type"],
            y=tissue_rows["stage"],
            mode="markers",
            marker=dict(
                color=apc.gray,
                size=tissue_rows["n_cells"].apply(size_scaling_function),
                opacity=0.3,
                line=dict(color=TISSUE_TYPE_PALETTE[tissue_type], width=1.5),
            ),
            zorder=1,
            showlegend=False,
            name="Total cells",
            hoverinfo="skip",
        )
        fig.add_trace(trace)


BUBBLE_PLOT_HOVERTEMPLATE = (
    "All <b>%{customdata[0]} cells</b><br>"
    + "at stage <b>%{customdata[1]}</b><br>"
    + "<b>% expressing:</b> %{customdata[2]:.2f} (%{customdata[3]} / %{customdata[4]} cells)<br>"
    + "<b>Mean expression:</b> %{customdata[5]:.2f}"
)

BUBBLE_PLOT_CUSTOMDATA = [
    "Tissue Type",
    "stage",
    "percent_expressing",
    "n_expressing_cells",
    "n_cells",
    "mean_expression",
]


def append_expressing_cell_bubbles(
    fig: go.Figure,
    input_id: str,
    data_collector: pd.DataFrame,
    size_scaling_function=lambda x: np.round(np.sqrt(x), 2),
):
    trace = go.Scatter(
        x=data_collector["Tissue Type"],
        y=data_collector["stage"],
        mode="markers",
        marker=dict(
            color=data_collector["mean_expression"],
            size=data_collector["n_expressing_cells"].apply(size_scaling_function),
            colorscale=GRAY_GRADIENT.to_plotly_colorscale(),
            opacity=1,
            colorbar=dict(
                title=f"<i>{input_id}</i>",
                title_font=dict(
                    size=16,
                    family=apc.style_defaults.DEFAULT_FONT_PLOTLY + "-SemiBold",
                ),
                tickfont=dict(family=apc.style_defaults.MONOSPACE_FONT_PLOTLY),
                thickness=15,
                outlinewidth=0,
            ),
        ),
        zorder=2,
        showlegend=False,
        hovertemplate=BUBBLE_PLOT_HOVERTEMPLATE,
        customdata=data_collector[BUBBLE_PLOT_CUSTOMDATA],
    )

    fig.add_trace(trace)


def _add_bubble_grid(
    fig: go.Figure,
    num_cols: int = len(TISSUE_TYPE_PALETTE),
    num_rows: int = len(CionaStage.ordered_stages()),
):
    fig.update_xaxes(
        tickson="boundaries",
        showgrid=True,
        gridcolor=apc.gray,
        scaleanchor="y",
        scaleratio=1,
        tickangle=-45,
        range=[-0.9, num_cols - 0.3],
        side="top",
    )
    fig.update_yaxes(
        autorange="reversed",
        tickson="boundaries",
        showgrid=True,
        gridcolor=apc.gray,
        scaleanchor="x",
        scaleratio=1,
        range=[-0.3, num_rows + 1],
    )


def _add_bubble_legend(
    fig: go.Figure,
    num_cols: int = len(TISSUE_TYPE_PALETTE),
    num_rows: int = len(CionaStage.ordered_stages()),
    back_circle_size=0.7,
    front_circle_size=0.2,
    front_circle_color=apc.steel,
):
    # Create white background rectangle for legend.
    fig.add_shape(
        type="rect",
        x0=num_cols - 2.5,
        y0=num_rows - 0.5,
        x1=num_cols - 0.5,
        y1=num_rows + 1,
        fillcolor=apc.white,
        line=dict(width=0),
    )

    # Add text annotations.
    fig.add_annotation(
        xanchor="center",
        xref="x",
        yref="y",
        x=num_cols - 1,
        y=num_rows + 0.5,
        text="Total cells",
        showarrow=False,
        font=dict(family=apc.style_defaults.DEFAULT_FONT_PLOTLY, size=10, color=apc.cloud),
    )
    fig.add_annotation(
        xanchor="right",
        yanchor="middle",
        align="right",
        xref="x",
        yref="y",
        x=num_cols - 1.2,
        y=num_rows,
        text="Expressing<br>cells",
        showarrow=False,
        font=dict(
            family=apc.style_defaults.DEFAULT_FONT_PLOTLY + "-SemiBold",
            size=10,
            color=front_circle_color,
        ),
    )

    fig.add_annotation(
        xanchor="left",
        yanchor="middle",
        align="left",
        xref="x",
        yref="y",
        x=num_cols - 2.4,
        y=num_rows + 0.5,
        text="Key",
        font=dict(
            family=apc.style_defaults.DEFAULT_FONT_PLOTLY + "-SemiBold",
            size=16,
        ),
        showarrow=False,
    )

    # Calculate how much inset relative to the square grid
    # the background circle needs to be, based on its size.
    back_circle_offset = (1 - back_circle_size) / 2

    back_circle_x0 = num_cols - 1.5 + back_circle_offset
    back_circle_y0 = num_rows - 0.5 + back_circle_offset

    back_circle_x1 = back_circle_x0 + back_circle_size
    back_circle_y1 = back_circle_y0 + back_circle_size

    fig.add_shape(
        type="circle",
        xref="x",
        yref="y",
        x0=back_circle_x0,
        y0=back_circle_y0,
        x1=back_circle_x1,
        y1=back_circle_y1,
        fillcolor=apc.gray,
        line=dict(
            color=apc.marine,
            width=1.5,
        ),
        opacity=0.2,
    )

    # Obtain the center position of the back circle; then subtract
    # the radius of the front circle.
    front_circle_x0 = (back_circle_x1 + back_circle_x0) / 2 - front_circle_size / 2
    front_circle_y0 = (back_circle_y1 + back_circle_y0) / 2 - front_circle_size / 2

    front_circle_x1 = front_circle_x0 + front_circle_size
    front_circle_y1 = front_circle_y0 + front_circle_size

    fig.add_shape(
        type="circle",
        xref="x",
        yref="y",
        x0=front_circle_x0,
        y0=front_circle_y0,
        x1=front_circle_x1,
        y1=front_circle_y1,
        fillcolor=front_circle_color,
        line=dict(
            color=apc.white,
            width=1,
        ),
    )


def plot_expression_bubbles(
    data_collector,
    input_id,
    width=620,
    height=800,
    size_scaling_function=lambda x: np.round(np.sqrt(x), 2) * 1.1,
    image_filepath: str | Path = None,
    html_filepath: str | Path = None,
):
    fig = go.Figure()

    append_total_cell_bubbles(fig, data_collector, size_scaling_function)
    append_expressing_cell_bubbles(fig, input_id, data_collector, size_scaling_function)

    fig.update_layout(
        height=height,
        width=width,
    )

    _add_bubble_grid(fig)
    _add_bubble_legend(fig)

    apc.plotly.set_axes_categorical(fig)

    if image_filepath:
        os.makedirs(os.path.dirname(image_filepath), exist_ok=True)
        fig.write_image(image_filepath)

    if html_filepath:
        os.makedirs(os.path.dirname(html_filepath), exist_ok=True)
        fig.write_html(html_filepath, config=create_save_fig_config(width=width, height=height))

    return fig


# These values correspond to the numerical developmental stage (e.g. S1, or 1 cell),
# as described in https://anatomypubs.onlinelibrary.wiley.com/doi/10.1002/dvdy.21188
BULK_RNA_SEQ_STAGES = [1, 8, 11, 12, 15, 21, 26]
BULK_RNA_SEQ_STAGES_AS_STRINGS = [f"S{stage}" for stage in BULK_RNA_SEQ_STAGES]

# These values corespond to the percent development at each stage,
# as described in https://anatomypubs.onlinelibrary.wiley.com/doi/10.1002/dvdy.21188
BULK_RNA_SEQ_PERCENT_DEVELOPMENT = [3, 23, 28, 32, 39, 57, 100]
BULK_RNA_SEQ_REPLICATES = [1, 2]

REPLICATE_COLORS = {"1": apc.sun, "2": apc.mustard}


class BulkRNASeqDataTypes(StrEnum):
    FPKM = "FPKM"  # FPKM: Fragments Per Kilobase Million.
    RLE = "RLE"  # Relative Log Expression.


def _extract_bulk_rna_seq_data(
    kh_id: str,
    data_dirpath: str | Path = BULK_RNA_SEQ_DATA_DIRPATH,
):
    data = pd.DataFrame()

    data_dirpath = Path(data_dirpath)

    for stage in BULK_RNA_SEQ_STAGES:
        for replicate in BULK_RNA_SEQ_REPLICATES:
            full_data = pd.read_csv(
                data_dirpath / f"Cirobu_RNA-Seq_stage{stage}_WT_Replicate{replicate}.csv",
                sep=";",
            )
            full_data.rename(columns={"Unnamed: 0": "kh_id"}, inplace=True)
            full_data["stage"] = f"S{stage}"
            full_data["percent_development"] = BULK_RNA_SEQ_PERCENT_DEVELOPMENT[
                BULK_RNA_SEQ_STAGES.index(stage)
            ]
            full_data["replicate"] = str(replicate)
            full_data = full_data.loc[full_data["kh_id"] == kh_id]

            data = pd.concat([data, full_data])

    data["stage"] = data["stage"].astype("category")
    data["stage"] = data["stage"].cat.reorder_categories(BULK_RNA_SEQ_STAGES_AS_STRINGS)
    data.sort_values(by="stage", inplace=True)

    return data


def plot_bulk_rna_seq_expression(
    input_id,
    input_id_type,
    mapper: IdentifierMapper,
    data_dirpath: str | Path = BULK_RNA_SEQ_DATA_DIRPATH,
    datatype=BulkRNASeqDataTypes.FPKM,
    width=400,
    height=250,
    image_filepath: str | Path | None = None,
    html_filepath: str | Path | None = None,
    spacing: Literal["stage", "percent_development"] = "stage",
    adjust_ylimits=False,
    override_kh_id: str | None = None,
    override_gene_symbol: str | None = None,
    replicate_colors: dict[str, str] = REPLICATE_COLORS,
):
    """Plot bulk RNA-seq expression data for a given KH ID across developmental stages.

    Args:
        input_id (str): Input ID of gene to plot
        input_id_type (CionaIDTypes): Type of input ID
        mapper (IdentifierMapper): Identifier mapper
        data_dirpath (str | Path): Path to data directory
        datatype (BulkRNASeqDataTypes): Expression data type to plot
        width (int): Plot width
        height (int): Plot height
        image_filepath (str | Path): Path to save image
        html_filepath (str | Path): Path to save HTML
        spacing (Literal["uniform", "scaled"]): Spacing of x-axis
        adjust_ylimits (bool): Whether to adjust y-limits
        override_kh_id (str | None): Override KH ID
        override_gene_symbol (str | None): Override gene symbol
        replicate_colors (dict[str, str]): Colors for replicates

    Returns:
        plotly.graph_objects.Figure: Expression plot
    """
    kh_id = input_id
    hgnc_gene_symbol = input_id

    try:
        all_ciona_ids = mapper.map_to_all(input_id, input_id_type)
        if not pd.isna(all_ciona_ids[CionaIDTypes.KH_ID]):
            kh_id = "KH2012:" + all_ciona_ids[CionaIDTypes.KH_ID]
            hgnc_gene_symbol = all_ciona_ids[CionaIDTypes.HGNC_GENE_SYMBOL]
        else:
            if override_kh_id is None:
                print(f"No KH ID found for {input_id}. Try override_kh_id.")
    except ValueError:
        pass

    if override_kh_id is not None:
        kh_id = "KH2012:" + override_kh_id

    if override_gene_symbol is not None:
        hgnc_gene_symbol = override_gene_symbol

    data = pd.DataFrame()

    data_dirpath = Path(data_dirpath)
    data = _extract_bulk_rna_seq_data(kh_id, data_dirpath)

    fig = px.line(
        data,
        x=spacing,
        y=datatype,
        color="replicate",
        color_discrete_map=replicate_colors,
        markers=True,
    )

    fig.update_layout(
        xaxis_title=spacing.title().replace("_", " "),
        yaxis_title=datatype,
        width=width,
        height=height,
        legend=dict(yanchor="top", y=1.02, xanchor="left", x=1),
        margin=dict(t=50),
    )

    title = f"Expression of {hgnc_gene_symbol} ({kh_id})<br>in bulk RNA-Seq"

    fig.add_annotation(
        text=title,
        xref="paper",
        yref="paper",
        x=0.5,
        xanchor="center",
        y=1.4,
        showarrow=False,
        align="center",
        font=dict(
            size=15,
            family=PLOTLY_TITLE_FONT,
        ),
    )

    if adjust_ylimits:
        fig.update_yaxes(tickfont=dict(size=10), range=[0, data[datatype].max() * 1.1])

    if spacing == "percent_development":
        fig.update_xaxes(
            ticktext=[f"S{i}" for i in BULK_RNA_SEQ_STAGES],
            tickvals=BULK_RNA_SEQ_PERCENT_DEVELOPMENT,
            tickfont=dict(size=10),
            tickangle=45,
        )
    elif spacing == "stage":
        pass

    if image_filepath is not None:
        fig.write_image(image_filepath)

    if html_filepath is not None:
        fig.write_html(html_filepath, config=create_save_fig_config(width, height))

    return fig


def _append_range_traces(
    stage_summaries: pd.DataFrame,
    fig: go.Figure,
):
    # Add a line for the maximum expression.
    trace0 = go.Scatter(
        x=stage_summaries["stage"],
        y=stage_summaries["stage_max"],
        mode="lines",
        line_color=apc.oat,
        showlegend=False,
    )
    fig.add_trace(trace0)

    # Add a line for the minimum expression, filled between the maximum and minimum.
    trace1 = go.Scatter(
        x=stage_summaries["stage"],
        y=stage_summaries["stage_min"],
        fill="tonexty",  # fill area between trace0 and trace1
        mode="lines",
        line_color=apc.oat,
        showlegend=False,
    )
    fig.add_trace(trace1)

    # Add a line for the mean expression.
    trace2 = go.Scatter(
        x=stage_summaries["stage"],
        y=stage_summaries["stage_mean"],
        mode="lines",
        line=dict(color=apc.canary),
        name="Mean",
    )
    fig.add_trace(trace2)


def _append_replicate_traces(
    data: pd.DataFrame,
    spacing: Literal["stage", "percent_development"],
    fig: go.Figure,
    replicate_colors: dict[str, str] = REPLICATE_COLORS,
    datatype: BulkRNASeqDataTypes = BulkRNASeqDataTypes.FPKM,
):
    for replicate in BULK_RNA_SEQ_REPLICATES:
        replicate_data = data[data["replicate"] == str(replicate)]

        trace1 = go.Scatter(
            x=replicate_data[spacing],
            y=replicate_data[datatype],
            marker=dict(color=replicate_colors[str(replicate)]),
            mode="markers",
            name=f"Replicate {replicate}",
        )
        fig.add_trace(trace1)


def plot_bulk_rna_seq_expression_range(
    input_id,
    input_id_type,
    mapper: IdentifierMapper,
    data_dirpath: str | Path = BULK_RNA_SEQ_DATA_DIRPATH,
    datatype=BulkRNASeqDataTypes.FPKM,
    width=450,
    height=250,
    image_filepath: str | Path | None = None,
    html_filepath: str | Path | None = None,
    spacing: Literal["stage", "percent_development"] = "stage",
    override_kh_id: str | None = None,
    override_gene_symbol: str | None = None,
    replicate_colors: dict[str, str] = REPLICATE_COLORS,
    tickangle: int = 45,
):
    """Plot bulk RNA-seq expression data for a given KH ID across developmental stages.

    Args:
        input_id (str): Input ID of gene to plot
        input_id_type (CionaIDTypes): Type of input ID
        mapper (IdentifierMapper): Identifier mapper
        data_dirpath (str | Path): Path to data directory
        datatype (BulkRNASeqDataTypes): Expression data type to plot
        width (int): Plot width
        height (int): Plot height
        image_filepath (str | Path): Path to save image
        html_filepath (str | Path): Path to save HTML
        spacing (Literal["uniform", "scaled"]): Spacing of x-axis
        adjust_ylimits (bool): Whether to adjust y-limits
        override_kh_id (str | None): Override KH ID
        override_gene_symbol (str | None): Override gene symbol

    Returns:
        plotly.graph_objects.Figure: Expression plot
    """
    kh_id = input_id
    hgnc_gene_symbol = input_id

    try:
        all_ciona_ids = mapper.map_to_all(input_id, input_id_type)
        if not pd.isna(all_ciona_ids[CionaIDTypes.KH_ID]):
            kh_id = "KH2012:" + all_ciona_ids[CionaIDTypes.KH_ID]
            hgnc_gene_symbol = all_ciona_ids[CionaIDTypes.HGNC_GENE_SYMBOL]
        else:
            if override_kh_id is None:
                print(f"No KH ID found for {input_id}. Try override_kh_id.")
    except ValueError:
        pass

    if override_kh_id is not None:
        kh_id = "KH2012:" + override_kh_id

    if override_gene_symbol is not None:
        hgnc_gene_symbol = override_gene_symbol

    data = _extract_bulk_rna_seq_data(kh_id, data_dirpath)

    fig = go.Figure()

    stage_summaries = (
        data.groupby("stage", observed=True)
        .agg(
            stage_max=(datatype, "max"),
            stage_min=(datatype, "min"),
            stage_mean=(datatype, "mean"),
        )
        .reset_index()
    )

    _append_range_traces(stage_summaries, fig)
    _append_replicate_traces(data, spacing, fig, replicate_colors, datatype)

    fig.update_layout(
        xaxis_title=spacing.title().replace("_", " "),
        yaxis_title=datatype,
        width=width,
        height=height,
        legend=dict(yanchor="top", y=1.02, xanchor="left", x=1),
        margin=dict(t=50),
    )

    title = f"Expression of {hgnc_gene_symbol} ({kh_id})<br>in bulk RNA-Seq"

    fig.add_annotation(
        text=title,
        xref="paper",
        yref="paper",
        x=0.5,
        xanchor="center",
        y=1.4,
        showarrow=False,
        align="center",
        font=dict(
            size=15,
            family=PLOTLY_TITLE_FONT,
        ),
    )

    fig.update_layout(legend_traceorder="normal")
    fig.update_xaxes(tickangle=tickangle)

    apc.plotly.set_yticklabel_monospaced(fig)
    apc.plotly.set_xticklabel_monospaced(fig)

    if image_filepath is not None:
        fig.write_image(image_filepath)

    if html_filepath is not None:
        fig.write_html(html_filepath, config=create_save_fig_config(width, height))

    return fig
