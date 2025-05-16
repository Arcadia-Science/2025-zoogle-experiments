import os
from pathlib import Path
from typing import Literal

import arcadia_pycolor as apc
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import scanpy as sc
from plotly.subplots import make_subplots
from tqdm import tqdm

from zoogletools.ciona.cluster_annotations import TISSUE_PALETTE_DICT
from zoogletools.ciona.constants import PIEKARZ_DATA_DIRPATH, CionaStage
from zoogletools.ciona.data_processing import load_ciona_scrnaseq_data
from zoogletools.ciona.identifier_mapping import CionaIDTypes, IdentifierMapper
from zoogletools.plotting import create_save_fig_config


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
        categories=TISSUE_PALETTE_DICT.keys(),
        ordered=True,
    )

    return umap_df.merge(cluster_annotations, on="seurat_clusters", how="left")


def append_cell_count_barchart(
    fig: go.Figure,
    ky_id: str,
    umap_df: pd.DataFrame,
    cluster_color_map: dict,
    row: int,
    col: int,
) -> go.Figure:
    cluster_stats = umap_df.groupby("seurat_clusters", observed=True).agg(
        total_cells=(ky_id, "count"),
        expressing_cells=(ky_id, lambda x: (x > 0).sum()),
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
    ky_id: str,
    umap_df: pd.DataFrame,
    cluster_color_map: dict,
    row: int,
    col: int,
) -> go.Figure:
    cluster_stats = umap_df.groupby("seurat_clusters", observed=True).agg(
        total_cells=(ky_id, "count"),
        expressing_cells=(ky_id, lambda x: (x > 0).sum()),
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
    ky_id: str,
    umap_df: pd.DataFrame,
    cluster_color_map: dict,
    row: int,
    col: int,
) -> go.Figure:
    umap_df_nonzero = umap_df[umap_df[ky_id] > 0]

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
                    y=cluster_data[ky_id],
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

    gene_index = np.where(adata.var.index == ky_id)[0][0]
    gene_expression = adata.X[:, gene_index].toarray().flatten()
    umap_df[ky_id] = gene_expression

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
        ky_id=ky_id,
        umap_df=umap_df,
        cluster_color_map=cluster_color_map,
        row=1,
        col=1,
    )

    append_expression_proportion_barchart(
        fig=fig,
        ky_id=ky_id,
        umap_df=umap_df,
        cluster_color_map=cluster_color_map,
        row=2,
        col=1,
    )

    append_expression_violin(
        fig=fig,
        ky_id=ky_id,
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
