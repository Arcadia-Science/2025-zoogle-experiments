import os
from enum import StrEnum
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
from zoogletools.ciona.constants import (
    CIONA_GENE_MODELS_DIRPATH,
    PIEKARZ_DATA_DIRPATH,
    ZOOGLE_RESULTS_DIRPATH,
    CionaStage,
)
from zoogletools.ciona.data_processing import load_ciona_scrnaseq_data
from zoogletools.plotting import create_save_fig_config


class CionaIDTypes(StrEnum):
    HGNC_GENE_SYMBOL = "hgnc_gene_symbol"
    NONREF_PROTEIN = "nonref_protein"
    KY_ID = "ky_id"
    KH_ID = "kh_id"


def extract_umap_coordinates(adata: sc.AnnData) -> pd.DataFrame:
    """
    Extract the UMAP coordinates from the scRNAseq data.
    """
    umap_embeddings = adata.obsm["X_umap"]
    umap_df = pd.DataFrame(umap_embeddings, columns=["UMAP_1", "UMAP_2"], index=adata.obs_names)
    umap_df = pd.concat([umap_df, adata.obs], axis=1)

    return umap_df


def _handle_multiple_hits(
    hits: pd.DataFrame,
    input_id_type: str,
    input_id: str,
    output_id_type: str,
) -> str:
    """Handle multiple hits for a given input ID.

    Args:
        hits: The hits dataframe.
        input_id_type: The type of the input ID.
        input_id: The input ID.

    Returns:
        The top hit.
    """
    if len(hits) == 0:
        raise ValueError(f"No hits found for {input_id_type} == {input_id}")
    elif len(hits) > 1:
        print(f"Multiple hits found for {input_id_type} == {input_id}")
        print(hits[[input_id_type, output_id_type]])
        top_hit = hits.iloc[0]
        print(f"Returning top hit: {top_hit[output_id_type]}")

    return hits[output_id_type].values[0]


def map_hgnc_to_ciona_uniprot(
    results_filepath: str,
    input_id: str,
    input_id_type: str = CionaIDTypes.HGNC_GENE_SYMBOL,
    output_id_type: str = CionaIDTypes.NONREF_PROTEIN,
):
    """Create a mapping between HGNC gene symbols and Salpingoeca UniProt IDs.

    Args:
        results_filepath: Path to Zoogle results TSV file containing HGNC gene symbols
            and Ciona protein IDs
        input_id: HGNC gene symbol or Ciona UniProt ID
        input_id_type: "hgnc_gene_symbol" or "nonref_protein"
        output_id_type: "hgnc_gene_symbol" or "nonref_protein"

    Returns:
        HGNC gene symbol or Ciona UniProt ID
    """
    zoogle_results = pd.read_csv(results_filepath, sep="\t")

    hits = zoogle_results[zoogle_results[input_id_type] == input_id]
    hits = hits.sort_values(by=["trait_dist"], ascending=False)

    return _handle_multiple_hits(
        hits=hits, input_id_type=input_id_type, input_id=input_id, output_id_type=output_id_type
    )


def map_uniprot_to_ky_id(
    input_filepath: str,
    input_id: str,
    input_id_type: str = CionaIDTypes.NONREF_PROTEIN,
    output_id_type: str = CionaIDTypes.KY_ID,
) -> dict:
    """Create a mapping between HGNC gene symbols and Ciona UniProt IDs.

    Args:
        input_filepath: Path to input TSV file containing HGNC gene symbols
            and Ciona protein IDs
        input_id: Ciona UniProt ID
        input_id_type: "nonref_protein"
        output_id_type: "ky_id"

    Returns:
        dict: Mapping from HGNC gene symbols to their corresponding Ciona UniProt IDs
    """
    input_df = pd.read_csv(input_filepath, sep="\t")

    hits = input_df[input_df[input_id_type] == input_id]

    return _handle_multiple_hits(
        hits=hits, input_id_type=input_id_type, input_id=input_id, output_id_type=output_id_type
    )


def map_ky_to_kh_id(
    input_filepath: str,
    input_id: str,
    input_id_type: str = CionaIDTypes.KY_ID,
    output_id_type: str = CionaIDTypes.KH_ID,
) -> str:
    """Map a Ciona KY ID to a Ciona KH ID.

    Args:
        input_filepath: Path to the input TSV file containing the KY ID to KH ID mapping.
        input_id: The KY ID to map.
        input_id_type: The type of the input ID.
        output_id_type: The type of the output ID.

    Returns:
        The KH ID.
    """
    input_df = pd.read_csv(input_filepath, sep="\t")

    hits = input_df[input_df[input_id_type] == input_id]

    return _handle_multiple_hits(
        hits=hits, input_id_type=input_id_type, input_id=input_id, output_id_type=output_id_type
    )


def map_input_id_to_all_ciona_ids(
    input_id: str,
    input_id_type: str = CionaIDTypes.HGNC_GENE_SYMBOL,
    zoogle_results_dirpath: str | Path = ZOOGLE_RESULTS_DIRPATH,
    ciona_gene_models_dirpath: str | Path = CIONA_GENE_MODELS_DIRPATH,
) -> dict:
    """Map an input ID to all Ciona IDs.

    Args:
        input_id: The input ID to map.
        input_id_type: The type of the input ID.
        zoogle_results_dirpath: The directory containing the Zoogle results.
        ciona_gene_models_dirpath: The directory containing the Ciona gene models.

    Returns:
        A dictionary containing the mapped IDs.
    """
    zoogle_results_filepath = (
        zoogle_results_dirpath / "per-nonref-species" / "Ciona-intestinalis.tsv"
    )

    ciona_uniprot_ky_map_filepath = ciona_gene_models_dirpath / "ciona_uniprot_ky_map.tsv"

    ciona_ky_kh_map_filepath = ciona_gene_models_dirpath / "ciona_ky_kh_map.tsv"

    if input_id_type == CionaIDTypes.HGNC_GENE_SYMBOL:
        hgnc_gene_symbol = input_id
        uniprot_id = map_hgnc_to_ciona_uniprot(
            results_filepath=zoogle_results_filepath,
            input_id=input_id,
            input_id_type=CionaIDTypes.HGNC_GENE_SYMBOL,
            output_id_type=CionaIDTypes.NONREF_PROTEIN,
        )
        ky_id = map_uniprot_to_ky_id(
            input_filepath=ciona_uniprot_ky_map_filepath,
            input_id=uniprot_id,
            input_id_type=CionaIDTypes.NONREF_PROTEIN,
            output_id_type=CionaIDTypes.KY_ID,
        )
        kh_id = map_ky_to_kh_id(
            input_filepath=ciona_ky_kh_map_filepath,
            input_id=ky_id,
            input_id_type=CionaIDTypes.KY_ID,
            output_id_type=CionaIDTypes.KH_ID,
        )
    elif input_id_type == CionaIDTypes.NONREF_PROTEIN:
        uniprot_id = input_id
        hgnc_gene_symbol = map_hgnc_to_ciona_uniprot(
            results_filepath=zoogle_results_filepath,
            input_id=uniprot_id,
            input_id_type=CionaIDTypes.NONREF_PROTEIN,
            output_id_type=CionaIDTypes.HGNC_GENE_SYMBOL,
        )
        ky_id = map_uniprot_to_ky_id(
            input_filepath=ciona_uniprot_ky_map_filepath,
            input_id=uniprot_id,
            input_id_type=CionaIDTypes.NONREF_PROTEIN,
            output_id_type=CionaIDTypes.KY_ID,
        )
        kh_id = map_ky_to_kh_id(
            input_filepath=ciona_ky_kh_map_filepath,
            input_id=ky_id,
            input_id_type=CionaIDTypes.KY_ID,
            output_id_type=CionaIDTypes.KH_ID,
        )
    elif input_id_type == CionaIDTypes.KY_ID:
        ky_id = input_id
        uniprot_id = map_uniprot_to_ky_id(
            input_filepath=ciona_uniprot_ky_map_filepath,
            input_id=ky_id,
            input_id_type=CionaIDTypes.KY_ID,
            output_id_type=CionaIDTypes.NONREF_PROTEIN,
        )
        hgnc_gene_symbol = map_hgnc_to_ciona_uniprot(
            results_filepath=zoogle_results_filepath,
            input_id=uniprot_id,
            input_id_type=CionaIDTypes.NONREF_PROTEIN,
            output_id_type=CionaIDTypes.HGNC_GENE_SYMBOL,
        )
        kh_id = map_ky_to_kh_id(
            input_filepath=ciona_ky_kh_map_filepath,
            input_id=ky_id,
            input_id_type=CionaIDTypes.KY_ID,
            output_id_type=CionaIDTypes.KH_ID,
        )
    elif input_id_type == CionaIDTypes.KH_ID:
        kh_id = input_id
        ky_id = map_ky_to_kh_id(
            input_filepath=ciona_ky_kh_map_filepath,
            input_id=kh_id,
            input_id_type=CionaIDTypes.KH_ID,
            output_id_type=CionaIDTypes.KY_ID,
        )
        uniprot_id = map_uniprot_to_ky_id(
            input_filepath=ciona_uniprot_ky_map_filepath,
            input_id=ky_id,
            input_id_type=CionaIDTypes.KY_ID,
            output_id_type=CionaIDTypes.NONREF_PROTEIN,
        )
        hgnc_gene_symbol = map_hgnc_to_ciona_uniprot(
            results_filepath=zoogle_results_filepath,
            input_id=uniprot_id,
            input_id_type=CionaIDTypes.NONREF_PROTEIN,
            output_id_type=CionaIDTypes.HGNC_GENE_SYMBOL,
        )
    else:
        raise ValueError(f"Invalid input_id_type: {input_id_type}")

    return {
        CionaIDTypes.HGNC_GENE_SYMBOL: hgnc_gene_symbol,
        CionaIDTypes.NONREF_PROTEIN: uniprot_id,
        CionaIDTypes.KY_ID: ky_id,
        CionaIDTypes.KH_ID: kh_id,
    }


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


def add_cell_count_barchart(
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
    input_id_type: str = CionaIDTypes.HGNC_GENE_SYMBOL,
    data_dirpath: str | Path = PIEKARZ_DATA_DIRPATH,
    zoogle_results_dirpath: str | Path = ZOOGLE_RESULTS_DIRPATH,
    ciona_gene_models_dirpath: str | Path = CIONA_GENE_MODELS_DIRPATH,
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

    all_ciona_ids = map_input_id_to_all_ciona_ids(
        input_id=input_id,
        input_id_type=input_id_type,
        zoogle_results_dirpath=zoogle_results_dirpath,
        ciona_gene_models_dirpath=ciona_gene_models_dirpath,
    )

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

    add_cell_count_barchart(
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
    input_id_type: str = CionaIDTypes.HGNC_GENE_SYMBOL,
    data_dirpath: str | Path = PIEKARZ_DATA_DIRPATH,
    zoogle_results_dirpath: str | Path = ZOOGLE_RESULTS_DIRPATH,
    ciona_gene_models_dirpath: str | Path = CIONA_GENE_MODELS_DIRPATH,
    color_mode: Literal["tissue", "cluster"] = "tissue",
    output_dirpath: str | Path = None,
) -> go.Figure:
    all_ciona_ids = map_input_id_to_all_ciona_ids(
        input_id=input_id,
        input_id_type=input_id_type,
        zoogle_results_dirpath=zoogle_results_dirpath,
        ciona_gene_models_dirpath=ciona_gene_models_dirpath,
    )

    uniprot_id = all_ciona_ids[CionaIDTypes.NONREF_PROTEIN]

    output_dirpath = Path(output_dirpath)

    for stage in tqdm(CionaStage.ordered_stages()):
        index = CionaStage.ordered_stages().index(stage) + 1
        index_stage = f"{index}_{stage}"

        plot_expression_violin(
            stage=stage,
            input_id=input_id,
            input_id_type=input_id_type,
            data_dirpath=data_dirpath,
            zoogle_results_dirpath=zoogle_results_dirpath,
            ciona_gene_models_dirpath=ciona_gene_models_dirpath,
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
