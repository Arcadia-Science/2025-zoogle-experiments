import itertools
from collections import defaultdict
from enum import StrEnum
from pathlib import Path

import arcadia_pycolor as apc
import numpy as np
import pandas as pd

import zoogletools as zt
from zoogletools.ciona.constants import CionaStage


def quantify_cluster_annotations(
    cell_clusters: pd.DataFrame,
    input_dirpath: Path = Path("../../data/Ciona_intestinalis_scRNAseq_data_Piekarz"),
    output_dirpath: Path = Path(
        "../../data/Ciona_intestinalis_scRNAseq_data_Piekarz/cluster_annotations"
    ),
) -> None:
    """Count cluster annotations from the Cao et al. data
        for each developmental stage in the Piekarz et al. data.

    Args:
        cell_clusters: DataFrame containing cell cluster annotations from Cao et al. 2019.
        scrnaseq_data_dir: Directory containing scRNAseq data from Piekarz et al 2024.
        output_dir: Directory to save cluster annotation files
    """
    for stage in CionaStage.ordered_stages():
        index = CionaStage.ordered_stages().index(stage) + 1

        adata = zt.ciona.data_processing.load_ciona_scrnaseq_data(stage, input_dirpath)

        replicates = adata.obs.rep.unique()
        replicate_cluster_annotations = pd.DataFrame()

        for replicate in replicates:
            cells_in_replicate = adata.obs[adata.obs["rep"] == replicate]
            annotations = cell_clusters[cell_clusters["stage_replicate"] == stage + "_" + replicate]
            annotations = annotations[["merging_barcode", "Tissue Type"]]

            cells_in_replicate = cells_in_replicate.merge(
                annotations, on="merging_barcode", how="left"
            )
            cells_in_replicate.replace(
                {
                    "Tissue Type": {
                        np.nan: "unannotated",
                        "muscle & heart": "muscle-heart",
                        "nervous system": "nervous-system",
                    }
                },
                inplace=True,
            )
            replicate_cluster_annotations = pd.concat(
                [
                    replicate_cluster_annotations,
                    cells_in_replicate,
                ]
            )

        stage_cluster_annotations = (
            replicate_cluster_annotations.groupby("seurat_clusters", observed=True)
            .agg(
                tissue_type_count=("Tissue Type", "value_counts"),
            )
            .reset_index(drop=False)
        )
        stage_cluster_annotations["stage"] = stage

        if not output_dirpath.exists():
            output_dirpath.mkdir(parents=True)

        stage_cluster_annotations.to_csv(
            output_dirpath / f"{index}_{stage}_cluster_annotations.tsv",
            sep="\t",
            index=False,
        )


def _get_top_cell_types(
    annotations: pd.DataFrame, n: int = 2, round_to: int = 2
) -> list[tuple[str, float]]:
    """Get the top n cell types and their fractions for a cluster.

    Args:
        annotations: DataFrame containing cell type annotations for a single cluster.
        n: Number of top cell types to return
        round_to: Number of decimal places to round the fractions

    Returns:
        List of tuples containing the top n cell types and their fractions
    """
    top_types = annotations.nlargest(n, "tissue_type_count")
    total_cells = annotations["tissue_type_count"].sum()
    return [
        (tissue, np.round(count / total_cells, round_to))
        for tissue, count in zip(
            top_types["Tissue Type"], top_types["tissue_type_count"], strict=True
        )
    ]


def _number_repeated_types(cluster_annotations: dict[str, str]) -> dict[str, int | float]:
    """Returns a dictionary of clusters and suffixes for cell types
        that appear in more than one cluster.

    This is based on the primary cell type of the cluster.

    For example, if there are two clusters with the annotation of "muscle",
    this function returns the suffixes "1" and "2" for those two clusters.
    If the cell type appears only once, the suffix is NaN.

    Args:
        cluster_annotations: Dictionary mapping clusters to their top cell types

    Returns:
        Dictionary mapping clusters to their suffixes
    """
    primary_types = cluster_annotations.values()
    type_counts = pd.Series(primary_types).value_counts()

    # This generates a defaultdict that counts the number of occurrences of each cell type.
    # The lambda function is used to initialize the counter for each cell type.
    # With each access, the counter is incremented by 1.
    type_counters = defaultdict(lambda: itertools.count(1))

    numbered_suffixes = defaultdict(lambda: np.nan)
    for cluster, annotation in cluster_annotations.items():
        # If the cell type appears in more than one cluster, add a number to it.
        if type_counts[annotation] > 1:
            # The number added to the cell type annotation increases
            # by 1 for each additional occurrence.
            numbered_suffixes[cluster] = next(type_counters[annotation])
        else:
            # If the cell type appears only once, keep it as is.
            numbered_suffixes[cluster] = np.nan

    return numbered_suffixes


class ClusterAnnotationsColumns(StrEnum):
    INDEX_STAGE = "index_stage"
    SEURAT_CLUSTERS = "seurat_clusters"
    TOP_CLUSTER_CELLTYPE = "top_cluster_celltype"
    TOP_CLUSTER_FRACTION = "top_cluster_fraction"
    SECOND_CLUSTER_CELLTYPE = "second_cluster_celltype"
    SECOND_CLUSTER_FRACTION = "second_cluster_fraction"
    TOP_CLUSTER_SUFFIX = "top_cluster_suffix"
    FORMATTED_CLUSTER_NAME = "formatted_cluster_name"
    CLUSTER_TISSUE_COLOR = "cluster_tissue_color"


def _format_cell_type_row(row: pd.Series) -> str:
    """Format a cell type annotation row into a string representation.

    If the top cluster fraction is more than double the second cluster fraction,
    only include the top cluster. Otherwise include both clusters joined by '+'.

    Adds numeric suffixes (e.g. '_1') for cell types that appear multiple times.
    """
    CAC = ClusterAnnotationsColumns

    # Format the top cluster (with a suffix if it appears multiple times).
    top_cluster = row[CAC.TOP_CLUSTER_CELLTYPE]
    if pd.notna(row[CAC.TOP_CLUSTER_SUFFIX]):
        top_cluster = f"{top_cluster}_{row[CAC.TOP_CLUSTER_SUFFIX]}"
    top_cluster = f"{top_cluster}({row[CAC.TOP_CLUSTER_FRACTION]})"

    # If second cluster is at least half as abundant, include it.
    if row[CAC.TOP_CLUSTER_FRACTION] <= 2 * row[CAC.SECOND_CLUSTER_FRACTION]:
        second_cluster = f"{row[CAC.SECOND_CLUSTER_CELLTYPE]}({row[CAC.SECOND_CLUSTER_FRACTION]})"
        return f"{top_cluster}+{second_cluster}"

    return top_cluster


def smooth_gradient_from_palette(palette: apc.Palette):
    return apc.Gradient(palette.name, palette.colors).interpolate_lightness()


TISSUE_PALETTE_DICT = {
    "epidermis": smooth_gradient_from_palette(apc.palettes.blue_shades),
    "nervous-system": smooth_gradient_from_palette(apc.palettes.purple_shades),
    "notochord": smooth_gradient_from_palette(apc.palettes.warm_gray_shades),
    "mesenchyme": smooth_gradient_from_palette(apc.palettes.red_shades),
    "muscle-heart": smooth_gradient_from_palette(apc.palettes.pink_shades),
    "endoderm": smooth_gradient_from_palette(apc.palettes.yellow_shades),
    "unannotated": smooth_gradient_from_palette(apc.palettes.cool_gray_shades),
    "germ": smooth_gradient_from_palette(apc.palettes.green_shades),
}


def _assign_cluster_tissue_colors(
    cluster_annotations: dict[str, str],
    tissue_palette: dict[str, apc.Gradient] = TISSUE_PALETTE_DICT,
) -> dict[str, str]:
    """Assign colors to clusters based on their top cell type.

    Args:
        cluster_annotations: Dictionary mapping clusters to their top cell types

    Returns:
        Dictionary mapping clusters to their colors
    """
    primary_types = cluster_annotations.values()
    type_counts = pd.Series(primary_types).value_counts()

    color_map = {}

    for tissue in type_counts.index:
        gradient = tissue_palette[tissue]
        colors = gradient.resample_as_palette(type_counts[tissue] + 2)

        # Get subset of dictionary with clusters that have the tissue as their top cell type.
        relevant_clusters = {k: v for k, v in cluster_annotations.items() if v == tissue}

        for cluster, color in zip(relevant_clusters.keys(), colors[1:], strict=False):
            color_map[cluster] = color

    return color_map


def process_quantified_cluster_annotations(
    annotation_filepath: str,
) -> tuple[pd.DataFrame, dict]:
    """Process cluster annotations from a TSV file and return formatted results.

    Args:
        annotation_filepath: Path to the TSV file containing cluster annotations

    Returns:
        tuple containing:
            - DataFrame with cluster statistics and formatted names
            - Dictionary mapping clusters to their top cell types
    """
    CAC = ClusterAnnotationsColumns

    cluster_annotations = pd.read_csv(annotation_filepath, sep="\t")
    stage = cluster_annotations["stage"].unique()[0]
    index = CionaStage.ordered_stages().index(stage) + 1

    clusters = []
    top_cluster_celltypes = []
    top_cluster_fractions = []
    second_cluster_celltypes = []
    second_cluster_fractions = []
    formatted_cluster_names = []

    clusters = cluster_annotations[CAC.SEURAT_CLUSTERS].unique()
    for cluster in clusters:
        cluster_data = cluster_annotations[cluster_annotations[CAC.SEURAT_CLUSTERS] == cluster]

        top_types = _get_top_cell_types(cluster_data, n=2)

        top_cluster_celltypes.append(top_types[0][0])
        top_cluster_fractions.append(top_types[0][1])

        if len(top_types) > 1:
            second_cluster_celltypes.append(top_types[1][0])
            second_cluster_fractions.append(top_types[1][1])
        else:
            second_cluster_celltypes.append(np.nan)
            second_cluster_fractions.append(np.nan)

    top_cluster_mapping = dict(zip(clusters, top_cluster_celltypes, strict=True))
    suffixes = _number_repeated_types(top_cluster_mapping)

    result = pd.DataFrame(
        {
            CAC.INDEX_STAGE: f"{index}_{stage}",
            CAC.SEURAT_CLUSTERS: clusters,
            CAC.TOP_CLUSTER_CELLTYPE: top_cluster_celltypes,
            CAC.TOP_CLUSTER_FRACTION: top_cluster_fractions,
            CAC.SECOND_CLUSTER_CELLTYPE: second_cluster_celltypes,
            CAC.SECOND_CLUSTER_FRACTION: second_cluster_fractions,
            CAC.TOP_CLUSTER_SUFFIX: [
                zt.utils.cast_numeric_id_as_string(value) for value in suffixes.values()
            ],
        }
    )

    cluster_colors = _assign_cluster_tissue_colors(top_cluster_mapping)
    result.insert(2, CAC.CLUSTER_TISSUE_COLOR, cluster_colors)

    formatted_cluster_names = result.apply(_format_cell_type_row, axis=1)
    result.insert(2, CAC.FORMATTED_CLUSTER_NAME, formatted_cluster_names)

    return result


def process_all_quantified_cluster_annotations(
    annotation_dirpath: str | Path, output_filepath: str | Path | None = None
) -> pd.DataFrame:
    """
    Compile cluster names across all developmental stages into a single dataframe,
        and optionally save the output to a TSV file.
    The output file includes the formatted cluster names and the tissue colors.

    Args:
        annotation_dir: Directory containing cluster annotation TSV files
        output_filepath: Optional path to save compiled cluster names

    Returns:
        DataFrame containing cluster names for all stages
    """
    annotation_dirpath = Path(annotation_dirpath)
    all_cluster_names = []

    for stage in CionaStage.ordered_stages():
        index = CionaStage.ordered_stages().index(stage) + 1

        cluster_annotations = process_quantified_cluster_annotations(
            annotation_dirpath / f"{index}_{stage}_cluster_annotations.tsv"
        )
        all_cluster_names.append(cluster_annotations)

    compiled_annotations = pd.concat(all_cluster_names, axis=0, ignore_index=True)

    if output_filepath:
        compiled_annotations.to_csv(output_filepath, sep="\t", index=False)

    return compiled_annotations
