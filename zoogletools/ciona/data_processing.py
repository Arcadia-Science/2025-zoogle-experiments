from pathlib import Path

import pandas as pd

from .constants import CAO_STAGE_NAME_TO_PIEKARZ_STAGE_NAME_MAP


def load_cell_clusters(data_dir: Path = Path("../../data/SCP454")) -> pd.DataFrame:
    """
    Load and process the cell cluster data from the Cao et al. dataset.

    Parameters:
        data_dir (Path): Path to the directory containing the cell cluster data

    Returns:
        pd.DataFrame: Processed cell cluster data with additional columns
            for stage and barcode information
    """
    cell_clusters = pd.read_csv(data_dir / "cluster/ciona10stage.cluster.upload.new.txt", sep="\t")
    # Drop the first row, which contains data type information.
    cell_clusters.drop(index=[0], inplace=True)

    # Add columns for stage and replicate information.
    # In the cell type data, the NAME column is of the form "<stage>.<replicate_number>_<barcode>".
    cell_clusters["Cao_stage"] = cell_clusters["NAME"].apply(
        lambda x: x.split("_")[0].split(".")[0]
    )
    cell_clusters["stage"] = cell_clusters["Cao_stage"].map(
        CAO_STAGE_NAME_TO_PIEKARZ_STAGE_NAME_MAP
    )
    cell_clusters["replicate"] = cell_clusters["NAME"].apply(
        lambda x: "rep" + x.split("_")[0].split(".")[1]
    )
    cell_clusters["stage_replicate"] = cell_clusters["stage"] + "_" + cell_clusters["replicate"]
    cell_clusters["barcode"] = cell_clusters["NAME"].apply(lambda x: x.split("_")[1])

    # Add a column for the merging barcode, which is the stage_replicate and barcode concatenated.
    # A column with identical values will be added to the scRNA-seq read count matrix.
    cell_clusters["merging_barcode"] = (
        cell_clusters["stage_replicate"] + "_" + cell_clusters["barcode"]
    )
    return cell_clusters


def _append_tech_replicate_to_barcode(row: pd.Series):
    if row["tech"] == "NA":
        return row["barcode"]
    return row["barcode"] + "-" + row["tech"].replace("tech", "")
