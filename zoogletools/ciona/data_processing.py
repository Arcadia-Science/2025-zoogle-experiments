from pathlib import Path

import pandas as pd
import scanpy as sc

from zoogletools.ciona.constants import (
    CAO_STAGE_NAME_TO_PIEKARZ_STAGE_NAME_MAP,
    CIONA_STAGE_CAO_TO_PIEKARZ_MAP,
    CionaStage,
)


def load_cell_clusters(data_dir: Path = Path("../../data/SCP454")) -> pd.DataFrame:
    """
    Load and process the cell cluster data from the Cao et al. dataset.

    Parameters:
        data_dir (Path): Path to the directory containing the cell cluster data

    Returns:
        pd.DataFrame: Processed cell cluster data with additional columns
            for stage and barcode information
    """
    cell_clusters = pd.read_csv(
        data_dir / "cluster" / "ciona10stage.cluster.upload.new.txt", sep="\t"
    )
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


def _append_tech_replicate_to_barcode(row: pd.Series) -> str:
    if row["tech"] == "NA":
        return row["barcode"]
    return row["barcode"] + "-" + row["tech"].replace("tech", "")


def load_ciona_scrnaseq_data(stage: CionaStage, data_dir: str | Path) -> sc.AnnData:
    """
    Load the scRNA-seq data for a given developmental stage.

    There are a number of inconsistencies in the Piekarz data, including:
      - inconsistent column names
      - incorrect replicate numbers
      - incorrect stage names
      - different cell barcode formats compared to Cao et al. for technical replicates

    This function takes care of all of these issues, returning a clean AnnData object
    for the correct stage.

    This function duplicates the similar function in the `disambiguation.py` module,
    but includes corrections for the stage and replicate mismatches that are
    illustrated in the `disambiguation.py` module.

    Parameters:
    stage (CionaStage): The developmental stage to load the data for.
    data_dir (str | Path): The directory containing the scRNA-seq data.

    Returns:
    adata (AnnData): The scRNA-seq data for the given developmental stage.
    """
    piekarz_stage_to_retrieve = CIONA_STAGE_CAO_TO_PIEKARZ_MAP[stage.value]
    index = CionaStage.ordered_stages().index(piekarz_stage_to_retrieve)
    filepath = (
        f"{data_dir}/{index + 1}_{piekarz_stage_to_retrieve}/final/"
        + f"{piekarz_stage_to_retrieve}.h5ad"
    )

    # The midTII stage file is named differently.
    if piekarz_stage_to_retrieve == CionaStage.MIDTII:
        filepath = (
            f"{data_dir}/{index + 1}_{piekarz_stage_to_retrieve}/final/"
            + f"{piekarz_stage_to_retrieve}_KY21.h5ad"
        )

    adata = sc.read_h5ad(filepath)

    # For some reason, the replicate column is named "larva" in the iniG stage.
    if piekarz_stage_to_retrieve == CionaStage.INIG:
        adata.obs.rename(columns={"larva": "rep"}, inplace=True)

    # The earTI stage file has an incorrect replicate number.
    if piekarz_stage_to_retrieve == CionaStage.EARTI:
        adata.obs["rep"] = adata.obs["rep"].cat.rename_categories({"rep2": "rep3"})

    # The latTI stage file has swapped replicate numbers.
    if piekarz_stage_to_retrieve == CionaStage.LATTI:
        adata.obs["rep"] = adata.obs["rep"].cat.rename_categories({"rep1": "rep2", "rep2": "rep1"})

    # The larva stage file has incorrect replicate numbers.
    # The Piekarz analysis seems to have omitted the original "rep1" file.
    if piekarz_stage_to_retrieve == CionaStage.LARVA:
        adata.obs["rep"] = adata.obs["rep"].cat.rename_categories({"rep1": "rep3", "rep2": "rep4"})

    # For some reason, the replicate column is named "replicate" in the midTII stage file.
    if piekarz_stage_to_retrieve == CionaStage.MIDTII:
        adata.obs.rename(columns={"replicate": "rep"}, inplace=True)
        adata.obs["rep"] = adata.obs["rep"].cat.rename_categories(
            {"midTII-1": "rep1", "midTII-2": "rep2"}
        )

    # The latTI and latTII files include technical replicates,
    # which are labeled in an additional "tech" column.
    # To be able to merge these with the Cao cell annotations,
    # we need to append "-1" and "-2" to the barcode.
    if piekarz_stage_to_retrieve in [CionaStage.LATTI, CionaStage.LATTI]:
        adata.obs["barcode"] = adata.obs.apply(_append_tech_replicate_to_barcode, axis=1)

    adata.obs["merging_barcode"] = (
        piekarz_stage_to_retrieve
        + "_"
        + adata.obs["rep"].astype(str)
        + "_"
        + adata.obs["barcode"].astype(str)
    )

    return adata
