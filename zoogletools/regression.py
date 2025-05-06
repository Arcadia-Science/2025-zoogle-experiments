import os
import pathlib

import pandas as pd
import scipy.stats
import tqdm

from zoogletools.dendropy_utils import get_species_distances


def calculate_protein_regressions(
    proteins_dirpath: str | pathlib.Path,
    tree_filepath: str | pathlib.Path,
    complexities_filepath: str | pathlib.Path,
    output_filepath: str | pathlib.Path,
) -> pd.DataFrame:
    """
    Calculate regression statistics for protein trait distances vs species distances.

    Args:
        proteins_dirpath: Path to directory containing per-protein data files
        tree_filepath: Path to phylogenetic tree file
        complexities_filepath: Path to species complexities file

    Returns:
        pandas.DataFrame: DataFrame containing regression statistics for each protein
    """
    collector = pd.DataFrame()

    distances = get_species_distances(
        tree_filepath=tree_filepath,
        complexities_filepath=complexities_filepath,
    )

    for filepath in tqdm(sorted(os.listdir(proteins_dirpath))):
        data = pd.read_csv(pathlib.Path(proteins_dirpath) / filepath, sep="\t")
        data = data.drop_duplicates(subset=["nonref_species"], keep="first")
        data = data.merge(distances, on="nonref_species", how="left")

        protein_id = filepath.split(".")[0]
        protein_symbol = data["hgnc_gene_symbol"].iloc[0]

        num_species = len(data)

        x = data["trait_dist"]
        y = data["species_dist"]

        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)

        collector = pd.concat(
            [
                collector,
                pd.DataFrame(
                    {
                        "protein_id": [protein_id],
                        "protein_symbol": [protein_symbol],
                        "slope": [slope],
                        "intercept": [intercept],
                        "r_value": [r_value],
                        "p_value": [p_value],
                        "std_err": [std_err],
                        "num_species": [num_species],
                    }
                ),
            ]
        )

        if output_filepath:
            collector.to_csv(output_filepath, index=False)

    return collector
