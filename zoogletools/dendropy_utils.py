import pathlib

import dendropy
import pandas as pd


def get_species_distances(
    tree_filepath: str | pathlib.Path,
    reference_species: str = "Homo-sapiens",
):
    """
    Get phylogenetic distances between reference species and all other species in tree

    Args:
        reference_species (str): Name of reference species to calculate distances from
        tree_filepath (str): Path to Newick tree file

    Returns:
        pandas.DataFrame: DataFrame with species names and distances from reference
    """
    tree = dendropy.Tree.get(path=tree_filepath, schema="newick")

    pdc = tree.phylogenetic_distance_matrix().as_data_table()
    pdc_df = pd.DataFrame(pdc._data)
    distances = pdc_df[reference_species].reset_index()
    distances.rename(
        columns={"index": "nonref_species", reference_species: "species_dist"}, inplace=True
    )

    return distances
