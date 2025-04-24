from pathlib import Path

import pandas as pd

from zoogletools.constants import (
    DATA_DIR,
)
from zoogletools.utils import (
    create_marrvel_url,
    create_omim_url,
    create_opentargets_url,
    create_orphanet_url,
    create_zoogle_url,
)

CLINVAR_DISEASES = DATA_DIR / "2025-03-20-disease-data.tsv"


def generate_minimum_rank_comparison_data(
    organism: str,
    data_dir: str | Path = DATA_DIR,
):
    data_dir = Path(data_dir)
    mouse_data = pd.read_csv(data_dir / f"per-nonref-species/{organism}.tsv", sep="\t")
    min_mouse_ranks = mouse_data.groupby("ref_protein")["portfolio_rank"].min().reset_index()
    min_mouse_ranks.rename(columns={"portfolio_rank": f"{organism}_min_rank"}, inplace=True)

    mouse_data = mouse_data.merge(min_mouse_ranks, on="ref_protein", how="left")

    results = mouse_data[["hgnc_gene_symbol", f"{organism}_min_rank", "orthogroup_homolog_count"]]
    results = results.rename(
        columns={"orthogroup_homolog_count": f"{organism}_orthogroup_homolog_count"}
    )
    results = results.drop_duplicates().reset_index(drop=True)
    return results


def filter_organism_proteins(
    target_organism: str,
    comparison_organisms: list[str] | None = None,
    max_portfolio_rank: int | None = None,
    max_pvalue_rowwise: float | None = None,
    max_pvalue_colwise: float | None = None,
    max_trait_distance: int | None = None,
    max_homolog_count: int | None = None,
    filter_monodisease: bool = True,
    filter_comparison_organism_ranks: bool = False,
    filter_comparsion_organism_homologs: bool = False,
    data_dir: str | Path = DATA_DIR,
) -> tuple[pd.DataFrame, dict]:
    filter_counts = dict()

    data_dir = Path(data_dir)

    data = pd.read_csv(data_dir / f"per-nonref-species/{target_organism}.tsv", sep="\t")
    filter_counts["total"] = len(data)

    if max_portfolio_rank:
        data = data[data["portfolio_rank"] <= max_portfolio_rank]
        filter_counts[f"portfolio_rank <= {max_portfolio_rank}"] = len(data)

    if max_trait_distance:
        data = data[data["trait_dist"] <= max_trait_distance]
        filter_counts[f"trait_distance <= {max_trait_distance}"] = len(data)

    if max_pvalue_colwise:
        data = data[data["pvalue_colwise"] <= max_pvalue_colwise]
        filter_counts[f"pvalue_colwise <= {max_pvalue_colwise}"] = len(data)

    if max_pvalue_rowwise:
        data = data[data["pvalue_rowwise"] <= max_pvalue_rowwise]
        filter_counts[f"pvalue_rowwise <= {max_pvalue_rowwise}"] = len(data)

    if max_homolog_count:
        data = data[data["orthogroup_homolog_count"] <= max_homolog_count]
        filter_counts[f"orthogroup_homolog_count <= {max_homolog_count}"] = len(data)

    clinvar_diseases = pd.read_csv(CLINVAR_DISEASES, sep="\t")
    data = data.merge(clinvar_diseases, how="left", left_on="hgnc_gene_symbol", right_on="symbol")
    data = data[data["symbol"].notna()]

    filter_counts["associated_with_disease"] = len(data)

    if filter_monodisease:
        data = data[data["disease_count"] == 1]
        filter_counts["monodisease"] = len(data)

    if comparison_organisms:
        for organism in comparison_organisms:
            comparison_data = generate_minimum_rank_comparison_data(organism)
            data = data.merge(
                comparison_data, how="left", left_on="hgnc_gene_symbol", right_on="hgnc_gene_symbol"
            )

    if filter_comparison_organism_ranks:
        for organism in comparison_organisms:
            data = data[data[f"{organism}_min_rank"] > data["portfolio_rank"]]
        filter_counts["comparison_organism_rank"] = len(data)

    if filter_comparsion_organism_homologs:
        for organism in comparison_organisms:
            data = data[
                data[f"{organism}_orthogroup_homolog_count"] >= data["orthogroup_homolog_count"]
            ]
        filter_counts["comparison_organism_homologs"] = len(data)

    return data, filter_counts


def create_final_filtered_sheet(
    input_path: str,
    output_path: str,
):
    """
    Joins two TSV files (genes and IDs) and adds columns with URLs to external resources.

    Args:
        genes_tsv_path (str): Path to the TSV file containing genes data
        ids_tsv_path (str): Path to the TSV file containing IDs data
        output_tsv_path (str): Path where the output TSV will be saved

    Returns:
        pd.DataFrame: The joined DataFrame with added URL columns
    """
    # Columns from the HGNC IDs file.
    omim_col = "omim_id"
    orphanet_col = "orphanet"
    uniprot_col = "uniprot_id"
    ensembl_col = "ensembl_gene_id"
    entrez_col = "entrez_id"

    # Columns from the genes file.
    gene_symbol_col = "hgnc_gene_symbol"

    input_data = pd.read_csv(input_path, sep="\t")

    input_data["zoogle_url"] = input_data.apply(
        lambda row: create_zoogle_url(row[gene_symbol_col], row[uniprot_col])
        if pd.notna(row[uniprot_col])
        else None,
        axis=1,
    )
    input_data["omim_url"] = input_data.apply(
        lambda row: create_omim_url(row[gene_symbol_col], row[omim_col])
        if pd.notna(row[omim_col])
        else None,
        axis=1,
    )
    input_data["orphanet_url"] = input_data.apply(
        lambda row: create_orphanet_url(row[gene_symbol_col], row[orphanet_col])
        if pd.notna(row[orphanet_col])
        else None,
        axis=1,
    )
    input_data["opentargets_url"] = input_data.apply(
        lambda row: create_opentargets_url(row[ensembl_col])
        if pd.notna(row[ensembl_col])
        else None,
        axis=1,
    )
    input_data["marrvel_url"] = input_data.apply(
        lambda row: create_marrvel_url(row[entrez_col]) if pd.notna(row[entrez_col]) else None,
        axis=1,
    )

    # Final columns to add to file.
    final_columns = {
        "hgnc_gene_symbol": "gene_symbol",
        "disease_names": "disease_names",
        "HighestPrevalenceClass": "orphanet_highest_prevalence_class",
        "PrettifiedPrevalences": "orphanet_prevalences",
        "ref_protein": "human_protein",
        "trait_dist": "trait_distance",
        "within_gene_family_percentile": "percentile",
        "pvalue_rowwise": "pvalue_rowwise",
        "pvalue_colwise": "pvalue_colwise",
        "nonref_protein": "organism_protein",
        "portfolio_rank": "portfolio_rank",
        "zoogle_url": "zoogle",
        "omim_url": "omim",
        "orphanet_url": "orphanet",
        "opentargets_url": "opentargets",
        "marrvel_url": "marrvel",
    }

    # Add additional columns (for organism-specific ranks and homolog counts)
    # that are not in the final_columns dictionary.
    additional_columns = [
        col
        for col in input_data.columns
        if "_min_rank" in col or "_orthogroup_homolog_count" in col
    ]
    final_df = input_data[list(final_columns.keys()) + additional_columns]

    final_df.rename(columns=final_columns, inplace=True)
    final_df.to_csv(output_path, sep="\t", index=False)

    return final_df
