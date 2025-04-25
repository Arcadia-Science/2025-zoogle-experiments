"""
This script is based on the `hgnc_dataset.py` script in the `2025-organism-selection-portal` repo.
Some additional comments are used to explain the changes made to the original script.
"""

import pathlib

import click
import pandas as pd

HGNC_REQUIRED_COLUMNS = [
    "hgnc_id",
    "symbol",
    "uniprot_ids",
    "name",
    "entrez_id",  # Additional ID used for linking ClinVar data
    "ensembl_gene_id",  # Additional ID used for generating OpenTargets links
    "omim_id",  # Additional ID used for generating OMIM links
    "orphanet",  # Additional ID used for generating Orphanet links
]


@click.command()
@click.option(
    "--input-filepath",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, path_type=pathlib.Path),
)
@click.option(
    "--output-filepath",
    required=True,
    type=click.Path(exists=False, file_okay=True, dir_okay=False, path_type=pathlib.Path),
)
def process(
    input_filepath: pathlib.Path,
    output_filepath: pathlib.Path,
) -> None:
    """
    Process the HGNC dataset into the form we need for the OS portal.

    Outputs a TSV file with one row for each unique combination of UniProt ID, primary gene symbol,
    and alias gene symbol in the HGNC dataset.

    The output file has the following columns:
    - hgnc_id: The HGNC ID.
    - uniprot_id: The UniProt ID associated with the HGNC entry.
    - symbol: The gene symbol (either primary or alias).
    - name: The gene name.
    - entrez_id: The Entrez ID associated with the HGNC entry.
    - ensembl_gene_id: The Ensembl gene ID associated with the HGNC entry.
    - omim_id: The OMIM ID associated with the HGNC entry.
    - orphanet: The Orphanet ID associated with the HGNC entry.
    """

    # The separator used in the HGNC dataset for list-like fields.
    sep = "|"

    hgnc_dataset = pd.read_csv(input_filepath, sep="\t")

    # Filter out locus types that we know are not relevant to the OS dataset.
    excluded_locus_types = ["RNA, long non-coding", "RNA, ribosomal"]
    hgnc_dataset["locus_type"] = hgnc_dataset["locus_type"].str.strip().str.lower()
    hgnc_dataset = hgnc_dataset[
        ~hgnc_dataset.locus_type.isin([s.strip().lower() for s in excluded_locus_types])
    ].copy()

    # Drop all columns except the ones we need.
    hgnc_dataset = hgnc_dataset[HGNC_REQUIRED_COLUMNS].copy()

    # A few gene symbols have multiple UniProt IDs, so we expand them into separate rows.
    hgnc_dataset["uniprot_id"] = hgnc_dataset["uniprot_ids"].str.split(sep)
    hgnc_dataset = hgnc_dataset.explode("uniprot_id").drop(columns=["uniprot_ids"])

    # Removed the alias symbol handling, since they're not needed for this application.

    hgnc_dataset.sort_values(by="symbol").to_csv(output_filepath, index=False, sep="\t")


if __name__ == "__main__":
    process()
