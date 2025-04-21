import locale
import warnings
from pathlib import Path

import click
import pandas as pd

from zoogletools.utils import join_non_na_values, join_unique_values

# Initialize locale settings
locale.setlocale(locale.LC_ALL, "")

VALMOY_DENOMINATOR = 100000

# The ranking of prevalence classes, from highest to lowest.
PREVALENCE_CLASS_RANKING = [
    ">1/1000",
    "6-9/10000",
    "1-5/10000",
    "1-9/100000",
    "1-9/1000000",
    "<1/1000000",
    "Unknown",
    "Notyetdocumented",
    pd.NA,
]


def _extract_swissprot(text: str) -> str | None:
    """
    Helper function to extract the SwissProt ID from a semicolon-separated string.

    Args:
        text: The text containing the SwissProt ID

    Returns:
        The extracted SwissProt ID or None if not found
    """
    if pd.isna(text):
        return None
    parts = text.split(";")
    for part in parts:
        if part.startswith("SwissProt:"):
            return part.split(":")[1]
    return None


def _format_number(num: float) -> str:
    """
    Helper function to format a number with a maximum of two decimal places.

    Args:
        num: The number to format

    Returns:
        The formatted number as a string
    """
    if pd.isna(num):
        return ""
    if num == 0:
        return "0"
    return locale.format_string("%.2f", num, grouping=True).rstrip("0").rstrip(".")


def _transform_row(row: pd.Series) -> str:
    """
    Transform a row of prevalence data into a prettified string.

    Args:
        row: A pandas Series containing prevalence data

    Returns:
        A formatted string representing the prevalence
    """
    case_or_family_to_plural = {"case": "cases", "family": "families"}

    result = ""

    if (
        row["PrevalenceQualification"] == "Value and class"
        and pd.notna(row["ValMoy"])
        and row["ValMoy"] > 0
    ):
        # Calculate the denominator: 100000 / val_moy
        denom = VALMOY_DENOMINATOR / row["ValMoy"]
        result = f"1/{_format_number(round(denom))}"

    elif row["PrevalenceType"] == "Cases/families" and pd.notna(row["ValMoy"]):
        suffix = (
            row["PrevalenceQualification"].lower()
            if pd.notna(row["PrevalenceQualification"])
            else ""
        )
        if suffix and row["ValMoy"] > 1 and suffix in case_or_family_to_plural:
            suffix = case_or_family_to_plural[suffix]
        result = f"{_format_number(row['ValMoy'])} {suffix}"

    else:
        result = row["PrevalenceClass"] if pd.notna(row["PrevalenceClass"]) else None

    geographic = row["PrevalenceGeographic"] if pd.notna(row["PrevalenceGeographic"]) else None

    if result is None and geographic is None:
        return "null (null)"
    elif result is None:
        return f"null ({geographic})"
    elif geographic is None:
        return f"{result} (null)"
    else:
        return f"{result} ({geographic})"


def _prettify_prevalence(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add a prettified prevalence column to the DataFrame.

    Args:
        df: DataFrame containing prevalence data

    Returns:
        DataFrame with an additional PrettifiedPrevalence column
    """
    result_df = df.copy()
    result_df["PrettifiedPrevalence"] = result_df.apply(_transform_row, axis=1)
    return result_df


def _get_highest_prevalence_class(prevalence_classes: pd.Series) -> str:
    """
    Returns the highest prevalence class from a series of values.

    Args:
        prevalence_classes: pandas Series of prevalence class values

    Returns:
        The highest ranking prevalence class value
    """
    # Filter out missing values
    valid_classes = [pc for pc in prevalence_classes if pd.notna(pc)]

    # If there are no valid classes, return NaN
    if not valid_classes:
        return pd.NA

    # Find the minimum index in the ranking list (highest rank)
    min_rank = float("inf")
    highest_class = None

    for pc in valid_classes:
        if pc in PREVALENCE_CLASS_RANKING:
            rank = PREVALENCE_CLASS_RANKING.index(pc)
            if rank < min_rank:
                min_rank = rank
                highest_class = pc

    return highest_class if highest_class is not None else valid_classes[0]


def process_orphanet_genes(
    genes_filepath: Path,
    prevalence_filepath: Path,
    output_filepath: Path | None = None,
) -> pd.DataFrame:
    """
    Process Orphanet genes and prevalence data.
    This function reads the Orphanet genes and prevalence CSV files,
    merges them, and extracts the UniProt IDs.

    Args:
        genes_filepath: Path to the Orphanet genes CSV file
        prevalence_filepath: Path to the Orphanet prevalence CSV file
        output_filepath: Optional path to save the processed output CSV file

    Returns:
        DataFrame containing the processed Orphanet data
    """
    orphanet_genes = pd.read_csv(genes_filepath)
    orphanet_prevalence = pd.read_csv(prevalence_filepath)
    orphanet_prevalence = _prettify_prevalence(orphanet_prevalence)

    # Link genes to prevalence data.
    orphanet_genes = orphanet_genes.merge(orphanet_prevalence, on="ORPHAcode", how="left")

    # Select relevant columns.
    orphanet_genes = orphanet_genes[
        ["ORPHAcode", "GeneSymbol", "ExternalReferences", "PrettifiedPrevalence", "PrevalenceClass"]
    ]

    # Extract SwissProt ID (which is a UniProt ID) from ExternalReferences.
    orphanet_genes["uniprot_id"] = orphanet_genes["ExternalReferences"].apply(_extract_swissprot)
    orphanet_genes = orphanet_genes[
        ["ORPHAcode", "GeneSymbol", "uniprot_id", "PrettifiedPrevalence", "PrevalenceClass"]
    ]

    # Aggregate prevalences based on UniProt ID.
    orphanet_prevalences = (
        orphanet_genes.groupby("uniprot_id")
        .agg(
            ORPHAcode=("ORPHAcode", join_unique_values),
            GeneSymbol=("GeneSymbol", "first"),
            PrettifiedPrevalences=("PrettifiedPrevalence", join_non_na_values),
            HighestPrevalenceClass=(
                "PrevalenceClass",
                _get_highest_prevalence_class,
            ),
        )
        .reset_index()
    )

    duplicate_genes = orphanet_prevalences["GeneSymbol"].duplicated(keep=False)

    if any(duplicate_genes):
        duplicate_gene_symbols = orphanet_prevalences.loc[duplicate_genes, "GeneSymbol"].unique()

        warning_message = (
            f"WARNING: Found {len(duplicate_gene_symbols)} gene symbols"
            + f"with multiple UniProt IDs: {', '.join(duplicate_gene_symbols)}"
        )
        warnings.warn(warning_message, stacklevel=1)

        duplicate_data = orphanet_prevalences[
            orphanet_prevalences["GeneSymbol"].isin(duplicate_gene_symbols)
        ]
        print("\nDuplicate gene data:")
        print(duplicate_data)

    if output_filepath:
        orphanet_prevalences.to_csv(output_filepath, index=False, sep="\t")
        click.echo(f"Results saved to {output_filepath}")

    return orphanet_prevalences


@click.command()
@click.option(
    "--genes-filepath",
    "-g",
    help="Path to the Orphanet genes CSV file",
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "--prevalence-filepath",
    "-p",
    help="Path to the Orphanet prevalence CSV file",
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "--output-filepath",
    "-o",
    help="Path to save the processed output CSV file",
    type=click.Path(path_type=Path),
)
@click.option(
    "--preview/--no-preview",
    default=False,
    help="Preview the first 5 rows of the processed data",
)
def process(
    genes_filepath: Path,
    prevalence_filepath: Path,
    output_filepath: Path,
    preview: bool,
) -> None:
    """Process Orphanet genes and prevalence data and extract UniProt IDs."""
    click.echo("Processing Orphanet genes and prevalence data...")
    result = process_orphanet_genes(genes_filepath, prevalence_filepath, output_filepath)

    if preview:
        click.echo("\nPreview of processed data:")
        click.echo(result.head(5))

    click.echo(f"Processed {len(result)} UniProt entries")


if __name__ == "__main__":
    process()
