import pandas as pd

ZOOGLE_URL_PREFIX = "https://zoogle.arcadiascience.com/search?gene="
OMIM_URL_PREFIX = "https://omim.org/entry/"
ORPHANET_URL_PREFIX = "https://www.orpha.net/en/disease/list/gene/"
OPENTARGETS_URL_PREFIX = "https://platform.opentargets.org/target/"


def create_zoogle_url(gene_symbol: str, uniprot_accession: str) -> str:
    """
    Creates a Zoogle URL based on gene symbol and UniProt accession.
    Example: https://zoogle.arcadiascience.com/search?gene=AAAS-Q9NRG9
    """
    if pd.isna(gene_symbol) or pd.isna(uniprot_accession):
        return None

    return f"{ZOOGLE_URL_PREFIX}{gene_symbol}-{uniprot_accession}"


def create_omim_url(gene_symbol: str, omim_id: str) -> str:
    """
    Creates an OMIM URL based on gene symbol and OMIM ID.
    Example: https://omim.org/entry/605378?search=aaas&highlight=aaas
    """
    if pd.isna(gene_symbol) or pd.isna(omim_id):
        return None

    gene_symbol_lower = gene_symbol.lower()

    return f"{OMIM_URL_PREFIX}{omim_id}?search={gene_symbol_lower}&highlight={gene_symbol_lower}"


def create_orphanet_url(gene_symbol: str, orphanet_id: float | int) -> str:
    """
    Creates an Orphanet URL based on gene symbol and Orphanet ID.
    Example: https://www.orpha.net/en/disease/list/gene/AAAS?orpha=117613
    """
    if pd.isna(gene_symbol) or pd.isna(orphanet_id):
        return None

    return f"{ORPHANET_URL_PREFIX}{gene_symbol}?orpha={int(orphanet_id)}"


def create_opentargets_url(ensembl_gene_id: str):
    """
    Creates an OpenTargets URL based on Ensembl gene ID.
    Example: https://www.targetvalidation.org/target/ENSG00000132155
    """
    if pd.isna(ensembl_gene_id):
        return None

    return f"{OPENTARGETS_URL_PREFIX}{ensembl_gene_id}"


def create_marrvel_url(entrez_id: str):
    """
    Creates a MARRVEL URL based on Entrez ID.
    Example: https://marrvel.org/variant/entrez_id
    """
    if pd.isna(entrez_id):
        return None

    return f"https://marrvel.org/human/gene/{int(entrez_id)}"


def check_many_to_many_mappings(
    df: pd.DataFrame, col1: str, col2: str, max_mappings: int = 10
) -> pd.Series | None:
    """
    Check for many-to-many mappings between two columns in a DataFrame.
    """
    # Group by col1 and count the unique values in col2 for each group
    mapping_counts = df.groupby(col1)[col2].nunique()

    # Find identifiers in col1 that map to multiple identifiers in col2
    multiple_mappings = mapping_counts[mapping_counts > 1]

    if len(multiple_mappings) > 0:
        print(
            f"Found {len(multiple_mappings)} identifiers in "
            + f"{col1} that map to multiple values in {col2}:"
        )

        print(f"Showing first {max_mappings} multi-mappings found:")

        # For each identifier with multiple mappings, show the actual mappings
        for identifier in multiple_mappings.index[:max_mappings]:
            mappings = df[df[col1] == identifier][col2].unique()
            print(f"  {identifier} maps to {len(mappings)} values: {list(mappings)}")

        return multiple_mappings
    else:
        print(f"No identifiers in {col1} map to multiple values in {col2}")
        return None


def join_non_na_values(series: pd.Series, sep: str = "; ") -> str:
    """
    Join non-NA values in a series into a string, separated by a separator.
    """
    non_na_values = [value for value in series if pd.notna(value)]
    string_values = [str(value) for value in non_na_values]
    return sep.join(string_values)


def join_unique_values(series: pd.Series, sep: str = "; ") -> str:
    """
    Join unique values in a series into a string, separated by a separator.
    """
    unique_values = series.dropna().unique()
    string_values = [str(value) for value in unique_values]
    return sep.join(string_values)


def cast_numeric_id_as_string(input: float | int | str) -> str:
    """
    Apply this function to a column to convert numeric IDs to strings,
    avoiding the default float representation.
    """
    if pd.notna(input):
        return str(int(input))
    else:
        return input
