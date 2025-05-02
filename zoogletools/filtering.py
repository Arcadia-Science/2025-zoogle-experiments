from abc import ABC, abstractmethod
from enum import Enum
from pathlib import Path

import pandas as pd

from zoogletools.utils import (
    create_marrvel_url,
    create_omim_url,
    create_opentargets_url,
    create_orphanet_url,
    create_zoogle_url,
)

# Column name constants
PORTFOLIO_RANK = "portfolio_rank"
TRAIT_DISTANCE = "trait_distance"
PVALUE_ROWWISE = "pvalue_rowwise"
PVALUE_COLWISE = "pvalue_colwise"
HOMOLOG_COUNT = "orthogroup_homolog_count"
DISEASE_ASSOCIATION = "associated_with_disease"
MONODISEASE = "monodisease"


def generate_minimum_rank_comparison_data(
    organism: str,
    data_dir: str | Path,
):
    data_dir = Path(data_dir)
    organism_data = pd.read_csv(data_dir / f"per-nonref-species/{organism}.tsv", sep="\t")
    min_organism_ranks = organism_data.groupby("ref_protein")["portfolio_rank"].min().reset_index()
    min_organism_ranks.rename(columns={"portfolio_rank": f"{organism}_min_rank"}, inplace=True)

    organism_data = organism_data.merge(min_organism_ranks, on="ref_protein", how="left")

    results = organism_data[
        ["hgnc_gene_symbol", f"{organism}_min_rank", "orthogroup_homolog_count"]
    ]
    results = results.rename(
        columns={"orthogroup_homolog_count": f"{organism}_orthogroup_homolog_count"}
    )
    results = results.drop_duplicates().reset_index(drop=True)
    return results


class Filter(ABC):
    """Base class for data filters."""

    name: str

    @abstractmethod
    def apply(self, data: pd.DataFrame) -> pd.DataFrame:
        """Apply the filter to the data.

        Args:
            data: Input DataFrame to filter

        Returns:
            Filtered DataFrame
        """
        pass


class FilteringPipeline:
    """Pipeline for applying a sequence of filters to data."""

    def __init__(self):
        self.filters = []
        self._last_run_stats: list[dict] | None = None

    def add_filter(self, filter: Filter):
        """Add a filter to the pipeline."""
        self.filters.append(filter)

    def apply(self, data: pd.DataFrame) -> pd.DataFrame:
        """Apply all filters in the pipeline to the data.

        Args:
            data: Input DataFrame to filter

        Returns:
            Filtered DataFrame
        """
        current_data = data
        stats = []

        for filter in self.filters:
            input_count = len(current_data)
            current_data = filter.apply(current_data)
            stats.append(
                {
                    "filter_name": filter.name,
                    "input_count": input_count,
                    "output_count": len(current_data),
                }
            )

        self._last_run_stats = stats
        return current_data

    def get_last_stats(self) -> dict:
        """Get statistics from the most recent pipeline run.

        Returns:
            dict: Statistics including total rows and filter counts
        """
        if self._last_run_stats is None:
            return {}

        stats = {"total_count": self._last_run_stats[0]["input_count"]}
        for filter_stats in self._last_run_stats:
            stats[f"{filter_stats['filter_name']}_count"] = filter_stats["output_count"]
        return stats


class PortfolioRankFilter(Filter):
    def __init__(self, max_rank: int):
        super().__init__()
        self.max_rank = max_rank
        self.name = PORTFOLIO_RANK

    def apply(self, data: pd.DataFrame) -> pd.DataFrame:
        filtered_data = data[data[PORTFOLIO_RANK] <= self.max_rank]
        return filtered_data


class TraitDistanceFilter(Filter):
    def __init__(self, max_distance: int):
        super().__init__()
        self.max_distance = max_distance
        self.name = TRAIT_DISTANCE

    def apply(self, data: pd.DataFrame) -> pd.DataFrame:
        filtered_data = data[data[TRAIT_DISTANCE] <= self.max_distance]
        return filtered_data


class PValueType(Enum):
    ROWWISE = "pvalue_rowwise"
    COLWISE = "pvalue_colwise"


class PValueFilter(Filter):
    def __init__(self, max_pvalue: float, pvalue_type: PValueType):
        super().__init__()
        if not isinstance(pvalue_type, PValueType):
            raise TypeError("pvalue_type must be a PValueType enum")
        self.max_pvalue = max_pvalue
        self.pvalue_type = pvalue_type
        self.name = pvalue_type.value

    def apply(self, data: pd.DataFrame) -> pd.DataFrame:
        filtered_data = data[data[self.pvalue_type.value] <= self.max_pvalue]
        return filtered_data


class HomologCountFilter(Filter):
    def __init__(self, max_count: int):
        super().__init__()
        self.max_count = max_count
        self.name = HOMOLOG_COUNT

    def apply(self, data: pd.DataFrame) -> pd.DataFrame:
        filtered_data = data[data[HOMOLOG_COUNT] <= self.max_count]
        return filtered_data


class DiseaseAssociationFilter(Filter):
    def __init__(self):
        super().__init__()
        self.name = DISEASE_ASSOCIATION

    def apply(self, data: pd.DataFrame) -> pd.DataFrame:
        # "symbol" is a column name only found in the processed disease data.
        # Rows missing a value in this column are not associated with a disease.
        filtered_data = data[data["symbol"].notna()]
        return filtered_data


class MonoDiseaseFilter(Filter):
    def __init__(self):
        super().__init__()
        self.name = MONODISEASE

    def apply(self, data: pd.DataFrame) -> pd.DataFrame:
        filtered_data = data[data["disease_count"] == 1]
        return filtered_data


default_filtering_pipeline = FilteringPipeline()
default_filtering_pipeline.add_filter(PValueFilter(max_pvalue=0.05, pvalue_type=PValueType.COLWISE))
default_filtering_pipeline.add_filter(DiseaseAssociationFilter())
default_filtering_pipeline.add_filter(HomologCountFilter(max_count=1))
default_filtering_pipeline.add_filter(MonoDiseaseFilter())


def filter_organism_proteins(
    target_organism: str,
    data_dirpath: str | Path,
    disease_data_filepath: str | Path,
    filtering_pipeline: FilteringPipeline = default_filtering_pipeline,
) -> tuple[pd.DataFrame, dict]:
    filter_counts = dict()

    data_dirpath = Path(data_dirpath)
    data_filepath = data_dirpath / f"per-nonref-species/{target_organism}.tsv"

    data = pd.read_csv(data_filepath, sep="\t")
    disease_data = pd.read_csv(disease_data_filepath, sep="\t")

    data = data.merge(disease_data, left_on="ref_protein", right_on="uniprot_id", how="left")

    data = filtering_pipeline.apply(data)
    filter_counts = filtering_pipeline.get_last_stats()

    return data, filter_counts


def _create_url_column(
    df: pd.DataFrame,
    url_creator: callable,
    *column_names: str,
    column_name: str,
) -> None:
    """Helper function to create URL columns in a DataFrame.

    Args:
        df: DataFrame to add URL column to
        url_creator: Function that creates the URL
        *column_names: Column names to pass to url_creator
        column_name: Name of the new URL column
    """
    df[column_name] = df.apply(
        lambda row: url_creator(*[row[col] for col in column_names])
        if all(pd.notna(row[col]) for col in column_names)
        else None,
        axis=1,
    )


def create_final_filtered_sheet(
    input_data: pd.DataFrame,
    output_path: str,
) -> pd.DataFrame:
    """
    Creates a final filtered DataFrame with added URL columns and saves it to a TSV file.

    Args:
        input_data: Input DataFrame containing gene/protein data
        output_path: Path where the output TSV will be saved

    Returns:
        pd.DataFrame: The processed DataFrame with added URL columns and renamed columns
    """
    omim_col = "omim_id"
    orphanet_col = "orphanet"
    uniprot_col = "uniprot_id"
    ensembl_col = "ensembl_gene_id"
    entrez_col = "entrez_id"
    gene_symbol_col = "hgnc_gene_symbol"

    result_data = input_data.copy()

    _create_url_column(
        result_data, create_zoogle_url, gene_symbol_col, uniprot_col, column_name="zoogle_url"
    )
    _create_url_column(
        result_data, create_omim_url, gene_symbol_col, omim_col, column_name="omim_url"
    )
    _create_url_column(
        result_data, create_orphanet_url, gene_symbol_col, orphanet_col, column_name="orphanet_url"
    )
    _create_url_column(
        result_data, create_opentargets_url, ensembl_col, column_name="opentargets_url"
    )
    _create_url_column(result_data, create_marrvel_url, entrez_col, column_name="marrvel_url")

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

    final_df = result_data[list(final_columns.keys())].copy()

    final_df.rename(columns=final_columns, inplace=True)
    final_df.sort_values(by="percentile", ascending=True, inplace=True)
    final_df.reset_index(drop=True, inplace=True)
    final_df.to_csv(output_path, sep="\t", index=False)

    return final_df
