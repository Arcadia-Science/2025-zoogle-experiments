from enum import StrEnum
from pathlib import Path
from typing import NamedTuple

import pandas as pd

from zoogletools.ciona.constants import (
    CIONA_GENE_MODELS_DIRPATH,
    ZOOGLE_RESULTS_DIRPATH,
)


class CionaIDTypes(StrEnum):
    HGNC_GENE_SYMBOL = "hgnc_gene_symbol"
    NONREF_PROTEIN = "nonref_protein"
    KY_ID = "ky_id"
    KH_ID = "kh_id"


class MappingConfig(NamedTuple):
    """Configuration for a direct mapping between two ID types."""

    from_type: CionaIDTypes
    to_type: CionaIDTypes
    hits_dict_name: str
    requires_intermediate: bool = False
    intermediate_type: CionaIDTypes | None = None


class IdentifierMapper:
    """A class to handle mapping between different Ciona identifier types.

    This class provides a unified interface for mapping between different identifier types
    used in our analysis, including HGNC gene symbols, UniProt IDs, KY IDs, and KH IDs.
    It handles loading and caching of mapping files, and provides both direct and
    multi-step mapping capabilities.
    """

    # Define all possible direct mappings
    DIRECT_MAPPINGS = [
        MappingConfig(
            CionaIDTypes.HGNC_GENE_SYMBOL,
            CionaIDTypes.NONREF_PROTEIN,
            "hgnc_to_uniprot_hits",
        ),
        MappingConfig(
            CionaIDTypes.NONREF_PROTEIN,
            CionaIDTypes.HGNC_GENE_SYMBOL,
            "uniprot_to_hgnc_hits",
        ),
        MappingConfig(
            CionaIDTypes.NONREF_PROTEIN,
            CionaIDTypes.KY_ID,
            "uniprot_to_ky_hits",
        ),
        MappingConfig(
            CionaIDTypes.KY_ID,
            CionaIDTypes.NONREF_PROTEIN,
            "ky_to_uniprot_hits",
        ),
        MappingConfig(
            CionaIDTypes.KY_ID,
            CionaIDTypes.KH_ID,
            "ky_to_kh_hits",
        ),
        MappingConfig(
            CionaIDTypes.KH_ID,
            CionaIDTypes.KY_ID,
            "kh_to_ky_hits",
        ),
    ]

    # Define multi-step mappings
    MULTI_STEP_MAPPINGS = [
        MappingConfig(
            CionaIDTypes.HGNC_GENE_SYMBOL,
            CionaIDTypes.KY_ID,
            "hgnc_to_uniprot_hits",
            requires_intermediate=True,
            intermediate_type=CionaIDTypes.NONREF_PROTEIN,
        ),
        MappingConfig(
            CionaIDTypes.HGNC_GENE_SYMBOL,
            CionaIDTypes.KH_ID,
            "hgnc_to_uniprot_hits",
            requires_intermediate=True,
            intermediate_type=CionaIDTypes.KY_ID,
        ),
        MappingConfig(
            CionaIDTypes.NONREF_PROTEIN,
            CionaIDTypes.KH_ID,
            "uniprot_to_ky_hits",
            requires_intermediate=True,
            intermediate_type=CionaIDTypes.KY_ID,
        ),
        MappingConfig(
            CionaIDTypes.KY_ID,
            CionaIDTypes.HGNC_GENE_SYMBOL,
            "ky_to_uniprot_hits",
            requires_intermediate=True,
            intermediate_type=CionaIDTypes.NONREF_PROTEIN,
        ),
        MappingConfig(
            CionaIDTypes.KH_ID,
            CionaIDTypes.HGNC_GENE_SYMBOL,
            "kh_to_ky_hits",
            requires_intermediate=True,
            intermediate_type=CionaIDTypes.KY_ID,
        ),
        MappingConfig(
            CionaIDTypes.KH_ID,
            CionaIDTypes.NONREF_PROTEIN,
            "kh_to_ky_hits",
            requires_intermediate=True,
            intermediate_type=CionaIDTypes.KY_ID,
        ),
    ]

    def __init__(
        self,
        zoogle_results_dirpath: str | Path = ZOOGLE_RESULTS_DIRPATH,
        ciona_gene_models_dirpath: str | Path = CIONA_GENE_MODELS_DIRPATH,
    ):
        """Initialize the IdentifierMapper with mapping files.

        Args:
            zoogle_results_dirpath: Path to the Zoogle results directory
            ciona_gene_models_dirpath: Path to the Ciona gene models directory
        """
        self.zoogle_results = pd.read_csv(
            zoogle_results_dirpath / "per-nonref-species" / "Ciona-intestinalis.tsv", sep="\t"
        )
        # Sort by trait_dist for consistent first-hit mapping
        self.zoogle_results = self.zoogle_results.sort_values(by=["trait_dist"], ascending=False)

        self.uniprot_ky_map = pd.read_csv(
            ciona_gene_models_dirpath / "ciona_uniprot_ky_map.tsv", sep="\t"
        )

        self.ky_kh_map = pd.read_csv(ciona_gene_models_dirpath / "ciona_ky_kh_map.tsv", sep="\t")

        self._build_mapping_dictionaries()

    def _build_mapping_dictionaries(self):
        """Build bidirectional mapping dictionaries for faster lookups."""
        # Store all hits for each mapping
        self.hgnc_to_uniprot_hits = {}
        self.uniprot_to_hgnc_hits = {}
        self.uniprot_to_ky_hits = {}
        self.ky_to_uniprot_hits = {}
        self.ky_to_kh_hits = {}
        self.kh_to_ky_hits = {}

        # Build HGNC to UniProt and vice versa mappings
        for _, row in self.zoogle_results.iterrows():
            hgnc = row[CionaIDTypes.HGNC_GENE_SYMBOL]
            uniprot = row[CionaIDTypes.NONREF_PROTEIN]

            if hgnc not in self.hgnc_to_uniprot_hits:
                self.hgnc_to_uniprot_hits[hgnc] = []
            self.hgnc_to_uniprot_hits[hgnc].append(row)

            if uniprot not in self.uniprot_to_hgnc_hits:
                self.uniprot_to_hgnc_hits[uniprot] = []
            self.uniprot_to_hgnc_hits[uniprot].append(row)

        # Build UniProt to KY and vice versa mappings
        for _, row in self.uniprot_ky_map.iterrows():
            uniprot = row[CionaIDTypes.NONREF_PROTEIN]
            ky = row[CionaIDTypes.KY_ID]

            if uniprot not in self.uniprot_to_ky_hits:
                self.uniprot_to_ky_hits[uniprot] = []
            self.uniprot_to_ky_hits[uniprot].append(row)

            if ky not in self.ky_to_uniprot_hits:
                self.ky_to_uniprot_hits[ky] = []
            self.ky_to_uniprot_hits[ky].append(row)

        # Build KY to KH and vice versa mappings
        for _, row in self.ky_kh_map.iterrows():
            ky = row[CionaIDTypes.KY_ID]
            kh = row[CionaIDTypes.KH_ID]

            if ky not in self.ky_to_kh_hits:
                self.ky_to_kh_hits[ky] = []
            self.ky_to_kh_hits[ky].append(row)

            if kh not in self.kh_to_ky_hits:
                self.kh_to_ky_hits[kh] = []
            self.kh_to_ky_hits[kh].append(row)

    def _handle_multiple_hits(
        self,
        hits: list[pd.Series],
        input_id_type: CionaIDTypes,
        input_id: str,
        output_id_type: CionaIDTypes,
    ) -> str:
        """Handle multiple hits for a given input ID.

        Args:
            hits: List of hit rows
            input_id_type: The type of the input ID
            input_id: The input ID
            output_id_type: The type of the output ID

        Returns:
            The top hit

        Raises:
            ValueError: If no hits are found
        """
        if not hits:
            raise ValueError(f"No hits found for {input_id_type} == {input_id}")
        elif len(hits) > 1:
            print(f"Multiple hits found for {input_id_type} == {input_id}")
            hits_df = pd.DataFrame(hits)
            print(hits_df[[input_id_type, output_id_type]])
            top_hit = hits[0]
            print(f"Returning top hit: {top_hit[output_id_type]}")

        return hits[0][output_id_type]

    def _find_mapping_config(
        self, from_type: CionaIDTypes, to_type: CionaIDTypes
    ) -> MappingConfig | None:
        """Find the mapping configuration for a given from_type and to_type.

        Args:
            from_type: The type of the input identifier
            to_type: The type of the output identifier

        Returns:
            The mapping configuration if found, None otherwise
        """
        for config in self.DIRECT_MAPPINGS + self.MULTI_STEP_MAPPINGS:
            if config.from_type == from_type and config.to_type == to_type:
                return config
        return None

    def map(
        self,
        input_id: str,
        from_type: CionaIDTypes,
        to_type: CionaIDTypes,
    ) -> str:
        """Map an identifier from one type to another.

        Args:
            input_id: The input identifier
            from_type: The type of the input identifier
            to_type: The type of the output identifier

        Returns:
            The mapped identifier

        Raises:
            ValueError: If no mapping is found
        """
        if from_type == to_type:
            return input_id

        config = self._find_mapping_config(from_type, to_type)
        if config is None:
            raise ValueError(f"Unsupported mapping from {from_type} to {to_type}")

        hits_dict = getattr(self, config.hits_dict_name)
        if input_id not in hits_dict:
            raise ValueError(f"No mapping found for {from_type}: {input_id}")

        if not config.requires_intermediate:
            return self._handle_multiple_hits(
                hits_dict[input_id],
                from_type,
                input_id,
                to_type,
            )

        # Handle multi-step mappings.
        intermediate_id = self.map(input_id, from_type, config.intermediate_type)
        return self.map(intermediate_id, config.intermediate_type, to_type)

    def map_to_all(
        self,
        input_id: str,
        input_type: CionaIDTypes,
    ) -> dict[CionaIDTypes, str]:
        """Map an identifier to all other types.

        Args:
            input_id: The input identifier
            input_type: The type of the input identifier

        Returns:
            A dictionary mapping each identifier type to its corresponding value

        Raises:
            ValueError: If no mapping is found
        """
        return {
            id_type: self.map(input_id, input_type, id_type)
            for id_type in CionaIDTypes
            if id_type != input_type
        } | {input_type: input_id}
