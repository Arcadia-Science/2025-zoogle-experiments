# Processing organism selection data to include p-values

## Data versions
In this repo, we don't use the dataset that was uploaded to Zenodo. Instead, we use the original organism selection dataset found in S3 as a starting point, which is the same dataset that had been used in the previous RAAS/ Compass repos.


## Dataset Descriptions
The current organism selection dataset uploaded to Zenodo contains the following columns:
- `gene_family`
- `nonref_species`
- `nonref_protein`
- `ref_protein`
- `phylo_dist`
- `trait_dist`
- `rank_trait_dist`

We add additional columns to the dataset when processing data for Zoogle, including:
- `within_gene_family_percentile`
- `organism_rank`
- `portfolio_rank` - this replaces `rank_trait_dist`, accounting for the 2 incorrect organisms.
- `orthogroup_homolog_rank`
- `orthogroup_homolog_count`
- `hgnc_gene_symbol`

To have more principled filtering of the dataset, we'd like to also include the row-wise and col-wise p-values that had been previously calculated for the dataset, explained in the [original raas repo](https://github.com/Arcadia-Science/raas-organism-prioritization?tab=readme-ov-file#key-outputs).

To do this, we'll need to download the RAAS dataset from S3. The S3 dataset includes the following columns, in addition to all the colunmns in the Zenodo data:
Columns we want to keep:
- `pvalue_rowwise`
- `pvalue_colwise`
Columns to drop (redundant with later processing):
- `associated_gene`
- `disease_name`
- `concept_id`
- `source_name`
- `source_id`
- `disease_mim`

## Downloading
Dataset downloading happens when running the `download_data.snakefile` workflow.

## Processing
To get the re-processed dataset with p-values, run the following commands from the top-level directory of the repo:
1. Remove unnecessary columns from the S3 dataset.
    ```{bash}
    python zoogletools/data_processing/clean_raas_dataset.py \
    --input-dirpath data/raas_final_protein_pair_summary_tables \
    --output-dirpath data/raas_cleaned
    ```

2. Use the slightly modified dataset processing scripts from the [`2025-organism-selection-portal`](https://github.com/Arcadia-Science/2025-organism-selection-portal) repo.
   1. Concatenate the cleaned dataset files.
    ```{bash}
    python zoogletools/data_processing/os_dataset.py concat \
    --dataset-dirpath data/raas_cleaned \
    --output-filepath data/YYYY-MM-DD-concatenated-os-dataset.tsv
    ```

    2. Process the HGNC dataset.
    ```{bash}
    python zoogletools/data_processing/hgnc_dataset.py \
    --input-filepath data/hgnc_complete_set.tsv \
    --output-filepath data/YYYY-MM-DD-hgnc-processed.tsv
    ```

    3. Process the ClinVar dataset.
    ```{bash}
    python zoogletools/data_processing/clinvar_dataset.py \
    --input-filepath data/gene_condition_source_id.tsv \
    --output-filepath data/YYYY-MM-DD-clinvar-processed.tsv
    ```

    4. Process the Orphanet dataset.
    ```{bash}
    python zoogletools/data_processing/orphanet_dataset.py \
    --genes-filepath data/orphanet/genes.csv \
    --prevalence-filepath data/orphanet/prevalence.csv \
    --output-filepath data/YYYY-MM-DD-orphanet-processed.tsv
    ```

    5. Merge the processed HGNC, ClinVar, and Orphanet datasets into a single TSV file.
    ```{bash}
    python zoogletools/data_processing/merge_disease_datasets.py \
    --hgnc-filepath data/YYYY-MM-DD-hgnc-processed.tsv \
    --orphanet-filepath data/YYYY-MM-DD-orphanet-processed.tsv \
    --clinvar-filepath data/YYYY-MM-DD-clinvar-processed.tsv \
    --output-filepath data/YYYY-MM-DD-merged-disease-datasets.tsv
    ```

    6. Process the concatenated OS dataset (from step 1 above) and merge it with the HGNC dataset.
    ```{bash}
    python zoogletools/data_processing/os_dataset.py process \
    --dataset-filepath data/YYYY-MM-DD-concatenated-os-dataset.tsv \
    --hgnc-dataset-filepath data/YYYY-MM-DD-hgnc-processed.tsv \
    --output-dirpath data/YYYY-MM-DD-os-portal-reprocessed
    ```

    7. Run a sanity check on the final OS dataset.
    ```{bash}
    python zoogletools/data_processing/os_dataset.py sanity-check \
    --dataset-filepath data/YYYY-MM-DD-concatenated-os-dataset.tsv \
    --output-dirpath data/YYYY-MM-DD-os-portal-reprocessed
    ```
