# 2025-zoogle-collabs

[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)

## Purpose

This repository contains the code for the 2025 Zoogle collabs. This includes code and utilities for the following analyses.
[TBI] indicates that the code has been implemented in another branch but needs polishing before it's ready for review.

### Organism selection data download and processing
1. Download organism selection data and human disease data from S3 and URLs. This includes data from the following sources:
    - Original RAAS dataset from back in 2024 [here](https://github.com/Arcadia-Science/raas-organism-prioritization). This is different from the published Zenodo dataset because it includes bootstrapped row- and column-wise p-values, which were used to generate the Discovery score in Compass.
    - Processed Orphanet data from Compass.
    - ClinVar gene-condition-source-id mapping.
    - HGNC gene set.
    - Species tree from speciesRAX.
1. Process downloaded data so that it's tidy and usable for filtering scripts.
1. Run filtering scripts.

### Salpingoeca rosetta RNA-Seq analysis.
1. Download processed RNA-Seq data from Figshare.
2. Generate bar plots and heatmaps of *Salpingoeca rosetta* gene expression.

### Ciona scRNA-Seq and RNA-Seq analysis.
1. [TBI] Download processed scRNA-Seq and RNA-Seq data from Figshare and other sources.
1. [TBI] Generate scatter plots and other plots of Ciona gene expression.

### Utilities
1. Utility to generate a scatter plot of best hits from each organism for a given gene, *a la* the organism selection publication.
2. [TBI] `dash` app for exploring the Ciona scRNA-Seq data.

## Installation and Setup

This repository uses conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/projects/miniconda/en/latest/). After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

```{bash}
TODO: Replace <NAME> with the name of your environment
mamba env create -n <NAME> --file envs/dev.yml
conda activate <NAME>
```

### `zoogletools`

This repository contains a Python package called `zoogletools`. This package contains the source code for the filtering pipeline and other utilities.

To install the package, run the following command from the top level of this repository.

```{bash}
pip install -e .
```

<details><summary>Developer Notes (click to expand/collapse)</summary>

1. Install your pre-commit hooks:

    ```{bash}
    pre-commit install
    ```

    This installs the pre-commit hooks defined in your config (`./.pre-commit-config.yaml`).

2. Export your conda environment before sharing:

    As your project develops, the number of dependencies in your environment may increase. Whenever you install new dependencies (using either `pip install` or `mamba install`), you should update the environment file using the following command.

    ```{bash}
    conda env export --no-builds > envs/dev.yml
    ```

    `--no-builds` removes build specification from the exported packages to increase portability between different platforms.
</details>

## Data

To download data for this analysis, we have developed a lightweight Snakemake workflow. Replace `<n>` with the number of cores you want to use.

```{bash}
snakemake --snakefile workflows/download_data.snakefile --use-conda --cores <n>
```

This will download the data from the URLs specified in the `workflows/download_config.yaml` file and save it to the `data` directory. This can take some time (30 minutes to 1 hour) to run, depending on the number of cores you use.

### Ciona scRNA-Seq data

Accessing the Ciona scRNA-Seq data requires a few manual steps.

1. Unzip the following files using your preferred method, e.g. the `unzip` command. Ensure that the files are unzipped into the same directory as the zipped files, with the same file names.
   - `data/Ciona_intestinalis_scRNAseq_data_Piekarz.zip`
   - `data/Cirobu_RNASeq.tar.gz`
   - `data/Ciona_gene_models/HT.KY21Gene.protein.2.fasta.zip`
   - `data/Ciona_gene_models/Ciona_intestinalis.faa.gz`

    Example commands are shown below.
    ```{bash}
    unzip data/Cirobu_RNASeq.tar.gz -d data/Cirobu_RNASeq
    unzip data/HT.KY21Gene.protein.2.fasta.zip -d data/HT.KY21Gene.protein.2.fasta
    unzip data/Ciona_intestinalis.faa.gz -d data/Ciona_intestinalis
    unzip data/Ciona_intestinalis_scRNAseq_data_Piekarz.zip -d data/Ciona_intestinalis_scRNAseq_data_Piekarz
    ```

2. Download the cell type annotation for the Ciona scRNA-Seq data. You need to create an account with the Broad Institute Single Cell Portal and download the data manually. The webpage for the dataset is [here](https://singlecell.broadinstitute.org/single_cell/study/SCP454/comprehensive-single-cell-transcriptome-lineages-of-a-proto-vertebrate). Navigate to the Download section and click the "Bulk download" button. A command will be generated with an authentication code, which you can run to download the data. An example command is shown below (`<auth_code>` should be replaced with the actual authentication code).

    ```{bash}
    cd data

    curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP454&auth_code=<auth_code>&directory=all&context=study"  -o cfg.txt; curl -K cfg.txt && rm cfg.txt
    ```

    This will put the data into a directory called `SCP454`, containing three directories: `cluster`, `expression`, and `metadata`. You may encounter a `curl: (2) no URL specified` error. If this happens, you can try the command again with the authentication code. The only essential file for downstream analyses is `SCP454/cluster/ciona10stage.cluster.upload.new.txt`, which contains the cell type annotation.

## Overview

### Description of the folder structure

This repository contains the following folders:
- `data`: Contains the data for the analysis.
- `notebooks`: Contains notebooks that frame the analyses.
- `workflows`: Contains the Snakemake workflow for downloading the data.
- `zoogletools`: A pip-installable Python package containing the source code for the filtering pipeline and other utilities.


### Methods

1. **Pre-processing.** After downloading the data, the data is preprocessed and filtered using scripts included in the `zoogletools` package. See [this README](zoogletools/data_processing/README.md) for more details.

2. **Filtering.** To run the filtering pipeline, refer to the [filtering notebook](notebooks/1_filtering.ipynb).

### Compute Specifications

The analysis was performed on a MacBook Pro with an M1 Pro processor and 32GB of RAM running macOS Sequoia 15.1.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md).

---
## For Developers

This section contains information for developers who are working off of this template. Please adjust or edit this section as appropriate when you're ready to share your repo.

### GitHub templates
This template uses GitHub templates to provide checklists when making new pull requests. These templates are stored in the [.github/](./.github/) directory.

### VSCode
This template includes recommendations to VSCode users for extensions, particularly the `ruff` linter. These recommendations are stored in `.vscode/extensions.json`. When you open the repository in VSCode, you should see a prompt to install the recommended extensions.

### `.gitignore`
This template uses a `.gitignore` file to prevent certain files from being committed to the repository.

### `pyproject.toml`
`pyproject.toml` is a configuration file to specify your project's metadata and to set the behavior of other tools such as linters, type checkers etc. You can learn more [here](https://packaging.python.org/en/latest/guides/writing-pyproject-toml/)

### Linting
This template automates linting and formatting using GitHub Actions and the `ruff` linter. When you push changes to your repository, GitHub will automatically run the linter and report any errors, blocking merges until they are resolved.
