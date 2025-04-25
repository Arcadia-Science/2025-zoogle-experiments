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
1. [TBI] Run filtering scripts.

### Salpingoeca rosetta RNA-Seq analysis.
1. [TBI] Download processed RNA-Seq data from Figshare.
1. [TBI] Generate scatter plots of S. rosetta gene expression.

### Ciona scRNA-Seq and RNA-Seq analysis.
1. [TBI] Download processed scRNA-Seq and RNA-Seq data from Figshare and other sources.
1. [TBI] Generate scatter plots and other plots of Ciona gene expression.

### Utilities
1. [TBI] Utility to generate a scatter plot of best hits from each organism for a given gene, *a la* the organism selection publication.
2. [TBI] `dash` app for exploring the Ciona scRNA-Seq data.

## Installation and Setup

This repository uses conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/projects/miniconda/en/latest/). After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

```{bash}
TODO: Replace <NAME> with the name of your environment
mamba env create -n <NAME> --file envs/dev.yml
conda activate <NAME>
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

To download data for this analysis, we have a Snakemake workflow. Replace `<n>` with the number of cores you want to use.

```{bash}
snakemake --snakefile workflows/download_data.snakefile --use-conda --cores <n>
```

This will download the data from the URLs specified in the `workflows/download_config.yaml` file and save it to the `data` directory.


## Overview

### Description of the folder structure

### Methods

TODO: Include a brief, step-wise overview of analyses performed.

> Example:
>
> 1.  Download scripts using `download.ipynb`.
> 2.  Preprocess using `./preprocessing.sh -a data/`
> 3.  Run Snakemake pipeline `snakemake --snakefile Snakefile`
> 4.  Generate figures using `pub/make_figures.ipynb`.

### Compute Specifications

TODO: Describe what compute resources were used to run the analysis. For example, you could list the operating system, number of cores, RAM, and storage space.

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
