import pathlib

import click
import pandas as pd

from zoogletools.utils import join_unique_values


@click.group()
def cli():
    pass


@cli.command()
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
    clinvar_dataset = pd.read_csv(input_filepath, sep="\t")

    # The `AssociatedGenes` column uses an HGNC gene symbol.
    aggregated = (
        clinvar_dataset.groupby("AssociatedGenes")
        .agg(
            DiseaseName=("DiseaseName", join_unique_values),
            DiseaseCount=("DiseaseName", "count"),
        )
        .reset_index()
    )

    aggregated.rename(
        columns={
            "AssociatedGenes": "associated_gene",
            "DiseaseName": "disease_names",
            "DiseaseCount": "disease_count",
        },
        inplace=True,
    )

    aggregated.to_csv(output_filepath, sep="\t", index=False)


if __name__ == "__main__":
    cli()
