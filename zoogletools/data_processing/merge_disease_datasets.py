import pathlib

import click
import pandas as pd


@click.command()
@click.option(
    "--hgnc-filepath",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, path_type=pathlib.Path),
)
@click.option(
    "--orphanet-filepath",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, path_type=pathlib.Path),
)
@click.option(
    "--clinvar-filepath",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, path_type=pathlib.Path),
)
@click.option(
    "--output-filepath",
    required=True,
    type=click.Path(exists=False, file_okay=True, dir_okay=False, path_type=pathlib.Path),
)
def process(
    hgnc_filepath: pathlib.Path,
    orphanet_filepath: pathlib.Path,
    clinvar_filepath: pathlib.Path,
    output_filepath: pathlib.Path,
) -> None:
    hgnc_dataset = pd.read_csv(hgnc_filepath, sep="\t")
    orphanet_dataset = pd.read_csv(orphanet_filepath, sep="\t")
    clinvar_dataset = pd.read_csv(clinvar_filepath, sep="\t")

    hgnc_with_orphanet = hgnc_dataset.merge(
        orphanet_dataset,
        how="left",
        left_on="uniprot_id",
        right_on="uniprot_id",
    )

    if len(hgnc_with_orphanet) > len(hgnc_dataset):
        raise ValueError(
            "Adding Orphanet data resulted in more rows than the original HGNC dataset."
        )

    hgnc_with_clinvar = hgnc_with_orphanet.merge(
        clinvar_dataset,
        how="left",
        left_on="symbol",
        right_on="associated_gene",
    )

    if len(hgnc_with_clinvar) > len(hgnc_with_orphanet):
        raise ValueError(
            "Adding ClinVar data resulted in more rows than the original HGNC dataset."
        )

    final_dataset = hgnc_with_clinvar.copy().drop("associated_gene", axis=1)

    final_dataset.sort_values(by="symbol").to_csv(output_filepath, index=False, sep="\t")


if __name__ == "__main__":
    process()
