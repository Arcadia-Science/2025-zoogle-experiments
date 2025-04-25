import os
import pathlib

import click
import pandas as pd
from tqdm import tqdm

# These columns are redundant with new columns that will be appended later.
DROPPED_COLUMNS = [
    "associated_gene",
    "disease_name",
    "concept_id",
    "source_name",
    "source_id",
    "disease_mim",
]


@click.command()
@click.option(
    "--input-dirpath",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Directory containing the input TSV files",
)
@click.option(
    "--output-dirpath",
    required=True,
    type=click.Path(file_okay=False, dir_okay=True),
    help="Directory where the cleaned TSV files will be saved",
)
def process(input_dirpath: pathlib.Path, output_dirpath: pathlib.Path):
    """
    Clean original RAAS dataset by removing unnecessary columns from TSV files.

    This tool processes all TSV files in the input directory and saves the cleaned
    versions to the output directory.
    """

    input_dirpath = pathlib.Path(input_dirpath)
    files = sorted(os.listdir(input_dirpath))

    output_dirpath = pathlib.Path(output_dirpath)
    output_dirpath.mkdir(parents=True, exist_ok=True)

    for file in tqdm(files, desc="Processing files", unit="file"):
        input_path = input_dirpath / file
        output_path = output_dirpath / file

        if not input_path.suffix == ".tsv":
            print(f"Skipping non-TSV file: {input_path}")
            continue

        df = pd.read_csv(input_path, sep="\t")
        df.drop(columns=DROPPED_COLUMNS, inplace=True)

        df.to_csv(output_path, sep="\t", index=False)


if __name__ == "__main__":
    process()
