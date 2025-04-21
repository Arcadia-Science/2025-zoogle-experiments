"""
Snakemake workflow to download files and directories from S3 and files from URLs.

Note: this Snakefile uses dynamic rule generation to create one rule
per file or directory to download. While this approach diverges
from Snakemake's typical filepath-pattern-matching paradigm,
it's used here as a convenient way to make the workflow configuration-driven
(i.e., all output files are specified in the config file and there are no input files).
"""

import os
from datetime import datetime
import json
from pathlib import Path


configfile: "workflows/download_config.yaml"


s3_env = "envs/s3_operations.yaml"
url_env = "envs/url_downloads.yaml"

data_dir = Path(config["data_dir"]).resolve()
log_dir = Path(config["log_dir"]).resolve()

S3_DIR_OUTPUTS = [data_dir / s3_dir["local_dir"] for s3_dir in config["s3_directories"]]
S3_FILE_OUTPUTS = [data_dir / s3_file["output"] for s3_file in config["s3_files"]]
URL_FILE_OUTPUTS = [data_dir / url_file["output"] for url_file in config["url_files"]]

os.makedirs(data_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)


rule all:
    input:
        *S3_DIR_OUTPUTS,
        *S3_FILE_OUTPUTS,
        *URL_FILE_OUTPUTS,
        data_dir / "metadata.json",


# These rules are programmatically generated, following the pattern in this StackOverflow answer:
# https://stackoverflow.com/questions/77244937

# Rules for S3 directories.
for s3_dir in config["s3_directories"]:

    def create_s3_dir_rule(s3_dir=s3_dir):

        rule:
            name:
                f"s3_dir_{s3_dir['local_dir'].replace('/', '_')}"
            output:
                directory(data_dir / s3_dir["local_dir"]),
            log:
                log_dir / f"{s3_dir['local_dir'].replace('/', '_')}.log",
            params:
                bucket=s3_dir["bucket"],
                prefix=s3_dir["prefix"],
            conda:
                s3_env
            shell:
                """
                mkdir -p {output}
                aws s3 sync s3://{params.bucket}/{params.prefix} {output} --no-progress 2> {log}
                """

    create_s3_dir_rule()


# Rules for S3 single files.
for s3_file in config["s3_files"]:

    def create_s3_file_rule(s3_file=s3_file):

        rule:
            name:
                f"s3_file_{s3_file['output'].replace('/', '_')}"
            output:
                data_dir / s3_file["output"],
            log:
                log_dir / f"{s3_file['output'].replace('/', '_')}.log",
            params:
                bucket=s3_file["bucket"],
                key=s3_file["key"],
            conda:
                s3_env
            shell:
                """
                mkdir -p $(dirname {output})
                aws s3 cp s3://{params.bucket}/{params.key} {output} 2> {log}
                """

    create_s3_file_rule()


# Rules for URL files.
for url_file in config["url_files"]:

    def create_url_file_rule(url_file=url_file):

        rule:
            name:
                f"url_file_{url_file['output'].replace('/', '_')}"
            output:
                data_dir / url_file["output"],
            log:
                log_dir / f"{url_file['output'].replace('/', '_')}.log",
            params:
                url=url_file["url"],
            conda:
                url_env
            shell:
                """
                mkdir -p $(dirname {output})
                curl -L {params.url} -o {output} 2> {log}
                """

    create_url_file_rule()


# Rule for metadata. This records the source, accessed date, and size of each file.
rule metadata:
    input:
        *URL_FILE_OUTPUTS,
    output:
        data_dir / "metadata.json",
    run:
        metadata = {
            "files": {
                os.path.basename(str(input_file)): {
                    "source": next(
                        url_file["url"]
                        for url_file in config["url_files"]
                        if str(data_dir / url_file["output"]) == str(input_file)
                    ),
                    "accessed": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "description": next(
                        url_file["description"]
                        for url_file in config["url_files"]
                        if str(data_dir / url_file["output"]) == str(input_file)
                    ),
                    "size_bytes": os.path.getsize(str(input_file)),
                }
                for input_file in input
            }
        }

        with open(str(output[0]), "w") as f:
            json.dump(metadata, f, indent=2)
