from pathlib import Path

DATA_DIR = Path("../data")
RESULTS_DIR = Path("../results")
CLINVAR_DISEASES = DATA_DIR / "gene_condition_source_id.tsv"
HGNC_IDS = DATA_DIR / "hgnc_complete_set.tsv"
ORPHANET_GENES = DATA_DIR / "orphanet/genes.csv"
ORPHANET_PREVALENCE = DATA_DIR / "orphanet/prevalence.csv"
ORPHANET_PREVALENCE_PROCESSED = DATA_DIR / "orphanet_prevalence_processed.csv"
