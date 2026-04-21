# RSV Global

This repository contains source code associated with the preprint [*Genomic surveillance in the UAE reveals the global origins and local diversification of RSV lineages*](https://doi.org/10.21203/rs.3.rs-9344370/v1).

The code in this repository covers two main analysis blocks:

- `F_protein_variants/`: RSV F-protein mutation calling, UAE-versus-public comparison, and RSV F EVE/EVEscape scoring helpers
- `RSV_phylo/`: contextual phylogenetic analysis and visualization utilities for RSV-A and RSV-B

## Source Code

### `F_protein_variants/`

- `01-alignment.sh`: F-gene alignment and consensus-generation workflow
- `02-compare_by_countries_new.py`: extracts and compares RSV F mutations between UAE and public datasets
- `03-mutation_frequency_map.py`: generates mutation-frequency summaries and heatmaps
- `run_qc.sh`: helper for QC checks
- `rsv_f_train_scoring_scripts/`: minimal RSV F EVE/EVEscape scoring subset
  - `make_eve_ready_msa.py`
  - `variants_to_mutations.py`
  - `compute_wcn.py`
  - `compute_evescape_scores.py`
  - `score_variants_evescape.py`
  - `escape_utils.py`

### `RSV_phylo/`

- `augur_pipelineA.sh` and `augur_pipelineB.sh`: subtype-specific phylogeny workflows
- `reads2vcf.sh`, `bam2fasta.sh`, and `maskLowCov.py`: sequence-processing helpers
- `augur2itol.py` and `plotIntros.py`: downstream tree annotation and visualization support
- `RSV_Epi_Analysis_Mar2026 1.Rmd`: analysis notebook for epidemiologic and phylogenetic summaries

## Dependencies And Tested Environment

The included tiny EVEscape example was checked on Linux with:

- Python 3.11.5
- numpy 2.4.4
- pandas 3.0.2

Other scripts in this repository also use:

- matplotlib
- scipy
- seaborn
- biopython

Optional Python packages:

- `statsmodels`
- `adjustText`

External tools used by the sequence-processing and phylogeny workflows include:

- `mafft`
- `minimap2`
- `samtools`
- `bcftools`
- `emboss`
- `snpEff`
- Nextstrain `augur` utilities where applicable

## Install

Create a Python environment and install the Python packages:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

If you want to run the alignment or phylogeny workflows, install the command-line tools separately. One reasonable setup is:

```bash
conda install -c conda-forge -c bioconda mafft minimap2 samtools bcftools emboss snpeff
```

## Tiny Example Input

A minimal example input is provided in `toy_examples/evescape_minimal/`. It contains:

- a toy RSV-A reference CDS
- two toy EVE evolutionary-index seed files
- a three-position toy WCN table
- a three-row toy dissimilarity table

## Exact Command To Run

From the repository root:

```bash
python F_protein_variants/rsv_f_train_scoring_scripts/compute_evescape_scores.py \
  --subtype RSVA_F \
  --evol_indices_pattern "toy_examples/evescape_minimal/evol_indices_seed*.csv" \
  --wcn_csv toy_examples/evescape_minimal/wcn.csv \
  --dissimilarity_csv toy_examples/evescape_minimal/dissimilarity_metrics.csv \
  --reference_cds_fasta toy_examples/evescape_minimal/reference_cds.fasta \
  --out_csv toy_examples/evescape_minimal/example_output.csv
```

This command writes `toy_examples/evescape_minimal/example_output.csv`.

## License

This repository is distributed under the MIT License. See `LICENSE` for details.
