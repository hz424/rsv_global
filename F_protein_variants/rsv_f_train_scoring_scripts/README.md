# Minimal Train + Scoring Scripts

This folder contains the minimal local RSV F scripts needed for manuscript submission.

A tiny runnable example for `compute_evescape_scores.py` is included in `toy_examples/evescape_minimal/`.

Included scripts:

- `make_eve_ready_msa.py`
  - prepare the focus-sequence-first MSA used before EVE training
- `variants_to_mutations.py`
  - convert protein FASTA variants to mutation lists for downstream scoring
- `compute_wcn.py`
  - compute structure-based WCN accessibility from an RSV F PDB; not restricted to 5UDD
- `compute_evescape_scores.py`
  - combine EVE evolutionary indices, WCN accessibility, and amino-acid dissimilarity into per-mutation EVEscape scores
- `score_variants_evescape.py`
  - aggregate per-mutation EVEscape scores to variant-level scores
- `escape_utils.py`
  - shared helpers used by `compute_evescape_scores.py` and `score_variants_evescape.py`

External EVE entrypoints are not duplicated here:

- `EVE/train_VAE.py`
- `EVE/compute_evol_indices.py`

Excluded from this minimal folder:

- `00_make_single_mut_table.py`
  - not needed when `EVE/compute_evol_indices.py` is run with `--computation_mode all_singles`
- `07_map_reference_to_pdb_5udd.py`
  - useful for alignment/QC, but not required to run the scoring pipeline
