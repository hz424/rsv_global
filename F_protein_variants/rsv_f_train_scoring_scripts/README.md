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
