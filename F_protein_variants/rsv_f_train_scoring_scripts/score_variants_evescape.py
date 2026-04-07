import argparse
import math
import os
import pandas as pd

from escape_utils import (
    VARIANT_METRIC_COLUMNS,
    build_score_lookup,
    parse_mutations_csv,
    summarize_variant_from_mutations,
)


def sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def main():
    ap = argparse.ArgumentParser(
        description="Score variant mutation lists by summing per-mutation EVEscape scores."
    )
    ap.add_argument("--single_mut_scores_csv", required=True,
                    help="CSV with columns wt,i,mut,evescape (like EVEscape outputs)")
    ap.add_argument("--variants_mutations_csv", required=True,
                    help="CSV produced by 01_variants_to_mutations.py")
    ap.add_argument("--out_csv", required=True)
    args = ap.parse_args()

    out_dir = os.path.dirname(os.path.abspath(args.out_csv))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    smm = pd.read_csv(args.single_mut_scores_csv)
    required_cols = {"wt", "i", "mut", "evescape"}
    missing = required_cols - set(smm.columns)
    if missing:
        raise ValueError(f"single_mut_scores_csv missing columns: {sorted(missing)}")

    smm, score_lookup = build_score_lookup(smm)
    smm["evescape_pos"] = smm["evescape"] - smm["evescape"].min()
    smm["evescape_sigmoid"] = smm["evescape"].apply(sigmoid)

    score_map_pos = dict(zip(smm["mutation"], smm["evescape_pos"]))
    score_map_sig = dict(zip(smm["mutation"], smm["evescape_sigmoid"]))

    vdf = pd.read_csv(args.variants_mutations_csv)
    if "mutations_csv" not in vdf.columns:
        raise ValueError("variants_mutations_csv must contain a 'mutations_csv' column")

    def score_row(muts_csv: str):
        muts = parse_mutations_csv(muts_csv)
        pos = sum(score_map_pos.get(m, 0.0) for m in muts)
        sig = sum(score_map_sig.get(m, 0.0) for m in muts)
        n_missing = sum(1 for m in muts if m not in score_map_pos)
        summary = summarize_variant_from_mutations(muts, score_lookup)
        return {
            "variant_score_pos": pos,
            "variant_score_sigmoid": sig,
            "n_unscored_muts": n_missing,
            **summary,
        }

    scores = vdf["mutations_csv"].apply(score_row)
    score_df = pd.DataFrame(scores.tolist(), index=vdf.index)
    vdf[score_df.columns] = score_df

    # Lightweight QC fields for publication-quality reporting.
    if "n_substitutions" in vdf.columns:
        vdf["n_scored_muts"] = vdf["n_substitutions"].astype(int) - vdf["n_unscored_muts"].astype(int)
        denom = vdf["n_substitutions"].astype(int).clip(lower=1)
        vdf["frac_unscored_muts"] = vdf["n_unscored_muts"].astype(int) / denom

    numeric_cols = [
        "variant_score_pos",
        "variant_score_sigmoid",
        "frac_unscored_muts",
        *VARIANT_METRIC_COLUMNS,
    ]
    for col in numeric_cols:
        if col in vdf.columns:
            vdf[col] = vdf[col].astype(float)

    vdf.to_csv(args.out_csv, index=False)


if __name__ == "__main__":
    main()
