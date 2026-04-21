#!/usr/bin/env python3
"""Compute EVEscape scores for RSV F protein.

Combines three components (EVEscape):
  1. Fitness        — EVE evolutionary indices (averaged across seeds)
  2. Accessibility  — WCN from PDB (negated, so higher = more accessible)
  3. Dissimilarity  — charge + Eisenberg-Weiss hydrophobicity difference

EVEscape = log σ(z_fitness/T1) + log σ(z_accessibility/T2) + log σ(z_dissimilarity/T3)
  where z = standardisation, σ = logistic, T1=1, T2=1, T3=2.

Uses only numpy + pandas (no sklearn).
"""
from __future__ import annotations
import argparse, glob, re, sys
from pathlib import Path
import numpy as np
import pandas as pd

from escape_utils import (
    accessible_amino_acids_by_distance,
    amino_acid_min_nt_distance,
    read_selected_reference_cds,
)

TEMPERATURES = {"fitness": 1.0, "surfacc": 1.0, "exchangability": 2.0}


def default_reference_cds_fasta() -> Path:
    for parent in Path(__file__).resolve().parents:
        candidate = parent / "inputs" / "references_RSV_AB" / "F_prot.ref.fasta"
        if candidate.exists():
            return candidate
    return Path(__file__).resolve().parents[2] / "inputs" / "references_RSV_AB" / "F_prot.ref.fasta"

def standardize(x):
    std = x.std()
    if std == 0 or not np.isfinite(std):
        return pd.Series(np.zeros(len(x)), index=x.index)
    return (x - x.mean()) / std

def logistic(x):
    return 1.0 / (1.0 + np.exp(-x))


def annotate_codon_distance(df: pd.DataFrame, subtype: str, reference_cds_fasta: Path) -> pd.DataFrame:
    cds = read_selected_reference_cds(reference_cds_fasta, subtype)
    if len(cds) % 3 != 0:
        raise ValueError(f"Reference CDS for {subtype} is not divisible by 3 (len={len(cds)})")

    site_meta = {}
    for pos in sorted(df["i"].astype(int).unique()):
        start = (pos - 1) * 3
        ref_codon = cds[start : start + 3]
        if len(ref_codon) != 3:
            raise ValueError(f"Reference CDS for {subtype} is too short for protein position {pos}")
        aa_min_dist = amino_acid_min_nt_distance(ref_codon)
        aa_by_dist = accessible_amino_acids_by_distance(ref_codon)
        one_nt_aas = aa_by_dist.get(1, [])
        site_meta[pos] = {
            "reference_codon": ref_codon,
            "aa_min_distance": aa_min_dist,
            "one_nt_accessible_aas": ",".join(one_nt_aas) if one_nt_aas else "",
            "n_one_nt_accessible_aas": len(one_nt_aas),
        }

    out = df.copy()
    out["reference_codon"] = out["i"].astype(int).map(lambda pos: site_meta[pos]["reference_codon"])
    out["one_nt_accessible_aas"] = out["i"].astype(int).map(lambda pos: site_meta[pos]["one_nt_accessible_aas"])
    out["n_one_nt_accessible_aas"] = out["i"].astype(int).map(lambda pos: site_meta[pos]["n_one_nt_accessible_aas"])
    out["min_nt_distance_from_reference_codon"] = out.apply(
        lambda row: site_meta[int(row["i"])]["aa_min_distance"].get(str(row["mut"]), np.nan),
        axis=1,
    )
    out["one_nt_accessible_from_reference"] = out["min_nt_distance_from_reference_codon"] == 1
    return out

def load_evol_indices(pattern, subtype):
    """Load and average evolutionary indices across multiple seed files."""
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matching {pattern}")
    dfs = []
    for f in files:
        df = pd.read_csv(f)
        df = df[df.mutations != "wt"].copy()
        dfs.append(df[["mutations", "evol_indices"]])
    merged = dfs[0].rename(columns={"evol_indices": "ei_0"})
    for k, d in enumerate(dfs[1:], 1):
        merged = merged.merge(d.rename(columns={"evol_indices": f"ei_{k}"}), on="mutations")
    ei_cols = [c for c in merged.columns if c.startswith("ei_")]
    merged["evol_indices"] = merged[ei_cols].mean(axis=1)
    merged["evol_indices_std"] = merged[ei_cols].std(axis=1)
    # Parse wt, position, mutant from mutation string like "M1A"
    merged["wt"]  = merged.mutations.str[0]
    merged["i"]   = merged.mutations.str[1:-1].astype(int)
    merged["mut"] = merged.mutations.str[-1]
    print(f"  Loaded {len(files)} seed files for {subtype}: {len(merged)} mutations")
    return merged[["wt", "i", "mut", "evol_indices", "evol_indices_std"]]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subtype", required=True, help="e.g. RSVA_F or RSVB_F")
    ap.add_argument("--evol_indices_pattern", required=True,
                    help="Glob pattern for evol indices CSVs, e.g. '.../*_seed*.csv'")
    ap.add_argument("--wcn_csv", required=True, help="WCN CSV from compute_wcn.py")
    ap.add_argument("--dissimilarity_csv",
                    default="EVEscape/data/aa_properties/dissimilarity_metrics.csv")
    ap.add_argument(
        "--reference_cds_fasta",
        default=str(default_reference_cds_fasta()),
        help="Subtype reference CDS FASTA containing KX858757.1 and KX858756.1.",
    )
    ap.add_argument("--out_csv", required=True)
    args = ap.parse_args()

    print(f"=== EVEscape scoring: {args.subtype} ===")

    # 1. Load EVE fitness (averaged across seeds)
    eve = load_evol_indices(args.evol_indices_pattern, args.subtype)

    # 2. Load WCN accessibility
    wcn = pd.read_csv(args.wcn_csv)
    print(f"  WCN: {len(wcn)} positions, {wcn.wcn_sc.notna().sum()} mapped from PDB")

    # 3. Load dissimilarity metrics
    props = pd.read_csv(args.dissimilarity_csv)
    # Standardize absolute differences (same as EVEscape sklearn StandardScaler)
    for col in ["eisenberg_weiss_diff", "charge_diff"]:
        vals = props[col].abs().values
        props[col + "_std"] = (vals - vals.mean()) / vals.std()
    props["charge_ew-hydro"] = props["eisenberg_weiss_diff_std"] + props["charge_diff_std"]

    # 4. Merge everything
    df = eve.copy()
    df = df.merge(
        wcn[["pos", "wcn_fill_r", "wcn_sc"]].rename(columns={"pos": "i"}),
        on="i",
        how="left",
    )
    df = df.merge(props[["wt", "mut", "charge_ew-hydro"]], on=["wt", "mut"], how="left")

    n_before = len(df)
    df["structure_mapped"] = df["wcn_sc"].notna()
    # Mean-impute any remaining NaN (matches EVEscape SimpleImputer behaviour)
    for col in ["evol_indices", "wcn_fill_r", "charge_ew-hydro"]:
        n_na = df[col].isna().sum()
        if n_na > 0:
            print(f"  Imputing {n_na} NaN in {col} with column mean")
            df[col] = df[col].fillna(df[col].mean())

    # 5. Compute EVEscape
    df["fitness_z"] = standardize(df["evol_indices"])
    df["accessibility_z"] = standardize(df["wcn_fill_r"])
    df["dissimilarity_z"] = standardize(df["charge_ew-hydro"])
    df["evescape"] = (
        np.log(logistic(df["fitness_z"] / TEMPERATURES["fitness"])) +
        np.log(logistic(df["accessibility_z"] / TEMPERATURES["surfacc"])) +
        np.log(logistic(df["dissimilarity_z"] / TEMPERATURES["exchangability"]))
    )
    df["evescape_composite"] = np.exp(df["evescape"])
    df["evescape_percentile"] = df["evescape"].rank(method="average", pct=True)
    df = annotate_codon_distance(df, args.subtype, Path(args.reference_cds_fasta))
    df["evescape_percentile_nt_distance_global"] = (
        df.groupby("min_nt_distance_from_reference_codon")["evescape"].rank(method="average", pct=True)
    )
    df["evescape_percentile_nt_distance_local"] = (
        df.groupby(["i", "min_nt_distance_from_reference_codon"])["evescape"].rank(method="average", pct=True)
    )
    df["evescape_percentile_one_nt_global"] = np.nan
    one_nt_mask = df["one_nt_accessible_from_reference"]
    df.loc[one_nt_mask, "evescape_percentile_one_nt_global"] = (
        df.loc[one_nt_mask, "evescape"].rank(method="average", pct=True)
    )
    df["evescape_percentile_one_nt_local"] = np.nan
    if one_nt_mask.any():
        df.loc[one_nt_mask, "evescape_percentile_one_nt_local"] = (
            df.loc[one_nt_mask].groupby("i")["evescape"].rank(method="average", pct=True)
        )

    # Rename for clarity (matches EVEscape convention)
    df = df.rename(columns={
        "evol_indices": "fitness_eve",
        "wcn_fill_r": "accessibility_wcn",
        "charge_ew-hydro": "dissimilarity_charge_hydro",
    })
    df["mutation"] = df["wt"] + df["i"].astype(int).astype(str) + df["mut"]

    out_cols = [
        "mutation",
        "wt",
        "i",
        "mut",
        "fitness_eve",
        "evol_indices_std",
        "fitness_z",
        "accessibility_wcn",
        "accessibility_z",
        "structure_mapped",
        "dissimilarity_charge_hydro",
        "dissimilarity_z",
        "evescape",
        "evescape_composite",
        "evescape_percentile",
        "reference_codon",
        "min_nt_distance_from_reference_codon",
        "one_nt_accessible_from_reference",
        "n_one_nt_accessible_aas",
        "one_nt_accessible_aas",
        "evescape_percentile_nt_distance_global",
        "evescape_percentile_nt_distance_local",
        "evescape_percentile_one_nt_global",
        "evescape_percentile_one_nt_local",
    ]
    df = df[out_cols].round(7)
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv, index=False)
    print(f"  Wrote {len(df)} rows → {args.out_csv}")
    print(f"  EVEscape range: [{df.evescape.min():.4f}, {df.evescape.max():.4f}]")
    print(
        "  Structure coverage:",
        f"{int(df.structure_mapped.sum()) // 19} positions mapped,"
        f" {df.structure_mapped.mean() * 100:.1f}% of single mutants directly mapped",
    )

if __name__ == "__main__":
    main()
