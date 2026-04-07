#!/usr/bin/env python3
"""Shared helpers for RSV F EVEscape analysis."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

AA20 = set("ACDEFGHIKLMNPQRSTVWY")
DNA_BASES = "ACGT"
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
REFERENCE_ACCESSION_BY_SUBTYPE = {"RSVA_F": "KX858757.1", "RSVB_F": "KX858756.1"}

ANTIGENIC_SITE_DATA = {
    "Site Ø": {"ranges": [(62, 96), (195, 227)], "color": "#FFB6C1"},
    "Site I": {"ranges": [(27, 45), (312, 318), (378, 389)], "color": "#ADD8E6"},
    "Site II": {"ranges": [(254, 277)], "color": "#90EE90"},
    "Site III": {"ranges": [(46, 54), (301, 311), (345, 352), (367, 378)], "color": "#FFD700"},
    "Site IV": {"ranges": [(422, 471)], "color": "#DA70D6"},
    "Site V": {"ranges": [(55, 61), (146, 194), (287, 300)], "color": "#FAA460"},
    "P27 (FP)": {"ranges": [(110, 136)], "color": "#778899"},
}

ANTIGENIC_SITES = {name: spec["ranges"] for name, spec in ANTIGENIC_SITE_DATA.items()}

SITE_COLORS = {name: spec["color"] for name, spec in ANTIGENIC_SITE_DATA.items()}
SITE_COLORS["None"] = "#8C8C8C"

COMPONENT_COLUMNS = [
    "fitness_z",
    "accessibility_z",
    "dissimilarity_z",
]

VARIANT_METRIC_COLUMNS = [
    "variant_score_raw_sum",
    "variant_score_raw_mean",
    "variant_score_raw_median",
    "variant_score_raw_max",
    "variant_score_composite_mean",
    "variant_score_composite_max",
    "variant_score_top3_composite_mean",
    "variant_score_mean_percentile",
    "variant_score_max_percentile",
    "variant_score_top3_mean_percentile",
    "variant_score_mean_percentile_nt_distance_global",
    "variant_score_max_percentile_nt_distance_global",
    "variant_score_top3_mean_percentile_nt_distance_global",
    "high_escape_mutation_count_90",
    "high_escape_mutation_fraction_90",
    "high_escape_mutation_count_95",
    "high_escape_mutation_fraction_95",
    "n_antigenic_site_mutations",
    "frac_antigenic_site_mutations",
    "n_structure_mapped_mutations",
    "frac_structure_mapped_mutations",
    "n_one_nt_accessible_mutations",
    "frac_one_nt_accessible_mutations",
    "mean_min_nt_distance_from_reference_codon",
]


def mutation_key(wt: str, pos: int, mut: str) -> str:
    return f"{wt}{int(pos)}{mut}"


def antigenic_sites_for_position(pos: int) -> str:
    hits = [site for site, ranges in ANTIGENIC_SITES.items() for lo, hi in ranges if lo <= pos <= hi]
    ordered_hits = list(dict.fromkeys(hits))
    return ",".join(ordered_hits) if ordered_hits else "None"


def primary_antigenic_site(pos: int) -> str:
    sites = antigenic_sites_for_position(pos)
    return sites.split(",")[0] if sites != "None" else "None"


def safe_float(value, default: float = 0.0) -> float:
    if value is None:
        return default
    try:
        if pd.isna(value):
            return default
    except TypeError:
        pass
    return float(value)


def parse_mutations_csv(muts_csv: object) -> list[str]:
    if muts_csv is None:
        return []
    if isinstance(muts_csv, float) and pd.isna(muts_csv):
        return []
    return [m.strip() for m in str(muts_csv).split(",") if m.strip()]


def read_selected_reference_cds(path: Path, subtype: str) -> str:
    accession = REFERENCE_ACCESSION_BY_SUBTYPE[subtype]
    header = None
    seq_parts: list[str] = []
    for line in Path(path).read_text().splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if header is not None and accession in header:
                return "".join(seq_parts).upper()
            header = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(line.strip())
    if header is not None and accession in header:
        return "".join(seq_parts).upper()
    raise ValueError(f"Could not find reference CDS for {subtype} ({accession}) in {path}")


def hamming_distance(a: str, b: str) -> int:
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))


def amino_acid_min_nt_distance(ref_codon: str) -> dict[str, int]:
    wt = CODON_TABLE.get(ref_codon)
    distances: dict[str, int] = {}
    for codon, aa in CODON_TABLE.items():
        if aa in {wt, "*"}:
            continue
        dist = hamming_distance(ref_codon, codon)
        if aa not in distances or dist < distances[aa]:
            distances[aa] = dist
    return distances


def accessible_amino_acids_by_distance(ref_codon: str) -> dict[int, list[str]]:
    dist_map = amino_acid_min_nt_distance(ref_codon)
    grouped: dict[int, list[str]] = {}
    for aa, dist in dist_map.items():
        grouped.setdefault(dist, []).append(aa)
    return {dist: sorted(aas) for dist, aas in grouped.items()}


def standardize(values: Iterable[float]) -> np.ndarray:
    arr = np.asarray(list(values), dtype=float)
    std = arr.std(ddof=0)
    if std == 0 or not np.isfinite(std):
        return np.zeros_like(arr)
    return (arr - arr.mean()) / std


def ensure_score_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "mutation" not in out.columns:
        out["mutation"] = out["wt"].astype(str) + out["i"].astype(int).astype(str) + out["mut"].astype(str)
    if "evescape_composite" not in out.columns:
        out["evescape_composite"] = np.exp(out["evescape"].astype(float))
    if "evescape_percentile" not in out.columns:
        out["evescape_percentile"] = out["evescape"].rank(method="average", pct=True)
    if "antigenic_sites" not in out.columns:
        out["antigenic_sites"] = out["i"].astype(int).map(antigenic_sites_for_position)
    if "primary_antigenic_site" not in out.columns:
        out["primary_antigenic_site"] = out["i"].astype(int).map(primary_antigenic_site)
    if "structure_mapped" not in out.columns:
        out["structure_mapped"] = True
    return out


def build_score_lookup(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, dict[str, object]]]:
    scores = ensure_score_columns(df)
    records = scores.to_dict(orient="records")
    return scores, {row["mutation"]: row for row in records}


def _mean_or_zero(values: list[float]) -> float:
    return float(np.mean(values)) if values else 0.0


def _median_or_zero(values: list[float]) -> float:
    return float(np.median(values)) if values else 0.0


def _max_or_zero(values: list[float]) -> float:
    return float(np.max(values)) if values else 0.0


def summarize_variant_from_mutations(mutations: list[str], lookup: dict[str, dict[str, object]]) -> dict[str, object]:
    matched = [lookup[m] for m in mutations if m in lookup]
    missing = [m for m in mutations if m not in lookup]

    raw_scores = [safe_float(r.get("evescape")) for r in matched]
    composites = [safe_float(r.get("evescape_composite")) for r in matched]
    percentiles = [safe_float(r.get("evescape_percentile")) for r in matched]
    percentiles_ntdist_global = [
        safe_float(r.get("evescape_percentile_nt_distance_global"))
        for r in matched
        if pd.notna(r.get("evescape_percentile_nt_distance_global"))
    ]
    mapped_flags = [bool(r.get("structure_mapped")) for r in matched]
    antigenic_flags = [r.get("antigenic_sites", "None") != "None" for r in matched]
    one_nt_flags = [bool(r.get("one_nt_accessible_from_reference")) for r in matched]
    min_nt_distances = [
        safe_float(r.get("min_nt_distance_from_reference_codon"))
        for r in matched
        if pd.notna(r.get("min_nt_distance_from_reference_codon"))
    ]

    ranked = sorted(matched, key=lambda row: safe_float(row.get("evescape")), reverse=True)
    top_mutations = [row["mutation"] for row in ranked[:3]]
    top3_composites = [safe_float(row.get("evescape_composite")) for row in ranked[:3]]
    top3_percentiles = [safe_float(row.get("evescape_percentile")) for row in ranked[:3]]
    top3_ntdist_percentiles = [
        safe_float(row.get("evescape_percentile_nt_distance_global"))
        for row in ranked[:3]
        if pd.notna(row.get("evescape_percentile_nt_distance_global"))
    ]
    touched_sites = sorted(
        {
            site
            for row in matched
            for site in str(row.get("antigenic_sites", "None")).split(",")
            if site and site != "None"
        }
    )

    summary = {
        "variant_score_raw_sum": float(np.sum(raw_scores)) if raw_scores else 0.0,
        "variant_score_raw_mean": _mean_or_zero(raw_scores),
        "variant_score_raw_median": _median_or_zero(raw_scores),
        "variant_score_raw_max": _max_or_zero(raw_scores),
        "variant_score_composite_mean": _mean_or_zero(composites),
        "variant_score_composite_max": _max_or_zero(composites),
        "variant_score_top3_composite_mean": _mean_or_zero(top3_composites),
        "variant_score_mean_percentile": _mean_or_zero(percentiles),
        "variant_score_max_percentile": _max_or_zero(percentiles),
        "variant_score_top3_mean_percentile": _mean_or_zero(top3_percentiles),
        "variant_score_mean_percentile_nt_distance_global": _mean_or_zero(percentiles_ntdist_global),
        "variant_score_max_percentile_nt_distance_global": _max_or_zero(percentiles_ntdist_global),
        "variant_score_top3_mean_percentile_nt_distance_global": _mean_or_zero(top3_ntdist_percentiles),
        "high_escape_mutation_count_90": int(sum(p >= 0.90 for p in percentiles)),
        "high_escape_mutation_fraction_90": _mean_or_zero([p >= 0.90 for p in percentiles]),
        "high_escape_mutation_count_95": int(sum(p >= 0.95 for p in percentiles)),
        "high_escape_mutation_fraction_95": _mean_or_zero([p >= 0.95 for p in percentiles]),
        "n_antigenic_site_mutations": int(sum(antigenic_flags)),
        "frac_antigenic_site_mutations": _mean_or_zero(antigenic_flags),
        "n_structure_mapped_mutations": int(sum(mapped_flags)),
        "frac_structure_mapped_mutations": _mean_or_zero(mapped_flags),
        "n_one_nt_accessible_mutations": int(sum(one_nt_flags)),
        "frac_one_nt_accessible_mutations": _mean_or_zero(one_nt_flags),
        "mean_min_nt_distance_from_reference_codon": _mean_or_zero(min_nt_distances),
        "top_escape_mutations_pipe": "|".join(top_mutations),
        "top_escape_mutation": top_mutations[0] if top_mutations else "",
        "antigenic_sites_touched": ",".join(touched_sites) if touched_sites else "None",
        "missing_mutations_csv": ",".join(missing),
    }
    return summary


def iter_fasta(path: Path):
    header = None
    seq_parts: list[str] = []
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts).upper()
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, "".join(seq_parts).upper()


def parse_msa_variants(msa_path: Path) -> pd.DataFrame:
    records = list(iter_fasta(msa_path))
    if not records:
        raise ValueError(f"No FASTA records in {msa_path}")

    ref_id, ref_seq = records[0]
    keep_cols = [idx for idx, aa in enumerate(ref_seq) if aa not in {"-", "."}]
    ref_positions = {}
    ref_pos = 0
    for idx, aa in enumerate(ref_seq):
        if aa in {"-", "."}:
            continue
        ref_pos += 1
        ref_positions[idx] = ref_pos

    rows = []
    for header, seq in records[1:]:
        muts = []
        x_ignored = 0
        for col in keep_cols:
            wt = ref_seq[col]
            alt = seq[col] if col < len(seq) else "-"
            if alt in {"-", "."} or alt == wt:
                continue
            if alt == "X":
                x_ignored += 1
                continue
            if wt not in AA20 or alt not in AA20:
                continue
            muts.append(mutation_key(wt, ref_positions[col], alt))
        rows.append(
            {
                "variant_id": header,
                "n_substitutions": len(muts),
                "has_X": bool(x_ignored > 0),
                "n_X_positions_ignored": int(x_ignored),
                "mutations_csv": ",".join(muts),
                "mutations_colon": ":".join(muts),
                "focus_reference_id": ref_id,
            }
        )
    return pd.DataFrame(rows)


def collect_mutation_counts(mutations_series: pd.Series) -> Counter:
    counts: Counter = Counter()
    for muts_csv in mutations_series:
        counts.update(parse_mutations_csv(muts_csv))
    return counts
