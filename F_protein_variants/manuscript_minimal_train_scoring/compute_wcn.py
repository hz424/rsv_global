#!/usr/bin/env python3
"""Compute Weighted Contact Number (WCN) from an RSV F PDB and map to reference positions.

Uses only stdlib + numpy + pandas (no BioPython). Follows the EVEscape convention:
  - WCN computed over the selected multimer (for example, chains A,B,C) using inverse-squared Cα distances
  - Only chain-A residues mapped to the 574-aa reference via Needleman-Wunsch
  - Missing positions filled with forward/backward fill average
  - Output negated so higher value = more accessible  (wcn_fill_r = -wcn_fill)
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd

AA3 = {"ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H",
       "ILE":"I","LYS":"K","LEU":"L","MET":"M","ASN":"N","PRO":"P","GLN":"Q",
       "ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y"}

def read_fasta(p):
    h, seq = None, []
    for ln in Path(p).read_text().splitlines():
        ln = ln.strip()
        if ln.startswith(">"):
            if h is not None: break
            h = ln[1:].split()[0]
        elif ln: seq.append(ln)
    return "".join(seq).upper()

def parse_pdb_atoms(pdb_path, chains):
    """Return list of dicts with chain, resseq, aa1, coord_ca, sc_coords per residue."""
    residues, seen = [], set()
    current = None
    for ln in Path(pdb_path).read_text().splitlines():
        if not ln.startswith("ATOM"): continue
        ch = ln[21].strip()
        if ch not in chains: continue
        rn = ln[17:20].strip()
        if rn not in AA3: continue
        resseq = int(ln[22:26])
        icode = ln[26].strip()
        key = (ch, resseq, icode)
        aname = ln[12:16].strip()
        x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
        if key not in seen:
            seen.add(key)
            current = {"chain": ch, "resseq": resseq, "aa1": AA3[rn],
                       "coord_ca": None, "sc_coords": []}
            residues.append(current)
        else:
            current = next(r for r in reversed(residues) if (r["chain"], r["resseq"]) == (ch, resseq))
        if aname == "CA":
            current["coord_ca"] = np.array([x, y, z])
        elif aname not in ("C", "N", "O"):
            current["sc_coords"].append(np.array([x, y, z]))
    # finalise sidechain centers
    for r in residues:
        if r["coord_ca"] is None:
            raise RuntimeError(f"Missing CA: chain {r['chain']} res {r['resseq']}")
        if len(r["sc_coords"]) == 0:
            r["sc_center"] = r["coord_ca"]
        else:
            r["sc_center"] = np.mean(r["sc_coords"], axis=0)
    return residues

def compute_wcn(residues):
    """Vectorised WCN using inverse-squared Cα distances."""
    ca = np.array([r["coord_ca"] for r in residues])       # (N,3)
    sc = np.array([r["sc_center"] for r in residues])       # (N,3)
    # pairwise squared distances
    diff_ca = ca[:, None, :] - ca[None, :, :]               # (N,N,3)
    d2_ca = (diff_ca**2).sum(axis=2)                         # (N,N)
    np.fill_diagonal(d2_ca, np.inf)
    diff_sc = sc[:, None, :] - sc[None, :, :]
    d2_sc = (diff_sc**2).sum(axis=2)
    np.fill_diagonal(d2_sc, np.inf)
    wcn_ca = (1.0 / d2_ca).sum(axis=1)
    wcn_sc = (1.0 / d2_sc).sum(axis=1)
    for i, r in enumerate(residues):
        r["wcn_ca"] = wcn_ca[i]
        r["wcn_sc"] = wcn_sc[i]
    return residues

def nw_align(a, b, match=2, mismatch=-1, gap=-2):
    n, m = len(a), len(b)
    S = np.zeros((n+1, m+1), dtype=np.int32)
    B = np.zeros((n+1, m+1), dtype=np.int8)
    S[1:,0] = np.arange(1, n+1)*gap; B[1:,0] = 2
    S[0,1:] = np.arange(1, m+1)*gap; B[0,1:] = 3
    for i in range(1, n+1):
        for j in range(1, m+1):
            s = match if a[i-1]==b[j-1] else mismatch
            vals = (S[i-1,j-1]+s, S[i-1,j]+gap, S[i,j-1]+gap)
            best = max(vals)
            S[i,j] = best
            B[i,j] = vals.index(best)+1
    i, j = n, m
    pairs = []
    while i>0 or j>0:
        d = B[i,j]
        if d==1: pairs.append((i,j)); i-=1; j-=1
        elif d==2: i-=1
        else: j-=1
    pairs.reverse()
    return pairs, S[n,m]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", default="adaptable_workflow/rsv_f/inputs/structures/5UDD.pdb")
    ap.add_argument("--trimer_chains", default="A,B,C")
    ap.add_argument("--target_chain", default="A")
    ap.add_argument("--reference_fasta", required=True)
    ap.add_argument("--out_csv", required=True)
    args = ap.parse_args()

    chains = [c.strip() for c in args.trimer_chains.split(",")]
    ref_seq = read_fasta(args.reference_fasta)
    print(f"Reference length: {len(ref_seq)}")

    # Parse PDB & compute WCN over full trimer
    residues = parse_pdb_atoms(args.pdb, set(chains))
    print(f"PDB residues (trimer): {len(residues)}")
    residues = compute_wcn(residues)

    # Extract target chain
    tgt = [r for r in residues if r["chain"] == args.target_chain]
    tgt_seq = "".join(r["aa1"] for r in tgt)
    print(f"Target chain {args.target_chain}: {len(tgt_seq)} residues")

    # Align reference → target chain
    pairs, score = nw_align(ref_seq, tgt_seq)
    ident = sum(1 for ri, pi in pairs if ref_seq[ri-1]==tgt_seq[pi-1])
    print(f"Alignment: {ident}/{len(pairs)} = {ident/len(pairs)*100:.1f}% identity")

    # Build mapping: ref_pos (1-based) → wcn values
    rows = []
    mapped_wcn = {}
    for ri, pi in pairs:
        r = tgt[pi-1]
        mapped_wcn[ri] = (r["wcn_ca"], r["wcn_sc"])

    # Create full-length table with forward/backward fill
    data = []
    for pos in range(1, len(ref_seq)+1):
        wt = ref_seq[pos-1]
        if pos in mapped_wcn:
            data.append({"pos": pos, "wt": wt, "wcn_ca": mapped_wcn[pos][0],
                         "wcn_sc": mapped_wcn[pos][1]})
        else:
            data.append({"pos": pos, "wt": wt, "wcn_ca": np.nan, "wcn_sc": np.nan})
    df = pd.DataFrame(data)
    df["wcn_bfil"] = df.wcn_sc.bfill()
    df["wcn_ffil"] = df.wcn_sc.ffill()
    df["wcn_fill"] = df[["wcn_ffil","wcn_bfil"]].mean(axis=1)
    df["wcn_fill_r"] = -df["wcn_fill"]   # EVEscape convention: negate
    df = df.drop(columns=["wcn_bfil","wcn_ffil"])

    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv, index=False)
    print(f"Wrote {len(df)} rows → {args.out_csv}")
    n_miss = df.wcn_sc.isna().sum()
    print(f"  Mapped: {len(df)-n_miss}, filled: {n_miss}")

if __name__ == "__main__":
    main()
