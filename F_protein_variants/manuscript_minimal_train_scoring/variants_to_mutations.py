import argparse
import os
import pandas as pd


AA20 = set("ACDEFGHIKLMNPQRSTVWY")


def read_single_fasta(path: str) -> str:
    n_headers = 0
    seq = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                n_headers += 1
                continue
            seq.append(line)
    if n_headers != 1:
        raise ValueError(
            f"Expected exactly 1 FASTA record in reference_fasta, found {n_headers}. "
            "Please provide a single-record protein FASTA as the reference."
        )
    return "".join(seq).strip().upper()


def read_multi_fasta(path: str):
    records = []
    header = None
    seq = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq).upper()))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line)
        if header is not None:
            records.append((header, "".join(seq).upper()))
    return records


def diff_to_mutations(reference: str, variant: str):
    if len(reference) != len(variant):
        raise ValueError(
            f"Length mismatch: reference={len(reference)} variant={len(variant)}. "
            "Please trim/align sequences to the same scope (e.g., ectodomain-only) "
            "or extend this script to support indels."
        )

    n_x_positions = 0
    if any(ch in variant for ch in ("B", "J")):
        raise ValueError(
            "Variant contains ambiguous residues 'B' or 'J'. Per the publication policy, "
            "drop these sequences from evaluation (and never include them in training)."
        )

    muts = []
    for idx, (wt, aa) in enumerate(zip(reference, variant), start=1):
        if wt == aa:
            continue
        if wt == "-" or aa == "-":
            continue
        if aa == "X":
            n_x_positions += 1
            continue
        if aa not in AA20 or wt not in AA20:
            raise ValueError(
                f"Non-AA20 character encountered at position {idx}: wt='{wt}' var='{aa}'. "
                "Please filter or clean sequences to AA20 (and optionally X for evaluation)."
            )
        muts.append(f"{wt}{idx}{aa}")
    return muts, n_x_positions


def main():
    ap = argparse.ArgumentParser(
        description="Convert RSV F variant FASTA sequences into substitution lists vs a reference."
    )
    ap.add_argument("--reference_fasta", required=True)
    ap.add_argument("--variants_fasta", required=True)
    ap.add_argument("--out_csv", required=True)
    args = ap.parse_args()

    out_dir = os.path.dirname(os.path.abspath(args.out_csv))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    ref = read_single_fasta(args.reference_fasta)
    records = read_multi_fasta(args.variants_fasta)

    rows = []
    for header, seq in records:
        muts, n_x = diff_to_mutations(ref, seq)
        rows.append({
            "variant_id": header,
            "n_substitutions": len(muts),
            "has_X": bool(n_x > 0),
            "n_X_positions_ignored": int(n_x),
            "mutations_csv": ",".join(muts),
            "mutations_colon": ":".join(muts),
        })

    df = pd.DataFrame(rows).sort_values(["n_substitutions", "variant_id"], ascending=[False, True])
    df.to_csv(args.out_csv, index=False)


if __name__ == "__main__":
    main()
