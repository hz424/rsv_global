#!/usr/bin/env python3
"""Make an EVE-compatible MSA by inserting a reference sequence as the focus sequence.

EVE's MSA_processing expects:
  * focus sequence is the FIRST record
  * focus header is of the form: >FOCUS_NAME/START-STOP

We already have a gapped MSA (e.g., MAFFT output). To avoid requiring external
alignment binaries at runtime, we create a gapped focus sequence by *projecting*
the ungapped reference sequence onto the gap-pattern of an "anchor" sequence
present in the MSA:

  - Choose an anchor record whose gap-stripped length equals len(reference)
    and that contains only AA20 and '-' (default; optionally allow 'X' for
    evaluation-only MSAs).
  - Construct the focus sequence by walking across the anchor:
      * if anchor has '-', focus has '-'
      * else consume the next residue from the reference

This keeps the alignment length identical to the input MSA, while ensuring the
focus sequence (after stripping '-') matches the reference exactly.
"""

from __future__ import annotations

import argparse
from pathlib import Path


AA20 = set("ACDEFGHIKLMNPQRSTVWY")


def read_single_fasta_seq(path: Path) -> str:
    text = path.read_text().splitlines()
    seq = []
    headers = 0
    for line in text:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            headers += 1
            continue
        seq.append(line)
    if headers != 1:
        raise ValueError(f"Expected 1 FASTA record in {path}, found {headers}")
    s = "".join(seq).upper()
    if "-" in s:
        raise ValueError(f"Reference FASTA must be ungapped; found '-' in {path}")
    bad = set(s) - AA20
    if bad:
        raise ValueError(f"Reference contains non-AA20 letters {sorted(bad)} in {path}")
    return s


def parse_fasta_records(fasta_text: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = None
    seq: list[str] = []
    for raw in fasta_text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq).upper()))
            header = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)
    if header is not None:
        records.append((header, "".join(seq).upper()))
    if not records:
        raise ValueError("No FASTA records parsed from MAFFT output")
    return records


def load_fasta_file(path: Path) -> list[tuple[str, str]]:
    text = path.read_text()
    return parse_fasta_records(text)


def pick_anchor(records: list[tuple[str, str]], ref: str, allow_X: bool = False) -> tuple[str, str, int]:
    """Pick an anchor record used to project reference gaps.

    Criteria:
      - default: only AA20/'-'
      - if allow_X: also allow 'X' in the anchor (for evaluation-only MSAs)
      - ungapped length equals len(ref)
      - prefer highest identity to ref (fewest mismatches), then fewest X, then fewest gaps
    """
    ref_len = len(ref)
    allowed = AA20 | {"-"} | ({"X"} if allow_X else set())
    best = None
    for h, s in records:
        s = s.upper()
        bad = set(s) - allowed
        if bad:
            continue
        ungapped = s.replace("-", "")
        if len(ungapped) != ref_len:
            continue
        # When allow_X=True, treat 'X' as missing/unknown, not a hard mismatch.
        mismatches = sum(1 for a, b in zip(ungapped, ref, strict=True) if a != "X" and a != b)
        n_x = ungapped.count("X")
        gaps = s.count("-")
        key = (mismatches, n_x, gaps)
        if best is None or key < best[0]:
            best = (key, h, s, mismatches)
    if best is None:
        raise ValueError(
            "Could not find an anchor record with the allowed alphabet and ungapped length equal to reference. "
            "This MSA may be heavily gapped/ambiguous, or the reference scope/subtype may not match."
        )
    (_key, h, s, mismatches) = best
    return h, s, mismatches


def project_reference_to_anchor_gaps(ref: str, anchor_gapped: str) -> str:
    ref_i = 0
    out = []
    for ch in anchor_gapped:
        if ch == "-":
            out.append("-")
        else:
            if ref_i >= len(ref):
                raise ValueError("Anchor has more non-gap positions than reference length")
            out.append(ref[ref_i])
            ref_i += 1
    if ref_i != len(ref):
        raise ValueError(
            f"Did not consume all reference residues while projecting to anchor gaps (consumed {ref_i}/{len(ref)})."
        )
    return "".join(out)


def write_fasta(records: list[tuple[str, str]], out_path: Path, wrap: int = 80) -> None:
    out_lines: list[str] = []
    for h, s in records:
        out_lines.append(">" + h)
        for i in range(0, len(s), wrap):
            out_lines.append(s[i : i + wrap])
    out_path.write_text("\n".join(out_lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--msa_in", required=True, help="Input gapped MSA (FASTA)")
    ap.add_argument("--reference_fasta", required=True, help="Ungapped reference protein FASTA (single record)")
    ap.add_argument("--focus_name", required=True, help="Focus name used in EVE header (e.g., RSVB_F)")
    ap.add_argument("--start", type=int, default=1)
    ap.add_argument("--stop", type=int, required=True)
    ap.add_argument("--out_msa", required=True, help="Output EVE-ready MSA (FASTA)")
    ap.add_argument(
        "--allow_X",
        action="store_true",
        help=(
            "Allow 'X' in the input MSA when selecting an anchor sequence for gap-projection. "
            "Use for evaluation-only MSAs that contain ambiguous residues. Note: EVE itself expects AA20+'-' "
            "and will not load an MSA containing X."
        ),
    )
    args = ap.parse_args()

    msa_in = Path(args.msa_in)
    reference_fasta = Path(args.reference_fasta)
    out_msa = Path(args.out_msa)

    ref = read_single_fasta_seq(reference_fasta)
    records = load_fasta_file(msa_in)

    anchor_header, anchor_seq, anchor_mismatches = pick_anchor(records, ref=ref, allow_X=bool(args.allow_X))
    ref_seq = project_reference_to_anchor_gaps(ref, anchor_seq)
    if ref_seq.replace("-", "") != ref:
        raise AssertionError("Internal error: projected focus does not match reference after stripping gaps")

    focus_header = f"{args.focus_name}/{args.start}-{args.stop}"

    # Output: focus first, then all original sequences in input order.
    out_records: list[tuple[str, str]] = [(focus_header, ref_seq)] + records

    out_msa.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(out_records, out_msa)

    print(f"Wrote EVE-ready MSA: {out_msa}")
    print(f"Focus header: >{focus_header}")
    print(
        f"Anchor used to project gaps: >{anchor_header} (gaps={anchor_seq.count('-')}, mismatches_to_ref={anchor_mismatches})"
    )


if __name__ == "__main__":
    main()
