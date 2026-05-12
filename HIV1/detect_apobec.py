#!/usr/bin/env python3
"""
Detect APOBEC-like HIV mutations from an aligned FASTA file.

Heuristic used:
- APOBEC3G-like: reference motif GG, observed AG (G->A at 2nd base)
- APOBEC3F-like: reference motif GA, observed AA (G->A at 2nd base)

The script compares each query sequence to a reference sequence and counts
G->A changes by context.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Set


VALID_BASES = set("ACGTN-")


@dataclass
class SeqRecord:
    header: str
    sequence: str

    @property
    def seq_id(self) -> str:
        return self.header.split()[0]


def read_fasta(path: Path) -> List[SeqRecord]:
    records: List[SeqRecord] = []
    header = None
    seq_chunks: List[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_chunks).upper()
                    records.append(SeqRecord(header=header, sequence=seq))
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if header is not None:
            seq = "".join(seq_chunks).upper()
            records.append(SeqRecord(header=header, sequence=seq))

    if not records:
        raise ValueError("No FASTA records found.")

    for rec in records:
        bad = set(rec.sequence) - VALID_BASES
        if bad:
            raise ValueError(
                f"Sequence '{rec.seq_id}' has invalid bases: {''.join(sorted(bad))}"
            )
    return records


def ensure_aligned(records: List[SeqRecord]) -> int:
    lengths = {len(r.sequence) for r in records}
    if len(lengths) != 1:
        raise ValueError(
            "Sequences are not the same length. Provide an aligned FASTA for mutation detection."
        )
    return lengths.pop()


def pick_reference(records: List[SeqRecord], reference_id: str | None) -> SeqRecord:
    if reference_id is None:
        return records[0]

    for rec in records:
        if rec.seq_id == reference_id:
            return rec
    available = ", ".join(r.seq_id for r in records[:10])
    raise ValueError(
        f"Reference id '{reference_id}' not found. Available ids include: {available}"
    )


def context_class(prev_base: str) -> str:
    if prev_base == "G":
        return "apobec3g_like"
    if prev_base == "A":
        return "apobec3f_like"
    return "other_ga"


def parse_position_spec(spec: str | None) -> Set[int]:
    """
    Parse 1-based position spec like: "67,103,184-190".
    Returns a set of 0-based nucleotide indices.
    """
    if not spec:
        return set()

    out: Set[int] = set()
    parts = [p.strip() for p in spec.split(",") if p.strip()]
    for part in parts:
        if "-" in part:
            a_str, b_str = part.split("-", 1)
            a = int(a_str)
            b = int(b_str)
            if a < 1 or b < 1 or b < a:
                raise ValueError(f"Invalid range in keep_pattern_for_drug: '{part}'")
            for pos in range(a, b + 1):
                out.add(pos - 1)
        else:
            pos = int(part)
            if pos < 1:
                raise ValueError(f"Invalid position in keep_pattern_for_drug: '{part}'")
            out.add(pos - 1)
    return out


def analyze_sequence(
    reference: str, query: str, keep_pattern_for_drug: Set[int]
) -> Dict[str, float]:
    total_compared = 0
    ga_total = 0
    apobec3g_like = 0
    apobec3f_like = 0
    other_ga = 0
    drug_pattern_ga = 0
    non_apobec_changes = 0

    for i in range(1, len(reference)):
        r_prev = reference[i - 1]
        r_base = reference[i]
        q_base = query[i]

        if r_base in {"-", "N"} or q_base in {"-", "N"}:
            continue

        total_compared += 1
        if q_base == r_base:
            continue

        if r_base == "G" and q_base == "A":
            ga_total += 1
            if i in keep_pattern_for_drug:
                drug_pattern_ga += 1
                continue
            cls = context_class(r_prev)
            if cls == "apobec3g_like":
                apobec3g_like += 1
            elif cls == "apobec3f_like":
                apobec3f_like += 1
            else:
                other_ga += 1
        else:
            non_apobec_changes += 1

    apobec_like = apobec3g_like + apobec3f_like
    mutation_total = ga_total + non_apobec_changes

    return {
        "positions_compared": total_compared,
        "all_mutations": mutation_total,
        "ga_mutations_total": ga_total,
        "apobec3g_like": apobec3g_like,
        "apobec3f_like": apobec3f_like,
        "apobec_like_total": apobec_like,
        "drug_pattern_ga": drug_pattern_ga,
        "other_ga": other_ga,
        "non_apobec_mutations": non_apobec_changes,
        "apobec_fraction_of_all_mutations": (
            apobec_like / mutation_total if mutation_total else 0.0
        ),
        "apobec_fraction_of_ga": (apobec_like / ga_total if ga_total else 0.0),
    }


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Detect APOBEC-like HIV mutation patterns from aligned FASTA."
    )
    parser.add_argument("fasta", type=Path, help="Aligned FASTA file.")
    parser.add_argument(
        "--reference-id",
        help="Sequence id to use as reference (default: first FASTA sequence).",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        help="Optional output CSV path for machine-readable results.",
    )
    parser.add_argument(
        "--keep-pattern-for-drug",
        default="",
        help=(
            "Comma-separated 1-based nucleotide positions/ranges to keep as "
            "drug-pattern G->A (excluded from APOBEC-like counts), e.g. 67,103,184-190"
        ),
    )
    args = parser.parse_args()

    try:
        records = read_fasta(args.fasta)
        ensure_aligned(records)
        ref = pick_reference(records, args.reference_id)
        keep_pattern_for_drug = parse_position_spec(args.keep_pattern_for_drug)
    except Exception as exc:  # noqa: BLE001
        print(f"Error: {exc}", file=sys.stderr)
        return 2

    rows: List[Dict[str, object]] = []
    print(f"Reference: {ref.seq_id}")
    print(
        "sequence_id\tpositions_compared\tall_mutations\tga_mutations_total\t"
        "apobec3g_like\tapobec3f_like\tapobec_like_total\tdrug_pattern_ga\tother_ga\t"
        "non_apobec_mutations\tapobec_fraction_of_all_mutations\t"
        "apobec_fraction_of_ga"
    )

    for rec in records:
        if rec.seq_id == ref.seq_id:
            continue

        result = analyze_sequence(ref.sequence, rec.sequence, keep_pattern_for_drug)
        row = {"sequence_id": rec.seq_id, **result}
        rows.append(row)
        print(
            f"{row['sequence_id']}\t{row['positions_compared']}\t{row['all_mutations']}\t"
            f"{row['ga_mutations_total']}\t{row['apobec3g_like']}\t"
            f"{row['apobec3f_like']}\t{row['apobec_like_total']}\t"
            f"{row['drug_pattern_ga']}\t{row['other_ga']}\t{row['non_apobec_mutations']}\t"
            f"{row['apobec_fraction_of_all_mutations']:.4f}\t"
            f"{row['apobec_fraction_of_ga']:.4f}"
        )

    if args.csv:
        write_csv(args.csv, rows)
        print(f"\nWrote CSV: {args.csv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
