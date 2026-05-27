from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path


GENETIC_CODE = {
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

IUPAC_BASES = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}

NT_GAP_SCORE = -2
AA_MATCH_SCORE = 2
AA_MISMATCH_SCORE = -1
AA_GAP_SCORE = -2


@dataclass
class LocalMatch:
    score: int
    ref_start: int
    ref_end: int
    target_start: int
    target_end: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find the HCV1NS3 nucleotide range inside na.fasta, translate both "
            "sequences in-frame, and write a codon-level amino-acid alignment."
        )
    )
    parser.add_argument(
        "--ns3",
        default="Genotype-Ref.fasta",
        help="Input HCV1NS3 nucleotide FASTA. Default: Genotype-Ref.fasta.",
    )
    parser.add_argument("--na", default="na.fasta", help="Input nucleotide FASTA.")
    parser.add_argument(
        "--out",
        default="ns3_codon_alignment.txt",
        help="Output text alignment file.",
    )
    parser.add_argument(
        "--width",
        type=int,
        default=20,
        help="Amino-acid columns per alignment block. Default: 20.",
    )
    return parser.parse_args()


def read_single_fasta(path: Path) -> tuple[str, str]:
    header = ""
    sequence_parts: list[str] = []

    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                raise ValueError(f"{path} contains more than one FASTA record")
            header = line[1:].strip()
        else:
            sequence_parts.append(line)

    if not header:
        raise ValueError(f"{path} does not contain a FASTA header")
    return header, "".join(sequence_parts).upper().replace("U", "T")


def nucleotide_score(left: str, right: str) -> int:
    left_bases = set(IUPAC_BASES.get(left, ""))
    right_bases = set(IUPAC_BASES.get(right, ""))
    if not left_bases or not right_bases:
        return -1
    if left == right and len(left_bases) == 1:
        return 2
    if left_bases & right_bases:
        return 1
    return -1


def local_match_range(reference: str, target: str) -> LocalMatch:
    previous_scores = [0] * (len(target) + 1)
    previous_ref_starts = [0] * (len(target) + 1)
    previous_target_starts = [0] * (len(target) + 1)
    best = LocalMatch(0, 0, 0, 0, 0)

    for ref_index, ref_base in enumerate(reference, start=1):
        current_scores = [0] * (len(target) + 1)
        current_ref_starts = [0] * (len(target) + 1)
        current_target_starts = [0] * (len(target) + 1)

        for target_index, target_base in enumerate(target, start=1):
            diagonal = previous_scores[target_index - 1] + nucleotide_score(
                ref_base, target_base
            )
            up = previous_scores[target_index] + NT_GAP_SCORE
            left = current_scores[target_index - 1] + NT_GAP_SCORE
            score = max(0, diagonal, up, left)

            if score == 0:
                continue
            if score == diagonal:
                if previous_scores[target_index - 1] == 0:
                    start_ref = ref_index
                    start_target = target_index
                else:
                    start_ref = previous_ref_starts[target_index - 1]
                    start_target = previous_target_starts[target_index - 1]
            elif score == up:
                start_ref = previous_ref_starts[target_index] or ref_index
                start_target = previous_target_starts[target_index] or target_index
            else:
                start_ref = current_ref_starts[target_index - 1] or ref_index
                start_target = current_target_starts[target_index - 1] or target_index

            current_scores[target_index] = score
            current_ref_starts[target_index] = start_ref
            current_target_starts[target_index] = start_target

            if score > best.score:
                best = LocalMatch(
                    score=score,
                    ref_start=start_ref,
                    ref_end=ref_index,
                    target_start=start_target,
                    target_end=target_index,
                )

        previous_scores = current_scores
        previous_ref_starts = current_ref_starts
        previous_target_starts = current_target_starts

    if best.score == 0:
        raise ValueError("Could not find a local nucleotide match")
    return best


def trim_to_codons(sequence: str) -> str:
    return sequence[: len(sequence) - (len(sequence) % 3)]


def extend_match_to_codons(
    match: LocalMatch, reference_length: int, target_length: int
) -> LocalMatch:
    ref_length = match.ref_end - match.ref_start + 1
    target_match_length = match.target_end - match.target_start + 1
    extra = max((-ref_length) % 3, (-target_match_length) % 3)
    if extra == 0:
        return match

    return LocalMatch(
        score=match.score,
        ref_start=match.ref_start,
        ref_end=min(reference_length, match.ref_end + extra),
        target_start=match.target_start,
        target_end=min(target_length, match.target_end + extra),
    )


def codons_from_dna(sequence: str) -> list[str]:
    coding_sequence = trim_to_codons(sequence.upper().replace("U", "T"))
    return [
        coding_sequence[index : index + 3]
        for index in range(0, len(coding_sequence), 3)
    ]


def translate_codon(codon: str) -> str:
    possible_codons = [""]
    for base in codon:
        possible_bases = IUPAC_BASES.get(base)
        if possible_bases is None:
            return "X"
        possible_codons = [
            prefix + possible_base
            for prefix in possible_codons
            for possible_base in possible_bases
        ]

    amino_acids = {GENETIC_CODE[possible_codon] for possible_codon in possible_codons}
    if len(amino_acids) == 1:
        return amino_acids.pop()
    return "X"


def translate_codons(codons: list[str]) -> str:
    return "".join(translate_codon(codon) for codon in codons)


def amino_acid_score(left: str, right: str) -> int:
    return AA_MATCH_SCORE if left == right else AA_MISMATCH_SCORE


def nw_score_row(left: str, right: str) -> list[int]:
    previous = [column * AA_GAP_SCORE for column in range(len(right) + 1)]

    for row_index, left_char in enumerate(left, start=1):
        current = [row_index * AA_GAP_SCORE]
        for column_index, right_char in enumerate(right, start=1):
            diagonal = previous[column_index - 1] + amino_acid_score(
                left_char, right_char
            )
            up = previous[column_index] + AA_GAP_SCORE
            left_score = current[column_index - 1] + AA_GAP_SCORE
            current.append(max(diagonal, up, left_score))
        previous = current

    return previous


def needleman_wunsch_small(left: str, right: str) -> tuple[str, str]:
    rows = len(left) + 1
    columns = len(right) + 1
    scores = [[0] * columns for _ in range(rows)]

    for row in range(1, rows):
        scores[row][0] = row * AA_GAP_SCORE
    for column in range(1, columns):
        scores[0][column] = column * AA_GAP_SCORE

    for row in range(1, rows):
        for column in range(1, columns):
            diagonal = scores[row - 1][column - 1] + amino_acid_score(
                left[row - 1], right[column - 1]
            )
            up = scores[row - 1][column] + AA_GAP_SCORE
            left_score = scores[row][column - 1] + AA_GAP_SCORE
            scores[row][column] = max(diagonal, up, left_score)

    aligned_left: list[str] = []
    aligned_right: list[str] = []
    row = len(left)
    column = len(right)

    while row > 0 or column > 0:
        if row > 0 and column > 0:
            diagonal = scores[row - 1][column - 1] + amino_acid_score(
                left[row - 1], right[column - 1]
            )
            if scores[row][column] == diagonal:
                aligned_left.append(left[row - 1])
                aligned_right.append(right[column - 1])
                row -= 1
                column -= 1
                continue

        if row > 0 and scores[row][column] == scores[row - 1][column] + AA_GAP_SCORE:
            aligned_left.append(left[row - 1])
            aligned_right.append("-")
            row -= 1
        else:
            aligned_left.append("-")
            aligned_right.append(right[column - 1])
            column -= 1

    return "".join(reversed(aligned_left)), "".join(reversed(aligned_right))


def hirschberg(left: str, right: str) -> tuple[str, str]:
    if not left:
        return "-" * len(right), right
    if not right:
        return left, "-" * len(left)
    if len(left) == 1 or len(right) == 1:
        return needleman_wunsch_small(left, right)

    midpoint = len(left) // 2
    left_score = nw_score_row(left[:midpoint], right)
    right_score = nw_score_row(left[midpoint:][::-1], right[::-1])
    split = max(
        range(len(right) + 1),
        key=lambda index: left_score[index] + right_score[len(right) - index],
    )

    aligned_left_a, aligned_right_a = hirschberg(left[:midpoint], right[:split])
    aligned_left_b, aligned_right_b = hirschberg(left[midpoint:], right[split:])
    return aligned_left_a + aligned_left_b, aligned_right_a + aligned_right_b


def aligned_codons(codons: list[str], aligned_amino_acids: str) -> list[str]:
    codon_index = 0
    output: list[str] = []

    for amino_acid in aligned_amino_acids:
        if amino_acid == "-":
            output.append("---")
            continue
        output.append(codons[codon_index])
        codon_index += 1

    return output


def marker_line(aligned_ref: str, aligned_target: str) -> str:
    markers: list[str] = []
    for ref_aa, target_aa in zip(aligned_ref, aligned_target):
        if ref_aa == "-" or target_aa == "-":
            markers.append(" ")
        elif ref_aa == target_aa:
            markers.append("|")
        else:
            markers.append(".")
    return "".join(markers)


def spaced_amino_acids(sequence: str) -> str:
    return " ".join(f" {amino_acid} " for amino_acid in sequence)


def spaced_markers(markers: str) -> str:
    return " ".join(f" {marker} " for marker in markers)


def alignment_stats(aligned_ref: str, aligned_target: str) -> tuple[int, int, int]:
    matches = 0
    mismatches = 0
    gaps = 0
    for ref_aa, target_aa in zip(aligned_ref, aligned_target):
        if ref_aa == "-" or target_aa == "-":
            gaps += 1
        elif ref_aa == target_aa:
            matches += 1
        else:
            mismatches += 1
    return matches, mismatches, gaps


def format_alignment(
    ns3_header: str,
    na_header: str,
    match: LocalMatch,
    ns3_sequence: str,
    matched_na_sequence: str,
    ns3_codons: list[str],
    na_codons: list[str],
    aligned_ns3: str,
    aligned_na: str,
    width: int,
) -> str:
    matches, mismatches, gaps = alignment_stats(aligned_ns3, aligned_na)
    identity = matches / len(aligned_ns3) * 100 if aligned_ns3 else 0.0
    ns3_aligned_codons = aligned_codons(ns3_codons, aligned_ns3)
    na_aligned_codons = aligned_codons(na_codons, aligned_na)
    markers = marker_line(aligned_ns3, aligned_na)

    lines = [
        "HCV1NS3 codon alignment summary",
        f"  HCV1NS3 FASTA: {ns3_header}",
        f"  na FASTA: {na_header}",
        f"  local nucleotide score: {match.score}",
        f"  HCV1NS3 matched range: {match.ref_start}..{match.ref_end}",
        f"  na.fasta matched range: {match.target_start}..{match.target_end}",
        f"  HCV1NS3 nucleotide length used: {len(ns3_sequence)} nt",
        f"  na nucleotide length used: {len(matched_na_sequence)} nt",
        f"  HCV1NS3 translated length: {len(ns3_codons)} aa",
        f"  na translated length: {len(na_codons)} aa",
        f"  amino-acid alignment length: {len(aligned_ns3)} aa",
        f"  matches: {matches}",
        f"  mismatches: {mismatches}",
        f"  amino-acid gaps: {gaps}",
        f"  identity: {identity:.2f}%",
        "",
    ]

    for start in range(0, len(aligned_ns3), width):
        end = start + width
        lines.append(f"NS3 codon {' '.join(ns3_aligned_codons[start:end])}")
        lines.append(f"NS3 aa    {spaced_amino_acids(aligned_ns3[start:end])}")
        lines.append(f"match     {spaced_markers(markers[start:end])}")
        lines.append(f"na aa     {spaced_amino_acids(aligned_na[start:end])}")
        lines.append(f"na codon  {' '.join(na_aligned_codons[start:end])}")
        lines.append("")

    return "\n".join(lines)


def entry() -> None:
    args = parse_args()
    ns3_header, ns3_sequence = read_single_fasta(Path(args.ns3))
    na_header, na_sequence = read_single_fasta(Path(args.na))

    raw_match = local_match_range(ns3_sequence, na_sequence)
    match = extend_match_to_codons(raw_match, len(ns3_sequence), len(na_sequence))
    ns3_matched = trim_to_codons(ns3_sequence[match.ref_start - 1 : match.ref_end])
    na_matched = trim_to_codons(na_sequence[match.target_start - 1 : match.target_end])

    ns3_codons = codons_from_dna(ns3_matched)
    na_codons = codons_from_dna(na_matched)
    ns3_protein = translate_codons(ns3_codons)
    na_protein = translate_codons(na_codons)
    aligned_ns3, aligned_na = hirschberg(ns3_protein, na_protein)

    output = format_alignment(
        ns3_header=ns3_header,
        na_header=na_header,
        match=match,
        ns3_sequence=ns3_matched,
        matched_na_sequence=na_matched,
        ns3_codons=ns3_codons,
        na_codons=na_codons,
        aligned_ns3=aligned_ns3,
        aligned_na=aligned_na,
        width=args.width,
    )
    Path(args.out).write_text(f"{output}\n")
    print(f"Wrote NS3 codon alignment to {args.out}")


if __name__ == "__main__":
    entry()
