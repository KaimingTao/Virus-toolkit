from __future__ import annotations

import argparse
from pathlib import Path


GENETIC_CODE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
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

MATCH_SCORE = 2
MISMATCH_SCORE = -1
GAP_SCORE = -2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Load nucleotide and protein FASTA files, translate the nucleotide "
            "sequence, globally align the translated protein to the input "
            "protein sequence, and save a text alignment."
        )
    )
    parser.add_argument("--na", default="na.fasta", help="Input nucleotide FASTA.")
    parser.add_argument("--aa", default="aa.fasta", help="Input protein FASTA.")
    parser.add_argument(
        "--out",
        default="translation_alignment.txt",
        help="Output text alignment file.",
    )
    parser.add_argument(
        "--width",
        type=int,
        default=25,
        help="Amino-acid columns per alignment block. Default: 25.",
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

    return header, "".join(sequence_parts).upper()


def translate_dna(sequence: str) -> str:
    normalized = sequence.upper().replace("U", "T")
    protein: list[str] = []

    for index in range(0, len(normalized) - 2, 3):
        codon = normalized[index : index + 3]
        protein.append(translate_codon(codon))

    return "".join(protein)


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


def codons_from_dna(sequence: str) -> list[str]:
    normalized = sequence.upper().replace("U", "T")
    return [
        normalized[index : index + 3]
        for index in range(0, len(normalized) - 2, 3)
    ]


def score_pair(left: str, right: str) -> int:
    return MATCH_SCORE if left == right else MISMATCH_SCORE


def nw_score_row(left: str, right: str) -> list[int]:
    previous = [column * GAP_SCORE for column in range(len(right) + 1)]

    for row_index, left_char in enumerate(left, start=1):
        current = [row_index * GAP_SCORE]
        for column_index, right_char in enumerate(right, start=1):
            diagonal = previous[column_index - 1] + score_pair(left_char, right_char)
            up = previous[column_index] + GAP_SCORE
            left_score = current[column_index - 1] + GAP_SCORE
            current.append(max(diagonal, up, left_score))
        previous = current

    return previous


def needleman_wunsch_small(left: str, right: str) -> tuple[str, str]:
    rows = len(left) + 1
    columns = len(right) + 1
    scores = [[0] * columns for _ in range(rows)]

    for row in range(1, rows):
        scores[row][0] = row * GAP_SCORE
    for column in range(1, columns):
        scores[0][column] = column * GAP_SCORE

    for row in range(1, rows):
        for column in range(1, columns):
            diagonal = scores[row - 1][column - 1] + score_pair(
                left[row - 1], right[column - 1]
            )
            up = scores[row - 1][column] + GAP_SCORE
            left_score = scores[row][column - 1] + GAP_SCORE
            scores[row][column] = max(diagonal, up, left_score)

    aligned_left: list[str] = []
    aligned_right: list[str] = []
    row = len(left)
    column = len(right)

    while row > 0 or column > 0:
        if row > 0 and column > 0:
            diagonal = scores[row - 1][column - 1] + score_pair(
                left[row - 1], right[column - 1]
            )
            if scores[row][column] == diagonal:
                aligned_left.append(left[row - 1])
                aligned_right.append(right[column - 1])
                row -= 1
                column -= 1
                continue

        if row > 0 and scores[row][column] == scores[row - 1][column] + GAP_SCORE:
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


def alignment_stats(aligned_left: str, aligned_right: str) -> dict[str, int]:
    matches = 0
    mismatches = 0
    gaps = 0

    for left_char, right_char in zip(aligned_left, aligned_right):
        if left_char == "-" or right_char == "-":
            gaps += 1
        elif left_char == right_char:
            matches += 1
        else:
            mismatches += 1

    return {
        "alignment_length": len(aligned_left),
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
    }


def match_line(aligned_left: str, aligned_right: str) -> str:
    symbols: list[str] = []
    for left_char, right_char in zip(aligned_left, aligned_right):
        if left_char == "-" or right_char == "-":
            symbols.append(" ")
        elif left_char == right_char:
            symbols.append("|")
        else:
            symbols.append(".")
    return "".join(symbols)


def aligned_codons(nucleotide_sequence: str, aligned_translation: str) -> list[str]:
    codons = codons_from_dna(nucleotide_sequence)
    codon_index = 0
    output: list[str] = []

    for amino_acid in aligned_translation:
        if amino_acid == "-":
            output.append("---")
            continue
        if codon_index >= len(codons):
            output.append("???")
            continue
        output.append(codons[codon_index])
        codon_index += 1

    return output


def format_amino_acid_row(sequence: str) -> str:
    return " ".join(f" {character} " for character in sequence)


def format_marker_row(markers: str) -> str:
    return " ".join(f" {marker} " for marker in markers)


def format_alignment(
    na_header: str,
    aa_header: str,
    nucleotide_sequence: str,
    translated_sequence: str,
    protein_sequence: str,
    aligned_translation: str,
    aligned_protein: str,
    width: int,
) -> str:
    stats = alignment_stats(aligned_translation, aligned_protein)
    identity = 0.0
    if stats["alignment_length"]:
        identity = stats["matches"] / stats["alignment_length"] * 100

    lines = [
        "Translation alignment summary",
        f"  nucleotide FASTA: {na_header}",
        f"  protein FASTA: {aa_header}",
        f"  nucleotide length: {len(nucleotide_sequence)} nt",
        f"  translated length: {len(translated_sequence)} aa",
        f"  protein length: {len(protein_sequence)} aa",
        f"  alignment length: {stats['alignment_length']} aa",
        f"  matches: {stats['matches']}",
        f"  mismatches: {stats['mismatches']}",
        f"  gaps: {stats['gaps']}",
        f"  identity: {identity:.2f}%",
        "",
    ]

    codons = aligned_codons(nucleotide_sequence, aligned_translation)
    markers = match_line(aligned_translation, aligned_protein)
    for start in range(0, len(aligned_translation), width):
        end = start + width
        lines.append(f"na codon   {' '.join(codons[start:end])}")
        lines.append(f"translated {format_amino_acid_row(aligned_translation[start:end])}")
        lines.append(f"match      {format_marker_row(markers[start:end])}")
        lines.append(f"protein    {format_amino_acid_row(aligned_protein[start:end])}")
        lines.append("")

    return "\n".join(lines)


def entry() -> None:
    args = parse_args()
    na_header, nucleotide_sequence = read_single_fasta(Path(args.na))
    aa_header, protein_sequence = read_single_fasta(Path(args.aa))
    translated_sequence = translate_dna(nucleotide_sequence)
    aligned_translation, aligned_protein = hirschberg(
        translated_sequence, protein_sequence
    )

    output = format_alignment(
        na_header=na_header,
        aa_header=aa_header,
        nucleotide_sequence=nucleotide_sequence,
        translated_sequence=translated_sequence,
        protein_sequence=protein_sequence,
        aligned_translation=aligned_translation,
        aligned_protein=aligned_protein,
        width=args.width,
    )
    Path(args.out).write_text(f"{output}\n")
    print(f"Wrote alignment to {args.out}")


if __name__ == "__main__":
    entry()
