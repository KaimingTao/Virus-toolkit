from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract a nucleotide FASTA from a GenBank record and write the "
            "annotated CDS /translation to a protein FASTA."
        )
    )
    parser.add_argument(
        "genbank",
        nargs="?",
        default="sequence.gb",
        help="Input GenBank file. Default: sequence.gb",
    )
    parser.add_argument(
        "--cds-index",
        type=int,
        default=1,
        help="1-based CDS index to use for nucleotide range and /translation. Default: 1.",
    )
    parser.add_argument("--na-out", default="na.fasta", help="Nucleotide FASTA output.")
    parser.add_argument("--aa-out", default="aa.fasta", help="Protein FASTA output.")
    return parser.parse_args()


def feature_range(feature: SeqFeature) -> tuple[int, int, bool]:
    start = int(feature.location.start) + 1
    end = int(feature.location.end)
    is_complement = feature.location.strand == -1
    return start, end, is_complement


def qualifier_text(feature: SeqFeature, key: str, default: str = "-") -> str:
    values = feature.qualifiers.get(key)
    if not values:
        return default
    return str(values[0])


def write_fasta(path: Path, header: str, sequence: str) -> None:
    wrapped = "\n".join(textwrap.wrap(sequence, 70))
    path.write_text(f">{header}\n{wrapped}\n")


def selected_cds_feature(record: SeqRecord, cds_index: int) -> SeqFeature:
    if cds_index < 1:
        raise ValueError("--cds-index must be 1 or greater")

    cds_features = [feature for feature in record.features if feature.type == "CDS"]
    if len(cds_features) < cds_index:
        raise ValueError(
            f"Could not find CDS index {cds_index}; this record has "
            f"{len(cds_features)} CDS feature(s)"
        )
    return cds_features[cds_index - 1]


def print_summary(
    genbank_path: Path,
    record: SeqRecord,
    cds: SeqFeature,
    cds_index: int,
    range_label: str,
    nucleotide_text: str,
    protein_sequence: str,
    na_out: str,
    aa_out: str,
) -> None:
    organism = record.annotations.get("organism", "-")
    accession = ",".join(record.annotations.get("accessions", [])) or "-"
    print("GenBank extraction summary")
    print(f"  input: {genbank_path}")
    print(f"  accession/version: {accession} / {record.id or '-'}")
    print(f"  definition: {record.description or '-'}")
    print(f"  organism: {organism}")
    print(f"  full sequence length: {len(record.seq)} nt")
    print(f"  selected nucleotide range from CDS: {range_label}")
    print(f"  selected nucleotide length: {len(nucleotide_text)} nt")
    print(f"  protein source: CDS {cds_index} /translation")
    print(f"  protein length written: {len(protein_sequence)} aa")
    print(f"  stop codons in /translation: {protein_sequence.count('*')}")
    print(f"  ambiguous amino acids X in /translation: {protein_sequence.count('X')}")
    print(f"  nucleotide output: {na_out}")
    print(f"  protein output: {aa_out}")
    print(f"  selected CDS location: {cds.location}")
    print(f"  selected CDS codon_start: {qualifier_text(cds, 'codon_start')}")
    print(f"  selected CDS product: {qualifier_text(cds, 'product')}")
    print(f"  selected CDS protein_id: {qualifier_text(cds, 'protein_id')}")


def main() -> None:
    args = parse_args()
    genbank_path = Path(args.genbank)
    record = SeqIO.read(genbank_path, "genbank")
    cds = selected_cds_feature(record, args.cds_index)
    start, end, is_complement = feature_range(cds)

    if end > len(record.seq):
        raise ValueError(f"Range end {end} is outside sequence length {len(record.seq)}")

    nucleotide_sequence = record.seq[start - 1 : end]
    if is_complement:
        nucleotide_sequence = nucleotide_sequence.reverse_complement()
    nucleotide_text = str(nucleotide_sequence).upper().replace("U", "T")

    translations = cds.qualifiers.get("translation")
    if not translations:
        raise ValueError(f"CDS index {args.cds_index} has no /translation qualifier")
    protein_sequence = str(translations[0])

    range_label = f"{start}..{end}"
    if is_complement:
        range_label = f"complement({range_label})"

    fasta_id = record.id or record.name or genbank_path.stem
    na_header = f"{fasta_id}|range={range_label}|length={len(nucleotide_text)}"
    aa_header = (
        f"{fasta_id}|cds_index={args.cds_index}|source=/translation"
        f"|length={len(protein_sequence)}"
    )

    write_fasta(Path(args.na_out), na_header, nucleotide_text)
    write_fasta(Path(args.aa_out), aa_header, protein_sequence)

    print_summary(
        genbank_path=genbank_path,
        record=record,
        cds=cds,
        cds_index=args.cds_index,
        range_label=range_label,
        nucleotide_text=nucleotide_text,
        protein_sequence=protein_sequence,
        na_out=args.na_out,
        aa_out=args.aa_out,
    )


if __name__ == "__main__":
    main()
