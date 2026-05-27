from __future__ import annotations

import argparse
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


@dataclass
class GenBankFeature:
    key: str
    location: str
    qualifiers: dict[str, list[str]] = field(default_factory=dict)


@dataclass
class GenBankRecord:
    locus: str = ""
    version: str = ""
    definition: str = ""
    accessions: list[str] = field(default_factory=list)
    organism: str = ""
    features: list[GenBankFeature] = field(default_factory=list)
    sequence: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse a GenBank file using only the Python standard library and "
            "extract a selected CDS to nucleotide/protein FASTA files."
        )
    )
    parser.add_argument(
        "genbank",
        nargs="?",
        default="KC248198.gb",
        help="Input GenBank file. Default: KC248198.gb",
    )
    parser.add_argument(
        "--cds-index",
        type=int,
        default=1,
        help="1-based CDS index to extract. Default: 1.",
    )
    parser.add_argument("--na-out", default="na.fasta", help="Nucleotide FASTA output.")
    parser.add_argument("--aa-out", default="aa.fasta", help="Protein FASTA output.")
    return parser.parse_args()


def parse_genbank(path: Path) -> GenBankRecord:
    record = GenBankRecord()
    section = "header"
    current_feature: GenBankFeature | None = None
    current_qualifier: str | None = None
    sequence_parts: list[str] = []

    for line in path.read_text().splitlines():
        keyword = line[:12].strip()

        if line.startswith("//"):
            break
        if keyword == "FEATURES":
            section = "features"
            current_feature = None
            current_qualifier = None
            continue
        if keyword == "ORIGIN":
            section = "origin"
            current_feature = None
            current_qualifier = None
            continue

        if section == "header":
            parse_header_line(record, line)
        elif section == "features":
            current_feature, current_qualifier = parse_feature_line(
                record, line, current_feature, current_qualifier
            )
        elif section == "origin":
            sequence_parts.append("".join(ch for ch in line if ch.isalpha()))

    record.sequence = "".join(sequence_parts).upper().replace("U", "T")
    return record


def parse_header_line(record: GenBankRecord, line: str) -> None:
    keyword = line[:12].strip()
    value = line[12:].strip()

    if keyword == "LOCUS":
        record.locus = value.split()[0] if value else ""
    elif keyword == "DEFINITION":
        record.definition = value
    elif keyword == "" and record.definition and not record.accessions:
        record.definition = f"{record.definition} {value}".strip()
    elif keyword == "ACCESSION":
        record.accessions = value.split()
    elif keyword == "VERSION":
        record.version = value.split()[0] if value else ""
    elif keyword == "ORGANISM":
        record.organism = value


def parse_feature_line(
    record: GenBankRecord,
    line: str,
    current_feature: GenBankFeature | None,
    current_qualifier: str | None,
) -> tuple[GenBankFeature | None, str | None]:
    if not line.startswith("     "):
        return current_feature, current_qualifier

    feature_key = line[5:21].strip()
    feature_text = line[21:].strip()

    if feature_key:
        current_feature = GenBankFeature(feature_key, feature_text)
        record.features.append(current_feature)
        return current_feature, None

    if current_feature is None:
        return current_feature, current_qualifier

    if feature_text.startswith("/"):
        key, value = parse_qualifier(feature_text)
        current_feature.qualifiers.setdefault(key, []).append(value)
        return current_feature, key

    if current_qualifier:
        values = current_feature.qualifiers[current_qualifier]
        values[-1] = append_qualifier_text(
            values[-1],
            feature_text,
            join_with_space=current_qualifier != "translation",
        )
    else:
        current_feature.location += feature_text

    return current_feature, current_qualifier


def parse_qualifier(text: str) -> tuple[str, str]:
    key_value = text[1:]
    if "=" not in key_value:
        return key_value, ""

    key, value = key_value.split("=", 1)
    return key, clean_qualifier_value(value)


def clean_qualifier_value(value: str) -> str:
    value = value.strip()
    if value.startswith('"'):
        value = value[1:]
    if value.endswith('"'):
        value = value[:-1]
    return value


def append_qualifier_text(existing: str, text: str, join_with_space: bool) -> str:
    text = clean_qualifier_value(text)
    if not text:
        return existing
    if not join_with_space:
        return f"{existing}{text}"
    if existing and not existing.endswith(" ") and not text.startswith(" "):
        return f"{existing} {text}"
    return f"{existing}{text}"


def selected_cds_feature(record: GenBankRecord, cds_index: int) -> GenBankFeature:
    if cds_index < 1:
        raise ValueError("--cds-index must be 1 or greater")

    cds_features = [feature for feature in record.features if feature.key == "CDS"]
    if len(cds_features) < cds_index:
        raise ValueError(
            f"Could not find CDS index {cds_index}; this record has "
            f"{len(cds_features)} CDS feature(s)"
        )
    return cds_features[cds_index - 1]


def feature_range(feature: GenBankFeature) -> tuple[int, int, bool]:
    location = feature.location.strip()
    is_complement = location.startswith("complement(") and location.endswith(")")
    if is_complement:
        location = location[len("complement(") : -1]

    start, end = parse_location_bounds(location)
    return start, end, is_complement


def parse_location_bounds(location: str) -> tuple[int, int]:
    if location.startswith("join(") and location.endswith(")"):
        parts = location[len("join(") : -1].split(",")
        bounds = [parse_location_bounds(part.strip()) for part in parts]
        return bounds[0][0], bounds[-1][1]

    if ".." in location:
        start_text, end_text = location.split("..", 1)
    else:
        start_text = end_text = location

    return parse_location_position(start_text), parse_location_position(end_text)


def parse_location_position(text: str) -> int:
    digits = "".join(ch for ch in text if ch.isdigit())
    if not digits:
        raise ValueError(f"Could not parse feature location position: {text}")
    return int(digits)


def coding_sequence(record: GenBankRecord, feature: GenBankFeature) -> str:
    location = feature.location.strip()
    is_complement = location.startswith("complement(") and location.endswith(")")
    if is_complement:
        location = location[len("complement(") : -1]

    sequence = "".join(extract_location_pieces(record.sequence, location))
    if is_complement:
        sequence = reverse_complement(sequence)
    return sequence.upper().replace("U", "T")


def extract_location_pieces(sequence: str, location: str) -> Iterable[str]:
    if location.startswith("join(") and location.endswith(")"):
        for part in location[len("join(") : -1].split(","):
            yield from extract_location_pieces(sequence, part.strip())
        return

    start, end = parse_location_bounds(location)
    if end > len(sequence):
        raise ValueError(f"Range end {end} is outside sequence length {len(sequence)}")
    yield sequence[start - 1 : end]


def reverse_complement(sequence: str) -> str:
    complement = str.maketrans(
        "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
        "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
    )
    return sequence.translate(complement)[::-1]


def qualifier_text(feature: GenBankFeature, key: str, default: str = "-") -> str:
    values = feature.qualifiers.get(key)
    if not values:
        return default
    return values[0]


def write_fasta(path: Path, header: str, sequence: str) -> None:
    wrapped = "\n".join(textwrap.wrap(sequence, 70))
    path.write_text(f">{header}\n{wrapped}\n")


def print_summary(
    genbank_path: Path,
    record: GenBankRecord,
    cds: GenBankFeature,
    cds_index: int,
    range_label: str,
    nucleotide_text: str,
    protein_sequence: str,
    na_out: str,
    aa_out: str,
) -> None:
    accession = ",".join(record.accessions) or "-"
    print("Pure-Python GenBank extraction summary")
    print(f"  input: {genbank_path}")
    print(f"  accession/version: {accession} / {record.version or '-'}")
    print(f"  definition: {record.definition or '-'}")
    print(f"  organism: {record.organism or '-'}")
    print(f"  full sequence length: {len(record.sequence)} nt")
    print(f"  selected nucleotide range from CDS: {range_label}")
    print(f"  selected nucleotide length: {len(nucleotide_text)} nt")
    print(f"  protein source: CDS {cds_index} /translation")
    print(f"  protein length written: {len(protein_sequence)} aa")
    print(f"  nucleotide output: {na_out}")
    print(f"  protein output: {aa_out}")
    print(f"  selected CDS location: {cds.location}")
    print(f"  selected CDS codon_start: {qualifier_text(cds, 'codon_start')}")
    print(f"  selected CDS product: {qualifier_text(cds, 'product')}")
    print(f"  selected CDS protein_id: {qualifier_text(cds, 'protein_id')}")


def entry() -> None:
    args = parse_args()
    genbank_path = Path(args.genbank)
    record = parse_genbank(genbank_path)
    cds = selected_cds_feature(record, args.cds_index)
    start, end, is_complement = feature_range(cds)

    nucleotide_text = coding_sequence(record, cds)
    translations = cds.qualifiers.get("translation")
    if not translations:
        raise ValueError(f"CDS index {args.cds_index} has no /translation qualifier")
    protein_sequence = translations[0]

    range_label = f"{start}..{end}"
    if is_complement:
        range_label = f"complement({range_label})"

    fasta_id = record.version or record.locus or genbank_path.stem
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
    entry()
