#!/usr/bin/env python3
"""Build docs/data/proteins.json from the team's screening TSVs.

Run from the repo root:
    python docs\\data\\build_from_screening.py

Reads every results/mp*_screening.tsv it can find and writes real BLAST
evidence into the website's schema. Also prints a per-split diagnostic report,
because the splits currently disagree about which reference database was used
and that has to be settled before anyone merges them.

What this fills in:  protein_id, annotation, hit_* fields, contamination_status,
                     confidence, flags.
What stays empty:    family, subfamily, domains  (InterProScan pending)
                     target_class, mechanism, priority_score, rank  (ranking pending)
                     category  (biology classification pending)
                     length_aa  (not in the screening file)

Nothing is invented. Every field is either read from the TSV or left null.

Standard library only.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

# --------------------------------------------------------------------- config

REPO = Path(__file__).resolve().parents[2]
RESULTS = REPO / "results"
OUT = REPO / "docs" / "data" / "proteins.json"
SAMPLE_BACKUP = REPO / "docs" / "data" / "proteins.sample.json"

# EDIT THIS once you know who ran which split. It drives the `owner` filter.
# Leave a split out and it becomes null rather than a wrong guess.
SPLIT_OWNERS: dict[int, str | None] = {
    1: None,        # not pushed to results/ yet
    2: None,        # unknown — git author is whoever committed, not whoever ran it
    3: "chiamin",   # confirmed: 6,488 rows matches my own run exactly
    4: None,        # unknown
}

# Assumed proteins per split, for hit-rate percentages only.
SPLIT_TOTALS = {1: 6622, 2: 6622, 3: 6622, 4: 6624}

# The screening script's classification vocabulary -> the website's enum.
CLASSIFICATION_MAP = {
    "KEEP_ARTHROPOD": "clean_insect",
    "NON_ARTHROPOD_EUKARYOTE": "non_arthropod_eukaryote",
    "AMBIGUOUS": "ambiguous",
    "POSSIBLE_BACTERIAL": "contaminant_bacterial",
    "POSSIBLE_FUNGAL": "contaminant_fungal",
    "NO_QUALIFYING_HIT": "no_qualifying_hit",
    "NO_HIT": "no_hit",
    "OTHER_PROKARYOTE": "contaminant_bacterial",
}

BLANKS = {"", "NA", "N/A", "-", "None", "none", "null", "nan", "NaN"}

# best_hit_description holds the raw UniProt FASTA header, e.g.
#   sp|P37469|DNAC_BACSU Replicative DNA helicase DnaC OS=Bacillus subtilis (strain 168) OX=224308 GN=dnaC PE=1 SV=2
# Judges should see "Replicative DNA helicase DnaC", not that whole line. Splitting
# it also recovers gene_symbol from GN=, which the screening file has no column for.
ACC_PREFIX = re.compile(r"^(?:sp|tr)\|[^|]+\|\S*\s+")
META_START = re.compile(r"\s+(?:OS|OX|GN|PE|SV)=")


def parse_uniprot_header(description: str | None) -> dict:
    """Split a UniProt FASTA header into protein name, gene symbol and organism."""
    if not description:
        return {}
    raw = description.strip()
    body = ACC_PREFIX.sub("", raw)
    out = {"hit_header": raw, "header_had_accession": body != raw}

    match = META_START.search(body)
    name = (body[: match.start()] if match else body).strip()
    out["annotation"] = name or None

    for key, field in (("OS", "organism"), ("GN", "gene_symbol")):
        found = re.search(rf"{key}=(.+?)(?=\s+(?:OS|OX|GN|PE|SV)=|$)", body)
        if found:
            out[field] = found.group(1).strip()
    return out

# ------------------------------------------------------------------ utilities

def clean(value: str | None) -> str | None:
    if value is None:
        return None
    value = value.strip()
    return None if value in BLANKS else value


def to_float(value: str | None) -> float | None:
    value = clean(value)
    if value is None:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def to_int(value: str | None) -> int | None:
    value = clean(value)
    if value is None:
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def accession_source(accession: str | None) -> str | None:
    """sp| = Swiss-Prot (curated), tr| = TrEMBL (unreviewed)."""
    if not accession:
        return None
    if accession.startswith("sp|"):
        return "swissprot"
    if accession.startswith("tr|"):
        return "trembl"
    return "unparsed"


def confidence_from_blast(identity: float | None, coverage: float | None) -> str:
    """Annotation confidence derived from alignment quality alone.

    This is the one derived field in the file. It is not a biological
    judgement about the gene, only about how safely we can copy the hit's
    label onto our protein. Stated explicitly in meta so it can be defended.
    """
    if identity is None or coverage is None:
        return "low"
    if identity >= 70 and coverage >= 80:
        return "high"
    if identity >= 40 and coverage >= 50:
        return "medium"
    return "low"


# -------------------------------------------------------------------- reading

def split_number(path: Path) -> int | None:
    stem = path.name
    for n in (1, 2, 3, 4):
        if stem.startswith(f"mp{n}_"):
            return n
    return None


def read_screening(path: Path, split: int, stats: dict) -> list[dict]:
    records = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        missing = {"gene_id", "classification"} - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"{path.name}: missing expected column(s) {sorted(missing)}\n"
                             f"  columns found: {reader.fieldnames}")

        for row in reader:
            gene_id = clean(row.get("gene_id"))
            if not gene_id:
                continue

            raw_class = (clean(row.get("classification")) or "").upper()
            status = CLASSIFICATION_MAP.get(raw_class)
            flags: list[str] = []
            if status is None:
                status = "ambiguous"
                flags.append("unmapped_classification")
                stats["unmapped"][raw_class] += 1

            accession = clean(row.get("best_hit_accession"))
            source = accession_source(accession)
            if source:
                flags.append(f"hit_{source}")
                stats["accession_source"][source] += 1
            else:
                stats["accession_source"]["none"] += 1

            identity = to_float(row.get("best_identity_pct"))
            coverage = to_float(row.get("best_query_coverage_pct"))
            evalue = to_float(row.get("best_evalue"))
            bitscore = to_float(row.get("best_bitscore"))

            if identity is not None and identity < 40:
                flags.append("low_identity")
            if coverage is not None and coverage < 50:
                flags.append("low_coverage")
            if status == "no_qualifying_hit":
                flags.append("no_qualifying_hit")

            description = clean(row.get("best_hit_description"))
            parsed = parse_uniprot_header(description)
            if parsed.get("header_had_accession") is not None:
                stats["header_format"][
                    "with_accession_prefix" if parsed["header_had_accession"] else "name_only"
                ] += 1

            species = clean(row.get("best_hit_species")) or parsed.get("organism")
            if species:
                stats["species"][species] += 1
            stats["classification"][status] += 1

            records.append({
                "protein_id": gene_id,
                "length_aa": None,
                "annotation": parsed.get("annotation") or "No qualifying hit",
                "gene_symbol": parsed.get("gene_symbol"),
                "family": None,
                "subfamily": None,
                "category": "unassigned",
                "owner": SPLIT_OWNERS.get(split),
                "split": split,
                "hit_subject_id": accession,
                "hit_subject_name": parsed.get("annotation"),
                "hit_header": parsed.get("hit_header"),
                "hit_organism": species,
                "hit_taxid": to_int(row.get("best_hit_taxid")),
                "hit_identity": identity,
                "hit_coverage": coverage,
                "hit_evalue": evalue,
                "hit_bitscore": bitscore,
                "domains": [],
                "annotation_source": "blast" if accession else "none",
                "contamination_status": status,
                "target_class": None,
                "mechanism": None,
                "priority_score": None,
                "confidence": confidence_from_blast(identity, coverage),
                "rank": None,
                "in_top10": False,
                "flags": flags,
                "notes": "",
                "total_blast_hits": to_int(row.get("total_blast_hits")),
                "qualifying_hits": to_int(row.get("qualifying_hits")),
                "hit_lineage": clean(row.get("best_hit_lineage")),
                "screen_decision": clean(row.get("decision")),
                "screen_reason": clean(row.get("reason")),
            })
    return records


# ------------------------------------------------------------------ reporting

def report(per_split: dict) -> bool:
    """Prints the diagnostic. Returns True if the splits look comparable."""
    print("\n" + "=" * 78)
    print("PER-SPLIT DIAGNOSTIC")
    print("=" * 78)

    shares = {}

    for split in sorted(per_split):
        stats = per_split[split]
        src = stats["accession_source"]
        sp, tr, unparsed = src["swissprot"], src["trembl"], src["unparsed"]
        qualifying = sp + tr + unparsed
        total = stats["rows"]
        owner = SPLIT_OWNERS.get(split) or "UNKNOWN"
        denom = SPLIT_TOTALS.get(split, total)

        print(f"\n  split {split}  ({owner})   {total:,} rows"
              f"   hit rate {total / denom * 100:.1f}% of {denom:,}")

        if qualifying:
            share = sp / qualifying * 100
            shares[split] = share
            print(f"    accessions   Swiss-Prot {sp:,} ({share:.1f}%)"
                  f"   TrEMBL {tr:,}   unparsed {unparsed:,}   blank {src['none']:,}")

        print("    screening    ", end="")
        print(",  ".join(f"{k} {v:,}" for k, v in stats["classification"].most_common()))

        top = stats["species"].most_common(5)
        if top:
            print("    top species  ", end="")
            print(",  ".join(f"{name} ({count:,})" for name, count in top))

        hf = stats["header_format"]
        if hf:
            print("    hit headers  ", end="")
            print(",  ".join(f"{k} {v:,}" for k, v in hf.most_common()))

        if stats["unmapped"]:
            print(f"    UNMAPPED classification values: {dict(stats['unmapped'])}")

    comparable = True
    if len(shares) > 1:
        spread = max(shares.values()) - min(shares.values())
        print("\n" + "-" * 78)
        if spread > 20:
            comparable = False
            print("  WARNING: the splits do not agree on reference database.")
            print(f"  Swiss-Prot share of qualifying hits ranges over {spread:.0f} percentage points:")
            for split in sorted(shares):
                print(f"      split {split}: {shares[split]:.1f}%")
            print("""
  The four splits are random slices of ONE proteome, so this ratio should be
  near-identical across them. A spread this wide means the splits were searched
  against different databases, or against differently-built versions of it.

  Why it matters: identity scores, e-values and taxonomic calls are only
  comparable within a single database. Ranking genes from different splits on
  one scale would silently favour whichever split had the deeper database.

  Compare the 'top species' lines above. A split dominated by Drosophila,
  human or mouse was searched against curated Swiss-Prot. A split dominated by
  Spodoptera, Bombyx, Plutella or Helicoverpa reached into TrEMBL, where the
  sequenced moth genomes live.

  Resolve before merging. Ask each person for the exact -db argument they used.""")
        else:
            print(f"  Swiss-Prot share is consistent across splits "
                  f"(spread {spread:.0f} points). Splits look comparable.")
    return comparable


# ----------------------------------------------------------------------- main

def main() -> int:
    if not RESULTS.is_dir():
        print(f"No results/ directory at {RESULTS}", file=sys.stderr)
        print("Run this from the repo root: python docs\\data\\build_from_screening.py", file=sys.stderr)
        return 1

    files = sorted(RESULTS.glob("mp*_screening.tsv"))
    if not files:
        print(f"No mp*_screening.tsv files in {RESULTS}", file=sys.stderr)
        return 1

    all_records: list[dict] = []
    per_split: dict[int, dict] = {}

    for path in files:
        split = split_number(path)
        if split is None:
            print(f"  skipping {path.name} — cannot tell which split it is")
            continue

        stats = {
            "rows": 0,
            "accession_source": Counter({"swissprot": 0, "trembl": 0, "unparsed": 0, "none": 0}),
            "classification": Counter(),
            "species": Counter(),
            "unmapped": defaultdict(int),
            "header_format": Counter(),
        }
        records = read_screening(path, split, stats)
        stats["rows"] = len(records)
        per_split[split] = stats
        all_records.extend(records)
        print(f"  read {path.name:24} split {split}   {len(records):,} proteins")

    # protein_id must be unique across the whole proteome
    seen: dict[str, int] = {}
    duplicates = []
    for rec in all_records:
        pid = rec["protein_id"]
        if pid in seen:
            duplicates.append((pid, seen[pid], rec["split"]))
        seen[pid] = rec["split"]
    if duplicates:
        print(f"\n  WARNING: {len(duplicates)} protein_id(s) appear in more than one split.")
        for pid, a, b in duplicates[:5]:
            print(f"      {pid}: splits {a} and {b}")
        print("      The chunks are supposed to be disjoint — worth telling the team.")

    comparable = report(per_split)

    # Keep the placeholder file so the site can be rolled back if needed.
    if OUT.exists() and not SAMPLE_BACKUP.exists():
        payload = json.loads(OUT.read_text(encoding="utf-8"))
        if payload.get("meta", {}).get("is_placeholder"):
            SAMPLE_BACKUP.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
            print(f"\n  kept the placeholder file as {SAMPLE_BACKUP.name}")

    unknown_owners = sorted(s for s in per_split if not SPLIT_OWNERS.get(s))

    doc = {
        "meta": {
            "schema_version": "1.2.0",
            "dataset": "proteins",
            "generated": date.today().isoformat(),
            "is_placeholder": False,
            "data_status": "REAL BLAST evidence from the team's screening files. "
                           "Annotation, family and ranking fields are still empty.",
            "record_count": len(all_records),
            "source_files": [p.name for p in files],
            "splits_included": sorted(per_split),
            "generated_by": "docs/data/build_from_screening.py",
            "empty_by_design": {
                "family, subfamily, domains": "waiting on InterProScan",
                "target_class, mechanism, priority_score, rank": "waiting on scoring and ranking",
                "category": "waiting on biology classification",
                "length_aa": "not present in the screening file",
            },
            "derived_fields": {
                "confidence": "From alignment quality only: high = identity >=70% and coverage >=80%; "
                              "medium = >=40% and >=50%; low = below that or missing. Describes how "
                              "safely the hit's label transfers, not the gene's importance.",
            },
            "open_questions": (
                ([] if comparable else
                 ["Splits disagree on reference database — see the Swiss-Prot share per split. "
                  "Resolve before ranking genes across splits."]) +
                ([f"Owner unknown for split(s) {unknown_owners} — set SPLIT_OWNERS in "
                  "build_from_screening.py" ] if unknown_owners else [])
            ),
        },
        "proteins": all_records,
    }

    OUT.write_text(json.dumps(doc, indent=2, ensure_ascii=False), encoding="utf-8")

    print("\n" + "=" * 78)
    print(f"  wrote {OUT.relative_to(REPO)}  —  {len(all_records):,} real proteins")
    if not comparable:
        print("  NOTE: written anyway so you can keep building, but the database")
        print("        inconsistency above is flagged inside the file's meta block.")
    print("  next: python docs\\data\\validate_data.py")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
