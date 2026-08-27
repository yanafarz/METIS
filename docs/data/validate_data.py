#!/usr/bin/env python3
"""Validate docs/data/*.json against the website schema.

Run from the repo root:
    python docs/data/validate_data.py

Exits 0 if the data is safe for the site to load, 1 if not.
Run this every time the data files change, especially the day the real
merged annotation table replaces the placeholders.

Standard library only. No install step.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

EXPECTED_MAJOR = 1

DATA_DIR = Path(__file__).resolve().parent

CATEGORY = {
    "detoxification", "chitin_cuticle", "digestion_insecticide_target",
    "hormone_immune", "other", "unassigned", "contaminant",
}
CONTAMINATION = {
    "clean_insect", "non_arthropod_eukaryote", "ambiguous",
    "contaminant_bacterial", "contaminant_fungal", "no_hit", "below_threshold",
}
TARGET_CLASS = {
    "lethal_rnai_target", "insecticide_binding_site",
    "resistance_breaker", "development_disruptor", "none",
}
CONFIDENCE = {"high", "medium", "low"}
ANNOTATION_SOURCE = {"blast", "interpro", "blast+interpro", "none"}
OWNER = {"chiamin", "ghaya", "syaf", "lyana"}

REQUIRED = [
    "protein_id", "annotation", "family", "category", "owner", "split",
    "hit_identity", "hit_coverage", "hit_evalue", "hit_bitscore", "domains",
    "contamination_status", "target_class", "confidence",
]

# Values that mean "missing" but are not null. These break sorting.
FAKE_NULLS = {"NA", "na", "N/A", "n/a", "-", "None", "null", "NaN", "nan"}

# Fields where the literal string "none" is a real enum value, not a missing value.
NONE_IS_VALID = {"target_class", "annotation_source"}

CORE_SHARED = [
    "annotation", "family", "subfamily", "category", "owner", "split",
    "length_aa", "target_class", "priority_score", "confidence",
    "hit_identity", "hit_coverage", "hit_bitscore",
]

errors: list[str] = []
warnings: list[str] = []


def err(msg: str) -> None:
    errors.append(msg)


def warn(msg: str) -> None:
    warnings.append(msg)


def load(name: str):
    path = DATA_DIR / name
    if not path.exists():
        err(f"{name}: file not found at {path}")
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        err(f"{name}: invalid JSON at line {exc.lineno}, column {exc.colno} — {exc.msg}")
        return None


def check_meta(name: str, doc: dict) -> None:
    meta = doc.get("meta")
    if not isinstance(meta, dict):
        err(f"{name}: missing 'meta' object")
        return
    version = meta.get("schema_version", "")
    major = version.split(".")[0] if version else ""
    if major != str(EXPECTED_MAJOR):
        err(f"{name}: schema_version is {version!r}, this validator expects {EXPECTED_MAJOR}.x")
    if meta.get("is_placeholder") is True:
        warn(f"{name}: still flagged as placeholder data")


def check_enum(where: str, field: str, value, allowed: set, nullable: bool = False) -> None:
    if value is None:
        if not nullable:
            err(f"{where}: {field} is null but must be one of {sorted(allowed)}")
        return
    if value not in allowed:
        err(f"{where}: {field}={value!r} is not a valid value. Allowed: {sorted(allowed)}")


def check_percent(where: str, field: str, value) -> None:
    if value is None:
        return
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        err(f"{where}: {field}={value!r} must be a number or null, not {type(value).__name__}")
        return
    if not 0 <= value <= 100:
        err(f"{where}: {field}={value} is outside 0-100")


def check_record(where: str, rec: dict) -> None:
    for field in REQUIRED:
        if field not in rec:
            err(f"{where}: missing required field '{field}'")

    pid = rec.get("protein_id")
    if not pid or not isinstance(pid, str):
        err(f"{where}: protein_id must be a non-empty string")

    for field, value in rec.items():
        if field in NONE_IS_VALID:
            continue
        if isinstance(value, str) and value.strip() in FAKE_NULLS:
            err(f"{where}: {field}={value!r} — use null for missing values, not a placeholder string")

    check_enum(where, "category", rec.get("category"), CATEGORY)
    check_enum(where, "owner", rec.get("owner"), OWNER)
    check_enum(where, "contamination_status", rec.get("contamination_status"), CONTAMINATION)
    check_enum(where, "target_class", rec.get("target_class"), TARGET_CLASS)
    check_enum(where, "confidence", rec.get("confidence"), CONFIDENCE)
    if "annotation_source" in rec:
        check_enum(where, "annotation_source", rec.get("annotation_source"), ANNOTATION_SOURCE)

    check_percent(where, "hit_identity", rec.get("hit_identity"))
    check_percent(where, "hit_coverage", rec.get("hit_coverage"))
    check_percent(where, "priority_score", rec.get("priority_score"))

    split = rec.get("split")
    if split is not None and split not in (1, 2, 3, 4):
        err(f"{where}: split={split!r} must be 1-4")

    domains = rec.get("domains")
    if not isinstance(domains, list):
        err(f"{where}: domains must be a list (use [] when there are none, never null)")
    else:
        for i, dom in enumerate(domains):
            if not isinstance(dom, dict) or not {"source", "id", "name"} <= set(dom):
                err(f"{where}: domains[{i}] needs keys source, id, name")

    flags = rec.get("flags", [])
    if not isinstance(flags, list):
        err(f"{where}: flags must be a list")

    rank, in_top10 = rec.get("rank"), rec.get("in_top10")
    if in_top10 is not None:
        if in_top10 and rank is None:
            err(f"{where}: in_top10 is true but rank is null")
        if not in_top10 and rank is not None:
            err(f"{where}: rank={rank} is set but in_top10 is false")

    ev = rec.get("hit_evalue")
    if ev is not None and (not isinstance(ev, (int, float)) or isinstance(ev, bool)):
        err(f"{where}: hit_evalue={ev!r} must be a number or null, not a string")

    if rec.get("contamination_status") == "no_hit" and rec.get("hit_identity") is not None:
        err(f"{where}: contamination_status is no_hit but hit_identity is set")


def main() -> int:
    proteins_doc = load("proteins.json")
    candidates_doc = load("candidates.json")
    qc_doc = load("qc_summary.json")

    if proteins_doc is None or candidates_doc is None or qc_doc is None:
        report()
        return 1

    for name, doc in (("proteins.json", proteins_doc),
                      ("candidates.json", candidates_doc),
                      ("qc_summary.json", qc_doc)):
        check_meta(name, doc)

    proteins = proteins_doc.get("proteins")
    candidates = candidates_doc.get("candidates")

    if not isinstance(proteins, list) or not proteins:
        err("proteins.json: 'proteins' must be a non-empty list")
        report()
        return 1
    if not isinstance(candidates, list) or not candidates:
        err("candidates.json: 'candidates' must be a non-empty list")
        report()
        return 1

    seen: dict[str, int] = {}
    for i, rec in enumerate(proteins):
        check_record(f"proteins[{i}]", rec)
        pid = rec.get("protein_id")
        if isinstance(pid, str):
            if pid in seen:
                err(f"proteins.json: duplicate protein_id {pid!r} at index {i} and {seen[pid]}")
            seen[pid] = i

    by_id = {r.get("protein_id"): r for r in proteins}

    ranks: list = []
    for i, cand in enumerate(candidates):
        where = f"candidates[{i}]"
        check_record(where, cand)

        pid = cand.get("protein_id")
        ranks.append(cand.get("rank"))

        partner = by_id.get(pid)
        if partner is None:
            err(f"{where}: protein_id {pid!r} is not in proteins.json — the detail modal will 404")
        else:
            for field in CORE_SHARED:
                if field in cand and field in partner and cand[field] != partner[field]:
                    err(f"{where}: {field} disagrees with proteins.json "
                        f"({cand[field]!r} vs {partner[field]!r})")
            if partner.get("rank") != cand.get("rank"):
                err(f"{where}: rank {cand.get('rank')!r} disagrees with "
                    f"proteins.json rank {partner.get('rank')!r}")

        for field in ("why_selected", "mechanism", "delivery_strategy", "next_step"):
            if not cand.get(field):
                warn(f"{where}: {field} is empty — this is what judges actually read")

        breakdown = cand.get("score_breakdown")
        score = cand.get("priority_score")
        if isinstance(breakdown, dict) and isinstance(score, (int, float)):
            total = round(sum(v for v in breakdown.values() if isinstance(v, (int, float))), 1)
            if abs(total - score) > 0.05:
                err(f"{where}: score_breakdown sums to {total} but priority_score is {score}")
            for key, value in breakdown.items():
                if isinstance(value, (int, float)) and not 0 <= value <= 20:
                    err(f"{where}: score_breakdown.{key}={value} is outside 0-20")

        spec = cand.get("specificity")
        if isinstance(spec, dict):
            status = spec.get("off_target_status")
            if status not in {"pending", "pass", "fail", None}:
                err(f"{where}: specificity.off_target_status={status!r} must be pending, pass or fail")
            if status == "pending":
                warn(f"{where}: off-target screen still pending")

    clean = [r for r in ranks if isinstance(r, int)]
    if sorted(clean) != list(range(1, len(clean) + 1)):
        err(f"candidates.json: ranks must be 1..{len(candidates)} with no gaps or repeats, got {sorted(clean)}")

    scored = [(c.get("rank"), c.get("priority_score")) for c in candidates
              if isinstance(c.get("rank"), int) and isinstance(c.get("priority_score"), (int, float))]
    scored.sort()
    for (r1, s1), (r2, s2) in zip(scored, scored[1:]):
        if s2 > s1:
            err(f"candidates.json: rank {r2} scores {s2} but rank {r1} only scores {s1} — ranking is inconsistent")

    for name, key in (("charts", dict), ("splits", list), ("proteome", dict)):
        if name not in qc_doc:
            err(f"qc_summary.json: missing '{name}'")

    for chart_key, chart in (qc_doc.get("charts") or {}).items():
        if "data" not in chart:
            err(f"qc_summary.json: charts.{chart_key} has no 'data'")
        if not chart.get("title"):
            warn(f"qc_summary.json: charts.{chart_key} has no title")

    report()
    print(f"\nChecked {len(proteins)} proteins and {len(candidates)} candidates.")
    return 1 if errors else 0


def report() -> None:
    if warnings:
        print(f"WARNINGS ({len(warnings)}) — not blocking:")
        for w in warnings:
            print(f"  ~ {w}")
        print()
    if errors:
        print(f"ERRORS ({len(errors)}) — fix before loading in the site:")
        for e in errors:
            print(f"  x {e}")
    else:
        print("PASS — data matches schema v1.x.")


if __name__ == "__main__":
    sys.exit(main())
