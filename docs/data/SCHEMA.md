# Website data schema v1.0.0

The website reads three JSON files from `docs/data/`. Nothing else. No 78 GB database, no BLAST output, no FASTA.

| File | What it holds | Used by |
|---|---|---|
| `candidates.json` | The Top-10 target genes with full justification | Top-10 evidence cards |
| `proteins.json` | Every annotated protein worth showing | Protein Explorer table, search, filters |
| `qc_summary.json` | Counts and chart-ready arrays | QC / contamination charts, headline stats |

All three currently contain **placeholder data**. Each has `meta.is_placeholder: true`. The site should render a visible "sample data" banner while that flag is true, so nobody demos fake results by accident.

---

## 1. Design rules

**Records are flat.** One JSON key per column. This is deliberate: whoever produces the merged annotation table will produce a CSV, and `df.to_json(orient="records")` gives exactly this shape. The only nested field is `domains`, which is genuinely a list.

**Missing means `null`, never `"NA"`, `"-"`, `""` or `"None"`.** The JS checks `value == null` in one place. Strings like `"NA"` will sort as text and break the table.

**Numbers are numbers.** `hit_identity: 78.4`, not `"78.4%"`. Formatting is the display layer's job.

**`protein_id` is the primary key** and must be identical across all files. Every entry in `candidates.json` must also exist in `proteins.json` with the same `protein_id`.

**`candidates.json` records are a superset of `proteins.json` records.** Same core field names, plus ranking and justification fields. This means one JS function renders the detail modal for both.

---

## 2. Core fields (both files)

| Field | Type | Notes |
|---|---|---|
| `protein_id` | string | e.g. `g10421.t1`. Primary key. Required, never null. |
| `length_aa` | int | Amino acid length. |
| `annotation` | string | Human-readable function. This is what search matches on. |
| `gene_symbol` | string \| null | e.g. `CYP6AE`. Null if unnamed. |
| `family` | string \| null | e.g. `Cytochrome P450`. |
| `subfamily` | string \| null | e.g. `CYP6AE`. Fills in after the classification step. |
| `category` | enum | See enums below. |
| `owner` | enum | `chiamin` \| `ghaya` \| `syaf` \| `lyana`. Who ran that split. |
| `split` | int | 1–4. |
| `hit_subject_id` | string \| null | UniProt accession of the best hit, e.g. `tr|A0A...|NAME`. |
| `hit_subject_name` | string \| null | Description of the best hit. |
| `hit_organism` | string \| null | Species of the best hit. |
| `hit_taxid` | int \| null | NCBI taxonomy ID. |
| `hit_identity` | float \| null | Percent, 0–100. |
| `hit_coverage` | float \| null | Percent, 0–100. Query coverage. |
| `hit_evalue` | float \| null | `0.0` legitimately occurs (below float precision) — display as `< 1e-180`. |
| `hit_bitscore` | float \| null | |
| `domains` | array | `[{ "source": "Pfam", "id": "PF00067", "name": "Cytochrome P450" }]`. Empty array if none. Never null. |
| `annotation_source` | enum | `blast` \| `interpro` \| `blast+interpro` \| `none`. |
| `contamination_status` | enum | See enums below. |
| `target_class` | enum | See enums below. |
| `mechanism` | string \| null | One or two sentences: how attacking this gene hurts the pest. |
| `priority_score` | float \| null | 0–100. Null if not scored. |
| `confidence` | enum | `high` \| `medium` \| `low`. |
| `rank` | int \| null | 1–10 for Top-10 members, null otherwise. |
| `in_top10` | bool | Must agree with `rank`. |
| `flags` | array of strings | Free-form warning tags, e.g. `low_coverage`, `too_conserved`. Drives the badge row. |
| `notes` | string | May be empty. |

## 3. Extra fields, `candidates.json` only

| Field | Type | Notes |
|---|---|---|
| `why_selected` | string | Two or three sentences a judge can read and immediately understand the pick. |
| `delivery_strategy` | string | e.g. `Oral dsRNA`, `Small-molecule modulator`. |
| `literature_support` | string | What is already published about knocking this out. |
| `score_breakdown` | object | Five keys, each 0–20, summing to `priority_score`. Drives the score bars on the card. |
| `specificity` | object | `off_target_status` (`pending` \| `pass` \| `fail`), `off_target_method`, `beneficial_conflicts` (array), `vertebrate_orthologue` (bool), `notes`. |
| `risks` | array of strings | Honest weaknesses. Judges reward these. |
| `next_step` | string | The experiment or check that would validate it. |

## 4. Enums

```
category              detoxification | chitin_cuticle | digestion_insecticide_target |
                      hormone_immune | other | unassigned | contaminant

contamination_status  clean_insect | non_arthropod_eukaryote | ambiguous |
                      contaminant_bacterial | contaminant_fungal |
                      no_hit | no_qualifying_hit

target_class          lethal_rnai_target | insecticide_binding_site |
                      resistance_breaker | development_disruptor | none

confidence            high | medium | low
annotation_source     blast | interpro | blast+interpro | none
owner                 chiamin | ghaya | syaf | lyana
```

Enum values are lowercase snake_case in the data. Display labels live in a lookup object in `script.js`, so the wording on screen can change without touching the data.

## 5. CSV to JSON conversion contract

If the merged table arrives as a CSV, these are the only rules needed to convert it:

- One CSV column per field name above, spelled identically.
- Empty cell → `null`. Do not write `NA`.
- `domains` column: semicolon-separated `SOURCE:ID:NAME`, e.g.
  `Pfam:PF02798:GST N-terminal;Pfam:PF00043:GST C-terminal`
- `flags` and `beneficial_conflicts`: semicolon-separated plain strings.
- `in_top10`: `TRUE` / `FALSE`.
- Then: `df.to_json("proteins.json", orient="records", indent=2)` and wrap it in the `{ "meta": {...}, "proteins": [...] }` envelope.

## 6. Versioning

`meta.schema_version` uses semver. Bump the minor version when adding a field, the major version when renaming or removing one. The loader logs a warning if the file's major version does not match the one it expects.

---

## 7. Message to send Lyana now

> Hi Lyana — for the website I need one merged annotation table (BLAST + InterPro + contamination + Top-10 ranking) as a single CSV. Two questions:
>
> 1. Who is producing it, and roughly when?
> 2. What exact column names will it have?
>
> I've already defined the schema the website expects and built the site against sample data in that shape, so if your columns match, plugging in the real table is a five-minute job. If they don't match I'll write a converter — that's fine, I just need to know the column names early. Schema is committed at `docs/data/SCHEMA.md`.
>
> The fields I'm assuming: `protein_id`, `annotation`, `gene_symbol`, `family`, `subfamily`, `category`, `owner`, `split`, `hit_subject_id`, `hit_organism`, `hit_identity`, `hit_coverage`, `hit_evalue`, `hit_bitscore`, `domains`, `contamination_status`, `target_class`, `mechanism`, `priority_score`, `confidence`, `rank`.
>
> Also: is anything other than `priority_score` deciding the Top-10 ranking? I want the website to show judges *why* each gene was picked, not just the order.
