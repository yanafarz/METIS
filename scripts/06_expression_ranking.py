import os
import re
import numpy as np
import pandas as pd
from pathlib import Path


# ==========================================================
# METISA PLANA
# RNA LIFE-STAGE EXPRESSION PRIORITIZATION
# ==========================================================

print("=" * 78)
print("METISA PLANA — RNA LIFE-STAGE EXPRESSION PRIORITIZATION")
print("=" * 78)


# ==========================================================
# HELPER FUNCTIONS
# ==========================================================

def ask_existing_file(prompt):

    while True:

        path = input(prompt).strip().strip('"')

        if os.path.isfile(path):
            return path

        print("\nERROR: File not found.")
        print("Please enter a valid file path.\n")


def ask_integer(prompt, default):

    while True:

        value = input(
            f"{prompt} [{default}]: "
        ).strip()

        if value == "":
            return default

        try:

            value = int(value)

            if value > 0:
                return value

        except ValueError:
            pass

        print("Please enter a positive integer.")


def ask_float(prompt, default):

    while True:

        value = input(
            f"{prompt} [{default}]: "
        ).strip()

        if value == "":
            return default

        try:

            value = float(value)

            if 0 <= value <= 100:
                return value

        except ValueError:
            pass

        print("Please enter a number between 0 and 100.")


def ask_stage(prompt):

    valid_stages = {
        "larva": "Larva",
        "larval": "Larva",
        "third-instar": "Larva",
        "third instar": "Larva",
        "adult": "Adult",
        "pupa": "Pupa",
        "egg": "Egg",
    }

    while True:

        value = input(
            prompt
        ).strip().lower()

        if value in valid_stages:
            return valid_stages[value]

        print(
            "\nPlease enter one of:"
        )

        print(
            "Larva, Adult, Pupa, Egg\n"
        )


def clean_gene_id(value):

    value = str(value).strip()

    # Remove transcript suffix:
    # g45.t1 -> g45
    # g45.t2 -> g45
    value = re.sub(
        r"\.t\d+$",
        "",
        value,
        flags=re.IGNORECASE
    )

    return value


# ==========================================================
# INPUT FILES
# ==========================================================

print("\n" + "=" * 78)
print("INPUT FILES")
print("=" * 78)


candidate_file = ask_existing_file(
    "Top-500 candidate TSV path: "
)


print(
    "\nHow many RNA quant.sf files do you want to use?"
)

print(
    "You may provide 1 file or multiple files."
)

print(
    "At least one Larva sample is required."
)


rna_file_count = ask_integer(
    "Number of RNA files",
    1
)


rna_inputs = []


for i in range(
    1,
    rna_file_count + 1
):

    print(
        f"\n--- RNA SAMPLE {i} ---"
    )

    rna_path = ask_existing_file(
        "quant.sf path: "
    )

    stage = ask_stage(
        "Stage "
        "(Larva/Adult/Pupa/Egg): "
    )

    rna_inputs.append({
        "file": rna_path,
        "stage": stage
    })


# ==========================================================
# CHECK FOR LARVAL SAMPLE
# ==========================================================

larval_inputs = [
    x
    for x in rna_inputs
    if x["stage"] == "Larva"
]


if not larval_inputs:

    raise ValueError(
        "\nERROR: At least one Larva RNA sample "
        "is required because this pipeline prioritizes "
        "larval-stage targets."
    )


# ==========================================================
# DETERMINE RNA ANALYSIS MODE
# ==========================================================

non_larval_inputs = [
    x
    for x in rna_inputs
    if x["stage"] != "Larva"
]


if non_larval_inputs:

    RNA_MODE = "MULTI_STAGE"

    print(
        "\nRNA analysis mode: MULTI-STAGE"
    )

    print(
        "Larval expression, larval specificity, "
        "and stage dominance will be calculated."
    )

else:

    RNA_MODE = "LARVA_ONLY"

    print(
        "\nRNA analysis mode: LARVA-ONLY"
    )

    print(
        "Only larval expression will be used "
        "for RNA prioritization."
    )

    print(
        "Missing Adult/Pupa/Egg data will NOT "
        "be treated as zero expression."
    )


# ==========================================================
# OUTPUT SETTINGS
# ==========================================================

print("\n" + "=" * 78)
print("OUTPUT SETTINGS")
print("=" * 78)


output_dir = str(
    Path(candidate_file).parent
)


output_basename = input(
    "Output file name: "
).strip()


if output_basename == "":
    output_basename = "metisa_plana_rna"


top_n = ask_integer(
    "Number of candidates to retain",
    100
)


rna_weight = ask_float(
    "RNA evidence weight (%)",
    15.0
)


annotation_weight = (
    100.0 - rna_weight
)


# ==========================================================
# OUTPUT PATHS
# ==========================================================

top_file = os.path.join(
    output_dir,
    output_basename + "_top.tsv"
)


all_file = os.path.join(
    output_dir,
    output_basename + "_all500.tsv"
)


summary_file = os.path.join(
    output_dir,
    output_basename + "_summary.txt"
)


excel_file = os.path.join(
    output_dir,
    output_basename + ".xlsx"
)


# ==========================================================
# DISPLAY SETTINGS
# ==========================================================

print("\n" + "=" * 78)
print("SETTINGS")
print("=" * 78)


print(
    f"Candidate file:       {candidate_file}"
)


print(
    f"RNA files:            {len(rna_inputs)}"
)


for i, item in enumerate(
    rna_inputs,
    start=1
):

    print(
        f"  RNA {i}: "
        f"{item['stage']} — "
        f"{item['file']}"
    )


print(
    f"RNA analysis mode:    {RNA_MODE}"
)


print(
    f"Output directory:     {output_dir}"
)


print(
    f"Output basename:      {output_basename}"
)


print(
    f"Candidates retained:  {top_n}"
)


print(
    f"RNA weight:           {rna_weight:.1f}%"
)


print(
    f"Annotation weight:    "
    f"{annotation_weight:.1f}%"
)


# ==========================================================
# READ CANDIDATES
# ==========================================================

print("\n" + "=" * 78)
print("READING TOP-500 CANDIDATES")
print("=" * 78)


candidates = pd.read_csv(
    candidate_file,
    sep="\t"
)


print(
    f"Rows loaded: "
    f"{len(candidates):,}"
)


if len(candidates) > 500:

    print(
        "Input contains more than 500 rows."
    )

    print(
        "Only the first 500 candidates "
        "will be used."
    )

    candidates = (
        candidates
        .head(500)
        .copy()
    )


elif len(candidates) < 500:

    print(
        "WARNING: Input contains fewer "
        "than 500 candidates."
    )


# ==========================================================
# IDENTIFY CANDIDATE ID
# ==========================================================

gene_column = candidates.columns[0]


print(
    f"Candidate ID column: "
    f"{gene_column}"
)


if "gene_id" in candidates.columns:

    candidates["gene_id"] = (
        candidates["gene_id"]
        .apply(clean_gene_id)
    )

else:

    candidates["gene_id"] = (
        candidates[gene_column]
        .apply(clean_gene_id)
    )


# ==========================================================
# CHECK FINAL SCORE
# ==========================================================

if "final_score" not in candidates.columns:

    raise ValueError(
        "\nERROR: 'final_score' column "
        "was not found in the candidate file."
    )


# ==========================================================
# READ ALL RNA FILES
# ==========================================================

print("\n" + "=" * 78)
print("READING RNA QUANTIFICATION")
print("=" * 78)


rna_tables = []


for index, item in enumerate(
    rna_inputs,
    start=1
):

    path = item["file"]
    stage = item["stage"]


    print(
        f"\nReading RNA {index}: "
        f"{stage}"
    )


    rna = pd.read_csv(
        path,
        sep="\t"
    )


    required_columns = {
        "Name",
        "TPM",
        "NumReads"
    }


    missing = (
        required_columns
        - set(rna.columns)
    )


    if missing:

        raise ValueError(
            f"\nERROR: RNA file is missing "
            f"required columns:\n"
            f"{path}\n"
            +
            ", ".join(
                sorted(missing)
            )
        )


    print(
        f"Records: {len(rna):,}"
    )


    # ------------------------------------------------------
    # CLEAN IDS
    # ------------------------------------------------------

    rna["gene_id"] = (
        rna["Name"]
        .apply(clean_gene_id)
    )


    rna["TPM"] = pd.to_numeric(
        rna["TPM"],
        errors="coerce"
    ).fillna(0)


    rna["NumReads"] = pd.to_numeric(
        rna["NumReads"],
        errors="coerce"
    ).fillna(0)


    # ------------------------------------------------------
    # KEEP REQUIRED DATA
    # ------------------------------------------------------

    rna_small = rna[
        [
            "gene_id",
            "TPM",
            "NumReads"
        ]
    ].copy()


    # ------------------------------------------------------
    # HANDLE DUPLICATE GENE IDS
    # ------------------------------------------------------

    rna_small = (
        rna_small
        .groupby(
            "gene_id",
            as_index=False
        )
        .agg({
            "TPM": "sum",
            "NumReads": "sum"
        })
    )


    # ------------------------------------------------------
    # RENAME SAMPLE-SPECIFIC COLUMNS
    # ------------------------------------------------------

    sample_number = index


    rna_small.rename(
        columns={
            "TPM":
                f"TPM_sample_{sample_number}",
            "NumReads":
                f"NumReads_sample_{sample_number}"
        },
        inplace=True
    )


    rna_small[
        f"stage_sample_{sample_number}"
    ] = stage


    rna_tables.append(
        rna_small
    )


# ==========================================================
# MERGE RNA SAMPLES
# ==========================================================

print("\n" + "=" * 78)
print("COMBINING RNA SAMPLES")
print("=" * 78)


rna_combined = rna_tables[0].copy()


for table in rna_tables[1:]:

    rna_combined = rna_combined.merge(
        table,
        on="gene_id",
        how="outer"
    )


print(
    f"Combined RNA genes: "
    f"{len(rna_combined):,}"
)


# ==========================================================
# FILL ONLY SAMPLE-SPECIFIC NUMERIC COLUMNS
# ==========================================================

sample_tpm_columns = [
    c
    for c in rna_combined.columns
    if c.startswith("TPM_sample_")
]


sample_read_columns = [
    c
    for c in rna_combined.columns
    if c.startswith("NumReads_sample_")
]


for column in (
    sample_tpm_columns
    +
    sample_read_columns
):

    rna_combined[column] = (
        pd.to_numeric(
            rna_combined[column],
            errors="coerce"
        )
        .fillna(0)
    )


# ==========================================================
# CALCULATE LIFE-STAGE EXPRESSION
# ==========================================================

print("\n" + "=" * 78)
print("CALCULATING LIFE-STAGE EXPRESSION")
print("=" * 78)


stage_sample_columns = {}


for index, item in enumerate(
    rna_inputs,
    start=1
):

    stage = item["stage"]

    tpm_column = (
        f"TPM_sample_{index}"
    )


    stage_sample_columns.setdefault(
        stage,
        []
    ).append(
        tpm_column
    )


# ==========================================================
# CALCULATE MEAN TPM PER STAGE
# ==========================================================

for stage, columns in (
    stage_sample_columns.items()
):

    valid_columns = [
        c
        for c in columns
        if c in rna_combined.columns
    ]


    if not valid_columns:
        continue


    rna_combined[
        f"{stage}_TPM"
    ] = (
        rna_combined[
            valid_columns
        ]
        .mean(axis=1)
    )


# ==========================================================
# ONLY CREATE STAGES THAT WERE ACTUALLY PROVIDED
# ==========================================================

provided_stages = set(
    stage_sample_columns.keys()
)


print(
    "Stages provided: "
    +
    ", ".join(
        sorted(provided_stages)
    )
)


# ==========================================================
# MERGE WITH TOP-500
# ==========================================================

print("\n" + "=" * 78)
print("MATCHING CANDIDATES TO RNA")
print("=" * 78)


results = candidates.merge(
    rna_combined,
    on="gene_id",
    how="left"
)


# ==========================================================
# RNA MATCH STATUS
# ==========================================================

sample_tpm_columns = [
    c
    for c in rna_combined.columns
    if c.startswith("TPM_sample_")
]


results["RNA_matched"] = (
    results[sample_tpm_columns]
    .notna()
    .any(axis=1)
)


matched_count = (
    results["RNA_matched"]
    .sum()
)


unmatched_count = (
    len(results)
    - matched_count
)


match_rate = (
    matched_count
    /
    len(results)
    *
    100
)


print(
    f"Matched:    {matched_count:,}"
)


print(
    f"Unmatched:  {unmatched_count:,}"
)


print(
    f"Match rate: {match_rate:.2f}%"
)


# ==========================================================
# ENSURE LARVA TPM EXISTS
# ==========================================================

if "Larva_TPM" not in results.columns:

    raise ValueError(
        "\nERROR: Larva TPM could not be calculated."
    )


results["Larva_TPM"] = pd.to_numeric(
    results["Larva_TPM"],
    errors="coerce"
).fillna(0)


# ==========================================================
# LARVAL EXPRESSION
# ==========================================================

results["Larva_log2_TPM"] = np.log2(
    results["Larva_TPM"] + 1
)


# ==========================================================
# NORMALIZE LARVAL EXPRESSION
# ==========================================================

max_larva_log = (
    results["Larva_log2_TPM"]
    .max()
)


if max_larva_log > 0:

    results["Larva_expression_score"] = (

        results["Larva_log2_TPM"]
        /
        max_larva_log
        *
        100
    )

else:

    results[
        "Larva_expression_score"
    ] = 0.0


# ==========================================================
# MULTI-STAGE ANALYSIS
# ==========================================================

if RNA_MODE == "MULTI_STAGE":

    # ------------------------------------------------------
    # CREATE NON-LARVAL STAGE MEAN
    # ONLY USING ACTUALLY PROVIDED STAGES
    # ------------------------------------------------------

    non_larval_stages = [
        stage
        for stage in [
            "Adult",
            "Pupa",
            "Egg"
        ]
        if stage in provided_stages
    ]


    non_larval_columns = [
        f"{stage}_TPM"
        for stage in non_larval_stages
        if f"{stage}_TPM" in results.columns
    ]


    if non_larval_columns:

        results[
            "Non_larval_mean_TPM"
        ] = (
            results[
                non_larval_columns
            ]
            .mean(axis=1)
        )


        results[
            "Non_larval_log2_TPM"
        ] = np.log2(
            results[
                "Non_larval_mean_TPM"
            ]
            + 1
        )


        # --------------------------------------------------
        # LARVAL ENRICHMENT
        # --------------------------------------------------
        #
        # +1 pseudocount prevents division by zero.
        #
        # Example:
        #
        # Larva = 150
        # Other stages = 2
        #
        # -> strong enrichment
        #
        # Larva = 150
        # Other stages = 100
        #
        # -> weaker enrichment
        #

        results[
            "Larva_vs_nonLarva_ratio"
        ] = (

            (
                results["Larva_TPM"]
                + 1
            )

            /

            (
                results[
                    "Non_larval_mean_TPM"
                ]
                + 1
            )
        )


        results[
            "Larva_vs_nonLarva_ratio"
        ] = (
            results[
                "Larva_vs_nonLarva_ratio"
            ]
            .clip(
                lower=0,
                upper=100
            )
        )


        # --------------------------------------------------
        # FIXED SCALE SPECIFICITY SCORE
        # --------------------------------------------------
        #
        # 1x enrichment  -> 0
        # 2x enrichment  -> ~20
        # 4x enrichment  -> ~40
        # 8x enrichment  -> ~60
        # 16x enrichment -> ~80
        # 32x enrichment -> 100
        #

        results[
            "Larva_specificity_score"
        ] = (

            np.log2(
                results[
                    "Larva_vs_nonLarva_ratio"
                ]
            )
            .clip(
                lower=0,
                upper=5
            )
            /
            5
            *
            100
        )


        # --------------------------------------------------
        # LIFE-STAGE LOG EXPRESSION
        # --------------------------------------------------

        for stage in provided_stages:

            column = f"{stage}_TPM"

            if column not in results.columns:
                continue

            results[
                f"{stage}_log2_TPM"
            ] = np.log2(
                results[column] + 1
            )


        # --------------------------------------------------
        # STAGE DOMINANCE
        # --------------------------------------------------

        stage_log_columns = [
            f"{stage}_log2_TPM"
            for stage in provided_stages
            if f"{stage}_log2_TPM"
            in results.columns
        ]


        if stage_log_columns:

            results[
                "Highest_expression_stage"
            ] = (
                results[
                    stage_log_columns
                ]
                .idxmax(axis=1)
                .str.replace(
                    "_log2_TPM",
                    "",
                    regex=False
                )
            )


            results[
                "Larva_stage_dominant"
            ] = (
                results[
                    "Highest_expression_stage"
                ]
                == "Larva"
            )

        else:

            results[
                "Highest_expression_stage"
            ] = "Unknown"

            results[
                "Larva_stage_dominant"
            ] = False


        # --------------------------------------------------
        # MULTI-STAGE RNA SCORE
        # --------------------------------------------------

        results["RNA_score"] = (

            results[
                "Larva_expression_score"
            ]
            * 0.60

            +

            results[
                "Larva_specificity_score"
            ]
            * 0.40
        )


        # --------------------------------------------------
        # SMALL LARVAL DOMINANCE BONUS
        # --------------------------------------------------

        results["RNA_score"] = (

            results["RNA_score"]

            +

            np.where(
                results[
                    "Larva_stage_dominant"
                ],
                5.0,
                0.0
            )
        )


        results["RNA_score"] = (
            results["RNA_score"]
            .clip(
                lower=0,
                upper=100
            )
        )


else:

    # ======================================================
    # LARVA-ONLY MODE
    # ======================================================
    #
    # IMPORTANT:
    #
    # Adult/Pupa/Egg are NOT available.
    #
    # Therefore they are NOT assigned TPM = 0.
    #
    # RNA score is based ONLY on larval expression.
    #

    results[
        "Non_larval_mean_TPM"
    ] = np.nan


    results[
        "Non_larval_log2_TPM"
    ] = np.nan


    results[
        "Larva_vs_nonLarva_ratio"
    ] = np.nan


    results[
        "Larva_specificity_score"
    ] = np.nan


    results[
        "Highest_expression_stage"
    ] = "Larva_only_data"


    results[
        "Larva_stage_dominant"
    ] = np.nan


    results["RNA_score"] = (
        results[
            "Larva_expression_score"
        ]
    )


# ==========================================================
# EXPRESSION CATEGORY
# ==========================================================

def expression_category(tpm):

    if tpm >= 100:
        return "Very High"

    elif tpm >= 10:
        return "High"

    elif tpm >= 1:
        return "Moderate"

    elif tpm > 0:
        return "Low"

    else:
        return "Not detected"


results["Larva_expression"] = (
    results["Larva_TPM"]
    .apply(expression_category)
)


# ==========================================================
# COMBINE CODE 1 + RNA
# ==========================================================

results[
    "annotation_score_before_RNA"
] = results["final_score"]


results["final_score_RNA"] = (

    results["final_score"]
    *
    annotation_weight
    /
    100.0

    +

    results["RNA_score"]
    *
    rna_weight
    /
    100.0
)


# ==========================================================
# RANK
# ==========================================================

print("\n" + "=" * 78)
print("RANKING RNA-PRIORITIZED CANDIDATES")
print("=" * 78)


sort_columns = [
    "final_score_RNA",
    "RNA_score",
    "Larva_expression_score",
    "Larva_TPM"
]


if RNA_MODE == "MULTI_STAGE":

    sort_columns.insert(
        2,
        "Larva_specificity_score"
    )


# Additional biological tie-breakers

for column in [
    "final_score",
    "homology_score",
    "taxonomy_score",
    "best_identity_pct",
    "best_query_coverage_pct",
    "best_bitscore"
]:

    if column in results.columns:

        sort_columns.append(
            column
        )


results = results.sort_values(
    by=sort_columns,
    ascending=[
        False
    ] * len(sort_columns),
    na_position="last"
).reset_index(
    drop=True
)


results["combined_rank"] = (
    np.arange(
        len(results)
    )
    +
    1
)


# ==========================================================
# TOP N
# ==========================================================

top_results = (
    results
    .head(top_n)
    .copy()
)


# ==========================================================
# REORDER IMPORTANT COLUMNS
# ==========================================================

priority_columns = [

    "gene_id",

    "RNA_rank",

    "final_score_RNA",

    "annotation_score_before_RNA",

    "final_score",

    "RNA_score",

    "Larva_expression_score",

    "Larva_specificity_score",

    "Larva_TPM",

    "Adult_TPM",

    "Pupa_TPM",

    "Egg_TPM",

    "Non_larval_mean_TPM",

    "Larva_vs_nonLarva_ratio",

    "Larva_expression",

    "Highest_expression_stage",

    "Larva_stage_dominant",

    "RNA_matched"
]


priority_columns = [
    x
    for x in priority_columns
    if x in top_results.columns
]


remaining_columns = [
    x
    for x in top_results.columns
    if x not in priority_columns
]


top_results = top_results[
    priority_columns
    +
    remaining_columns
]


# ==========================================================
# RNA SUMMARY STATISTICS
# ==========================================================

very_high = (
    results["Larva_TPM"] >= 100
).sum()


high = (
    (results["Larva_TPM"] >= 10)
    &
    (results["Larva_TPM"] < 100)
).sum()


moderate = (
    (results["Larva_TPM"] >= 1)
    &
    (results["Larva_TPM"] < 10)
).sum()


low = (
    (results["Larva_TPM"] > 0)
    &
    (results["Larva_TPM"] < 1)
).sum()


not_detected = (
    results["Larva_TPM"] == 0
).sum()


if RNA_MODE == "MULTI_STAGE":

    larva_dominant_count = (
        results[
            "Larva_stage_dominant"
        ]
        .fillna(False)
        .sum()
    )

else:

    larva_dominant_count = np.nan


# ==========================================================
# WRITE TSV
# ==========================================================

print("\n" + "=" * 78)
print("WRITING OUTPUT")
print("=" * 78)


top_results.to_csv(
    top_file,
    sep="\t",
    index=False
)


results.to_csv(
    all_file,
    sep="\t",
    index=False
)


print(
    f"Created: {top_file}"
)


print(
    f"Created: {all_file}"
)


# ==========================================================
# SUMMARY FILE
# ==========================================================

with open(
    summary_file,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "METISA PLANA — RNA LIFE-STAGE "
        "PRIORITIZATION\n"
    )

    f.write(
        "=" * 60
        +
        "\n\n"
    )


    f.write(
        "RNA ANALYSIS MODE\n"
    )

    f.write(
        "-" * 40
        +
        "\n"
    )

    f.write(
        f"{RNA_MODE}\n\n"
    )


    f.write(
        "RNA SAMPLES\n"
    )

    f.write(
        "-" * 40
        +
        "\n"
    )


    for i, item in enumerate(
        rna_inputs,
        start=1
    ):

        f.write(
            f"RNA {i}: "
            f"{item['stage']}\n"
        )

        f.write(
            f"File: "
            f"{item['file']}\n\n"
        )


    f.write(
        "PIPELINE\n"
    )

    f.write(
        "-" * 40
        +
        "\n"
    )

    f.write(
        "Input: Top-500 candidates "
        "from Code 1\n"
    )


    if RNA_MODE == "LARVA_ONLY":

        f.write(
            "RNA prioritization: "
            "Larval expression only\n"
        )

        f.write(
            "Specificity: Not calculated "
            "(non-larval data unavailable)\n\n"
        )

    else:

        f.write(
            "RNA prioritization: "
            "Larval expression + "
            "larval specificity + "
            "stage dominance\n\n"
        )


    f.write(
        "INPUT STATISTICS\n"
    )

    f.write(
        "-" * 40
        +
        "\n"
    )

    f.write(
        f"Input candidates: "
        f"{len(candidates)}\n"
    )

    f.write(
        f"RNA files: "
        f"{len(rna_inputs)}\n"
    )

    f.write(
        f"RNA genes combined: "
        f"{len(rna_combined):,}\n"
    )

    f.write(
        f"Matched candidates: "
        f"{matched_count}\n"
    )

    f.write(
        f"Unmatched candidates: "
        f"{unmatched_count}\n"
    )

    f.write(
        f"Match rate: "
        f"{match_rate:.2f}%\n\n"
    )


    f.write(
        "LARVAL EXPRESSION\n"
    )

    f.write(
        "-" * 40
        +
        "\n"
    )

    f.write(
        f"Very High (>=100 TPM): "
        f"{very_high}\n"
    )

    f.write(
        f"High (10-<100 TPM): "
        f"{high}\n"
    )

    f.write(
        f"Moderate (1-<10 TPM): "
        f"{moderate}\n"
    )

    f.write(
        f"Low (>0-<1 TPM): "
        f"{low}\n"
    )

    f.write(
        f"Not detected (0 TPM): "
        f"{not_detected}\n"
    )


    if RNA_MODE == "MULTI_STAGE":

        f.write(
            f"Larva highest-expression stage: "
            f"{larva_dominant_count}\n\n"
        )

    else:

        f.write(
            "Larva highest-expression stage: "
            "Not assessed\n\n"
        )


    f.write(
        "SCORING\n"
    )

    f.write(
        "-" * 40
        +
        "\n"
    )

    f.write(
        f"Annotation weight: "
        f"{annotation_weight:.1f}%\n"
    )

    f.write(
        f"RNA weight: "
        f"{rna_weight:.1f}%\n"
    )


    if RNA_MODE == "LARVA_ONLY":

        f.write(
            "RNA score composition:\n"
        )

        f.write(
            "  Larval expression: 100%\n"
        )

    else:

        f.write(
            "RNA score composition:\n"
        )

        f.write(
            "  Larval expression: 60%\n"
        )

        f.write(
            "  Larval specificity: 40%\n"
        )

        f.write(
            "  Larval stage dominance bonus: "
            "+5 points\n"
        )


    f.write(
        "\nOUTPUT\n"
    )

    f.write(
        "-" * 40
        +
        "\n"
    )

    f.write(
        f"Final candidates retained: "
        f"{len(top_results)}\n"
    )


print(
    f"Created: {summary_file}"
)


# ==========================================================
# EXCEL
# ==========================================================

try:

    with pd.ExcelWriter(
        excel_file,
        engine="openpyxl"
    ) as writer:


        # --------------------------------------------------
        # TOP CANDIDATES
        # --------------------------------------------------

        top_results.to_excel(
            writer,
            sheet_name="Top_Candidates",
            index=False
        )


        # --------------------------------------------------
        # ALL CANDIDATES
        # --------------------------------------------------

        results.to_excel(
            writer,
            sheet_name="All_500",
            index=False
        )


        # --------------------------------------------------
        # RNA STAGE DATA
        # --------------------------------------------------

        stage_summary_columns = [

            "gene_id",

            "Larva_TPM",
            "Adult_TPM",
            "Pupa_TPM",
            "Egg_TPM",

            "Non_larval_mean_TPM",

            "Larva_expression_score",

            "Larva_specificity_score",

            "Larva_vs_nonLarva_ratio",

            "RNA_score",

            "Highest_expression_stage",

            "Larva_stage_dominant"

        ]


        stage_summary_columns = [
            x
            for x in stage_summary_columns
            if x in results.columns
        ]


        results[
            stage_summary_columns
        ].to_excel(
            writer,
            sheet_name="RNA_Stage_Expression",
            index=False
        )


        # --------------------------------------------------
        # SUMMARY
        # --------------------------------------------------

        summary_table = pd.DataFrame({

            "Metric": [

                "RNA analysis mode",

                "Input candidates",

                "RNA files",

                "RNA genes combined",

                "Matched candidates",

                "Unmatched candidates",

                "Match rate (%)",

                "Very High Larval TPM",

                "High Larval TPM",

                "Moderate Larval TPM",

                "Low Larval TPM",

                "Not detected",

                "Larva highest-expression stage",

                "Annotation weight (%)",

                "RNA weight (%)",

                "Final candidates"

            ],

            "Value": [

                RNA_MODE,

                len(candidates),

                len(rna_inputs),

                len(rna_combined),

                matched_count,

                unmatched_count,

                match_rate,

                very_high,

                high,

                moderate,

                low,

                not_detected,

                (
                    larva_dominant_count
                    if RNA_MODE == "MULTI_STAGE"
                    else "Not assessed"
                ),

                annotation_weight,

                rna_weight,

                len(top_results)

            ]

        })


        summary_table.to_excel(
            writer,
            sheet_name="Summary",
            index=False
        )


    print(
        f"Created: {excel_file}"
    )


except Exception as e:

    print(
        "\nWARNING: Excel output "
        "could not be created."
    )

    print(e)


# ==========================================================
# SHOW TOP 20
# ==========================================================

print("\n" + "=" * 78)
print("TOP 20 AFTER RNA PRIORITIZATION")
print("=" * 78)


display_columns = [

    "gene_id",

    "combined_rank",

    "final_score_RNA",

    "final_score",

    "RNA_score",

    "Larva_expression_score",

    "Larva_specificity_score",

    "Larva_TPM",

    "Adult_TPM",

    "Pupa_TPM",

    "Egg_TPM",

    "Highest_expression_stage",

    "Larva_expression"

]


for column in [

    "best_hit_species",

    "functional_description",

    "target_class",

    "decision"

]:

    if column in top_results.columns:

        display_columns.append(
            column
        )


print(
    top_results[
        [
            c
            for c in display_columns
            if c in top_results.columns
        ]
    ]
    .head(20)
    .to_string(
        index=False
    )
)


# ==========================================================
# FINAL
# ==========================================================

print("\n" + "=" * 78)
print("DONE")
print("=" * 78)


print(
    f"Top {len(top_results)} "
    "candidates retained."
)


print(
    "\nRNA prioritization mode:"
)


if RNA_MODE == "LARVA_ONLY":

    print(
        "  Larva-only expression ranking"
    )

    print(
        "  High larval expression → higher priority"
    )

    print(
        "  Low/no larval expression → lower priority"
    )

else:

    print(
        "  1. Larval expression"
    )

    print(
        "  2. Larval specificity "
        "relative to available non-larval stages"
    )

    print(
        "  3. Whether larva is the "
        "highest-expression stage"
    )


print(
    "\nNext step: inspect the "
    "RNA-prioritized candidates "
    "before MAFFT/phylogenetic analysis."
)
