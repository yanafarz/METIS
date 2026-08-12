import os
import csv
import re
from collections import Counter


# ==========================================================
# CONFIGURATION
# ==========================================================

# Your chosen BLAST hit-quality thresholds
MIN_IDENTITY = 30.0
MIN_COVERAGE = 70.0
MAX_EVALUE = 0.05


# ==========================================================
# ORGANISM KEYWORDS
# ==========================================================

# IMPORTANT:
# These are screening keywords, NOT proof of contamination.
# A flagged sequence must be reviewed before removal.

BACTERIAL_KEYWORDS = [
    "bacter",
    "escherichia",
    "bacillus",
    "streptococcus",
    "staphylococcus",
    "pseudomonas",
    "salmonella",
    "enterobacter",
    "klebsiella",
    "vibrio",
    "lactobacillus",
    "clostridium",
    "mycobacterium",
    "cyanobacter",
    "actinobacter",
    "firmicutes",
    "proteobacter",
    "bacteroid",
    "spirochaet",
    "archaea",
    "archae",
]


FUNGAL_KEYWORDS = [
    "fung",
    "candida",
    "aspergillus",
    "saccharomyces",
    "schizosaccharomyces",
    "neurospora",
    "penicillium",
    "trichoderma",
    "cryptococcus",
    "kluyveromyces",
    "mucor",
    "rhizopus",
    "botrytis",
    "yeast",
]


# Broad arthropod/insect keywords.
# These indicate that a hit is probably from an arthropod.

ARTHROPOD_KEYWORDS = [
    "bombyx",
    "spodoptera",
    "manduca",
    "drosophila",
    "helicoverpa",
    "plutella",
    "trichoplusia",
    "epiphyas",
    "pieris",
    "papilio",
    "nasonia",
    "apis",
    "honeybee",
    "bee",
    "wasp",
    "ant",
    "formica",
    "camponotus",
    "harpegnathos",
    "tribolium",
    "tenebrio",
    "acyrthosiphon",
    "aphid",
    "bemisia",
    "aedes",
    "anopheles",
    "culex",
    "mosquito",
    "musca",
    "fly",
    "moth",
    "butterfly",
    "beetle",
    "locust",
    "grasshopper",
    "cricket",
    "tick",
    "mite",
    "spider",
    "arachnid",
    "arthropod",
    "insect",
    "lepidoptera",
    "diptera",
    "coleoptera",
    "hymenoptera",
    "hemiptera",
    "orthoptera",
    "odonata",
    "blattodea",
    "dermaptera",
    "thysanoptera",
]


# ==========================================================
# FILE FUNCTIONS
# ==========================================================

def expand_path(path):

    return os.path.abspath(
        os.path.expanduser(
            path.strip()
        )
    )


def check_file_exists(path):

    if not os.path.isfile(path):

        print("\nERROR: File not found:")
        print(path)

        return False

    return True


# ==========================================================
# READ BLAST RESULT
# ==========================================================

EXPECTED_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "qcovs",
    "evalue",
    "bitscore",
    "stitle",
]


def read_blast_file(path):

    print("\nReading BLAST result...")

    rows = []

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as file:

        reader = csv.reader(
            file,
            delimiter="\t"
        )

        for line_number, row in enumerate(
            reader,
            start=1
        ):

            # Skip empty lines
            if not row:
                continue

            # Skip manually added header
            if row[0].strip() == "qseqid":
                continue

            # BLAST should contain 10 columns
            if len(row) < 10:

                print(
                    f"WARNING: Skipping line {line_number}: "
                    f"only {len(row)} columns."
                )

                continue

            try:

                record = {

                    "qseqid": row[0].strip(),

                    "sseqid": row[1].strip(),

                    "pident": float(row[2]),

                    "length": int(float(row[3])),

                    "qlen": int(float(row[4])),

                    "slen": int(float(row[5])),

                    "qcovs": float(row[6]),

                    "evalue": float(row[7]),

                    "bitscore": float(row[8]),

                    "stitle": row[9].strip(),
                }

                rows.append(record)

            except ValueError:

                print(
                    f"WARNING: Could not parse line "
                    f"{line_number}."
                )

    print(
        f"BLAST hit rows read: {len(rows):,}"
    )

    return rows


# ==========================================================
# EXTRACT SPECIES
# ==========================================================

def extract_species(stitle):

    """
    Extract organism from UniProt-style descriptions.

    Example:

    sp|P12345|ABC_BOMMO Protein OS=Bombyx mori
    OX=...
    GN=...

    Returns:

    Bombyx mori
    """

    if not stitle:

        return "Unknown"

    # Standard UniProt OS= field

    match = re.search(
        r"\bOS=([^=]+?)(?=\s+(?:OX|GN|PE|SV|CC|KW|GO)=|$)",
        stitle
    )

    if match:

        species = match.group(1).strip()

        if species:

            return species

    # Alternative organism formats

    match = re.search(
        r"\borganism[=:]\s*([^,;]+)",
        stitle,
        flags=re.IGNORECASE
    )

    if match:

        species = match.group(1).strip()

        if species:

            return species

    return "Unknown"


# ==========================================================
# KEYWORD CHECK
# ==========================================================

def contains_keyword(text, keywords):

    text = text.lower()

    for keyword in keywords:

        if keyword.lower() in text:

            return True

    return False


# ==========================================================
# ORGANISM CLASSIFICATION
# ==========================================================

def classify_organism(species):

    text = species.lower()

    if text == "unknown":

        return "Unknown"

    # Bacteria first

    if contains_keyword(
        text,
        BACTERIAL_KEYWORDS
    ):

        return "Bacterial"

    # Fungi

    if contains_keyword(
        text,
        FUNGAL_KEYWORDS
    ):

        return "Fungal"

    # Arthropods

    if contains_keyword(
        text,
        ARTHROPOD_KEYWORDS
    ):

        return "Arthropod"

    # Everything else

    return "Other non-arthropod"


# ==========================================================
# HIT QUALITY FILTER
# ==========================================================

def passes_quality_filter(hit):

    """
    Determines whether a BLAST hit is strong enough
    to be used for contamination interpretation.

    Your chosen thresholds:

        Identity >= 30%
        Coverage >= 70%
        E-value <= 0.05
    """

    return (

        hit["pident"] >= MIN_IDENTITY

        and

        hit["qcovs"] >= MIN_COVERAGE

        and

        hit["evalue"] <= MAX_EVALUE

    )


# ==========================================================
# CHOOSE BEST QUALIFYING HIT
# ==========================================================

def select_best_hit(hits):

    qualifying = []

    for hit in hits:

        if passes_quality_filter(hit):

            qualifying.append(hit)

    if not qualifying:

        return None, []


    # Rank qualifying hits.
    #
    # Primary:
    #   lowest E-value
    #
    # Secondary:
    #   highest bitscore
    #
    # Tertiary:
    #   highest identity
    #
    # Fourth:
    #   highest coverage

    qualifying.sort(

        key=lambda x: (

            x["evalue"],

            -x["bitscore"],

            -x["pident"],

            -x["qcovs"],

        )

    )

    return (
        qualifying[0],
        qualifying
    )


# ==========================================================
# CONTAMINATION DECISION
# ==========================================================

def decide_action(classification):

    if classification == "Bacterial":

        return "FLAG_POSSIBLE_CONTAMINATION"

    if classification == "Fungal":

        return "FLAG_POSSIBLE_CONTAMINATION"

    if classification == "Other non-arthropod":

        return "FLAG_REVIEW"

    if classification == "Arthropod":

        return "KEEP"

    if classification == "Unknown":

        return "REVIEW_UNKNOWN"

    return "REVIEW"


# ==========================================================
# SCREEN ALL QUERY SEQUENCES
# ==========================================================

def screen_contamination(rows):

    print(
        "\nGrouping BLAST hits by query..."
    )

    grouped = {}

    for row in rows:

        query = row["qseqid"]

        if query not in grouped:

            grouped[query] = []

        grouped[query].append(row)


    print(
        f"Unique query sequences: "
        f"{len(grouped):,}"
    )


    results = []

    print(
        "\nPerforming contamination screening..."
    )


    for query_id, hits in grouped.items():

        best_hit, qualifying_hits = (
            select_best_hit(hits)
        )


        # --------------------------------------------------
        # NO QUALIFYING HIT
        # --------------------------------------------------

        if best_hit is None:

            results.append({

                "gene_id":
                    query_id,

                "best_hit_accession":
                    "",

                "best_hit_species":
                    "",

                "identity_pct":
                    "",

                "query_coverage_pct":
                    "",

                "evalue":
                    "",

                "bitscore":
                    "",

                "classification":
                    "No qualifying hit",

                "qc_decision":
                    "REVIEW_NO_QUALIFYING_HIT",

                "total_blast_hits":
                    len(hits),

                "qualifying_hits":
                    0,

                "best_hit_description":
                    "",

            })

            continue


        # --------------------------------------------------
        # CLASSIFY BEST HIT
        # --------------------------------------------------

        species = extract_species(
            best_hit["stitle"]
        )

        classification = classify_organism(
            species
        )

        action = decide_action(
            classification
        )


        # --------------------------------------------------
        # SAVE RESULT
        # --------------------------------------------------

        results.append({

            "gene_id":
                query_id,

            "best_hit_accession":
                best_hit["sseqid"],

            "best_hit_species":
                species,

            "identity_pct":
                best_hit["pident"],

            "query_coverage_pct":
                best_hit["qcovs"],

            "evalue":
                best_hit["evalue"],

            "bitscore":
                best_hit["bitscore"],

            "classification":
                classification,

            "qc_decision":
                action,

            "total_blast_hits":
                len(hits),

            "qualifying_hits":
                len(qualifying_hits),

            "best_hit_description":
                best_hit["stitle"],

        })


    return results


# ==========================================================
# SAVE DETAILED RESULTS
# ==========================================================

def save_detailed_results(
    results,
    output_path
):

    fields = [

        "gene_id",

        "best_hit_accession",

        "best_hit_species",

        "identity_pct",

        "query_coverage_pct",

        "evalue",

        "bitscore",

        "classification",

        "qc_decision",

        "total_blast_hits",

        "qualifying_hits",

        "best_hit_description",

    ]


    with open(
        output_path,
        "w",
        newline="",
        encoding="utf-8"
    ) as file:

        writer = csv.DictWriter(
            file,
            fieldnames=fields
        )

        writer.writeheader()

        writer.writerows(results)


    print(
        "\nDetailed screening result saved:"
    )

    print(output_path)


# ==========================================================
# SAVE FLAGGED RESULTS
# ==========================================================

def save_flagged_results(
    results,
    output_path
):

    flagged = [

        row

        for row in results

        if row["qc_decision"]

        in [

            "FLAG_POSSIBLE_CONTAMINATION",

            "FLAG_REVIEW",

        ]

    ]


    fields = [

        "gene_id",

        "best_hit_accession",

        "best_hit_species",

        "identity_pct",

        "query_coverage_pct",

        "evalue",

        "bitscore",

        "classification",

        "qc_decision",

        "total_blast_hits",

        "qualifying_hits",

        "best_hit_description",

    ]


    with open(
        output_path,
        "w",
        newline="",
        encoding="utf-8"
    ) as file:

        writer = csv.DictWriter(
            file,
            fieldnames=fields
        )

        writer.writeheader()

        writer.writerows(flagged)


    print(
        "\nFlagged sequences saved:"
    )

    print(output_path)

    return len(flagged)


# ==========================================================
# SAVE SUMMARY
# ==========================================================

def save_summary(
    results,
    output_path,
    blast_file
):

    classification_counts = Counter(

        row["classification"]

        for row in results

    )


    action_counts = Counter(

        row["qc_decision"]

        for row in results

    )


    total = len(results)


    with open(
        output_path,
        "w",
        encoding="utf-8"
    ) as file:

        file.write(
            "METISA PLANA BLAST CONTAMINATION SCREENING\n"
        )

        file.write(
            "=" * 65 + "\n\n"
        )


        file.write(
            "INPUT BLAST FILE\n"
        )

        file.write(
            "-" * 65 + "\n"
        )

        file.write(
            f"{blast_file}\n\n"
        )


        file.write(
            "CONTAMINATION SCREENING THRESHOLDS\n"
        )

        file.write(
            "-" * 65 + "\n"
        )

        file.write(
            f"Minimum identity: "
            f"{MIN_IDENTITY}%\n"
        )

        file.write(
            f"Minimum query coverage: "
            f"{MIN_COVERAGE}%\n"
        )

        file.write(
            f"Maximum E-value: "
            f"{MAX_EVALUE}\n\n"
        )


        file.write(
            "RESULTS\n"
        )

        file.write(
            "-" * 65 + "\n"
        )

        file.write(
            f"Queries screened: "
            f"{total:,}\n\n"
        )


        file.write(
            "ORGANISM CLASSIFICATION\n"
        )

        file.write(
            "-" * 65 + "\n"
        )


        categories = [

            "Arthropod",

            "Bacterial",

            "Fungal",

            "Other non-arthropod",

            "Unknown",

            "No qualifying hit",

        ]


        for category in categories:

            count = classification_counts.get(
                category,
                0
            )

            percentage = (

                count / total * 100

                if total > 0

                else 0

            )


            file.write(

                f"{category}: "
                f"{count:,} "
                f"({percentage:.2f}%)\n"

            )


        file.write(
            "\nQC DECISIONS\n"
        )

        file.write(
            "-" * 65 + "\n"
        )


        for action, count in action_counts.items():

            percentage = (

                count / total * 100

                if total > 0

                else 0

            )


            file.write(

                f"{action}: "
                f"{count:,} "
                f"({percentage:.2f}%)\n"

            )


        file.write(
            "\nIMPORTANT INTERPRETATION\n"
        )

        file.write(
            "-" * 65 + "\n"
        )

        file.write(
            "Bacterial and fungal classifications "
            "are possible contamination flags, "
            "not confirmed contamination.\n"
        )

        file.write(
            "Flagged sequences should be manually "
            "reviewed before removal.\n"
        )

        file.write(
            "The organism classification is based "
            "on the BLAST subject description and "
            "keyword screening.\n"
        )


    print(
        "\nSummary report saved:"
    )

    print(output_path)


# ==========================================================
# MAIN PROGRAM
# ==========================================================

def main():

    print(
        "\n"
        + "=" * 65
    )

    print(
        "       METISA PLANA CONTAMINATION SCREEN"
    )

    print(
        "=" * 65
    )


    # ------------------------------------------------------
    # INPUT
    # ------------------------------------------------------

    print(
        "\nSTEP 1: BLAST RESULT"
    )

    blast_file = expand_path(

        input(
            "Path to BLAST TSV result: "
        )

    )


    if not check_file_exists(
        blast_file
    ):

        return


    # ------------------------------------------------------
    # SHOW THRESHOLDS
    # ------------------------------------------------------

    print(
        "\nSTEP 2: SCREENING THRESHOLDS"
    )

    print(
        "\nUsing your selected thresholds:"
    )

    print(
        f"  Identity  >= {MIN_IDENTITY}%"
    )

    print(
        f"  Coverage  >= {MIN_COVERAGE}%"
    )

    print(
        f"  E-value   <= {MAX_EVALUE}"
    )


    # ------------------------------------------------------
    # READ BLAST
    # ------------------------------------------------------

    rows = read_blast_file(
        blast_file
    )


    if not rows:

        print(
            "\nERROR: No usable BLAST rows found."
        )

        return


    # ------------------------------------------------------
    # SCREEN
    # ------------------------------------------------------

    results = screen_contamination(
        rows
    )


    # ------------------------------------------------------
    # OUTPUT DIRECTORY
    # ------------------------------------------------------

    output_dir = os.path.join(

        os.path.dirname(
            blast_file
        ),

        "contamination_screen"

    )


    os.makedirs(
        output_dir,
        exist_ok=True
    )


    # ------------------------------------------------------
    # OUTPUT FILES
    # ------------------------------------------------------

    detailed_file = os.path.join(

        output_dir,

        "contamination_screening.csv"

    )


    flagged_file = os.path.join(

        output_dir,

        "flagged_possible_contamination.csv"

    )


    summary_file = os.path.join(

        output_dir,

        "contamination_summary.txt"

    )


    # ------------------------------------------------------
    # SAVE RESULTS
    # ------------------------------------------------------

    save_detailed_results(

        results,

        detailed_file

    )


    flagged_count = save_flagged_results(

        results,

        flagged_file

    )


    save_summary(

        results,

        summary_file,

        blast_file

    )


    # ------------------------------------------------------
    # TERMINAL SUMMARY
    # ------------------------------------------------------

    counts = Counter(

        row["classification"]

        for row in results

    )


    print(
        "\n"
        + "=" * 65
    )

    print(
        "CONTAMINATION SCREEN SUMMARY"
    )

    print(
        "=" * 65
    )


    print(
        f"Sequences screened: "
        f"{len(results):,}"
    )


    print(
        f"Arthropod: "
        f"{counts.get('Arthropod', 0):,}"
    )


    print(
        f"Bacterial: "
        f"{counts.get('Bacterial', 0):,}"
    )


    print(
        f"Fungal: "
        f"{counts.get('Fungal', 0):,}"
    )


    print(
        f"Other non-arthropod: "
        f"{counts.get('Other non-arthropod', 0):,}"
    )


    print(
        f"Unknown: "
        f"{counts.get('Unknown', 0):,}"
    )


    print(
        f"No qualifying hit: "
        f"{counts.get('No qualifying hit', 0):,}"
    )


    print(
        f"Flagged for review: "
        f"{flagged_count:,}"
    )


    print(
        "\nOutput folder:"
    )

    print(output_dir)


    print(
        "\nDone."
    )


# ==========================================================
# RUN
# ==========================================================

if __name__ == "__main__":

    main()