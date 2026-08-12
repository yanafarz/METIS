import os
import re
import csv
from collections import defaultdict, Counter

# ============================================================
# CONFIGURATION
# ============================================================

# BLAST thresholds requested by you
MIN_IDENTITY = 30.0
MIN_COVERAGE = 70.0
MAX_EVALUE = 0.05

# Taxonomic dominance required for classification
# Example: 0.70 means 70% or more of qualifying hits
DOMINANCE_THRESHOLD = 0.70

# Minimum number of qualifying hits required before
# making a strong taxonomic classification.
#
# If only 1 qualifying hit exists, we classify as REVIEW
# rather than trusting one hit too strongly.
MIN_HITS_FOR_STRONG_CLASSIFICATION = 2



# ============================================================
# CATEGORY DEFINITIONS
# ============================================================

CATEGORIES = {
    1: "KEEP_ARTHROPOD",
    2: "POSSIBLE_BACTERIAL",
    3: "POSSIBLE_FUNGAL",
    4: "OTHER_EUKARYOTE",
    5: "AMBIGUOUS",
    6: "NO_QUALIFYING_HIT",
    7: "UNKNOWN_TAXONOMY",
}


# ============================================================
# EXPECTED BLAST COLUMNS
# ============================================================

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


# ============================================================
# PATH FUNCTIONS
# ============================================================

def expand_path(path):
    return os.path.abspath(
        os.path.expanduser(
            path.strip().strip('"')
        )
    )


def check_file_exists(path, description):
    if not os.path.isfile(path):
        print("\nERROR: File not found:")
        print(f"{description}: {path}")
        return False

    return True


# ============================================================
# READ BLAST TSV
# ============================================================

def read_blast_file(path):

    print("\nReading BLAST TSV...")

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

            if not row:
                continue

            # Skip header
            if row[0].strip().lower() == "qseqid":
                continue

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
        f"BLAST rows read: {len(rows):,}"
    )

    return rows


# ============================================================
# BLAST QUALITY FILTER
# ============================================================

def passes_quality_filter(hit):

    return (
        hit["pident"] >= MIN_IDENTITY
        and
        hit["qcovs"] >= MIN_COVERAGE
        and
        hit["evalue"] <= MAX_EVALUE
    )


# ============================================================
# EXTRACT NCBI TAXON ID FROM UNIPROT DESCRIPTION
# ============================================================

def extract_taxid(stitle):

    if not stitle:
        return None

    # UniProt format:
    #
    # OS=Eumeta variegata OX=151549 GN=...
    #
    match = re.search(
        r"\bOX=(\d+)",
        stitle
    )

    if match:
        return int(match.group(1))

    return None


# ============================================================
# EXTRACT ORGANISM FROM STITLE
# ============================================================

def extract_species(stitle):

    if not stitle:
        return "Unknown"

    match = re.search(
        r"\bOS=([^=]+?)(?=\s+(?:OX|GN|PE|SV|CC|KW|GO)=|$)",
        stitle
    )

    if match:

        species = match.group(1).strip()

        if species:
            return species

    return "Unknown"


# ============================================================
# PARSE RANKEDLINEAGE.DMP LINE
# ============================================================

def parse_ranked_lineage_line(line):

    """
    NCBI rankedlineage.dmp normally contains:

    tax_id
    scientific_name
    species
    genus
    family
    order
    class
    phylum
    kingdom
    superkingdom
    """

    # Remove whitespace/newlines
    line = line.strip()

    if not line:
        return None

    # Most NCBI dump files use:
    #
    # value | value | value | ...
    #
    # with spaces around separators.

    parts = [
        x.strip()
        for x in line.split("|")
    ]

    if len(parts) < 10:
        return None

    try:
        taxid = int(parts[0])
    except ValueError:
        return None

    return {
        "taxid": taxid,
        "scientific_name": parts[1],
        "species": parts[2],
        "genus": parts[3],
        "family": parts[4],
        "order": parts[5],
        "class": parts[6],
        "phylum": parts[7],
        "kingdom": parts[8],
        "superkingdom": parts[9],
    }


# ============================================================
# LOAD ONLY TAXIDS NEEDED FROM RANKEDLINEAGE.DMP
# ============================================================

def load_required_taxonomy(
    rankedlineage_path,
    required_taxids
):

    print("\nLoading taxonomy from rankedlineage.dmp...")

    required_taxids = set(required_taxids)

    taxonomy = {}

    if not required_taxids:
        print("No TaxIDs found in BLAST results.")
        return taxonomy

    print(
        f"TaxIDs required: {len(required_taxids):,}"
    )

    with open(
        rankedlineage_path,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as file:

        for line in file:

            record = parse_ranked_lineage_line(line)

            if record is None:
                continue

            taxid = record["taxid"]

            if taxid in required_taxids:

                taxonomy[taxid] = record

                # Stop early if everything was found
                if len(taxonomy) == len(required_taxids):
                    break

    print(
        f"TaxIDs successfully found: "
        f"{len(taxonomy):,}/{len(required_taxids):,}"
    )

    return taxonomy


# ============================================================
# CLASSIFY TAXONOMY
# ============================================================

def classify_taxonomy(tax_record):

    if tax_record is None:

        return "UNKNOWN_TAXONOMY"

    superkingdom = (
        tax_record["superkingdom"]
        .strip()
        .lower()
    )

    kingdom = (
        tax_record["kingdom"]
        .strip()
        .lower()
    )

    phylum = (
        tax_record["phylum"]
        .strip()
        .lower()
    )

    # --------------------------------------------------------
    # BACTERIA
    # --------------------------------------------------------

    if superkingdom == "bacteria":

        return "BACTERIAL"

    # --------------------------------------------------------
    # ARCHAEA
    # --------------------------------------------------------

    if superkingdom == "archaea":

        return "OTHER_PROKARYOTE"

    # --------------------------------------------------------
    # FUNGI
    # --------------------------------------------------------

    if kingdom == "fungi":

        return "FUNGAL"

    # --------------------------------------------------------
    # ARTHROPODA
    # --------------------------------------------------------

    if phylum == "arthropoda":

        return "ARTHROPOD"

    # --------------------------------------------------------
    # OTHER EUKARYOTE
    # --------------------------------------------------------

    if superkingdom == "eukaryota":

        return "OTHER_EUKARYOTE"

    return "OTHER"


# ============================================================
# GROUP BLAST HITS
# ============================================================

def group_hits_by_query(rows):

    grouped = defaultdict(list)

    for row in rows:

        grouped[row["qseqid"]].append(row)

    return grouped


# ============================================================
# SELECT QUALIFYING HITS
# ============================================================

def get_qualifying_hits(hits):

    qualifying = [
        hit
        for hit in hits
        if passes_quality_filter(hit)
    ]

    # Sort strongest first
    qualifying.sort(
        key=lambda x: (
            x["evalue"],
            -x["bitscore"],
            -x["pident"],
            -x["qcovs"],
        )
    )

    return qualifying


# ============================================================
# CALCULATE TAXONOMIC COUNTS
# ============================================================

def calculate_taxonomic_support(
    qualifying_hits,
    taxonomy
):

    counts = Counter()

    classified_hits = []

    for hit in qualifying_hits:

        taxid = extract_taxid(
            hit["stitle"]
        )

        tax_record = (
            taxonomy.get(taxid)
            if taxid is not None
            else None
        )

        tax_class = classify_taxonomy(
            tax_record
        )

        counts[tax_class] += 1

        classified_hits.append({
            "hit": hit,
            "taxid": taxid,
            "taxonomy": tax_record,
            "tax_class": tax_class,
        })

    return counts, classified_hits


# ============================================================
# MAKE FINAL DECISION
# ============================================================

def make_decision(
    qualifying_hits,
    counts
):

    total = len(qualifying_hits)

    # --------------------------------------------------------
    # NO QUALIFYING HIT
    # --------------------------------------------------------

    if total == 0:

        return (
            "NO_QUALIFYING_HIT",
            "No BLAST hit passed all thresholds."
        )

    # --------------------------------------------------------
    # UNKNOWN TAXONOMY
    # --------------------------------------------------------

    known = (
        total
        - counts.get("UNKNOWN_TAXONOMY", 0)
    )

    if known == 0:

        return (
            "UNKNOWN_TAXONOMY",
            "Qualifying hits found, but taxonomy "
            "could not be resolved."
        )

    # --------------------------------------------------------
    # SINGLE HIT
    # --------------------------------------------------------

    if total < MIN_HITS_FOR_STRONG_CLASSIFICATION:

        top_class = (
            counts.most_common(1)[0][0]
        )

        if top_class == "BACTERIAL":

            return (
                "AMBIGUOUS",
                "Only one qualifying hit; "
                "bacterial hit requires review."
            )

        if top_class == "FUNGAL":

            return (
                "AMBIGUOUS",
                "Only one qualifying hit; "
                "fungal hit requires review."
            )

        if top_class == "ARTHROPOD":

            return (
                "AMBIGUOUS",
                "Only one qualifying hit; "
                "arthropod assignment requires review."
            )

        return (
            "AMBIGUOUS",
            "Too few qualifying hits for strong "
            "taxonomic classification."
        )

    # --------------------------------------------------------
    # CALCULATE FRACTIONS
    # --------------------------------------------------------

    arthropod_fraction = (
        counts.get("ARTHROPOD", 0) / total
    )

    bacterial_fraction = (
        counts.get("BACTERIAL", 0) / total
    )

    fungal_fraction = (
        counts.get("FUNGAL", 0) / total
    )

    other_euk_fraction = (
        counts.get("OTHER_EUKARYOTE", 0) / total
    )

    # --------------------------------------------------------
    # ARTHROPOD DOMINANT
    # --------------------------------------------------------

    if (
        arthropod_fraction
        >= DOMINANCE_THRESHOLD
    ):

        return (
            "KEEP_ARTHROPOD",
            f"Arthropod support "
            f"{arthropod_fraction:.1%} "
            f"of qualifying hits."
        )

    # --------------------------------------------------------
    # BACTERIAL DOMINANT
    # --------------------------------------------------------

    if (
        bacterial_fraction
        >= DOMINANCE_THRESHOLD
    ):

        return (
            "POSSIBLE_BACTERIAL",
            f"Bacterial support "
            f"{bacterial_fraction:.1%} "
            f"of qualifying hits."
        )

    # --------------------------------------------------------
    # FUNGAL DOMINANT
    # --------------------------------------------------------

    if (
        fungal_fraction
        >= DOMINANCE_THRESHOLD
    ):

        return (
            "POSSIBLE_FUNGAL",
            f"Fungal support "
            f"{fungal_fraction:.1%} "
            f"of qualifying hits."
        )

    # --------------------------------------------------------
    # OTHER EUKARYOTE DOMINANT
    # --------------------------------------------------------

    if (
        other_euk_fraction
        >= DOMINANCE_THRESHOLD
    ):

        return (
            "OTHER_EUKARYOTE",
            f"Other-eukaryote support "
            f"{other_euk_fraction:.1%} "
            f"of qualifying hits."
        )

    # --------------------------------------------------------
    # MIXED
    # --------------------------------------------------------

    return (
        "AMBIGUOUS",
        "No taxonomic group reached the "
        f"{DOMINANCE_THRESHOLD:.0%} dominance threshold."
    )


# ============================================================
# SCREEN ALL PROTEINS
# ============================================================

def screen_proteins(
    rows,
    taxonomy
):

    print("\nGrouping BLAST hits by protein...")

    grouped = group_hits_by_query(rows)

    print(
        f"Unique proteins: {len(grouped):,}"
    )

    results = []

    print("\nScreening contamination...")

    for query_id, hits in grouped.items():

        qualifying_hits = get_qualifying_hits(
            hits
        )

        counts, classified_hits = (
            calculate_taxonomic_support(
                qualifying_hits,
                taxonomy
            )
        )

        category, reason = make_decision(
            qualifying_hits,
            counts
        )

        # ----------------------------------------------------
        # BEST QUALIFYING HIT
        # ----------------------------------------------------

        if qualifying_hits:

            best = qualifying_hits[0]

            best_species = extract_species(
                best["stitle"]
            )

            best_taxid = extract_taxid(
                best["stitle"]
            )

            best_tax = taxonomy.get(
                best_taxid
            )

            if best_tax:

                best_lineage = (
                    best_tax["superkingdom"]
                    + "; "
                    + best_tax["kingdom"]
                    + "; "
                    + best_tax["phylum"]
                    + "; "
                    + best_tax["class"]
                    + "; "
                    + best_tax["order"]
                )

            else:

                best_lineage = ""

            best_accession = best["sseqid"]

            best_identity = best["pident"]

            best_coverage = best["qcovs"]

            best_evalue = best["evalue"]

            best_bitscore = best["bitscore"]

            best_description = best["stitle"]

        else:

            best_species = ""
            best_taxid = ""
            best_lineage = ""
            best_accession = ""
            best_identity = ""
            best_coverage = ""
            best_evalue = ""
            best_bitscore = ""
            best_description = ""

        total = len(qualifying_hits)

        results.append({

            "gene_id":
                query_id,

            "best_hit_accession":
                best_accession,

            "best_hit_species":
                best_species,

            "best_hit_taxid":
                best_taxid,

            "best_hit_lineage":
                best_lineage,

            "best_identity_pct":
                best_identity,

            "best_query_coverage_pct":
                best_coverage,

            "best_evalue":
                best_evalue,

            "best_bitscore":
                best_bitscore,

            "qualifying_hits":
                total,

            "arthropod_hits":
                counts.get(
                    "ARTHROPOD",
                    0
                ),

            "bacterial_hits":
                counts.get(
                    "BACTERIAL",
                    0
                ),

            "fungal_hits":
                counts.get(
                    "FUNGAL",
                    0
                ),

            "other_eukaryote_hits":
                counts.get(
                    "OTHER_EUKARYOTE",
                    0
                ),

            "other_prokaryote_hits":
                counts.get(
                    "OTHER_PROKARYOTE",
                    0
                ),

            "unknown_taxonomy_hits":
                counts.get(
                    "UNKNOWN_TAXONOMY",
                    0
                ),

            "arthropod_fraction":
                (
                    counts.get("ARTHROPOD", 0)
                    / total
                    if total > 0
                    else 0
                ),

            "bacterial_fraction":
                (
                    counts.get("BACTERIAL", 0)
                    / total
                    if total > 0
                    else 0
                ),

            "fungal_fraction":
                (
                    counts.get("FUNGAL", 0)
                    / total
                    if total > 0
                    else 0
                ),

            "classification":
                category,

            "decision":
                category,

            "reason":
                reason,

            "total_blast_hits":
                len(hits),

            "best_hit_description":
                best_description,
        })

    return results


# ============================================================
# WRITE TSV
# ============================================================

def save_tsv(
    rows,
    output_path
):

    if not rows:
        return

    fields = list(rows[0].keys())

    with open(
        output_path,
        "w",
        newline="",
        encoding="utf-8"
    ) as file:

        writer = csv.DictWriter(
            file,
            fieldnames=fields,
            delimiter="\t"
        )

        writer.writeheader()

        writer.writerows(rows)

    print(
        f"Saved TSV: {output_path}"
    )


# ============================================================
# WRITE EXCEL
# ============================================================

def save_excel(
    results,
    output_path
):

    try:

        from openpyxl import Workbook
        from openpyxl.styles import Font

    except ImportError:

        print(
            "\nWARNING: openpyxl is not installed."
        )

        print(
            "Install it with:"
        )

        print(
            "pip install openpyxl"
        )

        return

    workbook = Workbook()

    # --------------------------------------------------------
    # ALL RESULTS
    # --------------------------------------------------------

    sheet = workbook.active
    sheet.title = "All Results"

    fields = list(results[0].keys())

    for col, field in enumerate(
        fields,
        start=1
    ):

        cell = sheet.cell(
            row=1,
            column=col,
            value=field
        )

        cell.font = Font(
            bold=True
        )

    for row_number, record in enumerate(
        results,
        start=2
    ):

        for col, field in enumerate(
            fields,
            start=1
        ):

            sheet.cell(
                row=row_number,
                column=col,
                value=record[field]
            )

    sheet.freeze_panes = "A2"
    sheet.auto_filter.ref = sheet.dimensions

    # --------------------------------------------------------
    # CATEGORY SHEETS
    # --------------------------------------------------------

    categories = sorted(
        set(
            row["classification"]
            for row in results
        )
    )

    for category in categories:

        safe_name = category[:31]

        ws = workbook.create_sheet(
            title=safe_name
        )

        category_rows = [
            row
            for row in results
            if row["classification"]
            == category
        ]

        for col, field in enumerate(
            fields,
            start=1
        ):

            cell = ws.cell(
                row=1,
                column=col,
                value=field
            )

            cell.font = Font(
                bold=True
            )

        for row_number, record in enumerate(
            category_rows,
            start=2
        ):

            for col, field in enumerate(
                fields,
                start=1
            ):

                ws.cell(
                    row=row_number,
                    column=col,
                    value=record[field]
                )

        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions

    workbook.save(
        output_path
    )

    print(
        f"Saved Excel: {output_path}"
    )


# ============================================================
# SAVE CATEGORY TSV FILES
# ============================================================

def save_category_files(
    results,
    output_dir
):

    categories = sorted(
        set(
            row["classification"]
            for row in results
        )
    )

    for category in categories:

        category_rows = [
            row
            for row in results
            if row["classification"]
            == category
        ]

        filename = (
            category
            + ".tsv"
        )

        output_path = os.path.join(
            output_dir,
            filename
        )

        save_tsv(
            category_rows,
            output_path
        )


# ============================================================
# FASTA READER
# ============================================================

def read_fasta(path):

    sequences = {}

    current_id = None
    sequence_parts = []

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as file:

        for line in file:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                # Save previous sequence
                if current_id is not None:

                    sequences[current_id] = (
                        "".join(sequence_parts)
                    )

                header = line[1:]

                # First whitespace-separated token
                current_id = header.split()[0]

                sequence_parts = []

            else:

                sequence_parts.append(line)

    # Save last sequence
    if current_id is not None:

        sequences[current_id] = (
            "".join(sequence_parts)
        )

    return sequences


# ============================================================
# FASTA WRITER
# ============================================================

def write_fasta(
    sequence_ids,
    sequences,
    output_path
):

    written = 0
    missing = []

    with open(
        output_path,
        "w",
        encoding="utf-8"
    ) as file:

        for sequence_id in sequence_ids:

            if sequence_id not in sequences:

                missing.append(
                    sequence_id
                )

                continue

            file.write(
                ">"
                + sequence_id
                + "\n"
            )

            seq = sequences[
                sequence_id
            ]

            # 80 characters per line
            for i in range(
                0,
                len(seq),
                80
            ):

                file.write(
                    seq[i:i + 80]
                    + "\n"
                )

            written += 1

    print(
        f"\nFASTA sequences written: "
        f"{written:,}"
    )

    if missing:

        print(
            f"WARNING: {len(missing):,} "
            f"IDs were not found in FASTA."
        )

        missing_file = (
            output_path
            + ".missing.txt"
        )

        with open(
            missing_file,
            "w",
            encoding="utf-8"
        ) as file:

            for sequence_id in missing:

                file.write(
                    sequence_id
                    + "\n"
                )


# ============================================================
# CATEGORY SELECTION
# ============================================================

def ask_categories():

    print(
        "\n"
        + "=" * 70
    )

    print(
        "FASTA CATEGORY SELECTION"
    )

    print(
        "=" * 70
    )

    for number, category in CATEGORIES.items():

        print(
            f"{number} = {category}"
        )

    print(
        "\nYou can select multiple categories."
    )

    print(
        "Example: 1"
    )

    print(
        "Example: 1,5"
    )

    print(
        "Example: 2,3"
    )

    while True:

        answer = input(
            "\nEnter category numbers: "
        ).strip()

        if not answer:
            print(
                "Please enter at least one category."
            )
            continue

        try:

            numbers = [
                int(x.strip())
                for x in answer.split(",")
                if x.strip()
            ]

        except ValueError:

            print(
                "Invalid input. Use numbers such as 1,5."
            )

            continue

        invalid = [
            n
            for n in numbers
            if n not in CATEGORIES
        ]

        if invalid:

            print(
                "Invalid category number(s): "
                + ", ".join(
                    str(x)
                    for x in invalid
                )
            )

            continue

        # Remove duplicates while preserving order
        numbers = list(
            dict.fromkeys(numbers)
        )

        selected = [
            CATEGORIES[n]
            for n in numbers
        ]

        return numbers, selected


# ============================================================
# EXPORT SELECTED FASTA
# ============================================================

def export_selected_fasta(
    results,
    fasta_sequences,
    numbers,
    categories,
    output_dir
):

    selected_rows = [
        row
        for row in results
        if row["classification"]
        in categories
    ]

    selected_ids = [
        row["gene_id"]
        for row in selected_rows
    ]

    if not selected_ids:

        print(
            "\nNo sequences belong to the selected categories."
        )

        return

    if len(categories) == 1:

        category_text = categories[0]

    else:

        category_text = "_".join(
            categories
        )

    category_text = re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        category_text
    )

    output_path = os.path.join(
        output_dir,
        f"selected_{category_text}.fasta"
    )

    write_fasta(
        selected_ids,
        fasta_sequences,
        output_path
    )

    # Also save the corresponding screening table
    table_path = os.path.join(
        output_dir,
        f"selected_{category_text}.tsv"
    )

    save_tsv(
        selected_rows,
        table_path
    )

    print(
        "\nSelected categories:"
    )

    for category in categories:

        count = sum(
            1
            for row in results
            if row["classification"]
            == category
        )

        print(
            f"  {category}: {count:,}"
        )

    print(
        f"\nFASTA saved:"
    )

    print(
        output_path
    )


# ============================================================
# SUMMARY
# ============================================================

def print_summary(results):

    counts = Counter(
        row["classification"]
        for row in results
    )

    print(
        "\n"
        + "=" * 70
    )

    print(
        "CONTAMINATION SCREEN SUMMARY"
    )

    print(
        "=" * 70
    )

    print(
        f"Proteins screened: {len(results):,}"
    )

    for category in [
        "KEEP_ARTHROPOD",
        "POSSIBLE_BACTERIAL",
        "POSSIBLE_FUNGAL",
        "OTHER_EUKARYOTE",
        "AMBIGUOUS",
        "NO_QUALIFYING_HIT",
        "UNKNOWN_TAXONOMY",
    ]:

        count = counts.get(
            category,
            0
        )

        percentage = (
            count / len(results) * 100
            if results
            else 0
        )

        print(
            f"{category:25s} "
            f"{count:8,} "
            f"({percentage:6.2f}%)"
        )


# ============================================================
# MAIN
# ============================================================

def main():

    print(
        "\n"
        + "=" * 70
    )

    print(
        "METISA PLANA BLAST CONTAMINATION SCREEN"
    )

    print(
        "=" * 70
    )

    print(
        "\nThresholds:"
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

    print(
        f"  Taxonomic dominance >= "
        f"{DOMINANCE_THRESHOLD:.0%}"
    )

    # --------------------------------------------------------
    # INPUT BLAST
    # --------------------------------------------------------

    blast_file = expand_path(
        input(
            "\nPath to filtered BLAST TSV: "
        )
    )

    if not check_file_exists(
        blast_file,
        "BLAST TSV"
    ):

        return

    # --------------------------------------------------------
    # INPUT FASTA
    # --------------------------------------------------------

    fasta_file = expand_path(
        input(
            "\nPath to original protein FASTA: "
        )
    )

    if not check_file_exists(
        fasta_file,
        "Protein FASTA"
    ):

        return

    # --------------------------------------------------------
    # RANKED LINEAGE
    # --------------------------------------------------------

    rankedlineage_file = input(
        "\nPath to rankedlineage.dmp "
    ).strip()

    rankedlineage_file = expand_path(
        rankedlineage_file
    )

    if not check_file_exists(
        rankedlineage_file,
        "rankedlineage.dmp"
    ):

        return

    # --------------------------------------------------------
    # READ BLAST
    # --------------------------------------------------------

    rows = read_blast_file(
        blast_file
    )

    if not rows:

        print(
            "\nERROR: No BLAST rows found."
        )

        return

    # --------------------------------------------------------
    # GET TAXIDS FROM BLAST
    # --------------------------------------------------------

    print(
        "\nExtracting NCBI TaxIDs from BLAST stitle..."
    )

    required_taxids = set()

    for row in rows:

        taxid = extract_taxid(
            row["stitle"]
        )

        if taxid is not None:

            required_taxids.add(
                taxid
            )

    print(
        f"Unique TaxIDs found: "
        f"{len(required_taxids):,}"
    )

    # --------------------------------------------------------
    # LOAD TAXONOMY
    # --------------------------------------------------------

    taxonomy = load_required_taxonomy(
        rankedlineage_file,
        required_taxids
    )

    # --------------------------------------------------------
    # SCREEN
    # --------------------------------------------------------

    results = screen_proteins(
        rows,
        taxonomy
    )

    # --------------------------------------------------------
    # OUTPUT DIRECTORY
    # --------------------------------------------------------

    blast_dir = os.path.dirname(
        blast_file
    )

    output_dir = os.path.join(
        blast_dir,
        "contamination_screen"
    )

    os.makedirs(
        output_dir,
        exist_ok=True
    )

    fasta_output_dir = os.path.join(
        output_dir,
        "fasta_exports"
    )

    os.makedirs(
        fasta_output_dir,
        exist_ok=True
    )

    # --------------------------------------------------------
    # SAVE ALL TSV
    # --------------------------------------------------------

    all_tsv = os.path.join(
        output_dir,
        "screening_all.tsv"
    )

    save_tsv(
        results,
        all_tsv
    )

    # --------------------------------------------------------
    # SAVE CATEGORY TSVs
    # --------------------------------------------------------

    save_category_files(
        results,
        output_dir
    )

    # --------------------------------------------------------
    # SAVE EXCEL
    # --------------------------------------------------------

    excel_file = os.path.join(
        output_dir,
        "screening_all.xlsx"
    )

    save_excel(
        results,
        excel_file
    )

    # --------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------

    print_summary(
        results
    )

    # --------------------------------------------------------
    # READ ORIGINAL FASTA
    # --------------------------------------------------------

    print(
        "\nReading original protein FASTA..."
    )

    fasta_sequences = read_fasta(
        fasta_file
    )

    print(
        f"FASTA sequences loaded: "
        f"{len(fasta_sequences):,}"
    )

    # --------------------------------------------------------
    # ASK CATEGORY
    # --------------------------------------------------------

    numbers, categories = ask_categories()

    # --------------------------------------------------------
    # EXPORT FASTA
    # --------------------------------------------------------

    export_selected_fasta(
        results,
        fasta_sequences,
        numbers,
        categories,
        fasta_output_dir
    )

    # --------------------------------------------------------
    # DONE
    # --------------------------------------------------------

    print(
        "\n"
        + "=" * 70
    )

    print(
        "SCREENING COMPLETE"
    )

    print(
        "=" * 70
    )

    print(
        "\nResults folder:"
    )

    print(
        output_dir
    )

    print(
        "\nThe selected FASTA can now be used "
        "for downstream annotation such as InterProScan."
    )


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":
    main()
