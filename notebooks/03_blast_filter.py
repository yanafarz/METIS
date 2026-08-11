import os
import csv
import re
from collections import Counter

from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter


# ==========================================================
# CONFIGURATION
# ==========================================================

# BLAST hit-quality thresholds

MIN_IDENTITY = 30.0
MIN_COVERAGE = 70.0
MAX_EVALUE = 0.05


# ==========================================================
# EXPECTED BLAST COLUMNS
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


# ==========================================================
# PATH FUNCTIONS
# ==========================================================

def expand_path(path):

    return os.path.abspath(
        os.path.expanduser(
            path.strip().strip('"')
        )
    )


def check_file_exists(path):

    if not os.path.isfile(path):

        print("\nERROR: File not found:")
        print(path)

        return False

    return True


# ==========================================================
# ASK OUTPUT / RUN NAME
# ==========================================================

def ask_output_name(output_dir):

    print("\nOUTPUT / RUN NAME")

    print(
        "\nChoose a unique name for this screening run."
    )

    print(
        "Example:"
    )

    print(
        "metisa_blast_30id_70cov"
    )

    while True:

        name = input(
            "\nOutput/run name: "
        ).strip()

        if not name:

            print(
                "\nERROR: Output name cannot be empty."
            )

            continue

        # --------------------------------------------------
        # Remove characters unsafe for Windows filenames
        # --------------------------------------------------

        name = re.sub(
            r'[<>:"/\\|?*]',
            "_",
            name
        )

        # --------------------------------------------------
        # Remove accidental extension
        # --------------------------------------------------

        name = os.path.splitext(name)[0]

        if not name:

            print(
                "\nERROR: Invalid output name."
            )

            continue

        # --------------------------------------------------
        # Expected output files
        # --------------------------------------------------

        files = [

            os.path.join(
                output_dir,
                f"{name}.xlsx"
            ),

            os.path.join(
                output_dir,
                f"{name}_pass.fasta"
            ),

            os.path.join(
                output_dir,
                f"{name}_pass_ids.txt"
            ),

            os.path.join(
                output_dir,
                f"{name}_summary.txt"
            ),

        ]

        # --------------------------------------------------
        # Prevent overwriting an old run
        # --------------------------------------------------

        existing_files = [
            path
            for path in files
            if os.path.exists(path)
        ]

        if existing_files:

            print(
                "\nERROR: This output/run name already exists."
            )

            print(
                "\nExisting files:"
            )

            for path in existing_files:

                print(
                    f"  {os.path.basename(path)}"
                )

            print(
                "\nPlease choose a different output/run name."
            )

            continue

        return name


# ==========================================================
# READ BLAST RESULT
# ==========================================================

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

            # --------------------------------------------------
            # Skip empty lines
            # --------------------------------------------------

            if not row:
                continue

            if all(
                not cell.strip()
                for cell in row
            ):
                continue

            # --------------------------------------------------
            # Skip BLAST header
            # --------------------------------------------------

            if row[0].strip().lower() == "qseqid":
                continue

            # --------------------------------------------------
            # Check number of columns
            # --------------------------------------------------

            if len(row) < 10:

                print(
                    f"WARNING: Skipping line "
                    f"{line_number}: "
                    f"only {len(row)} columns found."
                )

                continue

            # --------------------------------------------------
            # Parse BLAST values
            # --------------------------------------------------

            try:

                record = {

                    "qseqid":
                        row[0].strip(),

                    "sseqid":
                        row[1].strip(),

                    "pident":
                        float(row[2]),

                    "length":
                        int(float(row[3])),

                    "qlen":
                        int(float(row[4])),

                    "slen":
                        int(float(row[5])),

                    "qcovs":
                        float(row[6]),

                    "evalue":
                        float(row[7]),

                    "bitscore":
                        float(row[8]),

                    "stitle":
                        row[9].strip(),

                }

                rows.append(record)

            except ValueError:

                print(
                    f"WARNING: Could not parse line "
                    f"{line_number}; skipping it."
                )

                continue

    print(
        f"BLAST hit rows read: "
        f"{len(rows):,}"
    )

    return rows


# ==========================================================
# READ PROTEIN FASTA
# ==========================================================

def read_fasta(path):

    print("\nReading protein FASTA...")

    sequences = {}

    current_id = None
    current_sequence = []

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as file:

        for line_number, line in enumerate(
            file,
            start=1
        ):

            line = line.strip()

            if not line:
                continue

            # --------------------------------------------------
            # FASTA header
            # --------------------------------------------------

            if line.startswith(">"):

                # Save previous sequence
                if current_id is not None:

                    sequences[current_id] = (
                        "".join(current_sequence)
                    )

                header = line[1:].strip()

                if not header:

                    print(
                        f"WARNING: Empty FASTA header "
                        f"at line {line_number}."
                    )

                    current_id = None
                    current_sequence = []

                    continue

                # FASTA ID = first word
                current_id = header.split()[0]

                current_sequence = []

            else:

                if current_id is None:
                    continue

                current_sequence.append(
                    line.replace(" ", "")
                )

    # ------------------------------------------------------
    # Save final sequence
    # ------------------------------------------------------

    if current_id is not None:

        sequences[current_id] = (
            "".join(current_sequence)
        )

    print(
        f"Protein sequences loaded: "
        f"{len(sequences):,}"
    )

    return sequences


# ==========================================================
# EXTRACT SPECIES
# ==========================================================

def extract_species(stitle):

    if not stitle:
        return "Unknown"

    # ------------------------------------------------------
    # UniProt OS= format
    # ------------------------------------------------------

    match = re.search(
        r"\bOS=([^=]+?)(?=\s+(?:OX|GN|PE|SV|CC|KW|GO)=|$)",
        stitle
    )

    if match:

        species = match.group(1).strip()

        if species:
            return species

    # ------------------------------------------------------
    # Alternative organism formats
    # ------------------------------------------------------

    match = re.search(
        r"\borganism[=:]\s*([^,;]+)",
        stitle,
        flags=re.IGNORECASE
    )

    if match:

        species = match.group(1).strip()

        if species:
            return species

    # ------------------------------------------------------
    # NCBI-style [species]
    # ------------------------------------------------------

    matches = re.findall(
        r"\[([^\[\]]+)\]",
        stitle
    )

    if matches:

        for value in reversed(matches):

            value = value.strip()

            if value:
                return value

    return "Unknown"


# ==========================================================
# EXTRACT NCBI TAXONOMY ID
# ==========================================================

def extract_taxid(stitle):

    if not stitle:
        return None

    match = re.search(
        r"\bOX=(\d+)",
        stitle
    )

    if match:
        return match.group(1)

    return None


# ==========================================================
# ASK FOR NCBI TAXONOMY FILE
# ==========================================================

def ask_taxonomy_path():

    print("\nSTEP 3: NCBI TAXONOMY")

    print(
        "\nYou need the NCBI rankedlineage.dmp file."
    )

    print(
        "Example:"
    )

    print(
        r"C:\BLAST\rankedlineage.dmp"
    )

    while True:

        taxonomy_file = expand_path(
            input(
                "\nPath to rankedlineage.dmp: "
            )
        )

        if not os.path.isfile(
            taxonomy_file
        ):

            print(
                "\nERROR: File not found:"
            )

            print(
                taxonomy_file
            )

            continue

        filename = os.path.basename(
            taxonomy_file
        ).lower()

        if filename != "rankedlineage.dmp":

            print(
                "\nWARNING: This does not appear "
                "to be rankedlineage.dmp."
            )

            answer = input(
                "\nUse this file anyway? (y/n): "
            ).strip().lower()

            if answer != "y":
                continue

        print(
            "\nNCBI taxonomy file detected:"
        )

        print(
            taxonomy_file
        )

        return taxonomy_file


# ==========================================================
# LOAD NCBI TAXONOMY
# ==========================================================

def load_ranked_lineage(taxonomy_file):

    print(
        "\nLoading NCBI ranked taxonomy..."
    )

    taxonomy = {}

    with open(
        taxonomy_file,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as file:

        for line in file:

            line = line.rstrip("\n")

            if not line:
                continue

            fields = [
                field.strip()
                for field in line.split("|")
            ]

            if len(fields) < 2:
                continue

            taxid = fields[0]

            taxonomy[taxid] = {

                "scientific_name":
                    fields[1]
                    if len(fields) > 1
                    else "",

                "genus":
                    fields[2]
                    if len(fields) > 2
                    else "",

                "family":
                    fields[3]
                    if len(fields) > 3
                    else "",

                "order":
                    fields[4]
                    if len(fields) > 4
                    else "",

                "class":
                    fields[5]
                    if len(fields) > 5
                    else "",

                "phylum":
                    fields[6]
                    if len(fields) > 6
                    else "",

                "kingdom":
                    fields[7]
                    if len(fields) > 7
                    else "",

                "superkingdom":
                    fields[8]
                    if len(fields) > 8
                    else "",

            }

    print(
        f"NCBI taxonomy records loaded: "
        f"{len(taxonomy):,}"
    )

    return taxonomy


# ==========================================================
# CLASSIFY TAXONOMY
# ==========================================================

def classify_taxonomy(taxonomy_record):

    if not taxonomy_record:
        return "Unknown"

    values = [

        taxonomy_record.get(
            "superkingdom",
            ""
        ),

        taxonomy_record.get(
            "kingdom",
            ""
        ),

        taxonomy_record.get(
            "phylum",
            ""
        ),

        taxonomy_record.get(
            "class",
            ""
        ),

        taxonomy_record.get(
            "order",
            ""
        ),

    ]

    text = " ".join(
        values
    ).lower()

    # ------------------------------------------------------
    # Bacteria
    # ------------------------------------------------------

    if (
        "bacteria" in text
        or
        "bacteri" in text
    ):

        return "Bacterial"

    # ------------------------------------------------------
    # Archaea
    # ------------------------------------------------------

    if (
        "archaea" in text
        or
        "archae" in text
    ):

        return "Archaeal"

    # ------------------------------------------------------
    # Fungi
    # ------------------------------------------------------

    if (
        "fungi" in text
        or
        "ascomycota" in text
        or
        "basidiomycota" in text
    ):

        return "Fungal"

    # ------------------------------------------------------
    # Arthropoda
    # ------------------------------------------------------

    if "arthropoda" in text:

        return "Arthropod"

    # ------------------------------------------------------
    # Plants
    # ------------------------------------------------------

    if (
        "viridiplantae" in text
        or
        "streptophyta" in text
        or
        "tracheophyta" in text
    ):

        return "Plant"

    # ------------------------------------------------------
    # Everything else
    # ------------------------------------------------------

    return "Other non-arthropod"


# ==========================================================
# BLAST QUALITY FILTER
# ==========================================================

def passes_quality_filter(hit):

    return (

        hit["pident"]
        >= MIN_IDENTITY

        and

        hit["qcovs"]
        >= MIN_COVERAGE

        and

        hit["evalue"]
        <= MAX_EVALUE

    )


# ==========================================================
# SELECT BEST QUALIFYING HIT
# ==========================================================

def select_best_hit(hits):

    qualifying = []

    for hit in hits:

        if passes_quality_filter(hit):

            qualifying.append(
                hit
            )

    if not qualifying:

        return None, []

    # ------------------------------------------------------
    # Ranking:
    #
    # 1. Lowest E-value
    # 2. Highest bitscore
    # 3. Highest identity
    # 4. Highest coverage
    # ------------------------------------------------------

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
# DECIDE ACTION
# ==========================================================

def decide_action(
    best_classification,
    arthropod_hits,
    bacterial_hits,
    fungal_hits,
    other_hits
):

    # ------------------------------------------------------
    # Unknown best classification
    # ------------------------------------------------------

    if best_classification == "Unknown":

        return "REVIEW_UNKNOWN"

    # ------------------------------------------------------
    # Bacterial dominance
    # ------------------------------------------------------

    if bacterial_hits > arthropod_hits:

        return "FLAG_POSSIBLE_BACTERIAL"

    # ------------------------------------------------------
    # Fungal dominance
    # ------------------------------------------------------

    if fungal_hits > arthropod_hits:

        return "FLAG_POSSIBLE_FUNGAL"

    # ------------------------------------------------------
    # Arthropod evidence
    # ------------------------------------------------------

    if arthropod_hits > 0:

        if (
            arthropod_hits
            >= bacterial_hits
            and
            arthropod_hits
            >= fungal_hits
            and
            arthropod_hits
            >= other_hits
        ):

            return "KEEP_ARTHROPOD_REVIEW"

    # ------------------------------------------------------
    # Mixed evidence
    # ------------------------------------------------------

    if (
        arthropod_hits > 0
        and
        other_hits > 0
    ):

        return "REVIEW_MIXED"

    # ------------------------------------------------------
    # Other non-arthropod
    # ------------------------------------------------------

    if other_hits > 0:

        return "REVIEW_NON_ARTHROPOD"

    return "REVIEW"


# ==========================================================
# SCREEN ALL QUERY SEQUENCES
# ==========================================================

def screen_contamination(
    rows,
    taxonomy
):

    print(
        "\nGrouping BLAST hits by query..."
    )

    grouped = {}

    # ------------------------------------------------------
    # Group hits
    # ------------------------------------------------------

    for row in rows:

        query = row["qseqid"]

        if query not in grouped:

            grouped[query] = []

        grouped[query].append(
            row
        )

    print(
        f"Unique query sequences: "
        f"{len(grouped):,}"
    )

    results = []

    print(
        "\nPerforming threshold filtering "
        "and taxonomy screening..."
    )

    # ------------------------------------------------------
    # Process every query
    # ------------------------------------------------------

    for query_id, hits in grouped.items():

        best_hit, qualifying_hits = (
            select_best_hit(hits)
        )

        # --------------------------------------------------
        # Count classifications
        # --------------------------------------------------

        arthropod_hits = 0
        bacterial_hits = 0
        fungal_hits = 0
        other_hits = 0
        unknown_hits = 0

        for hit in qualifying_hits:

            taxid = extract_taxid(
                hit["stitle"]
            )

            tax_record = (
                taxonomy.get(taxid)
                if taxid
                else None
            )

            classification = (
                classify_taxonomy(
                    tax_record
                )
            )

            if classification == "Arthropod":

                arthropod_hits += 1

            elif classification == "Bacterial":

                bacterial_hits += 1

            elif classification == "Fungal":

                fungal_hits += 1

            elif classification == "Unknown":

                unknown_hits += 1

            else:

                other_hits += 1

        # --------------------------------------------------
        # No qualifying hit
        # --------------------------------------------------

        if best_hit is None:

            results.append({

                "Gene_ID":
                    query_id,

                "Best_Hit":
                    "",

                "Best_Hit_Species":
                    "",

                "Best_Hit_Classification":
                    "No qualifying hit",

                "Best_Hit_Source":
                    "",

                "Best_Hit_Identity":
                    "",

                "Best_Hit_Coverage":
                    "",

                "Best_Hit_E_value":
                    "",

                "Total_Hits":
                    len(hits),

                "Qualifying_Hits":
                    0,

                "Arthropod_Hits":
                    0,

                "Bacterial_Hits":
                    0,

                "Fungal_Hits":
                    0,

                "Other_Non_Arthropod_Hits":
                    0,

                "Unknown_Hits":
                    0,

                "Action":
                    "REVIEW_NO_QUALIFYING_HIT",

            })

            continue

        # --------------------------------------------------
        # Best hit information
        # --------------------------------------------------

        species = extract_species(
            best_hit["stitle"]
        )

        taxid = extract_taxid(
            best_hit["stitle"]
        )

        tax_record = (
            taxonomy.get(taxid)
            if taxid
            else None
        )

        best_classification = (
            classify_taxonomy(
                tax_record
            )
        )

        # --------------------------------------------------
        # Determine source
        # --------------------------------------------------

        sseqid = best_hit["sseqid"]

        if sseqid.startswith("sp|"):

            source = "Swiss-Prot"

        elif sseqid.startswith("tr|"):

            source = "TrEMBL"

        elif sseqid.startswith("iso|"):

            source = "Isoform"

        else:

            source = "Other"

        # --------------------------------------------------
        # Action
        # --------------------------------------------------

        action = decide_action(

            best_classification,

            arthropod_hits,

            bacterial_hits,

            fungal_hits,

            other_hits

        )

        # --------------------------------------------------
        # Save result
        # --------------------------------------------------

        results.append({

            "Gene_ID":
                query_id,

            "Best_Hit":
                best_hit["sseqid"],

            "Best_Hit_Species":
                species,

            "Best_Hit_Classification":
                best_classification,

            "Best_Hit_Source":
                source,

            "Best_Hit_Identity":
                best_hit["pident"],

            "Best_Hit_Coverage":
                best_hit["qcovs"],

            "Best_Hit_E_value":
                best_hit["evalue"],

            "Total_Hits":
                len(hits),

            "Qualifying_Hits":
                len(qualifying_hits),

            "Arthropod_Hits":
                arthropod_hits,

            "Bacterial_Hits":
                bacterial_hits,

            "Fungal_Hits":
                fungal_hits,

            "Other_Non_Arthropod_Hits":
                other_hits,

            "Unknown_Hits":
                unknown_hits,

            "Action":
                action,

        })

    return results


# ==========================================================
# GET PASSING GENE IDS
# ==========================================================

def get_passing_gene_ids(results):

    return [
        row["Gene_ID"]
        for row in results
        if row["Qualifying_Hits"] > 0
    ]


# ==========================================================
# SAVE PASSING GENE IDS
# ==========================================================

def save_passing_ids(
    passing_ids,
    output_path
):

    print(
        "\nSaving passing gene IDs..."
    )

    with open(
        output_path,
        "w",
        encoding="utf-8"
    ) as file:

        for gene_id in passing_ids:

            file.write(
                gene_id + "\n"
            )

    print(
        f"Passing gene IDs saved: "
        f"{len(passing_ids):,}"
    )

    print(
        output_path
    )


# ==========================================================
# SAVE PASSING FASTA
# ==========================================================

def save_passing_fasta(
    passing_ids,
    sequences,
    output_path
):

    print(
        "\nCreating FASTA containing "
        "threshold-passing proteins..."
    )

    written = 0
    missing = []

    with open(
        output_path,
        "w",
        encoding="utf-8"
    ) as file:

        for gene_id in passing_ids:

            sequence = sequences.get(
                gene_id
            )

            if sequence is None:

                missing.append(
                    gene_id
                )

                continue

            file.write(
                f">{gene_id}\n"
            )

            # --------------------------------------------------
            # Write sequence in 80-character lines
            # --------------------------------------------------

            for start in range(
                0,
                len(sequence),
                80
            ):

                file.write(
                    sequence[
                        start:start + 80
                    ]
                    + "\n"
                )

            written += 1

    print(
        f"Passing FASTA sequences written: "
        f"{written:,}"
    )

    # ------------------------------------------------------
    # Save missing IDs if any
    # ------------------------------------------------------

    if missing:

        print(
            f"WARNING: {len(missing):,} "
            "passing IDs were not found "
            "in the FASTA."
        )

        base, ext = os.path.splitext(
            output_path
        )

        missing_file = (
            base + "_missing_ids.txt"
        )

        with open(
            missing_file,
            "w",
            encoding="utf-8"
        ) as file:

            for gene_id in missing:

                file.write(
                    gene_id + "\n"
                )

        print(
            "Missing IDs saved:"
        )

        print(
            missing_file
        )

    print(
        "\nPassing FASTA saved:"
    )

    print(
        output_path
    )


# ==========================================================
# SAVE EXCEL
# ==========================================================

def save_screening_excel(
    results,
    output_path
):

    print(
        "\nCreating Excel screening file..."
    )

    fieldnames = [

        "Gene_ID",

        "Best_Hit",

        "Best_Hit_Species",

        "Best_Hit_Classification",

        "Best_Hit_Source",

        "Best_Hit_Identity",

        "Best_Hit_Coverage",

        "Best_Hit_E_value",

        "Total_Hits",

        "Qualifying_Hits",

        "Arthropod_Hits",

        "Bacterial_Hits",

        "Fungal_Hits",

        "Other_Non_Arthropod_Hits",

        "Unknown_Hits",

        "Action",

    ]

    # ------------------------------------------------------
    # Create workbook
    # ------------------------------------------------------

    workbook = Workbook()

    worksheet = workbook.active

    worksheet.title = "Taxonomy Screening"

    # ------------------------------------------------------
    # Header
    # ------------------------------------------------------

    for column_number, fieldname in enumerate(
        fieldnames,
        start=1
    ):

        cell = worksheet.cell(
            row=1,
            column=column_number,
            value=fieldname
        )

        cell.font = Font(
            bold=True
        )

        cell.alignment = Alignment(
            horizontal="center",
            vertical="center"
        )

    # ------------------------------------------------------
    # Data
    # ------------------------------------------------------

    for row_number, result in enumerate(
        results,
        start=2
    ):

        for column_number, fieldname in enumerate(
            fieldnames,
            start=1
        ):

            worksheet.cell(
                row=row_number,
                column=column_number,
                value=result.get(
                    fieldname,
                    ""
                )
            )

    worksheet.freeze_panes = "A2"

    # ------------------------------------------------------
    # Excel filter
    # ------------------------------------------------------

    last_column = get_column_letter(
        len(fieldnames)
    )

    last_row = len(results) + 1

    worksheet.auto_filter.ref = (
        f"A1:{last_column}{last_row}"
    )

    # ------------------------------------------------------
    # Column widths
    # ------------------------------------------------------

    for column_number, fieldname in enumerate(
        fieldnames,
        start=1
    ):

        column_letter = get_column_letter(
            column_number
        )

        if fieldname == "Best_Hit":

            width = 35

        elif fieldname == "Best_Hit_Species":

            width = 25

        elif fieldname == "Action":

            width = 32

        elif fieldname == "Gene_ID":

            width = 18

        else:

            width = max(
                len(fieldname) + 2,
                15
            )

        worksheet.column_dimensions[
            column_letter
        ].width = width

    # ------------------------------------------------------
    # Number formats
    # ------------------------------------------------------

    identity_column = (
        fieldnames.index(
            "Best_Hit_Identity"
        ) + 1
    )

    coverage_column = (
        fieldnames.index(
            "Best_Hit_Coverage"
        ) + 1
    )

    evalue_column = (
        fieldnames.index(
            "Best_Hit_E_value"
        ) + 1
    )

    for row_number in range(
        2,
        last_row + 1
    ):

        worksheet.cell(
            row=row_number,
            column=identity_column
        ).number_format = "0.000"

        worksheet.cell(
            row=row_number,
            column=coverage_column
        ).number_format = "0.00"

        worksheet.cell(
            row=row_number,
            column=evalue_column
        ).number_format = "0.00E+00"

    # ------------------------------------------------------
    # Save
    # ------------------------------------------------------

    workbook.save(
        output_path
    )

    print(
        "\nExcel screening result saved:"
    )

    print(
        output_path
    )


# ==========================================================
# SAVE SUMMARY TXT
# ==========================================================

def save_summary(
    results,
    output_path
):

    classification_counts = Counter(

        row[
            "Best_Hit_Classification"
        ]

        for row in results

    )

    action_counts = Counter(

        row[
            "Action"
        ]

        for row in results

    )

    total = len(results)

    passing = sum(
        1
        for row in results
        if row["Qualifying_Hits"] > 0
    )

    no_passing = (
        total - passing
    )

    with open(
        output_path,
        "w",
        encoding="utf-8"
    ) as file:

        file.write(
            "METISA PLANA TAXONOMY SCREENING SUMMARY\n"
        )

        file.write(
            "=" * 60 + "\n\n"
        )

        file.write(
            "BLAST thresholds:\n"
        )

        file.write(
            f"Identity >= {MIN_IDENTITY}%\n"
        )

        file.write(
            f"Coverage >= {MIN_COVERAGE}%\n"
        )

        file.write(
            f"E-value <= {MAX_EVALUE}\n\n"
        )

        file.write(
            f"Sequences screened: "
            f"{total:,}\n"
        )

        file.write(
            f"Sequences passing threshold: "
            f"{passing:,}\n"
        )

        file.write(
            f"Sequences without qualifying hit: "
            f"{no_passing:,}\n\n"
        )

        file.write(
            "Best-hit classifications:\n"
        )

        for classification, count in (
            classification_counts.items()
        ):

            file.write(
                f"  {classification}: "
                f"{count:,}\n"
            )

        file.write(
            "\nActions:\n"
        )

        for action, count in (
            action_counts.items()
        ):

            file.write(
                f"  {action}: "
                f"{count:,}\n"
            )

    print(
        "\nScreening summary saved:"
    )

    print(
        output_path
    )


# ==========================================================
# TERMINAL SUMMARY
# ==========================================================

def print_summary(results):

    classification_counts = Counter(

        row[
            "Best_Hit_Classification"
        ]

        for row in results

    )

    action_counts = Counter(

        row[
            "Action"
        ]

        for row in results

    )

    total = len(results)

    passing = sum(
        1
        for row in results
        if row["Qualifying_Hits"] > 0
    )

    print(
        "\n"
        + "=" * 65
    )

    print(
        "TAXONOMY SCREEN SUMMARY"
    )

    print(
        "=" * 65
    )

    print(
        f"Sequences screened: "
        f"{total:,}"
    )

    print(
        f"Passed BLAST threshold: "
        f"{passing:,}"
    )

    print(
        f"Did not pass threshold: "
        f"{total - passing:,}"
    )

    print(
        f"\nArthropod best hits: "
        f"{classification_counts.get('Arthropod', 0):,}"
    )

    print(
        f"Bacterial best hits: "
        f"{classification_counts.get('Bacterial', 0):,}"
    )

    print(
        f"Fungal best hits: "
        f"{classification_counts.get('Fungal', 0):,}"
    )

    print(
        f"Archaeal best hits: "
        f"{classification_counts.get('Archaeal', 0):,}"
    )

    print(
        f"Plant best hits: "
        f"{classification_counts.get('Plant', 0):,}"
    )

    print(
        f"Other non-arthropod: "
        f"{classification_counts.get('Other non-arthropod', 0):,}"
    )

    print(
        f"Unknown: "
        f"{classification_counts.get('Unknown', 0):,}"
    )

    print(
        f"No qualifying hit: "
        f"{classification_counts.get('No qualifying hit', 0):,}"
    )

    print(
        "\nActions:"
    )

    for action, count in action_counts.items():

        print(
            f"  {action}: {count:,}"
        )


# ==========================================================
# MAIN PROGRAM
# ==========================================================

def main():

    print(
        "\n"
        + "=" * 65
    )

    print(
        "       METISA PLANA TAXONOMY SCREEN"
    )

    print(
        "=" * 65
    )

    # ======================================================
    # STEP 1: BLAST RESULT
    # ======================================================

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

    # ======================================================
    # STEP 2: ORIGINAL PROTEIN FASTA
    # ======================================================

    print(
        "\nSTEP 2: ORIGINAL PROTEIN FASTA"
    )

    fasta_file = expand_path(
        input(
            "Path to original protein FASTA: "
        )
    )

    if not check_file_exists(
        fasta_file
    ):

        return

    # ======================================================
    # STEP 3: NCBI TAXONOMY
    # ======================================================

    taxonomy_file = ask_taxonomy_path()

    if not taxonomy_file:

        return

    # ======================================================
    # STEP 4: LOAD TAXONOMY
    # ======================================================

    print(
        "\nSTEP 4: LOADING TAXONOMY"
    )

    taxonomy = load_ranked_lineage(
        taxonomy_file
    )

    if not taxonomy:

        print(
            "\nERROR: No taxonomy records loaded."
        )

        return

    # ======================================================
    # STEP 5: READ BLAST
    # ======================================================

    print(
        "\nSTEP 5: READING BLAST"
    )

    rows = read_blast_file(
        blast_file
    )

    if not rows:

        print(
            "\nERROR: No usable BLAST rows found."
        )

        return

    # ======================================================
    # STEP 6: READ FASTA
    # ======================================================

    print(
        "\nSTEP 6: READING FASTA"
    )

    sequences = read_fasta(
        fasta_file
    )

    if not sequences:

        print(
            "\nERROR: No protein sequences found."
        )

        return

    # ======================================================
    # STEP 7: SCREEN
    # ======================================================

    print(
        "\nSTEP 7: THRESHOLD + TAXONOMY SCREENING"
    )

    print(
        "\nUsing thresholds:"
    )

    print(
        f"  Identity >= {MIN_IDENTITY}%"
    )

    print(
        f"  Coverage >= {MIN_COVERAGE}%"
    )

    print(
        f"  E-value <= {MAX_EVALUE}"
    )

    results = screen_contamination(
        rows,
        taxonomy
    )

    if not results:

        print(
            "\nERROR: No screening results."
        )

        return

    # ======================================================
    # OUTPUT DIRECTORY
    # ======================================================

    output_dir = os.path.join(

        os.path.dirname(
            blast_file
        ),

        "taxonomy_screen"

    )

    os.makedirs(
        output_dir,
        exist_ok=True
    )

    # ======================================================
    # ASK FOR UNIQUE OUTPUT NAME
    # ======================================================

    output_name = ask_output_name(
        output_dir
    )

    # ======================================================
    # GET PASSING IDS
    # ======================================================

    passing_ids = get_passing_gene_ids(
        results
    )

    print(
        "\n"
        + "=" * 65
    )

    print(
        f"Genes passing BLAST threshold: "
        f"{len(passing_ids):,}"
    )

    print(
        "=" * 65
    )

    # ======================================================
    # OUTPUT FILES
    # ======================================================

    excel_file = os.path.join(
        output_dir,
        f"{output_name}.xlsx"
    )

    fasta_output = os.path.join(
        output_dir,
        f"{output_name}_pass.fasta"
    )

    ids_output = os.path.join(
        output_dir,
        f"{output_name}_pass_ids.txt"
    )

    summary_output = os.path.join(
        output_dir,
        f"{output_name}_summary.txt"
    )

    # ======================================================
    # SAVE EXCEL
    # ======================================================

    save_screening_excel(
        results,
        excel_file
    )

    # ======================================================
    # SAVE PASSING IDS
    # ======================================================

    save_passing_ids(
        passing_ids,
        ids_output
    )

    # ======================================================
    # SAVE PASSING FASTA
    # ======================================================

    save_passing_fasta(
        passing_ids,
        sequences,
        fasta_output
    )

    # ======================================================
    # SAVE SUMMARY
    # ======================================================

    save_summary(
        results,
        summary_output
    )

    # ======================================================
    # TERMINAL SUMMARY
    # ======================================================

    print_summary(
        results
    )

    # ======================================================
    # FINISH
    # ======================================================

    print(
        "\n"
        + "=" * 65
    )

    print(
        "DONE"
    )

    print(
        "=" * 65
    )

    print(
        "\nOutput folder:"
    )

    print(
        output_dir
    )

    print(
        "\nOutput files:"
    )

    print(
        f"  Excel:   {excel_file}"
    )

    print(
        f"  FASTA:   {fasta_output}"
    )

    print(
        f"  IDs:     {ids_output}"
    )

    print(
        f"  Summary: {summary_output}"
    )


# ==========================================================
# RUN
# ==========================================================

if __name__ == "__main__":

    main()
