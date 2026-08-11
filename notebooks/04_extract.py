
import os
import csv


# ==========================================================
# METISA PLANA FASTA EXTRACTION
# ==========================================================
#
# Purpose:
#
# 1. Read the ORIGINAL Metisa plana protein FASTA
# 2. Read taxonomy screening results (.xlsx or .csv)
# 3. Let the user choose an Action
# 4. Extract matching Gene_ID sequences
# 5. Write them into a new FASTA
# 6. Verify that all selected sequences were recovered
#
# Example workflow:
#
# original_proteins.fasta
#        +
# contamination_screening.xlsx
#        ↓
# KEEP_ARTHROPOD_REVIEW
#        ↓
# metisa_arthropod_candidates.fasta
#        ↓
# InterProScan
#
# ==========================================================


# ==========================================================
# PATH FUNCTIONS
# ==========================================================

def expand_path(path):
    """Convert user-entered path into an absolute path."""

    return os.path.abspath(
        os.path.expanduser(
            path.strip().strip('"')
        )
    )


def check_file_exists(path, description):

    if not os.path.isfile(path):

        print("\nERROR: File not found:")
        print(path)

        print(
            f"\nPlease check the path to the "
            f"{description}."
        )

        return False

    return True


# ==========================================================
# READ ORIGINAL FASTA
# ==========================================================

def read_fasta(fasta_path):

    print("\nReading original protein FASTA...")

    sequences = {}

    current_id = None
    current_header = None
    current_sequence = []

    record_count = 0

    with open(
        fasta_path,
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

                    sequence = "".join(
                        current_sequence
                    )

                    sequences[current_id] = {
                        "header": current_header,
                        "sequence": sequence
                    }

                    record_count += 1


                current_header = line[1:].strip()

                if not current_header:

                    print(
                        f"WARNING: Empty FASTA header "
                        f"at line {line_number}."
                    )

                    current_id = None
                    current_sequence = []

                    continue


                # FASTA ID = first item before whitespace
                current_id = (
                    current_header.split()[0]
                )

                current_sequence = []


            else:

                # Sequence before first FASTA header
                if current_id is None:

                    print(
                        f"WARNING: Sequence data found "
                        f"before FASTA header at line "
                        f"{line_number}."
                    )

                    continue


                current_sequence.append(
                    line
                )


    # ------------------------------------------------------
    # Save final sequence
    # ------------------------------------------------------

    if current_id is not None:

        sequence = "".join(
            current_sequence
        )

        sequences[current_id] = {
            "header": current_header,
            "sequence": sequence
        }

        record_count += 1


    print(
        f"FASTA records read: "
        f"{record_count:,}"
    )

    print(
        f"Unique FASTA IDs: "
        f"{len(sequences):,}"
    )


    if not sequences:

        print(
            "\nERROR: No FASTA sequences were found."
        )

        print(
            "Make sure you entered the ORIGINAL "
            "protein FASTA, not a BLAST TSV."
        )

        return {}


    return sequences


# ==========================================================
# READ SCREENING CSV
# ==========================================================

def read_csv_screening(csv_path):

    print(
        "\nReading screening CSV..."
    )

    rows = []

    with open(
        csv_path,
        "r",
        encoding="utf-8-sig",
        errors="replace",
        newline=""
    ) as file:

        reader = csv.DictReader(file)

        if reader.fieldnames is None:

            print(
                "\nERROR: CSV has no header."
            )

            return []


        # Remove whitespace from column names
        reader.fieldnames = [
            field.strip()
            for field in reader.fieldnames
        ]


        if "Gene_ID" not in reader.fieldnames:

            print(
                "\nERROR: 'Gene_ID' column "
                "was not found."
            )

            print(
                "\nColumns found:"
            )

            for field in reader.fieldnames:

                print(
                    f"  {field}"
                )

            return []


        if "Action" not in reader.fieldnames:

            print(
                "\nERROR: 'Action' column "
                "was not found."
            )

            return []


        for row_number, row in enumerate(
            reader,
            start=2
        ):

            gene_id = (
                row.get("Gene_ID", "")
                or ""
            ).strip()

            action = (
                row.get("Action", "")
                or ""
            ).strip()


            if not gene_id:

                print(
                    f"WARNING: Row {row_number} "
                    "has empty Gene_ID. "
                    "Skipping."
                )

                continue


            rows.append({
                "Gene_ID": gene_id,
                "Action": action
            })


    print(
        f"Screening rows loaded: "
        f"{len(rows):,}"
    )


    return rows


# ==========================================================
# READ SCREENING XLSX
# ==========================================================

def read_xlsx_screening(xlsx_path):

    print(
        "\nReading screening Excel file..."
    )


    # ------------------------------------------------------
    # Import pandas only when XLSX is used
    # ------------------------------------------------------

    try:

        import pandas as pd

    except ImportError:

        print(
            "\nERROR: pandas is required "
            "to read Excel files."
        )

        print(
            "\nInstall it with:"
        )

        print(
            "python -m pip install pandas openpyxl"
        )

        return []


    try:

        df = pd.read_excel(
            xlsx_path
        )

    except Exception as error:

        print(
            "\nERROR: Could not read Excel file."
        )

        print(
            f"Reason: {error}"
        )

        return []


    # ------------------------------------------------------
    # Clean column names
    # ------------------------------------------------------

    df.columns = [
        str(column).strip()
        for column in df.columns
    ]


    # ------------------------------------------------------
    # Check required columns
    # ------------------------------------------------------

    if "Gene_ID" not in df.columns:

        print(
            "\nERROR: 'Gene_ID' column "
            "was not found."
        )

        print(
            "\nColumns found:"
        )

        for column in df.columns:

            print(
                f"  {column}"
            )

        return []


    if "Action" not in df.columns:

        print(
            "\nERROR: 'Action' column "
            "was not found."
        )

        return []


    # ------------------------------------------------------
    # Convert to simple records
    # ------------------------------------------------------

    rows = []

    for row_number, row in df.iterrows():

        gene_id = str(
            row["Gene_ID"]
        ).strip()


        action = str(
            row["Action"]
        ).strip()


        if (
            not gene_id
            or
            gene_id.lower() == "nan"
        ):

            print(
                f"WARNING: Excel row "
                f"{row_number + 2} has "
                "empty Gene_ID. Skipping."
            )

            continue


        if (
            action.lower() == "nan"
        ):

            action = ""


        rows.append({
            "Gene_ID": gene_id,
            "Action": action
        })


    print(
        f"Screening rows loaded: "
        f"{len(rows):,}"
    )


    return rows


# ==========================================================
# READ SCREENING TABLE AUTOMATICALLY
# ==========================================================

def read_screening_file(path):

    extension = (
        os.path.splitext(path)[1]
        .lower()
    )


    if extension == ".csv":

        return read_csv_screening(
            path
        )


    if extension in [
        ".xlsx",
        ".xls"
    ]:

        return read_xlsx_screening(
            path
        )


    print(
        "\nERROR: Unsupported screening "
        "file type."
    )

    print(
        "\nAccepted formats:"
    )

    print(
        "  .xlsx"
    )

    print(
        "  .xls"
    )

    print(
        "  .csv"
    )

    return []


# ==========================================================
# SHOW ACTIONS
# ==========================================================

def choose_action(rows):

    action_counts = {}

    for row in rows:

        action = row["Action"]

        if not action:
            continue

        action_counts[action] = (
            action_counts.get(
                action,
                0
            ) + 1
        )


    actions = sorted(
        action_counts.keys()
    )


    if not actions:

        print(
            "\nERROR: No Action values "
            "were found."
        )

        return None


    print(
        "\nAvailable Actions:"
    )

    print(
        "-" * 70
    )


    for number, action in enumerate(
        actions,
        start=1
    ):

        print(
            f"{number}. "
            f"{action} "
            f"({action_counts[action]:,} sequences)"
        )


    print(
        "\nSelect the Action you want "
        "to extract."
    )


    while True:

        choice = input(
            "\nEnter action number or exact name: "
        ).strip()


        # --------------------------------------------------
        # Number
        # --------------------------------------------------

        if choice.isdigit():

            number = int(choice)

            if (
                1 <= number <= len(actions)
            ):

                return actions[
                    number - 1
                ]


            print(
                "\nERROR: Invalid number."
            )

            continue


        # --------------------------------------------------
        # Exact / case-insensitive
        # --------------------------------------------------

        for action in actions:

            if action.lower() == choice.lower():

                return action


        print(
            "\nERROR: Action not found."
        )


# ==========================================================
# SELECT GENE IDs
# ==========================================================

def select_genes(
    rows,
    selected_action
):

    selected = []

    for row in rows:

        if (
            row["Action"]
            ==
            selected_action
        ):

            selected.append(
                row["Gene_ID"]
            )


    # Remove duplicate Gene IDs
    unique_genes = list(
        dict.fromkeys(
            selected
        )
    )


    print(
        "\nSelection summary:"
    )

    print(
        f"Action: "
        f"{selected_action}"
    )

    print(
        f"Rows matching Action: "
        f"{len(selected):,}"
    )

    print(
        f"Unique Gene IDs: "
        f"{len(unique_genes):,}"
    )


    if (
        len(selected)
        !=
        len(unique_genes)
    ):

        print(
            f"Duplicate Gene IDs removed: "
            f"{len(selected) - len(unique_genes):,}"
        )


    return unique_genes


# ==========================================================
# MATCH GENES TO ORIGINAL FASTA
# ==========================================================

def match_genes_to_fasta(
    gene_ids,
    sequences
):

    found = []
    missing = []


    for gene_id in gene_ids:

        if gene_id in sequences:

            found.append(
                gene_id
            )

        else:

            missing.append(
                gene_id
            )


    return found, missing


# ==========================================================
# WRITE FASTA
# ==========================================================

def write_fasta(
    gene_ids,
    sequences,
    output_path
):

    print(
        "\nWriting selected protein sequences..."
    )


    with open(
        output_path,
        "w",
        encoding="utf-8"
    ) as file:

        for gene_id in gene_ids:

            record = sequences[
                gene_id
            ]


            # Keep original header
            file.write(
                f">{record['header']}\n"
            )


            sequence = record[
                "sequence"
            ]


            # Wrap sequence to 60 characters
            for start in range(
                0,
                len(sequence),
                60
            ):

                file.write(
                    sequence[
                        start:start + 60
                    ]
                )

                file.write("\n")


    print(
        "\nFASTA saved:"
    )

    print(
        output_path
    )


# ==========================================================
# SAVE MISSING IDs
# ==========================================================

def save_missing_ids(
    missing,
    output_path
):

    if not missing:
        return


    missing_path = (
        os.path.splitext(
            output_path
        )[0]
        +
        "_missing.txt"
    )


    with open(
        missing_path,
        "w",
        encoding="utf-8"
    ) as file:

        for gene_id in missing:

            file.write(
                gene_id + "\n"
            )


    print(
        "\nMissing Gene IDs saved:"
    )

    print(
        missing_path
    )


# ==========================================================
# VERIFY FASTA
# ==========================================================

def count_fasta_sequences(
    fasta_path
):

    count = 0

    with open(
        fasta_path,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as file:

        for line in file:

            if line.startswith(">"):

                count += 1


    return count


def verify_output(
    output_path,
    expected
):

    actual = count_fasta_sequences(
        output_path
    )


    print(
        "\nOutput verification:"
    )

    print(
        f"Expected sequences: "
        f"{expected:,}"
    )

    print(
        f"Output sequences:   "
        f"{actual:,}"
    )


    if actual == expected:

        print(
            "\nSUCCESS:"
        )

        print(
            "The number of sequences "
            "matches exactly."
        )

        return True


    print(
        "\nWARNING:"
    )

    print(
        "The expected and actual "
        "sequence counts do not match."
    )

    return False


# ==========================================================
# MAIN
# ==========================================================

def main():

    print(
        "\n"
        + "=" * 70
    )

    print(
        "       METISA PLANA FASTA EXTRACTION"
    )

    print(
        "=" * 70
    )


    # ======================================================
    # STEP 1 — ORIGINAL FASTA
    # ======================================================

    print(
        "\nSTEP 1: ORIGINAL PROTEIN FASTA"
    )

    print(
        "\nIMPORTANT:"
    )

    print(
        "Enter the ORIGINAL Metisa plana "
        "protein FASTA."
    )

    print(
        "Do NOT enter your BLAST TSV."
    )


    while True:

        fasta_path = expand_path(
            input(
                "\nPath to original protein FASTA: "
            )
        )


        if not check_file_exists(
            fasta_path,
            "original protein FASTA"
        ):

            continue


        extension = (
            os.path.splitext(
                fasta_path
            )[1]
            .lower()
        )


        if extension not in [
            ".fasta",
            ".faa",
            ".fa"
        ]:

            print(
                "\nWARNING:"
            )

            print(
                "This does not have a normal "
                "FASTA extension."
            )

            answer = input(
                "Use it anyway? (y/n): "
            ).strip().lower()


            if answer != "y":

                continue


        break


    # ======================================================
    # STEP 2 — SCREENING TABLE
    # ======================================================

    print(
        "\nSTEP 2: TAXONOMY SCREENING TABLE"
    )

    print(
        "\nAccepted formats:"
    )

    print(
        "  Excel: .xlsx / .xls"
    )

    print(
        "  CSV:   .csv"
    )


    while True:

        screening_path = expand_path(
            input(
                "\nPath to screening Excel/CSV: "
            )
        )


        if check_file_exists(
            screening_path,
            "screening Excel/CSV"
        ):

            extension = (
                os.path.splitext(
                    screening_path
                )[1]
                .lower()
            )


            if extension in [
                ".xlsx",
                ".xls",
                ".csv"
            ]:

                break


        print(
            "\nERROR: Please provide "
            "an .xlsx, .xls, or .csv file."
        )


    # ======================================================
    # STEP 3 — LOAD ORIGINAL FASTA
    # ======================================================

    print(
        "\nSTEP 3: LOADING ORIGINAL FASTA"
    )


    sequences = read_fasta(
        fasta_path
    )


    if not sequences:

        print(
            "\nExtraction stopped."
        )

        return


    # ======================================================
    # STEP 4 — LOAD SCREENING TABLE
    # ======================================================

    print(
        "\nSTEP 4: LOADING SCREENING TABLE"
    )


    rows = read_screening_file(
        screening_path
    )


    if not rows:

        print(
            "\nERROR: No usable screening "
            "rows were found."
        )

        return


    # ======================================================
    # STEP 5 — CHOOSE ACTION
    # ======================================================

    print(
        "\nSTEP 5: CHOOSE WHICH SEQUENCES "
        "TO EXTRACT"
    )


    selected_action = choose_action(
        rows
    )


    if not selected_action:

        return


    # ======================================================
    # STEP 6 — SELECT GENE IDs
    # ======================================================

    selected_genes = select_genes(
        rows,
        selected_action
    )


    if not selected_genes:

        print(
            "\nERROR: No Gene IDs matched "
            "the selected Action."
        )

        return


    # ======================================================
    # STEP 7 — MATCH AGAINST ORIGINAL FASTA
    # ======================================================

    print(
        "\nSTEP 6: MATCHING AGAINST "
        "ORIGINAL FASTA"
    )


    found, missing = match_genes_to_fasta(
        selected_genes,
        sequences
    )


    print(
        f"\nSelected Gene IDs:      "
        f"{len(selected_genes):,}"
    )

    print(
        f"Found in original FASTA: "
        f"{len(found):,}"
    )

    print(
        f"Missing from FASTA:     "
        f"{len(missing):,}"
    )


    # ======================================================
    # SAFETY CHECK
    # ======================================================

    if missing:

        print(
            "\nWARNING:"
        )

        print(
            "Some Gene IDs from the screening "
            "table were NOT found in the "
            "original FASTA."
        )

        print(
            "Those sequences will NOT be "
            "invented or substituted."
        )


        save_missing_ids(
            missing,
            os.path.join(
                os.path.dirname(
                    screening_path
                ),
                "extraction_missing.fasta"
            )
        )


        answer = input(
            "\nContinue using only sequences "
            "that were found? (y/n): "
        ).strip().lower()


        if answer != "y":

            print(
                "\nExtraction cancelled."
            )

            return


    if not found:

        print(
            "\nERROR: No selected Gene IDs "
            "were found in the original FASTA."
        )

        return


    # ======================================================
    # STEP 8 — OUTPUT FILENAME
    # ======================================================

    print(
        "\nSTEP 7: OUTPUT FASTA FILENAME"
    )


    while True:

        output_name = input(
            "\nOutput FASTA filename "
            "(example: metisa_arthropod_candidates.fasta): "
        ).strip().strip('"')


        if not output_name:

            print(
                "\nERROR: Filename cannot be empty."
            )

            continue


        if not output_name.lower().endswith(
            (
                ".fasta",
                ".faa",
                ".fa"
            )
        ):

            output_name += ".fasta"


        break


    output_path = os.path.join(
        os.path.dirname(
            screening_path
        ),
        output_name
    )


    # ======================================================
    # PREVENT ACCIDENTAL OVERWRITE
    # ======================================================

    if os.path.isfile(
        output_path
    ):

        print(
            "\nWARNING:"
        )

        print(
            "Output file already exists:"
        )

        print(
            output_path
        )


        answer = input(
            "\nOverwrite? (y/n): "
        ).strip().lower()


        if answer != "y":

            print(
                "\nExtraction cancelled."
            )

            return


    # ======================================================
    # STEP 9 — WRITE FASTA
    # ======================================================

    print(
        "\nSTEP 8: WRITING FASTA"
    )


    write_fasta(
        found,
        sequences,
        output_path
    )


    # ======================================================
    # STEP 10 — VERIFY
    # ======================================================

    print(
        "\nSTEP 9: VERIFYING OUTPUT"
    )


    verified = verify_output(
        output_path,
        len(found)
    )


    # ======================================================
    # FINAL SUMMARY
    # ======================================================

    print(
        "\n"
        + "=" * 70
    )

    print(
        "EXTRACTION COMPLETE"
    )

    print(
        "=" * 70
    )


    print(
        f"\nSelected Action:"
    )

    print(
        f"  {selected_action}"
    )


    print(
        f"\nSelected Gene IDs:"
        f" {len(selected_genes):,}"
    )


    print(
        f"Recovered from original FASTA:"
        f" {len(found):,}"
    )


    print(
        f"Missing:"
        f" {len(missing):,}"
    )


    print(
        "\nOutput FASTA:"
    )

    print(
        f"  {output_path}"
    )


    if verified:

        print(
            "\n✓ FASTA verification passed."
        )


    if not missing:

        print(
            "✓ Every selected Gene_ID was "
            "found in the original FASTA."
        )


    print(
        "\nThis FASTA is now ready to be "
        "used as input for InterProScan."
    )


    print(
        "\nDone."
    )


# ==========================================================
# RUN
# ==========================================================

if __name__ == "__main__":

    main()

