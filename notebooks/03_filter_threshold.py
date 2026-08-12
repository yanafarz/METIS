import os
import csv
import re


# ==========================================================
# CONFIGURATION
# ==========================================================

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
# ASK OUTPUT NAME
# ==========================================================

def ask_output_name(output_dir):

    print("\nOUTPUT FILE NAME")

    print(
        "\nEnter the output name WITHOUT an extension."
    )

    print(
        "Example:"
    )

    print(
        "metisa_blast_30id_70cov"
    )

    while True:

        name = input(
            "\nOutput name: "
        ).strip()

        if not name:

            print(
                "\nERROR: Output name cannot be empty."
            )

            continue

        # --------------------------------------------------
        # Remove unsafe Windows filename characters
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
        # Output files
        # --------------------------------------------------

        filtered_blast_file = os.path.join(
            output_dir,
            f"{name}_filtered_blast.tsv"
        )

        filtered_fasta_file = os.path.join(
            output_dir,
            f"{name}_filtered.fasta"
        )

        # --------------------------------------------------
        # Prevent accidental overwrite
        # --------------------------------------------------

        existing_files = []

        for path in [
            filtered_blast_file,
            filtered_fasta_file
        ]:

            if os.path.exists(path):

                existing_files.append(path)

        if existing_files:

            print(
                "\nERROR: These output files already exist:"
            )

            for path in existing_files:

                print(
                    f"  {os.path.basename(path)}"
                )

            print(
                "\nPlease choose a different output name."
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
            # Skip header
            # --------------------------------------------------

            if row[0].strip().lower() == "qseqid":

                continue

            # --------------------------------------------------
            # Check column count
            # --------------------------------------------------

            if len(row) < 10:

                print(
                    f"WARNING: Skipping line "
                    f"{line_number}: "
                    f"only {len(row)} columns found."
                )

                continue

            # --------------------------------------------------
            # Parse values
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

    print("\nReading original protein FASTA...")

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

                # --------------------------------------------------
                # FASTA ID = first whitespace-separated word
                # --------------------------------------------------

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
# FILTER BLAST HITS
# ==========================================================

def filter_blast_hits(rows):

    print(
        "\nApplying BLAST quality thresholds..."
    )

    filtered_rows = []

    for hit in rows:

        if passes_quality_filter(hit):

            filtered_rows.append(hit)

    removed = (
        len(rows)
        - len(filtered_rows)
    )

    print(
        f"\nOriginal BLAST hits: "
        f"{len(rows):,}"
    )

    print(
        f"Hits passing threshold: "
        f"{len(filtered_rows):,}"
    )

    print(
        f"Hits removed: "
        f"{removed:,}"
    )

    return filtered_rows


# ==========================================================
# GET UNIQUE PASSING GENE IDS
# ==========================================================

def get_passing_gene_ids(filtered_rows):

    passing_ids = []

    seen = set()

    for hit in filtered_rows:

        gene_id = hit["qseqid"]

        if gene_id not in seen:

            seen.add(gene_id)

            passing_ids.append(
                gene_id
            )

    return passing_ids


# ==========================================================
# SAVE FILTERED BLAST TSV
# ==========================================================

def save_filtered_blast(
    filtered_rows,
    output_path
):

    print(
        "\nSaving filtered BLAST TSV..."
    )

    with open(
        output_path,
        "w",
        encoding="utf-8",
        newline=""
    ) as file:

        writer = csv.writer(
            file,
            delimiter="\t",
            lineterminator="\n"
        )

        # --------------------------------------------------
        # Header
        # --------------------------------------------------

        writer.writerow(
            EXPECTED_COLUMNS
        )

        # --------------------------------------------------
        # Passing hits
        # --------------------------------------------------

        for hit in filtered_rows:

            writer.writerow([

                hit["qseqid"],

                hit["sseqid"],

                hit["pident"],

                hit["length"],

                hit["qlen"],

                hit["slen"],

                hit["qcovs"],

                hit["evalue"],

                hit["bitscore"],

                hit["stitle"],

            ])

    print(
        "\nFiltered BLAST TSV saved:"
    )

    print(
        output_path
    )


# ==========================================================
# SAVE FILTERED FASTA
# ==========================================================

def save_filtered_fasta(
    passing_ids,
    sequences,
    output_path
):

    print(
        "\nCreating FASTA containing "
        "only BLAST-passing proteins..."
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

            # --------------------------------------------------
            # FASTA header
            # --------------------------------------------------

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
        f"Filtered FASTA sequences written: "
        f"{written:,}"
    )

    # ------------------------------------------------------
    # Missing sequences
    # ------------------------------------------------------

    if missing:

        print(
            f"\nWARNING: {len(missing):,} "
            "passing BLAST IDs were not found "
            "in the original FASTA."
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
            "\nMissing IDs saved:"
        )

        print(
            missing_file
        )

    print(
        "\nFiltered FASTA saved:"
    )

    print(
        output_path
    )


# ==========================================================
# PRINT SUMMARY
# ==========================================================

def print_summary(
    original_rows,
    filtered_rows,
    passing_ids,
    sequences,
    output_dir,
    filtered_blast_file,
    filtered_fasta_file
):

    missing_count = 0

    for gene_id in passing_ids:

        if gene_id not in sequences:

            missing_count += 1

    print(
        "\n"
        + "=" * 65
    )

    print(
        "BLAST FILTER SUMMARY"
    )

    print(
        "=" * 65
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
        "\nResults:"
    )

    print(
        f"  Original BLAST hits:     "
        f"{len(original_rows):,}"
    )

    print(
        f"  Filtered BLAST hits:     "
        f"{len(filtered_rows):,}"
    )

    print(
        f"  Unique passing proteins: "
        f"{len(passing_ids):,}"
    )

    print(
        f"  FASTA sequences written:  "
        f"{len(passing_ids) - missing_count:,}"
    )

    if missing_count:

        print(
            f"  Missing FASTA sequences:  "
            f"{missing_count:,}"
        )

    print(
        "\nOutput files:"
    )

    print(
        f"  Filtered BLAST:"
    )

    print(
        f"    {filtered_blast_file}"
    )

    print(
        f"\n  Filtered FASTA:"
    )

    print(
        f"    {filtered_fasta_file}"
    )

    print(
        "\nOutput folder:"
    )

    print(
        f"  {output_dir}"
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
        "       METISA PLANA BLAST FILTER"
    )

    print(
        "=" * 65
    )

    print(
        "\nThis program performs ONLY BLAST threshold filtering."
    )

    print(
        "Taxonomy screening is NOT performed."
    )

    print(
        "InterPro annotation is NOT performed."
    )

    # ======================================================
    # STEP 1: BLAST RESULT
    # ======================================================

    print(
        "\nSTEP 1: BLAST RESULT"
    )

    blast_file = expand_path(
        input(
            "\nPath to BLAST TSV result: "
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
            "\nPath to original protein FASTA: "
        )
    )

    if not check_file_exists(
        fasta_file
    ):

        return

    # ======================================================
    # STEP 3: OUTPUT DIRECTORY
    # ======================================================
    #
    # SAME DIRECTORY AS BLAST INPUT
    # ======================================================

    output_dir = os.path.dirname(
        blast_file
    )

    # ======================================================
    # STEP 4: OUTPUT NAME
    # ======================================================

    output_name = ask_output_name(
        output_dir
    )

    # ======================================================
    # OUTPUT PATHS
    # ======================================================

    filtered_blast_file = os.path.join(
        output_dir,
        f"{output_name}_filtered_blast.tsv"
    )

    filtered_fasta_file = os.path.join(
        output_dir,
        f"{output_name}_filtered.fasta"
    )

    # ======================================================
    # STEP 5: READ BLAST
    # ======================================================

    print(
        "\nSTEP 3: READING BLAST"
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
    # STEP 6: READ ORIGINAL FASTA
    # ======================================================

    print(
        "\nSTEP 4: READING ORIGINAL FASTA"
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
    # STEP 7: FILTER BLAST
    # ======================================================

    print(
        "\nSTEP 5: BLAST QUALITY FILTER"
    )

    print(
        "\nThresholds:"
    )

    print(
        f"  Identity >= {MIN_IDENTITY}%"
    )

    print(
        f"  Coverage >= {MIN_COVERAGE}%"
    )

    print(
        f"  E-value  <= {MAX_EVALUE}"
    )

    filtered_rows = filter_blast_hits(
        rows
    )

    # ======================================================
    # SAVE FILTERED BLAST
    # ======================================================

    save_filtered_blast(
        filtered_rows,
        filtered_blast_file
    )

    # ======================================================
    # STOP IF NOTHING PASSED
    # ======================================================

    if not filtered_rows:

        print(
            "\nNO BLAST HITS PASSED THE THRESHOLD."
        )

        print(
            "\nNo filtered FASTA can be created."
        )

        return

    # ======================================================
    # STEP 6: GET UNIQUE PASSING PROTEINS
    # ======================================================

    print(
        "\nSTEP 6: SELECTING PASSING PROTEINS"
    )

    passing_ids = get_passing_gene_ids(
        filtered_rows
    )

    print(
        f"Unique proteins passing BLAST filter: "
        f"{len(passing_ids):,}"
    )

    # ======================================================
    # STEP 7: CREATE FILTERED FASTA
    # ======================================================

    print(
        "\nSTEP 7: CREATING FILTERED FASTA"
    )

    save_filtered_fasta(
        passing_ids,
        sequences,
        filtered_fasta_file
    )

    # ======================================================
    # SUMMARY
    # ======================================================

    print_summary(
        original_rows=rows,
        filtered_rows=filtered_rows,
        passing_ids=passing_ids,
        sequences=sequences,
        output_dir=output_dir,
        filtered_blast_file=filtered_blast_file,
        filtered_fasta_file=filtered_fasta_file
    )

    # ======================================================
    # DONE
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
        "\nThe filtered FASTA is ready for InterProScan."
    )


# ==========================================================
# RUN
# ==========================================================

if __name__ == "__main__":

    main()
