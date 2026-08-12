
import os
import statistics
import pandas as pd
import matplotlib.pyplot as plt


# ==========================================================
# METISA PLANA — CODE 1
# 01_input_qc.py
#
# Purpose:
#   Quality control for either:
#       1. DNA / nucleotide FASTA
#       2. Protein FASTA
#
# The script automatically detects the sequence type and
# applies the appropriate QC.
#
# Output:
#   - Excel QC report
#   - TSV QC summary
#   - FASTA length distribution PNG
# ==========================================================


# ==========================================================
# CONSTANTS
# ==========================================================

DNA_ALPHABET = set("ATGCN")

# Standard amino-acid letters + ambiguous X
PROTEIN_ALPHABET = set(
    "ACDEFGHIKLMNPQRSTVWYX"
)


# ==========================================================
# BASIC FUNCTIONS
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


def section(title):

    print(
        "\n"
        + "=" * 70
    )

    print(title)

    print(
        "=" * 70
    )


# ==========================================================
# READ FASTA
# ==========================================================

def read_fasta(input_file):

    sequences = []

    header = None
    sequence = ""

    try:

        with open(
            input_file,
            "r",
            encoding="utf-8",
            errors="replace"
        ) as file:

            for line in file:

                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):

                    # Save previous record
                    if header is not None:

                        sequences.append(
                            (
                                header,
                                sequence
                            )
                        )

                    header = line[1:].strip()
                    sequence = ""

                else:

                    if header is not None:

                        sequence += line.upper()

        # Save final record
        if header is not None:

            sequences.append(
                (
                    header,
                    sequence
                )
            )

    except Exception as error:

        print("\nERROR reading FASTA:")
        print(error)

        return None

    return sequences


# ==========================================================
# DETECT SEQUENCE TYPE
# ==========================================================

def detect_sequence_type(sequences):

    """
    Detect whether the FASTA is DNA/nucleotide or protein.

    Important:
    A sequence containing only A/T/G/C/N is treated as DNA.

    If characters outside the DNA alphabet are present,
    the script checks whether those characters are valid
    amino-acid characters.

    Terminal '*' is ignored during type detection because
    protein FASTA files may contain stop symbols.
    """

    all_letters = set()

    for header, sequence in sequences:

        cleaned = sequence.replace("*", "")

        all_letters.update(
            cleaned
        )

    # Completely empty sequence alphabet
    if not all_letters:

        return "UNKNOWN"

    # Pure DNA alphabet
    if all_letters.issubset(
        DNA_ALPHABET
    ):

        return "DNA"

    # Protein alphabet
    if all_letters.issubset(
        PROTEIN_ALPHABET
    ):

        return "PROTEIN"

    # If neither is cleanly valid
    return "UNKNOWN"


# ==========================================================
# DUPLICATE ID CHECK
# ==========================================================

def calculate_id_statistics(sequences):

    headers = [
        header
        for header, sequence in sequences
    ]

    sequence_ids = [
        header.split()[0]
        for header in headers
    ]

    counts = {}

    for sequence_id in sequence_ids:

        counts[sequence_id] = (
            counts.get(sequence_id, 0) + 1
        )

    duplicated_groups = {

        sequence_id: count

        for sequence_id, count
        in counts.items()

        if count > 1
    }

    duplicate_id_groups = len(
        duplicated_groups
    )

    duplicate_extra_records = sum(
        count - 1
        for count
        in duplicated_groups.values()
    )

    return (
        sequence_ids,
        duplicate_id_groups,
        duplicate_extra_records
    )


# ==========================================================
# COMMON FASTA QC
# ==========================================================

def calculate_common_qc(
    sequences
):

    lengths = [
        len(sequence)
        for header, sequence
        in sequences
    ]

    sequence_ids, duplicate_id_groups, duplicate_extra_records = (
        calculate_id_statistics(
            sequences
        )
    )

    empty_sequences = sum(

        len(sequence) == 0

        for header, sequence
        in sequences
    )

    results = {

        "Number of sequences":
            len(sequences),

        "Number of unique IDs":
            len(set(sequence_ids)),

        "Duplicate ID groups":
            duplicate_id_groups,

        "Extra duplicate ID records":
            duplicate_extra_records,

        "Empty sequences":
            empty_sequences,

    }

    return results


# ==========================================================
# DNA N50 / L50
# ==========================================================

def calculate_n50_l50(lengths):

    if not lengths:

        return 0, 0

    sorted_lengths = sorted(
        lengths,
        reverse=True
    )

    total_length = sum(
        sorted_lengths
    )

    half_total = (
        total_length / 2
    )

    cumulative = 0

    for index, length in enumerate(
        sorted_lengths,
        start=1
    ):

        cumulative += length

        if cumulative >= half_total:

            return (
                length,
                index
            )

    return (
        sorted_lengths[-1],
        len(sorted_lengths)
    )


def calculate_nxx(
    lengths,
    percentage
):

    if not lengths:

        return 0, 0

    sorted_lengths = sorted(
        lengths,
        reverse=True
    )

    total_length = sum(
        sorted_lengths
    )

    target = (
        total_length
        * percentage
        / 100
    )

    cumulative = 0

    for index, length in enumerate(
        sorted_lengths,
        start=1
    ):

        cumulative += length

        if cumulative >= target:

            return (
                length,
                index
            )

    return (
        sorted_lengths[-1],
        len(sorted_lengths)
    )


# ==========================================================
# DNA QC
# ==========================================================

def calculate_dna_qc(
    sequences
):

    section(
        "DNA / NUCLEOTIDE FASTA QC"
    )

    lengths = [
        len(sequence)
        for header, sequence
        in sequences
    ]

    common = calculate_common_qc(
        sequences
    )

    # ------------------------------------------------------
    # Invalid characters
    # ------------------------------------------------------

    invalid_counts = {}

    sequences_with_invalid = 0

    total_invalid = 0

    total_n = 0

    sequences_with_n = 0

    for header, sequence in sequences:

        sequence_invalid = False

        for character in sequence:

            if character not in DNA_ALPHABET:

                invalid_counts[
                    character
                ] = (
                    invalid_counts.get(
                        character,
                        0
                    ) + 1
                )

                total_invalid += 1
                sequence_invalid = True

        n_count = sequence.count(
            "N"
        )

        total_n += n_count

        if n_count > 0:

            sequences_with_n += 1

        if sequence_invalid:

            sequences_with_invalid += 1

    # ------------------------------------------------------
    # Length statistics
    # ------------------------------------------------------

    min_length = min(lengths)
    max_length = max(lengths)

    mean_length = statistics.mean(
        lengths
    )

    median_length = statistics.median(
        lengths
    )

    n50, l50 = calculate_n50_l50(
        lengths
    )

    n90, l90 = calculate_nxx(
        lengths,
        90
    )

    total_length = sum(
        lengths
    )

    n_content = (
        total_n
        / total_length
        * 100
        if total_length > 0
        else 0
    )

    results = common.copy()

    results.update({

        "Sequence type":
            "DNA / nucleotide",

        "Total sequence length (bp)":
            total_length,

        "Minimum sequence length (bp)":
            min_length,

        "Maximum sequence length (bp)":
            max_length,

        "Mean sequence length (bp)":
            round(
                mean_length,
                2
            ),

        "Median sequence length (bp)":
            median_length,

        "N50 (bp)":
            n50,

        "L50":
            l50,

        "N90 (bp)":
            n90,

        "L90":
            l90,

        "Total N bases":
            total_n,

        "N content (%)":
            round(
                n_content,
                4
            ),

        "Sequences containing N":
            sequences_with_n,

        "Total invalid characters":
            total_invalid,

        "Sequences containing invalid characters":
            sequences_with_invalid,

        "Invalid character breakdown":
            "; ".join(

                f"{character}: {count}"

                for character, count

                in sorted(
                    invalid_counts.items()
                )

            )

    })

    return results


# ==========================================================
# PROTEIN QC
# ==========================================================

def calculate_protein_qc(
    sequences
):

    section(
        "PROTEIN FASTA QC"
    )

    cleaned_records = []

    terminal_stop_count = 0

    internal_stop_count = 0

    sequences_with_x = 0

    total_x = 0

    invalid_counts = {}

    sequences_with_invalid = 0

    empty_after_cleanup = 0

    # ------------------------------------------------------
    # Process every protein
    # ------------------------------------------------------

    for header, raw_sequence in sequences:

        sequence = raw_sequence.upper()

        # --------------------------------------------------
        # Remove terminal stop symbol(s)
        #
        # Example:
        #
        # MAAKLL*
        #
        # becomes:
        #
        # MAAKLL
        #
        # --------------------------------------------------

        terminal_stops = 0

        while sequence.endswith("*"):

            sequence = sequence[:-1]

            terminal_stops += 1

        if terminal_stops > 0:

            terminal_stop_count += 1

        # --------------------------------------------------
        # Internal stop symbols
        # --------------------------------------------------

        if "*" in sequence:

            internal_stop_count += 1

        # --------------------------------------------------
        # X residues
        # --------------------------------------------------

        x_count = sequence.count(
            "X"
        )

        total_x += x_count

        if x_count > 0:

            sequences_with_x += 1

        # --------------------------------------------------
        # Invalid amino-acid characters
        # --------------------------------------------------

        sequence_invalid = False

        for character in sequence:

            if character not in PROTEIN_ALPHABET:

                invalid_counts[
                    character
                ] = (
                    invalid_counts.get(
                        character,
                        0
                    ) + 1
                )

                sequence_invalid = True

        if sequence_invalid:

            sequences_with_invalid += 1

        # --------------------------------------------------
        # Empty after terminal-stop removal
        # --------------------------------------------------

        if len(sequence) == 0:

            empty_after_cleanup += 1

        cleaned_records.append(
            (
                header,
                sequence
            )
        )

    # ------------------------------------------------------
    # Lengths AFTER terminal stop removal
    # ------------------------------------------------------

    lengths = [

        len(sequence)

        for header, sequence

        in cleaned_records
    ]

    non_empty_lengths = [

        length

        for length in lengths

        if length > 0

    ]

    if not non_empty_lengths:

        print(
            "\nERROR: No non-empty protein sequences remain."
        )

        return None

    # ------------------------------------------------------
    # Statistics
    # ------------------------------------------------------

    min_length = min(
        non_empty_lengths
    )

    max_length = max(
        non_empty_lengths
    )

    mean_length = statistics.mean(
        non_empty_lengths
    )

    median_length = statistics.median(
        non_empty_lengths
    )

    total_length = sum(
        non_empty_lengths
    )

    # ------------------------------------------------------
    # Short proteins
    # ------------------------------------------------------

    below_50 = sum(

        length < 50

        for length
        in non_empty_lengths
    )

    below_100 = sum(

        length < 100

        for length
        in non_empty_lengths
    )

    # ------------------------------------------------------
    # Invalid characters
    # ------------------------------------------------------

    total_invalid = sum(
        invalid_counts.values()
    )

    results = calculate_common_qc(
        sequences
    )

    results.update({

        "Sequence type":
            "PROTEIN",

        "Protein sequences after terminal-stop cleanup":
            len(non_empty_lengths),

        "Terminal stop-containing sequences":
            terminal_stop_count,

        "Internal stop-containing sequences":
            internal_stop_count,

        "Empty sequences after cleanup":
            empty_after_cleanup,

        "Total protein length (aa)":
            total_length,

        "Minimum protein length (aa)":
            min_length,

        "Maximum protein length (aa)":
            max_length,

        "Mean protein length (aa)":
            round(
                mean_length,
                2
            ),

        "Median protein length (aa)":
            median_length,

        "Proteins <50 aa":
            below_50,

        "Proteins <100 aa":
            below_100,

        "Proteins containing X":
            sequences_with_x,

        "Total X residues":
            total_x,

        "Total invalid amino-acid characters":
            total_invalid,

        "Sequences containing invalid amino-acid characters":
            sequences_with_invalid,

        "Invalid character breakdown":
            "; ".join(

                f"{character}: {count}"

                for character, count

                in sorted(
                    invalid_counts.items()
                )

            )

    })

    return (
        results,
        cleaned_records
    )


# ==========================================================
# OUTPUT FOLDER
# ==========================================================

def ask_output_folder():

    while True:

        folder = input(
            "\nEnter output folder path: "
        ).strip()

        folder = expand_path(
            folder
        )

        if not folder:

            print(
                "\nERROR: Output folder cannot be empty."
            )

            continue

        try:

            os.makedirs(
                folder,
                exist_ok=True
            )

        except Exception as error:

            print(
                "\nERROR creating output folder:"
            )

            print(error)

            continue

        return folder


# ==========================================================
# OUTPUT FILE NAME
# ==========================================================

def ask_output_filename(
    output_folder
):

    while True:

        filename = input(
            "\nEnter output filename "
            "(without extension): "
        ).strip()

        filename = filename.strip('"')

        if not filename:

            print(
                "\nERROR: Filename cannot be empty."
            )

            continue

        filename = os.path.splitext(
            filename
        )[0]

        return os.path.join(
            output_folder,
            filename
        )


# ==========================================================
# SAVE EXCEL
# ==========================================================

def save_excel_report(
    qc_results,
    output_base
):

    output_file = (
        output_base
        + ".xlsx"
    )

    qc_table = pd.DataFrame(

        list(
            qc_results.items()
        ),

        columns=[
            "QC Measurement",
            "Value"
        ]

    )

    try:

        with pd.ExcelWriter(
            output_file,
            engine="openpyxl"
        ) as writer:

            qc_table.to_excel(
                writer,
                sheet_name="QC Summary",
                index=False
            )

        print(
            "\nExcel report saved:"
        )

        print(
            output_file
        )

        return output_file

    except Exception as error:

        print(
            "\nERROR saving Excel:"
        )

        print(error)

        return None


# ==========================================================
# SAVE TSV
# ==========================================================

def save_tsv_report(
    qc_results,
    output_base
):

    output_file = (
        output_base
        + ".tsv"
    )

    qc_table = pd.DataFrame(

        list(
            qc_results.items()
        ),

        columns=[
            "QC_Measurement",
            "Value"
        ]

    )

    try:

        qc_table.to_csv(
            output_file,
            sep="\t",
            index=False
        )

        print(
            "\nTSV report saved:"
        )

        print(
            output_file
        )

        return output_file

    except Exception as error:

        print(
            "\nERROR saving TSV:"
        )

        print(error)

        return None


# ==========================================================
# SAVE LENGTH TABLE
# ==========================================================

def save_sequence_statistics(
    sequences,
    sequence_type,
    output_base
):

    output_file = (
        output_base
        + "_sequence_statistics.tsv"
    )

    rows = []

    for header, sequence in sequences:

        sequence_id = (
            header.split()[0]
        )

        rows.append({

            "Sequence_ID":
                sequence_id,

            "Header":
                header,

            "Length":
                len(sequence),

            "Sequence_Type":
                sequence_type

        })

    df = pd.DataFrame(
        rows
    )

    try:

        df.to_csv(
            output_file,
            sep="\t",
            index=False
        )

        print(
            "\nSequence statistics saved:"
        )

        print(
            output_file
        )

        return output_file

    except Exception as error:

        print(
            "\nERROR saving sequence statistics:"
        )

        print(error)

        return None


# ==========================================================
# SAVE LENGTH DISTRIBUTION
# ==========================================================

def save_length_plot(
    sequences,
    sequence_type,
    output_base
):

    output_file = (
        output_base
        + "_length_distribution.png"
    )

    lengths = [

        len(sequence)

        for header, sequence

        in sequences

    ]

    if not lengths:

        return None

    try:

        plt.figure(
            figsize=(10, 6)
        )

        plt.hist(
            lengths,
            bins=50
        )

        if sequence_type == "PROTEIN":

            plt.xlabel(
                "Protein length (aa)"
            )

            title = (
                "Protein Length Distribution"
            )

        else:

            plt.xlabel(
                "Sequence length (bp)"
            )

            title = (
                "DNA Sequence Length Distribution"
            )

        plt.ylabel(
            "Number of sequences"
        )

        plt.title(
            title
        )

        plt.tight_layout()

        plt.savefig(
            output_file,
            dpi=300
        )

        plt.close()

        print(
            "\nLength distribution plot saved:"
        )

        print(
            output_file
        )

        return output_file

    except Exception as error:

        print(
            "\nERROR creating length plot:"
        )

        print(error)

        return None


# ==========================================================
# MAIN
# ==========================================================

def main():

    section(
        "METISA PLANA — INPUT FASTA QUALITY CONTROL"
    )

    print(
        "\nThis script accepts either:"
    )

    print(
        "  1. DNA / nucleotide FASTA"
    )

    print(
        "  2. Predicted protein FASTA"
    )

    print(
        "\nThe sequence type will be detected automatically."
    )

    # ======================================================
    # STEP 1 — INPUT
    # ======================================================

    section(
        "STEP 1: INPUT FASTA"
    )

    input_file = expand_path(
        input(
            "\nEnter path to FASTA file: "
        )
    )

    if not check_file_exists(
        input_file
    ):

        return

    # ======================================================
    # STEP 2 — READ FASTA
    # ======================================================

    section(
        "STEP 2: READING FASTA"
    )

    print(
        "\nReading:"
    )

    print(
        input_file
    )

    sequences = read_fasta(
        input_file
    )

    if sequences is None:

        return

    print(
        "\nSequences loaded:",
        f"{len(sequences):,}"
    )

    if not sequences:

        print(
            "\nERROR: FASTA contains no sequences."
        )

        return

    # ======================================================
    # STEP 3 — DETECT TYPE
    # ======================================================

    section(
        "STEP 3: SEQUENCE TYPE DETECTION"
    )

    sequence_type = detect_sequence_type(
        sequences
    )

    print(
        "\nDetected sequence type:",
        sequence_type
    )

    if sequence_type == "UNKNOWN":

        print(
            "\nERROR:"
        )

        print(
            "The FASTA could not be confidently classified "
            "as DNA or protein."
        )

        print(
            "\nPlease inspect the FASTA manually."
        )

        return

    # ======================================================
    # STEP 4 — TYPE-SPECIFIC QC
    # ======================================================

    section(
        "STEP 4: TYPE-SPECIFIC QUALITY CONTROL"
    )

    cleaned_sequences = sequences

    if sequence_type == "DNA":

        qc_results = calculate_dna_qc(
            sequences
        )

    else:

        protein_result = calculate_protein_qc(
            sequences
        )

        if protein_result is None:

            return

        qc_results, cleaned_sequences = (
            protein_result
        )

    # ======================================================
    # STEP 5 — DISPLAY RESULTS
    # ======================================================

    section(
        "STEP 5: QC RESULTS"
    )

    for measurement, value in qc_results.items():

        print(
            f"{measurement}: {value}"
        )

    # ======================================================
    # STEP 6 — OUTPUT LOCATION
    # ======================================================

    section(
        "STEP 6: OUTPUT"
    )

    output_folder = ask_output_folder()

    output_base = ask_output_filename(
        output_folder
    )

    print(
        "\nOutput base:"
    )

    print(
        output_base
    )

    # ======================================================
    # STEP 7 — SAVE REPORTS
    # ======================================================

    section(
        "STEP 7: SAVING RESULTS"
    )

    excel_file = save_excel_report(
        qc_results,
        output_base
    )

    tsv_file = save_tsv_report(
        qc_results,
        output_base
    )

    # ======================================================
    # STEP 8 — SAVE PER-SEQUENCE STATISTICS
    # ======================================================

    sequence_stats_file = (
        save_sequence_statistics(
            cleaned_sequences,
            sequence_type,
            output_base
        )
    )

    # ======================================================
    # STEP 9 — SAVE PLOT
    # ======================================================

    plot_file = save_length_plot(
        cleaned_sequences,
        sequence_type,
        output_base
    )

    # ======================================================
    # FINAL SUMMARY
    # ======================================================

    section(
        "INPUT QC COMPLETE"
    )

    print(
        "\nSequence type:",
        sequence_type
    )

    print(
        "Sequences:",
        f"{len(sequences):,}"
    )

    print(
        "\nGenerated files:"
    )

    if excel_file:
        print(
            "  Excel:",
            excel_file
        )

    if tsv_file:
        print(
            "  TSV:",
            tsv_file
        )

    if sequence_stats_file:
        print(
            "  Sequence statistics:",
            sequence_stats_file
        )

    if plot_file:
        print(
            "  Length plot:",
            plot_file
        )

    print(
        "\nQC completed successfully."
    )


# ==========================================================
# RUN
# ==========================================================

if __name__ == "__main__":

    main()

