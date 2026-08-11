
import os
import statistics
import pandas as pd
import matplotlib.pyplot as plt


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
# BANNER
# ==========================================================

def banner():

    print(
        "\n"
        + "=" * 65
    )

    print(
        "              FASTA QUALITY CONTROL"
    )

    print(
        "=" * 65
    )


# ==========================================================
# READ FASTA
# ==========================================================

def read_fasta(input_file):

    sequences = []

    header = None
    sequence = ""

    with open(
        input_file,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as file:

        for line in file:

            line = line.strip()

            if line == "":
                continue

            if line.startswith(">"):

                # Save previous sequence
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

        # Save final sequence
        if header is not None:

            sequences.append(
                (
                    header,
                    sequence
                )
            )

    return sequences


# ==========================================================
# DETECT SEQUENCE TYPE
# ==========================================================

def detect_sequence_type(sequences):

    dna_letters = set("ATGCN")

    all_letters = set()

    for header, sequence in sequences:

        all_letters.update(
            sequence.replace("*", "")
        )

    if all_letters.issubset(dna_letters):

        return "DNA / GENE"

    else:

        return "PROTEIN"


# ==========================================================
# CALCULATE QC
# ==========================================================

def calculate_qc(sequences):

    if len(sequences) == 0:

        print(
            "ERROR: No FASTA sequences found."
        )

        return None

    headers = [

        header

        for header, sequence in sequences

    ]

    sequences_only = [

        sequence

        for header, sequence in sequences

    ]

    # ------------------------------------------------------
    # Empty sequences
    # ------------------------------------------------------

    empty_sequences = sum(

        len(sequence) == 0

        for sequence in sequences_only

    )

    # ------------------------------------------------------
    # Sequence lengths
    # ------------------------------------------------------

    lengths = [

        len(sequence)

        for sequence in sequences_only

    ]

    # ------------------------------------------------------
    # Sequence IDs
    # ------------------------------------------------------

    sequence_ids = [

        header.split()[0]

        for header in headers

    ]

    # ------------------------------------------------------
    # Duplicate IDs
    # ------------------------------------------------------

    unique_ids = set(
        sequence_ids
    )

    duplicate_ids = (

        len(sequence_ids)

        -

        len(unique_ids)

    )

    # ------------------------------------------------------
    # Length statistics
    # ------------------------------------------------------

    min_length = min(
        lengths
    )

    max_length = max(
        lengths
    )

    mean_length = statistics.mean(
        lengths
    )

    median_length = statistics.median(
        lengths
    )

    # ------------------------------------------------------
    # Short sequences
    # ------------------------------------------------------

    below_50 = sum(

        length < 50

        for length in lengths

    )

    # ------------------------------------------------------
    # X residues
    # ------------------------------------------------------

    containing_x = sum(

        "X" in sequence

        for sequence in sequences_only

    )

    # ------------------------------------------------------
    # Terminal stop
    # ------------------------------------------------------

    ending_stop = sum(

        sequence.endswith("*")

        for sequence in sequences_only

    )

    # ------------------------------------------------------
    # Internal stop
    # ------------------------------------------------------

    internal_stop = 0

    for sequence in sequences_only:

        if "*" in sequence[:-1]:

            internal_stop += 1

    # ------------------------------------------------------
    # Results
    # ------------------------------------------------------

    results = {

        "Number of sequences":
            len(sequences),

        "Number of unique IDs":
            len(unique_ids),

        "Duplicate IDs":
            duplicate_ids,

        "Empty sequences":
            empty_sequences,

        "Minimum length":
            min_length,

        "Maximum length":
            max_length,

        "Mean length":
            round(
                mean_length,
                2
            ),

        "Median length":
            median_length,

        "Sequences <50 residues":
            below_50,

        "Sequences containing X":
            containing_x,

        "Sequences ending with *":
            ending_stop,

        "Sequences with internal *":
            internal_stop

    }

    return results


# ==========================================================
# ASK OUTPUT FILE NAME
# ==========================================================

def ask_output_file(input_file):

    folder = os.path.dirname(
        input_file
    )

    print(
        "\nOUTPUT FILE"
    )

    print(
        "The report will be saved in:"
    )

    print(
        folder
    )

    while True:

        filename = input(
            "\nEnter Excel output filename "
            "(example: fasta_qc_report.xlsx): "
        ).strip()

        # --------------------------------------------------
        # Prevent empty filename
        # --------------------------------------------------

        if not filename:

            print(
                "\nERROR: Filename cannot be empty."
            )

            continue

        # --------------------------------------------------
        # Remove accidental quotes
        # --------------------------------------------------

        filename = filename.strip('"')

        # --------------------------------------------------
        # Add .xlsx automatically
        # --------------------------------------------------

        if not filename.lower().endswith(
            ".xlsx"
        ):

            filename += ".xlsx"

        # --------------------------------------------------
        # Prevent path input
        #
        # User should enter filename only.
        # --------------------------------------------------

        if (
            os.path.dirname(filename)
            or
            os.path.isabs(filename)
        ):

            print(
                "\nERROR: Enter a filename only."
            )

            print(
                "Do not enter a folder path."
            )

            continue

        output_file = os.path.join(
            folder,
            filename
        )

        # --------------------------------------------------
        # Existing file
        # --------------------------------------------------

        if os.path.exists(
            output_file
        ):

            print(
                "\nWARNING: File already exists:"
            )

            print(
                output_file
            )

            overwrite = input(
                "\nOverwrite this file? (y/n): "
            ).strip().lower()

            if overwrite != "y":

                print(
                    "\nPlease enter a different filename."
                )

                continue

        return output_file


# ==========================================================
# SAVE EXCEL REPORT
# ==========================================================

def save_excel_report(
    qc_results,
    output_file
):

    qc_table = pd.DataFrame(

        list(
            qc_results.items()
        ),

        columns=[
            "QC Measurement",
            "Value"
        ]

    )

    qc_table.to_excel(
        output_file,
        index=False
    )

    print(
        "\nQC report saved:"
    )

    print(
        output_file
    )


# ==========================================================
# MAIN PROGRAM
# ==========================================================

banner()


# ==========================================================
# INPUT FASTA
# ==========================================================

input_file = expand_path(

    input(
        "\nInput FASTA file "
        "(.fasta/.faa/.fna): "
    )

)


if not check_file_exists(
    input_file
):

    raise SystemExit


# ==========================================================
# READ FASTA
# ==========================================================

print(
    "\nReading FASTA..."
)


sequences = read_fasta(
    input_file
)


print(
    "Sequences loaded:",
    len(sequences)
)


if not sequences:

    print(
        "\nERROR: No FASTA sequences found."
    )

    raise SystemExit


# ==========================================================
# DETECT TYPE
# ==========================================================

sequence_type = detect_sequence_type(
    sequences
)


print(
    "Detected type:",
    sequence_type
)


# ==========================================================
# CALCULATE QC
# ==========================================================

qc_results = calculate_qc(
    sequences
)


if qc_results is None:

    raise SystemExit


# ==========================================================
# PRINT QC RESULTS
# ==========================================================

print(
    "\n"
    + "=" * 65
)

print(
    "                         QC RESULTS"
)

print(
    "=" * 65
)


for measurement, value in qc_results.items():

    print(
        f"{measurement}: {value}"
    )


# ==========================================================
# CREATE QC TABLE
# ==========================================================

qc_table = pd.DataFrame(

    list(
        qc_results.items()
    ),

    columns=[
        "QC Measurement",
        "Value"
    ]

)


# ==========================================================
# LENGTH HISTOGRAM
# ==========================================================

lengths = [

    len(sequence)

    for header, sequence in sequences

]


plt.figure(
    figsize=(10, 6)
)


plt.hist(
    lengths,
    bins=50
)


plt.xlabel(
    "Sequence length (aa / nt)"
)

plt.ylabel(
    "Number of sequences"
)

plt.title(
    f"FASTA Sequence Length Distribution — "
    f"{sequence_type}"
)


plt.tight_layout()

plt.show()


# ==========================================================
# ASK OUTPUT FILE NAME
# ==========================================================

output_file = ask_output_file(
    input_file
)


# ==========================================================
# SAVE EXCEL
# ==========================================================

save_excel_report(
    qc_results,
    output_file
)


# ==========================================================
# FINISH
# ==========================================================

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

