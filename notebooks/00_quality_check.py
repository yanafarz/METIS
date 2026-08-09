import os
import statistics
import pandas as pd
import matplotlib.pyplot as plt

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


def banner():

    print("""
===============================================================
                 FASTA INPUT QC TOOL
        Protein / Gene FASTA Quality Control
===============================================================
""")

def read_fasta(input_file):

    sequences = []

    header = None
    sequence = ""

    with open(
        input_file,
        "r",
        encoding="utf-8"
    ) as file:

        for line in file:

            line = line.strip()

            if line == "":
                continue

            if line.startswith(">"):

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

                sequence += line.upper()

        if header is not None:

            sequences.append(
                (
                    header,
                    sequence
                )
            )

    return sequences

def detect_sequence_type(sequences):

    dna_letters = set("ATGCN")

    all_letters = set()

    for header, sequence in sequences:

        all_letters.update(sequence.replace("*", ""))

    if all_letters.issubset(dna_letters):

        return "DNA / GENE"

    else:

        return "PROTEIN"

def calculate_qc(sequences):

    headers = [
        header
        for header, sequence in sequences
    ]

    sequences_only = [
        sequence
        for header, sequence in sequences
    ]

    lengths = [
        len(sequence)
        for sequence in sequences_only
    ]

    sequence_ids = [
        header.split()[0]
        for header in headers
    ]


    # ------------------------------------------------------
    # Duplicate IDs
    # ------------------------------------------------------

    unique_ids = set(sequence_ids)

    duplicate_ids = (
        len(sequence_ids) -
        len(unique_ids)
    )


    # ------------------------------------------------------
    # Length statistics
    # ------------------------------------------------------

    min_length = min(lengths)

    max_length = max(lengths)

    mean_length = statistics.mean(lengths)

    median_length = statistics.median(lengths)


    # ------------------------------------------------------
    # Short proteins
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

        "Minimum length":
            min_length,

        "Maximum length":
            max_length,

        "Mean length":
            round(mean_length, 2),

        "Median length":
            median_length,

        "Sequences <50 aa":
            below_50,

        "Sequences containing X":
            containing_x,

        "Sequences ending with *":
            ending_stop,

        "Sequences with internal *":
            internal_stop
    }


    return results

banner()


input_file = expand_path(
    input(
        "Input FASTA file (.fasta/.faa/.fna): "
    )
)


if not check_file_exists(input_file):

    raise SystemExit


print("\nReading FASTA...")


sequences = read_fasta(
    input_file
)


print(
    "Sequences loaded:",
    len(sequences)
)


sequence_type = detect_sequence_type(
    sequences
)


print(
    "Detected type:",
    sequence_type
)


qc_results = calculate_qc(
    sequences
)


print("\n===============================================================")
print("                         QC RESULTS")
print("===============================================================\n")


for measurement, value in qc_results.items():

    print(
        f"{measurement}: {value}"
    )

qc_table = pd.DataFrame(
    list(qc_results.items()),
    columns=[
        "QC Measurement",
        "Value"
    ]
)


qc_table

lengths = [
    len(sequence)
    for header, sequence in sequences
]


plt.figure(figsize=(10, 6))

plt.hist(
    lengths,
    bins=50
)

plt.xlabel("Sequence length (aa / nt)")
plt.ylabel("Number of sequences")
plt.title(
    f"FASTA Sequence Length Distribution — {sequence_type}"
)

plt.tight_layout()

plt.show()

folder = os.path.dirname(
    input_file
)


output_file = os.path.join(
    folder,
    "fasta_qc_report.xlsx"
)


qc_table.to_excel(
    output_file,
    index=False
)


print("\nQC report saved:")
print(output_file)

