import os
import statistics
import pandas as pd
import matplotlib.pyplot as plt

RAW_PROTEIN_FASTA = CLEANED_FASTA = QC_SUMMARY = None
SEQUENCE_STATISTICS = LENGTH_PLOT = None

sequences = cleaned_sequences = qc_results = None

PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYX")

def expand_path(path):
    return os.path.abspath(os.path.expanduser(path.strip().strip('"')))

def check_file_exists(path):
    if not os.path.isfile(path):
        print("\nERROR: File not found:")
        print(path)
        return False
    return True

def section(title):
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)

def read_fasta(input_file):
    sequences = []
    header = None
    sequence = ""

    try:
        with open(input_file, "r", encoding="utf-8", errors="replace") as file:
            for line in file:
                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):
                    if header is not None:
                        sequences.append((header, sequence))

                    header = line[1:].strip()
                    sequence = ""

                else:
                    if header is None:
                        raise ValueError(
                            "Sequence data was found before "
                            "the first FASTA header."
                        )

                    sequence += line.upper()

        if header is not None:
            sequences.append((header, sequence))

    except Exception as error:
        print("\nERROR reading FASTA:")
        print(error)
        return None

    return sequences

def calculate_id_statistics(sequences):
    sequence_ids = [header.split()[0] for header, sequence in sequences]

    counts = {}
    for sequence_id in sequence_ids:
        counts[sequence_id] = counts.get(sequence_id, 0) + 1

    duplicated_groups = {
        sequence_id: count
        for sequence_id, count in counts.items()
        if count > 1
    }

    duplicate_id_groups = len(duplicated_groups)
    duplicate_extra_records = sum(count - 1 for count in duplicated_groups.values())

    return sequence_ids, duplicate_id_groups, duplicate_extra_records

def calculate_protein_qc_and_prepare(sequences):
    cleaned_records = []
    terminal_stop_count = 0
    internal_stop_count = 0
    sequences_with_x = 0
    total_x = 0
    invalid_counts = {}
    sequences_with_invalid = 0
    empty_after_cleanup = 0

    sequence_ids, duplicate_id_groups, duplicate_extra_records = (
        calculate_id_statistics(sequences)
    )

    # Process every protein
    for header, raw_sequence in sequences:
        sequence = raw_sequence.upper().replace(" ", "").replace("\t", "")

        # Remove terminal stop characters
        terminal_stops = 0
        while sequence.endswith("*"):
            sequence = sequence[:-1]
            terminal_stops += 1

        if terminal_stops > 0:
            terminal_stop_count += 1

        # Report internal stops; do not remove them
        if "*" in sequence:
            internal_stop_count += 1

        # X residues
        x_count = sequence.count("X")
        total_x += x_count

        if x_count > 0:
            sequences_with_x += 1

        # Invalid amino-acid characters
        sequence_invalid = False

        for character in sequence:
            if character == "*":
                continue

            if character not in PROTEIN_ALPHABET:
                invalid_counts[character] = invalid_counts.get(character, 0) + 1
                sequence_invalid = True

        if sequence_invalid:
            sequences_with_invalid += 1

        # Terminal "*" already removed; internal "*" retained
        cleaned_sequence = sequence

        if len(cleaned_sequence) == 0:
            empty_after_cleanup += 1

        cleaned_records.append((header, cleaned_sequence))

    # Length statistics
    lengths = [
        len(sequence)
        for header, sequence in cleaned_records
        if len(sequence) > 0
    ]

    if not lengths:
        print("\nERROR: No non-empty protein sequences remain.")
        return None

    # Basic statistics
    min_length = min(lengths)
    max_length = max(lengths)
    mean_length = statistics.mean(lengths)
    median_length = statistics.median(lengths)
    total_length = sum(lengths)

    # Short proteins
    below_50 = sum(length < 50 for length in lengths)
    below_100 = sum(length < 100 for length in lengths)

    # Invalid characters
    total_invalid = sum(invalid_counts.values())

    # QC results
    results = {
        "Number of sequences": len(sequences),
        "Number of unique IDs": len(set(sequence_ids)),
        "Duplicate ID groups": duplicate_id_groups,
        "Extra duplicate ID records": duplicate_extra_records,
        "Empty input sequences": sum(
            len(sequence) == 0
            for header, sequence in sequences
        ),
        "Terminal stop-containing proteins": terminal_stop_count,
        "Internal stop-containing proteins": internal_stop_count,
        "Empty sequences after cleanup": empty_after_cleanup,
        "Proteins <50 aa": below_50,
        "Proteins <100 aa": below_100,
        "Proteins containing X": sequences_with_x,
        "Total X residues": total_x,
        "Total invalid amino-acid characters": total_invalid,
        "Proteins containing invalid amino-acid characters": sequences_with_invalid,
        "Invalid character breakdown": (
            "; ".join(
                f"{character}: {count}"
                for character, count in sorted(invalid_counts.items())
            )
            if invalid_counts else "None"
        ),
        "Total protein length (aa)": total_length,
        "Minimum protein length (aa)": min_length,
        "Maximum protein length (aa)": max_length,
        "Mean protein length (aa)": round(mean_length, 2),
        "Median protein length (aa)": median_length,
        "Proteins retained for downstream analysis": len(cleaned_records),
    }

    return results, cleaned_records

def save_fasta(sequences, output_file):
    try:
        with open(output_file, "w", encoding="utf-8") as file:
            for header, sequence in sequences:
                file.write(f">{header}\n")

                for i in range(0, len(sequence), 60):
                    file.write(sequence[i:i + 60] + "\n")

        print("\nCleaned FASTA saved:")
        print(output_file)

        return True

    except Exception as error:
        print("\nERROR saving FASTA:")
        print(error)
        return False

def save_qc_summary(qc_results, output_file):
    try:
        qc_table = pd.DataFrame(
            list(qc_results.items()),
            columns=["QC_Measurement", "Value"]
        )

        qc_table.to_csv(output_file, sep="\t", index=False)

        print("\nQC summary saved:")
        print(output_file)

        return True

    except Exception as error:
        print("\nERROR saving QC summary:")
        print(error)

        return False

def save_sequence_statistics(sequences, output_file):
    rows = []

    for header, sequence in sequences:
        sequence_id = header.split()[0]

        rows.append({
            "Protein_ID": sequence_id,
            "Header": header,
            "Length_aa": len(sequence),
            "Contains_X": "Yes" if "X" in sequence else "No",
            "Contains_internal_stop": "Yes" if "*" in sequence else "No"
        })

    df = pd.DataFrame(rows)

    try:
        df.to_csv(output_file, sep="\t", index=False)

        print("\nPer-protein statistics saved:")
        print(output_file)

        return True

    except Exception as error:
        print("\nERROR saving per-protein statistics:")
        print(error)

        return False

def save_length_plot(sequences, output_file):
    lengths = [len(sequence) for header, sequence in sequences]

    if not lengths:
        return False

    try:
        plt.figure(figsize=(10, 6))
        plt.hist(lengths, bins=50)

        plt.xlabel("Protein length (aa)")
        plt.ylabel("Number of proteins")
        plt.title("Metisa plana Protein Length Distribution")

        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()

        print("\nProtein length distribution saved:")
        print(output_file)

        return True

    except Exception as error:
        print("\nERROR creating length plot:")
        print(error)

        return False

def ask_output_name(input_file):
    folder = os.path.dirname(input_file)

    while True:
        filename = input(
            "\nEnter output filename (without extension): "
        ).strip().strip('"')

        if not filename:
            print("\nERROR: Filename cannot be empty.")
            continue

        filename = os.path.splitext(filename)[0]

        return (
            os.path.join(folder, filename + ".fasta"),
            os.path.join(folder, filename + "_QC.tsv"),
            os.path.join(folder, filename + "_protein_statistics.tsv"),
            os.path.join(folder, filename + "_length_distribution.png")
        )

section("METISA PLANA — INPUT QC + FASTA PREPARATION")
print("\nThis script checks and prepares a predicted")
print("protein FASTA for downstream analysis.")

section("STEP 1: INPUT FASTA")

RAW_PROTEIN_FASTA = expand_path(
    input("\nEnter path to raw protein FASTA: ")
)

if not check_file_exists(RAW_PROTEIN_FASTA):
    raise FileNotFoundError(RAW_PROTEIN_FASTA)

print("\nRaw protein FASTA:")
print(RAW_PROTEIN_FASTA)

section("STEP 2: READING FASTA")

sequences = read_fasta(RAW_PROTEIN_FASTA)

if sequences is None or not sequences:
    raise ValueError("FASTA contains no sequences.")

print("\nSequences loaded:", f"{len(sequences):,}")

section("STEP 3: PROTEIN QC + PREPARATION")

result = calculate_protein_qc_and_prepare(sequences)

if result is None:
    raise ValueError("No non-empty protein sequences remain.")

qc_results, cleaned_sequences = result

section("STEP 4: QC RESULTS")

for measurement, value in qc_results.items():
    print(f"{measurement}: {value}")

if qc_results["Internal stop-containing proteins"] > 0:
    print("\nWARNING: Internal stop-containing proteins were detected.")
    print("Review these proteins before downstream BLASTP analysis.")

section("STEP 5: OUTPUT")

CLEANED_FASTA, QC_SUMMARY, SEQUENCE_STATISTICS, LENGTH_PLOT = (
    ask_output_name(RAW_PROTEIN_FASTA)
)

section("STEP 6: SAVING CLEANED FASTA")

if not save_fasta(cleaned_sequences, CLEANED_FASTA):
    raise IOError("Failed to save cleaned FASTA.")

section("STEP 7: SAVING QC SUMMARY")

save_qc_summary(qc_results, QC_SUMMARY)

section("STEP 8: SAVING PROTEIN STATISTICS")

save_sequence_statistics(cleaned_sequences, SEQUENCE_STATISTICS)

section("STEP 9: SAVING LENGTH DISTRIBUTION")

save_length_plot(cleaned_sequences, LENGTH_PLOT)

section("INPUT QC + FASTA PREPARATION COMPLETE")

print("\nOriginal proteins:", f"{len(sequences):,}")
print("Prepared proteins:", f"{len(cleaned_sequences):,}")
print("Terminal '*' removed:",
      qc_results["Terminal stop-containing proteins"])
print("Internal '*' detected:",
      qc_results["Internal stop-containing proteins"])
print("Proteins <50 aa:",
      qc_results["Proteins <50 aa"])
print("Proteins containing X:",
      qc_results["Proteins containing X"])

print("\nPipeline variables saved:")
print("  RAW_PROTEIN_FASTA =", RAW_PROTEIN_FASTA)
print("  CLEANED_FASTA =", CLEANED_FASTA)
print("  QC_SUMMARY =", QC_SUMMARY)
print("  SEQUENCE_STATISTICS =", SEQUENCE_STATISTICS)
print("  LENGTH_PLOT =", LENGTH_PLOT)