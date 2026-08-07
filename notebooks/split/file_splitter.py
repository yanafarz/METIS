#!/usr/bin/env python3

import os
from pathlib import Path


def read_fasta(file_path):
    """
    Read FASTA file and return list of sequences.
    Each item: (header, sequence)
    """
    sequences = []

    header = None
    seq = []

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header:
                    sequences.append((header, "".join(seq)))

                header = line
                seq = []

            else:
                seq.append(line)

        # save last sequence
        if header:
            sequences.append((header, "".join(seq)))

    return sequences


def write_fasta(output_path, sequences):
    """
    Write sequences into FASTA file
    """

    with open(output_path, "w") as file:
        for header, sequence in sequences:
            file.write(header + "\n")

            # wrap sequence every 60 characters
            for i in range(0, len(sequence), 60):
                file.write(sequence[i:i+60] + "\n")


def split_evenly(sequences, number_of_files):
    """
    Split sequences equally
    """

    total = len(sequences)

    size = total // number_of_files
    remainder = total % number_of_files

    chunks = []

    start = 0

    for i in range(number_of_files):

        extra = 1 if i < remainder else 0

        end = start + size + extra

        chunks.append(sequences[start:end])

        start = end

    return chunks


def split_by_number(sequences, per_file):
    """
    Split by fixed number of sequences per file
    """

    chunks = []

    for i in range(0, len(sequences), per_file):
        chunks.append(sequences[i:i+per_file])

    return chunks


def main():

    print("=== FASTA SPLITTER ===")

    fasta_path = input(
        "Enter FASTA file path: "
    ).strip()

    fasta_path = Path(fasta_path)

    if not fasta_path.exists():
        print("ERROR: FASTA file not found")
        return


    print("\nReading FASTA...")
    sequences = read_fasta(fasta_path)

    total = len(sequences)

    print(f"\nTotal sequences found: {total}")


    print("\nChoose splitting method:")
    print("1. Split evenly by number of files")
    print("2. Split by number of sequences per file")


    choice = input(
        "Choose option (1/2): "
    ).strip()


    if choice == "1":

        number_of_files = int(
            input("How many files to create? ")
        )

        chunks = split_evenly(
            sequences,
            number_of_files
        )


    elif choice == "2":

        per_file = int(
            input(
                "How many sequences per file? "
            )
        )

        chunks = split_by_number(
            sequences,
            per_file
        )

    else:
        print("Invalid choice")
        return



    prefix = input(
        "Output file name prefix (without extension): "
    ).strip()


    if not prefix:
        prefix = fasta_path.stem



    output_folder = fasta_path.parent


    print("\nWriting files...")


    for i, chunk in enumerate(chunks, start=1):

        output_file = output_folder / (
            f"{prefix}_{i}.fasta"
        )

        write_fasta(
            output_file,
            chunk
        )

        print(
            f"{output_file.name}: {len(chunk)} sequences"
        )


    print("\nDONE")
    print(
        f"Created {len(chunks)} FASTA files"
    )


if __name__ == "__main__":
    main()
