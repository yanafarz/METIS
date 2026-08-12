
import os


# ==========================================================
# CODE 2 — FASTA PREPARATION
# ==========================================================

DNA_ALPHABET = set("ATGCN")

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


def banner():

    print(
        "\n"
        + "=" * 70
    )

    print(
        "             FASTA PREPARATION TOOL"
    )

    print(
        "       DNA / PROTEIN → BLAST-READY FASTA"
    )

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

                if header is None:

                    raise ValueError(
                        "Sequence data was found before "
                        "the first FASTA header."
                    )

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
# DETECT DNA / PROTEIN
# ==========================================================

def detect_sequence_type(sequences):

    if not sequences:

        return "UNKNOWN"

    all_letters = set()

    for header, sequence in sequences:

        all_letters.update(
            sequence.upper()
        )

    # Remove formatting characters
    all_letters.discard(" ")
    all_letters.discard("\t")
    all_letters.discard("*")

    if all_letters.issubset(
        DNA_ALPHABET
    ):

        return "DNA"

    if all_letters.issubset(
        PROTEIN_ALPHABET
    ):

        return "PROTEIN"

    return "UNKNOWN"


# ==========================================================
# PREPARE SEQUENCES
# ==========================================================

def prepare_sequences(
    sequences,
    sequence_type
):

    prepared = []

    terminal_stops_removed = 0

    internal_stops = 0

    invalid_characters = {}

    empty_sequences = 0

    duplicate_ids = 0

    sequence_ids = []


    # ------------------------------------------------------
    # Select valid alphabet
    # ------------------------------------------------------

    if sequence_type == "DNA":

        valid_alphabet = DNA_ALPHABET

    elif sequence_type == "PROTEIN":

        valid_alphabet = PROTEIN_ALPHABET

    else:

        raise ValueError(
            "Unable to determine whether the FASTA "
            "contains DNA or protein sequences."
        )


    # ------------------------------------------------------
    # Process every sequence
    # ------------------------------------------------------

    for header, sequence in sequences:

        sequence = sequence.upper()

        sequence_id = header.split()[0]

        sequence_ids.append(
            sequence_id
        )


        # --------------------------------------------------
        # Remove whitespace
        # --------------------------------------------------

        sequence = (
            sequence
            .replace(" ", "")
            .replace("\t", "")
        )


        # --------------------------------------------------
        # Protein terminal stop
        # --------------------------------------------------

        if sequence_type == "PROTEIN":

            if sequence.endswith("*"):

                sequence = sequence[:-1]

                terminal_stops_removed += 1


        # --------------------------------------------------
        # Detect internal stop
        # --------------------------------------------------

        if "*" in sequence:

            internal_stops += 1


        # --------------------------------------------------
        # Check invalid characters
        # --------------------------------------------------

        cleaned_sequence = []

        for character in sequence:

            if character in valid_alphabet:

                cleaned_sequence.append(
                    character
                )

            elif character == "*":

                # Do NOT silently remove internal stops.
                # Keep track of them so the user knows.
                continue

            else:

                invalid_characters[
                    character
                ] = (
                    invalid_characters.get(
                        character,
                        0
                    ) + 1
                )


        sequence = "".join(
            cleaned_sequence
        )


        # --------------------------------------------------
        # Empty sequence
        # --------------------------------------------------

        if not sequence:

            empty_sequences += 1


        prepared.append(
            (
                header,
                sequence
            )
        )


    # ------------------------------------------------------
    # Duplicate IDs
    # ------------------------------------------------------

    duplicate_ids = (
        len(sequence_ids)
        - len(set(sequence_ids))
    )


    return (
        prepared,
        terminal_stops_removed,
        internal_stops,
        invalid_characters,
        empty_sequences,
        duplicate_ids
    )


# ==========================================================
# SAVE FASTA
# ==========================================================

def save_fasta(
    sequences,
    output_file
):

    with open(
        output_file,
        "w",
        encoding="utf-8"
    ) as file:

        for header, sequence in sequences:

            file.write(
                ">" + header + "\n"
            )

            # FASTA line wrapping
            for i in range(
                0,
                len(sequence),
                60
            ):

                file.write(
                    sequence[i:i + 60]
                    + "\n"
                )


# ==========================================================
# SAVE SUMMARY
# ==========================================================

def save_summary(
    input_file,
    output_file,
    sequence_type,
    original_count,
    prepared_count,
    terminal_stops_removed,
    internal_stops,
    invalid_characters,
    empty_sequences,
    duplicate_ids,
    summary_file
):

    with open(
        summary_file,
        "w",
        encoding="utf-8"
    ) as file:

        file.write(
            "FASTA PREPARATION SUMMARY\n"
        )

        file.write(
            "=" * 60 + "\n\n"
        )

        file.write(
            f"Input file: {input_file}\n"
        )

        file.write(
            f"Output file: {output_file}\n"
        )

        file.write(
            f"Detected sequence type: {sequence_type}\n"
        )

        file.write(
            f"Original sequences: {original_count}\n"
        )

        file.write(
            f"Prepared sequences: {prepared_count}\n"
        )

        file.write(
            f"Terminal protein '*' removed: "
            f"{terminal_stops_removed}\n"
        )

        file.write(
            f"Sequences containing internal '*': "
            f"{internal_stops}\n"
        )

        file.write(
            f"Empty sequences after preparation: "
            f"{empty_sequences}\n"
        )

        file.write(
            f"Duplicate sequence IDs: "
            f"{duplicate_ids}\n"
        )

        file.write(
            "\nInvalid characters detected:\n"
        )

        if invalid_characters:

            for character, count in sorted(
                invalid_characters.items()
            ):

                file.write(
                    f"  {character}: {count}\n"
                )

        else:

            file.write(
                "  None\n"
            )

        file.write(
            "\nRecommended BLAST program:\n"
        )

        if sequence_type == "DNA":

            file.write(
                "  BLASTN\n"
            )

        elif sequence_type == "PROTEIN":

            file.write(
                "  BLASTP\n"
            )

        else:

            file.write(
                "  Cannot determine automatically\n"
            )


# ==========================================================
# OUTPUT NAME
# ==========================================================

def ask_output_name(input_file):

    folder = os.path.dirname(
        input_file
    )

    filename = input(
        "\nOutput filename "
        "(without extension): "
    ).strip()

    if not filename:

        filename = "prepared_sequences"

    filename = os.path.splitext(
        filename
    )[0]

    output_file = os.path.join(
        folder,
        filename + ".fasta"
    )

    summary_file = os.path.join(
        folder,
        filename + "_summary.txt"
    )

    return (
        output_file,
        summary_file
    )


# ==========================================================
# MAIN
# ==========================================================

def main():

    banner()


    # ======================================================
    # STEP 1 — INPUT
    # ======================================================

    input_file = expand_path(
        input(
            "\nInput FASTA file "
            "(.fasta/.faa/.fna): "
        )
    )


    if not check_file_exists(
        input_file
    ):

        return


    # ======================================================
    # STEP 2 — READ
    # ======================================================

    print(
        "\nReading FASTA..."
    )

    try:

        sequences = read_fasta(
            input_file
        )

    except Exception as error:

        print(
            "\nERROR reading FASTA:"
        )

        print(error)

        return


    if not sequences:

        print(
            "\nERROR: No FASTA sequences found."
        )

        return


    print(
        "Sequences found:",
        f"{len(sequences):,}"
    )


    # ======================================================
    # STEP 3 — DETECT TYPE
    # ======================================================

    sequence_type = detect_sequence_type(
        sequences
    )


    print(
        "Detected sequence type:",
        sequence_type
    )


    if sequence_type == "UNKNOWN":

        print(
            "\nERROR:"
        )

        print(
            "The sequence alphabet could not be "
            "confidently classified as DNA or protein."
        )

        print(
            "The input has not been modified."
        )

        return


    # ======================================================
    # STEP 4 — OUTPUT
    # ======================================================

    output_file, summary_file = (
        ask_output_name(
            input_file
        )
    )


    print(
        "\nPrepared FASTA:"
    )

    print(
        output_file
    )

    print(
        "\nPreparation summary:"
    )

    print(
        summary_file
    )


    # ======================================================
    # STEP 5 — PREPARE
    # ======================================================

    print(
        "\nPreparing sequences..."
    )


    (
        prepared,
        terminal_stops_removed,
        internal_stops,
        invalid_characters,
        empty_sequences,
        duplicate_ids
    ) = prepare_sequences(
        sequences,
        sequence_type
    )


    # ======================================================
    # STEP 6 — SAVE FASTA
    # ======================================================

    try:

        save_fasta(
            prepared,
            output_file
        )

    except Exception as error:

        print(
            "\nERROR saving FASTA:"
        )

        print(error)

        return


    # ======================================================
    # STEP 7 — SAVE SUMMARY
    # ======================================================

    try:

        save_summary(
            input_file,
            output_file,
            sequence_type,
            len(sequences),
            len(prepared),
            terminal_stops_removed,
            internal_stops,
            invalid_characters,
            empty_sequences,
            duplicate_ids,
            summary_file
        )

    except Exception as error:

        print(
            "\nERROR saving summary:"
        )

        print(error)

        return


    # ======================================================
    # STEP 8 — REPORT
    # ======================================================

    print(
        "\n"
        + "=" * 70
    )

    print(
        "                 PREPARATION COMPLETE"
    )

    print(
        "=" * 70
    )

    print(
        "\nSequence type:",
        sequence_type
    )

    print(
        "Original sequences:",
        f"{len(sequences):,}"
    )

    print(
        "Prepared sequences:",
        f"{len(prepared):,}"
    )

    print(
        "Terminal protein '*' removed:",
        terminal_stops_removed
    )

    print(
        "Sequences containing internal '*':",
        internal_stops
    )

    print(
        "Empty sequences:",
        empty_sequences
    )

    print(
        "Duplicate IDs:",
        duplicate_ids
    )

    print(
        "Invalid characters detected:",
        sum(
            invalid_characters.values()
        )
    )


    print(
        "\nPrepared FASTA:"
    )

    print(
        output_file
    )


    print(
        "\nSummary:"
    )

    print(
        summary_file
    )


    # ======================================================
    # BLAST RECOMMENDATION
    # ======================================================

    print(
        "\nRecommended downstream BLAST:"
    )

    if sequence_type == "PROTEIN":

        print(
            "BLASTP"
        )

    else:

        print(
            "BLASTN"
        )


    print(
        "\nNOTE:"
    )

    print(
        "This script does NOT predict genes,"
    )

    print(
        "does NOT translate DNA into protein,"
    )

    print(
        "and does NOT run BLAST."
    )

    print(
        "It only prepares the supplied sequences "
        "for downstream analysis."
    )


# ==========================================================
# RUN
# ==========================================================

if __name__ == "__main__":

    main()

