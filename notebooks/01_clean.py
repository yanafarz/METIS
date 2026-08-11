import os


# ==========================================================
# CONFIGURATION
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
                 FASTA CLEANING TOOL
          DNA / Protein FASTA Quality Control
===============================================================
""")


# ==========================================================
# DETECT SEQUENCE TYPE
# ==========================================================

def detect_sequence_type(sequence):

    letters = set(sequence)

    if letters.issubset(DNA_ALPHABET):

        return "DNA"

    else:

        return "PROTEIN"


# ==========================================================
# CLEAN FASTA
# ==========================================================

def clean_fasta(input_file, output_file):

    print("\nReading FASTA...")

    sequences = []

    header = None
    sequence = ""

    try:

        # --------------------------------------------------
        # READ FASTA
        # --------------------------------------------------

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

                    header = line
                    sequence = ""

                else:

                    sequence += line.upper()


            # Save final sequence
            if header is not None:

                sequences.append(
                    (
                        header,
                        sequence
                    )
                )


        print(
            "Sequences found:",
            len(sequences)
        )


        # --------------------------------------------------
        # DETECT FILE TYPE
        # --------------------------------------------------

        all_sequence = ""

        for header, sequence in sequences:

            all_sequence += sequence


        sequence_type = detect_sequence_type(
            all_sequence
        )


        print(
            "Detected type:",
            sequence_type
        )


        # --------------------------------------------------
        # CLEAN
        # --------------------------------------------------

        cleaned = []

        terminal_stops_removed = 0
        invalid_characters_removed = 0


        for header, sequence in sequences:


            # Remove spaces
            sequence = sequence.replace(
                " ",
                ""
            )


            # --------------------------------------------------
            # PROTEIN
            # --------------------------------------------------

            if sequence_type == "PROTEIN":

                # Remove ONLY terminal stop
                if sequence.endswith("*"):

                    sequence = sequence[:-1]

                    terminal_stops_removed += 1


                valid_alphabet = PROTEIN_ALPHABET


            # --------------------------------------------------
            # DNA
            # --------------------------------------------------

            else:

                valid_alphabet = DNA_ALPHABET


            # --------------------------------------------------
            # REMOVE INVALID CHARACTERS
            # --------------------------------------------------

            new_sequence = ""

            for letter in sequence:

                if letter in valid_alphabet:

                    new_sequence += letter

                else:

                    invalid_characters_removed += 1


            cleaned.append(
                (
                    header,
                    new_sequence
                )
            )

        # --------------------------------------------------
        # CHECK EMPTY SEQUENCES
        # --------------------------------------------------

        empty_after_cleaning = sum(
            len(sequence) == 0
            for header, sequence in cleaned
        )

        print(
            "Empty sequences after cleaning:",
            empty_after_cleaning
        )  


        # --------------------------------------------------
        # WRITE OUTPUT
        # --------------------------------------------------

        with open(
            output_file,
            "w",
            encoding="utf-8"
        ) as out:

            for header, sequence in cleaned:

                out.write(
                    header + "\n"
                )

                # FASTA wrap: 60 characters
                for i in range(
                    0,
                    len(sequence),
                    60
                ):

                    out.write(
                        sequence[i:i+60] + "\n"
                    )


        # --------------------------------------------------
        # SUMMARY
        # --------------------------------------------------

        print("\nCleaning completed.")

        print(
            "Sequences kept:",
            len(cleaned)
        )

        print(
            "Detected type:",
            sequence_type
        )

        print(
            "Terminal '*' removed:",
            terminal_stops_removed
        )

        print(
            "Other invalid characters removed:",
            invalid_characters_removed
        )


        if sequence_type == "PROTEIN":

            print(
                "X residues kept: YES"
            )


        print(
            "Short sequences removed: NO"
        )

        print(
            "Duplicate IDs removed: NO"
        )


        return True


    except Exception as error:

        print(
            "\nERROR:",
            error
        )

        return False


# ==========================================================
# MAIN PROGRAM
# ==========================================================

def main():

    banner()


    # ------------------------------------------------------
    # INPUT
    # ------------------------------------------------------

    input_file = expand_path(

        input(
            "Input FASTA file (.fasta/.faa/.fna): "
        )

    )


    if not check_file_exists(input_file):

        return


    # ------------------------------------------------------
    # OUTPUT
    # ------------------------------------------------------

    folder = os.path.dirname(
        input_file
    )


    filename = input(
        "\nOutput filename (without extension): "
    ).strip()


    if filename == "":

        filename = "cleaned_fasta"


    output_file = os.path.join(
        folder,
        filename + ".fasta"
    )


    print(
        "\nOutput location:"
    )

    print(
        output_file
    )


    # ------------------------------------------------------
    # CLEAN
    # ------------------------------------------------------

    success = clean_fasta(

        input_file,

        output_file

    )


    if success:

        print(
            "\nSUCCESS"
        )

        print(
            "Clean FASTA saved:"
        )

        print(
            output_file
        )

    else:

        print(
            "\nFAILED"
        )


# ==========================================================
# START
# ==========================================================

if __name__ == "__main__":

    main()
    
