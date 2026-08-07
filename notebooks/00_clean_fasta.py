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

        print(
            "\nERROR: File not found:"
        )

        print(path)

        return False

    return True



def banner():

    print("""
===============================================================
              FASTA CLEANING TOOL
     Protein / Gene FASTA Quality Control
===============================================================
""")


# ==========================================================
# DETECT SEQUENCE TYPE
# ==========================================================


def detect_type(sequence):


    letters = set(sequence)


    dna_count = len(
        letters.intersection(DNA_ALPHABET)
    )


    protein_count = len(
        letters.intersection(PROTEIN_ALPHABET)
    )



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
    seq = ""


    try:


        with open(
            input_file,
            "r",
            encoding="utf-8"
        ) as file:


            for line in file:


                line=line.strip()



                if line == "":

                    continue



                if line.startswith(">"):


                    if header is not None:


                        sequences.append(
                            (
                                header,
                                seq
                            )
                        )


                    header=line

                    seq=""



                else:

                    seq += line.upper()



            if header is not None:

                sequences.append(
                    (
                        header,
                        seq
                    )
                )




        print(
            "Sequences found:",
            len(sequences)
        )



        cleaned=[]

        total_removed=0



        for header,seq in sequences:



            seq = seq.replace(
                "*",
                ""
            )



            seq = seq.replace(
                " ",
                ""
            )



            seq_type = detect_type(seq)



            if seq_type=="DNA":


                valid = DNA_ALPHABET



            else:


                valid = PROTEIN_ALPHABET




            new_seq=""


            for letter in seq:


                if letter in valid:

                    new_seq += letter


                else:

                    total_removed += 1




            cleaned.append(
                (
                    header,
                    new_seq
                )
            )





        with open(
            output_file,
            "w",
            encoding="utf-8"
        ) as out:



            for header,seq in cleaned:


                out.write(
                    header+"\n"
                )


                # FASTA wrap 60 characters

                for i in range(
                    0,
                    len(seq),
                    60
                ):

                    out.write(
                        seq[i:i+60]+"\n"
                    )




        print("\nCleaning completed.")

        print(
            "Removed invalid characters:",
            total_removed
        )


        return True



    except Exception as error:


        print(
            "\nERROR:",
            error
        )


        return False



# ==========================================================
# MAIN
# ==========================================================


def main():


    banner()


    input_file = expand_path(

        input(
            "Input FASTA file (.fasta/.faa/.fna): "
        )

    )



    if not check_file_exists(input_file):

        return




    folder = os.path.dirname(
        input_file
    )



    filename=input(
        "\nOutput filename (without extension): "
    ).strip()



    if filename=="":

        filename="cleaned_fasta"



    output_file=os.path.join(

        folder,

        filename + ".fasta"

    )



    success = clean_fasta(

        input_file,

        output_file

    )



    if success:


        print(
            "\nOutput saved:"
        )

        print(
            output_file
        )


    else:


        print(
            "\nFAILED"
        )




if __name__=="__main__":

    main()