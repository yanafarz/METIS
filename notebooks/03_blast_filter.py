import os
import pandas as pd


# ==========================================================
# CONFIGURATION
# ==========================================================

MIN_IDENTITY = 30
MIN_COVERAGE = 70
MAX_EVALUE = 0.05


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
            f"\nERROR: File not found:\n{path}"
        )

        return False

    return True



def section(title):

    print("\n" + "=" * 65)
    print(title)
    print("=" * 65)



def banner():

    print("""
===============================================================
              BLAST RESULT FILTERING TOOL
        Identity / Coverage / E-value Filtering
===============================================================
""")


# ==========================================================
# FILTER FUNCTION
# ==========================================================


def filter_blast_result(
        input_file,
        output_file):


    print("\nReading BLAST result...")


    try:

        df = pd.read_csv(
            input_file,
            sep="\t"
        )


        print(
            "Original hits:",
            len(df)
        )


        filtered = df[

            (df["pident"] >= MIN_IDENTITY) &

            (df["qcovs"] >= MIN_COVERAGE) &

            (df["evalue"] <= MAX_EVALUE)

        ]


        print(
            "Filtered hits:",
            len(filtered)
        )



        filtered.to_csv(
            output_file,
            sep="\t",
            index=False
        )


        print(
            "\nFiltering completed."
        )

        return True



    except Exception as error:

        print(
            "\nERROR:",
            error
        )

        return False



# ==========================================================
# OUTPUT
# ==========================================================


def ask_output():


    filename = input(
        "\nOutput filename (without extension): "
    ).strip()


    if filename == "":

        filename = "blast_filtered_result"



    os.makedirs(
        "output",
        exist_ok=True
    )


    output_path = os.path.join(
        "output",
        filename + ".tsv"
    )


    print(
        "\nOutput location:",
        os.path.abspath(output_path)
    )


    return output_path



# ==========================================================
# MAIN PROGRAM
# ==========================================================


def main():


    banner()



    section(
        "STEP 1: INPUT BLAST RESULT"
    )


    input_file = expand_path(

        input(
            "BLAST result TSV path: "
        )

    )



    if not check_file_exists(input_file):

        return




    section(
        "STEP 2: FILTER SETTINGS"
    )


    print(
        """
Current filtering criteria:

Percentage identity >= 30%
Query coverage >= 70%
E-value <= 0.05
"""
    )



    section(
        "STEP 3: OUTPUT"
    )


    output_file = ask_output()



    success = filter_blast_result(

        input_file,

        output_file

    )



    if success:

        print(
            "\nSUCCESS"
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