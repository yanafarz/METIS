import os
import subprocess


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


def section(title):

    print("\n" + "=" * 65)
    print(title)
    print("=" * 65)


def banner():

    print("""
===============================================================
                INTERPROSCAN ANNOTATION TOOL
          Protein Domain / Family Annotation
===============================================================
""")


# ==========================================================
# COUNT FASTA SEQUENCES
# ==========================================================

def count_fasta_sequences(input_file):

    count = 0

    try:

        with open(
            input_file,
            "r",
            encoding="utf-8"
        ) as file:

            for line in file:

                if line.strip().startswith(">"):

                    count += 1

        return count

    except Exception as error:

        print("\nERROR reading FASTA:")
        print(error)

        return 0


# ==========================================================
# CONVERT WINDOWS PATH TO WSL PATH
# ==========================================================

def windows_to_wsl_path(windows_path):

    windows_path = os.path.abspath(
        windows_path
    )

    drive = windows_path[0].lower()

    remaining = windows_path[2:]

    remaining = remaining.replace(
        "\\",
        "/"
    )

    return "/mnt/" + drive + remaining


# ==========================================================
# CHECK INTERPROSCAN
# ==========================================================

def check_interproscan(interproscan_dir):

    interproscan_path = (
        interproscan_dir.rstrip("/") +
        "/interproscan.sh"
    )

    result = subprocess.run(
        [
            "wsl",
            "test",
            "-f",
            interproscan_path
        ],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:

        print(
            "\nERROR: InterProScan not found."
        )

        print(
            "Checked:"
        )

        print(
            interproscan_path
        )

        return False

    print(
        "\n✔ InterProScan found."
    )

    return True


# ==========================================================
# OUTPUT
# ==========================================================

def ask_output(input_file):

    folder = os.path.dirname(
        input_file
    )

    filename = input(
        "\nOutput filename (without extension): "
    ).strip()

    if filename == "":

        filename = "interproscan_result"

    output_file = os.path.join(
        folder,
        filename + ".tsv"
    )

    print(
        "\nOutput location:"
    )

    print(
        output_file
    )

    return output_file


# ==========================================================
# RUN INTERPROSCAN
# ==========================================================

def run_interproscan(
    input_file,
    output_file,
    interproscan_dir,
    sequence_count
):

    # ------------------------------------------------------
    # CONVERT PATHS
    # ------------------------------------------------------

    input_wsl = windows_to_wsl_path(
        input_file
    )

    output_wsl = windows_to_wsl_path(
        output_file
    )

    interproscan_path = (
        interproscan_dir.rstrip("/") +
        "/interproscan.sh"
    )


    # ------------------------------------------------------
    # DISPLAY INFORMATION
    # ------------------------------------------------------

    print(
        "\nInput FASTA:"
    )

    print(
        input_file
    )

    print(
        "\nNumber of sequences:"
    )

    print(
        sequence_count
    )

    print(
        "\nWSL input path:"
    )

    print(
        input_wsl
    )

    print(
        "\nWSL output path:"
    )

    print(
        output_wsl
    )

    print(
        "\nInterProScan:"
    )

    print(
        interproscan_path
    )


    # ------------------------------------------------------
    # BUILD COMMAND
    # ------------------------------------------------------

    command = [

        "wsl",

        "bash",

        interproscan_path,

        "-i",
        input_wsl,

        "-f",
        "tsv",

        "-o",
        output_wsl,

        "-goterms",

        "-iprlookup",

        "-pa"

    ]


    print(
        "\nInterProScan running..."
    )

    print(
        f"Processing {sequence_count:,} protein sequences."
    )

    print(
        "This may take some time."
    )


    print(
        "\nCommand:"
    )

    print(
        " ".join(command)
    )


    print(
        "\nStarting InterProScan...\n"
    )


    # ------------------------------------------------------
    # FORCE JAVA 11
    # ------------------------------------------------------

    env = os.environ.copy()

    env["JAVA_HOME"] = (
        "/usr/lib/jvm/java-11-openjdk-amd64"
    )

    env["PATH"] = (
        env["JAVA_HOME"] +
        "/bin:" +
        env["PATH"]
    )


    # ------------------------------------------------------
    # RUN INTERPROSCAN
    # ------------------------------------------------------

    try:

        result = subprocess.run(
            command,
            text=True,
            env=env
        )


        if result.returncode == 0:

            print(
                "\n✔ InterProScan completed successfully."
            )

            return True


        else:

            print(
                "\nERROR: InterProScan failed."
            )

            print(
                "Return code:",
                result.returncode
            )

            return False


    except Exception as error:

        print(
            "\nERROR:"
        )

        print(
            error
        )

        return False


# ==========================================================
# ADD INTERPROSCAN COLUMN HEADER
# ==========================================================

def add_interproscan_header(output_path):

    header = (
        "Protein Accession\t"
        "Sequence MD5 digest\t"
        "Sequence length\t"
        "Analysis\t"
        "Signature accession\t"
        "Signature description\t"
        "Start location\t"
        "Stop location\t"
        "Score\t"
        "Status\t"
        "Date\t"
        "InterPro accession\t"
        "InterPro description\t"
        "GO terms\t"
        "Pathways\n"
    )


    try:

        with open(
            output_path,
            "r",
            encoding="utf-8"
        ) as file:

            content = file.read()


        with open(
            output_path,
            "w",
            encoding="utf-8"
        ) as file:

            file.write(
                header
            )

            file.write(
                content
            )


        print(
            "\n✔ InterProScan column header added."
        )

        return True


    except Exception as error:

        print(
            "\nERROR adding column header:"
        )

        print(
            error
        )

        return False


# ==========================================================
# MAIN
# ==========================================================

def main():

    banner()


    # ------------------------------------------------------
    # STEP 1: INPUT
    # ------------------------------------------------------

    section(
        "STEP 1: INPUT PROTEIN FASTA"
    )


    input_file = expand_path(
        input(
            "Protein FASTA path: "
        )
    )


    if not check_file_exists(
        input_file
    ):

        return


    # ------------------------------------------------------
    # COUNT SEQUENCES
    # ------------------------------------------------------

    sequence_count = count_fasta_sequences(
        input_file
    )


    if sequence_count == 0:

        print(
            "\nERROR: No FASTA sequences found."
        )

        return


    print(
        "\n✔ FASTA detected."
    )

    print(
        f"✔ Sequences found: {sequence_count:,}"
    )


    # ------------------------------------------------------
    # STEP 2: INTERPROSCAN
    # ------------------------------------------------------

    section(
        "STEP 2: INTERPROSCAN"
    )


    interproscan_dir = input(
        "InterProScan directory path: "
    ).strip()


    if not check_interproscan(
        interproscan_dir
    ):

        return


    # ------------------------------------------------------
    # STEP 3: OUTPUT
    # ------------------------------------------------------

    section(
        "STEP 3: OUTPUT"
    )


    output_file = ask_output(
        input_file
    )


    # ------------------------------------------------------
    # STEP 4: RUN
    # ------------------------------------------------------

    section(
        "STEP 4: RUN INTERPROSCAN"
    )


    success = run_interproscan(

        input_file,

        output_file,

        interproscan_dir,

        sequence_count

    )


    # ------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------

    section(
        "INTERPROSCAN SUMMARY"
    )


    print(
        "INPUT :",
        input_file
    )


    print(
        "SEQUENCES:",
        f"{sequence_count:,}"
    )


    print(
        "INTERPROSCAN:",
        interproscan_dir
    )


    print(
        "OUTPUT:",
        output_file
    )


    if success:

        add_interproscan_header(
            output_file
        )

        print(
            "STATUS: SUCCESS"
        )


    else:

        print(
            "STATUS: FAILED"
        )


# ==========================================================
# START
# ==========================================================

if __name__ == "__main__":

    main()
