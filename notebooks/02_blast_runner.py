import os
import shutil
import subprocess
import re


# ==========================================================
# CONFIGURATION
# ==========================================================

BLAST_TYPES = {
    "blastp": "protein",
    "blastn": "nucleotide",
}

PROTEIN_DB_EXTS = [
    ".pin",
    ".psq",
    ".phr"
]

NUCLEOTIDE_DB_EXTS = [
    ".nin",
    ".nsq",
    ".nhr"
]

# BLAST output WITHOUT source.
# Source will be added by Python afterward.
ANNOTATION_OUTFMT = (
    "6 qseqid sseqid pident length "
    "qlen slen qcovs evalue bitscore stitle"
)

# Default number of hits retained per query
DEFAULT_MAX_TARGET_SEQS = "10"

# Default CPU threads
DEFAULT_THREADS = "4"


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

        print(f"\nERROR: File not found:")
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
                 UNIProt BLAST ANNOTATION
===============================================================
""")


# ==========================================================
# DATABASE CHECK
# ==========================================================

def find_db_type(db_root):

    # ------------------------------------------------------
    # Normal single-volume protein database
    # ------------------------------------------------------

    protein_files = [
        f"{db_root}{ext}"
        for ext in PROTEIN_DB_EXTS
    ]

    if all(os.path.isfile(path) for path in protein_files):

        return "protein"


    # ------------------------------------------------------
    # Normal single-volume nucleotide database
    # ------------------------------------------------------

    nucleotide_files = [
        f"{db_root}{ext}"
        for ext in NUCLEOTIDE_DB_EXTS
    ]

    if all(os.path.isfile(path) for path in nucleotide_files):

        return "nucleotide"


    # ------------------------------------------------------
    # Multi-volume protein database
    # ------------------------------------------------------

    protein_volume_files = [
        f"{db_root}.00{ext}"
        for ext in PROTEIN_DB_EXTS
    ]

    if all(os.path.isfile(path) for path in protein_volume_files):

        return "protein"


    # ------------------------------------------------------
    # Multi-volume nucleotide database
    # ------------------------------------------------------

    nucleotide_volume_files = [
        f"{db_root}.00{ext}"
        for ext in NUCLEOTIDE_DB_EXTS
    ]

    if all(os.path.isfile(path) for path in nucleotide_volume_files):

        return "nucleotide"


    # ------------------------------------------------------
    # BLAST alias database
    # ------------------------------------------------------

    if os.path.isfile(f"{db_root}.pal"):

        return "protein"


    if os.path.isfile(f"{db_root}.nal"):

        return "nucleotide"


    return None


def validate_db_root(db_root):

    if find_db_type(db_root) is None:

        print("\nERROR: BLAST database files not found.")

        print("\nExpected protein database:")
        print(db_root + ".pin/.psq/.phr")

        print("\nExpected nucleotide database:")
        print(db_root + ".nin/.nsq/.nhr")

        return False

    return True


# ==========================================================
# BLAST + FASTA VALIDATION
# ==========================================================

def validate_blast_input(
    blast_type,
    query_path,
    db_root
):

    db_type = find_db_type(db_root)


    # ------------------------------------------------------
    # DATABASE CHECK
    # ------------------------------------------------------

    if blast_type == "blastp":

        if db_type != "protein":

            print("""
ERROR:
blastp requires a protein BLAST database.

Required:
.pin
.psq
.phr
""")

            return False


    elif blast_type == "blastn":

        if db_type != "nucleotide":

            print("""
ERROR:
blastn requires a nucleotide BLAST database.

Required:
.nin
.nsq
.nhr
""")

            return False


    # ------------------------------------------------------
    # FASTA CHECK
    # ------------------------------------------------------

    print("\nChecking FASTA file...")


    try:

        with open(
            query_path,
            "r",
            encoding="utf-8"
        ) as file:

            first_line = file.readline()


            if not first_line.startswith(">"):

                print(
                    "ERROR: FASTA header missing."
                )

                return False


            sequence = ""


            for line in file:

                if not line.startswith(">"):

                    sequence += line.strip().upper()


        if len(sequence) == 0:

            print("ERROR: Empty FASTA.")

            return False


        dna_letters = set("ATGCN")

        sequence_letters = set(sequence)


        # --------------------------------------------------
        # Protein validation
        # --------------------------------------------------

        if blast_type == "blastp":

            if sequence_letters.issubset(dna_letters):

                print(
                    "ERROR: Protein BLAST received DNA."
                )

                return False


        # --------------------------------------------------
        # Nucleotide validation
        # --------------------------------------------------

        if blast_type == "blastn":

            if not sequence_letters.issubset(dna_letters):

                print(
                    "ERROR: Nucleotide BLAST received protein."
                )

                return False


    except Exception as error:

        print("ERROR:", error)

        return False


    print("✔ FASTA detected")

    return True


# ==========================================================
# FIND BLAST EXECUTABLE
# ==========================================================

def find_blast_executable(blast_type):

    program = shutil.which(blast_type)


    if program:

        return program


    possible = [

        r"C:\BLAST\blast-2.17.0+\bin",

        r"C:\BLAST\bin",

        r"C:\NCBI"

    ]


    for folder in possible:

        exe = os.path.join(
            folder,
            blast_type + ".exe"
        )


        if os.path.isfile(exe):

            return exe


    print(
        "\nERROR: BLAST executable not found."
    )

    return None


# ==========================================================
# BLAST TYPE
# ==========================================================

def ask_blast_type():

    print("""
Select BLAST method:

[1] blastp
    Protein → Protein

[2] blastn
    DNA → DNA
""")


    while True:

        choice = input(
            "Select: "
        ).strip()


        if choice == "":
            return "blastp"


        if choice == "1":
            return "blastp"


        if choice == "2":
            return "blastn"


        print("Invalid choice.")


# ==========================================================
# OUTPUT
# ==========================================================

def ask_output():

    filename = input(
        "\nOutput filename: "
    ).strip()


    if filename == "":

        filename = "blast_annotation_result"


    os.makedirs(
        "output",
        exist_ok=True
    )


    output_path = os.path.join(
        "output",
        filename + ".tsv"
    )


    print(
        "Output:",
        os.path.abspath(output_path)
    )


    return output_path


# ==========================================================
# BUILD BLAST COMMAND
# ==========================================================

def build_blast_command(
    blast_program,
    query_path,
    db_root,
    output_path,
    max_target_seqs,
    threads
):

    command = [

        blast_program,

        "-query",
        query_path,

        "-db",
        db_root,

        "-out",
        output_path,

        "-outfmt",
        ANNOTATION_OUTFMT,

        "-max_target_seqs",
        max_target_seqs,

        "-num_threads",
        threads

    ]


    return command


# ==========================================================
# DETERMINE UNIPROT SOURCE
# ==========================================================

def determine_uniprot_source(sseqid):

    """
    Determine whether a UniProt BLAST hit comes from:

        Swiss-Prot
        TrEMBL
        Isoform

    Examples:

        sp|A5D794|GAPD1_BOVIN
        -> Swiss-Prot

        tr|A0A123|PROTEIN_X
        -> TrEMBL

        sp|P12345-2|PROTEIN-2
        -> Isoform

        tr|A0A123-3|PROTEIN-3
        -> Isoform
    """


    if not sseqid:

        return "Unknown"


    # ------------------------------------------------------
    # Remove whitespace if any
    # ------------------------------------------------------

    sseqid = sseqid.strip()


    # ------------------------------------------------------
    # Split UniProt identifier
    #
    # Example:
    # sp|A5D794|GAPD1_BOVIN
    #
    # parts[0] = sp
    # parts[1] = A5D794
    # ------------------------------------------------------

    parts = sseqid.split("|")


    if len(parts) < 2:

        return "Unknown"


    database_type = parts[0].lower()

    accession = parts[1]


    # ------------------------------------------------------
    # Check for UniProt isoform accession
    #
    # Example:
    #
    # P12345-2
    # A0A123-3
    # ------------------------------------------------------

    if re.search(
        r"-\d+$",
        accession
    ):

        return "Isoform"


    # ------------------------------------------------------
    # Swiss-Prot
    # ------------------------------------------------------

    if database_type == "sp":

        return "Swiss-Prot"


    # ------------------------------------------------------
    # TrEMBL
    # ------------------------------------------------------

    if database_type == "tr":

        return "TrEMBL"


    # ------------------------------------------------------
    # Unknown source
    # ------------------------------------------------------

    return "Unknown"


# ==========================================================
# ADD HEADER + SOURCE
# ==========================================================

def add_header_and_source(output_path):

    header = (
        "qseqid\t"
        "sseqid\t"
        "pident\t"
        "length\t"
        "qlen\t"
        "slen\t"
        "qcovs\t"
        "evalue\t"
        "bitscore\t"
        "stitle\t"
        "source\n"
    )


    # ------------------------------------------------------
    # Read BLAST results
    # ------------------------------------------------------

    with open(
        output_path,
        "r",
        encoding="utf-8"
    ) as file:

        lines = file.readlines()


    new_lines = []


    # ------------------------------------------------------
    # Add source to every BLAST hit
    # ------------------------------------------------------

    for line in lines:

        line = line.rstrip("\n")


        if not line.strip():

            continue


        fields = line.split("\t")


        # Expected BLAST columns:
        #
        # 0 qseqid
        # 1 sseqid
        # 2 pident
        # 3 length
        # 4 qlen
        # 5 slen
        # 6 qcovs
        # 7 evalue
        # 8 bitscore
        # 9 stitle

        if len(fields) < 10:

            continue


        sseqid = fields[1]


        source = determine_uniprot_source(
            sseqid
        )


        fields.append(source)


        new_lines.append(
            "\t".join(fields) + "\n"
        )


    # ------------------------------------------------------
    # Rewrite output
    # ------------------------------------------------------

    with open(
        output_path,
        "w",
        encoding="utf-8"
    ) as file:

        file.write(header)

        file.writelines(new_lines)


# ==========================================================
# RUN BLAST
# ==========================================================

def run_blast(command):

    print("\nBLAST running...")

    print(
        subprocess.list2cmdline(command)
    )


    try:

        result = subprocess.run(
            command,
            capture_output=True,
            text=True
        )


        if result.returncode == 0:

            print(
                "\n✔ BLAST completed"
            )

            return True


        else:

            print(
                "\nBLAST failed"
            )

            print(
                result.stderr
            )

            return False


    except Exception as error:

        print(
            "\nERROR:",
            error
        )

        return False


# ==========================================================
# SUMMARY
# ==========================================================

def print_summary(
    blast_type,
    query,
    database,
    output,
    status
):

    print("""
===============================================================
                         SUMMARY
===============================================================
""")


    print(
        "BLAST    :",
        blast_type
    )

    print(
        "QUERY    :",
        query
    )

    print(
        "DATABASE :",
        database
    )

    print(
        "OUTPUT   :",
        output
    )

    print(
        "STATUS   :",
        status
    )


# ==========================================================
# MAIN PROGRAM
# ==========================================================

def main():

    banner()


    # ------------------------------------------------------
    # STEP 1
    # ------------------------------------------------------

    section(
        "STEP 1: BLAST TYPE"
    )


    blast_type = ask_blast_type()


    # ------------------------------------------------------
    # STEP 2
    # ------------------------------------------------------

    section(
        "STEP 2: INPUT"
    )


    query_path = expand_path(
        input(
            "Query FASTA path: "
        )
    )


    if not check_file_exists(
        query_path
    ):

        return


    db_root = expand_path(
        input(
            "BLAST database path: "
        )
    )


    if not validate_db_root(
        db_root
    ):

        return


    if not validate_blast_input(
        blast_type,
        query_path,
        db_root
    ):

        return


    # ------------------------------------------------------
    # STEP 3
    # ------------------------------------------------------

    section(
        "STEP 3: PARAMETERS"
    )


    max_target_seqs = input(
        "Maximum target sequences: "
    ).strip()


    if max_target_seqs == "":

        max_target_seqs = DEFAULT_MAX_TARGET_SEQS


    threads = input(
        "CPU threads: "
    ).strip()


    if threads == "":

        threads = DEFAULT_THREADS


    # ------------------------------------------------------
    # STEP 4
    # ------------------------------------------------------

    section(
        "STEP 4: OUTPUT"
    )


    output_path = ask_output()


    # ------------------------------------------------------
    # FIND BLAST
    # ------------------------------------------------------

    blast_program = find_blast_executable(
        blast_type
    )


    if blast_program is None:

        return


    # ------------------------------------------------------
    # BUILD COMMAND
    # ------------------------------------------------------

    command = build_blast_command(

        blast_program,

        query_path,

        db_root,

        output_path,

        max_target_seqs,

        threads

    )


    # ------------------------------------------------------
    # RUN BLAST
    # ------------------------------------------------------

    success = run_blast(
        command
    )


    # ------------------------------------------------------
    # ADD SOURCE COLUMN
    # ------------------------------------------------------

    if success:

        add_header_and_source(
            output_path
        )

        status = "SUCCESS"


    else:

        status = "FAILED"


    # ------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------

    print_summary(

        blast_type,

        query_path,

        db_root,

        output_path,

        status

    )


# ==========================================================
# RUN
# ==========================================================

if __name__ == "__main__":

    main()
