#import + configuration
import os
import shutil
import subprocess


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


ANNOTATION_OUTFMT = (
    "6 qseqid sseqid pident length "
    "qlen slen qcovs evalue bitscore stitle"
)

#basic func
def expand_path(path):

    return os.path.abspath(
        os.path.expanduser(
            path.strip()
        )
    )



def check_file_exists(path):

    if not os.path.isfile(path):

        print(f"\nERROR: File not found:\n{path}")

        return False

    return True



def section(title):

    print("\n" + "=" * 65)
    print(title)
    print("=" * 65)



def banner():

    print("""
===============================================================
              NCBI BLAST+ INTERACTIVE RUNNER
          Protein / Gene Homology Annotation Tool
===============================================================
""")

#database check
def find_db_type(db_root):
    # Normal single-volume protein database
    protein_files = [f"{db_root}{ext}" for ext in PROTEIN_DB_EXTS]
    if all(os.path.isfile(path) for path in protein_files):
        return "protein"

    # Normal single-volume nucleotide database
    nucleotide_files = [f"{db_root}{ext}" for ext in NUCLEOTIDE_DB_EXTS]
    if all(os.path.isfile(path) for path in nucleotide_files):
        return "nucleotide"

    # Multi-volume protein database
    protein_volume_files = [
        f"{db_root}.00{ext}" for ext in PROTEIN_DB_EXTS
    ]
    if all(os.path.isfile(path) for path in protein_volume_files):
        return "protein"

    # Multi-volume nucleotide database
    nucleotide_volume_files = [
        f"{db_root}.00{ext}" for ext in NUCLEOTIDE_DB_EXTS
    ]
    if all(os.path.isfile(path) for path in nucleotide_volume_files):
        return "nucleotide"

    # BLAST alias database (.pal = protein, .nal = nucleotide)
    if os.path.isfile(f"{db_root}.pal"):
        return "protein"

    if os.path.isfile(f"{db_root}.nal"):
        return "nucleotide"

    return None



def validate_db_root(db_root):

    if find_db_type(db_root) is None:

        print("\nERROR: BLAST database files not found.")

        print("Expected:")
        print(db_root + ".pin/.psq/.phr")
        print(db_root + ".nin/.nsq/.nhr")

        return False


    return True

#blast + fasta validation
def validate_blast_input(blast_type, query_path, db_root):


    db_type = find_db_type(db_root)


    # DATABASE CHECK

    if blast_type == "blastp":

        if db_type != "protein":

            print("""
ERROR:
blastp requires protein BLAST database.

Need:
.pin
.psq
.phr
""")

            return False



    elif blast_type == "blastn":

        if db_type != "nucleotide":

            print("""
ERROR:
blastn requires nucleotide BLAST database.

Need:
.nin
.nsq
.nhr
""")

            return False



    print("\nChecking FASTA file...")


    try:

        with open(query_path,"r",encoding="utf-8") as file:


            first_line = file.readline()


            if not first_line.startswith(">"):

                print("ERROR: FASTA header missing")

                return False



            sequence = ""


            for line in file:

                if not line.startswith(">"):

                    sequence += line.strip().upper()



        if len(sequence)==0:

            print("ERROR: Empty FASTA")

            return False



        dna_letters = set("ATGCN")

        sequence_letters = set(sequence)



        if blast_type=="blastp":

            if sequence_letters.issubset(dna_letters):

                print("ERROR: Protein BLAST received DNA")

                return False



        if blast_type=="blastn":

            if not sequence_letters.issubset(dna_letters):

                print("ERROR: Nucleotide BLAST received protein")

                return False



    except Exception as error:

        print("ERROR:",error)

        return False



    print("✔ FASTA detected")

    return True

#find blast executable
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



    print("\nERROR: BLAST executable not found.")

    return None

#user input
def ask_blast_type():

    print("""
Select BLAST method:

[1] blastp
    Protein → Protein

[2] blastn
    DNA → DNA
""")


    while True:

        choice = input("Select: ").strip()


        if choice == "":
            return "blastp"


        if choice == "1":

            return "blastp"


        if choice == "2":

            return "blastn"


        print("Invalid choice")




def ask_evalue():


    print("""
E-value:

[1] 1e-3  Exploratory
[2] 1e-5  Standard annotation
[3] 1e-10 High confidence
[4] Custom
""")


    while True:

        choice=input("Select: ").strip()


        if choice=="":
            return "1e-5"


        if choice=="1":
            return "1e-3"


        if choice=="2":
            return "1e-5"


        if choice=="3":
            return "1e-10"


        if choice=="4":

            value=input("Enter E-value: ")

            try:

                float(value)

                return value

            except:

                print("Invalid value")




def ask_output():

    filename=input(
        "\nOutput filename: "
    ).strip()


    if filename=="":
        filename="blast_annotation_result"



    os.makedirs(
        "output",
        exist_ok=True
    )


    output_path=os.path.join(
        "output",
        filename+".tsv"
    )


    print("Output:",os.path.abspath(output_path))


    return output_path

#buid blast command
def build_blast_command(
        blast_program,
        query_path,
        db_root,
        output_path,
        evalue,
        max_target_seqs,
        threads):


    command=[

        blast_program,

        "-query",
        query_path,

        "-db",
        db_root,

        "-out",
        output_path,

        "-outfmt",
        ANNOTATION_OUTFMT,

        "-evalue",
        evalue,

        "-max_target_seqs",
        max_target_seqs,

        "-num_threads",
        threads

    ]


    return command

#header + run blast
def add_output_header(output_path):


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
        "stitle\n"
    )


    with open(output_path,"r") as file:

        content=file.read()



    with open(output_path,"w") as file:

        file.write(header)

        file.write(content)




def run_blast(command):


    print("\nBLAST running...")
    
    print(
        subprocess.list2cmdline(command)
    )


    try:


        result=subprocess.run(
            command,
            capture_output=True,
            text=True
        )



        if result.returncode==0:

            print("\n✔ BLAST completed")

            return True



        else:

            print("\nBLAST failed")

            print(result.stderr)

            return False



    except Exception as error:

        print(error)

        return False

#summary
def print_summary(
        blast_type,
        query,
        database,
        output,
        status):


    print("""
================================================
BLAST SUMMARY
================================================
""")


    print("BLAST :",blast_type)

    print("QUERY :",query)

    print("DATABASE :",database)

    print("OUTPUT :",output)

    print("STATUS :",status)

#main program
def main():


    banner()


    section("STEP 1: BLAST TYPE")


    blast_type=ask_blast_type()



    section("STEP 2: INPUT")


    query_path=expand_path(
        input("Query FASTA path: ")
    )


    if not check_file_exists(query_path):

        return



    db_root=expand_path(
        input("BLAST database path: ")
    )



    if not validate_db_root(db_root):

        return



    if not validate_blast_input(
        blast_type,
        query_path,
        db_root
    ):

        return




    section("STEP 3: PARAMETERS")


    evalue=ask_evalue()


    max_target_seqs=input(
        "Maximum target sequences: "
    ).strip()


    if max_target_seqs=="":
        max_target_seqs="5"



    threads=input(
        "CPU threads: "
    ).strip()


    if threads=="":
        threads="4"




    section("STEP 4: OUTPUT")


    output_path=ask_output()



    blast_program=find_blast_executable(
        blast_type
    )



    if blast_program is None:

        return



    command=build_blast_command(

        blast_program,

        query_path,

        db_root,

        output_path,

        evalue,

        max_target_seqs,

        threads

    )



    success=run_blast(command)



    if success:

        add_output_header(output_path)

        status="SUCCESS"


    else:

        status="FAILED"



    print_summary(

        blast_type,

        query_path,

        db_root,

        output_path,

        status

    )

#run
main()