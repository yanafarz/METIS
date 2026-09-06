import os
import subprocess
import re

PROTEIN_DB_EXTS = [".pin", ".psq", ".phr"]

ANNOTATION_OUTFMT = (
    "6 qseqid sseqid pident length "
    "qlen slen qcov evalue bitscore stitle"
)

DEFAULT_MAX_TARGET_SEQS = "50"
DEFAULT_THREADS = "6"

def expand_path(path):
    return os.path.abspath(os.path.expanduser(path.strip().strip('"')))

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
                 METISA PLANA HOMOLOGY SEARCH
===============================================================
""")

def find_protein_db(db_root):

    # Standard single-volume protein database
    protein_files = [f"{db_root}{ext}" for ext in PROTEIN_DB_EXTS]

    if all(os.path.isfile(path) for path in protein_files):
        return True

    # Multi-volume protein database
    protein_volume_files = [
        f"{db_root}.00{ext}" for ext in PROTEIN_DB_EXTS
    ]

    if all(os.path.isfile(path) for path in protein_volume_files):
        return True

    # Protein alias database
    if os.path.isfile(f"{db_root}.pal"):
        return True

    return False

def validate_db_root(db_root):
    if not find_protein_db(db_root):
        print("\nERROR: Protein BLAST database files were not found.")
        print("\nExpected database:")
        print(db_root + ".pin/.psq/.phr")
        print("\nOr a BLAST protein alias database:")
        print(db_root + ".pal")
        return False

    print("\n✔ Protein BLAST database detected.")
    return True

def validate_protein_fasta(query_path):

    print("\nChecking protein FASTA...")

    # Standard amino-acid alphabet.
    #
    # X = unknown amino acid
    # * = stop codon
    protein_letters = set("ACDEFGHIKLMNPQRSTVWYX*")

    sequence_count = 0
    empty_sequences = 0
    invalid_sequences = 0
    internal_stop_sequences = 0

    current_sequence = []
    has_header = False

    try:

        with open(query_path, "r", encoding="utf-8", errors="replace") as file:

            for line in file:

                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):

                    # Process previous sequence.
                    if has_header:

                        sequence = "".join(current_sequence).upper()
                        sequence_count += 1

                        if not sequence:
                            empty_sequences += 1

                        else:

                            sequence_letters = set(sequence)

                            if not sequence_letters.issubset(protein_letters):
                                invalid_sequences += 1

                            # Internal stop detection.
                            # Terminal * is allowed because the
                            # previous QC script removes terminal *.
                            if "*" in sequence[:-1]:
                                internal_stop_sequences += 1

                    has_header = True
                    current_sequence = []

                else:

                    if not has_header:

                        print(
                            "\nERROR: Sequence data was found "
                            "before the first FASTA header."
                        )

                        return False

                    current_sequence.append(line.upper())

            if has_header:

                sequence = "".join(current_sequence).upper()
                sequence_count += 1

                if not sequence:
                    empty_sequences += 1

                else:

                    sequence_letters = set(sequence)

                    if not sequence_letters.issubset(protein_letters):
                        invalid_sequences += 1

                    if "*" in sequence[:-1]:
                        internal_stop_sequences += 1

        if sequence_count == 0:

            print("\nERROR: No FASTA sequences found.")
            return False

        if empty_sequences > 0:

            print(
                f"\nWARNING: {empty_sequences:,} "
                "empty sequence(s) detected."
            )

        if internal_stop_sequences > 0:

            print(
                f"\nWARNING: {internal_stop_sequences:,} "
                "protein(s) contain internal '*' stop characters."
            )

            print("Review these proteins before BLASTP.")

        if invalid_sequences > 0:

            print(
                f"\nERROR: {invalid_sequences:,} "
                "sequence(s) contain invalid protein characters."
            )

            print("Expected amino-acid characters:")
            print("ACDEFGHIKLMNPQRSTVWYX")

            return False

    except Exception as error:

        print("\nERROR checking FASTA:")
        print(error)

        return False

    print(
        f"✔ Protein FASTA detected: "
        f"{sequence_count:,} sequences"
    )

    if empty_sequences == 0:
        print("✔ No empty sequences detected.")

    if invalid_sequences == 0:
        print("✔ Protein characters are valid.")

    if internal_stop_sequences == 0:
        print("✔ No internal stop characters detected.")

    return True

def find_blastp_executable():
    print("\nEnter the full path to blastp.exe")
    executable = expand_path(input("Path to blastp.exe: ").strip().strip('"'))

    if not os.path.isfile(executable):
        print("\nERROR: blastp.exe not found:")
        print(executable)
        return None

    if os.path.basename(executable).lower() != "blastp.exe":
        print("\nERROR: Selected file is not blastp.exe.")
        return None

    try:
        result = subprocess.run(
            [executable, "-version"],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            print("\nERROR: blastp.exe was found but could not be executed.")
            print(result.stderr)
            return None

        print("\n✔ BLASTP executable found:")
        print(executable)
        print("\n✔ BLASTP is working.")

        version_output = result.stdout.strip().splitlines()

        if version_output:
            print("Version:", version_output[0])

        return executable

    except Exception as error:
        print("\nERROR testing BLASTP:")
        print(error)
        return None

def ask_output(query_path):
    folder = os.path.dirname(query_path)

    filename = input(
        "\nOutput filename (without extension): "
    ).strip().strip('"')

    if not filename:
        filename = "blast_annotation_result"

    filename = os.path.splitext(filename)[0]
    output_path = os.path.join(folder, filename + ".tsv")

    print("\nOutput:")
    print(output_path)

    return output_path

def build_blast_command(
    blast_program, query_path, db_root,
    output_path, max_target_seqs, threads
):
    return [
        blast_program, "-query", query_path,
        "-db", db_root, "-out", output_path,
        "-outfmt", ANNOTATION_OUTFMT,
        "-max_target_seqs", max_target_seqs,
        "-num_threads", threads
    ]

def determine_uniprot_source(sseqid):
    if not sseqid:
        return "Unknown"

    sseqid = sseqid.strip()
    parts = sseqid.split("|")

    if len(parts) < 2:
        return "Unknown"

    database_type = parts[0].lower()
    accession = parts[1]

    # Isoform
    if re.search(r"-\d+$", accession):
        return "Isoform"

    # Swiss-Prot
    if database_type == "sp":
        return "Swiss-Prot"

    # TrEMBL
    if database_type == "tr":
        return "TrEMBL"

    return "Unknown"

def add_header_and_source(output_path):
    header = (
        "qseqid\tsseqid\tpident\tlength\tqlen\tslen\t"
        "qcov\tevalue\tbitscore\tstitle\tsource\n"
    )

    try:
        with open(output_path, "r", encoding="utf-8") as file:
            lines = file.readlines()

        new_lines = []

        for line in lines:
            line = line.rstrip("\n")

            if not line.strip():
                continue

            fields = line.split("\t")

            if len(fields) < 10:
                continue

            source = determine_uniprot_source(fields[1])
            fields.append(source)
            new_lines.append("\t".join(fields) + "\n")

        with open(output_path, "w", encoding="utf-8") as file:
            file.write(header)
            file.writelines(new_lines)

        return True

    except Exception as error:
        print("\nERROR adding output header:")
        print(error)
        return False

def run_blast(command):
    print("\nBLASTP running...")
    print(subprocess.list2cmdline(command))

    try:
        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode == 0:
            print("\n✔ BLASTP completed.")
            return True

        print("\nBLASTP failed.")

        if result.stderr:
            print(result.stderr)

        return False

    except Exception as error:
        print("\nERROR running BLASTP:")
        print(error)
        return False

def print_summary(query, database, output, status):

    print("""
===============================================================
                         SUMMARY
===============================================================
""")

    print("BLAST    : blastp")
    print("QUERY    :", query)
    print("DATABASE :", database)
    print("OUTPUT   :", output)
    print("STATUS   :", status)

def main():

    banner()

    section("STEP 1: INPUT")

    query_path = expand_path(
        input("Query protein FASTA path: ")
    )

    if not check_file_exists(query_path):
        return

    db_root = expand_path(
        input("Protein BLAST database path: ")
    )

    if not validate_db_root(db_root):
        return

    if not validate_protein_fasta(query_path):
        return

    section("STEP 2: PARAMETERS")

    while True:
        max_target_seqs = input(
            "Maximum target sequences per query [50]: "
        ).strip()

        if not max_target_seqs:
            max_target_seqs = DEFAULT_MAX_TARGET_SEQS
            break

        if max_target_seqs.isdigit() and int(max_target_seqs) > 0:
            break

        print("Please enter a positive whole number.")

    while True:
        threads = input("CPU threads [4]: ").strip()

        if not threads:
            threads = DEFAULT_THREADS
            break

        if threads.isdigit() and int(threads) > 0:
            break

        print("Please enter a positive whole number.")

    section("STEP 3: OUTPUT")
    output_path = ask_output(query_path)

    section("STEP 4: BLASTP")
    blast_program = find_blastp_executable()

    if blast_program is None:
        return

    command = build_blast_command(
        blast_program,
        query_path,
        db_root,
        output_path,
        max_target_seqs,
        threads
    )

    section("STEP 5: RUNNING BLASTP")
    success = run_blast(command)

    if not success:
        print_summary(
            query_path,
            db_root,
            output_path,
            "FAILED"
        )
        return

    section("STEP 6: PROCESSING RESULTS")
    processed = add_header_and_source(output_path)

    if not processed:
        print_summary(
            query_path,
            db_root,
            output_path,
            "FAILED DURING RESULT PROCESSING"
        )
        return

    print_summary(
        query_path,
        db_root,
        output_path,
        "SUCCESS"
    )

    print(
        "\nBLASTP result is ready for downstream "
        "filtering and top-hit ranking."
    )

if __name__ == "__main__":
    main()