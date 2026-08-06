import os
import shutil
import subprocess

BLAST_TYPES = {
    "blastn": "nucleotide",
    "blastp": "protein",
    "blastx": "protein",
    "tblastn": "nucleotide",
    "tblastx": "nucleotide",
}

PROTEIN_DB_EXTS = [".pin", ".psq", ".phr"]
NUCLEOTIDE_DB_EXTS = [".nin", ".nsq", ".nhr"]


def expand_path(value):
    return os.path.abspath(os.path.expanduser(value.strip()))


def prompt_text(prompt, default=None, allow_empty=False):
    if default:
        prompt = f"{prompt} [{default}]: "
    else:
        prompt = f"{prompt}: "

    while True:
        value = input(prompt).strip()
        if value == "" and default is not None:
            return default
        if value == "" and allow_empty:
            return ""
        if value:
            return value
        print("Please enter a value or press Enter to accept the default.")


def prompt_choice(prompt, choices, default=None):
    choices_text = ", ".join(choices)
    while True:
        value = prompt_text(f"{prompt} ({choices_text})", default=default)
        if value in choices:
            return value
        print(f"Invalid choice. Please select one of: {choices_text}")


def prompt_number(prompt, default, cast, positive=True):
    while True:
        value = prompt_text(prompt, default=default)
        try:
            number = cast(value)
            if positive and number <= 0:
                raise ValueError
            return str(number)
        except ValueError:
            print(f"Please enter a valid {'positive ' if positive else ''}{cast.__name__} value.")


def check_file_exists(path):
    if not os.path.isfile(path):
        print(f"ERROR: File not found: {path}")
        return False
    return True


def find_db_type(db_root):
    protein_files = [f"{db_root}{ext}" for ext in PROTEIN_DB_EXTS]
    nucleotide_files = [f"{db_root}{ext}" for ext in NUCLEOTIDE_DB_EXTS]

    if all(os.path.isfile(path) for path in protein_files):
        return "protein"
    if all(os.path.isfile(path) for path in nucleotide_files):
        return "nucleotide"
    return None


def validate_db_root(db_root):
    if not db_root:
        print("ERROR: BLAST database base name cannot be empty.")
        return False
    db_root = expand_path(db_root)
    if find_db_type(db_root) is None:
        print("ERROR: BLAST database files were not found for:")
        print(f"  {db_root}{PROTEIN_DB_EXTS[0]} / {db_root}{PROTEIN_DB_EXTS[1]} / {db_root}{PROTEIN_DB_EXTS[2]}")
        print(f"  or {db_root}{NUCLEOTIDE_DB_EXTS[0]} / {db_root}{NUCLEOTIDE_DB_EXTS[1]} / {db_root}{NUCLEOTIDE_DB_EXTS[2]}")
        return False
    return True


def find_blast_executable(blast_type):
    program = shutil.which(blast_type)
    if program:
        return program
    print(f"ERROR: BLAST executable '{blast_type}' was not found in PATH.")
    print("Please install BLAST+ or make sure the folder containing the BLAST executables is in your PATH.")
    return None


def warn_on_type_mismatch(selected_type, db_type):
    expected = BLAST_TYPES[selected_type]
    if expected != db_type:
        print("WARNING: BLAST type and database type may not match.")
        print(f"  Selected BLAST type: {selected_type} ({expected} database expected)")
        print(f"  Detected database type: {db_type}")
        print("This may still run, but results can be invalid if the database type is wrong.")


def build_blast_command(blast_type, query_path, db_root, output_path, outfmt, evalue, max_target_seqs, num_threads):
    command = [
        blast_type,
        "-query",
        query_path,
        "-db",
        db_root,
        "-out",
        output_path,
        "-outfmt",
        outfmt,
        "-evalue",
        evalue,
        "-max_target_seqs",
        max_target_seqs,
        "-num_threads",
        num_threads,
    ]
    return command


def run_blast(command):
    print("\nExecuting command:")
    print(subprocess.list2cmdline(command))
    try:
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode == 0:
            print("\nBLAST finished successfully.")
            if result.stdout:
                print("--- BLAST output ---")
                print(result.stdout)
            return True, result
        print("\nBLAST failed with return code", result.returncode)
        if result.stderr:
            print("--- BLAST error output ---")
            print(result.stderr.strip())
        return False, result
    except FileNotFoundError:
        print("ERROR: Failed to launch BLAST. The executable was not found.")
        return False, None
    except subprocess.SubprocessError as error:
        print("ERROR: BLAST execution failed:", error)
        return False, None


def ensure_output_directory_exists(output_path):
    directory = os.path.dirname(output_path)
    if not directory:
        return True
    if os.path.isdir(directory):
        return True
    print(f"ERROR: Output directory does not exist: {directory}")
    return False


def main():
    print("=== Interactive NCBI BLAST Runner ===")
    print("This tool builds and runs a BLAST+ command with simple prompts.")
    print("Press Enter to reuse a previous value where shown.")

    previous_inputs = {
        "blast_type": None,
        "query_path": None,
        "db_root": None,
        "output_path": None,
        "evalue": "1e-5",
        "outfmt": "6",
        "max_target_seqs": "5",
        "num_threads": "2",
    }

    while True:
        blast_type = prompt_choice(
            "Choose BLAST type",
            list(BLAST_TYPES.keys()),
            default=previous_inputs["blast_type"] or "blastp",
        )
        previous_inputs["blast_type"] = blast_type

        query_path = prompt_text(
            "Path to query FASTA file",
            default=previous_inputs["query_path"],
        )
        query_path = expand_path(query_path)
        previous_inputs["query_path"] = query_path
        if not check_file_exists(query_path):
            continue

        db_root = prompt_text(
            "Path to BLAST database base name (without extension)",
            default=previous_inputs["db_root"],
        )
        db_root = expand_path(db_root)
        previous_inputs["db_root"] = db_root
        if not validate_db_root(db_root):
            continue

        output_path = prompt_text(
            "Output file path",
            default=previous_inputs["output_path"],
        )
        output_path = expand_path(output_path)
        previous_inputs["output_path"] = output_path
        if not ensure_output_directory_exists(output_path):
            continue

        evalue = prompt_text("E-value cutoff", default=previous_inputs["evalue"])
        try:
            float(evalue)
        except ValueError:
            print("ERROR: Please enter a valid numeric E-value, for example 1e-5.")
            continue
        previous_inputs["evalue"] = evalue

        outfmt = prompt_text("Output format", default=previous_inputs["outfmt"])
        if not outfmt:
            print("ERROR: Output format cannot be empty.")
            continue
        previous_inputs["outfmt"] = outfmt

        max_target_seqs = prompt_number(
            "Max target sequences",
            default=previous_inputs["max_target_seqs"],
            cast=int,
            positive=True,
        )
        previous_inputs["max_target_seqs"] = max_target_seqs

        num_threads = prompt_number(
            "Number of threads",
            default=previous_inputs["num_threads"],
            cast=int,
            positive=True,
        )
        previous_inputs["num_threads"] = num_threads

        blast_executable = find_blast_executable(blast_type)
        if blast_executable is None:
            status = "FAILED"
            print_summary(
                blast_type,
                query_path,
                db_root,
                output_path,
                status,
            )
            return

        db_type = find_db_type(db_root)
        if db_type:
            warn_on_type_mismatch(blast_type, db_type)

        command = build_blast_command(
            blast_executable,
            query_path,
            db_root,
            output_path,
            outfmt,
            evalue,
            max_target_seqs,
            num_threads,
        )

        success, _ = run_blast(command)
        status = "SUCCESS" if success else "FAILED"
        print_summary(
            blast_type,
            query_path,
            db_root,
            output_path,
            status,
        )
        return


def print_summary(blast_type, query_path, db_root, output_path, status):
    print("\n=== BLAST RUN SUMMARY ===")
    print(f"Type: {blast_type}")
    print(f"Query: {query_path}")
    print(f"Database: {db_root}")
    print(f"Output: {output_path}")
    print(f"Status: {status}")


if __name__ == "__main__":
    main()
