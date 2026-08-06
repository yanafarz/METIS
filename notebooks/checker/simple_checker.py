#!/usr/bin/env python3
"""BLAST installation checker and database manager."""

import os
import shutil
import subprocess
import urllib.request
from pathlib import Path

BLAST_COMMANDS = ["blastn", "blastp", "makeblastdb", "blastdbcmd"]
COMMON_PATHS = [
    r"C:\Program Files\NCBI",
    r"C:\Program Files (x86)\NCBI",
    r"C:\BLAST",
]
DB_EXTENSIONS = {
    "protein": {"pin", "psq", "phr"},
    "nucleotide": {"nin", "nsq", "nhr"},
}


def find_executable(cmd):
    path = shutil.which(cmd)
    if path:
        return os.path.abspath(path)
    for base in COMMON_PATHS:
        if os.path.isdir(base):
            for root, _, files in os.walk(base):
                if f"{cmd}.exe" in files:
                    return os.path.abspath(os.path.join(root, f"{cmd}.exe"))
    return None


def get_blast_version(executable_path):
    try:
        result = subprocess.run(
            [executable_path, "-version"],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode == 0:
            first_line = result.stdout.splitlines()[0].strip()
            return first_line
    except (OSError, subprocess.SubprocessError):
        pass
    return "Unknown"


def print_installation_summary():
    print("=== BLAST INSTALLATION CHECK ===")
    blast_paths = {cmd: find_executable(cmd) for cmd in BLAST_COMMANDS}

    installed = any(blast_paths.values())
    print(f"BLAST installed: {'YES' if installed else 'NO'}")

    version = "Not available"
    if blast_paths.get("blastn"):
        version = get_blast_version(blast_paths["blastn"])
    elif blast_paths.get("blastp"):
        version = get_blast_version(blast_paths["blastp"])
    elif blast_paths.get("makeblastdb"):
        version = get_blast_version(blast_paths["makeblastdb"])

    print(f"Version: {version}")
    print("Location:")
    for cmd in BLAST_COMMANDS:
        path = blast_paths.get(cmd)
        if path:
            print(f"  {cmd}: {path}")
        else:
            print(f"  {cmd}: Not found")
    print()


def prompt_yes_no(message):
    while True:
        try:
            answer = input(message).strip().lower()
        except EOFError:
            return False
        if answer in {"y", "yes"}:
            return True
        if answer in {"n", "no"}:
            return False
        print("Please enter Y or N.")


def check_internet_connection():
    try:
        request = urllib.request.Request(
            "https://ftp.ncbi.nlm.nih.gov/",
            headers={"User-Agent": "Mozilla/5.0"},
        )
        with urllib.request.urlopen(request, timeout=10) as response:
            return response.status < 400
    except Exception:
        return False


def ensure_download_prerequisites():
    tools = {cmd: find_executable(cmd) for cmd in ["blastdbcmd", "makeblastdb"]}
    missing = [name for name, path in tools.items() if not path]
    if missing:
        print("Missing required BLAST tools:")
        for name in missing:
            print(f"  - {name}")
        print("Install BLAST+ and try again.")
        return None, False

    if not check_internet_connection():
        print("No internet connection detected.")
        print("Please connect to the internet before downloading databases.")
        return tools, False

    return tools, True


def prompt_folder_path(prompt_text, default_name=None):
    try:
        user_input = input(prompt_text).strip()
    except EOFError:
        return None

    if not user_input:
        if default_name:
            return str(Path(default_name).expanduser())
        return None

    return user_input


def resolve_download_path(base_path, default_file_name):
    path = Path(base_path).expanduser()
    if path.suffix:
        return path
    if path.exists() and path.is_dir():
        return path / default_file_name
    if not path.suffix:
        return path.with_suffix(".fasta")
    return path


def download_with_progress(url, destination_path):
    destination = Path(destination_path).expanduser()
    destination.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading: {url}")
    print(f"Saving to: {destination}")

    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        total_size = response.headers.get("Content-Length")
        total_size = int(total_size) if total_size and total_size.isdigit() else None
        downloaded = 0
        chunk_size = 1024 * 1024
        with open(destination, "wb") as handle:
            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                handle.write(chunk)
                downloaded += len(chunk)
                if total_size:
                    percent = int(downloaded * 100 / total_size)
                    print(f"Progress: {percent}% ({downloaded}/{total_size} bytes)", end="\r")
                else:
                    print(f"Progress: {downloaded} bytes downloaded", end="\r")
        print()
    print("Download complete.")
    return destination


def extract_gzip_file(gzip_path, output_path):
    output = Path(output_path).expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)
    if shutil.which("gzip"):
        try:
            with open(output, "wb") as dst:
                result = subprocess.run(
                    ["gzip", "-dc", str(gzip_path)],
                    stdout=dst,
                    stderr=subprocess.PIPE,
                    check=False,
                    text=False,
                )
            if result.returncode != 0:
                raise RuntimeError(result.stderr.decode("utf-8", errors="ignore"))
            return output
        except OSError as exc:
            print(f"Failed to decompress with gzip: {exc}")
    print("gzip is not available. Please extract the downloaded file manually before creating the BLAST database.")
    return None


def run_command(command, description):
    print(f"Running: {' '.join(command)}")
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False)
    except OSError as exc:
        print(f"Failed to start {description}: {exc}")
        return False, ""

    if result.stdout:
        print(result.stdout.strip())
    if result.stderr:
        print(result.stderr.strip())

    if result.returncode != 0:
        print(f"{description} failed.")
        return False, result.stderr
    return True, result.stdout


def get_database_info(database_name):
    blastdbcmd = find_executable("blastdbcmd")
    if not blastdbcmd:
        return None

    command = [blastdbcmd, "-info", "-db", database_name]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        return None

    info_lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    sequence_count = "Unknown"
    for line in info_lines:
        if "sequences" in line.lower():
            parts = line.split()
            for index, part in enumerate(parts):
                if part.lower() == "sequences;":
                    if index > 0:
                        sequence_count = parts[index - 1]
                        break
            if sequence_count != "Unknown":
                break

    return {"raw_info": info_lines, "sequence_count": sequence_count}


def print_database_summary(database_path, database_type):
    database_name = str(database_path)
    info = get_database_info(database_name)
    print("Database created successfully.")
    print(f"Database name: {database_path.name}")
    print(f"Location: {database_path}")
    print(f"Type: {database_type}")
    if info:
        print(f"Number of sequences: {info['sequence_count']}")
        if info["raw_info"]:
            print("Database info:")
            for line in info["raw_info"]:
                print(f"  {line}")
    else:
        print("Number of sequences: Unknown")
    print()


def download_ncbi_database(title, url, fasta_name, database_type):
    print(f"\n=== {title} ===")
    tools, ready = ensure_download_prerequisites()
    if not ready:
        return

    if title.lower().startswith("ncbi nr"):
        warning = "This database may require hundreds of GB. Continue? (Y/N)"
        if not prompt_yes_no(warning):
            print("Download cancelled.")
            return

    download_folder = prompt_folder_path(
        "Enter a folder to save the downloaded files: ",
        "./blast_databases",
    )
    if not download_folder:
        print("No folder selected. Cancelled.")
        return

    folder = Path(download_folder).expanduser()
    folder.mkdir(parents=True, exist_ok=True)

    archive_path = folder / f"{fasta_name}.gz"
    downloaded_file = download_with_progress(url, archive_path)
    fasta_path = folder / f"{fasta_name}.fasta"
    extract_gzip_file(downloaded_file, fasta_path)

    output_base = folder / fasta_name
    makeblastdb = find_executable("makeblastdb")
    command = [
        makeblastdb,
        "-in",
        str(fasta_path),
        "-dbtype",
        "prot",
        "-parse_seqids",
        "-out",
        str(output_base),
    ]
    success, _ = run_command(command, "makeblastdb")
    if not success:
        print("Database creation failed.")
        return

    print_database_summary(output_base, database_type)


def download_uniprot_database(database_name, url, output_file_name):
    print(f"\n=== {database_name} ===")
    tools, ready = ensure_download_prerequisites()
    if not ready:
        return

    save_path = prompt_folder_path(
        "Enter a file path to save the FASTA file (for example: C:/blastdb/uniprot_sprot.fasta): ",
        "./blast_databases/uniprot_sprot.fasta",
    )
    if not save_path:
        print("No file path selected. Cancelled.")
        return

    output_path = resolve_download_path(save_path, output_file_name)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading FASTA from UniProt REST API")
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            with open(output_path, "wb") as handle:
                shutil.copyfileobj(response, handle)
    except Exception as exc:
        print(f"Download failed: {exc}")
        return

    makeblastdb = find_executable("makeblastdb")
    output_base = output_path.with_suffix("")
    command = [
        makeblastdb,
        "-in",
        str(output_path),
        "-dbtype",
        "prot",
        "-parse_seqids",
        "-out",
        str(output_base),
    ]
    success, _ = run_command(command, "makeblastdb")
    if not success:
        print("Database creation failed.")
        return

    print_database_summary(output_base, "Protein")


def create_custom_blast_database():
    print("\n=== CREATE BLAST DATABASE FROM LOCAL FASTA ===")
    tools, ready = ensure_download_prerequisites()
    if not ready:
        return

    fasta_path_input = prompt_folder_path("Enter the FASTA file path: ")
    if not fasta_path_input:
        print("No FASTA file path provided. Cancelled.")
        return

    fasta_path = Path(fasta_path_input).expanduser()
    if not fasta_path.exists() or not fasta_path.is_file():
        print(f"FASTA file not found: {fasta_path}")
        return

    database_name = prompt_folder_path("Enter the database name (for example: mydb): ")
    if not database_name:
        print("No database name provided. Cancelled.")
        return

    db_type = input("Choose database type: protein or nucleotide: ").strip().lower()
    if db_type not in {"protein", "nucleotide"}:
        print("Invalid database type. Please enter protein or nucleotide.")
        return

    db_type_flag = "prot" if db_type == "protein" else "nucl"
    output_base = Path(database_name).expanduser()
    output_base.parent.mkdir(parents=True, exist_ok=True)

    makeblastdb = find_executable("makeblastdb")
    command = [
        makeblastdb,
        "-in",
        str(fasta_path),
        "-dbtype",
        db_type_flag,
        "-parse_seqids",
        "-out",
        str(output_base),
    ]
    success, _ = run_command(command, "makeblastdb")
    if not success:
        print("Database creation failed.")
        return

    print_database_summary(output_base, "Protein" if db_type_flag == "prot" else "Nucleotide")


def scan_existing_databases():
    print("\n=== SCAN EXISTING DATABASES ===")
    folder_input = prompt_folder_path(
        "Enter a folder to scan for BLAST databases (press Enter to scan a common location): ",
        ".",
    )
    if not folder_input:
        scan_root = Path(".").expanduser()
    else:
        scan_root = Path(folder_input).expanduser()

    if not scan_root.exists() or not scan_root.is_dir():
        print(f"Folder not found: {scan_root}")
        return

    found = []
    for root, _, files in os.walk(scan_root):
        for file_name in files:
            extension = Path(file_name).suffix.lower().lstrip(".")
            if extension in DB_EXTENSIONS["protein"] or extension in DB_EXTENSIONS["nucleotide"]:
                base_name = Path(file_name).stem
                full_path = os.path.join(root, file_name)
                found.append((base_name, full_path))

    if not found:
        print("No BLAST databases were found.")
        return

    print("Found databases:")
    for base_name, full_path in sorted(set(found)):
        print(f"  {base_name}: {full_path}")
    print()


def show_menu():
    print("=== DATABASE MANAGER ===")
    print("1. Download NCBI NR protein database")
    print("2. Download NCBI RefSeq protein database")
    print("3. Download UniProt Swiss-Prot protein database")
    print("4. Download UniProt complete protein database")
    print("5. Create BLAST database from local FASTA file")
    print("6. Scan existing databases")
    print("7. Exit")
    print()


def main():
    print_installation_summary()
    while True:
        show_menu()
        try:
            choice = input("Choose an option (1-7): ").strip()
        except EOFError:
            print("Goodbye.")
            break

        if choice == "1":
            download_ncbi_database(
                "Download NCBI NR protein database",
                "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz",
                "nr",
                "Protein",
            )
        elif choice == "2":
            download_ncbi_database(
                "Download NCBI RefSeq protein database",
                "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/refseq_protein.gz",
                "refseq_protein",
                "Protein",
            )
        elif choice == "3":
            download_uniprot_database(
                "Download UniProt Swiss-Prot protein database",
                "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28reviewed%3Atrue%29",
                "uniprot_swissprot.fasta",
            )
        elif choice == "4":
            download_uniprot_database(
                "Download UniProt complete protein database",
                "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%2A%29",
                "uniprot_complete.fasta",
            )
        elif choice == "5":
            create_custom_blast_database()
        elif choice == "6":
            scan_existing_databases()
        elif choice == "7":
            print("Goodbye.")
            break
        else:
            print("Invalid choice. Please enter a number from 1 to 7.")


if __name__ == "__main__":
    main()
