import os
import shutil
import subprocess
import sys
from pathlib import Path


WINDOWS = sys.platform.startswith("win")
DEFAULT_INSTALL_DIR = r"C:\\BLAST"
COMMON_BLAST_PATHS = [
    r"C:\\Program Files\\NCBI\\blast\\bin",
    r"C:\\Program Files\\NCBI\\blast-2.14.0+\\bin",
    r"C:\\Program Files\\NCBI\\blast-2.13.0+\\bin",
    r"C:\\Program Files\\NCBI\\blast-2.12.0+\\bin",
    r"C:\\Program Files\\NCBI\\blast\\bin",
    r"C:\\Program Files (x86)\\NCBI\\blast\\bin",
]


def run_cmd(command: str, capture_output: bool = True):
    """Run a command and return a completed process."""
    return subprocess.run(
        command,
        shell=True,
        capture_output=capture_output,
        text=True,
    )


def is_windows():
    if not WINDOWS:
        print("This script is intended for Windows systems only.")
        return False
    return True


def find_blast_executable():
    """Try to locate blastn or other BLAST binaries already on the system."""
    exe = shutil.which("blastn")
    if exe:
        return exe

    for base in COMMON_BLAST_PATHS:
        possible = [
            os.path.join(base, "blastn.exe"),
            os.path.join(base, "blastn"),
            os.path.join(base, "update_blastdb.pl"),
        ]
        for path in possible:
            if os.path.exists(path):
                return path
    return None


def install_blast_via_winget():
    """Install BLAST using winget when available."""
    print("Checking for winget...")
    winget = shutil.which("winget")
    if not winget:
        print("winget is not installed or not on PATH.")
        print("Please install App Installer / winget from Microsoft, then rerun this script.")
        return False

    package_ids = [
        "NCBI.NCBIblast",
        "NCBI.BLAST",
        "NCBI.NCBIBLAST",
        "NCBI.Blast",
    ]

    for pkg in package_ids:
        print(f"Trying to install BLAST with package ID: {pkg}")
        result = run_cmd(
            f'winget install --id "{pkg}" --accept-source-agreements --accept-package-agreements --silent'
        )
        if result.returncode == 0:
            print(f"BLAST installed successfully using {pkg}.")
            return True
        print(f"Package failed: {pkg}")
        if result.stderr:
            print(result.stderr.strip())

    return False


def install_blast_via_choco():
    """Fallback: try installing BLAST with Chocolatey."""
    print("Checking for choco...")
    choco = shutil.which("choco")
    if not choco:
        print("Chocolatey was not found.")
        return False

    result = run_cmd(
        'choco install ncbi-blast -y --no-progress',
    )
    if result.returncode == 0:
        print("BLAST installed successfully with Chocolatey.")
        return True
    if result.stderr:
        print(result.stderr.strip())
    return False


def install_blast():
    blast_path = find_blast_executable()
    if blast_path:
        print(f"BLAST is already installed: {blast_path}")
        return True

    print("BLAST was not found on this machine.")
    choice = input("Install BLAST now? [Y/N]: ").strip().lower()
    if choice not in {"y", "yes"}:
        print("BLAST install skipped.")
        return False

    if install_blast_via_winget():
        return True

    if install_blast_via_choco():
        return True

    print("Automatic installation failed. Please install BLAST manually and rerun the script.")
    print("Recommended: install NCBI BLAST and then run this script again.")
    return False


def find_blastdb_update_script():
    """Look for update_blastdb.pl, which is the official way to download BLAST databases."""
    for root in COMMON_BLAST_PATHS:
        script = os.path.join(root, "update_blastdb.pl")
        if os.path.exists(script):
            return script

    for base, _, files in os.walk(r"C:\\"):
        if "update_blastdb.pl" in files:
            return os.path.join(base, "update_blastdb.pl")

    return None


def install_blast_database(db_name: str, install_dir: str):
    """Use the installed BLAST tools to download a common NCBI database."""
    update_script = find_blastdb_update_script()

    if update_script:
        print(f"Found BLAST database updater at: {update_script}")
        out_dir = install_dir
        os.makedirs(out_dir, exist_ok=True)

        # This command downloads the database into the target directory.
        command = f'"{update_script}" -decompress -verbose -nt {db_name} -outdir "{out_dir}"'
        # The above line may not be valid for every BLAST installation.
        # The more standard invocation is below, which also works for common DBs.
        command = f'"{update_script}" --decompress --verbose "{db_name}" -outdir "{out_dir}"'

        print(f"Running: {command}")
        result = run_cmd(command)
        if result.returncode == 0:
            print(f"BLAST database '{db_name}' installed successfully in {out_dir}")
            return True
        if result.stderr:
            print(result.stderr.strip())
        if result.stdout:
            print(result.stdout.strip())
        print("The updater command failed. Trying a fallback method.")

    # Fallback: download a single database file directly from NCBI for common DBs.
    print("Trying direct download fallback...")
    os.makedirs(install_dir, exist_ok=True)

    db_urls = {
        "swissprot": "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz",
        "nr": "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz",
        "nt": "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz",
        "refseq_rna": "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/refseq_rna.gz",
    }

    chosen_url = db_urls.get(db_name.lower())
    if not chosen_url:
        print(f"Database '{db_name}' is not in the predefined download list.")
        return False

    archive_path = os.path.join(install_dir, os.path.basename(chosen_url))
    try:
        print(f"Downloading {chosen_url}")
        import urllib.request

        urllib.request.urlretrieve(chosen_url, archive_path)
        print(f"Downloaded {archive_path}")

        # Decompress with gzip if needed.
        import gzip
        target_path = archive_path.replace(".gz", "")
        with gzip.open(archive_path, "rb") as gz_in:
            with open(target_path, "wb") as out:
                out.write(gz_in.read())

        print(f"Database unpacked to: {target_path}")
        return True
    except Exception as exc:
        print(f"Direct download failed: {exc}")
        return False


def prompt_for_database():
    print("\nAvailable database choices:")
    print("  1) nt")
    print("  2) nr")
    print("  3) swissprot")
    print("  4) refseq_rna")
    print("  5) custom database name")

    answer = input("Choose a database [1-5]: ").strip()
    mapping = {
        "1": "nt",
        "2": "nr",
        "3": "swissprot",
        "4": "refseq_rna",
        "5": "custom",
    }

    db_name = mapping.get(answer)
    if not db_name:
        return "nt"

    if db_name == "custom":
        return input("Enter the BLAST database name or exact NCBI database key: ").strip() or "nt"

    return db_name


def main():
    if not is_windows():
        sys.exit(1)

    print("Windows BLAST installer")
    print("=" * 25)

    if not install_blast():
        print("BLAST could not be installed. Exiting.")
        sys.exit(1)

    install_dir = input(f"Where should the BLAST database be stored? [default: {DEFAULT_INSTALL_DIR}]: ").strip() or DEFAULT_INSTALL_DIR
    install_dir = install_dir.strip('"')
    os.makedirs(install_dir, exist_ok=True)

    db_name = prompt_for_database()

    print(f"Installing BLAST database '{db_name}' into {install_dir}...")
    ok = install_blast_database(db_name, install_dir)

    if ok:
        print("\nInstallation complete.")
        print("You can test the install with: blastn -help")
        print(f"Your database directory is: {install_dir}")
    else:
        print("\nThe BLAST installation finished, but the database download did not complete successfully.")
        print("You can download a database manually from NCBI and save it into the folder above.")


if __name__ == "__main__":
    main()
