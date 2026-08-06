import os
import shutil
import subprocess
import sys
import urllib.request
import gzip
from pathlib import Path


# ==========================
# SETTINGS
# ==========================

WINDOWS = sys.platform.startswith("win")

BLAST_HOME = Path(r"C:\BLAST")

BLAST_BIN_PATHS = [

    r"C:\BLAST\bin",

    r"C:\Program Files\NCBI\blast\bin",

    r"C:\Program Files (x86)\NCBI\blast\bin",

]


DATABASE_FOLDER = BLAST_HOME / "db"



# ==========================
# COMMAND RUNNER
# ==========================

def run_cmd(command):

    return subprocess.run(
        command,
        shell=True,
        capture_output=True,
        text=True
    )



# ==========================
# WINDOWS CHECK
# ==========================

def check_windows():

    if not WINDOWS:

        print("❌ This script only supports Windows")

        return False


    print("✅ Windows detected")

    return True



# ==========================
# FIND BLAST
# ==========================

def find_blast():

    print("\n🔍 Checking BLAST installation...")


    # PATH check

    blastp = shutil.which("blastp")


    if blastp:

        print("✅ BLAST found in PATH:")
        print(blastp)

        return blastp



    # Folder check

    for folder in BLAST_BIN_PATHS:


        exe = Path(folder) / "blastp.exe"


        if exe.exists():

            print("✅ BLAST found:")
            print(exe)

            return str(exe)



    print("\n❌ BLAST not found")

    print("Expected location:")
    print(r"C:\BLAST\bin")


    return None



# ==========================
# INSTALL BLAST
# ==========================

def install_blast_winget():


    print("\nChecking winget...")


    if not shutil.which("winget"):

        print("❌ winget not available")

        return False



    packages = [

        "NCBI.NCBIblast",
        "NCBI.BLAST"

    ]



    for pkg in packages:


        print("Trying:", pkg)


        result = run_cmd(
            f'winget install --id "{pkg}" '
            '--accept-package-agreements '
            '--accept-source-agreements'
        )


        if result.returncode == 0:

            print("✅ BLAST installed")

            return True



    return False




def install_blast_choco():


    print("\nChecking Chocolatey...")


    if not shutil.which("choco"):

        print("❌ Chocolatey not found")

        return False



    result = run_cmd(
        "choco install ncbi-blast -y"
    )


    if result.returncode == 0:

        print("✅ BLAST installed")

        return True


    return False




def install_blast():


    blast = find_blast()


    if blast:

        return True



    print("\nBLAST missing")

    answer = input(
        "Install BLAST now? (Y/N): "
    ).lower()



    if answer not in ["y","yes"]:

        return False



    if install_blast_winget():

        return True



    if install_blast_choco():

        return True



    print(
        "❌ Automatic installation failed"
    )

    print(
        "Install manually from NCBI"
    )

    return False



# ==========================
# CREATE DATABASE FOLDER
# ==========================

def create_database_folder():


    print("\nChecking database folder")


    DATABASE_FOLDER.mkdir(
        parents=True,
        exist_ok=True
    )


    print(
        "Database location:"
    )

    print(
        DATABASE_FOLDER
    )



# ==========================
# CHOOSE DATABASE
# ==========================

def choose_database():


    print("\nAvailable databases")

    print(
        """
1. swissprot
2. nr
3. nt
4. refseq_protein
"""
    )


    choice = input(
        "Choose database: "
    )



    databases = {

        "1":"swissprot",

        "2":"nr",

        "3":"nt",

        "4":"refseq_protein"

    }



    if choice not in databases:

        print("❌ Invalid choice")

        return None



    return databases[choice]



# ==========================
# DOWNLOAD DATABASE
# ==========================

def download_database(db):


    print(
        "\nDownloading:",
        db
    )


    urls = {


    "swissprot":
    "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz",


    "nr":
    "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz",


    "nt":
    "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz",


    "refseq_protein":
    "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/refseq_protein.gz"

    }



    if db not in urls:

        print("❌ Database unavailable")

        return False



    archive = DATABASE_FOLDER / f"{db}.gz"


    print(
        "Downloading file..."
    )


    try:


        urllib.request.urlretrieve(
            urls[db],
            archive
        )


        print(
            "✅ Download complete"
        )



    except Exception as e:


        print(
            "❌ Download failed"
        )

        print(e)

        return False



    return True



# ==========================
# CHECK DATABASE
# ==========================

def check_database(db):


    print(
        "\nChecking database..."
    )


    database = DATABASE_FOLDER / db



    result = run_cmd(

        [
            "blastdbcmd",
            "-db",
            str(database),
            "-info"

        ]

    )


    print(result.stdout)

    print(result.stderr)




# ==========================
# MAIN
# ==========================

def main():


    print(
        "======================"
    )

    print(
        "NCBI BLAST INSTALLER"
    )

    print(
        "======================"
    )



    if not check_windows():

        return



    if not install_blast():

        return



    create_database_folder()



    db = choose_database()



    if db is None:

        return



    if download_database(db):


        print(
            "\nDatabase downloaded"
        )

        print(
            "Location:"
        )

        print(
            DATABASE_FOLDER
        )



if __name__ == "__main__":

    main()