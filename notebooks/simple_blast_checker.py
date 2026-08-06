#!/usr/bin/env python3
"""Standalone Windows script to verify NCBI BLAST+ and database files."""

import os
import subprocess
import shutil
import sys

BLAST_COMMANDS = ["blastn", "blastp", "makeblastdb"]
COMMON_PATHS = [
    r"C:\Program Files\NCBI",
    r"C:\Program Files (x86)\NCBI",
    r"C:\NCBI",
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


def scan_databases(base_path):
    found = {
        "protein": set(),
        "nucleotide": set(),
    }
    for root, _, files in os.walk(base_path):
        for file_name in files:
            extension = os.path.splitext(file_name)[1].lower().lstrip('.')
            if extension in DB_EXTENSIONS["protein"]:
                found["protein"].add(root)
            elif extension in DB_EXTENSIONS["nucleotide"]:
                found["nucleotide"].add(root)
    return found


def format_location(paths):
    return paths if paths else ["Not found"]


def prompt_database_folder():
    try:
        user_input = input("Enter a folder path to scan for BLAST databases (press Enter to skip): ").strip()
    except EOFError:
        return None
    return user_input or None


def main():
    print("=== NCBI BLAST CHECK ===")

    blast_paths = {}
    for cmd in BLAST_COMMANDS:
        exe_path = find_executable(cmd)
        blast_paths[cmd] = exe_path

    installed = any(blast_paths.values())
    print(f"BLAST installed: {'YES' if installed else 'NO'}")

    version = "Not available"
    location_lines = []

    if installed:
        if blast_paths.get("blastn"):
            version = get_blast_version(blast_paths["blastn"])
        elif blast_paths.get("blastp"):
            version = get_blast_version(blast_paths["blastp"])
        elif blast_paths.get("makeblastdb"):
            version = get_blast_version(blast_paths["makeblastdb"])

    print(f"Version: {version}")

    for cmd in BLAST_COMMANDS:
        path = blast_paths.get(cmd)
        if path:
            location_lines.append(f"{cmd}: {path}")
        else:
            location_lines.append(f"{cmd}: Not found")

    print("Location:")
    for line in location_lines:
        print(f"  {line}")

    print("\n=== DATABASE CHECK ===")
    folder = prompt_database_folder()
    if folder:
        if not os.path.isdir(folder):
            print(f"Database folder path is invalid or inaccessible: {folder}")
            print("Database found: NO")
            return
        found = scan_databases(folder)
    else:
        found = {"protein": set(), "nucleotide": set()}
        for base in COMMON_PATHS:
            if os.path.isdir(base):
                dbs = scan_databases(base)
                found["protein"].update(dbs["protein"])
                found["nucleotide"].update(dbs["nucleotide"])

    database_type = "None"
    database_locations = []
    if found["protein"]:
        database_type = "Protein"
        database_locations.extend(sorted(found["protein"]))
    if found["nucleotide"]:
        database_type = "Nucleotide" if database_type == "None" else "Protein + Nucleotide"
        database_locations.extend(sorted(found["nucleotide"]))

    database_found = bool(database_locations)
    print(f"Database found: {'YES' if database_found else 'NO'}")
    print("Database location:")
    if database_locations:
        for loc in sorted(set(database_locations)):
            print(f"  {loc}")
    else:
        print("  None")
    print(f"Database type: {database_type}")


if __name__ == "__main__":
    main()
