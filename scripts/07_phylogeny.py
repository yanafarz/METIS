import logging
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


DEFAULT_HOMOLOGS_PER_TREE = 10
MIN_SEQUENCE_LENGTH = 30

MAFFT_TIMEOUT = 1800      # 30 minutes
TRIMAL_TIMEOUT = 1800     # 30 minutes
IQTREE_TIMEOUT = 7200     # 2 hours
BLASTDBCMD_TIMEOUT = 600   # 10 minutes

UFBBOOT_REPLICATES = 1000
SH_ALRT_REPLICATES = 1000


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)

logger = logging.getLogger(__name__)


def clean_path(path_string: str) -> Path:
    value = path_string.strip().strip('"')

    # Expand ~
    if value.startswith("~"):
        value = str(Path.home()) + value[1:]

    match = re.match(r"^/([a-zA-Z])/(.*)$", value)

    if match:
        drive = match.group(1).upper()
        remainder = match.group(2)
        value = f"{drive}:/{remainder}"

    return Path(value)


def ask_existing_file(prompt: str) -> Path:
    while True:
        value = input(prompt).strip()

        if not value:
            print("ERROR: Path cannot be empty.")
            continue

        path = clean_path(value)

        if not path.is_file():
            print("\nERROR: File not found:")
            print(f"  {path}\n")
            continue

        return path

def ask_output_folder() -> Path:
    """Ask for output folder and create it if necessary."""

    while True:
        value = input(
            "Path to OUTPUT folder (will be created if necessary): "
        ).strip()

        if not value:
            print("ERROR: Output folder path cannot be empty.")
            continue

        path = clean_path(value)

        try:
            path.mkdir(parents=True, exist_ok=True)
            return path
        except Exception as e:
            print(f"\nERROR: Cannot create output folder:\n  {e}\n")


def ask_positive_integer(prompt: str, default: int) -> int:
    """
    Ask for a positive integer.
    """

    while True:
        value = input(f"{prompt} [{default}]: ").strip()

        if not value:
            return default

        try:
            number = int(value)

            if number < 1:
                raise ValueError

            return number

        except ValueError:
            print("ERROR: Please enter a positive integer.")

def windows_to_gitbash_path(path: Path) -> str:
    """
    Convert a Windows path to Git Bash format.
    """

    path = Path(path)
    drive = path.drive

    if not drive:
        return path.as_posix()

    drive_letter = drive[0].lower()
    remainder = path.as_posix()[2:]

    return f"/{drive_letter}{remainder}"

def find_git_bash() -> Path | None:
    """Find Git Bash from the system PATH."""

    bash = shutil.which("bash")

    if bash and Path(bash).is_file():
        return Path(bash)

    return None
    
def test_executable(executable: Path, version_arguments=None) -> bool:
    args = version_arguments or ["--version"]

    print("\n" + "=" * 70)
    print("EXECUTABLE DIAGNOSTIC")
    print("=" * 70)

    # 1 — PATH
    print(f"\n[1/5] Checking path:\n      {executable}")
    path = executable.absolute()

    if not path.is_file():
        print("      FAILED: File does not exist.")
        return False

    print("      OK")

    # 2 — FILE
    try:
        size = path.stat().st_size
        print(
            f"\n[2/5] File check:"
            f"\n      Name: {path.name}"
            f"\n      Size: {size:,} bytes"
            f"\n      OK"
        )
    except Exception as e:
        print(f"      FAILED: {e}")
        return False

    # 3 — EXECUTE
    command = [str(path)] + args
    print(f"\n[3/5] Running:\n      {command}")

    try:
        result = subprocess.run(
            command, capture_output=True, text=True, timeout=30
        )

    except FileNotFoundError:
        print("      FAILED: Executable could not be launched.")
        return False

    except PermissionError:
        print("      FAILED: Permission denied.")
        return False

    except subprocess.TimeoutExpired:
        print("      FAILED: Timeout.")
        return False

    except Exception as e:
        print(f"      FAILED: {type(e).__name__}: {e}")
        return False

    # 4 — RESULT
    print(f"\n[4/5] Result: {result.returncode}")

    if result.stdout.strip():
        print("      STDOUT:")
        print("      " + result.stdout.strip().replace("\n", "\n      "))

    if result.stderr.strip():
        print("      STDERR:")
        print("      " + result.stderr.strip().replace("\n", "\n      "))

    # 5 — DIAGNOSIS
    print("\n[5/5] FINAL DIAGNOSIS")

    if result.returncode == 0:
        print("      EXECUTABLE WORKS.")
        print("=" * 70)
        return True

    error = result.stderr.lower()

    if "usage" in error or "version" in error:
        print("      EXECUTABLE LAUNCHED.")
        print("      Argument/version check returned an error.")
        print("=" * 70)
        return True

    print("      EXECUTABLE FAILED.")
    print("=" * 70)
    return False

def test_trimal_executable(
    trimal: Path
) -> tuple[bool, Path | None]:
    """Test Windows trimAl through Git Bash."""

    print("\nTesting trimAl through Git Bash...")
    print(f"  Executable: {trimal}")

    if not trimal.is_file():
        print("  trimAl: FAILED — file does not exist.")
        return False, None

    git_bash = find_git_bash()

    if git_bash is None:
        print("  trimAl: FAILED — Git Bash not found.")
        return False, None

    print(f"  Git Bash: {git_bash}")

    trimal_dir = windows_to_gitbash_path(trimal.parent)
    command = f'cd "{trimal_dir}" && ./trimal.exe --version'

    try:
        result = subprocess.run(
            [str(git_bash), "-lc", command],
            capture_output=True,
            text=True,
            timeout=30
        )

    except subprocess.TimeoutExpired:
        print("  trimAl: FAILED — timeout.")
        return False, None

    except Exception as e:
        print(f"  trimAl: FAILED — {e}")
        return False, None

    output = (
        (result.stdout or "") +
        (result.stderr or "")
    ).strip()

    if result.returncode == 0:
        print("  trimAl: WORKING")

        if output:
            print(output)

        return True, git_bash

    print(f"  trimAl: FAILED — return code {result.returncode}")

    if output:
        print(output)

    return False, None


def test_trimal_wsl(wsl_path: str) -> bool:
    print("=" * 70)
    print("WSL trimAl EXECUTABLE DIAGNOSTIC")
    print("=" * 70)

    print("\n[1/4] Checking WSL and trimAl file...")

    try:
        result = subprocess.run(
            [
                "wsl.exe",
                "sh",
                "-c",
                f'test -f "{wsl_path}" && test -x "{wsl_path}"'
            ],
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode != 0:
            print("      FAILED: WSL or trimAl file check.")
            return False

        print("      WSL: WORKING")
        print("      trimAl: EXISTS + EXECUTABLE")

    except Exception as e:
        print(f"      FAILED: {e}")
        return False

    print("\n[2/4] Checking trimAl dependencies...")

    result = subprocess.run(
        [
            "wsl.exe",
            "sh",
            "-c",
            f'ldd "{wsl_path}" 2>&1'
        ],
        capture_output=True,
        text=True,
        timeout=30
    )

    if "not found" in result.stdout.lower():
        print("      FAILED: Missing dependency.")
        print(result.stdout)
        return False

    print("      Dependencies: OK")

    print("\n[3/4] Checking WSL environment...")

    result = subprocess.run(
        [
            "wsl.exe",
            "sh",
            "-c",
            "uname -m"
        ],
        capture_output=True,
        text=True,
        timeout=30
    )

    print(f"      Architecture: {result.stdout.strip()}")

    print("\n[4/4] Running trimAl...")

    try:
        result = subprocess.run(
            [
                "wsl.exe",
                "sh",
                "-c",
                f'"{wsl_path}" -h'
            ],
            capture_output=True,
            text=True,
            timeout=30
        )

    except subprocess.TimeoutExpired:
        print("      FAILED: trimAl timed out.")
        return False

    output = (
        (result.stdout or "") +
        (result.stderr or "")
    ).strip()

    if "trimal" not in output.lower():
        print("      FAILED: No recognizable trimAl output.")
        print(output[:1000])
        return False

    print("      trimAl: WORKING")
    print(f"      Return code: {result.returncode}")

    print("\n" + "=" * 70)
    print("FINAL DIAGNOSIS")
    print("trimAl WORKS INSIDE WSL.")
    print("=" * 70)

    return True

def ask_executable(
    tool_name: str
) -> tuple[Path, Path | None]:

    while True:

        value = input(
            f"Path to {tool_name} executable: "
        ).strip()

        if not value:

            print(
                "ERROR: Path cannot be empty."
            )

            continue

        path = clean_path(value)

        if tool_name.lower() == "trimal":

            original_value = value.strip().strip('"')

            if original_value.startswith("/"):

                print()
                print("=" * 70)
                print("WSL trimAl path detected")
                print("=" * 70)

                print()
                print("Path:")
                print(f"  {original_value}")

                print()
                print("Testing trimAl INSIDE WSL...")
                print()

                success = test_trimal_wsl(original_value)

                if success:

                    print()
                    print("=" * 70)
                    print("FINAL DIAGNOSIS")
                    print()
                    print("        trimAl WORKS INSIDE WSL.")
                    print("=" * 70)

                    return Path(original_value), None

                print()
                print("=" * 70)
                print("FINAL DIAGNOSIS")
                print()
                print("        trimAl FAILED INSIDE WSL.")
                print("=" * 70)

                print()
                print(
                    "The supplied WSL trimAl executable "
                    "could not be run."
                )

                print(
                    "Please provide another path."
                )

                print()

                continue

        if not path.is_file():
            print("\nERROR: Executable not found:")
            print(f"  {path}\n")
            continue

        if tool_name.lower() == "trimal":
            success, git_bash = test_trimal_executable(path)

            if success:
                return path, git_bash

            print("\nThe supplied trimAl executable could not be run.")
            print("Please provide another path.\n")
            continue

        print(f"\nTesting {tool_name}...")
        print(f"  Executable: {path}")

        if tool_name.lower() == "blastdbcmd":
            test_ok = test_executable(path, ["-version"])
        else:
            test_ok = test_executable(path)

        if test_ok:
            print(f"  {tool_name}: WORKING")
            return path, None

        print(f"  {tool_name}: FAILED")
        print(f"\nThe supplied {tool_name} executable could not be run.")
        print("Please provide another path.\n")

def detect_blast_database(db_root: Path):
    print("\nChecking BLAST protein database...")
    print(f"  Database root: {db_root}")

    single_pin = Path(str(db_root) + ".pin")
    single_phr = Path(str(db_root) + ".phr")
    single_psq = Path(str(db_root) + ".psq")

    if single_pin.is_file() and single_phr.is_file() and single_psq.is_file():
        print("\n  COMPLETE SINGLE-VOLUME DATABASE FOUND")
        print(f"    {single_pin.name}")
        print(f"    {single_phr.name}")
        print(f"    {single_psq.name}")

        return {
            "type": "single",
            "volumes": [db_root]
        }

    parent = db_root.parent
    prefix = db_root.name

    volumes = []

    if parent.exists():
        pattern = re.compile(
            re.escape(prefix) + r"\.(\d+)\.pin$",
            re.IGNORECASE
        )

        for pin_file in parent.iterdir():
            match = pattern.match(pin_file.name)

            if not match:
                continue

            volume_number = int(match.group(1))

            volume_root = (
                parent
                / f"{prefix}.{volume_number:02d}"
            )

            phr = Path(str(volume_root) + ".phr")
            psq = Path(str(volume_root) + ".psq")

            if phr.is_file() and psq.is_file():
                volumes.append(
                    (volume_number, volume_root)
                )

    volumes.sort(key=lambda x: x[0])

    if volumes:
        print()
        print("  SPLIT BLAST PROTEIN DATABASE DETECTED")
        print(f"  Number of volumes: {len(volumes)}")

        for number, root in volumes:
            print(
                f"    Volume {number:02d}: "
                f"{root.name}"
            )

        print()

        return {
            "type": "split",
            "volumes": [
                root
                for number, root in volumes
            ]
        }

    possible_files = (
        list(parent.glob(prefix + ".*"))
        if parent.exists()
        else []
    )

    if possible_files:
        print()
        print(
            "WARNING: Database files were found "
            "but the database appears incomplete."
        )

        print()

        for file in sorted(possible_files):
            print(f"  {file.name}")

    else:
        print()
        print("No BLAST database files found.")

    return None

def ask_blast_database():

    while True:

        value = input(
            "Path to BLAST protein database ROOT "
            "(without .pin/.phr/.psq): "
        ).strip()

        if not value:

            print(
                "ERROR: Database path "
                "cannot be empty."
            )

            continue

        db_root = clean_path(value)

        database_info = (
            detect_blast_database(
                db_root
            )
        )

        if database_info is not None:

            return (
                db_root,
                database_info
            )

        print()
        print(
            "ERROR: Could not detect a "
            "complete protein BLAST database."
        )

        print()

def load_fasta(fasta_file: Path) -> dict:
    """
    Load FASTA into:

        {header: sequence}
    """

    sequences = {}
    current_header = None
    current_sequence = []

    with open(
        fasta_file,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = "".join(current_sequence)

                current_header = line[1:].strip()
                current_sequence = []

            else:
                current_sequence.append(line)

        if current_header is not None:
            sequences[current_header] = "".join(current_sequence)

    return sequences


def get_fasta_id(header: str) -> str:
    return header.split()[0]


def clean_protein_sequence(sequence: str) -> str:
    """
    Clean protein sequence.

    Removes whitespace and terminal *.
    """

    sequence = (
        str(sequence)
        .replace(" ", "")
        .replace("\r", "")
        .replace("\n", "")
        .strip()
        .rstrip("*")
        .upper()
    )

    return sequence


def build_fasta_lookup(fasta_file: Path):
    """
    Build lookup using both:

        full header
        first FASTA ID
    """

    raw = load_fasta(fasta_file)
    lookup = {}

    for header, sequence in raw.items():
        fasta_id = get_fasta_id(header)
        sequence = clean_protein_sequence(sequence)

        lookup[fasta_id] = sequence
        lookup[header] = sequence

    logger.info(
        f"Loaded {len(raw)} Metisa "
        f"protein sequences"
    )

    return lookup

def load_blast_tsv(blast_tsv: Path) -> pd.DataFrame:
    logger.info(f"Loading BLAST TSV: {blast_tsv}")

    df = pd.read_csv(blast_tsv, sep="\t", dtype=str)

    df.columns = [str(c).strip() for c in df.columns]

    required = ["qseqid", "sseqid"]
    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(
            "BLAST TSV is missing required "
            f"columns: {missing}"
        )

    for column in ["pident", "length", "qlen", "slen", "qcovs", "evalue", "bitscore"]:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")

    df["qseqid"] = (
        df["qseqid"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    df["sseqid"] = (
        df["sseqid"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    df = df[
        (df["qseqid"] != "")
        &
        (df["sseqid"] != "")
    ].copy()

    logger.info(
        f"Loaded {len(df)} BLAST rows"
    )

    logger.info(
        f"Unique Metisa proteins: "
        f"{df['qseqid'].nunique()}"
    )

    return df

def extract_accession(
    sseqid: str
) -> str:

    if pd.isna(sseqid):

        return ""

    value = str(
        sseqid
    ).strip()

    if not value:

        return ""

    if "|" in value:

        parts = value.split("|")

        if len(parts) >= 2:

            accession = (
                parts[1]
                .strip()
            )

            if accession:

                return accession

    return value.split()[0]

def select_top_homologs(
    blast_df: pd.DataFrame,
    max_hits: int = DEFAULT_HOMOLOGS_PER_TREE
) -> pd.DataFrame:

    df = blast_df.copy()

    df["accession"] = (
        df["sseqid"]
        .apply(extract_accession)
    )

    df = df[
        df["accession"].str.len() > 0
    ].copy()

    df = df.drop_duplicates(
        subset=[
            "qseqid",
            "accession"
        ],
        keep="first"
    )

    sort_columns = []
    ascending = []

    if "bitscore" in df.columns:

        sort_columns.append(
            "bitscore"
        )

        ascending.append(False)

    if "qcovs" in df.columns:

        sort_columns.append(
            "qcovs"
        )

        ascending.append(False)

    if "pident" in df.columns:

        sort_columns.append(
            "pident"
        )

        ascending.append(False)

    if "evalue" in df.columns:

        sort_columns.append(
            "evalue"
        )

        ascending.append(True)

    if sort_columns:

        df = df.sort_values(
            sort_columns,
            ascending=ascending,
            na_position="last"
        )

    selected = (
        df
        .groupby(
            "qseqid",
            sort=False,
            group_keys=False
        )
        .head(max_hits)
        .reset_index(drop=True)
    )

    selected["homolog_rank"] = (
        selected
        .groupby("qseqid")
        .cumcount()
        + 1
    )

    return selected

def parse_fasta_records(fasta_text: str):
    """
    Parse FASTA text into:

        [(header, sequence), ...]
    """

    records = []
    current_header = None
    current_sequence = []

    for line in fasta_text.splitlines():
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            if current_header is not None:
                sequence = clean_protein_sequence(
                    "".join(current_sequence)
                )

                records.append(
                    (current_header, sequence)
                )

            current_header = line[1:].strip()
            current_sequence = []

        else:
            current_sequence.append(line)

    if current_header is not None:
        sequence = clean_protein_sequence(
            "".join(current_sequence)
        )

        records.append(
            (current_header, sequence)
        )

    return records

def run_blastdbcmd_batch(
    blastdbcmd: Path,
    db_root: Path,
    accessions: list[str],
    output_fasta: Path
) -> bool:

    if not accessions:

        return False

    accession_file = (
        output_fasta.parent
        /
        (
            output_fasta.stem
            +
            "_accessions.txt"
        )
    )

    try:

        accession_file.write_text(
            "\n".join(
                accessions
            )
            +
            "\n",
            encoding="utf-8"
        )

        command = [

            str(blastdbcmd),

            "-db",
            str(db_root),

            "-entry_batch",
            str(accession_file),

            "-outfmt",
            "%f",

            "-out",
            str(output_fasta)
        ]

        logger.info(
            f"Retrieving {len(accessions)} "
            f"sequences with blastdbcmd"
        )

        result = subprocess.run(

            command,

            cwd=str(
                blastdbcmd.parent
            ),

            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,

            text=True,

            timeout=BLASTDBCMD_TIMEOUT
        )

        if result.returncode != 0:

            logger.error(
                "blastdbcmd failed:"
            )

            if result.stderr:

                logger.error(
                    result.stderr.strip()
                )

            return False

        if not output_fasta.is_file():

            logger.error(
                "blastdbcmd did not "
                "produce an output FASTA."
            )

            return False

        if output_fasta.stat().st_size == 0:

            logger.error(
                "blastdbcmd produced "
                "an empty FASTA."
            )

            return False

        return True

    except subprocess.TimeoutExpired:

        logger.error(
            "blastdbcmd timed out."
        )

        return False

    except Exception as e:

        logger.error(
            f"blastdbcmd error: {e}"
        )

        return False

    finally:

        accession_file.unlink(
            missing_ok=True
        )

def safe_filename(
    value: str
) -> str:

    value = str(value)

    value = re.sub(
        r'[<>:"/\\|?*]',
        "_",
        value
    )

    value = value.strip()

    if not value:

        value = "unknown"

    return value

def retrieve_homolog_sequences(
    blastdbcmd: Path,
    db_root: Path,
    candidate_df: pd.DataFrame,
    homolog_dir: Path
):


    homolog_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    retrieval_records = []

    for gene_id, group in candidate_df.groupby(
        "qseqid",
        sort=False
    ):

        logger.info("")
        logger.info(
            f"{gene_id}: retrieving "
            f"{len(group)} homologs"
        )

        safe_gene = safe_filename(
            gene_id
        )

        output_fasta = (
            homolog_dir
            /
            f"{safe_gene}_homologs.fasta"
        )

        accessions = (
            group["accession"]
            .dropna()
            .astype(str)
            .tolist()
        )

        accessions = list(
            dict.fromkeys(
                accessions
            )
        )

        if not accessions:

            logger.warning(
                f"{gene_id}: no accessions"
            )

            continue

        success = run_blastdbcmd_batch(
            blastdbcmd,
            db_root,
            accessions,
            output_fasta
        )

        if not success:

            for accession in accessions:

                retrieval_records.append({

                    "qseqid":
                        gene_id,

                    "accession":
                        accession,

                    "retrieved":
                        False
                })

            continue

        try:

            text = (
                output_fasta
                .read_text(
                    encoding="utf-8",
                    errors="replace"
                )
            )

            records = parse_fasta_records(
                text
            )

        except Exception as e:

            logger.error(
                f"{gene_id}: could not "
                f"parse retrieved FASTA: {e}"
            )

            records = []

        retrieved_accessions = set()

        for header, sequence in records:

            accession = extract_accession(
                header
            )

            if accession:

                retrieved_accessions.add(
                    accession
                )

                retrieval_records.append({

                    "qseqid":
                        gene_id,

                    "accession":
                        accession,

                    "retrieved":
                        True,

                    "sequence_length":
                        len(sequence)
                })

        for accession in accessions:

            if accession not in retrieved_accessions:

                retrieval_records.append({

                    "qseqid":
                        gene_id,

                    "accession":
                        accession,

                    "retrieved":
                        False,

                    "sequence_length":
                        0
                })

        logger.info(
            f"{gene_id}: retrieved "
            f"{len(records)} sequences"
        )

    return pd.DataFrame(
        retrieval_records
    )

def filter_homolog_fasta(
    homolog_fasta: Path,
    min_length: int = MIN_SEQUENCE_LENGTH
):


    if not homolog_fasta.is_file():

        return "", 0

    text = (
        homolog_fasta
        .read_text(
            encoding="utf-8",
            errors="replace"
        )
    )

    records = parse_fasta_records(
        text
    )

    valid_records = []

    for header, sequence in records:

        if len(sequence) < min_length:

            logger.warning(
                f"Removing short homolog "
                f"{header} "
                f"({len(sequence)} aa)"
            )

            continue

        allowed = set(
            "ACDEFGHIKLMNPQRSTVWYBXZJUO"
        )

        invalid = set(sequence) - allowed

        if invalid:

            logger.warning(
                f"Removing sequence "
                f"{header}: invalid "
                f"characters {invalid}"
            )

            continue

        valid_records.append(
            (
                header,
                sequence
            )
        )

    if not valid_records:

        return "", 0

    output_lines = []

    for header, sequence in valid_records:

        output_lines.append(
            f">{header}"
        )

        # Wrap sequence at 80 characters

        for i in range(
            0,
            len(sequence),
            80
        ):

            output_lines.append(
                sequence[
                    i:i + 80
                ]
            )

    return (
        "\n".join(output_lines)
        +
        "\n",
        len(valid_records)
    )

def build_phylo_fasta(
    gene_id: str,
    homolog_fasta: Path,
    metisa_lookup: dict,
    protein_id: str,
    output_fasta: Path
) -> tuple[bool, int]:

    if not homolog_fasta.is_file():

        return False, 0

    metisa_sequence = None

    if protein_id in metisa_lookup:

        metisa_sequence = (
            metisa_lookup[
                protein_id
            ]
        )

    elif gene_id in metisa_lookup:

        metisa_sequence = (
            metisa_lookup[
                gene_id
            ]
        )

    if not metisa_sequence:

        logger.warning(
            f"{gene_id}: Metisa protein "
            f"{protein_id} not found."
        )

        return False, 0

    metisa_sequence = clean_protein_sequence(
        metisa_sequence
    )

    if len(metisa_sequence) < MIN_SEQUENCE_LENGTH:

        logger.warning(
            f"{gene_id}: Metisa sequence "
            f"is too short "
            f"({len(metisa_sequence)} aa)."
        )

        return False, 0

    homolog_text, homolog_count = (
        filter_homolog_fasta(
            homolog_fasta
        )
    )

    if homolog_count == 0:

        logger.warning(
            f"{gene_id}: no valid homolog "
            f"sequences available."
        )

        return False, 0

    with open(
        output_fasta,
        "w",
        encoding="utf-8"
    ) as handle:

        handle.write(
            f">METISA_PLANA|{protein_id}\n"
        )

        for i in range(
            0,
            len(metisa_sequence),
            80
        ):

            handle.write(
                metisa_sequence[
                    i:i + 80
                ]
                +
                "\n"
            )

        handle.write(
            homolog_text
        )

    total_sequences = (
        homolog_count + 1
    )

    logger.info(
        f"{gene_id}: phylogenetic FASTA "
        f"contains {total_sequences} sequences "
        f"(1 Metisa + {homolog_count} homologs)"
    )

    return True, total_sequences

def run_mafft(
    mafft: Path,
    input_fasta: Path,
    output_alignment: Path
) -> bool:

    command = [

        str(mafft),

        "--auto",

        str(input_fasta)
    ]

    logger.info(
        f"Running MAFFT: "
        f"{input_fasta.name}"
    )

    try:

        with open(
            output_alignment,
            "w",
            encoding="utf-8"
        ) as output_handle:

            result = subprocess.run(

                command,

                cwd=str(
                    mafft.parent
                ),

                stdout=output_handle,

                stderr=subprocess.PIPE,

                text=True,

                timeout=MAFFT_TIMEOUT
            )

        if result.returncode != 0:

            logger.error(
                "MAFFT failed:"
            )

            logger.error(
                result.stderr[-5000:]
            )

            return False

        if not output_alignment.is_file():

            return False

        if output_alignment.stat().st_size == 0:

            return False

        return True

    except subprocess.TimeoutExpired:

        logger.error(
            "MAFFT timed out."
        )

        return False

    except Exception as e:

        logger.error(
            f"MAFFT error: {e}"
        )

        return False

def windows_to_wsl_path(path: Path) -> str:
    """
    Convert a Windows path to a WSL /mnt/<drive>/... path.
    """

    path_str = str(path)

    # Normalize Windows backslashes
    path_str = path_str.replace("\\", "/")

    if len(path_str) >= 2 and path_str[1] == ":":
        drive = path_str[0].lower()
        remainder = path_str[2:].lstrip("/")

        return f"/mnt/{drive}/{remainder}"

    return path_str

def run_trimal(
    trimal: Path,
    git_bash: Path | None,
    input_alignment: Path,
    output_alignment: Path
) -> bool:
    """
    Run trimAl.

    Supports:
        - WSL trimAl
        - Native Windows trimAl

    WSL trimAl is executed through:
        wsl <trimal_path>

    Git Bash is not required.
    """

    logger.info(
        f"Running trimAl: "
        f"{input_alignment.name}"
    )

    # ========================================================
    # DETERMINE IF THIS IS THE WSL TRIMAL
    # ========================================================

    trimal_str = str(trimal)

    if (
        trimal_str.startswith("/")
        or trimal_str.startswith("\\home\\")
        or trimal_str.startswith("\\usr\\")
    ):

        # ----------------------------------------------------
        # Convert Python's path representation back to WSL
        # ----------------------------------------------------

        trimal_wsl = trimal_str.replace("\\", "/")

        # Make sure it starts with /
        if not trimal_wsl.startswith("/"):
            trimal_wsl = "/" + trimal_wsl

        # ----------------------------------------------------
        # Convert Windows alignment paths to WSL paths
        # ----------------------------------------------------

        input_wsl = windows_to_wsl_path(
            input_alignment
        )

        output_wsl = windows_to_wsl_path(
            output_alignment
        )

        command = [
            "wsl",
            trimal_wsl,
            "-in",
            input_wsl,
            "-out",
            output_wsl,
            "-automated1"
        ]

        logger.info(
            f"Using WSL trimAl: {trimal_wsl}"
        )

    # ========================================================
    # NATIVE WINDOWS TRIMAL
    # ========================================================

    else:

        command = [
            str(trimal),
            "-in",
            str(input_alignment),
            "-out",
            str(output_alignment),
            "-automated1"
        ]

        logger.info(
            f"Using Windows trimAl: {trimal}"
        )

    # ========================================================
    # RUN TRIMAL
    # ========================================================

    try:

        result = subprocess.run(

            command,

            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,

            text=True,

            timeout=TRIMAL_TIMEOUT
        )

        # ----------------------------------------------------
        # Check exit status
        # ----------------------------------------------------

        if result.returncode != 0:

            logger.error(
                "trimAl failed."
            )

            if result.stdout:

                logger.error(
                    result.stdout
                )

            if result.stderr:

                logger.error(
                    result.stderr
                )

            return False

        # ----------------------------------------------------
        # Check output file
        # ----------------------------------------------------

        if not output_alignment.is_file():

            logger.error(
                "trimAl did not produce "
                "the output alignment."
            )

            return False

        if output_alignment.stat().st_size == 0:

            logger.error(
                "trimAl produced an "
                "empty alignment."
            )

            return False

        logger.info(
            f"trimAl completed successfully: "
            f"{output_alignment.name}"
        )

        return True

    except subprocess.TimeoutExpired:

        logger.error(
            "trimAl timed out."
        )

        return False

    except FileNotFoundError:

        logger.error(
            f"trimAl executable not found: "
            f"{trimal}"
        )

        return False

    except Exception as e:

        logger.error(
            f"trimAl error: {e}"
        )

        return False


# ============================================================
# IQ-TREE
# ============================================================

def run_iqtree(
    iqtree: Path,
    input_alignment: Path,
    output_prefix: Path
) -> bool:

    command = [

        str(iqtree),

        "-s",
        str(input_alignment),

        "-m",
        "MFP",

        "-B",
        str(UFBBOOT_REPLICATES),

        "-alrt",
        str(SH_ALRT_REPLICATES),

        "-nt",
        "AUTO",

        "-pre",
        str(output_prefix)
    ]

    logger.info(
        f"Running IQ-TREE: "
        f"{input_alignment.name}"
    )

    logger.info(
        "Model: ModelFinder Plus (MFP)"
    )

    logger.info(
        f"UFBoot: {UFBBOOT_REPLICATES}"
    )

    logger.info(
        f"SH-aLRT: {SH_ALRT_REPLICATES}"
    )

    try:

        result = subprocess.run(

            command,

            cwd=str(
                iqtree.parent
            ),

            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,

            text=True,

            timeout=IQTREE_TIMEOUT
        )

        if result.returncode != 0:

            logger.error(
                "IQ-TREE failed:"
            )

            if result.stdout:

                logger.error(
                    result.stdout[-5000:]
                )

            if result.stderr:

                logger.error(
                    result.stderr[-5000:]
                )

            return False

        treefile = Path(
            str(output_prefix)
            +
            ".treefile"
        )

        iqtree_report = Path(
            str(output_prefix)
            +
            ".iqtree"
        )

        if not treefile.is_file():

            logger.error(
                "IQ-TREE finished but "
                ".treefile was not produced."
            )

            return False

        if treefile.stat().st_size == 0:

            logger.error(
                "IQ-TREE produced "
                "an empty treefile."
            )

            return False

        if not iqtree_report.is_file():

            logger.warning(
                "IQ-TREE tree was produced "
                "but .iqtree report was not found."
            )

        return True

    except subprocess.TimeoutExpired:

        logger.error(
            "IQ-TREE timed out."
        )

        return False

    except Exception as e:

        logger.error(
            f"IQ-TREE error: {e}"
        )

        return False


# ============================================================
# READ IQ-TREE REPORT
# ============================================================

def extract_iqtree_model(
    iqtree_report: Path
) -> str:
    """
    Try to extract the selected model
    from the IQ-TREE report.
    """

    if not iqtree_report.is_file():

        return ""

    try:

        text = (
            iqtree_report
            .read_text(
                encoding="utf-8",
                errors="replace"
            )
        )

    except Exception:

        return ""

    patterns = [

        r"Best-fit model:\s*(\S+)",

        r"Best-fit model according to BIC:\s*(\S+)",

        r"Best-fit model according to AICc:\s*(\S+)",

        r"Model of substitution:\s*(\S+)"
    ]

    for pattern in patterns:

        match = re.search(
            pattern,
            text,
            flags=re.IGNORECASE
        )

        if match:

            return match.group(1)

    return ""


# ============================================================
# TREE SUPPORT
# ============================================================

def extract_support_values(
    treefile: Path
):
    """
    Extract support values from IQ-TREE
    Newick labels.

    IQ-TREE commonly writes branch support
    as:

        SH-aLRT/UFBoot

    e.g.

        98.5/100

    This function returns both components
    when available.
    """

    if not treefile.is_file():

        return {
            "ufboot_values": [],
            "sh_alrt_values": []
        }

    try:

        tree = (
            treefile
            .read_text(
                encoding="utf-8",
                errors="replace"
            )
            .strip()
        )

    except Exception:

        return {
            "ufboot_values": [],
            "sh_alrt_values": []
        }

    # --------------------------------------------------------
    # Match support labels such as:
    #
    # 98/100
    # 98.5/100
    # --------------------------------------------------------

    pair_matches = re.findall(
        r"\)(\d+(?:\.\d+)?)\/"
        r"(\d+(?:\.\d+)?)"
        r"(?=:|,|\))",
        tree
    )

    sh_alrt_values = []
    ufboot_values = []

    for sh_value, ufboot_value in pair_matches:

        try:

            sh_alrt_values.append(
                float(sh_value)
            )

            ufboot_values.append(
                float(ufboot_value)
            )

        except ValueError:

            continue

    # --------------------------------------------------------
    # If only one support value exists
    # --------------------------------------------------------

    if not pair_matches:

        single_matches = re.findall(
            r"\)(\d+(?:\.\d+)?)(?=:|,|\))",
            tree
        )

        for value in single_matches:

            try:

                ufboot_values.append(
                    float(value)
                )

            except ValueError:

                continue

    return {
        "ufboot_values":
            ufboot_values,

        "sh_alrt_values":
            sh_alrt_values
    }


def summarize_values(values):
    if not values:
        return {"mean": None, "minimum": None, "maximum": None, "nodes": 0}

    return {
        "mean": sum(values) / len(values),
        "minimum": min(values),
        "maximum": max(values),
        "nodes": len(values)
    }


# ============================================================
# PROCESS ONE CANDIDATE
# ============================================================

def process_candidate(
    gene_id: str,
    protein_id: str,
    homolog_fasta: Path,
    metisa_lookup: dict,
    mafft: Path,
    trimal: Path,
    git_bash: Path,
    iqtree: Path,
    phylo_fasta_dir: Path,
    alignments_dir: Path,
    trees_dir: Path
):

    logger.info("")
    logger.info(
        "=" * 75
    )

    logger.info(
        f"PROCESSING CANDIDATE: {gene_id}"
    )

    logger.info(
        "=" * 75
    )

    result = {

        "gene_id":
            gene_id,

        "protein_id":
            protein_id,

        "success":
            False,

        "sequence_count":
            0,

        "alignment_length":
            None,

        "trimmed_alignment_length":
            None,

        "model":
            "",

        "ufboot_mean_pct":
            None,

        "ufboot_min_pct":
            None,

        "ufboot_max_pct":
            None,

        "ufboot_nodes":
            0,

        "sh_alrt_mean_pct":
            None,

        "sh_alrt_min_pct":
            None,

        "sh_alrt_max_pct":
            None,

        "sh_alrt_nodes":
            0,

        "error":
            ""
    }

    safe_gene = safe_filename(
        gene_id
    )

    combined_fasta = (
        phylo_fasta_dir
        /
        f"{safe_gene}_phylo.fasta"
    )

    alignment = (
        alignments_dir
        /
        f"{safe_gene}.aln"
    )

    trimmed = (
        alignments_dir
        /
        f"{safe_gene}.trimmed.aln"
    )

    tree_prefix = (
        trees_dir
        /
        safe_gene
    )

    treefile = Path(
        str(tree_prefix)
        +
        ".treefile"
    )

    iqtree_report = Path(
        str(tree_prefix)
        +
        ".iqtree"
    )

    fasta_success, sequence_count = (
        build_phylo_fasta(

            gene_id,

            homolog_fasta,

            metisa_lookup,

            protein_id,

            combined_fasta
        )
    )

    if not fasta_success:

        result["error"] = (
            "Could not build "
            "phylogenetic FASTA"
        )

        return result

    result[
        "sequence_count"
    ] = sequence_count


    if not alignment.is_file():

        if not run_mafft(

            mafft,

            combined_fasta,

            alignment

        ):

            result["error"] = (
                "MAFFT failed"
            )

            return result

    else:

        logger.info(
            f"MAFFT output already exists: "
            f"{alignment.name}"
        )

    try:

        alignment_text = (
            alignment
            .read_text(
                encoding="utf-8",
                errors="replace"
            )
        )

        alignment_records = (
            parse_fasta_records(
                alignment_text
            )
        )

        if alignment_records:

            result[
                "alignment_length"
            ] = len(
                alignment_records[0][1]
            )

    except Exception:

        pass

    if not trimmed.is_file():

        if not run_trimal(

            trimal,

            git_bash,

            alignment,

            trimmed

        ):

            result["error"] = (
                "trimAl failed"
            )

            return result

    else:

        logger.info(
            f"trimAl output already exists: "
            f"{trimmed.name}"
        )

    try:

        trimmed_text = (
            trimmed
            .read_text(
                encoding="utf-8",
                errors="replace"
            )
        )

        trimmed_records = (
            parse_fasta_records(
                trimmed_text
            )
        )

        if trimmed_records:

            result[
                "trimmed_alignment_length"
            ] = len(
                trimmed_records[0][1]
            )

    except Exception:

        pass

    if not treefile.is_file():

        if not run_iqtree(

            iqtree,

            trimmed,

            tree_prefix

        ):

            result["error"] = (
                "IQ-TREE failed"
            )

            return result

    else:

        logger.info(
            f"IQ-TREE tree already exists: "
            f"{treefile.name}"
        )

    result[
        "model"
    ] = extract_iqtree_model(
        iqtree_report
    )

    support = extract_support_values(
        treefile
    )

    ufboot_summary = summarize_values(
        support["ufboot_values"]
    )

    sh_alrt_summary = summarize_values(
        support["sh_alrt_values"]
    )

    result[
        "ufboot_mean_pct"
    ] = ufboot_summary["mean"]

    result[
        "ufboot_min_pct"
    ] = ufboot_summary["minimum"]

    result[
        "ufboot_max_pct"
    ] = ufboot_summary["maximum"]

    result[
        "ufboot_nodes"
    ] = ufboot_summary["nodes"]

    result[
        "sh_alrt_mean_pct"
    ] = sh_alrt_summary["mean"]

    result[
        "sh_alrt_min_pct"
    ] = sh_alrt_summary["minimum"]

    result[
        "sh_alrt_max_pct"
    ] = sh_alrt_summary["maximum"]

    result[
        "sh_alrt_nodes"
    ] = sh_alrt_summary["nodes"]

    result["success"] = True

    logger.info(
        f"{gene_id}: phylogenetic analysis complete"
    )

    if result[
        "model"
    ]:

        logger.info(
            f"{gene_id}: model = "
            f"{result['model']}"
        )

    if result[
        "ufboot_mean_pct"
    ] is not None:

        logger.info(
            f"{gene_id}: mean UFBoot = "
            f"{result['ufboot_mean_pct']:.2f}%"
        )

    if result[
        "sh_alrt_mean_pct"
    ] is not None:

        logger.info(
            f"{gene_id}: mean SH-aLRT = "
            f"{result['sh_alrt_mean_pct']:.2f}%"
        )

    return result

def load_interpro(
    interpro_file: Path
) -> pd.DataFrame:

    try:

        df = pd.read_csv(

            interpro_file,

            sep="\t",

            header=None,

            dtype=str,

            comment="#"
        )

        logger.info(
            f"Loaded InterProScan data: "
            f"{len(df)} rows"
        )

        return df

    except Exception as e:

        logger.warning(
            f"Could not parse InterProScan TSV: "
            f"{e}"
        )

        return pd.DataFrame()

def get_interpro_hits(interpro_df, protein_id):

    if interpro_df.empty:
        return ""

    matches = []

    for _, row in interpro_df.iterrows():

        values = [
            str(x) for x in row.tolist()
            if pd.notna(x)
        ]

        if not values:
            continue

        first = values[0]

        if first == protein_id or first.startswith(protein_id + "."):

            if len(values) >= 12:

                matches.append(
                    f"{values[4]}: {values[5]}"
                )

                if values[11] and values[11] != "-":
                    matches.append(
                        f"InterPro={values[11]}"
                    )

    return "; ".join(dict.fromkeys(matches))

def build_annotation_summary(
    top_hits: pd.DataFrame,
    results_df: pd.DataFrame,
    interpro_df: pd.DataFrame,
    metisa_lookup: dict
) -> pd.DataFrame:

    rows = []

    candidate_ids = (
        top_hits[
            "qseqid"
        ]
        .drop_duplicates()
        .tolist()
    )

    for gene_id in candidate_ids:

        group = top_hits[
            top_hits["qseqid"]
            ==
            gene_id
        ]

        if group.empty:

            continue

        first = group.iloc[0]

        protein_id = gene_id

        if (
            "protein_id" in group.columns
            and
            pd.notna(
                first.get(
                    "protein_id"
                )
            )
        ):

            protein_id = str(
                first[
                    "protein_id"
                ]
            )

        interpro_hits = (
            get_interpro_hits(
                interpro_df,
                protein_id
            )
        )

        row = {

            "gene_id":
                gene_id,

            "protein_id":
                protein_id,

            "top_hit_sseqid":
                first.get(
                    "sseqid",
                    ""
                ),

            "top_hit_accession":
                first.get(
                    "accession",
                    ""
                ),

            "top_hit_identity_pct":
                first.get(
                    "pident",
                    ""
                ),

            "top_hit_query_coverage_pct":
                first.get(
                    "qcovs",
                    ""
                ),

            "top_hit_evalue":
                first.get(
                    "evalue",
                    ""
                ),

            "top_hit_bitscore":
                first.get(
                    "bitscore",
                    ""
                ),

            "top_hit_description":
                first.get(
                    "stitle",
                    ""
                ),

            "interpro_evidence":
                interpro_hits
        }

        if not results_df.empty:

            matching = results_df[
                results_df[
                    "gene_id"
                ].astype(str)
                ==
                str(gene_id)
            ]

        else:

            matching = pd.DataFrame()

        if not matching.empty:

            phylo = matching.iloc[0]

            row[
                "phylo_success"
            ] = phylo.get(
                "success",
                False
            )

            row[
                "phylo_sequence_count"
            ] = phylo.get(
                "sequence_count",
                ""
            )

            row[
                "alignment_length"
            ] = phylo.get(
                "alignment_length",
                ""
            )

            row[
                "trimmed_alignment_length"
            ] = phylo.get(
                "trimmed_alignment_length",
                ""
            )

            row[
                "iqtree_model"
            ] = phylo.get(
                "model",
                ""
            )

            row[
                "ufboot_mean_pct"
            ] = phylo.get(
                "ufboot_mean_pct",
                ""
            )

            row[
                "ufboot_min_pct"
            ] = phylo.get(
                "ufboot_min_pct",
                ""
            )

            row[
                "ufboot_max_pct"
            ] = phylo.get(
                "ufboot_max_pct",
                ""
            )

            row[
                "ufboot_nodes"
            ] = phylo.get(
                "ufboot_nodes",
                0
            )

            row[
                "sh_alrt_mean_pct"
            ] = phylo.get(
                "sh_alrt_mean_pct",
                ""
            )

            row[
                "sh_alrt_min_pct"
            ] = phylo.get(
                "sh_alrt_min_pct",
                ""
            )

            row[
                "sh_alrt_max_pct"
            ] = phylo.get(
                "sh_alrt_max_pct",
                ""
            )

            row[
                "sh_alrt_nodes"
            ] = phylo.get(
                "sh_alrt_nodes",
                0
            )

            row[
                "phylo_error"
            ] = phylo.get(
                "error",
                ""
            )

        rows.append(
            row
        )

    return pd.DataFrame(
        rows
    )

def main():

    print()
    print(
        "=" * 80
    )

    print(
        "METISA PLANA — SCRIPT 2 v2"
    )

    print(
        "MAFFT + trimAl + IQ-TREE"
    )

    print(
        "=" * 80
    )

    print()

    print(
        "Workflow:"
    )

    print(
        "TOP-10 RNA candidates → "
        "RAW BLAST matching candidates → "
        "Top-N homologs → blastdbcmd → "
        "MAFFT → trimAl → IQ-TREE"
    )

    print()

    print(
        "Default homologs per tree: "
        f"{DEFAULT_HOMOLOGS_PER_TREE}"
    )

    print()

    print(
        "=" * 80
    )

    print(
        "STEP 1 — INPUT FILES"
    )

    print(
        "=" * 80
    )

    print()

    rna_candidates_tsv = ask_existing_file(
        "Path to TOP RNAi CANDIDATES TSV: "
    )

    blast_tsv = ask_existing_file(
        "Path to RAW BLAST TSV: "
    )

    metisa_fasta = ask_existing_file(
        "Path to Metisa plana PROTEIN FASTA: "
    )

    interpro_tsv = ask_existing_file(
        "Path to InterProScan TSV: "
    )

    print()
    print(
        "=" * 80
    )

    print(
        "STEP 2 — PHYLOGENETIC HOMOLOG COUNT"
    )

    print(
        "=" * 80
    )

    print()

    print(
        "Your BLAST file contains the top 10 hits."
    )

    print(
        "The script will select the best N "
        "homologs for each Metisa protein."
    )

    print()

    homologs_per_tree = ask_positive_integer(
        "Number of homologs per tree",
        DEFAULT_HOMOLOGS_PER_TREE
    )

    print()
    print(
        "=" * 80
    )

    print(
        "STEP 3 — BLAST PROTEIN DATABASE"
    )

    print(
        "=" * 80
    )

    db_root, db_info = (
        ask_blast_database()
    )

    print()
    print(
        "=" * 80
    )

    print(
        "STEP 4 — SOFTWARE EXECUTABLES"
    )

    print(
        "=" * 80
    )

    print()

    mafft, _ = ask_executable(
        "MAFFT"
    )

    trimal, git_bash = ask_executable(
        "trimAl"
    )

    iqtree, _ = ask_executable(
        "IQ-TREE"
    )

    blastdbcmd, _ = ask_executable(
        "BLASTDBCMD"
    )

    print()
    print(
        "=" * 80
    )

    print(
        "STEP 5 — OUTPUT"
    )

    print(
        "=" * 80
    )

    print()

    output_folder = (
        ask_output_folder()
    )

    homolog_dir = (
        output_folder
        /
        "homolog_sequences"
    )

    phylo_fasta_dir = (
        output_folder
        /
        "phylo_fastas"
    )

    alignments_dir = (
        output_folder
        /
        "alignments"
    )

    trees_dir = (
        output_folder
        /
        "trees"
    )

    work_dir = (
        output_folder
        /
        "temp"
    )

    for directory in [

        homolog_dir,

        phylo_fasta_dir,

        alignments_dir,

        trees_dir,

        work_dir

    ]:

        directory.mkdir(
            parents=True,
            exist_ok=True
        )

    print()

    print(
        "=" * 80
    )

    print(
        "LOADING TOP RNA CANDIDATES"
    )

    print(
        "=" * 80
    )

    try:

        rna_candidates_df = pd.read_csv(
            rna_candidates_tsv,
            sep="\t",
            dtype=str
        )

    except Exception as e:

        logger.error(
            f"Cannot load RNA candidates TSV: {e}"
        )

        sys.exit(1)

    if "gene_id" not in rna_candidates_df.columns:

        logger.error(
            "RNA candidates TSV must contain "
            "a 'gene_id' column."
        )

        sys.exit(1)

    rna_candidate_ids = (
        rna_candidates_df[
            "gene_id"
        ]
        .dropna()
        .astype(str)
        .str.strip()
        .drop_duplicates()
        .tolist()
    )

    logger.info(
        f"RNA candidates loaded: "
        f"{len(rna_candidate_ids)}"
    )

    print("\n" + "=" * 80)
    print("LOADING BLAST DATA")
    print("=" * 80)

    try:
        blast_df = load_blast_tsv(blast_tsv)
    except Exception as e:
        logger.error(f"Cannot load BLAST TSV: {e}")
        sys.exit(1)

    print("\nLoading Metisa protein FASTA...")
    metisa_lookup = build_fasta_lookup(metisa_fasta)

    print("\nLoading InterProScan...")
    interpro_df = load_interpro(interpro_tsv)

    print("\n" + "=" * 80)
    print("SELECTING TOP HOMOLOGS FOR RNA CANDIDATES")
    print("=" * 80)

    blast_df["_gene_id"] = (
        blast_df["qseqid"]
        .astype(str)
        .str.strip()
        .str.replace(r"\.t\d+$", "", regex=True)
    )

    blast_candidate_df = blast_df[
        blast_df["_gene_id"].isin(rna_candidate_ids)
    ].copy()

    logger.info(
        f"BLAST rows belonging to TOP RNA candidates: "
        f"{len(blast_candidate_df)}"
    )

    blast_candidate_ids = (
        blast_candidate_df["qseqid"]
        .drop_duplicates()
        .tolist()
    )

    logger.info(
        f"RNA candidates found in BLAST: "
        f"{len(blast_candidate_ids)}"
    )

    top_hits = select_top_homologs(
        blast_candidate_df,
        max_hits=homologs_per_tree
    )

    logger.info(
        f"Selected {len(top_hits)} "
        f"homolog records"
    )

    candidate_ids = (
        top_hits[
            "qseqid"
        ]
        .drop_duplicates()
        .tolist()
    )

    logger.info(
        f"Candidates with selected homologs: "
        f"{len(candidate_ids)}"
    )

    missing_candidates = [
        gene_id
        for gene_id in rna_candidate_ids
        if gene_id not in blast_candidate_ids
    ]

    if missing_candidates:

        logger.warning(
            f"{len(missing_candidates)} RNA candidates "
            f"were not found in the BLAST TSV."
        )

        logger.warning(
            "Missing RNA candidate IDs: "
            +
            ", ".join(missing_candidates)
        )

    selected_tsv = output_folder / "selected_homologs.tsv"
    top_hits.to_csv(selected_tsv, sep="\t", index=False)

    print("\nSelected homolog table saved:")
    print(f"  {selected_tsv}")

    print("\n" + "=" * 80)
    print("RETRIEVING HOMOLOG SEQUENCES")
    print("=" * 80)

    retrieval_df = retrieve_homolog_sequences(
        blastdbcmd,
        db_root,
        top_hits,
        homolog_dir
    )

    retrieval_tsv = output_folder / "homolog_retrieval.tsv"
    retrieval_df.to_csv(retrieval_tsv, sep="\t", index=False)

    print("\n" + "=" * 80)
    print("PHYLOGENETIC ANALYSIS")
    print("=" * 80)

    results = []

    for index, gene_id in enumerate(candidate_ids, start=1):
        print("\n" + "-" * 80)
        print(f"Candidate {index}/{len(candidate_ids)}: {gene_id}")
        print("-" * 80)

        group = top_hits[top_hits["qseqid"] == gene_id]

        if group.empty:
            continue

        protein_id = gene_id

        if "protein_id" in group.columns:

            candidate_ids_from_column = (

                group[
                    "protein_id"
                ]

                .dropna()

                .astype(str)

                .tolist()
            )

            if candidate_ids_from_column:

                protein_id = (
                    candidate_ids_from_column[0]
                )

        homolog_fasta = homolog_dir / f"{safe_filename(gene_id)}_homologs.fasta"

        if not homolog_fasta.is_file():
            logger.warning(
                f"{gene_id}: homolog FASTA not found. Skipping."
            )

            results.append({
                "gene_id": gene_id,
                "protein_id": protein_id,
                "success": False,
                "sequence_count": 0,
                "error": "No homolog FASTA"
            })

            continue

        result = process_candidate(
            gene_id=gene_id,
            protein_id=protein_id,
            homolog_fasta=homolog_fasta,
            metisa_lookup=metisa_lookup,
            mafft=mafft,
            trimal=trimal,
            git_bash=git_bash,
            iqtree=iqtree,
            phylo_fasta_dir=phylo_fasta_dir,
            alignments_dir=alignments_dir,
            trees_dir=trees_dir
        )

        results.append(result)

    results_df = pd.DataFrame(results)
    results_tsv = output_folder / "phylo_results.tsv"
    results_df.to_csv(results_tsv, sep="\t", index=False)

    print()
    print("=" * 80)
    print("BUILDING ANNOTATION SUMMARY")
    print("=" * 80)

    annotation_df = build_annotation_summary(
        top_hits,
        results_df,
        interpro_df,
        metisa_lookup
    )

    annotation_tsv = output_folder / "phylo_annotation_summary.tsv"
    annotation_df.to_csv(annotation_tsv, sep="\t", index=False)

    successful = 0
    if not results_df.empty:
        successful = int((results_df["success"] == True).sum())

    failed = len(candidate_ids) - successful

    print("\n" + "=" * 80)
    print("SCRIPT 2 v2 COMPLETE")
    print("=" * 80 + "\n")

    print(f"Input BLAST rows: {len(blast_df)}")
    print(f"RNA candidates processed: {len(candidate_ids)}")
    print(f"Homologs per tree: {homologs_per_tree}")
    print(f"Selected homolog rows: {len(top_hits)}")
    print(f"Successful phylogenies: {successful}")
    print(f"Failed phylogenies: {failed}\n")

    print("BLAST database:")
    if db_info["type"] == "split":
        print("  Type: SPLIT")
        print(f"  Volumes: {len(db_info['volumes'])}")
    else:
        print("  Type: SINGLE")

    print("\nPhylogenetic settings:")
    print("  MAFFT: --auto")
    print("  trimAl: -automated1")
    print("  IQ-TREE model: MFP")
    print(f"  UFBoot: {UFBBOOT_REPLICATES}")
    print(f"  SH-aLRT: {SH_ALRT_REPLICATES}")
    print("  Threads: AUTO\n")

    print("Output files:")
    print(f"  Selected homologs:\n    {selected_tsv}\n")
    print(f"  Homolog retrieval:\n    {retrieval_tsv}\n")
    print(f"  Phylogenetic results:\n    {results_tsv}\n")
    print(f"  Annotation summary:\n    {annotation_tsv}\n")

    print("Output directories:")
    print(f"  Homolog sequences:\n    {homolog_dir}")
    print(f"  Phylogenetic FASTA:\n    {phylo_fasta_dir}")
    print(f"  Alignments:\n    {alignments_dir}")
    print(f"  Trees:\n    {trees_dir}\n")

    print("NEXT STEP:")
    print(
        "Use the phylogenetic results together "
        "with BLAST and InterPro evidence "
        "for candidate prioritization."
    )

    print("\n" + "=" * 80)

if __name__ == "__main__":
    try:
        main()

    except KeyboardInterrupt:
        print("\nProcess cancelled by user.")
        sys.exit(1)

    except Exception as e:
        print("\nFATAL ERROR:")
        print(f"  {e}")
        logger.exception("Full error:")
        sys.exit(1)