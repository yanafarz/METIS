{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "77234fa3",
   "metadata": {},
   "outputs": [],
   "source": [
    "import os\n",
    "import shutil\n",
    "import subprocess\n",
    "import sys\n",
    "from pathlib import Path"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "a8f19e52",
   "metadata": {},
   "outputs": [],
   "source": [
    "WINDOWS = sys.platform.startswith(\"win\")\n",
    "\n",
    "# Required BLAST + DB location\n",
    "DEFAULT_INSTALL_DIR = r\"C:\\BLAST\"\n",
    "\n",
    "# Common BLAST install locations (bin folder)\n",
    "COMMON_BLAST_PATHS = [\n",
    "    r\"C:\\BLAST\\bin\",  # add this FIRST (your requirement)\n",
    "    r\"C:\\Program Files\\NCBI\\blast\\bin\",\n",
    "    r\"C:\\Program Files\\NCBI\\blast-2.14.0+\\bin\",\n",
    "    r\"C:\\Program Files\\NCBI\\blast-2.13.0+\\bin\",\n",
    "    r\"C:\\Program Files\\NCBI\\blast-2.12.0+\\bin\",\n",
    "    r\"C:\\Program Files (x86)\\NCBI\\blast\\bin\",\n",
    "]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "989d90bb",
   "metadata": {},
   "outputs": [],
   "source": [
    "def run_cmd(command, capture_output=True):\n",
    "    \"\"\"Run a command and return a completed process.\"\"\"\n",
    "    return subprocess.run(\n",
    "        command,\n",
    "        shell=isinstance(command, str),\n",
    "        capture_output=capture_output,\n",
    "        text=True,\n",
    "    )"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "41ce6f94",
   "metadata": {},
   "outputs": [],
   "source": [
    "def is_windows():\n",
    "    if not WINDOWS:\n",
    "        print(\"❌ ERROR: This pipeline only supports Windows.\")\n",
    "        print(\"👉 Please run on a Windows machine.\")\n",
    "        return False\n",
    "    return True"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "fba1f8dd",
   "metadata": {},
   "outputs": [],
   "source": [
    "def find_blast():\n",
    "\n",
    "    print(\"\\n🔍 Checking BLAST installation...\")\n",
    "\n",
    "    # First: check if already in PATH\n",
    "    blastp = shutil.which(\"blastp\")\n",
    "\n",
    "    if blastp:\n",
    "        print(f\"✅ BLAST found in PATH: {blastp}\")\n",
    "        return blastp\n",
    "\n",
    "    # Second: check common install folders\n",
    "    for path in COMMON_BLAST_PATHS:\n",
    "        exe = Path(path) / \"blastp.exe\"\n",
    "        if exe.exists():\n",
    "            print(f\"✅ BLAST found at: {exe}\")\n",
    "            return str(exe)\n",
    "\n",
    "    # If not found\n",
    "    print(\"\\n❌ ERROR: BLAST not found.\")\n",
    "    print(\"👉 Please install NCBI BLAST+\")\n",
    "    print(\"👉 Expected location: C:\\\\BLAST\\\\bin\")\n",
    "    print(\"👉 Download from: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download\")\n",
    "\n",
    "    return None"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "c05db8df",
   "metadata": {},
   "outputs": [],
   "source": [
    "def check_blast_folder():\n",
    "\n",
    "    print(\"\\n📁 Checking C:\\\\BLAST folder...\")\n",
    "\n",
    "    folder = Path(DEFAULT_INSTALL_DIR)\n",
    "\n",
    "    if not folder.exists():\n",
    "        print(\"❌ ERROR: C:\\\\BLAST folder not found\")\n",
    "        print(\"👉 Please create C:\\\\BLAST and place your database there\")\n",
    "        return False\n",
    "\n",
    "    print(\"✅ C:\\\\BLAST folder exists\")\n",
    "    return True"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "id": "0ab7bae3",
   "metadata": {},
   "outputs": [],
   "source": [
    "def find_databases():\n",
    "\n",
    "    print(\"\\n🔍 Checking BLAST databases...\")\n",
    "\n",
    "    folder = Path(DEFAULT_INSTALL_DIR)\n",
    "\n",
    "    db_files = [\".pin\", \".psq\", \".phr\"]\n",
    "    databases = set()\n",
    "\n",
    "    for file in folder.iterdir():\n",
    "        if file.suffix.lower() in db_files:\n",
    "            databases.add(file.stem)\n",
    "\n",
    "    if not databases:\n",
    "        print(\"❌ ERROR: No BLAST protein database found\")\n",
    "        print(\"👉 Place database files (.pin/.psq/.phr) in C:\\\\BLAST\")\n",
    "        return None\n",
    "\n",
    "    print(\"✅ Found databases:\")\n",
    "    for db in databases:\n",
    "        print(f\"   - {db}\")\n",
    "\n",
    "    return list(databases)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "9932698c",
   "metadata": {},
   "outputs": [
    {
     "ename": "NameError",
     "evalue": "name 'is_windows' is not defined",
     "output_type": "error",
     "traceback": [
      "\u001b[31m---------------------------------------------------------------------------\u001b[39m",
      "\u001b[31mNameError\u001b[39m                                 Traceback (most recent call last)",
      "\u001b[36mCell\u001b[39m\u001b[36m \u001b[39m\u001b[32mIn[1]\u001b[39m\u001b[32m, line 1\u001b[39m\n\u001b[32m----> \u001b[39m\u001b[32m1\u001b[39m \u001b[38;5;28;01mif\u001b[39;00m \u001b[38;5;28;01mnot\u001b[39;00m is_windows():\n\u001b[32m      2\u001b[39m     sys.exit()\n\u001b[32m      3\u001b[39m \n\u001b[32m      4\u001b[39m print(\u001b[33m\"\\n==============================\"\u001b[39m)\n",
      "\u001b[31mNameError\u001b[39m: name 'is_windows' is not defined"
     ]
    }
   ],
   "source": [
    "if not is_windows():\n",
    "    sys.exit()\n",
    "\n",
    "print(\"\\n==============================\")\n",
    "print(\"BLAST INSTALLATION CHECK\")\n",
    "print(\"==============================\")\n",
    "\n",
    "\n",
    "blast = find_blast()\n",
    "\n",
    "\n",
    "if not blast:\n",
    "\n",
    "    print(\"\\n⚠️ BLAST is not installed.\")\n",
    "    print(\"Starting installation process...\")\n",
    "\n",
    "    installed = install_blast()\n",
    "\n",
    "    if not installed:\n",
    "        print(\"\\n❌ BLAST installation failed.\")\n",
    "        sys.exit()\n",
    "\n",
    "    else:\n",
    "        print(\"\\n✅ BLAST installation completed.\")\n",
    "\n",
    "\n",
    "else:\n",
    "\n",
    "    print(\"\\n✅ BLAST already available:\")\n",
    "    print(blast)\n",
    "\n",
    "\n",
    "print(\"\\n🎉 BLAST program check completed.\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "id": "191bafe6",
   "metadata": {},
   "outputs": [],
   "source": [
    "def install_blast_via_winget():\n",
    "\n",
    "    print(\"\\n🔍 Checking winget...\")\n",
    "\n",
    "    winget = shutil.which(\"winget\")\n",
    "\n",
    "    if not winget:\n",
    "        print(\"❌ winget not found\")\n",
    "        print(\"👉 Install Microsoft App Installer first\")\n",
    "        return False\n",
    "\n",
    "\n",
    "    packages = [\n",
    "        \"NCBI.NCBIblast\",\n",
    "        \"NCBI.BLAST\",\n",
    "        \"NCBI.NCBIBLAST\",\n",
    "        \"NCBI.Blast\"\n",
    "    ]\n",
    "\n",
    "\n",
    "    for package in packages:\n",
    "\n",
    "        print(f\"\\nTrying package: {package}\")\n",
    "\n",
    "        result = run_cmd(\n",
    "            f'winget install --id \"{package}\" '\n",
    "            '--accept-source-agreements '\n",
    "            '--accept-package-agreements '\n",
    "            '--silent'\n",
    "        )\n",
    "\n",
    "\n",
    "        if result.returncode == 0:\n",
    "\n",
    "            print(\"✅ BLAST installed successfully\")\n",
    "            return True\n",
    "\n",
    "\n",
    "        else:\n",
    "\n",
    "            print(\"❌ Package failed\")\n",
    "\n",
    "            if result.stderr:\n",
    "                print(result.stderr)\n",
    "\n",
    "\n",
    "    print(\"\\n❌ Automatic BLAST installation failed\")\n",
    "    return False"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "id": "cdbfb66f",
   "metadata": {},
   "outputs": [],
   "source": [
    "def install_blast_via_choco():\n",
    "\n",
    "    print(\"\\n🔍 Checking Chocolatey...\")\n",
    "\n",
    "\n",
    "    choco = shutil.which(\"choco\")\n",
    "\n",
    "\n",
    "    if not choco:\n",
    "        print(\"❌ Chocolatey not found\")\n",
    "        return False\n",
    "\n",
    "\n",
    "    result = run_cmd(\n",
    "        \"choco install ncbi-blast -y\"\n",
    "    )\n",
    "\n",
    "\n",
    "    if result.returncode == 0:\n",
    "\n",
    "        print(\"✅ BLAST installed using Chocolatey\")\n",
    "        return True\n",
    "\n",
    "\n",
    "    print(\"❌ Chocolatey installation failed\")\n",
    "\n",
    "    return False"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "id": "7a424423",
   "metadata": {},
   "outputs": [],
   "source": [
    "def install_blast():\n",
    "\n",
    "    print(\"\\n==========================\")\n",
    "    print(\"Checking BLAST installation\")\n",
    "    print(\"==========================\")\n",
    "\n",
    "\n",
    "    blast = find_blast()\n",
    "\n",
    "\n",
    "    if blast:\n",
    "\n",
    "        print(\"✅ BLAST already installed\")\n",
    "        print(blast)\n",
    "        return True\n",
    "\n",
    "\n",
    "\n",
    "    print(\"\\n❌ BLAST not found\")\n",
    "    print(\"Expected location:\")\n",
    "    print(r\"C:\\BLAST\\bin\")\n",
    "\n",
    "\n",
    "    choice = input(\n",
    "        \"\\nInstall BLAST now? (Y/N): \"\n",
    "    ).lower()\n",
    "\n",
    "\n",
    "\n",
    "    if choice not in [\"y\",\"yes\"]:\n",
    "\n",
    "        print(\"Installation cancelled\")\n",
    "        return False\n",
    "\n",
    "\n",
    "\n",
    "    if install_blast_via_winget():\n",
    "\n",
    "        return True\n",
    "\n",
    "\n",
    "\n",
    "    if install_blast_via_choco():\n",
    "\n",
    "        return True\n",
    "\n",
    "\n",
    "\n",
    "    print(\"\\n❌ Could not install BLAST\")\n",
    "    print(\"Please install manually from NCBI\")\n",
    "\n",
    "    return False"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "id": "c024753b",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "\n",
      "==========================\n",
      "Checking BLAST installation\n",
      "==========================\n",
      "\n",
      "🔍 Checking BLAST installation...\n",
      "✅ BLAST found in PATH: C:\\Users\\nurly\\OneDrive\\Documents\\UKM\\BioHackathon\\blast-2.17.0+\\bin\\blastp.EXE\n",
      "✅ BLAST already installed\n",
      "C:\\Users\\nurly\\OneDrive\\Documents\\UKM\\BioHackathon\\blast-2.17.0+\\bin\\blastp.EXE\n",
      "\n",
      "✅ BLAST installation check passed\n"
     ]
    }
   ],
   "source": [
    "blast_ready = install_blast()\n",
    "\n",
    "if blast_ready:\n",
    "    print(\"\\n✅ BLAST installation check passed\")\n",
    "else:\n",
    "    print(\"\\n❌ BLAST installation failed\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 13,
   "id": "ccc65a68",
   "metadata": {},
   "outputs": [],
   "source": [
    "def create_database_folder():\n",
    "\n",
    "    db_folder = Path(DEFAULT_INSTALL_DIR)\n",
    "\n",
    "    print(\"\\nChecking database location...\")\n",
    "    print(db_folder)\n",
    "\n",
    "    if not db_folder.exists():\n",
    "        db_folder.mkdir(parents=True)\n",
    "        print(\"✅ Created C:\\\\BLAST\")\n",
    "\n",
    "    else:\n",
    "        print(\"✅ C:\\\\BLAST already exists\")\n",
    "\n",
    "    return db_folder"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 14,
   "id": "c8b050db",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "\n",
      "Checking database location...\n",
      "C:\\BLAST\n",
      "✅ C:\\BLAST already exists\n"
     ]
    }
   ],
   "source": [
    "database_folder = create_database_folder()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "id": "68efc312",
   "metadata": {},
   "outputs": [],
   "source": [
    "def choose_database():\n",
    "\n",
    "    print(\"\\nAvailable BLAST databases\")\n",
    "    print(\"=========================\")\n",
    "    print(\"1. swissprot\")\n",
    "    print(\"2. nr\")\n",
    "    print(\"3. nt\")\n",
    "    print(\"4. refseq_protein\")\n",
    "\n",
    "    choice = input(\"\\nChoose database: \")\n",
    "\n",
    "    database_map = {\n",
    "        \"1\": \"swissprot\",\n",
    "        \"2\": \"nr\",\n",
    "        \"3\": \"nt\",\n",
    "        \"4\": \"refseq_protein\"\n",
    "    }\n",
    "\n",
    "    if choice not in database_map:\n",
    "        print(\"❌ Invalid choice\")\n",
    "        return None\n",
    "\n",
    "    return database_map[choice]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "id": "aa26803a",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "\n",
      "Available BLAST databases\n",
      "=========================\n",
      "1. swissprot\n",
      "2. nr\n",
      "3. nt\n",
      "4. refseq_protein\n",
      "❌ Invalid choice\n",
      "\n",
      "Selected database:\n",
      "None\n"
     ]
    }
   ],
   "source": [
    "db_name = choose_database()\n",
    "\n",
    "print(\"\\nSelected database:\")\n",
    "print(db_name)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "id": "c693b0fb",
   "metadata": {},
   "outputs": [],
   "source": [
    "def find_update_blastdb():\n",
    "\n",
    "    print(\"\\nSearching for update_blastdb.pl\")\n",
    "\n",
    "    for root in COMMON_BLAST_PATHS:\n",
    "\n",
    "        file = Path(root) / \"update_blastdb.pl\"\n",
    "\n",
    "        if file.exists():\n",
    "            print(\"✅ Found:\")\n",
    "            print(file)\n",
    "            return str(file)\n",
    "\n",
    "\n",
    "    print(\"❌ update_blastdb.pl not found\")\n",
    "    print(\"Database download may fail\")\n",
    "\n",
    "    return None"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 18,
   "id": "14ad8558",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "\n",
      "Searching for update_blastdb.pl\n",
      "❌ update_blastdb.pl not found\n",
      "Database download may fail\n"
     ]
    }
   ],
   "source": [
    "update_script = find_update_blastdb()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 19,
   "id": "9dedc46f",
   "metadata": {},
   "outputs": [],
   "source": [
    "def download_database(db_name):\n",
    "\n",
    "    update_script = find_update_blastdb()\n",
    "\n",
    "    if update_script is None:\n",
    "        print(\"❌ Cannot download database\")\n",
    "        return False\n",
    "\n",
    "\n",
    "    command = (\n",
    "        f'\"{update_script}\" '\n",
    "        f'--decompress '\n",
    "        f'--verbose '\n",
    "        f'{db_name} '\n",
    "        f'-outdir \"{DEFAULT_INSTALL_DIR}\"'\n",
    "    )\n",
    "\n",
    "\n",
    "    print(\"\\nRunning:\")\n",
    "    print(command)\n",
    "\n",
    "\n",
    "    result = run_cmd(command)\n",
    "\n",
    "\n",
    "    print(\"\\nSTDOUT:\")\n",
    "    print(result.stdout)\n",
    "\n",
    "\n",
    "    print(\"\\nSTDERR:\")\n",
    "    print(result.stderr)\n",
    "\n",
    "\n",
    "    if result.returncode == 0:\n",
    "\n",
    "        print(\"\\n✅ Database download complete\")\n",
    "        return True\n",
    "\n",
    "\n",
    "    else:\n",
    "\n",
    "        print(\"\\n❌ Database download failed\")\n",
    "        return False"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 20,
   "id": "9b58098e",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "\n",
      "Searching for update_blastdb.pl\n",
      "❌ update_blastdb.pl not found\n",
      "Database download may fail\n",
      "❌ Cannot download database\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "False"
      ]
     },
     "execution_count": 20,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "download_database(db_name)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 21,
   "id": "b5aed5be",
   "metadata": {},
   "outputs": [],
   "source": [
    "def check_database():\n",
    "\n",
    "    folder = Path(DEFAULT_INSTALL_DIR)\n",
    "\n",
    "    print(\"\\nChecking databases in:\")\n",
    "    print(folder)\n",
    "\n",
    "\n",
    "    db_files = [\n",
    "        \".pin\",\n",
    "        \".psq\",\n",
    "        \".phr\"\n",
    "    ]\n",
    "\n",
    "\n",
    "    databases=set()\n",
    "\n",
    "\n",
    "    for file in folder.iterdir():\n",
    "\n",
    "        if file.suffix in db_files:\n",
    "            databases.add(file.stem)\n",
    "\n",
    "\n",
    "    if databases:\n",
    "\n",
    "        print(\"\\n✅ Found databases:\")\n",
    "\n",
    "        for db in databases:\n",
    "            print(\n",
    "                folder / db\n",
    "            )\n",
    "\n",
    "    else:\n",
    "\n",
    "        print(\"\\n❌ No BLAST database detected\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 22,
   "id": "0a5a8038",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "\n",
      "Checking databases in:\n",
      "C:\\BLAST\n",
      "\n",
      "✅ Found databases:\n",
      "C:\\BLAST\\swissprot\n"
     ]
    }
   ],
   "source": [
    "check_database()"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.14.3"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
