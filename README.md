# Metisa_plana_target_mining

for blast can install from here: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

for uniprot database here: [https://www.uniprot.org/help/downloads](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/)
 (make sure down all; both swiss and trembl also isoform)

when do blast make 10 max target seq

             pipeline:
     
     quality check [00] -> clean file [01] -> blastp [02] -> threshold [03]
                           |
     contamination screening -> (bacteria etc) -> domain annotation ->  gene family/subfamily assignment
                           |
     candidate target selection -> Scoring candidates -> top-10 table -> biological justification


The pipeline combines:

- FASTA quality control
- FASTA preparation
- BLASTP against UniProt
- Taxonomic contamination screening
- HMMER/Pfam domain annotation
- InterProScan annotation
- Combined annotation
- Target-family classification
- Candidate extraction and scoring

---

# 1. Requirements

## Operating system

Recommended:

- Windows 10/11
- WSL2
- Ubuntu Linux distribution

The pipeline can run Python scripts from Windows PowerShell,
while tools such as HMMER and InterProScan can run inside WSL.

---

# 2. Python

Install Python 3.11 or newer.

Check:

```powershell
python --version
