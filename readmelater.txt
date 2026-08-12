Here is a much cleaner version you can use as README.md:

# Metisa_plana_Target_Mining

## Automated Bioinformatics Pipeline for Metisa plana Protein Target Mining

This project provides a reproducible bioinformatics pipeline for identifying,
annotating, screening, and prioritizing potential molecular targets from
predicted **Metisa plana** proteins.

The pipeline combines sequence similarity, protein-domain evidence,
functional annotation, taxonomy screening, and target-family classification
to produce a ranked list of candidate proteins for further biological
investigation.

---

# Pipeline Overview

The overall workflow is:

```text
                    Metisa plana Protein FASTA
                              |
                              v
                    [00] Quality Control
                              |
                              v
                    [01] FASTA Preparation
                              |
                              v
                    [02] BLASTP vs UniProt
                              |
                              v
                    [03] BLAST Thresholding
                              |
                              v
                    Contamination Screening
                    Bacteria / Fungi / etc.
                              |
                              v
                      Domain Annotation
                        HMMER + Pfam
                              |
                              v
                       InterProScan
                              |
                              v
                     Combined Annotation
                              |
                              v
                    Gene Family / Subfamily
                           Assignment
                              |
                              v
                      Candidate Selection
                              |
                              v
                        Candidate Scoring
                              |
                              v
                         Top-10 Table
                              |
                              v
                 Biological Justification

What the Pipeline Combines

The final annotation is based on multiple independent evidence sources:

1. Sequence similarity

BLASTP against UniProt

Provides:

Best protein homolog
Protein description
Organism/species
Sequence identity
Query coverage
E-value
Bitscore
Multiple homologous hits

2. Taxonomic evidence

Used to identify possible non-target-origin proteins, including:

Bacteria
Fungi
Plants
Other non-metazoan organisms

Taxonomic screening is used as an evidence layer and should not
automatically remove a protein based on a single weak BLAST hit.

3. Protein-domain evidence

HMMER + Pfam-A

Provides:

Pfam family
Pfam accession
Domain boundaries
Domain significance
Multiple domains within a protein

HMM-based domain detection can identify conserved protein families even
when the best BLAST hit is weak or poorly annotated.

4. Functional/domain annotation

InterProScan

Provides:

InterPro entries
Protein families
Domains
Functional signatures
GO terms
Pathways
Additional protein evidence
5. Integrated annotation

BLASTP, Pfam/HMMER, and InterProScan evidence are combined to assign:

Target class
Gene family
Functional annotation
Annotation confidence
Candidate evidence
Candidate ranking
1. System Requirements
Operating System

Recommended:

Windows 10 or Windows 11
WSL2
Ubuntu Linux distribution

The main Python pipeline can be executed from Windows PowerShell.

Linux-based tools such as:

HMMER
InterProScan
Pfam

can be executed inside WSL.

2. Python

Python 3.11 or newer is recommended.

Check the installed version:

python --version

Example:

Python 3.12.x

The pipeline does not require Python packages to install the bioinformatics
software itself.

3. Python Environment

A virtual environment is recommended.

From the project directory:

python -m venv .venv

Activate it:

.venv\Scripts\Activate.ps1

Install the required Python packages:

python -m pip install -r requirements.txt
4. Python Dependencies

The project uses:

pandas
matplotlib
openpyxl

These are listed in:

requirements.txt

Install with:

python -m pip install -r requirements.txt

Python standard-library modules such as:

os
subprocess
statistics

do not need to be included in requirements.txt.

5. WSL2

WSL2 is required for Linux-based tools used by the pipeline.

Check WSL:

wsl --status

List installed distributions:

wsl -l -v

Example:

NAME            STATE           VERSION
Ubuntu-22.04    Running         2

The pipeline should not assume that every computer uses
Ubuntu-22.04.

Users should determine their installed WSL distribution with:

wsl -l -v

Scripts should use the detected/selected distribution rather than
hard-coding a specific user's Ubuntu installation.

6. NCBI BLAST+

BLASTP is used to identify homologous proteins in UniProt.

Official NCBI BLAST+ download:

https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Download the appropriate BLAST+ package for the operating system.

After installation, verify:

blastp -version

Example:

blastp: 2.x.x+

The BLAST executable must be available to the pipeline either through
the system PATH or through the configured BLAST installation directory.

7. UniProt Reference Database

BLASTP uses UniProt protein sequences as the primary reference database.

UniProt downloads:

https://www.uniprot.org/help/downloads

The reference dataset should include:

Swiss-Prot reviewed proteins
TrEMBL unreviewed proteins
Isoform sequences where required by the database construction

The purpose is to provide broad protein coverage while retaining the
higher-quality reviewed Swiss-Prot annotations.

The Swiss-Prot and TrEMBL datasets should be combined into the reference
protein FASTA used to create the local BLAST database.

Recommended UniProt Reference

The conceptual reference set is:

UniProt
├── Swiss-Prot
│   └── Reviewed proteins
│
└── TrEMBL
    └── Unreviewed proteins

Isoform sequences should also be included when constructing the reference
dataset if isoform-level matches are desired.

8. Building the Local UniProt BLAST Database

After obtaining the UniProt protein FASTA:

makeblastdb `
    -in uniprot_combined.fasta `
    -dbtype prot `
    -parse_seqids `
    -out uniprot_all

The resulting database will contain files such as:

uniprot_all.pin
uniprot_all.psq
uniprot_all.phr

The database can be stored, for example, as:

C:\BLAST\db\uniprot_all

The pipeline scripts should ask the user for the database location rather
than assuming a specific computer path.

9. BLASTP Settings

The primary homology search uses:

BLASTP

Recommended settings for the candidate-mining stage:

Maximum target sequences: 10

This is important because the pipeline needs multiple homologous proteins
rather than relying only on a single best hit.

The BLAST output should retain information such as:

Query protein ID
Subject accession
Percentage identity
Alignment length
Query length
Subject length
Query coverage
E-value
Bitscore
Subject description

Example output fields:

qseqid
sseqid
pident
length
qlen
slen
qcovs
evalue
bitscore
stitle
10. HMMER

HMMER is used for protein-domain detection.

Check installation inside WSL:

hmmscan -h

or:

hmmsearch -h

The pipeline uses HMMER to search protein sequences against Pfam HMM
profiles.

11. Pfam-A

Pfam-A is used as the HMM library for protein-domain annotation.

Example location:

~/pfam/Pfam-A.hmm

After downloading:

hmmpress ~/pfam/Pfam-A.hmm

This creates:

Pfam-A.hmm
Pfam-A.hmm.h3f
Pfam-A.hmm.h3i
Pfam-A.hmm.h3m
Pfam-A.hmm.h3p

The complete Pfam-A library should be retained.

Do not restrict the HMM library to only chitin, detoxification, or
insecticide-related families before the general annotation stage.

Target-specific filtering is performed later.

This allows unexpected but biologically relevant protein families to be
detected.

12. InterProScan

InterProScan provides additional protein-domain, family, and functional
annotation.

It can provide:

InterPro entries
Protein families
Domains
Functional signatures
GO terms
Pathways

InterProScan is installed separately from Python.

The Python script asks the user to provide the InterProScan directory.

Example:

/home/<username>/interproscan/

The executable should be:

interproscan.sh

Check inside WSL:

cd ~/interproscan
./interproscan.sh --version
13. Java

InterProScan requires Java.

The current pipeline expects Java 11.

Check:

java -version

If Java 11 is not installed:

sudo apt update
sudo apt install openjdk-11-jdk

Then:

java -version

The pipeline should verify the Java installation rather than assuming
that every computer has Java installed in exactly the same location.

14. Input Data

The primary input is a predicted protein FASTA file from the
Metisa plana genome annotation.

Example:

data/metisa_plana_proteins.fasta

Example FASTA:

>g1.t1
MXXXXXX...
>g1.t2
MXXXXXX...

Each protein should have a unique sequence identifier.

15. Pipeline Steps
[00] Input Quality Control

Script:

01_input_qc.py

Purpose:

Detect DNA/protein FASTA
Count sequences
Check duplicate IDs
Check empty sequences
Calculate sequence-length statistics
Detect invalid characters
Calculate protein statistics
Generate QC reports
Generate length distribution plot

Outputs include:

QC summary
Sequence statistics
Excel report
TSV report
Length distribution plot
[01] FASTA Preparation

Purpose:

Standardize FASTA sequences
Remove terminal protein stop symbols
Detect internal stop symbols
Detect invalid characters
Check duplicate IDs
Generate BLAST-ready FASTA

This step does not:

Predict genes
Translate DNA
Run BLAST
[02] BLASTP

Purpose:

Identify homologous proteins using the local UniProt database.

Reference:

UniProt Swiss-Prot + TrEMBL

Recommended:

Maximum target sequences = 10

The output is a TSV containing multiple homologous hits per query protein.

[03] BLAST Thresholding

BLAST hits are filtered according to sequence similarity evidence.

The pipeline should consider:

E-value
Percentage identity
Query coverage
Bitscore

A strong hit should have both:

high identity
+
high query coverage
+
low E-value

A low E-value alone does not necessarily indicate a strong functional
annotation.

Contamination Screening

Taxonomic screening is performed against BLAST evidence.

Potential non-target organisms include:

Bacteria
Fungi
Plants
Other non-metazoans

The purpose is to identify proteins whose strongest evidence is inconsistent
with the expected Metisa plana biological origin.

Important:

A single weak bacterial or fungal hit should not automatically result in
protein removal.

Taxonomic evidence should be combined with:

identity
coverage
E-value
number of supporting hits
taxonomic consistency
Pfam domains
InterPro annotation
Domain Annotation
HMMER + Pfam

Protein sequences are searched against Pfam-A HMM profiles.

This provides domain-level evidence that is independent of the BLAST
best-hit description.

Example:

Protein
  |
  +-- PFxxxxx
  |
  +-- PFyyyyy
  |
  +-- PFzzzzz

Multiple domains should be retained.

InterProScan Annotation

InterProScan provides additional:

Domain evidence
Family evidence
Signature evidence
GO terms
Pathways
InterPro classifications

InterProScan is combined with BLASTP and Pfam/HMMER rather than replacing
them.

Combined Annotation

The evidence from:

BLASTP
   +
HMMER / Pfam
   +
InterProScan
   +
Taxonomy

is combined into a unified annotation table.

The combined table can contain:

Protein_ID
Protein_Length
BLAST_Best_Hit
BLAST_Identity
BLAST_Query_Coverage
BLAST_Evalue
BLAST_Bitscore
Pfam_Domains
InterPro_Entries
GO_Terms
Pathways
Taxonomic_Evidence
Functional_Annotation
Gene_Family
Target_Class
Annotation_Confidence
Gene Family / Subfamily Assignment

Gene-family assignment should not rely exclusively on the BLAST best hit.

The assignment should consider:

BLAST homologs
BLAST identity
Query coverage
E-value
Pfam domains
InterPro families
Conserved signatures
Taxonomic consistency
Functional descriptions

The objective is to distinguish:

Protein
   |
   v
Gene family
   |
   v
Subfamily / clade
   |
   v
Functional role
   |
   v
Target class
Candidate Target Classes

The pipeline can prioritize classes such as:

Detoxification
Cytochrome P450
Glutathione S-transferase
Carboxylesterase
UDP-glycosyltransferase
ABC transporter
Other detoxification-associated proteins
Chitin / Cuticle
Chitin synthase
Chitinase
Chitin deacetylase
Chitin-binding proteins
Cuticle-associated proteins
Insecticide / Physiological Targets

Examples include:

Acetylcholinesterase
Nicotinic acetylcholine receptor
GABA receptor
Voltage-gated sodium channel
Ryanodine receptor
Acetyl-CoA carboxylase
Ecdysone receptor
Other validated insect molecular targets

The target-class list can be expanded without changing the upstream
annotation pipeline.

Candidate Selection

Candidates are selected after annotation and contamination screening.

Candidate selection should consider:

Target-class membership
Gene-family assignment
Domain evidence
BLAST evidence
InterPro evidence
Taxonomic consistency
Functional annotation confidence
Candidate Scoring

Candidates are ranked using a multi-evidence scoring framework.

Example categories:

Evidence category                 Weight
------------------------------------------------
Essential / critical role          30%
Known insect target                25%
Strong annotation evidence         20%
Larval expression support          15%
Species specificity                10%

The weights can be modified according to the biological objective.

Top-10 Candidate Table

The final output should contain the highest-ranked candidates.

Example:

Rank	Protein_ID	Gene Family	Target Class	Score	Confidence
1	...	...	...	...	...
2	...	...	...	...	...
3	...	...	...	...	...

The final table should also retain the evidence supporting each candidate.

Biological Justification

Each final candidate should receive a short biological justification.

The justification should explain:

What protein family it belongs to
Which domains support the assignment
What BLAST evidence supports the annotation
What InterPro/Pfam evidence supports the annotation
Whether the protein is relevant to an insect biological process
Whether contamination evidence is present
Why the protein is suitable for further target investigation

The goal is to produce an evidence-based candidate list rather than
simply ranking proteins by BLAST score.

Recommended Final Output

The final candidate table should contain at least:

Rank
Protein_ID
Protein_Length_AA
Target_Class
Gene_Family
Subfamily_or_Clade
Best_Hit_Accession
Best_Hit_Species
Pfam_Domains
InterPro_Domains
Functional_Annotation
Evalue
Identity_Pct
Query_Coverage_Pct
Taxonomic_Evidence
Expression_Support
Life_Stage_Relevance
Offtarget_Check
Annotation_Confidence
Candidate_Score
Justification
Reproducibility

The pipeline is designed so that another user can reproduce the analysis
without using the original computer's paths.

The scripts should request or configure:

Input FASTA
BLAST database
Pfam-A HMM
InterProScan directory
Output directory
WSL distribution where required

User-specific paths such as:

C:\Users\nurly\...
/home/nurly/...

should not be hard-coded into the pipeline.

External Software Summary
Software / Database	Purpose	Environment
Python	Pipeline scripts	Windows
pandas	Data processing	Python
matplotlib	QC plots	Python
openpyxl	Excel output	Python
NCBI BLAST+	BLASTP	Windows
UniProt	Protein reference database	Local
WSL2	Linux environment	Windows
Ubuntu	Linux tools	WSL
HMMER	HMM searches	WSL
Pfam-A	Protein HMM library	WSL
InterProScan	Domain/family annotation	WSL
Java 11	InterProScan runtime	WSL
Important Notes
BLAST

The BLAST reference database should be updated periodically because
UniProt is continuously updated.

Record the UniProt release/date used for each analysis.

Pfam

Record the Pfam release/date used for each analysis.

InterProScan

Record the InterProScan version and database version.

BLAST settings

The pipeline uses:

Maximum target sequences = 10

unless explicitly changed for a particular analysis.

Reproducibility

For every final analysis, record:

Input FASTA
BLAST version
UniProt release
Pfam release
HMMER version
InterProScan version
Java version
Threshold parameters
Candidate scoring weights

This allows the analysis to be reproduced later.

Project Structure

Recommended:

Metisa_plana_Target_Mining/
│
├── README.md
├── SETUP.md
├── requirements.txt
│
├── scripts/
│   ├── 01_input_qc.py
│   ├── 02_fasta_preparation.py
│   ├── 03_blastp.py
│   ├── 04_blast_threshold.py
│   ├── 05_hmmer.py
│   ├── 06_interproscan.py
│   ├── 07_combine_annotations.py
│   ├── 08_target_classification.py
│   ├── 09_candidate_scoring.py
│   └── ...
│
├── data/
│
├── databases/
│
└── output/

Large reference databases should normally be stored outside the Git
repository or managed separately rather than committed to Git.

Citation / Software Sources

NCBI BLAST+:

https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

UniProt downloads:

https://www.uniprot.org/help/downloads

Pfam:

https://www.ebi.ac.uk/interpro/

InterPro:

https://www.ebi.ac.uk/interpro/

HMMER:

https://www.ebi.ac.uk/Tools/hmmer/
