# Track 2 — Gene Mining for Pesticide Targets in *Metisa plana*

## Overview

This repository contains the computational workflow developed by **Team CSGL** for **i-Biohackathon 2026 Track 2: Gene Mining for Pesticide Targets in *Metisa plana***.

The aim was to identify and prioritise the **top 10 candidate molecular targets** from the *M. plana* predicted proteome by integrating multiple lines of biological evidence, including sequence quality, protein homology, taxonomic screening, functional annotation, developmental-stage RNA expression, phylogenetic support, and target relevance.

The workflow progressively reduces the search space from **26,490 predicted proteins → 500 candidates → 100 RNA-supported candidates → 10 final targets**.

---

## Website
https://gxxyu.github.io/-metisa-target-portal/

## Workflow

```text
*M. plana* Predicted Proteome
            │
            ▼
   Protein QC & Cleaning
            │
            ▼
       BLASTP Homology
            │
            ▼
 Taxonomic & Contamination Screening
            │
            ▼
   InterProScan / Pfam Annotation
            │
            ▼
     Candidate Mining
            │
            ▼
   Top 500 Candidates
            │
            ▼
 Developmental RNA-seq Analysis
            │
            ▼
   Top 100 RNA-supported
            │
            ▼
   MAFFT → trimAl → IQ-TREE
            │
            ▼
    Integrated Scoring
            │
            ▼
       Top 10 Targets
```

---

## Project Structure

```text
team_CSGL/
├── README.md
├── environment.yml
│
├── notebook/
│   └── analysis.ipynb
│
├── scripts/
│   ├── 01_input_qc.ipynb
│   ├── 02_blast_homology.ipynb
│   ├── 03_filter_threshold.ipynb
│   ├── 04_interproscan.ipynb
│   ├── 05_candidate_mining.ipynb
│   ├── 06_expression_ranking.ipynb
│   ├── 07_phylogeny.ipynb
│   └── 08_final_ranking.ipynb
│
├── results/
│   └── top10_candidates.xlsx
│
└── report/
    └── report.pdf
```

---

## Input Data

The primary input was the organiser-provided predicted *M. plana* proteome:

```text
proteins_longest_isoform.fasta
```

The dataset contains **26,490 predicted proteins**, with one representative protein per gene locus.

---

## Analysis Components

### 1. Protein Quality Control

* Sequence length assessment
* Ambiguous residue detection
* Stop-character detection
* Terminal `*` removal
* FASTA preparation for downstream analysis

### 2. BLASTP Homology Analysis

Protein sequences were searched against reference protein databases to identify homologous proteins and obtain sequence-similarity and taxonomic evidence.

### 3. Contamination Screening

BLAST hits were examined for taxonomic consistency, with particular attention to bacterial and fungal matches that could represent non-target or contaminant-associated sequences.

### 4. Functional Annotation

InterProScan, Pfam and related signatures were used to obtain independent functional and domain-level evidence for candidate proteins.

### 5. Candidate Mining

Homology, taxonomy, functional annotation and target relevance were integrated to produce a **500-candidate pool** for downstream prioritisation.

### 6. Developmental RNA Expression

RNA-seq datasets representing **egg, third-instar larva, pupa and adult stages** were quantified to identify candidates with relevant developmental-stage expression. The **top 100 RNA-supported candidates** were retained.

### 7. Phylogenetic Analysis

Candidate proteins and reference sequences were analysed using:

* **MAFFT** — multiple sequence alignment
* **trimAl** — alignment trimming
* **IQ-TREE** — phylogenetic inference and branch-support analysis

### 8. Final Ranking

The final candidates were ranked using integrated evidence from:

| Evidence              |  Weight |
| --------------------- | ------: |
| RNA expression        | **30%** |
| Homology              | **20%** |
| Functional annotation | **20%** |
| Phylogenetic support  | **15%** |
| Target relevance      | **15%** |

The **10 highest-ranked distinct candidates** were selected as the final target set.

---

## Key Dataset Reduction

| Stage                      |              Number |
| -------------------------- | ------------------: |
| Initial predicted proteome | **26,490 proteins** |
| Candidate-mining pool      |  **500 candidates** |
| RNA-supported pool         |  **100 candidates** |
| Final target set           |   **10 candidates** |

---

## Software

The workflow uses the following major bioinformatics tools:

* **Python**
* **NCBI BLAST+**
* **InterProScan**
* **Salmon**
* **MAFFT**
* **trimAl**
* **IQ-TREE**

Python package dependencies and versions are provided in:

```text
environment.yml
```

---

## Reproducibility

Create the conda environment:

```bash
conda env create -f environment.yml
```

Activate it:

```bash
conda activate <environment-name>
```

Then open:

```text
notebook/analysis.ipynb
```

The notebook documents the main workflow, parameters, filtering criteria, scoring procedures and analysis steps used to generate the submitted results.

Individual analysis notebooks are also provided in:

```text
scripts/
```

They can be followed sequentially from **01 → 08** when reproducing the workflow.

> **Note:** External databases and third-party software may produce different results if database versions, software versions or search parameters differ from those used in the original analysis.

---

## Main Outputs

### Final candidates

```text
results/top10_candidates.xlsx
```

Contains the final ranked Top 10 candidate targets and their integrated evidence.

### Complete analysis

```text
notebook/analysis.ipynb
```

Contains the consolidated computational workflow used for the project.

### Report

```text
report/report.pdf
```

Contains the methodology, results, discussion and interpretation of the candidate-target prioritisation.

---

## Interpretation

The Top 10 candidates represent **computationally prioritised pesticide targets**, rather than experimentally validated targets. The ranking provides a focused starting point for subsequent experimental investigation, including RNAi, biochemical and pest-control validation.

---

## Team

**Team CSGL**
i-Biohackathon 2026 — Track 2: Bio-Analyst

**Project:** Gene Mining for Pesticide Targets in *Metisa plana*
