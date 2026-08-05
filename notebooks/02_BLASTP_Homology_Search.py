# Environment setup
from __future__ import annotations

import os
import sys
import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

print(f"Python version: {sys.version}")

# Resolve notebook root in a portable way
current_dir = Path.cwd().resolve()
candidate_roots = [current_dir, current_dir.parent]
project_root = None
for root in candidate_roots:
    if (root / 'data').exists() or (root / '../data').exists():
        project_root = root
        break
if project_root is None:
    project_root = current_dir

# Prefer the relative paths used by the workflow, but allow fallback if the repo is laid out slightly differently
data_dir = project_root / 'data'
results_dir = project_root / 'results'
if not data_dir.exists() and (project_root.parent / 'data').exists():
    data_dir = project_root.parent / 'data'
if not results_dir.exists() and (project_root.parent / 'results').exists():
    results_dir = project_root.parent / 'results'

data_dir.mkdir(parents=True, exist_ok=True)
results_dir.mkdir(parents=True, exist_ok=True)

print(f"Project root: {project_root}")
print(f"Data directory: {data_dir}")
print(f"Results directory: {results_dir}")

# Check BLAST+ availability and display version
blastp_path = shutil.which('blastp')
if blastp_path is None:
    raise FileNotFoundError('BLAST+ executable `blastp` was not found in PATH. Please install BLAST+ or activate the correct environment.')

print(f"blastp executable: {blastp_path}")

blast_version = subprocess.run(
    [blastp_path, '-version'],
    capture_output=True,
    text=True,
    check=True,
)
print(blast_version.stdout.strip())

# FASTA quality assessment
from collections import defaultdict

def parse_fasta(fasta_path: Path):
    sequences = []
    current_id = None
    current_seq = []

    with fasta_path.open('r', encoding='utf-8') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    sequences.append((current_id, ''.join(current_seq)))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        sequences.append((current_id, ''.join(current_seq)))
    return sequences

query_candidates = [
    project_root / '../data/metisa_plana_proteins.fasta',
    project_root / 'data' / 'metisa_plana_proteins.fasta',
    project_root.parent / 'data' / 'metisa_plana_proteins.fasta',
]
query_file = None
for candidate in query_candidates:
    resolved = candidate.resolve()
    if resolved.exists():
        query_file = resolved
        break

if query_file is None:
    raise FileNotFoundError(
        'Metisa plana protein FASTA file was not found. Expected one of: ' + ', '.join(str(c) for c in query_candidates)
    )

seqs = parse_fasta(query_file)
lengths = [len(seq) for _, seq in seqs]
valid = []
invalid = []
for seq_id, seq in seqs:
    aa_chars = set('ACDEFGHIKLMNPQRSTVWY*X-')
    if not seq or any(char not in aa_chars for char in seq):
        invalid.append(seq_id)
    else:
        valid.append(seq_id)

print(f"Input FASTA file: {query_file}")
print(f"Total protein sequences: {len(seqs)}")
print(f"Valid sequences: {len(valid)}")
print(f"Invalid or missing sequences: {len(invalid)}")
print(f"Sequence length min: {min(lengths)} aa")
print(f"Sequence length max: {max(lengths)} aa")
print(f"Average sequence length: {sum(lengths)/len(lengths):.2f} aa")
print('Example sequence IDs:')
for seq_id, _ in seqs[:5]:
    print(f"  - {seq_id}")

if invalid:
   print()
   print(f"Warning: invalid sequence IDs found: {invalid[:10]}")

# Determine protein database for BLASTP
def find_blast_db():
    candidates = []
    env_db = os.environ.get('BLASTDB')
    if env_db:
        candidates.append(Path(env_db))
    candidates.extend([
        project_root / 'data' / 'swissprot',
        project_root / 'data' / 'uniprot_sprot',
        project_root / 'data' / 'ref_db',
        project_root / 'ref_db',
        project_root.parent / 'data' / 'swissprot',
        project_root.parent / 'data' / 'uniprot_sprot',
        Path('/usr/local/share/ncbi/blast/db/swissprot'),
        Path('/usr/local/share/ncbi/blast/db/uniprot_sprot'),
        Path('/databases/blast/swissprot'),
        Path('/databases/blast/uniprot_sprot'),
    ])

    for candidate in candidates:
        if candidate.exists():
            return candidate

    # Fallback: search for common SwissProt-like filenames in the workspace
    for root in [project_root, project_root.parent]:
        for path in root.rglob('*swissprot*'):
            if path.is_file():
                return path.with_suffix('')
        for path in root.rglob('*uniprot*'):
            if path.is_file():
                return path.with_suffix('')

    return None

blast_db = find_blast_db()
if blast_db is None:
    raise FileNotFoundError('No reference protein BLAST database was found. Please add a SwissProt or other protein database to the project data directory or set BLASTDB.')

print(f"BLAST database: {blast_db}")

blast_output = results_dir / 'blast_results.tsv'
blast_command = [
    blastp_path,
    '-query', str(query_file),
    '-db', str(blast_db),
    '-evalue', '1e-5',
    '-max_target_seqs', '5',
    '-outfmt', '6 qseqid sseqid pident length qcovs evalue bitscore stitle',
    '-out', str(blast_output),
]

print('Running BLASTP...')
result = subprocess.run(blast_command, capture_output=True, text=True)

print(f"Return code: {result.returncode}")
if result.stdout.strip():
    print(result.stdout.strip()[:500])
if result.stderr.strip():
    print('BLAST STDERR:')
    print(result.stderr.strip()[:1000])

if not blast_output.exists():
    raise FileNotFoundError(f"BLAST output file was not created: {blast_output}")

print(f"BLAST results saved to: {blast_output}")

# Load and process BLAST output
blast_cols = [
    'qseqid', 'sseqid', 'pident', 'length', 'qcovs', 'evalue', 'bitscore', 'stitle'
]

blast_df = pd.read_csv(blast_output, sep='\t', header=None, names=blast_cols, comment='#')
if blast_df.empty:
    raise ValueError(f"No BLAST hits were found in {blast_output}. Check the database and search parameters.")

blast_df['pident'] = pd.to_numeric(blast_df['pident'], errors='coerce')
blast_df['qcovs'] = pd.to_numeric(blast_df['qcovs'], errors='coerce')
blast_df['length'] = pd.to_numeric(blast_df['length'], errors='coerce')
blast_df['evalue'] = pd.to_numeric(blast_df['evalue'], errors='coerce')
blast_df['bitscore'] = pd.to_numeric(blast_df['bitscore'], errors='coerce')

strong_mask = (blast_df['pident'] >= 50) & (blast_df['qcovs'] >= 50) & (blast_df['evalue'] <= 1e-5)
moderate_mask = (blast_df['pident'] >= 30) & (blast_df['qcovs'] >= 70) & (blast_df['evalue'] <= 1e-10)
candidate_mask = strong_mask | moderate_mask

blast_df['classification'] = 'weak'
blast_df.loc[strong_mask, 'classification'] = 'strong'
blast_df.loc[moderate_mask, 'classification'] = 'moderate'

filtered_df = blast_df.loc[candidate_mask].copy()
filtered_df = filtered_df.sort_values(['classification', 'evalue', 'bitscore'], ascending=[True, True, False])

strong_hits = int(strong_mask.sum())
moderate_hits = int(moderate_mask.sum())
filtered_candidates = int(filtered_df.shape[0])

print(f"Total BLAST hits: {len(blast_df)}")
print(f"Strong hits: {strong_hits}")
print(f"Moderate hits: {moderate_hits}")
print(f"Filtered candidate hits: {filtered_candidates}")
print()
print(blast_df.head(10).to_string(index=False))

summary = pd.DataFrame({
    'category': ['total', 'strong', 'moderate', 'filtered_candidates'],
    'count': [len(blast_df), strong_hits, moderate_hits, filtered_candidates]
})
print('\nSummary table:')
print(summary.to_string(index=False))

filtered_output = results_dir / 'filtered_blast_candidates.tsv'
filtered_df.to_csv(filtered_output, sep='\t', index=False)
print(f"Saved filtered BLAST candidates: {filtered_output}")

# Plot BLAST quality metrics
plt.style.use('seaborn-v0_8-whitegrid')

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

axes[0].hist(blast_df['pident'].dropna(), bins=30, color='#2E86C1', edgecolor='black')
axes[0].set_title('Distribution of Sequence Identity (%)')
axes[0].set_xlabel('Identity (%)')
axes[0].set_ylabel('Count')

evalues = blast_df['evalue'].dropna()
axes[1].hist(evalues.map(lambda x: max(x, 1e-300)).map(lambda x: -__import__('math').log10(x)), bins=30, color='#28B463', edgecolor='black')
axes[1].set_title('Distribution of BLAST E-values')
axes[1].set_xlabel('-log10(E-value)')
axes[1].set_ylabel('Count')

hits_per_protein = blast_df.groupby('qseqid').size()
axes[2].hist(hits_per_protein.values, bins=20, color='#D35400', edgecolor='black')
axes[2].set_title('Number of Hits per Protein')
axes[2].set_xlabel('Hits per query protein')
axes[2].set_ylabel('Number of proteins')

plt.tight_layout()
plt.show()

# Prepare a candidate list for downstream annotation
candidate_df = filtered_df.copy()
candidate_df['query_length'] = candidate_df['qseqid'].map({seq_id: len(seq) for seq_id, seq in seqs})
candidate_df = candidate_df[['qseqid', 'query_length', 'sseqid', 'stitle', 'pident', 'qcovs', 'evalue', 'bitscore', 'classification']].copy()
candidate_df = candidate_df.sort_values(['classification', 'evalue', 'bitscore'], ascending=[True, True, False])
candidate_df.columns = ['query_id', 'query_length', 'subject_id', 'subject_title', 'identity_pct', 'query_coverage_pct', 'evalue', 'bitscore', 'classification']

candidate_export = results_dir / 'candidates_for_annotation.tsv'
candidate_df.to_csv(candidate_export, sep='\t', index=False)

print(f"Candidate list saved to: {candidate_export}")
print(candidate_df.head(20).to_string(index=False))