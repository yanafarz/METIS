import pandas as pd
import numpy as np
import re
from pathlib import Path

TOP_N = 10

WEIGHT_RNA = 30
WEIGHT_HOMOLOGY = 20
WEIGHT_ANNOTATION = 20
WEIGHT_PHYLO = 15
WEIGHT_TARGET = 15

def ask_file(label):

    while True:

        value = input(
            f"\nEnter {label} path:\n> "
        ).strip().strip('"')

        path = Path(value)

        if path.is_file():
            return path

        print(
            f"ERROR: File not found:\n{path}"
        )


def ask_output_directory():

    value = input(
        "\nEnter OUTPUT DIRECTORY "
        "(press Enter to use the same folder as the RNA file):\n> "
    ).strip().strip('"')

    if not value:
        return None

    path = Path(value)

    path.mkdir(
        parents=True,
        exist_ok=True
    )

    return path

def load_tsv(path, label):

    print(
        f"\nLoading {label}..."
    )

    try:

        df = pd.read_csv(
            path,
            sep="\t",
            dtype=str,
            low_memory=False
        )

    except Exception as e:

        raise RuntimeError(
            f"Could not read {label}:\n{e}"
        )

    df.columns = [
        str(c).strip()
        for c in df.columns
    ]

    print(
        f"{label} rows loaded: "
        f"{len(df):,}"
    )

    return df

def require_columns(
    df,
    columns,
    label
):

    missing = [
        c
        for c in columns
        if c not in df.columns
    ]

    if missing:

        raise ValueError(
            f"\n{label} is missing required columns:\n"
            + "\n".join(
                f"  - {x}"
                for x in missing
            )
        )

def numeric(df, column):

    if column not in df.columns:

        return pd.Series(
            np.nan,
            index=df.index
        )

    return pd.to_numeric(
        df[column],
        errors="coerce"
    )

def normalize_gene_id(value):

    if pd.isna(value):
        return ""

    value = str(value).strip()

    match = re.match(
        r"^(g\d+)",
        value,
        flags=re.IGNORECASE
    )

    if match:
        return match.group(1)

    return value

def clean_text(value):

    if pd.isna(value):
        return ""

    value = str(value).strip()

    if value.lower() in {
        "nan",
        "none",
        "null",
        "-"
    }:
        return ""

    return value

def extract_accession(value):

    value = clean_text(value)

    if not value:
        return ""

    parts = value.split("|")

    if len(parts) >= 2:

        accession = parts[1].strip()

        if accession:
            return accession

    return value

def extract_species(description):

    description = clean_text(
        description
    )

    if not description:
        return ""

    match = re.search(
        r"OS=(.*?)(?:\s+OX=|\s+GN=|\s+PE=|\s+SV=|$)",
        description
    )

    if match:
        return match.group(1).strip()

    return ""

def extract_function(description):

    description = clean_text(
        description
    )

    if not description:
        return ""

    description = re.sub(
        r"^(?:sp|tr)\|[^|]+\|[^ ]+\s+",
        "",
        description,
        flags=re.IGNORECASE
    )

    description = re.split(
        r"\s+OS=",
        description,
        maxsplit=1
    )[0]

    return description.strip()

def extract_pfam(interpro_rows):

    if interpro_rows.empty:
        return ""

    values = []

    if "Signature accession" in interpro_rows.columns:

        for value in interpro_rows[
            "Signature accession"
        ]:

            value = clean_text(value)

            if re.match(
                r"^PF\d+$",
                value
            ):

                values.append(value)

    return ";".join(
        dict.fromkeys(values)
    )

def extract_interpro(interpro_rows):

    if interpro_rows.empty:
        return ""

    values = []

    if "InterPro accession" in interpro_rows.columns:

        for value in interpro_rows[
            "InterPro accession"
        ]:

            value = clean_text(value)

            if re.match(
                r"^IPR\d+$",
                value
            ):

                values.append(value)

    return ";".join(
        dict.fromkeys(values)
    )

def extract_interpro_descriptions(
    interpro_rows
):

    if interpro_rows.empty:
        return ""

    values = []

    if "InterPro description" in interpro_rows.columns:

        for value in interpro_rows[
            "InterPro description"
        ]:

            value = clean_text(value)

            if value:
                values.append(value)

    return "; ".join(
        dict.fromkeys(values)
    )

def determine_gene_family(
    functional_annotation,
    interpro_description,
    pfam
):

    text = (
        clean_text(functional_annotation)
        + " "
        + clean_text(interpro_description)
        + " "
        + clean_text(pfam)
    )

    lower = text.lower()

    rules = [

        ("chitin synthase", "Chitin synthase"),

        ("trypsin", "Trypsin / Peptidase S1"),

        ("chymotrypsin", "Chymotrypsin / Peptidase S1"),

        ("serine protease", "Serine protease"),

        ("collagenase", "Collagenase / Serine protease"),

        ("pancreatic triacylglycerol lipase",
         "Triacylglycerol lipase"),

        ("triacylglycerol lipase",
         "Triacylglycerol lipase"),

        ("lipase", "Lipase"),

        ("nucleoside diphosphate kinase",
         "Nucleoside diphosphate kinase"),

        ("mitochondrial-processing peptidase",
         "Mitochondrial processing peptidase"),

        ("mitochondrial processing peptidase",
         "Mitochondrial processing peptidase"),

        ("presequence protease",
         "Mitochondrial presequence protease"),

        ("serine protease inhibitor",
         "Serine protease inhibitor"),

        ("protease inhibitor",
         "Protease inhibitor"),

        ("sh3 domain",
         "SH3-domain protein"),

        ("immunoglobulin-like",
         "Immunoglobulin-like protein"),

        ("leucine-rich repeat",
         "Leucine-rich repeat protein"),

        ("ef-hand",
         "EF-hand protein"),

        ("transcription factor",
         "Transcription factor"),
    ]

    for keyword, family in rules:

        if keyword in lower:
            return family

    ipr = clean_text(
        interpro_description
    )

    if ipr:

        first = ipr.split(";")[0].strip()

        if first:
            return first[:120]

    annotation = clean_text(
        functional_annotation
    )

    if annotation:
        return annotation[:120]

    return ""

def determine_target_class(
    functional_annotation,
    gene_family,
    interpro_description,
    existing_target_class
):

    annotation = (
        clean_text(functional_annotation)
        + " "
        + clean_text(gene_family)
        + " "
        + clean_text(interpro_description)
    )

    lower = annotation.lower()

    # CHITIN / STRUCTURAL

    if (
        "chitin synthase" in lower
        or "chitin pathway" in lower
    ):

        return "STRUCTURAL/CHITIN"

    # DIGESTION

    if any(
        x in lower
        for x in [
            "trypsin",
            "chymotrypsin",
            "serine protease",
            "collagenase",
            "lipase",
            "peptidase"
        ]
    ):

        if (
            "mitochondrial" not in lower
            and "presequence" not in lower
        ):

            return "DIGESTION"

    # TARGET SITE

    if any(
        x in lower
        for x in [
            "ion channel",
            "ligand-gated",
            "voltage-gated",
            "acetylcholine receptor",
            "nicotinic acetylcholine receptor",
            "gaba receptor",
            "gaba-gated",
            "glutamate-gated",
            "ryanodine receptor"
        ]
    ):

        return "TARGET-SITE"

    # SIGNALLING

    if any(
        x in lower
        for x in [
            "signaling",
            "signalling",
            "signal transduction",
            "sh3 domain"
        ]
    ):

        return "SIGNALING"

    existing = clean_text(
        existing_target_class
    )

    if existing:
        return existing.upper()

    return ""

def calculate_phylo_score(row):

    if not row["phylo_success"]:
        return 0.0

    values = []

    uf = row.get(
        "ufboot_mean_pct",
        np.nan
    )

    sh = row.get(
        "sh_alrt_mean_pct",
        np.nan
    )

    if pd.notna(uf):
        values.append(float(uf))

    if pd.notna(sh):
        values.append(float(sh))

    if not values:
        return 0.0

    return min(
        100.0,
        np.mean(values)
    )

def calculate_homology_score(
    identity,
    coverage,
    evalue
):

    if pd.isna(identity):
        identity = 0

    if pd.isna(coverage):
        coverage = 0

    if pd.isna(evalue):
        evalue = 1.0

    identity_score = min(
        100,
        max(
            0,
            (identity - 20)
            / 60
            * 100
        )
    )

    coverage_score = min(
        100,
        max(
            0,
            coverage
        )
    )

    if evalue <= 1e-50:
        evalue_score = 100

    elif evalue <= 1e-20:
        evalue_score = 90

    elif evalue <= 1e-10:
        evalue_score = 80

    elif evalue <= 1e-5:
        evalue_score = 70

    elif evalue <= 0.05:
        evalue_score = 50

    else:
        evalue_score = 0

    return (
        identity_score * 0.40
        + coverage_score * 0.40
        + evalue_score * 0.20
    )

def calculate_rna_score(row):

    score = row.get(
        "RNA_score_numeric",
        np.nan
    )

    if pd.notna(score):

        return min(
            100,
            max(
                0,
                float(score)
            )
        )

    rank = row.get(
        "RNA_rank_numeric",
        np.nan
    )

    if pd.isna(rank):
        return 0

    if rank <= 10:
        return 100

    if rank <= 25:
        return 90

    if rank <= 50:
        return 80

    if rank <= 75:
        return 70

    return 60

def calculate_stage_expression_score(
    tpm,
    all_stage_tpm
):

    if pd.isna(tpm):
        return 0.0

    tpm = max(
        0.0,
        float(tpm)
    )

    if tpm <= 0:
        return 0.0

    other_values = [
        float(x)
        for x in all_stage_tpm
        if pd.notna(x)
        and float(x) >= 0
    ]

    if not other_values:
        return 0.0

    maximum = max(
        other_values
    )

    if maximum <= 0:
        return 0.0

    # Relative expression score.
    # 100 = highest expression among stages.

    score = (
        tpm / maximum
    ) * 100

    return round(
        min(100, score),
        2
    )

def calculate_annotation_score(
    interpro_count,
    pfam_count,
    functional_annotation
):

    score = 0

    if interpro_count > 0:
        score += 50

    if pfam_count > 0:
        score += 30

    if clean_text(
        functional_annotation
    ):

        score += 20

    return min(
        100,
        score
    )

def calculate_target_score(
    target_class,
    functional_annotation
):

    target = clean_text(
        target_class
    ).upper()

    annotation = clean_text(
        functional_annotation
    ).lower()

    if target == "STRUCTURAL/CHITIN":
        return 100

    if target == "TARGET-SITE":
        return 95

    if target == "DIGESTION":

        if any(
            x in annotation
            for x in [
                "trypsin",
                "chymotrypsin",
                "lipase",
                "collagenase",
                "serine protease"
            ]
        ):

            return 90

        return 70

    if target == "SIGNALING":
        return 75

    if target:
        return 50

    return 0

def create_functional_annotation(
    functional_annotation,
    gene_family,
    target_class
):

    source = clean_text(
        functional_annotation
    )

    family = clean_text(
        gene_family
    )

    target = clean_text(
        target_class
    ).upper()

    if not source:
        source = family

    if not source:
        return "Uncharacterized protein with limited functional annotation evidence"

    # Remove database-style extras.

    source = re.sub(
        r"\([^)]*\)",
        "",
        source
    )

    source = re.sub(
        r"\s+",
        " ",
        source
    ).strip()

    # Remove EC numbers and similar trailing information.

    source = re.sub(
        r"\s+EC=\S+.*$",
        "",
        source,
        flags=re.IGNORECASE
    )

    words = source.split()

    # If already within requested range.

    if 7 <= len(words) <= 15:
        return source

    # Build a concise biological description.

    if "chitin synthase" in source.lower():

        text = (
            "Chitin synthase involved in insect cuticle formation "
            "and developmental moulting"
        )

    elif "trypsin" in source.lower():

        text = (
            "Trypsin-like serine protease involved in larval "
            "dietary protein digestion"
        )

    elif "chymotrypsin" in source.lower():

        text = (
            "Chymotrypsin-like serine protease involved in "
            "larval digestive protein hydrolysis"
        )

    elif "lipase" in source.lower():

        text = (
            "Lipase involved in dietary lipid hydrolysis "
            "during larval feeding"
        )

    elif "collagenase" in source.lower():

        text = (
            "Collagenase-like serine protease associated with "
            "protein digestion during larval feeding"
        )

    elif "protease inhibitor" in source.lower():

        text = (
            "Protease inhibitor regulating proteolytic activity "
            "within insect physiological processes"
        )

    elif target == "TARGET-SITE":

        text = (
            f"{family or 'Receptor-like protein'} associated with "
            "insect neuronal signalling and physiological regulation"
        )

    elif target == "SIGNALING":

        text = (
            f"{family or 'Signalling protein'} involved in "
            "cellular signal transduction and developmental regulation"
        )

    else:

        # Use first 15 words if possible.

        text = " ".join(
            words[:15]
        )

    # Enforce 7–15 words.

    words = text.split()

    if len(words) < 7:

        additions = [
            "in insect biological processes",
            "during normal insect development",
            "with predicted biological activity"
        ]

        for addition in additions:

            words.extend(
                addition.split()
            )

            if len(words) >= 7:
                break

    return " ".join(
        words[:15]
    )

def determine_phylogenetic_group(row):

    group = clean_text(
        row.get(
            "phylogenetic_group",
            ""
        )
    )

    uf = row.get(
        "ufboot_mean_pct",
        np.nan
    )

    sh = row.get(
        "sh_alrt_mean_pct",
        np.nan
    )

    if group:

        if pd.notna(uf):

            return (
                f"{group} "
                f"(UFBoot {float(uf):.0f}%)"
            )

        return group

    if not row.get(
        "phylo_success",
        False
    ):

        return "Not resolved"

    if pd.notna(uf):

        if pd.notna(sh):

            return (
                "Supported phylogenetic placement "
                f"(UFBoot {float(uf):.0f}%; "
                f"SH-aLRT {float(sh):.0f}%)"
            )

        return (
            "Supported phylogenetic placement "
            f"(UFBoot {float(uf):.0f}%)"
        )

    return "Phylogenetic placement obtained"

def phylo_evidence(row):

    if not row["phylo_success"]:
        return ""

    model = clean_text(
        row.get(
            "iqtree_model",
            ""
        )
    )

    uf = row.get(
        "ufboot_mean_pct",
        np.nan
    )

    sh = row.get(
        "sh_alrt_mean_pct",
        np.nan
    )

    text = (
        "MAFFT + trimAl + IQ-TREE"
    )

    if model:
        text += (
            f"; model={model}"
        )

    if pd.notna(uf):

        text += (
            f"; UFBoot={float(uf):.2f}%"
        )

    if pd.notna(sh):

        text += (
            f"; SH-aLRT={float(sh):.2f}%"
        )

    return text

def generate_justification(
    gene_family,
    functional_annotation,
    target_class,
    larva_tpm,
    adult_tpm,
    pupa_tpm,
    egg_tpm,
    phylogenetic_group
):

    family = clean_text(
        gene_family
    )

    annotation = clean_text(
        functional_annotation
    )

    target = clean_text(
        target_class
    ).upper()

    stages = {
        "larval": larva_tpm,
        "adult": adult_tpm,
        "pupal": pupa_tpm,
        "egg": egg_tpm
    }

    valid = {
        k: float(v)
        for k, v in stages.items()
        if pd.notna(v)
    }

    larva = (
        float(larva_tpm)
        if pd.notna(larva_tpm)
        else 0.0
    )

    adult = (
        float(adult_tpm)
        if pd.notna(adult_tpm)
        else 0.0
    )

    pupa = (
        float(pupa_tpm)
        if pd.notna(pupa_tpm)
        else 0.0
    )

    egg = (
        float(egg_tpm)
        if pd.notna(egg_tpm)
        else 0.0
    )

    non_larval = [
        adult,
        pupa,
        egg
    ]

    max_non_larval = max(
        non_larval
    )

    if larva > max_non_larval * 2:

        expression_sentence = (
            f"RNA expression is strongly larva-biased "
            f"(Larva TPM {larva:.1f} versus Adult {adult:.1f}, "
            f"Pupa {pupa:.1f}, and Egg {egg:.1f})."
        )

    elif larva > max_non_larval:

        expression_sentence = (
            f"RNA expression is highest in larvae "
            f"(Larva TPM {larva:.1f}; Adult {adult:.1f}, "
            f"Pupa {pupa:.1f}, Egg {egg:.1f})."
        )

    else:

        expression_sentence = (
            f"RNA expression is detectable across stages "
            f"(Larva TPM {larva:.1f}, Adult {adult:.1f}, "
            f"Pupa {pupa:.1f}, Egg {egg:.1f})."
        )

    if (
        "chitin synthase" in family.lower()
        or "chitin synthase" in annotation.lower()
    ):

        return (
            "Cuticular chitin synthase supports chitin production "
            "required for insect cuticle formation and moulting. "
            f"{expression_sentence} "
            "Reducing chitin synthesis could compromise cuticle "
            "integrity and interfere with successful moulting. "
            "The combination of functional annotation, developmental "
            "relevance and phylogenetic support makes this a strong "
            "RNAi candidate."
        )

    if (
        "trypsin" in family.lower()
        or "trypsin" in annotation.lower()
    ):

        return (
            "Trypsin-like serine proteases contribute to dietary "
            "protein digestion in feeding insect larvae. "
            f"{expression_sentence} "
            "Knockdown could reduce proteolytic capacity and limit "
            "amino-acid acquisition during feeding. This provides "
            "a plausible mechanism for impaired larval growth "
            "and development."
        )

    if (
        "chymotrypsin" in family.lower()
        or "chymotrypsin" in annotation.lower()
    ):

        return (
            "Chymotrypsin-like serine proteases contribute to "
            "hydrolysis of dietary proteins in the larval gut. "
            f"{expression_sentence} "
            "RNAi-mediated reduction could decrease digestive "
            "capacity and nutrient acquisition. This makes the "
            "candidate biologically relevant for targeting the "
            "feeding stage of M. plana."
        )

    if (
        "lipase" in family.lower()
        or "lipase" in annotation.lower()
    ):

        return (
            "Lipases support hydrolysis of dietary lipids during "
            "larval feeding and nutrient acquisition. "
            f"{expression_sentence} "
            "Knockdown could reduce the availability of fatty acids "
            "needed for energy production and cellular processes. "
            "This provides a plausible digestive vulnerability "
            "for RNAi-based control."
        )

    if "collagenase" in annotation.lower():

        return (
            "The candidate is annotated as a collagenase-like "
            "serine protease with predicted proteolytic activity. "
            f"{expression_sentence} "
            "Knockdown could reduce digestive proteolysis and "
            "nutrient acquisition in feeding larvae. The predicted "
            "mechanism is plausible, although the specific substrate "
            "and phenotype require experimental validation."
        )

    if target == "TARGET-SITE":

        return (
            f"The candidate is annotated as {annotation or family}, "
            "consistent with a receptor or ion-channel function. "
            f"{expression_sentence} "
            "RNAi disruption could alter neuronal signalling or "
            "other physiological processes controlled by the target. "
            "The exact phenotype depends on the protein's biological "
            "role and requires experimental validation."
        )

    if target == "SIGNALING":

        return (
            f"The candidate is associated with {annotation or family} "
            "and may contribute to cellular signal transduction. "
            f"{expression_sentence} "
            "Knockdown could disrupt signalling processes required "
            "for normal development or physiological function. "
            "The specific downstream phenotype cannot be established "
            "from sequence evidence alone."
        )

    if "protease inhibitor" in family.lower():

        return (
            "The candidate encodes a protease inhibitor that may "
            "regulate proteolytic activity in the insect. "
            f"{expression_sentence} "
            "Knockdown could disturb the balance between proteases "
            "and their inhibitors, potentially affecting digestion "
            "or other physiological processes. Functional validation "
            "is needed to establish the relevant phenotype."
        )

    if target == "DIGESTION":

        return (
            f"The candidate is associated with {annotation or family} "
            "and is predicted to contribute to digestive function. "
            f"{expression_sentence} "
            "Knockdown could reduce the corresponding digestive "
            "activity and limit nutrient acquisition during feeding. "
            "This provides a plausible RNAi mechanism, although "
            "functional redundancy should be considered."
        )

    return (
        f"The candidate is annotated as {annotation or family or 'an uncharacterized protein'} "
        "with supporting sequence and functional evidence. "
        f"{expression_sentence} "
        f"The phylogenetic analysis provides {phylogenetic_group.lower() if phylogenetic_group else 'additional evolutionary support'}. "
        "The biological effect of RNAi remains to be confirmed experimentally."
    )

print("\n")
print("=" * 80)
print("FINAL TOP-10 CANDIDATE RANKING")
print("=" * 80)

print(
    "\nThe Top 100 RNA file is the authoritative candidate pool."
    "\nBLAST, InterPro and phylogeny provide supporting evidence."
)

rna_path = ask_file("TOP 100 RNA TSV")
blast_path = ask_file("TOP 500 BLAST-UniProt TSV")
interpro_path = ask_file("InterPro TSV")
phylo_path = ask_file("Phylogenetic annotation summary TSV")

output_dir = ask_output_directory()

if output_dir is None:
    output_dir = rna_path.parent

output_dir.mkdir(parents=True, exist_ok=True)

rna = load_tsv(rna_path, "Top 100 RNA")
blast = load_tsv(blast_path, "BLAST")
interpro = load_tsv(interpro_path, "InterPro")
phylo = load_tsv(
    phylo_path,
    "Phylogenetic annotation summary"
)

require_columns(
    rna,
    [
        "gene_id",
        "RNA_score",
        "Larva_TPM",
        "Adult_TPM",
        "Pupa_TPM",
        "Egg_TPM",
        "Larva_expression"
    ],
    "Top 100 RNA file"
)

require_columns(
    blast,
    [
        "gene_id",
        "hit_rank",
        "sseqid",
        "pident",
        "qcov",
        "evalue",
        "bitscore",
        "stitle"
    ],
    "BLAST file"
)

require_columns(
    interpro,
    [
        "Protein Accession",
        "Signature accession",
        "InterPro accession"
    ],
    "InterPro file"
)

require_columns(
    phylo,
    [
        "gene_id",
        "protein_id",
        "phylo_success",
        "phylo_sequence_count",
        "alignment_length",
        "trimmed_alignment_length",
        "iqtree_model",
        "ufboot_mean_pct",
        "sh_alrt_mean_pct"
    ],
    "Phylogenetic annotation summary"
)

if "phylogenetic_group" in phylo.columns:
    print("\n✓ phylogenetic_group detected.")
    print("  Named phylogenetic groups will be used.")
else:
    print("\nNOTE: phylogenetic_group column not found.")
    print(
        "      UFBoot/SH-aLRT support will be reported,"
        "\n      but no named clade will be invented."
    )

rna["gene_key"] = rna["gene_id"].map(normalize_gene_id)
blast["gene_key"] = blast["gene_id"].map(normalize_gene_id)
interpro["gene_key"] = (
    interpro["Protein Accession"].map(normalize_gene_id)
)
phylo["gene_key"] = phylo["gene_id"].map(normalize_gene_id)

rna = rna[
    rna["gene_key"].astype(str).str.strip() != ""
].copy()

rna = rna.drop_duplicates(
    subset=["gene_key"],
    keep="first"
).copy()

print(
    f"\nAuthoritative candidate pool: "
    f"{len(rna):,} genes"
)

candidate_genes = set(rna["gene_key"])

blast = blast[
    blast["gene_key"].isin(candidate_genes)
].copy()

interpro = interpro[
    interpro["gene_key"].isin(candidate_genes)
].copy()

phylo = phylo[
    phylo["gene_key"].isin(candidate_genes)
].copy()

rna["RNA_score_numeric"] = numeric(rna, "RNA_score")
rna["Larva_TPM_numeric"] = numeric(rna, "Larva_TPM")
rna["Adult_TPM_numeric"] = numeric(rna, "Adult_TPM")
rna["Pupa_TPM_numeric"] = numeric(rna, "Pupa_TPM")
rna["Egg_TPM_numeric"] = numeric(rna, "Egg_TPM")

blast["hit_rank_numeric"] = numeric(blast, "hit_rank")
blast["pident_numeric"] = numeric(blast, "pident")
blast["qcov_numeric"] = numeric(blast, "qcov")
blast["evalue_numeric"] = numeric(blast, "evalue")
blast["bitscore_numeric"] = numeric(blast, "bitscore")

phylo["ufboot_mean_pct"] = numeric(
    phylo,
    "ufboot_mean_pct"
)

phylo["sh_alrt_mean_pct"] = numeric(
    phylo,
    "sh_alrt_mean_pct"
)

blast = blast.sort_values(
    by=[
        "gene_key",
        "hit_rank_numeric",
        "evalue_numeric",
        "bitscore_numeric",
        "pident_numeric",
        "qcov_numeric"
    ],
    ascending=[
        True, True, True, False, False, False
    ],
    na_position="last"
)

best_blast = (
    blast
    .drop_duplicates(
        subset=["gene_key"],
        keep="first"
    )
    .copy()
)

print(f"\nBest BLAST hits selected: {len(best_blast):,}")

old_columns = [
    "best_hit_sseqid", "best_hit_accession",
    "best_hit_species", "best_hit_description",
    "best_identity_pct", "best_query_coverage_pct",
    "best_evalue", "best_bitscore",
    "functional_annotation", "gene_family",
    "subfamily_or_clade", "phylogenetic_group",
    "pfam_domains", "interpro_domains",
    "target_class", "offtarget_check"
]

rna = rna.drop(
    columns=[c for c in old_columns if c in rna.columns],
    errors="ignore"
)

best_blast_small = best_blast[
    [
        "gene_key",
        "sseqid",
        "pident_numeric",
        "qcov_numeric",
        "evalue_numeric",
        "bitscore_numeric",
        "stitle"
    ]
].copy()

best_blast_small = best_blast_small.rename(
    columns={
        "sseqid": "best_hit_sseqid",
        "pident_numeric": "best_identity_pct",
        "qcov_numeric": "best_query_coverage_pct",
        "evalue_numeric": "best_evalue",
        "bitscore_numeric": "best_bitscore",
        "stitle": "best_hit_description"
    }
)

merged = rna.merge(
    best_blast_small,
    on="gene_key",
    how="left",
    validate="one_to_one"
)

phylo_columns = [
    "gene_key",
    "protein_id",
    "phylo_success",
    "phylo_sequence_count",
    "alignment_length",
    "trimmed_alignment_length",
    "iqtree_model",
    "ufboot_mean_pct",
    "sh_alrt_mean_pct"
]

if "phylogenetic_group" in phylo.columns:

    phylo_columns.append(
        "phylogenetic_group"
    )


phylo_small = (
    phylo[
        phylo_columns
    ]
    .drop_duplicates(
        subset=["gene_key"],
        keep="first"
    )
    .copy()
)


merged = merged.merge(
    phylo_small,
    on="gene_key",
    how="left",
    validate="one_to_one"
)

interpro_records = {}

for gene, group in interpro.groupby(
    "gene_key"
):

    pfam = extract_pfam(
        group
    )

    ipr = extract_interpro(
        group
    )

    ipr_desc = extract_interpro_descriptions(
        group
    )

    interpro_records[gene] = {

        "pfam":
            pfam,

        "interpro":
            ipr,

        "interpro_desc":
            ipr_desc,

        "interpro_count":
            len([
                x
                for x in ipr.split(";")
                if x
            ]),

        "pfam_count":
            len([
                x
                for x in pfam.split(";")
                if x
            ])
    }

rows = []

for _, row in merged.iterrows():

    gene = row["gene_key"]

    ipr_info = interpro_records.get(
        gene,
        {
            "pfam": "",
            "interpro": "",
            "interpro_desc": "",
            "interpro_count": 0,
            "pfam_count": 0
        }
    )

    raw_function = extract_function(
        row.get("best_hit_description", "")
    )

    gene_family = determine_gene_family(
        raw_function,
        ipr_info["interpro_desc"],
        ipr_info["pfam"]
    )

    existing_target = row.get("target_class", "")

    target_class = determine_target_class(
        raw_function,
        gene_family,
        ipr_info["interpro_desc"],
        existing_target
    )

    functional_annotation = create_functional_annotation(
        raw_function,
        gene_family,
        target_class
    )

    identity = row.get("best_identity_pct", np.nan)
    coverage = row.get("best_query_coverage_pct", np.nan)
    evalue = row.get("best_evalue", np.nan)
    bitscore = row.get("best_bitscore", np.nan)

    homology_score = calculate_homology_score(
        identity,
        coverage,
        evalue
    )

    annotation_score = calculate_annotation_score(
        ipr_info["interpro_count"],
        ipr_info["pfam_count"],
        functional_annotation
    )

    phylo_success = (
        str(row.get("phylo_success", ""))
        .strip()
        .lower()
        in {"true", "1", "yes", "success"}
    )

    phylo_score = calculate_phylo_score(
        {
            "phylo_success": phylo_success,
            "ufboot_mean_pct": row.get(
                "ufboot_mean_pct",
                np.nan
            ),
            "sh_alrt_mean_pct": row.get(
                "sh_alrt_mean_pct",
                np.nan
            )
        }
    )

    phylogenetic_group = determine_phylogenetic_group(
        {
            "phylo_success": phylo_success,
            "phylogenetic_group": row.get(
                "phylogenetic_group",
                ""
            ),
            "ufboot_mean_pct": row.get(
                "ufboot_mean_pct",
                np.nan
            ),
            "sh_alrt_mean_pct": row.get(
                "sh_alrt_mean_pct",
                np.nan
            )
        }
    )

    phylo_text = phylo_evidence(
        {
            "phylo_success": phylo_success,
            "iqtree_model": row.get(
                "iqtree_model",
                ""
            ),
            "ufboot_mean_pct": row.get(
                "ufboot_mean_pct",
                np.nan
            ),
            "sh_alrt_mean_pct": row.get(
                "sh_alrt_mean_pct",
                np.nan
            )
        }
    )

    rna_score = calculate_rna_score(row)

    larva_tpm = row.get("Larva_TPM_numeric", np.nan)
    adult_tpm = row.get("Adult_TPM_numeric", np.nan)
    pupa_tpm = row.get("Pupa_TPM_numeric", np.nan)
    egg_tpm = row.get("Egg_TPM_numeric", np.nan)

    all_stage_tpm = [
        larva_tpm,
        adult_tpm,
        pupa_tpm,
        egg_tpm
    ]

    larva_expression_score = calculate_stage_expression_score(
        larva_tpm,
        all_stage_tpm
    )

    adult_expression_score = calculate_stage_expression_score(
        adult_tpm,
        all_stage_tpm
    )

    pupa_expression_score = calculate_stage_expression_score(
        pupa_tpm,
        all_stage_tpm
    )

    egg_expression_score = calculate_stage_expression_score(
        egg_tpm,
        all_stage_tpm
    )

    target_score = calculate_target_score(
        target_class,
        functional_annotation
    )

    final_score = (
        rna_score * WEIGHT_RNA / 100
        + homology_score * WEIGHT_HOMOLOGY / 100
        + annotation_score * WEIGHT_ANNOTATION / 100
        + phylo_score * WEIGHT_PHYLO / 100
        + target_score * WEIGHT_TARGET / 100
    )

    best_species = extract_species(
        row.get("best_hit_description", "")
    )

    protein_id = clean_text(
        row.get("protein_id", "")
    )

    if not protein_id:
        protein_id = clean_text(
            row.get("gene_id", "")
        )

    larval_expression = clean_text(
        row.get("Larva_expression", "")
    )

    justification = generate_justification(
        gene_family,
        functional_annotation,
        target_class,
        larva_tpm,
        adult_tpm,
        pupa_tpm,
        egg_tpm,
        phylogenetic_group
    )

    rows.append({
        "_gene_key": gene,

        "gene_id": clean_text(
            row.get("gene_id", "")
        ),

        "protein_id": protein_id,

        "protein_length_aa": clean_text(
            row.get("sequence_length", "")
        ),

        "gene_family": gene_family,
        "phylogenetic_group": phylogenetic_group,

        "best_hit_accession": extract_accession(
            row.get("best_hit_sseqid", "")
        ),

        "best_hit_species": best_species,

        "pfam_domains": ipr_info["pfam"],
        "interpro_domains": ipr_info["interpro"],

        "functional_annotation": functional_annotation,
        "target_class": target_class,

        "evalue": (
            evalue if pd.notna(evalue) else ""
        ),

        "identity_pct": (
            identity if pd.notna(identity) else ""
        ),

        "query_coverage_pct": (
            coverage if pd.notna(coverage) else ""
        ),

        "phylogenetic_evidence": phylo_text,

        "Larva_TPM": (
            larva_tpm if pd.notna(larva_tpm) else ""
        ),

        "Adult_TPM": (
            adult_tpm if pd.notna(adult_tpm) else ""
        ),

        "Pupa_TPM": (
            pupa_tpm if pd.notna(pupa_tpm) else ""
        ),

        "Egg_TPM": (
            egg_tpm if pd.notna(egg_tpm) else ""
        ),

        "larva_expression_score":
            larva_expression_score,

        "adult_expression_score":
            adult_expression_score,

        "pupa_expression_score":
            pupa_expression_score,

        "egg_expression_score":
            egg_expression_score,

        "Larva_expression": larval_expression,

        "rnai_or_chemistry": "RNAi",

        "justification": justification,

        # INTERNAL

        "_rna_score": rna_score,
        "_homology_score": homology_score,
        "_annotation_score": annotation_score,
        "_phylo_score": phylo_score,
        "_target_score": target_score,
        "_final_score": final_score,

        "_rna_rank": row.get(
            "RNA_rank_numeric",
            np.nan
        )
    })

results = pd.DataFrame(rows)

results = results.sort_values(
    by=[
        "_final_score",
        "_rna_score",
        "_homology_score",
        "_annotation_score",
        "_phylo_score",
        "_target_score",
        "larva_expression_score"
    ],

    ascending=[
        False,
        False,
        False,
        False,
        False,
        False,
        False
    ],

    na_position="last"
).reset_index(
    drop=True
)

results.insert(
    0,
    "rank",
    range(
        1,
        len(results) + 1
    )
)

top10 = results.head(
    TOP_N
).copy()

internal_columns = [
    c
    for c in results.columns
    if c.startswith("_")
]

results_public = results.drop(
    columns=internal_columns
)

top10_public = top10.drop(
    columns=internal_columns
)

final_columns = [

    "rank",
    "gene_id",
    "protein_id",
    "protein_length_aa",

    "gene_family",
    "phylogenetic_group",

    "best_hit_accession",
    "best_hit_species",

    "pfam_domains",
    "interpro_domains",

    "functional_annotation",
    "target_class",

    "evalue",
    "identity_pct",
    "query_coverage_pct",

    "phylogenetic_evidence",

    "Larva_TPM",
    "Adult_TPM",
    "Pupa_TPM",
    "Egg_TPM",

    "larva_expression_score",
    "adult_expression_score",
    "pupa_expression_score",
    "egg_expression_score",

    "Larva_expression",

    "rnai_or_chemistry",

    "justification"
]


top10_public = top10_public[
    final_columns
]

results_public = results_public[
    final_columns
]

top10_tsv = (
    output_dir
    / "FINAL_TOP10_CANDIDATES.tsv"
)

top10_csv = (
    output_dir
    / "FINAL_TOP10_CANDIDATES.csv"
)

full_tsv = (
    output_dir
    / "RANKED_TOP100_CANDIDATES.tsv"
)

top10_public.to_csv(
    top10_tsv,
    sep="\t",
    index=False
)

top10_public.to_csv(
    top10_csv,
    index=False
)

results_public.to_csv(
    full_tsv,
    sep="\t",
    index=False
)

print("\n" + "=" * 80)
print("FINAL TOP 10")
print("=" * 80)

display_columns = [
    "rank",
    "gene_id",
    "gene_family",
    "phylogenetic_group",
    "target_class"
]

print(
    top10_public[
        display_columns
    ].to_string(
        index=False
    )
)

print("\n" + "=" * 80)
print("TOP 10 SCORE BREAKDOWN")
print("=" * 80)

score_display = results[
    [
        "rank",
        "gene_id",
        "_rna_score",
        "_homology_score",
        "_annotation_score",
        "_phylo_score",
        "_target_score",
        "_final_score"
    ]
].head(
    TOP_N
).copy()

score_display.columns = [
    "rank",
    "gene_id",
    "RNA_score",
    "Homology_score",
    "Annotation_score",
    "Phylo_score",
    "Target_score",
    "Final_score"
]

print(
    score_display.to_string(
        index=False,
        float_format=lambda x: f"{x:.2f}"
    )
)

print("\n" + "=" * 80)
print("OUTPUT COMPLETE")
print("=" * 80)

print(f"\nAuthoritative candidate pool: {len(rna)} genes")
print(f"Final ranked candidates: {len(results_public)}")
print(f"Final Top 10 candidates: {len(top10_public)}")

print("\nTop 10 TSV:")
print(top10_tsv)

print("\nTop 10 CSV:")
print(top10_csv)

print("\nFull ranked Top 100 TSV:")
print(full_tsv)

print("\n" + "=" * 80)
print("FINAL VALIDATION")
print("=" * 80)

# Candidate count
assert len(rna) == len(results_public), (
    f"RNA candidate count ({len(rna)}) does not match "
    f"ranked candidate count ({len(results_public)})"
)

# Top 10
assert len(top10_public) == 10, (
    f"Expected 10 final candidates, found {len(top10_public)}"
)

# Unique genes
assert results_public["gene_id"].nunique() == len(results_public), (
    "Duplicate gene IDs detected."
)

# Unique ranks
assert results_public["rank"].nunique() == len(results_public), (
    "Duplicate ranks detected."
)

# Rank sequence
assert list(results_public["rank"]) == list(
    range(1, len(results_public) + 1)
), "Rank sequence is invalid."

# Functional annotation length
annotation_word_counts = (
    results_public["functional_annotation"]
    .fillna("")
    .apply(lambda x: len(str(x).split()))
)

assert annotation_word_counts.between(7, 15).all(), (
    "Functional annotation contains entries outside the 7–15 word range."
)

# BLAST coverage
blast_coverage = merged["best_hit_sseqid"].notna().sum()

print(f"\n✓ RNA candidate pool = {len(rna)}")
print(f"✓ Ranked candidates = {len(results_public)}")
print(f"✓ Final Top 10 = {len(top10_public)}")
print(f"✓ Unique gene IDs = {results_public['gene_id'].nunique()}")
print(f"✓ Candidates with BLAST evidence = {blast_coverage}")
print("✓ BLAST hit_rank=1 used as best hit")
print("✓ InterPro/Pfam used as annotation evidence")
print("✓ MAFFT + trimAl + IQ-TREE used as phylogenetic evidence")
print("✓ UFBoot and SH-aLRT retained")
print("✓ Named phylogenetic groups used only when provided")
print("✓ No unsupported CHS-A/CHS-B assignments invented")
print("✓ Stage-specific RNA expression retained")
print("✓ Off-target column removed")
print("✓ Functional annotations constrained to 7–15 words")
print("✓ Justifications generated as natural 3–4 sentence explanations")
print("✓ Exactly 10 candidates exported")
print("\nDone.")