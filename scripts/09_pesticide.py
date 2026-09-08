
# ============================================================
# METISA PLANA — PESTICIDE / COMPOUND EVIDENCE + TOP 100 / TOP 10
# ============================================================
#
# PURPOSE:
#   Combine:
#       1. Top-500 pesticide / compound evidence mapping
#       2. Top-100 pesticide evidence addition
#       3. Top-10 pesticide evidence addition
#
# INPUT:
#   1. Top-500 candidate TSV
#   2. Top-100 candidate TSV
#   3. Top-10 candidate TSV
#
# OUTPUT:
#
#   top500_pesticide_evidence.tsv
#   top500_pesticide_evidence_summary.tsv
#   top500_target_discovery.tsv
#   top500_target_discovery_raw.tsv
#
#   top100_with_pesticide.tsv
#
#   top10_with_pesticide.csv
#   top10_with_pesticide.tsv
#
# IMPORTANT:
#   - Top 100 and Top 10 are NOT derived from Top 500.
#   - They are supplied separately because they may contain
#     different ranking/results.
#   - Pesticide evidence is generated from the Top 500 and
#     then mapped to Top 100 / Top 10 using gene_id.
#   - Original ranking, RNA, phylogeny, scoring and
#     justification columns are preserved.
#   - Pesticide evidence is literature-backed and generally
#     based on insect homologues or target families.
#   - It does NOT prove that the corresponding Metisa plana
#     protein is itself an experimentally validated pesticide
#     target.
#
# ============================================================

from pathlib import Path

import pandas as pd
import re
import sys


# ============================================================
# VERIFIED EVIDENCE LIBRARY
# ============================================================
#
# These are literature-backed evidence records.
#
# They are NOT scoring rules.
#
# Evidence levels:
#
# DIRECT_TARGET
#     Experimental evidence connects the target protein/family
#     to an insecticidal compound or toxin.
#
# INSECTICIDAL_INHIBITOR
#     Experimental inhibitor affects insect growth,
#     development or mortality, but the exact M. plana
#     protein is not experimentally validated.
#
# SYNERGIST_RESISTANCE
#     Mainly relevant to insecticide detoxification,
#     resistance or synergism rather than direct target
#     validation.
#
# Evidence scope:
#
# Metisa_plana_direct
#     Exact M. plana protein experimentally validated.
#
# Insect_homologue
#     Evidence demonstrated using homologous proteins
#     in other insect species.
#
# Target_family
#     Evidence supports the broader target/enzyme family.
#
# Resistance_synergist
#     Evidence concerns resistance or insecticide synergism.
#
# ============================================================

EVIDENCE_LIBRARY = [

    # --------------------------------------------------------
    # EV01 — TRYPSIN / CHYMOTRYPSIN
    # --------------------------------------------------------

    {
        "evidence_id":
            "EV01",

        "target_family":
            "Trypsin / chymotrypsin serine proteases",

        "match":
            r"\btrypsin\b|\bchymotrypsin\b|serine protease",

        "compound":
            "Bowman-Birk proteinase inhibitors",

        "compound_class":
            "Proteinase inhibitor",

        "compound_type":
            "Biological / experimental insecticidal inhibitor",

        "mechanism":
            "Inhibits insect digestive trypsin/chymotrypsin-like proteases",

        "evidence_level":
            "INSECTICIDAL_INHIBITOR",

        "evidence_scope":
            "Insect_homologue",

        "evidence_species":
            "Manduca sexta; Anagasta kuehniella; "
            "Diatraea saccharalis; Heliothis virescens",

        "evidence_summary":
            "A Bowman-Birk proteinase inhibitor inhibited larval "
            "midgut trypsin-like proteinase activity in Manduca sexta "
            "and reduced larval growth and development.",

        "reference":
            "PMID 20146519; PMID 26330217",

        "source_url":
            "https://pubmed.ncbi.nlm.nih.gov/20146519/ ; "
            "https://pubmed.ncbi.nlm.nih.gov/26330217/",
    },


    # --------------------------------------------------------
    # EV02 — LIPASE
    # --------------------------------------------------------

    {
        "evidence_id":
            "EV02",

        "target_family":
            "Midgut lipases / triacylglycerol lipases",

        "match":
            r"\blipase\b",

        "compound":
            "Tetrahydrolipstatin (THL; orlistat)",

        "compound_class":
            "Lipase inhibitor",

        "compound_type":
            "Experimental insecticidal / growth-disrupting inhibitor",

        "mechanism":
            "Inhibits midgut lipase activity",

        "evidence_level":
            "INSECTICIDAL_INHIBITOR",

        "evidence_scope":
            "Insect_homologue",

        "evidence_species":
            "Epiphyas postvittana",

        "evidence_summary":
            "Tetrahydrolipstatin (THL) reduced midgut lipase activity "
            "in Epiphyas postvittana and was associated with reduced "
            "larval growth, development and pupation.",

        "reference":
            "PMID 21910995; DOI 10.1016/j.jinsphys.2011.08.018",

        "source_url":
            "https://pubmed.ncbi.nlm.nih.gov/21910995/",
    },


    # --------------------------------------------------------
    # EV03 — AMINOPEPTIDASE N
    # --------------------------------------------------------

    {
        "evidence_id":
            "EV03",

        "target_family":
            "Aminopeptidase N (APN)",

        "match":
            r"aminopeptidase\s+n\b|aminopeptidase n-type",

        "compound":
            "Bacillus thuringiensis Cry toxins",

        "compound_class":
            "Bt insecticidal crystal toxin",

        "compound_type":
            "Biological insecticide",

        "mechanism":
            "APN functions as a midgut receptor/binding protein "
            "for specific Cry toxins in several insects",

        "evidence_level":
            "DIRECT_TARGET",

        "evidence_scope":
            "Insect_homologue",

        "evidence_species":
            "Manduca sexta; Bombyx mori; Lymantria dispar; "
            "Spodoptera litura; Aedes aegypti",

        "evidence_summary":
            "APN proteins have been experimentally identified as "
            "receptors or toxin-binding proteins for specific Bt Cry "
            "toxins in several insect species.",

        "reference":
            "PMID 7908713; PMID 7629076; PMID 8580914; "
            "PMID 12200317; PMID 24128608",

        "source_url":
            "https://pubmed.ncbi.nlm.nih.gov/7908713/ ; "
            "https://pubmed.ncbi.nlm.nih.gov/7629076/ ; "
            "https://pubmed.ncbi.nlm.nih.gov/8580914/ ; "
            "https://pubmed.ncbi.nlm.nih.gov/12200317/ ; "
            "https://pubmed.ncbi.nlm.nih.gov/24128608/",
    },


    # --------------------------------------------------------
    # EV04 — CHITINASE
    # --------------------------------------------------------

    {
        "evidence_id":
            "EV04",

        "target_family":
            "Insect chitinases",

        "match":
            r"\bchitinase\b",

        "compound":
            "Allosamidin",

        "compound_class":
            "Chitinase inhibitor",

        "compound_type":
            "Experimental insecticidal compound",

        "mechanism":
            "Competitive inhibition of insect family-18 "
            "chitinases; disrupts moulting",

        "evidence_level":
            "INSECTICIDAL_INHIBITOR",

        "evidence_scope":
            "Target_family",

        "evidence_species":
            "Bombyx mori; Leucania separata; Lucilia cuprina; "
            "Tineola bisselliella; Myzus persicae",

        "evidence_summary":
            "Allosamidin specifically inhibits insect chitinase "
            "activity and has been reported to interfere with "
            "moulting and produce insecticidal effects.",

        "reference":
            "PMID 3570982; DOI 10.1271/bbb1961.51.471; "
            "PMID 25486024",

        "source_url":
            "https://pubmed.ncbi.nlm.nih.gov/3570982/ ; "
            "https://www.jstage.jst.go.jp/article/bbb1961/51/2/51_2_471/_article ; "
            "https://pubmed.ncbi.nlm.nih.gov/25486024/",
    },


    # --------------------------------------------------------
    # EV05 — ACE / PEPTIDYL-DIPEPTIDASE A
    # --------------------------------------------------------

    {
        "evidence_id":
            "EV05",

        "target_family":
            "Peptidyl-dipeptidase A / ACE-like peptidases",

        "match":
            r"peptidyl[- ]dipeptidase\s+a|"
            r"angiotensin[- ]converting|ace[- ]like",

        "compound":
            "Captopril; fosinopril/fosinoprilat",

        "compound_class":
            "ACE inhibitor",

        "compound_type":
            "Experimental larvicide",

        "mechanism":
            "Inhibits larval ACE-like peptidyl-dipeptidase activity",

        "evidence_level":
            "INSECTICIDAL_INHIBITOR",

        "evidence_scope":
            "Insect_homologue",

        "evidence_species":
            "Aedes aegypti; Anopheles gambiae",

        "evidence_summary":
            "ACE inhibitors inhibited larval peptidyl-dipeptidase "
            "activity; captopril and fosinopril showed larvicidal activity.",

        "reference":
            "PMID 28345667",

        "source_url":
            "https://pubmed.ncbi.nlm.nih.gov/28345667/",
    },


    # --------------------------------------------------------
    # EV06 — NICOTINIC ACETYLCHOLINE RECEPTORS
    # --------------------------------------------------------

    {
        "evidence_id":
            "EV06",

        "target_family":
            "Nicotinic acetylcholine receptors (nAChRs)",

        "match":
            r"nicotinic acetylcholine receptor|\bnachr\b",

        "compound":
            "Neonicotinoids (e.g. imidacloprid, clothianidin, "
            "thiamethoxam); spinosyns; sulfoximines; "
            "related nAChR insecticides",

        "compound_class":
            "nAChR-active insecticides",

        "compound_type":
            "Commercial insecticides",

        "mechanism":
            "Agonist/modulatory action at insect "
            "nicotinic acetylcholine receptors",

        "evidence_level":
            "DIRECT_TARGET",

        "evidence_scope":
            "Target_family",

        "evidence_species":
            "Multiple insect species",

        "evidence_summary":
            "Insect nAChRs are established molecular targets "
            "of several commercial insecticide classes including "
            "neonicotinoids and spinosyns.",

        "reference":
            "PMID 17216290; PMID 11698101; PMID 28176635; "
            "IRAC Mode of Action",

        "source_url":
            "https://pubmed.ncbi.nlm.nih.gov/17216290/ ; "
            "https://pubmed.ncbi.nlm.nih.gov/11698101/ ; "
            "https://pubmed.ncbi.nlm.nih.gov/28176635/ ; "
            "https://irac-online.org/mode-of-action/",
    },


    # --------------------------------------------------------
    # EV07 — GST
    # --------------------------------------------------------

    {
        "evidence_id":
            "EV07",

        "target_family":
            "Glutathione S-transferases (GSTs)",

        "match":
            r"glutathione s-transferase|\bgst\b",

        "compound":
            "Taxifolin; quercetin; diethyl maleate; "
            "S-hexyl glutathione",

        "compound_class":
            "GST inhibitor / insecticide synergist",

        "compound_type":
            "Resistance-management / synergist evidence",

        "mechanism":
            "Inhibition of GST-mediated detoxification can "
            "increase insecticide susceptibility",

        "evidence_level":
            "SYNERGIST_RESISTANCE",

        "evidence_scope":
            "Resistance_synergist",

        "evidence_species":
            "Leptinotarsa decemlineata; Plutella xylostella; "
            "other insects",

        "evidence_summary":
            "GSTs are detoxification enzymes associated with "
            "insecticide resistance; GST inhibitors have been "
            "investigated as insecticide synergists.",

        "reference":
            "PMID 25270601; PMID 30207467; PMID 34301351",

        "source_url":
            "https://pubmed.ncbi.nlm.nih.gov/25270601/ ; "
            "https://pubmed.ncbi.nlm.nih.gov/30207467/ ; "
            "https://pubmed.ncbi.nlm.nih.gov/34301351/",
    },


    # --------------------------------------------------------
    # EV08 — CYTOCHROME P450
    # --------------------------------------------------------

    {
        "evidence_id":
            "EV08",

        "target_family":
            "Cytochrome P450 monooxygenases",

        "match":
            r"cytochrome p450|p450 superfamily",

        "compound":
            "Piperonyl butoxide (PBO)",

        "compound_class":
            "P450 inhibitor / insecticide synergist",

        "compound_type":
            "Resistance-management / synergist",

        "mechanism":
            "Inhibits cytochrome P450-mediated "
            "xenobiotic metabolism",

        "evidence_level":
            "SYNERGIST_RESISTANCE",

        "evidence_scope":
            "Resistance_synergist",

        "evidence_species":
            "Drosophila melanogaster and other insects",

        "evidence_summary":
            "PBO is an established insecticide synergist "
            "that inhibits insect cytochrome P450 activity.",

        "reference":
            "PMID 17514638",

        "source_url":
            "https://pubmed.ncbi.nlm.nih.gov/17514638/",
    },

]


# ============================================================
# REQUIRED TOP-500 COLUMNS
# ============================================================

REQUIRED_TOP500_COLUMNS = [

    "gene_id",
    "rank",
    "functional_description",
    "best_hit_description",
    "target_class",
    "interpro_accessions",
    "interpro_descriptions",
    "pfam",
    "panther",
    "go_terms",
    "functional_keywords",

]


# ============================================================
# EXPECTED TOP-100 / TOP-10 COLUMNS
# ============================================================

EXPECTED_RANKING_COLUMNS = [

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
    "justification",

]


# ============================================================
# BASIC TEXT CLEANING
# ============================================================

def clean_text(value):

    if pd.isna(value):
        return ""

    return str(value).strip()


def normalize_spaces(text):

    return re.sub(
        r"\s+",
        " ",
        clean_text(text)
    ).strip()


# ============================================================
# GENE ID CLEANING
# ============================================================

def clean_gene_id(value):

    if pd.isna(value):
        return ""

    value = str(value).strip()

    # Remove accidental whitespace
    value = value.replace(
        " ",
        ""
    )

    # Convert:
    #
    # g13783.t1
    #
    # into:
    #
    # g13783
    value = re.sub(
        r"\.t\d+$",
        "",
        value,
        flags=re.IGNORECASE
    )

    return value


# ============================================================
# SPLIT MULTI-VALUE ANNOTATIONS
# ============================================================

def split_multi(value):

    text = clean_text(value)

    if not text:
        return []

    parts = re.split(
        r"\s*[;|]\s*",
        text
    )

    return [
        normalize_spaces(x)
        for x in parts
        if normalize_spaces(x)
    ]


# ============================================================
# FIND EVIDENCE MATCHES
# ============================================================

def evidence_matches(row):

    # --------------------------------------------------------
    # IMPORTANT:
    #
    # Evidence matching is based ONLY on:
    #
    #   functional_description
    #   best_hit_description
    #   target_class
    #
    # GO / Pfam / InterPro are NOT used alone to assign
    # pesticide evidence.
    # --------------------------------------------------------

    searchable_text = " | ".join([

        clean_text(
            row.get(
                "functional_description",
                ""
            )
        ),

        clean_text(
            row.get(
                "best_hit_description",
                ""
            )
        ),

        clean_text(
            row.get(
                "target_class",
                ""
            )
        ),

    ])

    matches = []

    for evidence in EVIDENCE_LIBRARY:

        pattern = evidence["match"]

        try:

            if re.search(
                pattern,
                searchable_text,
                flags=re.IGNORECASE
            ):

                matches.append(
                    evidence
                )

        except re.error as error:

            print(
                f"WARNING: Invalid evidence pattern "
                f"{evidence['evidence_id']}: {error}"
            )

    return matches


# ============================================================
# CREATE PESTICIDE EVIDENCE STRING
# ============================================================

def make_evidence_string(matches):

    if not matches:

        return (
            "No verified pesticide/compound evidence mapped"
        )

    evidence_strings = []

    seen = set()

    for evidence in matches:

        text = (
            f'{evidence["target_family"]}: '
            f'{evidence["compound"]} '
            f'[{evidence["evidence_level"]}]'
        )

        if text not in seen:

            evidence_strings.append(
                text
            )

            seen.add(
                text
            )

    return " || ".join(
        evidence_strings
    )


# ============================================================
# CREATE EVIDENCE LEVEL STRING
# ============================================================

def make_evidence_level_string(matches):

    if not matches:
        return "NONE"

    levels = []

    for evidence in matches:

        level = evidence[
            "evidence_level"
        ]

        if level not in levels:

            levels.append(
                level
            )

    return "; ".join(
        sorted(levels)
    )


# ============================================================
# CREATE EVIDENCE SCOPE STRING
# ============================================================

def make_evidence_scope_string(matches):

    if not matches:
        return "NONE"

    scopes = []

    for evidence in matches:

        scope = evidence[
            "evidence_scope"
        ]

        if scope not in scopes:

            scopes.append(
                scope
            )

    return "; ".join(
        sorted(scopes)
    )


# ============================================================
# BUILD TOP-500 EVIDENCE TABLE
# ============================================================

def build_top500_evidence(df):

    output = df.copy()

    pesticide_evidence = []

    pesticide_levels = []

    pesticide_scopes = []

    for _, row in df.iterrows():

        matches = evidence_matches(
            row
        )

        # ----------------------------------------------------
        # Evidence text
        # ----------------------------------------------------

        pesticide_evidence.append(
            make_evidence_string(
                matches
            )
        )

        # ----------------------------------------------------
        # Evidence level
        # ----------------------------------------------------

        pesticide_levels.append(
            make_evidence_level_string(
                matches
            )
        )

        # ----------------------------------------------------
        # Evidence scope
        # ----------------------------------------------------

        pesticide_scopes.append(
            make_evidence_scope_string(
                matches
            )
        )

    output[
        "pesticide_evidence"
    ] = pesticide_evidence

    output[
        "pesticide_evidence_level"
    ] = pesticide_levels

    output[
        "pesticide_evidence_scope"
    ] = pesticide_scopes

    return output


# ============================================================
# BUILD PESTICIDE EVIDENCE SUMMARY
# ============================================================

def build_pesticide_summary(df):

    rows = []

    for _, row in df.iterrows():

        matches = evidence_matches(
            row
        )

        for evidence in matches:

            rows.append({

                "evidence_id":
                    evidence["evidence_id"],

                "target_family":
                    evidence["target_family"],

                "compound":
                    evidence["compound"],

                "compound_class":
                    evidence["compound_class"],

                "compound_type":
                    evidence["compound_type"],

                "mechanism":
                    evidence["mechanism"],

                "evidence_level":
                    evidence["evidence_level"],

                "evidence_scope":
                    evidence["evidence_scope"],

                "evidence_species":
                    evidence["evidence_species"],

                "gene_id":
                    clean_text(
                        row.get(
                            "gene_id",
                            ""
                        )
                    ),

                "rank":
                    clean_text(
                        row.get(
                            "rank",
                            ""
                        )
                    ),

                "functional_description":
                    clean_text(
                        row.get(
                            "functional_description",
                            ""
                        )
                    ),

                "best_hit_description":
                    clean_text(
                        row.get(
                            "best_hit_description",
                            ""
                        )
                    ),

                "target_class":
                    clean_text(
                        row.get(
                            "target_class",
                            ""
                        )
                    ),

                "evidence_summary":
                    evidence["evidence_summary"],

                "reference":
                    evidence["reference"],

                "source_url":
                    evidence["source_url"],

            })


    raw_evidence = pd.DataFrame(
        rows
    )


    # --------------------------------------------------------
    # No matches
    # --------------------------------------------------------

    if raw_evidence.empty:

        return pd.DataFrame(
            columns=[

                "evidence_id",
                "target_family",
                "compound",
                "compound_class",
                "compound_type",
                "mechanism",
                "evidence_level",
                "evidence_scope",
                "evidence_species",
                "candidate_count",
                "genes",
                "ranks",
                "evidence_summary",
                "reference",
                "source_url",

            ]
        )


    # --------------------------------------------------------
    # Group evidence
    # --------------------------------------------------------

    summary = (

        raw_evidence

        .groupby(
            [

                "evidence_id",
                "target_family",
                "compound",
                "compound_class",
                "compound_type",
                "mechanism",
                "evidence_level",
                "evidence_scope",
                "evidence_species",
                "evidence_summary",
                "reference",
                "source_url",

            ],
            dropna=False
        )

        .agg(

            candidate_count=(
                "gene_id",
                "nunique"
            ),

            genes=(
                "gene_id",
                lambda x:
                    "; ".join(
                        sorted(
                            set(
                                x
                            )
                        )
                    )
            ),

            ranks=(
                "rank",
                lambda x:
                    "; ".join(
                        sorted(
                            set(
                                map(
                                    str,
                                    x
                                )
                            )
                        )
                    )
            ),

        )

        .reset_index()
    )


    summary = summary.sort_values(

        [

            "evidence_level",
            "candidate_count",
            "target_family",

        ],

        ascending=[

            True,
            False,
            True,

        ]
    )

    return summary


# ============================================================
# ANNOTATION INVENTORY
# ============================================================

def annotation_inventory(df):

    rows = []

    annotation_columns = [

        "target_class",
        "functional_description",
        "best_hit_description",
        "interpro_accessions",
        "interpro_descriptions",
        "pfam",
        "panther",
        "go_terms",
        "functional_keywords",

    ]


    for _, row in df.iterrows():

        gene_id = clean_text(
            row.get(
                "gene_id",
                ""
            )
        )

        rank = clean_text(
            row.get(
                "rank",
                ""
            )
        )


        for column in annotation_columns:

            value = clean_text(
                row.get(
                    column,
                    ""
                )
            )

            if not value:
                continue


            # ------------------------------------------------
            # Main descriptive annotations
            # ------------------------------------------------

            if column in {

                "target_class",
                "functional_description",
                "best_hit_description",

            }:

                values = [

                    normalize_spaces(
                        value
                    )

                ]

            else:

                values = split_multi(
                    value
                )


            for information in values:

                rows.append({

                    "gene_id":
                        gene_id,

                    "rank":
                        rank,

                    "information_type":
                        column,

                    "information":
                        information,

                })


    raw = pd.DataFrame(
        rows
    )


    # --------------------------------------------------------
    # Empty result
    # --------------------------------------------------------

    if raw.empty:

        empty_summary = pd.DataFrame(

            columns=[

                "information_type",
                "information",
                "candidate_count",
                "genes",

            ]

        )

        return raw, empty_summary


    # --------------------------------------------------------
    # Group and count annotations
    # --------------------------------------------------------

    summary = (

        raw

        .groupby(
            [

                "information_type",
                "information",

            ],
            dropna=False
        )

        .agg(

            candidate_count=(
                "gene_id",
                "nunique"
            ),

            genes=(
                "gene_id",
                lambda x:
                    "; ".join(
                        sorted(
                            set(
                                x
                            )
                        )
                    )
            ),

        )

        .reset_index()
    )


    summary = summary.sort_values(

        [

            "information_type",
            "candidate_count",
            "information",

        ],

        ascending=[

            True,
            False,
            True,

        ]
    )


    return raw, summary


# ============================================================
# READ TABLE
# ============================================================

def read_table(path):

    path = Path(
        path
    )


    if not path.exists():

        raise FileNotFoundError(
            f"File not found:\n{path}"
        )


    suffix = path.suffix.lower()


    if suffix == ".csv":

        return pd.read_csv(

            path,
            dtype=str,
            keep_default_na=False,

        )


    if suffix in {

        ".tsv",
        ".txt",

    }:

        return pd.read_csv(

            path,
            sep="\t",
            dtype=str,
            keep_default_na=False,

        )


    # --------------------------------------------------------
    # Fallback
    # --------------------------------------------------------

    return pd.read_csv(

        path,
        sep=None,
        engine="python",
        dtype=str,
        keep_default_na=False,

    )


# ============================================================
# VALIDATE COLUMNS
# ============================================================

def validate_columns(
    df,
    required_columns,
    file_label
):

    missing = [

        column

        for column
        in required_columns

        if column not in df.columns

    ]


    if missing:

        print(
            f"\nERROR: {file_label} is missing "
            f"required columns:"
        )

        for column in missing:

            print(
                f"  - {column}"
            )


        print(
            f"\nAvailable columns in {file_label}:"
        )

        for column in df.columns:

            print(
                f"  - {column}"
            )

        return False


    return True


# ============================================================
# CREATE GENE EVIDENCE MAP
# ============================================================

def create_gene_evidence_map(
    pesticide_df
):

    # --------------------------------------------------------
    # Work on a copy
    # --------------------------------------------------------

    working = pesticide_df.copy()


    working[
        "_merge_gene_id"
    ] = (

        working[
            "gene_id"
        ]

        .apply(
            clean_gene_id
        )

    )


    # --------------------------------------------------------
    # Combine duplicate gene evidence
    # --------------------------------------------------------

    gene_map = {}


    for gene_id, group in (

        working.groupby(
            "_merge_gene_id"
        )

    ):

        evidence_values = []

        for value in group[
            "pesticide_evidence"
        ]:

            value = clean_text(
                value
            )

            if (

                value

                and

                value
                !=
                "No verified pesticide/compound evidence mapped"

            ):

                if value not in evidence_values:

                    evidence_values.append(
                        value
                    )


        level_values = []

        for value in group[
            "pesticide_evidence_level"
        ]:

            value = clean_text(
                value
            )

            if (

                value
                and
                value != "NONE"

            ):

                for level in value.split(";"):

                    level = level.strip()

                    if (
                        level
                        and
                        level not in level_values
                    ):

                        level_values.append(
                            level
                        )


        scope_values = []

        for value in group[
            "pesticide_evidence_scope"
        ]:

            value = clean_text(
                value
            )

            if (

                value
                and
                value != "NONE"

            ):

                for scope in value.split(";"):

                    scope = scope.strip()

                    if (
                        scope
                        and
                        scope not in scope_values
                    ):

                        scope_values.append(
                            scope
                        )


        if evidence_values:

            gene_map[
                gene_id
            ] = {

                "pesticide_evidence":
                    " || ".join(
                        evidence_values
                    ),

                "pesticide_evidence_level":
                    "; ".join(
                        sorted(
                            level_values
                        )
                    ),

                "pesticide_evidence_scope":
                    "; ".join(
                        sorted(
                            scope_values
                        )
                    ),

            }


    return gene_map


# ============================================================
# ADD PESTICIDE EVIDENCE TO TOP 100 / TOP 10
# ============================================================

def add_pesticide_evidence_to_ranking(
    ranking_df,
    gene_evidence_map
):

    output = ranking_df.copy()


    # --------------------------------------------------------
    # Temporary normalized gene ID
    # --------------------------------------------------------

    output[
        "_merge_gene_id"
    ] = (

        output[
            "gene_id"
        ]

        .apply(
            clean_gene_id
        )

    )


    # --------------------------------------------------------
    # Evidence maps
    # --------------------------------------------------------

    output[
        "pesticide_evidence"
    ] = (

        output[
            "_merge_gene_id"
        ]

        .map(
            lambda gene_id:
                gene_evidence_map.get(
                    gene_id,
                    {}
                ).get(
                    "pesticide_evidence",
                    "No verified pesticide/compound evidence mapped"
                )
        )

    )


    output[
        "pesticide_evidence_level"
    ] = (

        output[
            "_merge_gene_id"
        ]

        .map(
            lambda gene_id:
                gene_evidence_map.get(
                    gene_id,
                    {}
                ).get(
                    "pesticide_evidence_level",
                    "NONE"
                )
        )

    )


    output[
        "pesticide_evidence_scope"
    ] = (

        output[
            "_merge_gene_id"
        ]

        .map(
            lambda gene_id:
                gene_evidence_map.get(
                    gene_id,
                    {}
                ).get(
                    "pesticide_evidence_scope",
                    "NONE"
                )
        )

    )


    # --------------------------------------------------------
    # Remove temporary column
    # --------------------------------------------------------

    output = output.drop(

        columns=[

            "_merge_gene_id"

        ]

    )


    # --------------------------------------------------------
    # Force exact final column order
    # --------------------------------------------------------

    final_columns = (

        EXPECTED_RANKING_COLUMNS

        +

        [

            "pesticide_evidence",
            "pesticide_evidence_level",
            "pesticide_evidence_scope",

        ]

    )


    output = output[
        final_columns
    ]


    return output


# ============================================================
# MAIN
# ============================================================

def main():

    print(
        "=" * 75
    )

    print(
        "METISA PLANA — PESTICIDE / COMPOUND EVIDENCE"
    )

    print(
        "TOP 500 + TOP 100 + TOP 10"
    )

    print(
        "=" * 75
    )


    # ========================================================
    # INPUT 1 — TOP 500
    # ========================================================

    print(
        "\nINPUT 1 OF 3 — TOP 500"
    )

    top500_path = input(
        "Mining candidates path:\n> "
    ).strip().strip('"')


    if not top500_path:

        print(
            "\nERROR: No Top 500 path supplied."
        )

        sys.exit(1)


    top500_path = Path(
        top500_path
    )


    if not top500_path.exists():

        print(
            f"\nERROR: Top 500 file not found:\n"
            f"{top500_path}"
        )

        sys.exit(1)


    # ========================================================
    # INPUT 2 — TOP 100
    # ========================================================

    print(
        "\nINPUT 2 OF 3 — TOP 100"
    )

    print(
        "This is entered separately because your Top 100 "
        "is not necessarily the same subset/file as Top 500."
    )

    top100_path = input(
        "Enter path to TOP 100 TSV:\n> "
    ).strip().strip('"')


    if not top100_path:

        print(
            "\nERROR: No Top 100 path supplied."
        )

        sys.exit(1)


    top100_path = Path(
        top100_path
    )


    if not top100_path.exists():

        print(
            f"\nERROR: Top 100 file not found:\n"
            f"{top100_path}"
        )

        sys.exit(1)


    # ========================================================
    # INPUT 3 — TOP 10
    # ========================================================

    print(
        "\nINPUT 3 OF 3 — TOP 10"
    )

    print(
        "This is entered separately because your Top 10 "
        "is not necessarily the same subset/file as Top 500."
    )

    top10_path = input(
        "Enter path to TOP 10 TSV:\n> "
    ).strip().strip('"')


    if not top10_path:

        print(
            "\nERROR: No Top 10 path supplied."
        )

        sys.exit(1)


    top10_path = Path(
        top10_path
    )


    if not top10_path.exists():

        print(
            f"\nERROR: Top 10 file not found:\n"
            f"{top10_path}"
        )

        sys.exit(1)


    # ========================================================
    # READ TOP 500
    # ========================================================

    print(
        "\n" + "-" * 75
    )

    print(
        "READING TOP 500"
    )

    print(
        "-" * 75
    )


    try:

        top500 = read_table(
            top500_path
        )

    except Exception as error:

        print(
            f"ERROR reading Top 500:\n{error}"
        )

        sys.exit(1)


    print(
        f"Rows: {len(top500):,}"
    )

    print(
        f"Columns: {len(top500.columns):,}"
    )


    # ========================================================
    # VALIDATE TOP 500
    # ========================================================

    if not validate_columns(

        top500,
        REQUIRED_TOP500_COLUMNS,
        "Top 500"

    ):

        sys.exit(1)


    # ========================================================
    # READ TOP 100
    # ========================================================

    print(
        "\n" + "-" * 75
    )

    print(
        "READING TOP 100"
    )

    print(
        "-" * 75
    )


    try:

        top100 = read_table(
            top100_path
        )

    except Exception as error:

        print(
            f"ERROR reading Top 100:\n{error}"
        )

        sys.exit(1)


    print(
        f"Rows: {len(top100):,}"
    )

    print(
        f"Columns: {len(top100.columns):,}"
    )


    # ========================================================
    # VALIDATE TOP 100
    # ========================================================

    if not validate_columns(

        top100,
        EXPECTED_RANKING_COLUMNS,
        "Top 100"

    ):

        sys.exit(1)


    # ========================================================
    # READ TOP 10
    # ========================================================

    print(
        "\n" + "-" * 75
    )

    print(
        "READING TOP 10"
    )

    print(
        "-" * 75
    )


    try:

        top10 = read_table(
            top10_path
        )

    except Exception as error:

        print(
            f"ERROR reading Top 10:\n{error}"
        )

        sys.exit(1)


    print(
        f"Rows: {len(top10):,}"
    )

    print(
        f"Columns: {len(top10.columns):,}"
    )


    # ========================================================
    # VALIDATE TOP 10
    # ========================================================

    if not validate_columns(

        top10,
        EXPECTED_RANKING_COLUMNS,
        "Top 10"

    ):

        sys.exit(1)


    # ========================================================
    # BUILD TOP 500 PESTICIDE EVIDENCE
    # ========================================================

    print(
        "\n" + "=" * 75
    )

    print(
        "STEP 1 — MAPPING VERIFIED PESTICIDE / COMPOUND EVIDENCE"
    )

    print(
        "=" * 75
    )


    top500_evidence = build_top500_evidence(
        top500
    )


    # ========================================================
    # BUILD PESTICIDE SUMMARY
    # ========================================================

    print(
        "\nBuilding pesticide evidence summary..."
    )


    pesticide_summary = build_pesticide_summary(
        top500
    )


    # ========================================================
    # BUILD TARGET DISCOVERY
    # ========================================================

    print(
        "Building target discovery tables..."
    )


    target_raw, target_summary = (
        annotation_inventory(
            top500
        )
    )


    # ========================================================
    # CREATE GENE → PESTICIDE EVIDENCE MAP
    # ========================================================

    print(
        "Creating gene-level pesticide evidence map..."
    )


    gene_evidence_map = (
        create_gene_evidence_map(
            top500_evidence
        )
    )


    print(
        f"Genes with mapped pesticide/compound evidence: "
        f"{len(gene_evidence_map):,}"
    )


    # ========================================================
    # ADD EVIDENCE TO TOP 100
    # ========================================================

    print(
        "\nAdding pesticide evidence to Top 100..."
    )


    top100_out = (
        add_pesticide_evidence_to_ranking(
            top100,
            gene_evidence_map
        )
    )


    # ========================================================
    # ADD EVIDENCE TO TOP 10
    # ========================================================

    print(
        "Adding pesticide evidence to Top 10..."
    )


    top10_out = (
        add_pesticide_evidence_to_ranking(
            top10,
            gene_evidence_map
        )
    )


    # ========================================================
    # MATCHING CHECK
    # ========================================================

    no_evidence_text = (
        "No verified pesticide/compound evidence mapped"
    )


    top100_matched = (

        top100_out[
            "pesticide_evidence"
        ]

        !=

        no_evidence_text

    ).sum()


    top10_matched = (

        top10_out[
            "pesticide_evidence"
        ]

        !=

        no_evidence_text

    ).sum()


    print(
        "\n" + "-" * 75
    )

    print(
        "MATCHING RESULTS"
    )

    print(
        "-" * 75
    )


    print(
        f"Top 100 evidence matches: "
        f"{top100_matched}/{len(top100_out)}"
    )


    print(
        f"Top 10 evidence matches: "
        f"{top10_matched}/{len(top10_out)}"
    )


    # ========================================================
    # OUTPUT DIRECTORY
    # ========================================================

    print(
        "\n" + "=" * 75
    )

    print(
        "OUTPUT DIRECTORY"
    )

    print(
        "=" * 75
    )


    output_path = input(

        "\nEnter output folder path "
        "(press Enter = same folder as Top 500):\n> "

    ).strip().strip('"')


    if not output_path:

        output_dir = (
            top500_path.parent
        )

    else:

        output_dir = Path(
            output_path
        )


    output_dir.mkdir(

        parents=True,
        exist_ok=True

    )


    print(
        f"\nOutput folder:\n{output_dir}"
    )


    # ========================================================
    # SAVE TOP 500 PESTICIDE EVIDENCE — TSV ONLY
    # ========================================================

    print(
        "\nSaving Top 500 pesticide evidence..."
    )


    top500_evidence_tsv = (

        output_dir /
        "top500_pesticide_evidence.tsv"

    )


    top500_evidence.to_csv(

        top500_evidence_tsv,
        sep="\t",
        index=False,
        encoding="utf-8-sig"

    )


    # ========================================================
    # SAVE TOP 500 PESTICIDE SUMMARY — TSV ONLY
    # ========================================================

    top500_summary_tsv = (

        output_dir /
        "top500_pesticide_evidence_summary.tsv"

    )


    pesticide_summary.to_csv(

        top500_summary_tsv,
        sep="\t",
        index=False,
        encoding="utf-8-sig"

    )


    # ========================================================
    # SAVE TARGET DISCOVERY — TSV ONLY
    # ========================================================

    top500_discovery_tsv = (

        output_dir /
        "top500_target_discovery.tsv"

    )


    target_summary.to_csv(

        top500_discovery_tsv,
        sep="\t",
        index=False,
        encoding="utf-8-sig"

    )


    # ========================================================
    # SAVE RAW TARGET DISCOVERY — TSV ONLY
    # ========================================================

    top500_raw_tsv = (

        output_dir /
        "top500_target_discovery_raw.tsv"

    )


    target_raw.to_csv(

        top500_raw_tsv,
        sep="\t",
        index=False,
        encoding="utf-8-sig"

    )


    # ========================================================
    # SAVE TOP 100 — TSV ONLY
    # ========================================================

    top100_tsv = (

        output_dir /
        "top100_with_pesticide.tsv"

    )


    top100_out.to_csv(

        top100_tsv,
        sep="\t",
        index=False,
        encoding="utf-8-sig"

    )


    # ========================================================
    # SAVE TOP 10 — CSV
    # ========================================================
    #
    # THIS IS THE MAIN FINAL SUBMISSION TABLE.
    #

    top10_csv = (

        output_dir /
        "top10_with_pesticide.csv"

    )


    top10_out.to_csv(

        top10_csv,
        index=False,
        encoding="utf-8-sig"

    )


    # ========================================================
    # SAVE TOP 10 — TSV
    # ========================================================

    top10_tsv = (

        output_dir /
        "top10_with_pesticide.tsv"

    )


    top10_out.to_csv(

        top10_tsv,
        sep="\t",
        index=False,
        encoding="utf-8-sig"

    )


    # ========================================================
    # FINAL OUTPUT CHECK
    # ========================================================

    print(
        "\n" + "=" * 75
    )

    print(
        "COMPLETED SUCCESSFULLY"
    )

    print(
        "=" * 75
    )


    print(
        "\nOUTPUT FILES:"
    )


    output_files = [

        top500_evidence_tsv,
        top500_summary_tsv,
        top500_discovery_tsv,
        top500_raw_tsv,
        top100_tsv,
        top10_csv,
        top10_tsv,

    ]


    for file in output_files:

        print(
            f"  ✓ {file.name}"
        )


    # ========================================================
    # FINAL COUNTS
    # ========================================================

    print(
        "\n" + "-" * 75
    )

    print(
        "FINAL COUNTS"
    )

    print(
        "-" * 75
    )


    print(
        f"Top 500: "
        f"{len(top500_evidence):,} candidates"
    )


    print(
        f"Top 100: "
        f"{len(top100_out):,} candidates"
    )


    print(
        f"Top 10: "
        f"{len(top10_out):,} candidates"
    )


    # ========================================================
    # TOP 10 PREVIEW
    # ========================================================

    print(
        "\n" + "=" * 75
    )

    print(
        "TOP 10 PESTICIDE EVIDENCE"
    )

    print(
        "=" * 75
    )


    preview_columns = [

        "rank",
        "gene_id",
        "gene_family",
        "target_class",
        "pesticide_evidence",
        "pesticide_evidence_level",
        "pesticide_evidence_scope",

    ]


    print(

        top10_out[
            preview_columns
        ]

        .to_string(
            index=False
        )

    )


    # ========================================================
    # FINAL NOTE
    # ========================================================

    print(
        "\n" + "=" * 75
    )

    print(
        "IMPORTANT"
    )

    print(
        "=" * 75
    )


    print(
        "Pesticide evidence is literature-backed and generally "
        "based on insect homologues or target families."
    )


    print(
        "It does NOT prove that the corresponding Metisa plana "
        "protein is itself an experimentally validated pesticide target."
    )


    print(
        "\nEvidence scope distinguishes direct M. plana validation "
        "from insect-homologue, target-family and resistance/synergist evidence."
    )


    print(
        "\nTop 100 and Top 10 were read independently and "
        "were NOT generated by slicing the Top 500."
    )


    print(
        "\nThe final Top 10 submission table is:"
    )


    print(
        f"  {top10_csv}"
    )


    print(
        "\nDone."
    )


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":

    main()



