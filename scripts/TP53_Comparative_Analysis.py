#!/usr/bin/env python3

"""
TP53 Comparative Analysis
=========================

Comparative sequence analysis of human and elephant TP53-related sequences.

Workflow
--------
1. Locate repository root
2. Load TP53_clean.fasta from data/processed/
3. Validate and clean protein sequences
4. Identify human TP53 reference
5. Calculate sequence lengths and amino-acid composition
6. Calculate pairwise global identity to human TP53
7. Map human TP53 hotspot positions using alignment-aware coordinates
8. Calculate hotspot conservation
9. Generate sequence-level feature table
10. Generate exploratory clustering
11. Save reproducible CSV outputs
12. Save hotspot mapping and summary tables

Important
---------
This is an exploratory comparative bioinformatics workflow.

It does NOT provide:
- clinical predictions
- cancer-risk predictions
- diagnostic predictions
- proof of elephant cancer resistance
- causal biological inference

The six predefined TP53 hotspots are:
R175, G245, R248, R249, R273 and R282.

Input
-----
data/processed/TP53_clean.fasta

Optional input:
data/processed/TP53_all_sequences.fasta

Outputs
-------
results/tp53_sequence_features.csv
results/tp53_hotspot_mapping.csv
results/tp53_comparative_features.csv
results/tp53_summary.json

figures/tp53_identity_to_human.png
figures/tp53_hotspot_conservation.png
figures/tp53_feature_clustering.png
"""

# ============================================================
# 1. STANDARD LIBRARY
# ============================================================

from pathlib import Path
import json
import sys
import warnings

warnings.filterwarnings("ignore")


# ============================================================
# 2. THIRD-PARTY LIBRARIES
# ============================================================

try:
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    from Bio import SeqIO
    from Bio.Align import PairwiseAligner

    from sklearn.preprocessing import StandardScaler
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score

except ImportError as exc:
    print("\nERROR: Missing required Python package.")
    print(f"Details: {exc}")
    print("\nInstall dependencies with:")
    print("pip install biopython numpy pandas matplotlib seaborn scikit-learn")
    sys.exit(1)


# ============================================================
# 3. REPOSITORY PATHS
# ============================================================

def find_repository_root():
    """
    Locate the repository root.

    The script is expected to live in:

        repository/scripts/TP53_Comparative_Analysis.py

    Therefore the repository root is normally the parent directory
    of the scripts directory.
    """

    script_path = Path(__file__).resolve()

    candidates = [
        script_path.parent.parent,
        Path.cwd().resolve(),
        Path.cwd().resolve().parent,
    ]

    for candidate in candidates:

        if (
            (candidate / "data").exists()
            and (
                (candidate / ".git").exists()
                or (candidate / "README.md").exists()
                or (candidate / "notebooks").exists()
                or (candidate / "scripts").exists()
            )
        ):
            return candidate

    # Fallback
    return script_path.parent.parent


ROOT = find_repository_root()

DATA_DIR = ROOT / "data"
PROCESSED_DIR = DATA_DIR / "processed"
DATABASE_DIR = DATA_DIR / "Database"

RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"

RESULTS_DIR.mkdir(
    parents=True,
    exist_ok=True
)

FIGURES_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 4. INPUT FILES
# ============================================================

CLEAN_FASTA = PROCESSED_DIR / "TP53_clean.fasta"
ALL_FASTA = PROCESSED_DIR / "TP53_all_sequences.fasta"

# Existing processed pair files are optional.
PAIR_FASTA = PROCESSED_DIR / "human_elephant_tp53_pair.fasta"


# ============================================================
# 5. TP53 HOTSPOTS
# ============================================================

HOTSPOTS = {
    175: "R175",
    245: "G245",
    248: "R248",
    249: "R249",
    273: "R273",
    282: "R282",
}


# ============================================================
# 6. VALID AMINO ACIDS
# ============================================================

VALID_AMINO_ACIDS = set(
    "ACDEFGHIKLMNPQRSTVWY"
)


# ============================================================
# 7. PRINT HEADER
# ============================================================

print("=" * 70)
print("TP53 COMPARATIVE ANALYSIS")
print("=" * 70)

print(f"\nRepository root:")
print(ROOT)

print("\nInput directory:")
print(PROCESSED_DIR)

print("\nOutput directory:")
print(RESULTS_DIR)


# ============================================================
# 8. CHECK INPUT
# ============================================================

if not CLEAN_FASTA.exists():

    print("\nERROR: Required input file was not found:")
    print(CLEAN_FASTA)

    print("\nExpected repository structure:")
    print("data/")
    print("└── processed/")
    print("    └── TP53_clean.fasta")

    print(
        "\nThe 41 MB elephant_proteome.fasta is NOT required "
        "for this analysis."
    )

    sys.exit(1)


# ============================================================
# 9. FASTA LOADING
# ============================================================

def load_fasta(path):
    """
    Load FASTA records and verify that records exist.
    """

    records = list(
        SeqIO.parse(
            path,
            "fasta"
        )
    )

    if not records:
        raise ValueError(
            f"No FASTA records found in: {path}"
        )

    return records


records = load_fasta(
    CLEAN_FASTA
)

print(
    f"\nLoaded sequences: {len(records)}"
)


# ============================================================
# 10. CLEAN SEQUENCES
# ============================================================

def clean_sequence(sequence):
    """
    Normalize a protein sequence.

    - uppercase
    - remove whitespace
    """

    sequence = str(sequence)
    sequence = sequence.upper()
    sequence = "".join(sequence.split())

    return sequence


cleaned_records = []

invalid_records = []

seen_ids = set()

seen_sequences = set()


for record in records:

    sequence = clean_sequence(
        record.seq
    )

    # --------------------------------------------------------
    # Validate sequence
    # --------------------------------------------------------

    invalid_characters = sorted(
        set(sequence) - VALID_AMINO_ACIDS
    )

    if invalid_characters:

        invalid_records.append(
            {
                "id": record.id,
                "invalid_characters": ",".join(
                    invalid_characters
                )
            }
        )

        continue

    # --------------------------------------------------------
    # Normalize ID
    # --------------------------------------------------------

    clean_id = (
        record.id
        .replace("|", "_")
        .replace(" ", "_")
    )

    # Keep original ID uniqueness
    if clean_id in seen_ids:
        continue

    # Avoid exact duplicate sequences
    if sequence in seen_sequences:
        continue

    seen_ids.add(clean_id)
    seen_sequences.add(sequence)

    cleaned_records.append(
        {
            "id": clean_id,
            "description": record.description,
            "sequence": sequence
        }
    )


print(
    f"Valid unique sequences retained: "
    f"{len(cleaned_records)}"
)

print(
    f"Sequences excluded: "
    f"{len(invalid_records)}"
)


# ============================================================
# 11. HUMAN TP53 REFERENCE DETECTION
# ============================================================

def identify_human_reference(sequence_records):
    """
    Identify the human TP53 reference.

    The previous workflow identified the UniProt P04637 record as:

        sp_P04637_P53_HUMAN

    Several fallback patterns are included.
    """

    preferred_patterns = [
        "P04637",
        "P53_HUMAN",
        "TP53_HUMAN",
        "HUMAN_TP53",
        "HUMAN"
    ]

    for pattern in preferred_patterns:

        for record in sequence_records:

            searchable = (
                record["id"] + " " +
                record["description"]
            ).upper()

            if pattern.upper() in searchable:

                return record

    return None


human_record = identify_human_reference(
    cleaned_records
)


if human_record is None:

    print(
        "\nERROR: Human TP53 reference could not "
        "be identified automatically."
    )

    print("\nAvailable sequence IDs:")

    for record in cleaned_records:
        print(" -", record["id"])

    print(
        "\nThe script will not guess a human reference."
    )

    sys.exit(1)


human_id = human_record["id"]
human_sequence = human_record["sequence"]


print("\nHuman TP53 reference:")
print("ID:", human_id)
print("Length:", len(human_sequence))


# ============================================================
# 12. VERIFY HUMAN TP53 LENGTH
# ============================================================

EXPECTED_HUMAN_LENGTH = 393

if len(human_sequence) != EXPECTED_HUMAN_LENGTH:

    print(
        "\nWARNING:"
        "\nThe detected human TP53 sequence has length "
        f"{len(human_sequence)}, not the expected "
        f"{EXPECTED_HUMAN_LENGTH} residues."
    )

    print(
        "Verify the human reference before interpreting "
        "hotspot coordinates."
    )

else:

    print(
        "\nHuman TP53 length check: PASS "
        f"({EXPECTED_HUMAN_LENGTH} aa)"
    )


# ============================================================
# 13. VERIFY HUMAN HOTSPOTS
# ============================================================

print("\nHuman TP53 hotspot reference residues:")

for position, label in HOTSPOTS.items():

    if position > len(human_sequence):

        print(
            f"{label}: ERROR — position outside sequence"
        )

        continue

    residue = human_sequence[
        position - 1
    ]

    expected_residue = label[0]

    status = (
        "PASS"
        if residue == expected_residue
        else "WARNING"
    )

    print(
        f"{label}: "
        f"reference={residue} "
        f"expected={expected_residue} "
        f"[{status}]"
    )


# ============================================================
# 14. SEQUENCE SUMMARY
# ============================================================

def amino_acid_composition(sequence):

    length = len(sequence)

    composition = {}

    for amino_acid in sorted(
        VALID_AMINO_ACIDS
    ):

        count = sequence.count(
            amino_acid
        )

        composition[
            f"frac_{amino_acid}"
        ] = (
            count / length
            if length > 0
            else 0.0
        )

    return composition


sequence_rows = []

for record in cleaned_records:

    sequence = record["sequence"]

    row = {
        "id": record["id"],
        "length": len(sequence),
        "is_human_reference": (
            record["id"] == human_id
        )
    }

    row.update(
        amino_acid_composition(
            sequence
        )
    )

    sequence_rows.append(row)


sequence_df = pd.DataFrame(
    sequence_rows
)


# ============================================================
# 15. PAIRWISE ALIGNER
# ============================================================

aligner = PairwiseAligner()

aligner.mode = "global"

# Identity-oriented scoring.
# This is used to calculate positional correspondence
# and exact residue identity.
aligner.match_score = 1.0
aligner.mismatch_score = 0.0

aligner.open_gap_score = -1.0
aligner.extend_gap_score = -0.1


# ============================================================
# 16. ALIGNMENT-AWARE IDENTITY
# ============================================================

def calculate_pairwise_identity(
    reference,
    query
):
    """
    Calculate exact amino-acid identity after global alignment.

    Identity =
        identical aligned residues /
        aligned non-gap residue pairs
    """

    alignment = aligner.align(
        reference,
        query
    )[0]

    aligned_reference = alignment[0]
    aligned_query = alignment[1]

    matches = 0
    comparable_positions = 0

    for ref_residue, query_residue in zip(
        aligned_reference,
        aligned_query
    ):

        if (
            ref_residue != "-"
            and query_residue != "-"
        ):

            comparable_positions += 1

            if (
                ref_residue
                == query_residue
            ):

                matches += 1

    if comparable_positions == 0:

        return np.nan

    return (
        matches
        / comparable_positions
    )


# ============================================================
# 17. MAP HUMAN REFERENCE POSITION
# ============================================================

def map_reference_position(
    reference,
    query,
    reference_position
):
    """
    Map a human-reference residue position onto
    a query sequence using a global alignment.

    The reference counter advances only when a
    non-gap residue occurs in the aligned reference.

    This prevents direct indexing errors caused
    by alignment gaps.
    """

    alignment = aligner.align(
        reference,
        query
    )[0]

    aligned_reference = alignment[0]
    aligned_query = alignment[1]

    reference_counter = 0

    for (
        ref_residue,
        query_residue
    ) in zip(
        aligned_reference,
        aligned_query
    ):

        if ref_residue != "-":

            reference_counter += 1

        if (
            reference_counter
            == reference_position
        ):

            return query_residue

    return None


# ============================================================
# 18. CALCULATE SEQUENCE IDENTITIES
# ============================================================

similarity_rows = []

print(
    "\nCalculating pairwise identity "
    "relative to human TP53..."
)


for record in cleaned_records:

    identity = calculate_pairwise_identity(
        human_sequence,
        record["sequence"]
    )

    similarity_rows.append(
        {
            "id": record["id"],
            "identity_to_human": identity,
            "identity_to_human_percent": (
                identity * 100
                if not np.isnan(identity)
                else np.nan
            )
        }
    )


similarity_df = pd.DataFrame(
    similarity_rows
)


# ============================================================
# 19. HOTSPOT MAPPING
# ============================================================

hotspot_rows = []

print(
    "\nMapping TP53 hotspots..."
)


for record in cleaned_records:

    sequence = record["sequence"]

    row = {
        "id": record["id"]
    }

    conserved_values = []

    for position, label in HOTSPOTS.items():

        query_residue = map_reference_position(
            human_sequence,
            sequence,
            position
        )

        human_residue = human_sequence[
            position - 1
        ]

        conserved = (
            query_residue is not None
            and query_residue
            == human_residue
        )

        row[
            f"{label}_human"
        ] = human_residue

        row[
            f"{label}_query"
        ] = (
            query_residue
            if query_residue is not None
            else "-"
        )

        row[
            f"{label}_conserved"
        ] = bool(conserved)

        conserved_values.append(
            int(conserved)
        )

    row[
        "hotspot_conservation_fraction"
    ] = (
        np.mean(
            conserved_values
        )
        if conserved_values
        else np.nan
    )

    row[
        "hotspot_conservation_percent"
    ] = (
        row[
            "hotspot_conservation_fraction"
        ] * 100
    )

    hotspot_rows.append(
        row
    )


hotspot_df = pd.DataFrame(
    hotspot_rows
)


# ============================================================
# 20. MERGE FEATURES
# ============================================================

feature_df = (
    sequence_df
    .merge(
        similarity_df,
        on="id",
        how="left"
    )
    .merge(
        hotspot_df,
        on="id",
        how="left"
    )
)


# ============================================================
# 21. ADD SEQUENCE CATEGORY
# ============================================================

def classify_sequence(record):

    text = (
        record["id"]
        + " "
        + record["description"]
    ).lower()

    if (
        "human" in text
        or "p53_human" in text
        or "p04637" in text
    ):

        return "Human"

    if (
        "retrogene" in text
        or "retro" in text
        or "rtg" in text
    ):

        return "Elephant_retrogene"

    if (
        "elephas" in text
        or "elephant" in text
    ):

        if (
            "african" in text
            or "loxodonta" in text
        ):

            return "African_elephant"

        if (
            "asian" in text
            or "maximus" in text
        ):

            return "Asian_elephant"

        return "Elephant"

    return "Other"


classification_lookup = {}

for record in cleaned_records:

    classification_lookup[
        record["id"]
    ] = classify_sequence(
        record
    )


feature_df["sequence_category"] = (
    feature_df["id"]
    .map(classification_lookup)
)


# ============================================================
# 22. SORT BY HUMAN IDENTITY
# ============================================================

feature_df = feature_df.sort_values(
    by="identity_to_human",
    ascending=False
).reset_index(
    drop=True
)


# ============================================================
# 23. SAVE BASIC SEQUENCE FEATURES
# ============================================================

sequence_features_output = (
    RESULTS_DIR
    / "tp53_sequence_features.csv"
)

sequence_df.to_csv(
    sequence_features_output,
    index=False
)


# ============================================================
# 24. SAVE HOTSPOT MAPPING
# ============================================================

hotspot_output = (
    RESULTS_DIR
    / "tp53_hotspot_mapping.csv"
)

hotspot_df.to_csv(
    hotspot_output,
    index=False
)


# ============================================================
# 25. SAVE COMPLETE COMPARATIVE FEATURES
# ============================================================

comparative_output = (
    RESULTS_DIR
    / "tp53_comparative_features.csv"
)

feature_df.to_csv(
    comparative_output,
    index=False
)


print("\nSaved tables:")
print(" -", sequence_features_output)
print(" -", hotspot_output)
print(" -", comparative_output)


# ============================================================
# 26. PRINT HOTSPOT SUMMARY
# ============================================================

print("\n" + "=" * 70)
print("HOTSPOT CONSERVATION SUMMARY")
print("=" * 70)

for position, label in HOTSPOTS.items():

    column = (
        f"{label}_conserved"
    )

    if column not in hotspot_df.columns:
        continue

    conservation_rate = (
        hotspot_df[column]
        .mean()
        * 100
    )

    print(
        f"{label}: "
        f"{conservation_rate:.2f}% "
        f"of analyzed sequences conserved"
    )


# ============================================================
# 27. PRINT SEQUENCE IDENTITY SUMMARY
# ============================================================

print("\n" + "=" * 70)
print("SEQUENCE IDENTITY SUMMARY")
print("=" * 70)

identity_summary = (
    feature_df[
        [
            "id",
            "sequence_category",
            "length",
            "identity_to_human_percent",
            "hotspot_conservation_percent"
        ]
    ]
)

print(
    identity_summary.to_string(
        index=False
    )
)


# ============================================================
# 28. EXPLORATORY CLUSTERING
# ============================================================

CLUSTER_FEATURES = [
    "length",
    "identity_to_human",
    "hotspot_conservation_fraction"
]


cluster_df = feature_df.copy()

available_cluster_features = [
    column
    for column in CLUSTER_FEATURES
    if column in cluster_df.columns
]


cluster_data = (
    cluster_df[
        available_cluster_features
    ]
    .replace(
        [np.inf, -np.inf],
        np.nan
    )
)


# Remove rows with incomplete
# clustering features.

valid_cluster_rows = (
    cluster_data
    .dropna()
)


if len(valid_cluster_rows) >= 3:

    scaler = StandardScaler()

    X_scaled = scaler.fit_transform(
        valid_cluster_rows
    )

    # --------------------------------------------------------
    # Determine a conservative cluster count.
    # --------------------------------------------------------

    max_k = min(
        4,
        len(valid_cluster_rows) - 1
    )

    best_k = None
    best_score = -np.inf

    for k in range(
        2,
        max_k + 1
    ):

        model = KMeans(
            n_clusters=k,
            random_state=42,
            n_init=20
        )

        labels = model.fit_predict(
            X_scaled
        )

        # Silhouette requires at least
        # two unique clusters.

        if len(
            np.unique(labels)
        ) < 2:

            continue

        score = silhouette_score(
            X_scaled,
            labels
        )

        if score > best_score:

            best_score = score
            best_k = k

    if best_k is not None:

        final_model = KMeans(
            n_clusters=best_k,
            random_state=42,
            n_init=20
        )

        final_labels = (
            final_model
            .fit_predict(X_scaled)
        )

        cluster_df.loc[
            valid_cluster_rows.index,
            "cluster"
        ] = final_labels

        cluster_df[
            "cluster"
        ] = cluster_df[
            "cluster"
        ].astype(
            "Int64"
        )

        cluster_df[
            "silhouette_score"
        ] = np.nan

        cluster_df.loc[
            valid_cluster_rows.index,
            "silhouette_score"
        ] = best_score

        print(
            "\nExploratory clustering:"
        )

        print(
            f"Selected clusters: {best_k}"
        )

        print(
            f"Silhouette score: "
            f"{best_score:.4f}"
        )

    else:

        print(
            "\nClustering skipped: "
            "no valid cluster solution."
        )

else:

    cluster_df["cluster"] = pd.NA

    cluster_df[
        "silhouette_score"
    ] = np.nan

    print(
        "\nClustering skipped: "
        "insufficient valid observations."
    )


# ============================================================
# 29. SAVE CLUSTERED DATASET
# ============================================================

cluster_output = (
    RESULTS_DIR
    / "tp53_comparative_features_clustered.csv"
)

cluster_df.to_csv(
    cluster_output,
    index=False
)

print(
    "\nSaved clustered dataset:"
)

print(
    cluster_output
)


# ============================================================
# 30. FIGURE 1 — HUMAN IDENTITY
# ============================================================

plot_df = feature_df.copy()

plot_df = plot_df.sort_values(
    "identity_to_human_percent",
    ascending=True
)


plt.figure(
    figsize=(
        10,
        max(
            5,
            len(plot_df) * 0.35
        )
    )
)

plt.barh(
    plot_df["id"],
    plot_df[
        "identity_to_human_percent"
    ]
)

plt.xlabel(
    "Identity to human TP53 (%)"
)

plt.ylabel(
    "Sequence"
)

plt.title(
    "Sequence Identity Relative to Human TP53"
)

plt.tight_layout()

identity_figure = (
    FIGURES_DIR
    / "tp53_identity_to_human.png"
)

plt.savefig(
    identity_figure,
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 31. FIGURE 2 — HOTSPOT CONSERVATION
# ============================================================

conservation_columns = [
    f"{label}_conserved"
    for label in HOTSPOTS.values()
]


heatmap_data = (
    hotspot_df[
        [
            "id"
        ]
        + conservation_columns
    ]
    .set_index("id")
    .astype(int)
)


heatmap_data.columns = [
    label
    for label in HOTSPOTS.values()
]


plt.figure(
    figsize=(
        10,
        max(
            5,
            len(heatmap_data) * 0.35
        )
    )
)

sns.heatmap(
    heatmap_data,
    vmin=0,
    vmax=1,
    annot=True,
    fmt="d",
    linewidths=0.5,
    cbar_kws={
        "label": "Conserved"
    }
)

plt.xlabel(
    "Human TP53 hotspot"
)

plt.ylabel(
    "Sequence"
)

plt.title(
    "TP53 Hotspot Conservation"
)

plt.tight_layout()

hotspot_figure = (
    FIGURES_DIR
    / "tp53_hotspot_conservation.png"
)

plt.savefig(
    hotspot_figure,
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 32. FIGURE 3 — EXPLORATORY CLUSTERING
# ============================================================

if (
    "cluster" in cluster_df.columns
    and
    cluster_df["cluster"].notna().sum() >= 2
):

    plt.figure(
        figsize=(9, 6)
    )

    sns.scatterplot(
        data=cluster_df,
        x="identity_to_human_percent",
        y="hotspot_conservation_percent",
        hue="cluster",
        style="sequence_category",
        s=100
    )

    plt.xlabel(
        "Identity to human TP53 (%)"
    )

    plt.ylabel(
        "Hotspot conservation (%)"
    )

    plt.title(
        "Exploratory TP53 Sequence Clustering"
    )

    plt.tight_layout()

    cluster_figure = (
        FIGURES_DIR
        / "tp53_feature_clustering.png"
    )

    plt.savefig(
        cluster_figure,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

else:

    print(
        "\nClustering figure skipped."
    )


# ============================================================
# 33. SUMMARY JSON
# ============================================================

summary = {

    "project": (
        "Elephant TP53 Hotspot Mapping"
    ),

    "analysis": (
        "TP53 Comparative Sequence Analysis"
    ),

    "input_file": str(
        CLEAN_FASTA.relative_to(ROOT)
    ),

    "number_of_input_records": len(
        records
    ),

    "number_of_valid_unique_sequences": len(
        cleaned_records
    ),

    "number_of_excluded_sequences": len(
        invalid_records
    ),

    "human_reference": human_id,

    "human_reference_length": len(
        human_sequence
    ),

    "expected_human_reference_length": (
        EXPECTED_HUMAN_LENGTH
    ),

    "hotspots": {
        str(position): label
        for position, label
        in HOTSPOTS.items()
    },

    "hotspot_conservation": {},

    "outputs": {
        "sequence_features": str(
            sequence_features_output.relative_to(
                ROOT
            )
        ),
        "hotspot_mapping": str(
            hotspot_output.relative_to(
                ROOT
            )
        ),
        "comparative_features": str(
            comparative_output.relative_to(
                ROOT
            )
        ),
        "clustered_features": str(
            cluster_output.relative_to(
                ROOT
            )
        )
    },

    "interpretation": (
        "Exploratory comparative sequence analysis; "
        "not a clinical or validated predictive model."
    )
}


for position, label in HOTSPOTS.items():

    column = (
        f"{label}_conserved"
    )

    if column in hotspot_df.columns:

        summary[
            "hotspot_conservation"
        ][label] = (
            float(
                hotspot_df[column].mean()
            )
        )


summary_output = (
    RESULTS_DIR
    / "tp53_summary.json"
)


with open(
    summary_output,
    "w",
    encoding="utf-8"
) as handle:

    json.dump(
        summary,
        handle,
        indent=2
    )


# ============================================================
# 34. FINAL REPORT
# ============================================================

print("\n" + "=" * 70)
print("ANALYSIS COMPLETED")
print("=" * 70)

print("\nInput:")
print(
    CLEAN_FASTA.relative_to(ROOT)
)

print("\nReference:")
print(human_id)

print(
    f"Human TP53 length: "
    f"{len(human_sequence)} aa"
)

print(
    f"Sequences analyzed: "
    f"{len(cleaned_records)}"
)

print("\nResults:")
print(
    sequence_features_output.relative_to(ROOT)
)

print(
    hotspot_output.relative_to(ROOT)
)

print(
    comparative_output.relative_to(ROOT)
)

print(
    cluster_output.relative_to(ROOT)
)

print(
    summary_output.relative_to(ROOT)
)

print("\nFigures:")
print(
    identity_figure.relative_to(ROOT)
)

print(
    hotspot_figure.relative_to(ROOT)
)

if (
    "cluster_figure" in locals()
):
    print(
        cluster_figure.relative_to(ROOT)
    )

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
