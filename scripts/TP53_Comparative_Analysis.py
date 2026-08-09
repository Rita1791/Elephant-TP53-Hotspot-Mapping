from pathlib import Path

code = r'''#!/usr/bin/env python3
"""
TP53 Comparative Analysis
=========================

Reproducible comparative sequence analysis of human and elephant
TP53-related protein sequences.

This workflow:
1. Locates the repository root.
2. Loads data/processed/TP53_clean.fasta.
3. Validates and normalizes protein sequences.
4. Identifies the human TP53 reference (UniProt P04637).
5. Calculates sequence length and amino-acid composition.
6. Calculates global pairwise amino-acid identity to human TP53.
7. Maps predefined human TP53 hotspot coordinates through global alignments.
8. Separates exact residue identity from broader biological interpretation.
9. Calculates comparative hotspot identity statistics excluding the human reference.
10. Classifies sequences using accession/description metadata.
11. Performs optional exploratory K-means clustering on sequence-level features.
12. Generates publication-ready PNG figures.
13. Writes CSV and JSON outputs documenting the analysis.

Important scientific scope
--------------------------
This is an exploratory comparative bioinformatics workflow.

It does NOT establish:
- clinical utility,
- diagnostic or prognostic value,
- individual cancer risk,
- causal mechanisms of elephant cancer resistance,
- proof of Peto's paradox,
- functional equivalence of elephant TP53-related sequences,
- or experimental validation.

"Hotspot conservation" in this repository means EXACT AMINO-ACID
IDENTITY at an alignment-mapped human TP53 hotspot position. It is
not a substitute for a formal evolutionary conservation score.

Human TP53 hotspots used:
R175, G245, R248, R249, R273, R282.

Primary input
-------------
data/processed/TP53_clean.fasta

Primary outputs
--------------
results/tp53_sequence_features.csv
results/tp53_hotspot_mapping.csv
results/tp53_hotspot_identity_summary.csv
results/tp53_comparative_features.csv
results/tp53_comparative_features_clustered.csv
results/tp53_excluded_sequences.csv
results/tp53_summary.json

Figures
-------
figures/tp53_identity_to_human.png
figures/tp53_hotspot_conservation.png
figures/tp53_feature_clustering.png

Run
---
python scripts/TP53_Comparative_Analysis.py
"""

from __future__ import annotations

from pathlib import Path
import json
import re
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler


# ============================================================
# 1. CONFIGURATION
# ============================================================

PROJECT_NAME = "Elephant TP53 Hotspot Mapping"
ANALYSIS_NAME = "TP53 Comparative Sequence Analysis"

# Human canonical TP53 protein: UniProt P04637, 393 aa.
HUMAN_REFERENCE_ACCESSION = "P04637"
EXPECTED_HUMAN_LENGTH = 393

# Human TP53 cancer hotspot positions used by this project.
HOTSPOTS = {
    175: "R175",
    245: "G245",
    248: "R248",
    249: "R249",
    273: "R273",
    282: "R282",
}

VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

# Alignment parameters are explicit so the workflow is reproducible.
ALIGNMENT_MATCH_SCORE = 1.0
ALIGNMENT_MISMATCH_SCORE = -1.0
ALIGNMENT_GAP_OPEN_SCORE = -2.0
ALIGNMENT_GAP_EXTEND_SCORE = -0.5

# Exploratory clustering settings.
CLUSTER_RANDOM_STATE = 42
CLUSTER_MAX_K = 4
MIN_CLUSTER_OBSERVATIONS = 4


# ============================================================
# 2. REPOSITORY PATHS
# ============================================================

def find_repository_root() -> Path:
    """Locate the repository root from the script location or CWD."""
    script_path = Path(__file__).resolve()

    candidates = [
        script_path.parent.parent,
        Path.cwd().resolve(),
        Path.cwd().resolve().parent,
    ]

    for candidate in candidates:
        if (
            (candidate / "data").is_dir()
            and (
                (candidate / ".git").exists()
                or (candidate / "README.md").exists()
                or (candidate / "scripts").is_dir()
            )
        ):
            return candidate

    raise RuntimeError(
        "Could not identify the repository root. "
        "Run this script from the repository or place it in scripts/."
    )


ROOT = find_repository_root()

DATA_DIR = ROOT / "data"
PROCESSED_DIR = DATA_DIR / "processed"

RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

CLEAN_FASTA = PROCESSED_DIR / "TP53_clean.fasta"


# ============================================================
# 3. OUTPUT PATHS
# ============================================================

SEQUENCE_FEATURES_OUTPUT = RESULTS_DIR / "tp53_sequence_features.csv"
HOTSPOT_OUTPUT = RESULTS_DIR / "tp53_hotspot_mapping.csv"
HOTSPOT_SUMMARY_OUTPUT = RESULTS_DIR / "tp53_hotspot_identity_summary.csv"
COMPARATIVE_OUTPUT = RESULTS_DIR / "tp53_comparative_features.csv"
CLUSTER_OUTPUT = RESULTS_DIR / "tp53_comparative_features_clustered.csv"
EXCLUDED_OUTPUT = RESULTS_DIR / "tp53_excluded_sequences.csv"
SUMMARY_OUTPUT = RESULTS_DIR / "tp53_summary.json"

IDENTITY_FIGURE = FIGURES_DIR / "tp53_identity_to_human.png"
HOTSPOT_FIGURE = FIGURES_DIR / "tp53_hotspot_conservation.png"
CLUSTER_FIGURE = FIGURES_DIR / "tp53_feature_clustering.png"


# ============================================================
# 4. UTILITY FUNCTIONS
# ============================================================

def print_header(title: str) -> None:
    """Print a consistent terminal section header."""
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def relative_path(path: Path) -> str:
    """Return a repository-relative path when possible."""
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def clean_sequence(sequence: str) -> str:
    """Normalize a protein sequence to uppercase without whitespace."""
    return "".join(str(sequence).upper().split())


def extract_accession(record_id: str) -> str | None:
    """
    Extract a common accession-like token from a FASTA identifier.

    Examples:
        sp_P04637_P53_HUMAN -> P04637
        XP_049714738.1       -> XP_049714738.1
    """
    match = re.search(
        r"(?:sp_)?([A-Z]{1,4}_?\d+(?:\.\d+)?)",
        record_id.upper(),
    )
    return match.group(1) if match else None


def amino_acid_composition(sequence: str) -> dict[str, float]:
    """Calculate amino-acid fractions for the 20 standard residues."""
    length = len(sequence)
    return {
        f"frac_{aa}": (
            sequence.count(aa) / length if length else np.nan
        )
        for aa in sorted(VALID_AMINO_ACIDS)
    }


def classify_sequence(record: dict) -> str:
    """
    Assign a broad sequence category from accession/description metadata.

    Classification is metadata-based and should not be interpreted as
    an independently verified taxonomic annotation.
    """
    text = (
        f'{record["id"]} {record["description"]}'
    ).lower()

    if (
        "p04637" in text
        or "p53_human" in text
        or "tp53_human" in text
        or "human tp53" in text
    ):
        return "Human"

    if (
        "retrogene" in text
        or "retro" in text
        or "rtg" in text
    ):
        return "Elephant_retrogene"

    if (
        "loxodonta" in text
        or "african elephant" in text
    ):
        return "African_elephant"

    if (
        "elephas" in text
        or "asian elephant" in text
    ):
        return "Asian_elephant"

    if "elephant" in text:
        return "Elephant"

    return "Other"


# ============================================================
# 5. FASTA VALIDATION AND LOADING
# ============================================================

def load_and_validate_fasta(path: Path) -> tuple[list[dict], list[dict], dict]:
    """
    Load FASTA records.

    Records with invalid amino-acid symbols are excluded and documented.
    Duplicate sequence IDs are excluded because they are ambiguous.
    Exact duplicate sequences with different IDs are retained because
    they may represent distinct biological records/accessions.
    """
    if not path.exists():
        raise FileNotFoundError(
            f"Required input file was not found: {path}"
        )

    raw_records = list(SeqIO.parse(path, "fasta"))

    if not raw_records:
        raise ValueError(f"No FASTA records found in: {path}")

    cleaned_records: list[dict] = []
    excluded_records: list[dict] = []
    seen_ids: set[str] = set()

    for record in raw_records:
        sequence = clean_sequence(record.seq)

        invalid = sorted(set(sequence) - VALID_AMINO_ACIDS)

        clean_id = (
            record.id
            .replace("|", "_")
            .replace(" ", "_")
        )

        if clean_id in seen_ids:
            excluded_records.append(
                {
                    "id": clean_id,
                    "reason": "duplicate_sequence_id",
                    "details": "A record with the same normalized ID was already retained.",
                }
            )
            continue

        seen_ids.add(clean_id)

        if invalid:
            excluded_records.append(
                {
                    "id": clean_id,
                    "reason": "invalid_amino_acid_symbols",
                    "details": ",".join(invalid),
                }
            )
            continue

        if not sequence:
            excluded_records.append(
                {
                    "id": clean_id,
                    "reason": "empty_sequence",
                    "details": "Sequence contained no residues.",
                }
            )
            continue

        cleaned_records.append(
            {
                "id": clean_id,
                "description": record.description,
                "sequence": sequence,
            }
        )

    metadata = {
        "input_records": len(raw_records),
        "retained_records": len(cleaned_records),
        "excluded_records": len(excluded_records),
    }

    return cleaned_records, excluded_records, metadata


# ============================================================
# 6. HUMAN REFERENCE IDENTIFICATION
# ============================================================

def identify_human_reference(records: list[dict]) -> dict:
    """
    Identify the canonical human TP53 reference.

    P04637 is required as the preferred identifier. A conservative
    fallback checks explicit human TP53 naming rather than the broad
    word 'human' alone.
    """
    preferred_patterns = [
        "P04637",
        "sp_P04637_P53_HUMAN",
        "P53_HUMAN",
        "TP53_HUMAN",
    ]

    for pattern in preferred_patterns:
        matches = [
            record
            for record in records
            if pattern.lower() in (
                f'{record["id"]} {record["description"]}'
            ).lower()
        ]

        if matches:
            return matches[0]

    raise RuntimeError(
        "Canonical human TP53 reference P04637 could not be identified. "
        "The script will not silently guess another human sequence."
    )


# ============================================================
# 7. HUMAN REFERENCE VALIDATION
# ============================================================

def validate_human_reference(human_record: dict) -> dict:
    """Validate expected human TP53 length and hotspot residues."""
    sequence = human_record["sequence"]

    validation = {
        "accession": HUMAN_REFERENCE_ACCESSION,
        "sequence_id": human_record["id"],
        "length": len(sequence),
        "expected_length": EXPECTED_HUMAN_LENGTH,
        "length_pass": len(sequence) == EXPECTED_HUMAN_LENGTH,
        "hotspots": {},
    }

    for position, label in HOTSPOTS.items():
        if position > len(sequence):
            validation["hotspots"][label] = {
                "position": position,
                "observed": None,
                "expected": label[0],
                "pass": False,
                "status": "outside_sequence",
            }
            continue

        observed = sequence[position - 1]
        expected = label[0]

        validation["hotspots"][label] = {
            "position": position,
            "observed": observed,
            "expected": expected,
            "pass": observed == expected,
            "status": "PASS" if observed == expected else "WARNING",
        }

    return validation


# ============================================================
# 8. GLOBAL ALIGNMENT
# ============================================================

def create_aligner() -> PairwiseAligner:
    """Create the explicit global alignment configuration."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = ALIGNMENT_MATCH_SCORE
    aligner.mismatch_score = ALIGNMENT_MISMATCH_SCORE
    aligner.open_gap_score = ALIGNMENT_GAP_OPEN_SCORE
    aligner.extend_gap_score = ALIGNMENT_GAP_EXTEND_SCORE
    return aligner


ALIGNER = create_aligner()


def get_best_alignment(reference: str, query: str):
    """Return the highest-scoring global alignment."""
    return ALIGNER.align(reference, query)[0]


def calculate_pairwise_identity(
    reference: str,
    query: str,
) -> tuple[float, int, int]:
    """
    Calculate exact amino-acid identity after global alignment.

    Denominator:
        aligned positions where neither sequence contains a gap.

    Returns:
        identity_fraction, identical_positions, comparable_positions
    """
    alignment = get_best_alignment(reference, query)

    aligned_reference = alignment[0]
    aligned_query = alignment[1]

    matches = 0
    comparable = 0

    for ref_residue, query_residue in zip(
        aligned_reference,
        aligned_query,
    ):
        if ref_residue == "-" or query_residue == "-":
            continue

        comparable += 1

        if ref_residue == query_residue:
            matches += 1

    if comparable == 0:
        return np.nan, 0, 0

    return matches / comparable, matches, comparable


def map_reference_position(
    reference: str,
    query: str,
    reference_position: int,
) -> dict:
    """
    Map a 1-based human TP53 residue position through a global alignment.

    Returns a structured result distinguishing:
        mapped residue,
        alignment gap,
        and unmapped coordinate.
    """
    alignment = get_best_alignment(reference, query)

    aligned_reference = alignment[0]
    aligned_query = alignment[1]

    reference_counter = 0

    for ref_residue, query_residue in zip(
        aligned_reference,
        aligned_query,
    ):
        if ref_residue != "-":
            reference_counter += 1

            if reference_counter == reference_position:
                if query_residue == "-":
                    return {
                        "status": "gap",
                        "residue": "-",
                    }

                return {
                    "status": "mapped",
                    "residue": query_residue,
                }

    return {
        "status": "unmapped",
        "residue": None,
    }


# ============================================================
# 9. SEQUENCE-LEVEL FEATURES
# ============================================================

def build_sequence_features(
    records: list[dict],
    human_id: str,
) -> pd.DataFrame:
    """Build sequence length and amino-acid composition features."""
    rows = []

    for record in records:
        sequence = record["sequence"]

        row = {
            "id": record["id"],
            "accession": extract_accession(record["id"]),
            "length": len(sequence),
            "is_human_reference": record["id"] == human_id,
            "sequence_category": classify_sequence(record),
        }

        row.update(amino_acid_composition(sequence))
        rows.append(row)

    return pd.DataFrame(rows)


# ============================================================
# 10. PAIRWISE IDENTITY FEATURES
# ============================================================

def build_identity_features(
    records: list[dict],
    human_sequence: str,
) -> pd.DataFrame:
    """Calculate global identity of every sequence to human TP53."""
    rows = []

    for record in records:
        identity, matches, comparable = calculate_pairwise_identity(
            human_sequence,
            record["sequence"],
        )

        rows.append(
            {
                "id": record["id"],
                "identity_to_human": identity,
                "identity_to_human_percent": (
                    identity * 100
                    if not np.isnan(identity)
                    else np.nan
                ),
                "identical_aligned_residues": matches,
                "comparable_aligned_residues": comparable,
            }
        )

    return pd.DataFrame(rows)


# ============================================================
# 11. HOTSPOT MAPPING
# ============================================================

def build_hotspot_mapping(
    records: list[dict],
    human_id: str,
    human_sequence: str,
) -> pd.DataFrame:
    """
    Map all predefined human TP53 hotspots onto every sequence.

    Exact residue identity is recorded independently from mapping status.
    """
    rows = []

    for record in records:
        sequence = record["sequence"]

        row = {
            "id": record["id"],
            "is_human_reference": record["id"] == human_id,
        }

        conserved_flags = []
        evaluable_flags = []

        for position, label in HOTSPOTS.items():
            human_residue = human_sequence[position - 1]

            mapped = map_reference_position(
                human_sequence,
                sequence,
                position,
            )

            query_residue = mapped["residue"]
            status = mapped["status"]

            is_conserved = (
                status == "mapped"
                and query_residue == human_residue
            )

            is_evaluable = status == "mapped"

            row[f"{label}_human"] = human_residue
            row[f"{label}_query"] = (
                query_residue if query_residue is not None else "-"
            )
            row[f"{label}_status"] = status
            row[f"{label}_conserved"] = bool(is_conserved)

            conserved_flags.append(int(is_conserved))
            evaluable_flags.append(int(is_evaluable))

        evaluable_count = int(sum(evaluable_flags))
        conserved_count = int(sum(conserved_flags))

        row["hotspots_total"] = len(HOTSPOTS)
        row["hotspots_evaluable"] = evaluable_count
        row["hotspots_conserved"] = conserved_count
        row["hotspot_conservation_fraction"] = (
            conserved_count / evaluable_count
            if evaluable_count
            else np.nan
        )
        row["hotspot_conservation_percent"] = (
            row["hotspot_conservation_fraction"] * 100
            if not np.isnan(row["hotspot_conservation_fraction"])
            else np.nan
        )

        rows.append(row)

    return pd.DataFrame(rows)


# ============================================================
# 12. COMPARATIVE HOTSPOT SUMMARY
# ============================================================

def build_hotspot_identity_summary(
    hotspot_df: pd.DataFrame,
    human_id: str,
    human_sequence: str,
) -> pd.DataFrame:
    """
    Summarize hotspot identity across comparative sequences.

    The human reference is explicitly excluded from the denominator.
    Gaps and unmapped positions are not counted as exact residue identity.
    """
    comparative_df = hotspot_df[
        hotspot_df["id"] != human_id
    ].copy()

    rows = []

    for position, label in HOTSPOTS.items():
        residue_column = f"{label}_query"
        status_column = f"{label}_status"
        conserved_column = f"{label}_conserved"

        mapped_mask = comparative_df[status_column] == "mapped"
        gap_mask = comparative_df[status_column] == "gap"
        unmapped_mask = comparative_df[status_column] == "unmapped"

        mapped_count = int(mapped_mask.sum())
        conserved_count = int(
            comparative_df.loc[
                mapped_mask,
                conserved_column,
            ].sum()
        )

        substituted_count = mapped_count - conserved_count

        gap_count = int(gap_mask.sum())
        unmapped_count = int(unmapped_mask.sum())

        rows.append(
            {
                "hotspot": label,
                "human_position": position,
                "human_residue": human_sequence[position - 1],
                "comparative_sequences": len(comparative_df),
                "mapped": mapped_count,
                "conserved": conserved_count,
                "substituted": substituted_count,
                "gaps": gap_count,
                "unmapped": unmapped_count,
                "exact_residue_identity_percent": (
                    round(
                        conserved_count / mapped_count * 100,
                        2,
                    )
                    if mapped_count
                    else np.nan
                ),
            }
        )

    return pd.DataFrame(rows)


# ============================================================
# 13. MERGE COMPLETE FEATURE TABLE
# ============================================================

def build_complete_feature_table(
    sequence_df: pd.DataFrame,
    identity_df: pd.DataFrame,
    hotspot_df: pd.DataFrame,
) -> pd.DataFrame:
    """Merge all sequence-level, identity, and hotspot features."""
    feature_df = (
        sequence_df
        .merge(identity_df, on="id", how="left")
        .merge(hotspot_df, on=["id"], how="left")
    )

    return feature_df.sort_values(
        by="identity_to_human",
        ascending=False,
        na_position="last",
    ).reset_index(drop=True)


# ============================================================
# 14. EXPLORATORY CLUSTERING
# ============================================================

def exploratory_clustering(
    feature_df: pd.DataFrame,
) -> tuple[pd.DataFrame, dict]:
    """
    Perform exploratory K-means clustering on comparative sequences.

    The human reference is excluded from clustering because it is the
    reference anchor and would otherwise be treated as another biological
    observation.

    Features:
        length
        identity_to_human
        hotspot_conservation_fraction

    Clustering is descriptive/exploratory and is not a predictive model.
    """
    result = feature_df.copy()

    result["cluster"] = pd.Series(
        pd.NA,
        index=result.index,
        dtype="Int64",
    )
    result["silhouette_score"] = np.nan

    comparative_mask = ~result["is_human_reference"].astype(bool)

    cluster_features = [
        "length",
        "identity_to_human",
        "hotspot_conservation_fraction",
    ]

    valid = (
        result.loc[
            comparative_mask,
            cluster_features,
        ]
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
    )

    metadata = {
        "performed": False,
        "features": cluster_features,
        "observations": int(len(valid)),
        "selected_k": None,
        "silhouette_score": None,
    }

    if len(valid) < MIN_CLUSTER_OBSERVATIONS:
        return result, metadata

    max_k = min(CLUSTER_MAX_K, len(valid) - 1)

    if max_k < 2:
        return result, metadata

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(valid)

    best_k = None
    best_score = -np.inf

    for k in range(2, max_k + 1):
        model = KMeans(
            n_clusters=k,
            random_state=CLUSTER_RANDOM_STATE,
            n_init=20,
        )

        labels = model.fit_predict(X_scaled)

        if len(np.unique(labels)) < 2:
            continue

        score = silhouette_score(X_scaled, labels)

        if score > best_score:
            best_score = score
            best_k = k

    if best_k is None:
        return result, metadata

    final_model = KMeans(
        n_clusters=best_k,
        random_state=CLUSTER_RANDOM_STATE,
        n_init=20,
    )

    final_labels = final_model.fit_predict(X_scaled)

    result.loc[
        valid.index,
        "cluster",
    ] = final_labels

    result.loc[
        valid.index,
        "silhouette_score",
    ] = best_score

    result["cluster"] = result["cluster"].astype("Int64")

    metadata.update(
        {
            "performed": True,
            "selected_k": int(best_k),
            "silhouette_score": float(best_score),
        }
    )

    return result, metadata


# ============================================================
# 15. FIGURES
# ============================================================

def configure_plot_style() -> None:
    """Apply a clean, consistent figure style."""
    sns.set_theme(
        style="whitegrid",
        context="notebook",
    )


def plot_identity(feature_df: pd.DataFrame) -> None:
    """Plot global sequence identity relative to human TP53."""
    plot_df = feature_df.sort_values(
        "identity_to_human_percent",
        ascending=True,
    )

    height = max(5, len(plot_df) * 0.35)

    plt.figure(figsize=(10, height))

    sns.barplot(
        data=plot_df,
        x="identity_to_human_percent",
        y="id",
        hue="sequence_category",
        dodge=False,
        legend=True,
    )

    plt.xlabel("Global amino-acid identity to human TP53 (%)")
    plt.ylabel("Sequence")
    plt.title(
        "TP53 Sequence Identity Relative to Human Reference"
    )
    plt.xlim(0, 100)
    plt.tight_layout()

    plt.savefig(
        IDENTITY_FIGURE,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


def plot_hotspot_conservation(
    hotspot_df: pd.DataFrame,
    human_id: str,
) -> None:
    """Plot exact residue identity at mapped TP53 hotspots."""
    comparative_df = hotspot_df[
        hotspot_df["id"] != human_id
    ].copy()

    if comparative_df.empty:
        return

    conservation_columns = [
        f"{label}_conserved"
        for label in HOTSPOTS.values()
    ]

    heatmap_data = (
        comparative_df[
            ["id"] + conservation_columns
        ]
        .set_index("id")
        .astype(int)
    )

    heatmap_data.columns = list(HOTSPOTS.values())

    plt.figure(
        figsize=(
            10,
            max(5, len(heatmap_data) * 0.35),
        )
    )

    sns.heatmap(
        heatmap_data,
        vmin=0,
        vmax=1,
        annot=True,
        fmt="d",
        linewidths=0.5,
        cbar_kws={"label": "Exact residue identity"},
    )

    plt.xlabel("Human TP53 hotspot")
    plt.ylabel("Comparative sequence")
    plt.title(
        "Exact Amino-Acid Identity at Human TP53 Hotspots"
    )
    plt.tight_layout()

    plt.savefig(
        HOTSPOT_FIGURE,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


def plot_clusters(cluster_df: pd.DataFrame) -> None:
    """Plot exploratory clustering of comparative sequences."""
    if "cluster" not in cluster_df.columns:
        return

    plot_df = cluster_df[
        (~cluster_df["is_human_reference"].astype(bool))
        & cluster_df["cluster"].notna()
    ].copy()

    if len(plot_df) < 2:
        return

    plt.figure(figsize=(9, 6))

    sns.scatterplot(
        data=plot_df,
        x="identity_to_human_percent",
        y="hotspot_conservation_percent",
        hue="cluster",
        style="sequence_category",
        s=100,
    )

    plt.xlabel("Global identity to human TP53 (%)")
    plt.ylabel("Exact hotspot residue identity (%)")
    plt.title(
        "Exploratory Clustering of TP53-Related Sequences"
    )
    plt.tight_layout()

    plt.savefig(
        CLUSTER_FIGURE,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


# ============================================================
# 16. JSON-SAFE SUMMARY
# ============================================================

def clean_for_json(value):
    """Convert NumPy/Pandas values to JSON-compatible Python values."""
    if isinstance(value, dict):
        return {
            str(k): clean_for_json(v)
            for k, v in value.items()
        }

    if isinstance(value, list):
        return [clean_for_json(v) for v in value]

    if isinstance(value, (np.integer,)):
        return int(value)

    if isinstance(value, (np.floating,)):
        return (
            None
            if np.isnan(value)
            else float(value)
        )

    if isinstance(value, (np.bool_,)):
        return bool(value)

    return value


def build_summary(
    input_metadata: dict,
    human_record: dict,
    human_validation: dict,
    hotspot_summary_df: pd.DataFrame,
    clustering_metadata: dict,
    excluded_records: list[dict],
) -> dict:
    """Build a machine-readable run summary."""
    hotspot_summary = {}

    for _, row in hotspot_summary_df.iterrows():
        hotspot_summary[row["hotspot"]] = {
            "human_position": int(row["human_position"]),
            "human_residue": row["human_residue"],
            "comparative_sequences": int(
                row["comparative_sequences"]
            ),
            "mapped": int(row["mapped"]),
            "conserved": int(row["conserved"]),
            "substituted": int(row["substituted"]),
            "gaps": int(row["gaps"]),
            "unmapped": int(row["unmapped"]),
            "exact_residue_identity_percent": (
                None
                if pd.isna(
                    row["exact_residue_identity_percent"]
                )
                else float(
                    row["exact_residue_identity_percent"]
                )
            ),
        }

    summary = {
        "project": PROJECT_NAME,
        "analysis": ANALYSIS_NAME,
        "analysis_timestamp_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "input": {
            "file": relative_path(CLEAN_FASTA),
            **input_metadata,
        },
        "reference": {
            "accession": HUMAN_REFERENCE_ACCESSION,
            "sequence_id": human_record["id"],
            "length": len(human_record["sequence"]),
            "expected_length": EXPECTED_HUMAN_LENGTH,
            "validation": human_validation,
        },
        "hotspots": {
            str(position): label
            for position, label in HOTSPOTS.items()
        },
        "hotspot_identity_summary": hotspot_summary,
        "clustering": clustering_metadata,
        "excluded_records": excluded_records,
        "outputs": {
            "sequence_features": relative_path(
                SEQUENCE_FEATURES_OUTPUT
            ),
            "hotspot_mapping": relative_path(
                HOTSPOT_OUTPUT
            ),
            "hotspot_identity_summary": relative_path(
                HOTSPOT_SUMMARY_OUTPUT
            ),
            "comparative_features": relative_path(
                COMPARATIVE_OUTPUT
            ),
            "clustered_features": relative_path(
                CLUSTER_OUTPUT
            ),
            "excluded_sequences": relative_path(
                EXCLUDED_OUTPUT
            ),
            "summary": relative_path(
                SUMMARY_OUTPUT
            ),
        },
        "figures": {
            "identity": relative_path(IDENTITY_FIGURE),
            "hotspot_identity": relative_path(HOTSPOT_FIGURE),
            "clustering": relative_path(CLUSTER_FIGURE),
        },
        "interpretation": {
            "hotspot_metric": (
                "Exact amino-acid identity at an "
                "alignment-mapped human TP53 hotspot position."
            ),
            "comparative_denominator": (
                "Human TP53 reference is excluded; "
                "only mapped comparative residues contribute."
            ),
            "scope": (
                "Exploratory comparative sequence analysis; "
                "not a clinical or validated predictive model."
            ),
        },
    }

    return clean_for_json(summary)


# ============================================================
# 17. MAIN WORKFLOW
# ============================================================

def main() -> int:
    print_header("TP53 COMPARATIVE ANALYSIS")

    print(f"Repository root : {ROOT}")
    print(f"Input           : {relative_path(CLEAN_FASTA)}")
    print(f"Results         : {relative_path(RESULTS_DIR)}")
    print(f"Figures         : {relative_path(FIGURES_DIR)}")

    # --------------------------------------------------------
    # Load and validate input
    # --------------------------------------------------------
    print_header("1. INPUT VALIDATION")

    records, excluded_records, input_metadata = (
        load_and_validate_fasta(CLEAN_FASTA)
    )

    print(
        f"Input records          : "
        f"{input_metadata['input_records']}"
    )
    print(
        f"Valid records retained : "
        f"{input_metadata['retained_records']}"
    )
    print(
        f"Records excluded       : "
        f"{input_metadata['excluded_records']}"
    )

    if not records:
        raise RuntimeError(
            "No valid protein sequences remain after validation."
        )

    # Save exclusions even when empty for auditability.
    pd.DataFrame(
        excluded_records,
        columns=["id", "reason", "details"],
    ).to_csv(
        EXCLUDED_OUTPUT,
        index=False,
    )

    # --------------------------------------------------------
    # Human reference
    # --------------------------------------------------------
    print_header("2. HUMAN TP53 REFERENCE")

    human_record = identify_human_reference(records)
    human_id = human_record["id"]
    human_sequence = human_record["sequence"]

    print(f"Reference ID : {human_id}")
    print(f"Length       : {len(human_sequence)} aa")

    human_validation = validate_human_reference(
        human_record
    )

    if not human_validation["length_pass"]:
        raise RuntimeError(
            "The detected P04637 sequence is not 393 aa. "
            "Verify data/processed/TP53_clean.fasta before "
            "interpreting hotspot coordinates."
        )

    for label, details in human_validation["hotspots"].items():
        print(
            f"{label}: observed={details['observed']} "
            f"expected={details['expected']} "
            f"[{details['status']}]"
        )

        if not details["pass"]:
            raise RuntimeError(
                f"Human TP53 hotspot validation failed for {label}. "
                "Do not proceed with hotspot interpretation."
            )

    # --------------------------------------------------------
    # Sequence features
    # --------------------------------------------------------
    print_header("3. SEQUENCE FEATURES")

    sequence_df = build_sequence_features(
        records,
        human_id,
    )

    print(
        f"Sequence feature rows: {len(sequence_df)}"
    )

    # --------------------------------------------------------
    # Pairwise identity
    # --------------------------------------------------------
    print_header("4. GLOBAL PAIRWISE IDENTITY")

    identity_df = build_identity_features(
        records,
        human_sequence,
    )

    # --------------------------------------------------------
    # Hotspot mapping
    # --------------------------------------------------------
    print_header("5. HOTSPOT MAPPING")

    hotspot_df = build_hotspot_mapping(
        records,
        human_id,
        human_sequence,
    )

    # --------------------------------------------------------
    # Comparative hotspot statistics
    # --------------------------------------------------------
    print_header(
        "6. COMPARATIVE HOTSPOT IDENTITY"
    )

    hotspot_summary_df = build_hotspot_identity_summary(
        hotspot_df,
        human_id,
        human_sequence,
    )

    print(
        hotspot_summary_df.to_string(index=False)
    )

    # --------------------------------------------------------
    # Complete feature table
    # --------------------------------------------------------
    print_header("7. COMPLETE FEATURE TABLE")

    feature_df = build_complete_feature_table(
        sequence_df,
        identity_df,
        hotspot_df,
    )

    # --------------------------------------------------------
    # Save primary tables
    # --------------------------------------------------------
    print_header("8. SAVE RESULTS")

    sequence_df.to_csv(
        SEQUENCE_FEATURES_OUTPUT,
        index=False,
    )

    hotspot_df.to_csv(
        HOTSPOT_OUTPUT,
        index=False,
    )

    hotspot_summary_df.to_csv(
        HOTSPOT_SUMMARY_OUTPUT,
        index=False,
    )

    feature_df.to_csv(
        COMPARATIVE_OUTPUT,
        index=False,
    )

    print("Saved:")
    print(f" - {relative_path(SEQUENCE_FEATURES_OUTPUT)}")
    print(f" - {relative_path(HOTSPOT_OUTPUT)}")
    print(f" - {relative_path(HOTSPOT_SUMMARY_OUTPUT)}")
    print(f" - {relative_path(COMPARATIVE_OUTPUT)}")
    print(f" - {relative_path(EXCLUDED_OUTPUT)}")

    # --------------------------------------------------------
    # Exploratory clustering
    # --------------------------------------------------------
    print_header("9. EXPLORATORY CLUSTERING")

    cluster_df, clustering_metadata = (
        exploratory_clustering(feature_df)
    )

    cluster_df.to_csv(
        CLUSTER_OUTPUT,
        index=False,
    )

    if clustering_metadata["performed"]:
        print(
            f"Selected k          : "
            f"{clustering_metadata['selected_k']}"
        )
        print(
            f"Silhouette score     : "
            f"{clustering_metadata['silhouette_score']:.4f}"
        )
        print(
            "Human reference      : excluded from clustering"
        )
    else:
        print(
            "Clustering was not performed because the "
            "comparative dataset was insufficient."
        )

    # --------------------------------------------------------
    # Figures
    # --------------------------------------------------------
    print_header("10. FIGURES")

    configure_plot_style()

    plot_identity(feature_df)
    plot_hotspot_conservation(
        hotspot_df,
        human_id,
    )
    plot_clusters(cluster_df)

    print(f" - {relative_path(IDENTITY_FIGURE)}")
    print(f" - {relative_path(HOTSPOT_FIGURE)}")

    if CLUSTER_FIGURE.exists():
        print(f" - {relative_path(CLUSTER_FIGURE)}")

    # --------------------------------------------------------
    # Machine-readable summary
    # --------------------------------------------------------
    print_header("11. RUN SUMMARY")

    summary = build_summary(
        input_metadata=input_metadata,
        human_record=human_record,
        human_validation=human_validation,
        hotspot_summary_df=hotspot_summary_df,
        clustering_metadata=clustering_metadata,
        excluded_records=excluded_records,
    )

    SUMMARY_OUTPUT.write_text(
        json.dumps(
            summary,
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    print(
        f"Summary JSON: {relative_path(SUMMARY_OUTPUT)}"
    )

    # --------------------------------------------------------
    # Final report
    # --------------------------------------------------------
    print_header("ANALYSIS COMPLETED")

    print(
        f"Human reference        : {human_id}"
    )
    print(
        f"Human TP53 length      : {len(human_sequence)} aa"
    )
    print(
        f"Valid sequences         : {len(records)}"
    )
    print(
        f"Comparative sequences   : "
        f"{len(records) - 1}"
    )
    print(
        "Hotspot metric          : "
        "exact alignment-mapped residue identity"
    )
    print(
        "Human reference in hotspot denominator: NO"
    )

    print("\nPrimary outputs:")
    for path in [
        SEQUENCE_FEATURES_OUTPUT,
        HOTSPOT_OUTPUT,
        HOTSPOT_SUMMARY_OUTPUT,
        COMPARATIVE_OUTPUT,
        CLUSTER_OUTPUT,
        EXCLUDED_OUTPUT,
        SUMMARY_OUTPUT,
    ]:
        print(f" - {relative_path(path)}")

    print("\nFigures:")
    for path in [
        IDENTITY_FIGURE,
        HOTSPOT_FIGURE,
        CLUSTER_FIGURE,
    ]:
        if path.exists():
            print(f" - {relative_path(path)}")

    print("\nDONE.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nAnalysis interrupted by user.")
        raise SystemExit(130)
    except Exception as exc:
        print("\nERROR")
        print("=" * 72)
        print(str(exc))
        print("=" * 72)
        raise SystemExit(1)
'''

out = Path("/mnt/data/TP53_Comparative_Analysis_improved.py")
out.write_text(code, encoding="utf-8")

# Syntax validation without executing the biological workflow.
compile(code, str(out), "exec")

print(f"Created: {out}")
print(f"Lines: {len(code.splitlines())}")
print("Syntax check: PASSED")
