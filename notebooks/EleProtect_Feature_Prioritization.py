from pathlib import Path
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from Bio import SeqIO

print("Libraries loaded successfully.")

def find_repository_root():

    current = Path.cwd().resolve()

    candidates = [
        current,
        current.parent,
        current.parent.parent,
        current.parent.parent.parent
    ]

    for path in candidates:

        if (
            (path / "data").exists()
            and (path / "README.md").exists()
        ):

            return path

    raise FileNotFoundError(
        "Repository root could not be located."
    )


ROOT = find_repository_root()

INPUT = (
    ROOT
    / "data"
    / "processed"
    / "TP53_clean.fasta"
)

RESULTS_DIR = (
    ROOT / "results"
)

FIGURES_DIR = (
    ROOT / "figures"
)

RESULTS_DIR.mkdir(
    exist_ok=True
)

FIGURES_DIR.mkdir(
    exist_ok=True
)

print("Repository:", ROOT)
print("Input:", INPUT)

if not INPUT.exists():

    raise FileNotFoundError(
        f"Input FASTA not found: {INPUT}"
    )

records = list(
    SeqIO.parse(
        INPUT,
        "fasta"
    )
)

print(
    "Sequences loaded:",
    len(records)
)

VALID_AA = set(
    "ACDEFGHIKLMNPQRSTVWY"
)


def sequence_features(
    sequence
):

    sequence = (
        str(sequence)
        .upper()
        .replace(" ", "")
        .replace("\n", "")
    )

    length = len(sequence)

    if length == 0:

        raise ValueError(
            "Sequence cannot be empty."
        )

    hydrophobic = set(
        "AVILMFWY"
    )

    charged = set(
        "DEKR"
    )

    aromatic = set(
        "FWY"
    )

    polar = set(
        "STNQ"
    )

    return {

        "length":
            length,

        "hydrophobic_fraction":
            sum(
                aa in hydrophobic
                for aa in sequence
            ) / length,

        "charged_fraction":
            sum(
                aa in charged
                for aa in sequence
            ) / length,

        "aromatic_fraction":
            sum(
                aa in aromatic
                for aa in sequence
            ) / length,

        "polar_fraction":
            sum(
                aa in polar
                for aa in sequence
            ) / length,

        "proline_fraction":
            sequence.count("P")
            / length,

        "glycine_fraction":
            sequence.count("G")
            / length,

        "cysteine_fraction":
            sequence.count("C")
            / length,

        "lysine_fraction":
            sequence.count("K")
            / length,

        "arginine_fraction":
            sequence.count("R")
            / length
    }

rows = []

for record in records:

    sequence = (
        str(record.seq)
        .upper()
        .replace(" ", "")
        .replace("\n", "")
    )

    invalid = (
        set(sequence)
        - VALID_AA
    )

    if invalid:

        print(
            f"Skipping {record.id}: "
            f"invalid residues = "
            f"{sorted(invalid)}"
        )

        continue

    row = {
        "id": record.id
    }

    row.update(
        sequence_features(
            sequence
        )
    )

    rows.append(row)


feature_df = pd.DataFrame(
    rows
)

feature_df

numeric_features = [
    column
    for column in feature_df.columns
    if column != "id"
]

feature_summary = (
    feature_df[
        numeric_features
    ]
    .describe()
    .T
)

feature_summary

plt.figure(
    figsize=(12, 7)
)

sns.boxplot(
    data=feature_df[
        numeric_features
    ]
)

plt.xticks(
    rotation=45,
    ha="right"
)

plt.ylabel(
    "Feature value"
)

plt.title(
    "Distribution of Sequence-Derived Features"
)

plt.tight_layout()

plt.savefig(
    FIGURES_DIR
    / "eleprotect_feature_distributions.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

correlation = (
    feature_df[
        numeric_features
    ]
    .corr()
)

plt.figure(
    figsize=(10, 8)
)

sns.heatmap(
    correlation,
    annot=True,
    fmt=".2f",
    center=0,
    linewidths=0.5
)

plt.title(
    "Correlation Between Sequence-Derived Features"
)

plt.tight_layout()

plt.savefig(
    FIGURES_DIR
    / "eleprotect_feature_correlation.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

from sklearn.preprocessing import MinMaxScaler

scaler = MinMaxScaler()

scaled_values = (
    scaler.fit_transform(
        feature_df[
            numeric_features
        ]
    )
)

scaled_df = pd.DataFrame(
    scaled_values,
    columns=numeric_features
)

feature_df[
    "exploratory_feature_score"
] = (
    scaled_df.mean(
        axis=1
    )
)

ranking = (
    feature_df[
        [
            "id",
            "exploratory_feature_score"
        ]
    ]
    .sort_values(
        "exploratory_feature_score",
        ascending=False
    )
)

ranking

output = (
    RESULTS_DIR
    / "eleprotect_sequence_features.csv"
)

feature_df.to_csv(
    output,
    index=False
)

print(
    "Saved:",
    output
)

