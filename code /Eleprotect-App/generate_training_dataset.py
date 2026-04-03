# generate_training_dataset.py

import pandas as pd
import numpy as np
import os

from utils import build_model_features

np.random.seed(42)

# =========================================================
# Synthetic dataset size (increase for stronger model)
# =========================================================

N_SAMPLES = 800   # <-- increase to 1000 if you want stronger ML

# =========================================================
# Generate synthetic hotspot feature distributions
# =========================================================

rows = []

for i in range(N_SAMPLES):

    conservation = np.random.uniform(0.2, 1.0)
    mutation_freq = np.random.uniform(0.01, 0.08)
    mean_blosum = np.random.normal(3, 2)
    std_blosum = abs(np.random.normal(1.5, 0.8))
    damaging_fraction = np.random.uniform(0, 1)
    hotspot_count = np.random.randint(3, 7)
    weighted_burden = mutation_freq * abs(mean_blosum)

    retrogene_variability = 1 - conservation
    sequence_identity = conservation

    # Synthetic biological risk target
    target_score = (
        0.6 * damaging_fraction +
        0.4 * weighted_burden +
        0.3 * mutation_freq -
        0.5 * conservation
    )

    rows.append({
        "conservation_score": conservation,
        "mutation_frequency": mutation_freq,
        "retrogene_variability": retrogene_variability,
        "sequence_identity": sequence_identity,
        "mean_blosum": mean_blosum,
        "std_blosum": std_blosum,
        "damaging_fraction": damaging_fraction,
        "hotspot_count": hotspot_count,
        "weighted_mutation_burden": weighted_burden,
        "target_score": target_score
    })

df = pd.DataFrame(rows)

# =========================================================
# Save dataset
# =========================================================

df.to_csv("tp53_expanded_training_dataset.csv", index=False)

print("Expanded synthetic training dataset created successfully.")
print("Samples:", len(df))