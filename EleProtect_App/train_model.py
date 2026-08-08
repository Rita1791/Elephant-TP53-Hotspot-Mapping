# train_model.py

import pandas as pd
import numpy as np
import os
import joblib

from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import cross_val_score

# =========================================================
# CONFIGURATION
# =========================================================

FEATURE_CSV = "tp53_expanded_training_dataset.csv"

if not os.path.exists(FEATURE_CSV):
    raise FileNotFoundError(
        "Training dataset not found. Run generate_training_dataset.py first."
    )

# =========================================================
# LOAD DATA
# =========================================================

df = pd.read_csv(FEATURE_CSV)

if df.empty:
    raise ValueError("Training dataset is empty.")

# =========================================================
# DEFINE FEATURE COLUMNS (STRICT ORDER)
# =========================================================

FEATURE_COLUMNS = [
    "conservation_score",
    "mutation_frequency",
    "retrogene_variability",
    "sequence_identity",
    "mean_blosum",
    "std_blosum",
    "damaging_fraction",
    "hotspot_count",
    "weighted_mutation_burden"
]

missing = [c for c in FEATURE_COLUMNS if c not in df.columns]
if missing:
    raise ValueError(f"Missing required feature columns: {missing}")

X = df[FEATURE_COLUMNS].fillna(0)

# =========================================================
# TARGET VARIABLE
# =========================================================

if "target_score" in df.columns:
    y = df["target_score"]
else:
    print("No explicit target_score found. Creating composite exploratory target.")
    y = (
        0.6 * (1 - df["conservation_score"])
        + 0.4 * df["mutation_frequency"]
    )

# =========================================================
# MODEL DEFINITION (BOOSTING > RANDOM FOREST)
# =========================================================

model = GradientBoostingRegressor(
    n_estimators=500,
    learning_rate=0.03,
    max_depth=4,
    random_state=42
)

# =========================================================
# TRAIN MODEL
# =========================================================

model.fit(X, y)

# =========================================================
# CROSS-VALIDATION
# =========================================================

cv_scores = cross_val_score(model, X, y, cv=5)

print("Cross-validation R2 scores:", cv_scores)
print("Mean CV R2:", np.mean(cv_scores))

# =========================================================
# SAVE MODEL + FEATURE ORDER
# =========================================================

os.makedirs("models", exist_ok=True)

joblib.dump(model, "models/eleprotect_model.joblib")
joblib.dump(FEATURE_COLUMNS, "models/feature_columns.joblib")

print("Model and feature column order saved successfully.")