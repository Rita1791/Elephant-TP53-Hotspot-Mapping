# model.py

import joblib
import os
import pandas as pd
import numpy as np

MODEL_PATH = "models/eleprotect_model.joblib"
FEATURE_COL_PATH = "models/feature_columns.joblib"


# =========================================================
# ADVANCED FALLBACK SCORING
# =========================================================

def fallback_score(df: pd.DataFrame) -> float:
    """
    Advanced descriptive composite score if ML model unavailable.
    Uses expanded feature set.
    """

    row = df.iloc[0]

    weights = {
        "conservation_score": -0.5,
        "mutation_frequency": 0.4,
        "retrogene_variability": 0.3,
        "sequence_identity": -0.4,
        "mean_blosum": -0.3,
        "std_blosum": 0.2,
        "damaging_fraction": 0.5,
        "hotspot_count": 0.1,
        "weighted_mutation_burden": 0.6
    }

    score = 0.0
    for k, w in weights.items():
        score += w * float(row.get(k, 0))

    return float(score)


# =========================================================
# MAIN PREDICTION FUNCTION
# =========================================================

def predict_score(df_features: pd.DataFrame) -> float:
    """
    Predict using trained model if available.
    Falls back to composite scoring if needed.
    """

    # -----------------------------
    # Attempt ML prediction
    # -----------------------------
    try:
        if os.path.exists(MODEL_PATH) and os.path.exists(FEATURE_COL_PATH):

            model = joblib.load(MODEL_PATH)
            feature_columns = joblib.load(FEATURE_COL_PATH)

            # Ensure feature order consistency
            df = df_features.reindex(columns=feature_columns, fill_value=0.0)

            X = df.astype(float)

            preds = model.predict(X)

            if hasattr(preds, "__len__"):
                return float(np.mean(preds))
            return float(preds)

    except Exception as e:
        print("ML prediction failed. Using fallback:", e)

    # -----------------------------
    # Fallback composite score
    # -----------------------------
    return fallback_score(df_features)