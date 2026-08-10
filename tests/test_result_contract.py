import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
HOTSPOTS = {"R175": 175, "G245": 245, "R248": 248, "R249": 249, "R273": 273, "R282": 282}


def test_all_primary_outputs_exist_and_are_nonempty():
    expected = [
        RESULTS / "tp53_sequence_features.csv",
        RESULTS / "tp53_hotspot_mapping.csv",
        RESULTS / "tp53_hotspot_identity_summary.csv",
        RESULTS / "tp53_comparative_features.csv",
        RESULTS / "tp53_comparative_features_clustered.csv",
        RESULTS / "tp53_excluded_sequences.csv",
        RESULTS / "tp53_summary.json",
        FIGURES / "tp53_identity_to_human.png",
        FIGURES / "tp53_hotspot_conservation.png",
    ]
    for path in expected:
        assert path.is_file(), f"Missing output: {path.relative_to(ROOT)}"
        assert path.stat().st_size > 0, f"Empty output: {path.relative_to(ROOT)}"


def test_hotspot_summary_excludes_human_from_comparative_denominator():
    frame = pd.read_csv(RESULTS / "tp53_hotspot_identity_summary.csv")
    assert dict(zip(frame["hotspot"], frame["human_position"])) == HOTSPOTS
    assert (frame["mapped"] + frame["gaps"] + frame["unmapped"] == frame["comparative_sequences"]).all()
    assert (frame["conserved"] <= frame["mapped"]).all()
    assert frame["exact_residue_identity_percent"].dropna().between(0, 100).all()


def test_mapping_contains_exact_status_and_identity_fields():
    frame = pd.read_csv(RESULTS / "tp53_hotspot_mapping.csv")
    assert frame["is_human_reference"].sum() == 1
    for label in HOTSPOTS:
        assert {f"{label}_human", f"{label}_query", f"{label}_status", f"{label}_conserved"} <= set(frame.columns)
        assert set(frame[f"{label}_status"]) <= {"mapped", "gap", "unmapped"}


def test_summary_states_research_not_clinical_scope():
    payload = json.loads((RESULTS / "tp53_summary.json").read_text(encoding="utf-8"))
    scope = payload["interpretation"]["scope"].lower()
    assert "not a clinical" in scope

