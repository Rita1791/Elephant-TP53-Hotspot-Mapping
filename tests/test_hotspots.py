"""
Tests for TP53 hotspot mapping and comparative hotspot statistics.
"""

import pandas as pd

from scripts.TP53_Comparative_Analysis import (
    HOTSPOTS,
    build_hotspot_identity_summary,
    build_hotspot_mapping,
)


def make_reference_sequence():
    """Create a minimal test sequence containing all six hotspot residues."""

    sequence = ["A"] * 393

    sequence[174] = "R"
    sequence[244] = "G"
    sequence[247] = "R"
    sequence[248] = "R"
    sequence[272] = "R"
    sequence[281] = "R"

    return "".join(sequence)


def test_hotspot_mapping_contains_all_hotspots():
    """Every predefined hotspot should appear in the output."""

    reference = make_reference_sequence()

    records = [
        {
            "id": "sp_P04637_P53_HUMAN",
            "description": "TP53_HUMAN",
            "sequence": reference,
        },
        {
            "id": "ELEPHANT_001",
            "description": "Elephant TP53",
            "sequence": reference,
        },
    ]

    result = build_hotspot_mapping(
        records,
        "sp_P04637_P53_HUMAN",
        reference,
    )

    assert len(result) == 2

    row = result[
        result["id"] == "ELEPHANT_001"
    ].iloc[0]

    for label in HOTSPOTS.values():
        assert f"{label}_human" in result.columns
        assert f"{label}_query" in result.columns
        assert f"{label}_status" in result.columns
        assert f"{label}_conserved" in result.columns

        assert row[f"{label}_conserved"] is True


def test_human_reference_is_excluded_from_comparative_summary():
    """
    The human reference must not contribute to the comparative denominator.
    """

    reference = make_reference_sequence()

    records = [
        {
            "id": "sp_P04637_P53_HUMAN",
            "description": "TP53_HUMAN",
            "sequence": reference,
        },
        {
            "id": "ELEPHANT_001",
            "description": "Elephant TP53",
            "sequence": reference,
        },
    ]

    hotspot_df = build_hotspot_mapping(
        records,
        "sp_P04637_P53_HUMAN",
        reference,
    )

    summary = build_hotspot_identity_summary(
        hotspot_df,
        "sp_P04637_P53_HUMAN",
        reference,
    )

    assert all(
        summary["comparative_sequences"] == 1
    )


def test_conserved_and_substituted_counts():
    """
    Create one identical and one substituted comparative sequence
    and verify the hotspot accounting.
    """

    reference = make_reference_sequence()

    substituted = list(reference)

    # R175 -> K
    substituted[174] = "K"

    records = [
        {
            "id": "sp_P04637_P53_HUMAN",
            "description": "TP53_HUMAN",
            "sequence": reference,
        },
        {
            "id": "ELEPHANT_CONSERVED",
            "description": "Elephant TP53 conserved",
            "sequence": reference,
        },
        {
            "id": "ELEPHANT_SUBSTITUTED",
            "description": "Elephant TP53 substituted",
            "sequence": "".join(substituted),
        },
    ]

    hotspot_df = build_hotspot_mapping(
        records,
        "sp_P04637_P53_HUMAN",
        reference,
    )

    summary = build_hotspot_identity_summary(
        hotspot_df,
        "sp_P04637_P53_HUMAN",
        reference,
    )

    r175 = summary[
        summary["hotspot"] == "R175"
    ].iloc[0]

    assert r175["comparative_sequences"] == 2
    assert r175["mapped"] == 2
    assert r175["conserved"] == 1
    assert r175["substituted"] == 1
    assert r175["gaps"] == 0
    assert r175["unmapped"] == 0
    assert r175["exact_residue_identity_percent"] == 50.0


def test_all_six_hotspots_are_present_in_summary():
    """The summary must contain exactly the six defined hotspots."""

    reference = make_reference_sequence()

    records = [
        {
            "id": "sp_P04637_P53_HUMAN",
            "description": "TP53_HUMAN",
            "sequence": reference,
        },
        {
            "id": "ELEPHANT_001",
            "description": "Elephant TP53",
            "sequence": reference,
        },
    ]

    hotspot_df = build_hotspot_mapping(
        records,
        "sp_P04637_P53_HUMAN",
        reference,
    )

    summary = build_hotspot_identity_summary(
        hotspot_df,
        "sp_P04637_P53_HUMAN",
        reference,
    )

    assert set(summary["hotspot"]) == set(
        HOTSPOTS.values()
    )

    assert len(summary) == 6
