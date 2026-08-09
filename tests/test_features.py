"""
Tests for sequence feature extraction and metadata classification.
"""

import numpy as np

from scripts.TP53_Comparative_Analysis import (
    amino_acid_composition,
    build_sequence_features,
    extract_accession,
    classify_sequence,
)


def test_amino_acid_composition():
    sequence = "ACDEFGHIKLMNPQRSTVWY"

    composition = amino_acid_composition(sequence)

    assert np.isclose(
        sum(composition.values()),
        1.0,
    )


def test_amino_acid_composition_all_a():
    sequence = "AAAA"

    composition = amino_acid_composition(sequence)

    assert composition["frac_A"] == 1.0

    for key, value in composition.items():
        if key != "frac_A":
            assert value == 0.0


def test_uniprot_accession():
    assert (
        extract_accession("sp_P04637_P53_HUMAN")
        == "P04637"
    )


def test_refseq_accession():
    assert (
        extract_accession("XP_049714738.1")
        == "XP_049714738.1"
    )


def test_human_classification():
    record = {
        "id": "sp_P04637_P53_HUMAN",
        "description": "TP53_HUMAN",
    }

    assert classify_sequence(record) == "Human"


def test_retrogene_classification():
    record = {
        "id": "ELEPHANT_RTG_001",
        "description": "Elephant TP53 retrogene",
    }

    assert (
        classify_sequence(record)
        == "Elephant_retrogene"
    )


def test_sequence_features():
    records = [
        {
            "id": "sp_P04637_P53_HUMAN",
            "description": "TP53_HUMAN",
            "sequence": "ACDEFGHIKLMNPQRSTVWY",
        }
    ]

    result = build_sequence_features(
        records,
        "sp_P04637_P53_HUMAN",
    )

    assert result.iloc[0]["length"] == 20
    assert result.iloc[0]["is_human_reference"] is True
    assert "frac_A" in result.columns
    assert "frac_Y" in result.columns
