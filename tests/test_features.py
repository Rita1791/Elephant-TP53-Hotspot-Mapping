"""
Tests for sequence feature extraction and metadata classification.
"""

import numpy as np

from scripts.TP53_Comparative_Analysis import (
    amino_acid_composition,
    build_sequence_features,
    classify_sequence,
    extract_accession,
)


def test_amino_acid_composition_sums_to_one():
    """Fractions of the 20 standard amino acids should sum to 1."""

    sequence = "ACDEFGHIKLMNPQRSTVWY"

    composition = amino_acid_composition(sequence)

    total = sum(composition.values())

    assert np.isclose(total, 1.0)


def test_amino_acid_composition_counts_correctly():
    """Each amino acid should have the expected fraction."""

    sequence = "AAAA"

    composition = amino_acid_composition(sequence)

    assert composition["frac_A"] == 1.0

    for amino_acid, fraction in composition.items():
        if amino_acid != "frac_A":
            assert fraction == 0.0


def test_accession_extraction_uniprot():
    """UniProt P04637 should be extracted correctly."""

    accession = extract_accession(
        "sp_P04637_P53_HUMAN"
    )

    assert accession == "P04637"


def test_accession_extraction_refseq():
    """RefSeq-style accession should be retained."""

    accession = extract_accession(
        "XP_049714738.1"
    )

    assert accession == "XP_049714738.1"


def test_human_sequence_classification():
    """Human TP53 should be classified as Human."""

    record = {
        "id": "sp_P04637_P53_HUMAN",
        "description": "TP53_HUMAN",
    }

    assert classify_sequence(record) == "Human"


def test_elephant_sequence_classification():
    """Elephant metadata should produce an elephant category."""

    record = {
        "id": "XP_ELEPHANT_001",
        "description": "Loxodonta africana TP53 protein",
    }

    category = classify_sequence(record)

    assert category in {
        "African_elephant",
        "Elephant",
    }


def test_retrogene_classification():
    """Retrogene metadata should be recognized."""

    record = {
        "id": "ELEPHANT_RTG_001",
        "description": "Elephant TP53 retrogene",
    }

    assert classify_sequence(record) == "Elephant_retrogene"


def test_sequence_feature_table_has_expected_columns():
    """Feature generation should preserve core metadata."""

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

    expected_columns = {
        "id",
        "accession",
        "length",
        "is_human_reference",
        "sequence_category",
        "frac_A",
        "frac_C",
        "frac_D",
        "frac_E",
        "frac_F",
        "frac_G",
        "frac_H",
        "frac_I",
        "frac_K",
        "frac_L",
        "frac_M",
        "frac_N",
        "frac_P",
        "frac_Q",
        "frac_R",
        "frac_S",
        "frac_T",
        "frac_V",
        "frac_W",
        "frac_Y",
    }

    assert expected_columns.issubset(
        set(result.columns)
    )

    assert result.iloc[0]["length"] == 20
