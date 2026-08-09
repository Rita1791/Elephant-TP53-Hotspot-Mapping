"""
Tests for global pairwise alignment and reference-coordinate mapping.
"""

import numpy as np

from scripts.TP53_Comparative_Analysis import (
    calculate_pairwise_identity,
    map_reference_position,
)


def test_identical_sequences():
    sequence = "ACDEFGHIKLMNPQRSTVWY"

    identity, matches, comparable = calculate_pairwise_identity(
        sequence,
        sequence,
    )

    assert identity == 1.0
    assert matches == len(sequence)
    assert comparable == len(sequence)


def test_single_substitution():
    reference = "ACDEFG"
    query = "ACDFFG"

    identity, matches, comparable = calculate_pairwise_identity(
        reference,
        query,
    )

    assert matches == 5
    assert comparable == 6
    assert np.isclose(identity, 5 / 6)


def test_position_mapping():
    reference = "ACDEFG"
    query = "ACDEFG"

    result = map_reference_position(
        reference,
        query,
        3,
    )

    assert result["status"] == "mapped"
    assert result["residue"] == "D"


def test_one_based_coordinate_system():
    reference = "ABCDE"
    query = "ABCDE"

    result = map_reference_position(
        reference,
        query,
        1,
    )

    assert result["residue"] == "A"


def test_last_coordinate():
    reference = "ABCDE"
    query = "ABCDE"

    result = map_reference_position(
        reference,
        query,
        5,
    )

    assert result["residue"] == "E"
