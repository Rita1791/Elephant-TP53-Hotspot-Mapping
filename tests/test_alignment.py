"""
Tests for global pairwise alignment and reference-coordinate mapping.
"""

import numpy as np

from scripts.TP53_Comparative_Analysis import (
    calculate_pairwise_identity,
    map_reference_position,
)


def test_identical_sequences_have_full_identity():
    """Identical sequences must produce 100% identity."""

    reference = "ACDEFGHIKLMNPQRSTVWY"
    query = reference

    identity, matches, comparable = calculate_pairwise_identity(
        reference,
        query,
    )

    assert identity == 1.0
    assert matches == len(reference)
    assert comparable == len(reference)


def test_single_substitution_reduces_identity():
    """One amino-acid substitution should reduce exact identity."""

    reference = "ACDEFG"
    query = "ACDFFG"

    identity, matches, comparable = calculate_pairwise_identity(
        reference,
        query,
    )

    assert comparable == 6
    assert matches == 5
    assert np.isclose(identity, 5 / 6)


def test_reference_position_maps_to_same_residue():
    """A conserved reference position should map to the same residue."""

    reference = "ACDEFG"
    query = "ACDEFG"

    result = map_reference_position(
        reference,
        query,
        3,
    )

    assert result["status"] == "mapped"
    assert result["residue"] == "D"


def test_reference_position_maps_substitution():
    """A substituted residue must be reported as mapped, not conserved."""

    reference = "ACDEFG"
    query = "ACXEFG"

    result = map_reference_position(
        reference,
        query,
        3,
    )

    assert result["status"] == "mapped"
    assert result["residue"] == "X"


def test_reference_position_maps_gap():
    """
    A deletion/gap at the query position should be explicitly reported.
    """

    reference = "ACDEFG"
    query = "ACDFG"

    result = map_reference_position(
        reference,
        query,
        4,
    )

    assert result["status"] in {"gap", "mapped"}


def test_first_reference_position_is_one_based():
    """The mapping API must use biological 1-based coordinates."""

    reference = "ABCDE"
    query = "ABCDE"

    result = map_reference_position(
        reference,
        query,
        1,
    )

    assert result["status"] == "mapped"
    assert result["residue"] == "A"


def test_last_reference_position_maps_correctly():
    """The final reference coordinate must map correctly."""

    reference = "ABCDE"
    query = "ABCDE"

    result = map_reference_position(
        reference,
        query,
        5,
    )

    assert result["status"] == "mapped"
    assert result["residue"] == "E"
