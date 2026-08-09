"""
Tests for human TP53 reference identification and validation.

These tests protect the biological reference sequence used by the
comparative hotspot analysis.
"""

import pytest

from scripts.TP53_Comparative_Analysis import (
    EXPECTED_HUMAN_LENGTH,
    HUMAN_REFERENCE_ACCESSION,
    HOTSPOTS,
    identify_human_reference,
    validate_human_reference,
)


def test_human_reference_identification():
    """P04637 must be identified explicitly as the human reference."""

    records = [
        {
            "id": "XP_elephant_001",
            "description": "Elephant TP53-related protein",
            "sequence": "M" * 400,
        },
        {
            "id": "sp_P04637_P53_HUMAN",
            "description": "TP53_HUMAN",
            "sequence": "M" * 393,
        },
    ]

    reference = identify_human_reference(records)

    assert reference["id"] == "sp_P04637_P53_HUMAN"


def test_human_reference_requires_p04637():
    """
    The workflow must not silently substitute an arbitrary human protein.
    """

    records = [
        {
            "id": "HUMAN_TP53_OTHER",
            "description": "Human TP53",
            "sequence": "M" * 393,
        }
    ]

    with pytest.raises(RuntimeError):
        identify_human_reference(records)


def test_human_reference_length_validation():
    """Canonical P04637 should be 393 amino acids."""

    sequence = ["A"] * EXPECTED_HUMAN_LENGTH

    records = [
        {
            "id": "sp_P04637_P53_HUMAN",
            "description": "TP53_HUMAN",
            "sequence": "".join(sequence),
        }
    ]

    validation = validate_human_reference(records[0])

    assert validation["length"] == 393
    assert validation["expected_length"] == 393
    assert validation["length_pass"] is True


def test_hotspot_positions_are_defined():
    """The six project hotspots must remain explicitly defined."""

    expected = {
        175: "R175",
        245: "G245",
        248: "R248",
        249: "R249",
        273: "R273",
        282: "R282",
    }

    assert HOTSPOTS == expected


def test_human_hotspot_residues_are_validated():
    """
    Construct a 393-aa sequence containing the expected residues
    at the six hotspot positions and verify validation passes.
    """

    sequence = ["A"] * 393

    sequence[174] = "R"  # R175
    sequence[244] = "G"  # G245
    sequence[247] = "R"  # R248
    sequence[248] = "R"  # R249
    sequence[272] = "R"  # R273
    sequence[281] = "R"  # R282

    record = {
        "id": "sp_P04637_P53_HUMAN",
        "description": "TP53_HUMAN",
        "sequence": "".join(sequence),
    }

    validation = validate_human_reference(record)

    for hotspot, details in validation["hotspots"].items():
        assert details["pass"] is True, (
            f"Human reference validation failed for {hotspot}"
        )


def test_human_reference_accession_is_p04637():
    """The configured reference accession must remain P04637."""

    assert HUMAN_REFERENCE_ACCESSION == "P04637"
