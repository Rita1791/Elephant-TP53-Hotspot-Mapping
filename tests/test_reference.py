"""
Tests for human TP53 reference identification and validation.

These tests protect the biological reference sequence used by the
comparative hotspot analysis.
"""

import pytest

from scripts.TP53_Comparative_Analysis import (
    HUMAN_REFERENCE_ACCESSION,
    EXPECTED_HUMAN_LENGTH,
    HOTSPOTS,
    identify_human_reference,
    validate_human_reference,
)


def test_reference_accession():
    assert HUMAN_REFERENCE_ACCESSION == "P04637"


def test_reference_length():
    assert EXPECTED_HUMAN_LENGTH == 393


def test_hotspot_definition():
    assert HOTSPOTS == {
        175: "R175",
        245: "G245",
        248: "R248",
        249: "R249",
        273: "R273",
        282: "R282",
    }


def test_identify_human_reference():
    records = [
        {
            "id": "XP_ELEPHANT_001",
            "description": "Elephant TP53",
            "sequence": "A" * 400,
        },
        {
            "id": "sp_P04637_P53_HUMAN",
            "description": "TP53_HUMAN",
            "sequence": "A" * 393,
        },
    ]

    result = identify_human_reference(records)

    assert result["id"] == "sp_P04637_P53_HUMAN"


def test_reference_validation():
    sequence = ["A"] * 393

    sequence[174] = "R"
    sequence[244] = "G"
    sequence[247] = "R"
    sequence[248] = "R"
    sequence[272] = "R"
    sequence[281] = "R"

    record = {
        "id": "sp_P04637_P53_HUMAN",
        "description": "TP53_HUMAN",
        "sequence": "".join(sequence),
    }

    validation = validate_human_reference(record)

    assert validation["length"] == 393
    assert validation["length_pass"] is True

    for hotspot in validation["hotspots"].values():
        assert hotspot["pass"] is True

def test_human_reference_accession_is_p04637():
    """The configured reference accession must remain P04637."""

    assert HUMAN_REFERENCE_ACCESSION == "P04637"
