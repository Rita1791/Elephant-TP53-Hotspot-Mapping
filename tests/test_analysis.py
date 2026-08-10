import math

from scripts.TP53_Comparative_Analysis import (
    calculate_pairwise_identity,
    clean_sequence,
    map_reference_position,
)


def test_clean_sequence_normalizes_case_and_whitespace():
    assert clean_sequence(" ac d\nE ") == "ACDE"


def test_pairwise_identity_is_exact_for_identical_sequences():
    identity, matches, comparable = calculate_pairwise_identity("ACDE", "ACDE")
    assert math.isclose(identity, 1.0)
    assert matches == 4
    assert comparable == 4


def test_reference_coordinate_mapping_handles_query_deletion():
    assert map_reference_position("ACD", "AD", 1) == {"status": "mapped", "residue": "A"}
    assert map_reference_position("ACD", "AD", 2) == {"status": "gap", "residue": "-"}
    assert map_reference_position("ACD", "AD", 3) == {"status": "mapped", "residue": "D"}


def test_out_of_range_reference_coordinate_is_unmapped():
    assert map_reference_position("ACD", "ACD", 4) == {
        "status": "unmapped",
        "residue": None,
    }

