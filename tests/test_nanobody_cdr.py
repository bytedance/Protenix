from protenix.nanobody.cdr import annotate_cdrh3


def test_cdrh3_uses_metadata_bounds():
    annotation = annotate_cdrh3(metadata_start="105", metadata_end="117")
    assert annotation.start == 105
    assert annotation.end == 117
    assert annotation.source == "metadata"


def test_cdrh3_missing_bounds_behavior():
    annotation = annotate_cdrh3(sequence=None, metadata_start="", metadata_end="")
    assert annotation.start is None
    assert annotation.end is None
    assert annotation.source == "missing"


def test_cdrh3_heuristic_sequence():
    seq = "Q" * 75 + "CARDRSTYY" + "WGQG" + "Q" * 10
    annotation = annotate_cdrh3(sequence=seq, use_optional_numbering=False)
    assert annotation.source == "heuristic"
    assert annotation.start is not None
    assert annotation.end is not None
    assert annotation.start <= annotation.end
