"""Truth denominators, strand/locus checks and full-span candidate recall."""

import json

import pytest

from genocache.evaluate import evaluate_alignments, evaluate_retrieval, load_truth


def write_jsonl(path, rows):
    path.write_text("".join(json.dumps(row) + "\n" for row in rows))
    return path


def truth_row(name, start=1000, *, strand="+", chrom="chr1", stratum="unique",
              read_length=200, reference_span=200, origins=None):
    return {"read_id": name, "read_length": read_length, "reference_span": reference_span,
            "stratum": stratum,
            "origins": origins if origins is not None else [{"chrom": chrom, "start": start, "strand": strand}]}


def anchor_rows(name, *, origin=1000, strand="+", read_length=200, chrom="chr1"):
    rows = []
    for qs in (0, 100):
        qe = qs + 50
        oriented_query = qs if strand == "+" else read_length - qe
        rows.append({"read_id": name, "read_length": read_length,
                     "query_start": qs, "query_end": qe, "seed_id": str(qs), "chrom": chrom,
                     "ref_start": origin + oriented_query, "ref_end": origin + oriented_query + 50,
                     "strand": strand, "similarity": 0.9})
    return rows


def write_sam(path, records):
    lines = ["@HD\tVN:1.6\tSO:unsorted", "@SQ\tSN:chr1\tLN:5000", "@SQ\tSN:chr2\tLN:5000"]
    for name, chrom, start, flag, mapq in records:
        mapped = not flag & 4
        lines.append("\t".join(map(str, [name, flag, chrom if mapped else "*",
                                            start + 1 if mapped else 0, mapq if mapped else 0,
                                            "100M" if mapped else "*", "*", 0, 0, "A" * 100, "*"])))
    path.write_text("\n".join(lines) + "\n")
    return path


def test_mapping_rate_is_distinct_from_truth_recall_and_missing_reads_remain(tmp_path):
    truth = [
        truth_row("correct", 100), truth_row("wrong_pos", 100, stratum="repeats"),
        truth_row("wrong_strand", 300, stratum="repeats"),
        truth_row("unmapped", 400, stratum="hard"), truth_row("missing", 500, stratum="hard"),
        truth_row("negative", origins=[], stratum="negative"),
        truth_row("ambiguous", stratum="repeats", origins=[
            {"chrom": "chr1", "start": 700, "strand": "+"},
            {"chrom": "chr2", "start": 800, "strand": "-"},
        ]),
        truth_row("unknown_mapq", 900, stratum="quality_unknown"),
    ]
    sam = write_sam(tmp_path / "reads.sam", [
        ("correct", "chr1", 100, 0, 60), ("wrong_pos", "chr1", 1000, 0, 60),
        ("wrong_strand", "chr1", 300, 16, 60), ("unmapped", "*", 0, 4, 0),
        ("negative", "chr2", 600, 0, 60), ("ambiguous", "chr2", 800, 16, 60),
        ("unknown_mapq", "chr1", 1000, 0, 255),
    ])
    report = evaluate_alignments(sam, write_jsonl(tmp_path / "truth.jsonl", truth), tolerance=0)
    overall = report["overall"]
    assert overall["reads"] == 8
    assert overall["reads_with_known_origins"] == 7
    assert overall["correct"] == 2
    assert overall["recall"] == pytest.approx(2 / 7)
    assert overall["mapped"] == 6
    assert overall["mapping_rate"] == 0.75
    assert overall["missing_primary_records"] == 1
    assert overall["false_mapped_known_unmappable"] == 1
    assert report["strata"]["hard"]["recall"] == 0
    assert report["strata"]["negative"]["recall"] is None
    assert not report["cigar_or_variant_accuracy_evaluated"]
    assert report["mapq_ge_20"]["reads"] == 5
    assert report["mapq_ge_20"]["wrong"] == 3


def test_per_stratum_uncertainty_is_a_wilson_record_not_a_scalar(tmp_path):
    truth = [truth_row("yes", 100, stratum="easy"), truth_row("no", 200, stratum="hard")]
    sam = write_sam(tmp_path / "reads.sam", [("yes", "chr1", 100, 0, 60)])
    report = evaluate_alignments(sam, write_jsonl(tmp_path / "truth.jsonl", truth), tolerance=0)
    for stratum, successes in (("easy", 1), ("hard", 0)):
        interval = report["strata"][stratum]["recall_wilson_95"]
        assert isinstance(interval, dict)
        assert interval["confidence"] == 0.95
        assert interval["successes"] == successes
        assert interval["trials"] == 1
        assert 0 <= interval["lower"] < interval["upper"] <= 1
    assert report["strata"]["hard"]["recall_wilson_95"]["upper"] > 0.5


def test_secondary_correct_locus_does_not_rescue_wrong_primary_origin(tmp_path):
    truth = write_jsonl(tmp_path / "truth.jsonl", [truth_row("read", 100)])
    sam = write_sam(tmp_path / "reads.sam", [
        ("read", "chr1", 1000, 0, 60), ("read", "chr1", 100, 256, 20),
        ("read", "chr1", 100, 2048, 20),
    ])
    report = evaluate_alignments(sam, truth, tolerance=0)
    assert report["overall"]["correct"] == 0
    assert report["overall"]["secondary_records"] == 1
    assert report["overall"]["supplementary_records"] == 1


def test_alignment_origin_evaluation_does_not_require_unmeasured_read_lengths(tmp_path):
    truth = write_jsonl(tmp_path / "truth.jsonl", [{"read_id": "read", "origins": [
        {"chrom": "chr1", "start": 100, "strand": "+"}]}])
    sam = write_sam(tmp_path / "reads.sam", [("read", "chr1", 100, 0, 60)])
    assert evaluate_alignments(sam, truth, tolerance=0)["overall"]["correct"] == 1


def test_retrieval_requires_full_span_and_correct_locus_and_strand(tmp_path):
    truth = [
        truth_row("correct"), truth_row("missing", stratum="hard"),
        truth_row("wrong_locus", stratum="repeats"), truth_row("wrong_strand", stratum="repeats"),
        truth_row("span_not_contained", reference_span=240, stratum="indels"),
        truth_row("negative", origins=[], stratum="negative"),
    ]
    anchors = (anchor_rows("correct") + anchor_rows("wrong_locus", origin=1500)
               + anchor_rows("wrong_strand", strand="-") + anchor_rows("span_not_contained")
               + anchor_rows("negative"))
    report = evaluate_retrieval(
        write_jsonl(tmp_path / "anchors.jsonl", anchors),
        write_jsonl(tmp_path / "truth.jsonl", truth), {"chr1": 5000}, positional_slop=0,
    )
    assert report["overall"]["reads_with_known_origins"] == 5
    assert report["overall"]["correct"] == 1
    assert report["overall"]["recall"] == 0.2
    by_name = {row["read_id"]: row for row in report["per_read"]}
    assert by_name["missing"]["candidate_count"] == 0
    assert not by_name["span_not_contained"]["candidate_contains_origin"]
    assert by_name["negative"]["candidate_count"] == 1
    assert not report["authoritative_alignment"]


def test_empty_anchors_are_counted_as_misses_for_every_known_origin(tmp_path):
    truth = [truth_row("one"), truth_row("two", stratum="hard"), truth_row("negative", origins=[])]
    report = evaluate_retrieval(
        write_jsonl(tmp_path / "anchors.jsonl", []), write_jsonl(tmp_path / "truth.jsonl", truth),
        {"chr1": 5000}, positional_slop=0,
    )
    assert report["overall"]["reads"] == 3
    assert report["overall"]["reads_with_known_origins"] == 2
    assert report["overall"]["recall"] == 0
    assert all(row["candidate_count"] == 0 for row in report["per_read"])


@pytest.mark.parametrize("field,value", [("read_length", None), ("read_length", 0),
                                         ("read_length", True), ("reference_span", None),
                                         ("reference_span", 0), ("reference_span", -1)])
def test_retrieval_rejects_missing_or_invalid_truth_span_metadata(tmp_path, field, value):
    truth = truth_row("read")
    if value is None:
        truth.pop(field)
    else:
        truth[field] = value
    with pytest.raises(ValueError, match="read_length|reference_span"):
        evaluate_retrieval(write_jsonl(tmp_path / "anchors.jsonl", anchor_rows("read")),
                           write_jsonl(tmp_path / "truth.jsonl", [truth]), {"chr1": 5000})


def test_retrieval_null_origin_requires_read_length_but_not_a_reference_span(tmp_path):
    truth = truth_row("negative", origins=[])
    truth.pop("reference_span")
    report = evaluate_retrieval(write_jsonl(tmp_path / "anchors.jsonl", []),
                                write_jsonl(tmp_path / "truth.jsonl", [truth]), {"chr1": 5000})
    assert report["overall"]["recall"] is None
    assert report["overall"]["recall_wilson_95"] is None


def test_anchor_read_length_must_match_independent_truth(tmp_path):
    with pytest.raises(ValueError, match="read_length"):
        evaluate_retrieval(
            write_jsonl(tmp_path / "anchors.jsonl", anchor_rows("read", read_length=1000)),
            write_jsonl(tmp_path / "truth.jsonl", [truth_row("read", reference_span=500)]),
            {"chr1": 5000}, positional_slop=0,
        )


def test_zero_query_interval_is_invalid_input_not_a_retrieval_miss(tmp_path):
    anchors = anchor_rows("read")
    anchors[0]["query_end"] = anchors[0]["query_start"]
    with pytest.raises(ValueError, match="anchor|interval"):
        evaluate_retrieval(write_jsonl(tmp_path / "anchors.jsonl", anchors),
                           write_jsonl(tmp_path / "truth.jsonl", [truth_row("read")]),
                           {"chr1": 5000})


@pytest.mark.parametrize("bad", [
    {"read_id": "read", "origins": [None]},
    {"read_id": "read", "origins": [{"chrom": "", "start": 0, "strand": "+"}]},
    {"read_id": "read", "origins": [], "stratum": []},
    {"read_id": "read", "origins": [], "read_length": -1},
    {"read_id": "read", "origins": [], "reference_span": False},
])
def test_malformed_truth_fails_with_a_readable_validation_error(tmp_path, bad):
    with pytest.raises(ValueError):
        load_truth(write_jsonl(tmp_path / "truth.jsonl", [bad]))


@pytest.mark.parametrize("truth", [truth_row("read", chrom="unknown"),
                                  truth_row("read", start=4900, reference_span=200)])
def test_truth_origins_must_fit_the_declared_reference(tmp_path, truth):
    with pytest.raises(ValueError, match="contig|reference"):
        evaluate_retrieval(write_jsonl(tmp_path / "anchors.jsonl", []),
                           write_jsonl(tmp_path / "truth.jsonl", [truth]), {"chr1": 5000})


def test_unknown_read_ids_and_duplicate_primaries_are_errors(tmp_path):
    truth = write_jsonl(tmp_path / "truth.jsonl", [truth_row("read", 100)])
    for name, records in [
        ("unknown", [("other", "chr1", 100, 0, 60)]),
        ("duplicate", [("read", "chr1", 100, 0, 60)] * 2),
    ]:
        with pytest.raises(ValueError):
            evaluate_alignments(write_sam(tmp_path / f"{name}.sam", records), truth)
