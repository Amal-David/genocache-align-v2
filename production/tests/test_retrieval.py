"""Independent numeric, provenance and publication checks for embedding retrieval."""

import hashlib
import json

import faiss
import numpy as np
import pytest

from genocache import retrieval


REFERENCE_SHA = "a" * 64
ENCODER_SHA = "b" * 64


def write_jsonl(path, records):
    path.write_text("".join(json.dumps(record) + "\n" for record in records))
    return path


@pytest.fixture
def packed(tmp_path):
    vectors = np.array([[3, 0, 0, 0], [0, 4, 0, 0], [0, 0, -5, 0], [1, 2, 3, 4]], dtype=np.float64)
    rows = [
        {"chrom": "chrB", "start": 101, "end": 201, "strand": "+"},
        {"chrom": "chrA", "start": 501, "end": 601, "strand": "-"},
        {"chrom": "chrB", "start": 1001, "end": 1101, "strand": "-"},
        {"chrom": "chrA", "start": 1501, "end": 1601, "strand": "+"},
    ]
    source = tmp_path / "source.npy"
    np.save(source, vectors)
    windows = write_jsonl(tmp_path / "windows.jsonl", rows)
    pack = tmp_path / "pack"
    manifest = retrieval.build_pack(
        source, windows, pack, reference_sha256=REFERENCE_SHA,
        encoder_sha256=ENCODER_SHA, batch_rows=2,
    )
    return {"pack": pack, "vectors": vectors, "windows": rows, "source": source,
            "windows_path": windows, "manifest": manifest}


def resign_manifest(pack, changed_files=(), **updates):
    """Create a self-consistent envelope around bad data to test semantic validation."""
    path = pack / "manifest.json"
    manifest = json.loads(path.read_text())
    manifest.update(updates)
    for name in changed_files:
        manifest["files"][name] = retrieval.file_sha256(pack / name)
    identity = {key: value for key, value in manifest.items()
                if key not in {"pack_id", "build_seconds"}}
    manifest["pack_id"] = hashlib.sha256(
        json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    path.write_text(json.dumps(manifest))


def query_metadata(count):
    return [{"read_id": "query", "read_length": 1000,
             "query_start": i * 100, "query_end": i * 100 + 100,
             "seed_id": str(i), "query_strand": "+"} for i in range(count)]


def search_args(packed, query_path, metadata_path, output, *extra):
    return ["search", "--pack", str(packed["pack"]), "--queries", str(query_path),
            "--metadata", str(metadata_path), "--encoder-sha256", ENCODER_SHA,
            "--output", str(output), "--batch-size", "1", "--k", "2", *extra]


def test_normalization_and_row_identity_round_trip(packed):
    index = retrieval.VectorIndex(packed["pack"])
    expected = packed["vectors"] / np.sqrt(np.sum(packed["vectors"] ** 2, axis=1, keepdims=True))
    np.testing.assert_allclose(index.vectors, expected, rtol=1e-6, atol=1e-7)
    assert index.vectors.dtype == np.dtype("float32")
    scores, ids = index.search(packed["vectors"], k=1)
    np.testing.assert_array_equal(ids[:, 0], np.arange(4))
    anchors = list(index.anchors(scores, ids, query_metadata(4)))
    assert [(a["chrom"], a["ref_start"], a["ref_end"], a["strand"]) for a in anchors] == [
        (row["chrom"], row["start"], row["end"], row["strand"]) for row in packed["windows"]
    ]
    assert packed["manifest"]["source_vectors_sha256"] == retrieval.file_sha256(packed["source"])
    assert packed["manifest"]["source_windows_sha256"] == retrieval.file_sha256(packed["windows_path"])
    assert not packed["manifest"]["authoritative_alignment"]


def test_exact_and_flat_search_agree_with_independent_float64_cosine_oracle(packed):
    index = retrieval.VectorIndex(packed["pack"])
    rng = np.random.default_rng(231)
    queries = rng.normal(size=(19, 4))
    reference = packed["vectors"]
    cosine = (queries @ reference.T) / (
        np.sqrt(np.sum(queries ** 2, axis=1))[:, None]
        * np.sqrt(np.sum(reference ** 2, axis=1))[None, :]
    )
    expected_ids = np.argsort(-cosine, axis=1)[:, :3]
    expected_scores = np.take_along_axis(cosine, expected_ids, axis=1)
    for scores, ids in (index.search(queries, 3), index.exact_search(queries, 3, block_rows=2)):
        np.testing.assert_array_equal(ids, expected_ids)
        np.testing.assert_allclose(scores, expected_scores, atol=2e-6)


@pytest.mark.parametrize("bad", [np.zeros((1, 4)), np.full((1, 4), np.nan),
                                 np.full((1, 4), np.inf), np.ones(4),
                                 np.ones((1, 4), dtype=np.complex64) * (1 + 1j)])
def test_invalid_embeddings_are_rejected(bad):
    with pytest.raises(ValueError):
        retrieval.normalize(bad)


def test_query_dimension_and_empty_result_contracts(packed):
    index = retrieval.VectorIndex(packed["pack"])
    with pytest.raises(ValueError, match="dimension"):
        index.search(np.ones((1, 3)), k=2)
    scores, ids = index.search(np.empty((0, 4), dtype=np.float32), k=2)
    assert scores.shape == ids.shape == (0, 2)
    assert list(index.anchors(scores, ids, [])) == []
    assert index.search(np.ones((1, 4)), k=100)[0].shape == (1, 4)


def test_reverse_complement_orientation_is_composed_without_changing_forward_offsets(packed):
    index = retrieval.VectorIndex(packed["pack"])
    metadata = query_metadata(2)
    metadata[1]["query_strand"] = "-"
    scores = np.ones((2, 2), dtype=np.float32)
    ids = np.array([[0, 1], [0, 1]], dtype=np.int64)
    anchors = list(index.anchors(scores, ids, metadata))
    assert [a["strand"] for a in anchors] == ["+", "-", "-", "+"]
    assert [(a["query_start"], a["query_end"]) for a in anchors] == [(0, 100)] * 2 + [(100, 200)] * 2


@pytest.mark.parametrize("scores,ids", [
    (np.ones((1, 2)), np.array([[0]], dtype=np.int64)),
    (np.ones(1), np.array([0], dtype=np.int64)),
    (np.ones((1, 1)), np.array([[0.5]])),
    (np.ones((1, 1)), np.array([[True]])),
    (np.ones((1, 1)), np.array([[-2]], dtype=np.int64)),
    (np.ones((1, 1)), np.array([[4]], dtype=np.int64)),
    (np.full((1, 1), np.nan), np.array([[0]], dtype=np.int64)),
])
def test_invalid_search_result_shapes_and_row_ids_fail_explicitly(packed, scores, ids):
    index = retrieval.VectorIndex(packed["pack"])
    with pytest.raises(ValueError):
        list(index.anchors(scores, ids, query_metadata(1)))


def test_only_faiss_minus_one_sentinel_is_skipped(packed):
    index = retrieval.VectorIndex(packed["pack"])
    assert list(index.anchors(np.array([[-np.inf]]), np.array([[-1]]), query_metadata(1))) == []


@pytest.mark.parametrize("field,value", [("query_end", 0), ("query_strand", []),
                                         ("read_length", None), ("seed_id", "")])
def test_invalid_query_metadata_fails_explicitly(packed, field, value):
    index = retrieval.VectorIndex(packed["pack"])
    metadata = query_metadata(1)
    metadata[0][field] = value
    with pytest.raises(ValueError):
        list(index.anchors(np.ones((1, 1)), np.array([[0]]), metadata))


@pytest.mark.parametrize("filename", ["vectors.npy", "positions.npy", "search.faiss", "contigs.json"])
def test_corrupted_pack_member_is_rejected_before_use(packed, filename):
    with (packed["pack"] / filename).open("ab") as stream:
        stream.write(b"corrupt")
    with pytest.raises(ValueError, match="checksum"):
        retrieval.VectorIndex(packed["pack"])


@pytest.mark.parametrize("field", ["reference_sha256", "encoder_sha256"])
def test_manifest_provenance_cannot_change_without_invalidating_pack_identity(packed, field):
    manifest_path = packed["pack"] / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest[field] = "c" * 64
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="identity|pack_id"):
        retrieval.VectorIndex(packed["pack"])


@pytest.mark.parametrize("field,value", [("chrom", 5), ("start", 201), ("strand", 0)])
def test_self_consistent_manifest_cannot_hide_invalid_position_rows(packed, field, value):
    positions = np.load(packed["pack"] / "positions.npy", mmap_mode="r+")
    positions[field][0] = value
    positions.flush()
    del positions
    resign_manifest(packed["pack"], ["positions.npy"])
    with pytest.raises(ValueError, match="position|coordinate|strand|contig"):
        retrieval.VectorIndex(packed["pack"])


@pytest.mark.parametrize("replacement", [np.zeros((4, 4), dtype=np.float32),
                                        np.ones((4, 4), dtype=np.float32),
                                        np.full((4, 4), np.nan, dtype=np.float32)])
def test_pack_vectors_must_remain_normalized_and_finite(packed, replacement):
    np.save(packed["pack"] / "vectors.npy", replacement)
    resign_manifest(packed["pack"], ["vectors.npy"])
    with pytest.raises(ValueError, match="vector|normalized|finite"):
        retrieval.VectorIndex(packed["pack"])


def test_pack_contig_dictionary_requires_unique_valid_names(packed):
    (packed["pack"] / "contigs.json").write_text(json.dumps(["chrA", "chrA"]))
    resign_manifest(packed["pack"], ["contigs.json"])
    with pytest.raises(ValueError, match="contig"):
        retrieval.VectorIndex(packed["pack"])


@pytest.mark.parametrize("row_delta", [-1, 1])
def test_window_row_count_failure_does_not_publish_a_partial_pack(tmp_path, row_delta):
    vector_file = tmp_path / "vectors.npy"
    np.save(vector_file, np.eye(3, dtype=np.float32))
    rows = [{"chrom": "chr1", "start": i * 10, "end": i * 10 + 10} for i in range(3 + row_delta)]
    windows = write_jsonl(tmp_path / "windows.jsonl", rows)
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="row"):
        retrieval.build_pack(vector_file, windows, output,
                             reference_sha256=REFERENCE_SHA, encoder_sha256=ENCODER_SHA)
    assert not output.exists()
    assert not list(tmp_path.glob(".embedding-*"))


def test_vector_source_mutation_during_build_invalidates_provenance(packed, tmp_path, monkeypatch):
    normalize = retrieval.normalize
    mutated = False

    def mutate_after_first_batch(batch):
        nonlocal mutated
        result = normalize(batch)
        if not mutated:
            values = np.load(packed["source"], mmap_mode="r+")
            values[0, 1] = 0.7
            values.flush()
            mutated = True
        return result

    monkeypatch.setattr(retrieval, "normalize", mutate_after_first_batch)
    output = tmp_path / "mutated-pack"
    with pytest.raises((ValueError, RuntimeError), match="changed|mutat"):
        retrieval.build_pack(packed["source"], packed["windows_path"], output,
                             reference_sha256=REFERENCE_SHA, encoder_sha256=ENCODER_SHA,
                             batch_rows=2)
    assert not output.exists()


@pytest.mark.parametrize("nlist,pq_m", [(1, 2), (1_000_000, 2), (1, 3)])
def test_invalid_ivfpq_training_configuration_is_rejected_before_native_allocation(
    packed, tmp_path, monkeypatch, nlist, pq_m
):
    def must_not_construct(*args, **kwargs):
        raise AssertionError("Native IVF allocation must follow feasibility checks")

    monkeypatch.setattr(faiss, "IndexIVFPQ", must_not_construct)
    output = tmp_path / "ivf"
    with pytest.raises(ValueError, match="training|divisible"):
        retrieval.build_pack(packed["source"], packed["windows_path"], output,
                             reference_sha256=REFERENCE_SHA, encoder_sha256=ENCODER_SHA,
                             kind="ivfpq", nlist=nlist, pq_m=pq_m)
    assert not output.exists()


@pytest.mark.parametrize("metadata_count", [1, 3])
def test_missing_or_extra_query_metadata_preserves_preexisting_output(packed, tmp_path, metadata_count):
    queries = tmp_path / "queries.npy"
    np.save(queries, packed["vectors"][:2])
    metadata = write_jsonl(tmp_path / "queries.jsonl", query_metadata(metadata_count))
    output = tmp_path / "hits.jsonl"
    output.write_text("previous complete output\n")
    with pytest.raises(SystemExit) as raised:
        retrieval.main(search_args(packed, queries, metadata, output))
    assert raised.value.code == 2
    assert output.read_text() == "previous complete output\n"
    assert not list(tmp_path.glob(".hits-*"))


def test_empty_exact_audit_does_not_publish_output_before_failure(packed, tmp_path):
    queries = tmp_path / "queries.npy"
    np.save(queries, np.empty((0, 4), dtype=np.float32))
    metadata = write_jsonl(tmp_path / "queries.jsonl", [])
    output = tmp_path / "hits.jsonl"
    output.write_text("previous complete output\n")
    with pytest.raises(SystemExit) as raised:
        retrieval.main(search_args(packed, queries, metadata, output, "--audit-exact", "1"))
    assert raised.value.code == 2
    assert output.read_text() == "previous complete output\n"


def test_failed_optional_exact_audit_does_not_publish_output(packed, tmp_path, monkeypatch):
    queries = tmp_path / "queries.npy"
    np.save(queries, packed["vectors"][:2])
    metadata = write_jsonl(tmp_path / "queries.jsonl", query_metadata(2))
    output = tmp_path / "hits.jsonl"
    output.write_text("previous complete output\n")

    def fail_oracle(*args, **kwargs):
        raise RuntimeError("controlled exact oracle failure")

    monkeypatch.setattr(retrieval.VectorIndex, "exact_search", fail_oracle)
    with pytest.raises(SystemExit) as raised:
        retrieval.main(search_args(packed, queries, metadata, output, "--audit-exact", "1"))
    assert raised.value.code == 2
    assert output.read_text() == "previous complete output\n"


def test_successful_batched_search_publishes_all_hits_and_separate_exact_audit(packed, tmp_path, capsys):
    queries = tmp_path / "queries.npy"
    np.save(queries, np.array([[3, 1, 0.2, 0.1], [0.2, 3, 1, 0.1]]))
    metadata = write_jsonl(tmp_path / "queries.jsonl", query_metadata(2))
    output = tmp_path / "hits.jsonl"
    assert retrieval.main(search_args(packed, queries, metadata, output, "--audit-exact", "2")) == 0
    report = json.loads(capsys.readouterr().out)
    assert len(list(retrieval.json_lines(output))) == report["hits"] == 4
    assert report["exact_vector_audit"]["top_k_overlap"] == 1
    assert report["exact_vector_audit"]["biological_locus_recall"] is None
    assert not report["authoritative_alignment"]
