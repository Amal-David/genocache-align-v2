"""Batched embedding retrieval for experiments, never an alignment-output authority.

An encoder is an explicit input artifact. This module does not silently replace
missing learned weights with random or handcrafted vectors. Reference positions
are columnar, search uses normalized inner product, and query orientations travel
with every hit. Dense vectors remain on disk for an exact-vector recall audit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
import time
from typing import Iterator

import numpy as np

from .io import file_signature


POSITION_DTYPE = np.dtype(
    [("chrom", "<u4"), ("start", "<u8"), ("end", "<u8"), ("strand", "i1")]
)
SCHEMA = "genocache.embedding-pack.v1"
PACK_FILES = {"vectors.npy", "positions.npy", "search.faiss", "contigs.json"}
MAX_COORDINATE = (1 << 63) - 1


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_lines(path: Path) -> Iterator[dict]:
    with Path(path).open(encoding="utf-8") as stream:
        for number, line in enumerate(stream, 1):
            if line.strip():
                record = json.loads(line)
                if not isinstance(record, dict):
                    raise ValueError(f"JSONL row {number} must be an object")
                yield record


def normalize(vectors: np.ndarray) -> np.ndarray:
    raw = np.asarray(vectors)
    if raw.ndim != 2 or raw.dtype.kind not in "iuf":
        raise ValueError("embeddings must be a real numeric two-dimensional array")
    with np.errstate(over="ignore", invalid="ignore"):
        values = np.asarray(raw, dtype=np.float32, order="C")
    if not np.isfinite(values).all():
        raise ValueError("embeddings must be a finite two-dimensional array")
    with np.errstate(over="ignore", invalid="ignore"):
        norm = np.linalg.norm(values, axis=1, keepdims=True)
    if np.any(norm <= 1e-12) or not np.isfinite(norm).all():
        raise ValueError("zero or overflowing embedding vectors are invalid")
    return np.ascontiguousarray(values / norm, dtype=np.float32)


def _sha(value: str, name: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")
    return value


def _integer(value, name: str, minimum: int = 0, maximum: int = MAX_COORDINATE) -> int:
    if (
        isinstance(value, bool) or not isinstance(value, (int, np.integer))
        or not minimum <= value <= maximum
    ):
        raise ValueError(f"{name} must be an integer between {minimum} and {maximum}")
    return int(value)


def _pack_identity(manifest: dict) -> str:
    payload = {key: value for key, value in manifest.items()
               if key not in {"pack_id", "build_seconds"}}
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def _validate_manifest(manifest: dict) -> None:
    if not isinstance(manifest, dict) or manifest.get("schema") != SCHEMA:
        raise ValueError("unsupported embedding-pack schema")
    if not isinstance(manifest.get("files"), dict) or set(manifest["files"]) != PACK_FILES:
        raise ValueError("invalid embedding-pack file inventory")
    for field in ("pack_id", "reference_sha256", "encoder_sha256",
                  "source_vectors_sha256", "source_windows_sha256"):
        _sha(manifest.get(field), field)
    for name, digest in manifest["files"].items():
        _sha(digest, name)
    if _pack_identity(manifest) != manifest["pack_id"]:
        raise ValueError("embedding-pack manifest identity does not match pack_id")
    _integer(manifest.get("rows"), "pack rows", 1)
    _integer(manifest.get("dimension"), "pack dimension", 1, (1 << 31) - 1)
    if manifest.get("authoritative_alignment") is not False:
        raise ValueError("embedding packs cannot claim authoritative alignments")
    config = manifest.get("config")
    if (
        not isinstance(config, dict) or config.get("kind") not in ("flat", "ivfpq")
        or config.get("metric") != "normalized_inner_product"
    ):
        raise ValueError("unsupported embedding-pack index configuration")
    if config["kind"] == "ivfpq":
        _integer(config.get("nlist"), "nlist", 1, (1 << 31) - 1)
        _integer(config.get("pq_m"), "pq_m", 1, manifest["dimension"])
        if config.get("pq_bits") != 8 or manifest["dimension"] % config["pq_m"]:
            raise ValueError("invalid embedding-pack product quantization configuration")


def _assert_sources_unchanged(signatures: dict[Path, tuple], stage: str) -> None:
    if any(file_signature(path) != before for path, before in signatures.items()):
        raise ValueError(f"embedding input changed during {stage}; retry with immutable input files")


def build_pack(
    vectors: Path,
    windows: Path,
    output: Path,
    *,
    reference_sha256: str,
    encoder_sha256: str,
    kind: str = "flat",
    nlist: int = 1024,
    pq_m: int = 16,
    batch_rows: int = 8192,
    seed: int = 17,
) -> dict:
    """Import an existing encoder's vectors without materializing the whole matrix.

    windows JSONL has chrom,start,end,strand, in exactly the vector row order.
    Coordinates are zero-based, half-open and always on the forward reference.
    strand identifies which sequence orientation the encoder saw.
    """
    import faiss

    reference_sha256 = _sha(reference_sha256, "reference_sha256")
    encoder_sha256 = _sha(encoder_sha256, "encoder_sha256")
    if kind not in {"flat", "ivfpq"}:
        raise ValueError("kind must be flat or ivfpq")
    batch_rows = _integer(batch_rows, "batch_rows", 1)
    vectors, windows = Path(vectors), Path(windows)
    source_signatures = {path: file_signature(path) for path in (vectors, windows)}
    values = np.load(vectors, mmap_mode="r", allow_pickle=False)
    if values.ndim != 2 or not all(values.shape):
        raise ValueError("vectors must be a nonempty N x D .npy array")
    n, dim = values.shape
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".embedding-", dir=output.parent))
    started = time.perf_counter()
    try:
        normalized = np.lib.format.open_memmap(
            temporary / "vectors.npy", mode="w+", dtype="<f4", shape=(n, dim)
        )
        for offset in range(0, n, batch_rows):
            normalized[offset : offset + batch_rows] = normalize(values[offset : offset + batch_rows])
        normalized.flush()
        positions = np.lib.format.open_memmap(
            temporary / "positions.npy", mode="w+", dtype=POSITION_DTYPE, shape=(n,)
        )
        contigs: list[str] = []
        contig_ids: dict[str, int] = {}
        count = 0
        for row, item in enumerate(json_lines(windows)):
            if row >= n:
                raise ValueError("more window rows than embedding vectors")
            chrom = item.get("chrom")
            if not isinstance(chrom, str) or not chrom or any(c.isspace() for c in chrom):
                raise ValueError("chrom must be a nonempty name without whitespace")
            start = _integer(item.get("start"), "start")
            end = _integer(item.get("end"), "end", 1)
            if end <= start:
                raise ValueError("window end must be greater than start")
            strand = item.get("strand", "+")
            if not isinstance(strand, str) or strand not in {"+", "-"}:
                raise ValueError("window strand must be + or -")
            if chrom not in contig_ids:
                contig_ids[chrom] = len(contigs)
                contigs.append(chrom)
            positions[row] = (contig_ids[chrom], start, end, 1 if strand == "+" else -1)
            count += 1
        if count != n:
            raise ValueError("window row count differs from embedding vector count")
        positions.flush()
        if kind == "flat":
            index = faiss.IndexFlatIP(dim)
            config = {"kind": kind, "metric": "normalized_inner_product"}
        else:
            nlist = _integer(nlist, "nlist", 1)
            pq_m = _integer(pq_m, "pq_m", 1)
            seed = _integer(seed, "training seed", 0, (1 << 31) - 1)
            if dim % pq_m:
                raise ValueError("embedding dimension must be divisible by pq_m")
            sample_size = min(n, 262144)
            if sample_size < max(39 * nlist, 39 * 256):
                raise ValueError("insufficient representative training rows for this IVF-PQ config")
            rng = np.random.default_rng(seed)
            sample_ids = np.sort(rng.choice(n, sample_size, replace=False))
            training = np.ascontiguousarray(normalized[sample_ids], dtype=np.float32)
            quantizer = faiss.IndexFlatIP(dim)
            index = faiss.IndexIVFPQ(quantizer, dim, nlist, pq_m, 8, faiss.METRIC_INNER_PRODUCT)
            index.cp.seed = seed
            index.pq.cp.seed = seed
            index.train(training)
            del training
            config = {
                "kind": kind,
                "metric": "normalized_inner_product",
                "nlist": nlist,
                "pq_m": pq_m,
                "pq_bits": 8,
                "training_rows": sample_size,
                "training_seed": seed,
            }
        for offset in range(0, n, batch_rows):
            index.add(np.ascontiguousarray(normalized[offset : offset + batch_rows]))
        faiss.write_index(index, str(temporary / "search.faiss"))
        (temporary / "contigs.json").write_text(json.dumps(contigs), encoding="utf-8")
        del normalized, positions
        source_vectors_sha256 = file_sha256(vectors)
        source_windows_sha256 = file_sha256(windows)
        _assert_sources_unchanged(source_signatures, "index construction")
        manifest = {
            "schema": SCHEMA,
            "rows": n,
            "dimension": dim,
            "reference_sha256": reference_sha256,
            "encoder_sha256": encoder_sha256,
            "source_vectors_sha256": source_vectors_sha256,
            "source_windows_sha256": source_windows_sha256,
            "faiss_version": faiss.__version__,
            "config": config,
            "files": {
                name: file_sha256(temporary / name)
                for name in ("vectors.npy", "positions.npy", "search.faiss", "contigs.json")
            },
            "authoritative_alignment": False,
        }
        manifest["pack_id"] = _pack_identity(manifest)
        manifest["build_seconds"] = time.perf_counter() - started
        (temporary / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        os.rename(temporary, output)
        return manifest
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


class VectorIndex:
    """A loaded index reused for all batches in an experiment."""

    def __init__(self, pack: Path, *, device: str = "cpu", verify: bool = True):
        import faiss

        self.path = Path(pack)
        self.manifest = json.loads((self.path / "manifest.json").read_text())
        m = self.manifest
        _validate_manifest(m)
        if verify:
            for name, digest in m["files"].items():
                if file_sha256(self.path / name) != digest:
                    raise ValueError(f"embedding-pack checksum mismatch: {name}")
        self.positions = np.load(self.path / "positions.npy", mmap_mode="r", allow_pickle=False)
        self.vectors = np.load(self.path / "vectors.npy", mmap_mode="r", allow_pickle=False)
        self.contigs = json.loads((self.path / "contigs.json").read_text())
        if (
            not isinstance(self.contigs, list) or not self.contigs
            or any(not isinstance(name, str) or not name or any(c.isspace() for c in name)
                   for name in self.contigs)
            or len(set(self.contigs)) != len(self.contigs)
        ):
            raise ValueError("invalid or duplicate embedding-pack contig names")
        self.index = faiss.read_index(str(self.path / "search.faiss"))
        if (
            self.positions.shape != (m["rows"],)
            or self.positions.dtype != POSITION_DTYPE
            or self.vectors.shape != (m["rows"], m["dimension"])
            or self.vectors.dtype != np.dtype("<f4")
            or self.index.ntotal != m["rows"]
            or self.index.d != m["dimension"]
            or self.index.metric_type != faiss.METRIC_INNER_PRODUCT
        ):
            raise ValueError("embedding-pack shape/metric mismatch")
        if m["config"]["kind"] == "flat":
            if not isinstance(self.index, faiss.IndexFlat):
                raise ValueError("embedding-pack index kind differs from manifest")
        elif (
            not isinstance(self.index, faiss.IndexIVFPQ)
            or self.index.nlist != m["config"]["nlist"]
            or self.index.pq.M != m["config"]["pq_m"]
            or self.index.pq.nbits != m["config"]["pq_bits"]
            or not self.index.is_trained
        ):
            raise ValueError("embedding-pack IVF-PQ configuration mismatch")
        for offset in range(0, m["rows"], 8192):
            positions = self.positions[offset : offset + 8192]
            if (
                np.any(positions["chrom"] >= len(self.contigs))
                or np.any(positions["end"] <= positions["start"])
                or np.any(positions["end"] > MAX_COORDINATE)
                or np.any((positions["strand"] != 1) & (positions["strand"] != -1))
            ):
                raise ValueError("embedding-pack position coordinate/strand/contig mismatch")
            if verify:
                vectors = self.vectors[offset : offset + 8192]
                with np.errstate(over="ignore", invalid="ignore"):
                    norm = np.linalg.norm(vectors, axis=1)
                if not np.isfinite(vectors).all() or not np.allclose(norm, 1, rtol=1e-5, atol=1e-6):
                    raise ValueError("embedding-pack vectors must be finite and normalized")
        if device == "cuda":
            if not hasattr(faiss, "StandardGpuResources"):
                raise RuntimeError("CUDA requested but installed FAISS has no GPU support")
            self.resources = faiss.StandardGpuResources()
            self.index = faiss.index_cpu_to_gpu(self.resources, 0, self.index)
        elif device != "cpu":
            raise ValueError("device must be cpu or cuda")
        self.device = device

    def search(self, queries: np.ndarray, k: int = 32, nprobe: int = 16):
        queries = normalize(queries)
        if queries.shape[1] != self.manifest["dimension"]:
            raise ValueError("query embedding dimension does not match reference encoder")
        k = min(_integer(k, "k", 1), self.manifest["rows"])
        nprobe = _integer(nprobe, "nprobe", 1)
        if self.manifest["config"]["kind"] == "ivfpq":
            self.index.nprobe = min(nprobe, self.manifest["config"]["nlist"])
        return self.index.search(queries, k)

    def exact_search(self, queries: np.ndarray, k: int, block_rows: int = 16384):
        """Blocked dense-vector oracle; explicitly unsuitable as a genome-scale default."""
        q = normalize(queries)
        if q.shape[1] != self.manifest["dimension"]:
            raise ValueError("query dimension mismatch")
        block_rows = _integer(block_rows, "block_rows", 1)
        k = min(_integer(k, "k", 1), self.manifest["rows"])
        best_scores = np.full((len(q), 0), -np.inf, dtype=np.float32)
        best_ids = np.empty((len(q), 0), dtype=np.int64)
        for start in range(0, len(self.vectors), block_rows):
            end = min(start + block_rows, len(self.vectors))
            scores = q @ self.vectors[start:end].T
            ids = np.broadcast_to(np.arange(start, end), scores.shape)
            scores = np.concatenate((best_scores, scores), axis=1)
            ids = np.concatenate((best_ids, ids), axis=1)
            keep = min(k, scores.shape[1])
            selected = np.argpartition(-scores, keep - 1, axis=1)[:, :keep]
            best_scores = np.take_along_axis(scores, selected, axis=1)
            best_ids = np.take_along_axis(ids, selected, axis=1)
        order = np.argsort(-best_scores, axis=1, kind="stable")
        return np.take_along_axis(best_scores, order, 1), np.take_along_axis(best_ids, order, 1)

    def anchors(self, scores, ids, query_rows: list[dict]) -> Iterator[dict]:
        """Adapt hits while retaining original-forward-read query coordinates."""
        scores, ids = np.asarray(scores), np.asarray(ids)
        if (
            scores.ndim != 2 or ids.ndim != 2 or scores.shape != ids.shape
            or scores.dtype.kind not in "iuf" or ids.dtype.kind not in "iu"
        ):
            raise ValueError("search scores/IDs must have matching two-dimensional numeric shapes and integer IDs")
        if len(query_rows) != len(ids):
            raise ValueError("query metadata count does not match search result count")
        for metadata, distances, neighbors in zip(query_rows, scores, ids, strict=True):
            if not isinstance(metadata, dict):
                raise ValueError("query metadata rows must be objects")
            read_length = _integer(metadata.get("read_length"), "read_length", 1)
            qs = _integer(metadata.get("query_start"), "query_start")
            qe = _integer(metadata.get("query_end"), "query_end", 1)
            if not qs < qe <= read_length:
                raise ValueError("invalid query interval")
            orientation = metadata.get("query_strand", "+")
            if not isinstance(orientation, str) or orientation not in {"+", "-"}:
                raise ValueError("query_strand must be + or -")
            read_id, seed_id = metadata.get("read_id"), metadata.get("seed_id")
            if (
                not isinstance(read_id, str) or not read_id.strip()
                or not isinstance(seed_id, str) or not seed_id.strip()
            ):
                raise ValueError("read_id and seed_id must be nonempty strings")
            for score, row in zip(distances, neighbors, strict=True):
                if row == -1:
                    continue  # FAISS uses -1 when no neighbor was found.
                if row < -1 or row >= len(self.positions) or not math.isfinite(float(score)):
                    raise ValueError("invalid FAISS result")
                item = self.positions[row]
                strand = (1 if orientation == "+" else -1) * int(item["strand"])
                yield {
                    "read_id": read_id,
                    "read_length": read_length,
                    "query_start": qs,
                    "query_end": qe,
                    "seed_id": seed_id,
                    "chrom": self.contigs[int(item["chrom"])],
                    "ref_start": int(item["start"]),
                    "ref_end": int(item["end"]),
                    "strand": "+" if strand > 0 else "-",
                    "similarity": float(score),
                }


def _search_pack(args) -> dict:
    _integer(args.batch_size, "batch_size", 1)
    _integer(args.k, "k", 1)
    _integer(args.nprobe, "nprobe", 1)
    _integer(args.audit_exact, "audit_exact", 0, 128)
    started = time.perf_counter()
    index = VectorIndex(args.pack, device=args.device)
    load_seconds = time.perf_counter() - started
    if args.encoder_sha256 != index.manifest["encoder_sha256"]:
        raise ValueError("query encoder SHA-256 differs from reference encoder")
    source_signatures = {path: file_signature(path) for path in (args.queries, args.metadata)}
    output_path = args.output.resolve()
    protected = [*source_signatures, *(args.pack / name for name in PACK_FILES | {"manifest.json"})]
    if output_path in {path.resolve() for path in protected}:
        raise ValueError("search output must not overwrite query inputs or embedding-pack members")
    queries = np.load(args.queries, mmap_mode="r", allow_pickle=False)
    if queries.ndim != 2 or queries.shape[1] != index.manifest["dimension"]:
        raise ValueError("query array shape mismatch")
    normalize(queries[:0])  # Validate numeric dtype even when the query set is empty.
    if args.audit_exact and not len(queries):
        raise ValueError("cannot audit an empty query set")
    rows = iter(json_lines(args.metadata))
    batch_seconds = 0.0
    hit_count = 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=".hits-", dir=args.output.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as out:
            for start in range(0, len(queries), args.batch_size):
                q = queries[start : start + args.batch_size]
                metadata = []
                for _ in range(len(q)):
                    try:
                        metadata.append(next(rows))
                    except StopIteration as exc:
                        raise ValueError("missing query metadata rows") from exc
                t = time.perf_counter()
                scores, ids = index.search(q, args.k, args.nprobe)
                batch_seconds += time.perf_counter() - t
                for anchor in index.anchors(scores, ids, metadata):
                    out.write(json.dumps(anchor, separators=(",", ":"), allow_nan=False) + "\n")
                    hit_count += 1
            if next(rows, None) is not None:
                raise ValueError("extra query metadata rows")
            out.flush()
            os.fsync(out.fileno())
        result = {
            "pack_id": index.manifest["pack_id"],
            "query_vectors_sha256": file_sha256(args.queries),
            "query_metadata_sha256": file_sha256(args.metadata),
            "queries": len(queries), "hits": hit_count, "device": args.device,
            "index_load_and_verify_seconds": load_seconds,
            "normalize_and_search_seconds": batch_seconds,
            "encoding_seconds": None,
            "authoritative_alignment": False,
        }
        if args.audit_exact:
            count = min(args.audit_exact, len(queries))
            selected = np.linspace(0, len(queries) - 1, count, dtype=int)
            subset = queries[selected]
            t = time.perf_counter()
            _, approx = index.search(subset, args.k, args.nprobe)
            _, exact = index.exact_search(subset, args.k)
            overlap = sum(len(set(a) & set(b)) for a, b in zip(approx, exact, strict=True))
            result["exact_vector_audit"] = {
                "queries": count, "top_k_overlap": overlap / exact.size,
                "seconds": time.perf_counter() - t,
                "biological_locus_recall": None,
                "ties_note": "Equal-score vector ties can reduce ID-overlap despite equivalent scores.",
            }
        _assert_sources_unchanged(source_signatures, "retrieval")
        result["wall_seconds"] = time.perf_counter() - started
        json.dumps(result, allow_nan=False)
        # Publish only after metadata, optional audit and input provenance all pass.
        os.replace(temporary_name, args.output)
        return result
    finally:
        rows.close()
        Path(temporary_name).unlink(missing_ok=True)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="command", required=True)
    build = subs.add_parser("build")
    build.add_argument("--vectors", type=Path, required=True)
    build.add_argument("--windows", type=Path, required=True)
    build.add_argument("--output", type=Path, required=True)
    build.add_argument("--reference-sha256", required=True)
    build.add_argument("--encoder-sha256", required=True)
    build.add_argument("--kind", choices=["flat", "ivfpq"], default="flat")
    build.add_argument("--nlist", type=int, default=1024)
    build.add_argument("--pq-m", type=int, default=16)
    search = subs.add_parser("search")
    search.add_argument("--pack", type=Path, required=True)
    search.add_argument("--queries", type=Path, required=True)
    search.add_argument("--metadata", type=Path, required=True)
    search.add_argument("--encoder-sha256", required=True)
    search.add_argument("--output", type=Path, required=True)
    search.add_argument("--k", type=int, default=32)
    search.add_argument("--nprobe", type=int, default=16)
    search.add_argument("--batch-size", type=int, default=256)
    search.add_argument("--device", choices=["cpu", "cuda"], default="cpu")
    search.add_argument("--audit-exact", type=int, default=0, help="Dense oracle query sample; 0 disables")
    args = parser.parse_args(argv)
    try:
        if args.command == "build":
            result = build_pack(
                args.vectors, args.windows, args.output,
                reference_sha256=args.reference_sha256,
                encoder_sha256=args.encoder_sha256,
                kind=args.kind, nlist=args.nlist, pq_m=args.pq_m,
            )
        else:
            result = _search_pack(args)
        print(json.dumps(result, indent=2, allow_nan=False))
        return 0
    except (ValueError, OSError, RuntimeError, ImportError) as exc:
        parser.exit(2, f"genocache-retrieval: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
