#!/usr/bin/env python3
"""Reproducible synthetic integration benchmark, not a human-genome claim."""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import platform
import random
import subprocess
import time

import pysam

from genocache import Engine
from genocache.evaluate import evaluate_alignments
from genocache.retrieval import file_sha256


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def make_panel(directory: Path):
    directory.mkdir(parents=True, exist_ok=True)
    rng = random.Random(20261007)
    def dna(length):
        return "".join(rng.choices("ACGT", k=length))
    references = {"chr1": dna(400000), "chr2": dna(400000)}
    repeat = dna(12000)
    copies = [("chr1", 100000), ("chr2", 150000), ("chr2", 300000)]
    for chrom, start in copies:
        ref = references[chrom]
        references[chrom] = ref[:start] + repeat + ref[start + len(repeat):]
    def mutate(sequence, substitution, insertion, deletion):
        out = []
        for base in sequence:
            if rng.random() < deletion:
                continue
            if rng.random() < substitution:
                base = rng.choice([x for x in "ACGT" if x != base])
            out.append(base)
            if rng.random() < insertion:
                out.append(rng.choice("ACGT"))
        return "".join(out)
    records = []
    truth = []
    for stratum in ("unique_exact", "substitutions_5pct", "mixed_indels", "repeat_copies"):
        for i in range(32):
            length = (1200, 2400, 4800, 8000)[i % 4]
            chrom = "chr1" if i % 2 else "chr2"
            if stratum == "repeat_copies":
                length = min(length, 6000)
                offset = 1000 + (i * 131) % (10000 - length)
                chrom, start = copies[i % len(copies)]
                start += offset
                origins = [{"chrom": c, "start": p + offset} for c, p in copies]
            else:
                while True:
                    start = rng.randrange(1000, 390000 - length)
                    if not any(c == chrom and start < p + 12000 and start + length > p for c, p in copies):
                        break
                origins = [{"chrom": chrom, "start": start}]
            read = references[chrom][start:start + length]
            if stratum == "substitutions_5pct":
                read = mutate(read, 0.05, 0, 0)
            elif stratum == "mixed_indels":
                read = mutate(read, 0.03, 0.025, 0.025)
            strand = "-" if i % 2 else "+"
            if strand == "-":
                read = reverse_complement(read)
            name = f"{stratum}_{i:03d}"
            records.append((name, read))
            truth.append({
                "read_id": name, "read_length": len(read), "reference_span": length,
                "stratum": stratum,
                "origins": [{**origin, "strand": strand} for origin in origins],
            })
    for i in range(8):
        name, read = f"unmapped_{i:03d}", dna(1800)
        records.append((name, read))
        truth.append({
            "read_id": name, "read_length": len(read),
            "stratum": "known_unmappable", "origins": [],
        })
    reference = directory / "reference.fa"
    reference.write_text("".join(f">{name}\n{seq}\n" for name, seq in references.items()))
    content = "".join(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n" for name, seq in records).encode()
    reads = directory / "reads.fastq.gz"
    reads.write_bytes(gzip.compress(content, mtime=0))
    truth_path = directory / "truth.jsonl"
    truth_path.write_text("".join(json.dumps(item) + "\n" for item in truth))
    (directory / "contigs.json").write_text(json.dumps({c: len(s) for c, s in references.items()}))
    return reference, reads, truth_path


def record_multiset(path):
    with pysam.AlignmentFile(str(path)) as stream:
        return Counter(record.to_string() for record in stream.fetch(until_eof=True))


def run(directory: Path, binary: str):
    if directory.exists() and any(directory.iterdir()):
        raise ValueError("choose a new empty work directory to measure cold execution honestly")
    directory.mkdir(parents=True, exist_ok=True)
    reference, reads, truth = make_panel(directory / "panel")
    engine = Engine(directory / "cache", binary)
    pack = engine.build_index(reference, threads=2)
    pack_warm = engine.build_index(reference, threads=2)
    cold = engine.align(pack["reference_id"], reads, threads=2, output_format="bam")
    warm = engine.align(pack["reference_id"], reads, threads=2, output_format="bam")
    native_sam = directory / "native.sam"
    native_bam = directory / "native.bam"
    phase = time.perf_counter()
    with native_sam.open("wb") as out:
        subprocess.run(
            [binary, "-a", "-x", "map-ont", "-t", "2", pack["index_path"], str(reads)],
            stdout=out, stderr=subprocess.PIPE, check=True,
        )
    native_map = time.perf_counter() - phase
    phase = time.perf_counter()
    pysam.sort("-@", "2", "-m", "256M", "-o", str(native_bam), str(native_sam))
    pysam.index("-c", str(native_bam))
    native_sort = time.perf_counter() - phase
    same_records = record_multiset(native_bam) == record_multiset(Path(cold["output_path"]))
    if not same_records:
        raise AssertionError("native BAM and GenoCache BAM record multisets differ")
    evaluation = evaluate_alignments(Path(cold["output_path"]), truth, tolerance=100)
    result = {
        "schema": "genocache.local-integration-benchmark.v1",
        "scope": "synthetic 800 kb reference; not a real human genome or a learned-model benchmark",
        "random_seed": 20261007,
        "platform": platform.platform(),
        "python": platform.python_version(),
        "minimap2_version": subprocess.check_output([binary, "--version"], text=True).strip(),
        "minimap2_sha256": file_sha256(Path(binary)),
        "pysam_version": pysam.__version__,
        "inputs": {
            "reference_sha256": file_sha256(reference),
            "reads_sha256": file_sha256(reads),
            "truth_sha256": file_sha256(truth),
            **cold["input_stats"],
        },
        "threads": 2,
        "native_equivalent_output": {
            "mapping_seconds": native_map,
            "sort_and_index_seconds": native_sort,
            "total_seconds": native_map + native_sort,
            "record_multiset_identical": same_records,
            "timing_note": "Already-built .mmi; no GenoCache input snapshot/checksum/validation overhead.",
        },
        "genocache_index_cold": pack["metrics"],
        "genocache_index_warm": pack_warm["metrics"],
        "genocache_job_cold": cold["metrics"],
        "genocache_job_warm": warm["metrics"],
        "exact_repeat_cache_hit": warm["cache_hit"],
        "same_job_id": cold["job_id"] == warm["job_id"],
        "warm_artifact_sha256_matches": file_sha256(Path(warm["output_path"])) == file_sha256(Path(cold["output_path"])),
        "mapping_evaluation": {k: v for k, v in evaluation.items() if k != "per_read"},
        "interpretation": {
            "algorithmic_speedup_claim": False,
            "cold_ratio_native_over_genocache": (native_map + native_sort) / cold["metrics"]["request_wall_seconds"],
            "exact_repeat_ratio_cold_over_warm": cold["metrics"]["request_wall_seconds"] / warm["metrics"]["request_wall_seconds"],
            "caution": "Warm hit replays the identical input batch; it does not measure novel-read throughput or real-workload hit rate.",
        },
    }
    (directory / "per-read-evaluation.json").write_text(json.dumps(evaluation, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--minimap2", required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.workdir, args.minimap2)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
