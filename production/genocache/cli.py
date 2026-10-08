"""JSON-only stdout command-line interface; no read sequences are printed."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

from . import __version__
from .engine import Engine, SUPPORTED_PRESETS
from .io import EngineError


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Verified minimap2 batch alignment with exact-input caching")
    parser.add_argument("--version", action="version", version=f"genocache {__version__}")
    parser.add_argument("--root", type=Path, default=Path(os.environ.get("GENOCACHE_ROOT", "./genocache-data")),
                        help="Local worker cache directory; never a shared/network volume")
    parser.add_argument("--minimap2", default=os.environ.get("MINIMAP2", "minimap2"), help="minimap2 executable")
    commands = parser.add_subparsers(dest="command", required=True)
    index = commands.add_parser("index", help="Build or verify an immutable reference pack")
    index.add_argument("--reference", type=Path, required=True)
    index.add_argument("--preset", choices=sorted(SUPPORTED_PRESETS), default="map-ont")
    index.add_argument("--threads", type=int, default=1)
    index.add_argument("--timeout", type=float)
    align = commands.add_parser("align", help="Align one FASTA/FASTQ batch or verify an exact cached job")
    align.add_argument("--reference-id", required=True)
    align.add_argument("--reads", type=Path, required=True)
    align.add_argument("--threads", type=int, default=1)
    align.add_argument("--timeout", type=float)
    align.add_argument("--output-format", choices=["bam", "sam"], default="bam")
    inspect = commands.add_parser("inspect", help="Verify and inspect a reference pack")
    inspect.add_argument("--reference-id", required=True)
    commands.add_parser("doctor", help="Check native binary, dependencies, and local engine directory")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        engine = Engine(args.root, minimap2=args.minimap2)
        if args.command == "index":
            result = engine.build_index(args.reference, preset=args.preset, threads=args.threads, timeout=args.timeout)
        elif args.command == "align":
            result = engine.align(args.reference_id, args.reads, threads=args.threads,
                                  timeout=args.timeout, output_format=args.output_format)
        elif args.command == "inspect":
            result = engine.inspect_reference(args.reference_id)
        else:
            binary, identity = engine._binary_identity()
            result = {"ok": True, "engine_version": __version__, "root": str(engine.root),
                      "minimap2": {"path": str(binary), **identity},
                      "dependencies": engine._dependency_identity(),
                      "presets": sorted(SUPPORTED_PRESETS), "index_residency": "disk; loaded once per cold batch"}
        print(json.dumps(result, sort_keys=True, allow_nan=False))
        return 0
    except EngineError as exc:
        print(json.dumps({"ok": False, "error": type(exc).__name__, "message": str(exc)}, sort_keys=True), file=sys.stderr)
        return 2
    except OSError as exc:
        print(json.dumps({"ok": False, "error": "FilesystemError", "message": exc.strerror or "Filesystem operation failed"}), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
