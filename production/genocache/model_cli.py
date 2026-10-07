"""JSON interface to the mathematical model; missing observations stay unknown."""
import argparse
import json
from pathlib import Path

from .model import analyze_model


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    try:
        report = analyze_model(json.loads(args.config.read_text()))
        text = json.dumps(report, indent=2, allow_nan=False) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(text)
        print(text, end="")
        return 0
    except (OSError, ValueError) as exc:
        parser.exit(2, f"genocache-model: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
