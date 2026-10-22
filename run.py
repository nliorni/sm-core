#!/usr/bin/env python3

"""
sm-core suite dispatcher: run any pipeline from the repository root.

    ./run.py --list                                  # the pipelines
    ./run.py rnaseq --list-workflows                 # one pipeline's targets
    ./run.py rnaseq -w all -c my_config.yaml -q 32   # = pipelines/rnaseq/run.py ...

Everything after the pipeline name is passed to that pipeline's own run.py
unchanged, so both spellings are equivalent.
"""

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PIPELINES = ROOT / "pipelines"


def available():
    return sorted(p.name for p in PIPELINES.iterdir() if (p / "run.py").is_file())


def main() -> int:
    args = sys.argv[1:]
    if not args or args[0] in ("-h", "--help", "--list"):
        print(__doc__)
        print("Pipelines:")
        for name in available():
            print(f"  {name}")
        return 0
    name, rest = args[0], args[1:]
    launcher = PIPELINES / name / "run.py"
    if not launcher.is_file():
        sys.exit(f"Unknown pipeline '{name}'. Available: {', '.join(available())}")
    os.execv(sys.executable, [sys.executable, str(launcher), *rest])


if __name__ == "__main__":
    raise SystemExit(main())
