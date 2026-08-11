#!/usr/bin/env python3

"""
Bootstrap a new pipeline from this sm-core template.

Usage:
    ./bootstrap.py <PipelineName>

Replaces the {{PIPELINE_NAME}} placeholder in run.py, README.md, and
PIPELINE.md, then prints a checklist of what to edit next.
"""

import sys
from pathlib import Path

PLACEHOLDER = "{{PIPELINE_NAME}}"
FILES = ["run.py", "README.md", "PIPELINE.md"]


def main() -> int:
    if len(sys.argv) != 2:
        sys.exit("Usage: ./bootstrap.py <PipelineName>")

    name = sys.argv[1]
    root = Path(__file__).resolve().parent

    for fname in FILES:
        path = root / fname
        if not path.exists():
            continue
        text = path.read_text()
        if PLACEHOLDER not in text:
            continue
        path.write_text(text.replace(PLACEHOLDER, name))
        print(f"updated {fname}")

    print()
    print(f"Bootstrapped '{name}'. Next steps (see PIPELINE.md for detail):")
    print("  1. Edit WORKFLOWS in run.py with your real workflow targets.")
    print("  2. Replace workflow/rules/example.smk (+ scripts/, envs/) with your rules.")
    print("  3. Edit config_template.yaml for your project's config shape.")
    print("  4. Update tests/tutorial/ with real (or still-toy) demo data.")
    print("  5. rm -rf .git && git init   # start this pipeline's own history")
    print("  6. Delete this bootstrap.py once you no longer need it (optional).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
