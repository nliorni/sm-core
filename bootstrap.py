#!/usr/bin/env python3

"""
Start a new pipeline of the sm-core suite from templates/pipeline/.

Usage:
    ./bootstrap.py <name> "<one-line description>"

    ./bootstrap.py chipseq "ChIP-seq: FASTQ to peaks"

Copies templates/pipeline/ to pipelines/<name>/, fills in the
{{PIPELINE_NAME}} / {{PIPELINE_DESCRIPTION}} placeholders, and prints what to
edit next. The new pipeline already runs: its toy `example` target works
through the shared launcher (smcore/launcher.py) and helpers
(smcore/common.smk) exactly as every other pipeline of the suite does.
"""

import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TEMPLATE = ROOT / "templates" / "pipeline"
PLACEHOLDERS = ("{{PIPELINE_NAME}}", "{{PIPELINE_DESCRIPTION}}")
TEXT_SUFFIXES = {".py", ".smk", ".yaml", ".yml", ".md", ".txt", ""}


def main() -> int:
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    name, description = sys.argv[1], sys.argv[2]
    if not re.fullmatch(r"[a-z][a-z0-9_-]*", name):
        sys.exit(f"Pipeline name '{name}' must be lowercase letters, digits, '-' or '_'.")

    target = ROOT / "pipelines" / name
    if target.exists():
        sys.exit(f"{target} already exists.")

    shutil.copytree(
        TEMPLATE, target,
        ignore=shutil.ignore_patterns(".snakemake", ".smk_locks", "logs", "__pycache__"),
    )
    for path in target.rglob("*"):
        if not path.is_file() or path.suffix not in TEXT_SUFFIXES:
            continue
        text = path.read_text()
        if any(p in text for p in PLACEHOLDERS):
            text = text.replace(PLACEHOLDERS[0], name).replace(PLACEHOLDERS[1], description)
            path.write_text(text)
    for leftover in (target / "tests" / "tutorial" / "output").glob("*"):
        if leftover.name != ".gitkeep":
            shutil.rmtree(leftover) if leftover.is_dir() else leftover.unlink()

    print(f"Created pipelines/{name}/. It runs already:")
    print(f"  cd pipelines/{name} && ./run.py -w example -c tests/tutorial/config_tutorial.yaml -q 2 --no-conda")
    print()
    print("Next steps (see docs/conventions.md):")
    print("  1. Edit `workflows` in run.py with your real targets.")
    print("  2. Replace workflow/rules/example.smk (+ scripts/, envs/) with real rules.")
    print("  3. Edit config/config.yaml for the project config shape (and add a schema).")
    print("  4. Replace tests/tutorial/ with a real dataset + truth set, and add tests/ci/.")
    print("  5. Add the pipeline to the table in the top-level README.md and to CI.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
