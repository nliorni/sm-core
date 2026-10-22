#!/usr/bin/env python3

"""
Launcher for sm-core/{{PIPELINE_NAME}}.

    ./run.py --list-workflows
    ./run.py -w example -c tests/tutorial/config_tutorial.yaml -q 2 -n
    ./run.py -w example -c tests/tutorial/config_tutorial.yaml -q 2

All flags: ./run.py --help (shared by every sm-core pipeline; see
smcore/launcher.py).
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))  # the sm-core repository root

from smcore.launcher import Pipeline, main  # noqa: E402

PIPELINE = Pipeline(
    name="{{PIPELINE_NAME}}",
    version="0.1.0",
    description="{{PIPELINE_DESCRIPTION}}",
    directory=HERE,
    # >>> EDIT PER PIPELINE <<<
    # Workflow target -> one-line description. Every key must be a rule in
    # the Snakefile. "example" is the toy demo shipped with the template.
    workflows={
        "example": "Toy demo rule chain (see tests/tutorial/) -- replace with real targets.",
    },
)

if __name__ == "__main__":
    raise SystemExit(main(PIPELINE))
