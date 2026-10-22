# Tutorial

A tiny, dependency-free smoke test proving the template runs end-to-end
right after `./bootstrap.py` -- no bioinformatics tools, conda envs, or
downloads required. `data/` holds four small synthetic text files standing in
for paired-end reads; `count_lines`/`summarize` (workflow/rules/example.smk)
count and aggregate them.

Run from this pipeline's directory:

```bash
./run.py -w example -c tests/tutorial/config_tutorial.yaml -q 2 -n --no-conda   # dry-run
./run.py -w example -c tests/tutorial/config_tutorial.yaml -q 2 --no-conda
cat tests/tutorial/output/summary.tsv
```

Expected:

```
sample	r1_lines	r2_lines
SAMPLE_A	4	4
SAMPLE_B	2	2
```

`tests/tutorial/output/` is gitignored -- delete it freely to re-run from
scratch.
