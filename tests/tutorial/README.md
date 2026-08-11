# Tutorial

A tiny, dependency-free smoke test proving the template runs end-to-end
right after cloning -- no bioinformatics tools, conda envs, or downloads
required. `data/` holds four small synthetic text files standing in for
paired-end reads; `count_lines`/`summarize` (workflow/rules/example.smk)
count and aggregate them.

Run from the repo root:

```bash
# Dry-run first
./run.py -w example -c tests/tutorial/config_tutorial.yaml -q 2 -n

# Then for real
./run.py -w example -c tests/tutorial/config_tutorial.yaml -q 2
```

Check the result:

```bash
cat tests/tutorial/output/summary.tsv
```

Expected:

```
sample	r1_lines	r2_lines
SAMPLE_A	4	4
SAMPLE_B	2	2
```

`tests/tutorial/output/` is gitignored -- delete it freely to re-run from
scratch (`rm -rf tests/tutorial/output`).
