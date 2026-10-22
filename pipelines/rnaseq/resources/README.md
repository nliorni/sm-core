# resources

| Path | What it is |
|---|---|
| `get_reference.sh` | Downloads a matched genome + annotation (GENCODE human/mouse, or any Ensembl species), optionally the ERCC92 spike-ins, verifies them against the publisher's checksums, and prints the `reference:` block for the config. `--help` for options. |

Genomes, annotations and the indices built from them are **not** kept here. They are far
too large and are shared between projects. [docs/resources.md](../docs/resources.md) says
where they come from and why those files.
