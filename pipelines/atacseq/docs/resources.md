# Getting the resources

| File | Config key | Source |
|---|---|---|
| Genome FASTA | `reference.genome` | GENCODE primary assembly, via `../rnaseq/resources/get_reference.sh` |
| Annotation GTF | `reference.gtf` | the matching GENCODE GTF, same script |
| Blacklist | `reference.blacklist` | ENCODE blacklist v2 (Amemiya et al. 2019), via `resources/get_atac_resources.sh` |
| somalier sites | `qc.somalier_sites` | somalier's hg38 sites, same script (optional) |

## Genome and annotation

atacseq uses the same genome and GTF as rnaseq, so a project with both assays shares one
download. The reasons for GENCODE's primary assembly are in
[rnaseq's resources.md](../../rnaseq/docs/resources.md). In short: no alternative
haplotypes, so no reads lost to multimapping between a chromosome and its alternates. The
GTF is used only for TSSs and peak annotation, so any release on the same assembly works.

## Blacklist and somalier sites

```bash
resources/get_atac_resources.sh --build hg38 --dest /refs/GRCh38_gencode50
resources/get_atac_resources.sh --build mm10 --dest /refs/GRCm38 --ensembl   # Ensembl names
```

The blacklist exists for hg38, hg19 and mm10 (and fly and worm). For GRCm39 there is no
official v2 list yet. Lifting mm10's over with UCSC liftOver is the usual stopgap, or run
without one (`reference.blacklist: null`), accepting a few artefact peaks at satellite
repeats. somalier's sites are published for hg38 (and GRCh37). For other genomes, skip
the identity check.

## What gets built (`-w reference`)

Under `reference.dir`: `genome.fa` (+ `.fai`), `chrom.sizes`, `mito_name.txt`,
`nuclear.bed`, `effective_genome_size.txt`, `genes.tsv`, `tss.bed`, `blacklist.bed` and
`bowtie2/`. Point several projects at one `reference.dir` to build the index once. Use a
directory separate from rnaseq's (e.g. `<ref>/sm-core-atacseq`), since the two pipelines
build different files.

## Memory and disk (human)

| Step | Memory | Time | Disk |
|---|---|---|---|
| bowtie2-build, GRCh38, 16 threads | ~16 GB | 1–2 h | 4 GB index + 3 GB genome copy |
| bowtie2, 40 M pairs, 16 threads | < 8 GB | ~1–2 h | |
| Picard MarkDuplicates | 16 GB (heap 12.8 GB) | ~30 min | |
| final BAM | | | ~1 GB per 10 M fragments |
