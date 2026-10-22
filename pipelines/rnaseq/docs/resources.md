# Getting the resources

Everything sm-core/rnaseq reads but does not build itself: a genome, its gene annotation,
and, optionally, spike-in sequences. That's all. Every index (STAR, Salmon, the
transcriptome, the gene tables) is **built by the pipeline** from these two or three files
(`-w reference`). So nothing else needs downloading, and no count can silently come from a
different annotation than the one the matrices are labelled with.

All URLs below were checked on 2026-10-22. Sizes are download sizes.

## At a glance

| Resource | Config key | Size (human) | Needed by |
|---|---|---|---|
| Genome FASTA, primary assembly | `reference.genome` | 0.85 GB gz (3.1 GB) | everything |
| Gene annotation GTF, same source and release | `reference.gtf` | 0.12 GB gz (4.4 GB) | everything |
| Spike-in FASTA (e.g. ERCC92) | `reference.spikein_fasta` | 30 KB | optional |
| *Built by the pipeline:* STAR + Salmon indices, transcriptome, gene tables | `reference.dir` | 56 GB | built once |

Plan for about **60 GB per genome + annotation**: 1 GB of downloads and 56 GB of indices,
which are built once and shared between projects.

**The one rule that matters: genome and GTF from the same source and release.** The GTF's
sequence names must be the genome's (`chr1` in GENCODE and UCSC, `1` in Ensembl), and its
transcripts must be modelled on that same assembly. Mixing sources is the most common way
to get a quietly wrong count matrix. The pipeline catches the loud half of the problem. GTF
records on sequences that are not in the genome are dropped and counted in
`logs/reference/prepare_annotation.log`, and a GTF whose names match nothing stops the
run. But it cannot detect an annotation from a different *patch* of the same assembly.
Download both from the same release directory, which is what the script below does.

## 1. The quick way: `get_reference.sh`

```bash
cd pipelines/rnaseq
resources/get_reference.sh --dest /path/to/references/GRCh38_gencode50              # human
resources/get_reference.sh --species mouse --dest /path/to/references/GRCm39_gencodeM39
resources/get_reference.sh --source ensembl --species danio_rerio --dest /path/to/references/GRCz11_ens116
```

For each one it:

1. picks the right two files: the primary assembly, and the comprehensive annotation from
   the same release (sections 2–4 explain the choice);
2. downloads them, resuming if interrupted, and verifies each against the publisher's own
   checksums (GENCODE `MD5SUMS`, Ensembl `CHECKSUMS`). A file that fails is deleted, not
   kept;
3. checks that every sequence the GTF names exists in the genome;
4. writes a `MANIFEST.tsv` (file, URL, MD5, date) for provenance;
5. prints the `reference:` block to paste into the project config.

```
Done. Put this in the project config (SECTION 2):

reference:
  genome: "/path/to/references/GRCh38_gencode50/GRCh38.primary_assembly.genome.fa.gz"
  gtf: "/path/to/references/GRCh38_gencode50/gencode.v50.primary_assembly.annotation.gtf.gz"
  # Build the indices once, here, and share them between projects:
  dir: "/path/to/references/GRCh38_gencode50/sm-core-rnaseq"
```

Add `--spikein ercc` for the ERCC92 spike-ins (section 5), `--release N` for another
release, and `--check` to resolve the URLs without downloading. The script needs only
`curl`, `gzip`, `md5sum`, `sum` and `unzip`. `--help` lists every option.

Then [build the indices](#6-building-the-indices) once.

## 2. Human: GRCh38, GENCODE (recommended)

GENCODE is Ensembl's annotation with UCSC-style sequence names. It's what GTEx, ENCODE,
TCGA and nf-core use. Release **50** (= Ensembl 116) at the time of writing; check
<https://www.gencodegenes.org/human/> for newer ones.

```bash
REL=50
BASE=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${REL}
curl -O ${BASE}/GRCh38.primary_assembly.genome.fa.gz                     # 0.85 GB
curl -O ${BASE}/gencode.v${REL}.primary_assembly.annotation.gtf.gz        # 0.12 GB
curl -O ${BASE}/MD5SUMS && md5sum -c MD5SUMS --ignore-missing
```

Use `ftp://`, not `https://`: EBI serves the same files about ten times faster over plain
FTP.

**Which files, and why.**

- `GRCh38.primary_assembly.genome.fa.gz`: chromosomes plus unlocalised and unplaced
  scaffolds, with **no alternate haplotypes**. Do not use `GRCh38.p14.genome.fa.gz` (the
  "ALL" file). Its ~260 alt contigs duplicate large stretches of chr6 (HLA), chr17, chr19
  and others, so reads there become multimappers: STAR drops them from its gene counts,
  and Salmon's decoys absorb them.
- `gencode.vNN.primary_assembly.annotation.gtf.gz`: the comprehensive annotation on the
  same sequences. The plain `gencode.vNN.annotation.gtf.gz` ("CHR") covers the reference
  chromosomes only. It also works, since all its records exist in the primary assembly,
  but it misses the few genes on scaffolds.
- **Comprehensive, not `basic`.** Salmon needs every transcript a read may have come from
  in order to apportion reads correctly. A missing isoform does not remove its reads; it
  pushes them onto the isoforms that are present.

## 3. Mouse: GRCm39, GENCODE

Release **M39** at the time of writing (<https://www.gencodegenes.org/mouse/>). The same
choices apply as for human.

```bash
REL=M39
BASE=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${REL}
curl -O ${BASE}/GRCm39.primary_assembly.genome.fa.gz                     # 0.77 GB
curl -O ${BASE}/gencode.v${REL}.primary_assembly.annotation.gtf.gz        # 0.09 GB
curl -O ${BASE}/MD5SUMS && md5sum -c MD5SUMS --ignore-missing
```

## 4. Any other species: Ensembl

Ensembl publishes a genome and a matching GTF for every species it annotates: the
unmasked `primary_assembly` for assemblies that have alt haplotypes, `toplevel` for those
that don't. Release **116** at the time of writing. Ensembl Genomes
(<https://ensemblgenomes.org>) covers plants, fungi, protists and invertebrates the same
way.

```bash
REL=116; SP=danio_rerio; Sp=Danio_rerio; ASM=GRCz11      # adjust all four
curl -O ftp://ftp.ensembl.org/pub/release-${REL}/fasta/${SP}/dna/${Sp}.${ASM}.dna.primary_assembly.fa.gz
curl -O ftp://ftp.ensembl.org/pub/release-${REL}/gtf/${SP}/${Sp}.${ASM}.${REL}.gtf.gz
```

Take the plain `dna` file, never `dna_sm` (soft-masked) or `dna_rm` (repeat-masked). And
take the `<release>.gtf.gz` GTF, not the `.chr.` or `.abinitio.` variants. The script
makes both choices for you (`--source ensembl --species <name>`).

Ensembl names chromosomes `1`, `2`, … `MT`, and its biotype attribute is `gene_biotype`
instead of GENCODE's `gene_type`. The pipeline reads either, so the matrices carry biotypes
and the count-composition QC works the same way. The yeast tutorial runs on an Ensembl GTF.

## 5. Spike-ins: ERCC

The ERCC ExFold / RNA Spike-In Mixes (Thermo Fisher 4456739 / 4456740) are 92 synthetic
transcripts at known concentrations. Mix 1 and Mix 2 differ by known ratios (4, 1, 2/3,
1/2, in four subgroups), which makes them a built-in truth set for fold changes.

```bash
resources/get_reference.sh --spikein ercc --dest ...      # or by hand:
curl -O https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip && unzip ERCC92.zip ERCC92.fa
curl -O https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt   # concentrations
```

```yaml
reference:
  spikein_fasta: "/path/to/references/ERCC92.fa"
```

Only the FASTA is needed. The pipeline appends it to the genome and writes one
gene/transcript/exon per sequence itself (`gene_type "spikein"`), so Thermo Fisher's
`ERCC92.gtf` is not used. The spike-ins then appear in both branches' count matrices like
any other gene, and their share of each library shows in the count-composition QC. The
same works for any other spike-in set (SIRVs, sequins): one sequence per transcript.

A different spike-in set is a different reference: give it its own `reference.dir`.

## 6. Building the indices

```bash
./run.py -w reference -c config_my_project.yaml -q 16
```

This builds under `reference.dir` (default `<output_dir>/reference/`). Every other target
builds whatever it is missing too, so this step is optional, but running it once on its
own is the clean way to share a reference (see below).

| File | Built by | Used by | Human GRCh38 + GENCODE 50 |
|---|---|---|---|
| `genome.fa` (+ `.fai`) | `prepare_genome` | everything below | 3.1 GB |
| `genes.gtf` | `prepare_annotation` (filtered, + spike-ins) | STAR index, gffread | 4.4 GB |
| `tx2gene.tsv` | `prepare_annotation` | both count matrices | 40 MB |
| `transcripts.fa` | gffread | Salmon index | 1.4 GB |
| `genes.collapsed.gtf` | `collapse_annotation` | RNA-SeQC (star) | 0.1 GB |
| `star/` | STAR `genomeGenerate` | star branch | 29 GB on disk; 38 GB RAM; 1 h 54 min |
| `salmon/` | `salmon index`, decoy-aware | salmon branch | 18 GB on disk; 20 GB RAM; 1 h 11 min |

Figures measured on this pipeline's human tutorial: GENCODE 50 + ERCC92, 16 threads each,
on a 2011-era AMD Opteron 6272 node. Expect roughly half the time on a current CPU. Each
run's own numbers are in `benchmarks/reference/`. The reference directory takes
**56 GB** in all.

### Sharing one build between projects

Point `reference.dir` at a shared directory in every project config that uses the same
genome + GTF (+ spike-ins), as the script's suggested `dir:` does. The first project builds
it, and the rest find it built. Two caveats:

- Never let two projects *build* it at the same time. Run `-w reference` once on its own
  first.
- A different GTF, release or spike-in set is a different reference. Give it its own
  directory.

### Prebuilt indices

`reference.star_index` / `reference.salmon_index` use an index built elsewhere instead.
It must come from exactly the genome + GTF (+ spike-ins) configured, because the count
matrices are annotated from that GTF. A STAR index built with another annotation counts
other genes, and a Salmon index with other transcript IDs cannot be mapped to genes. The
STAR and Salmon versions must also match the pinned ones (`workflow/envs/`). Both refuse
or misread indices from other major versions. When in doubt, let the pipeline build them.

## 7. Software

There is nothing to download by hand. Every tool comes from bioconda/conda-forge at an
exact pinned version (`workflow/envs/*.yaml`) and is built by Snakemake on first use. See
[installation.md](installation.md).

## Memory and disk at a glance (human)

| Step | RAM | Disk |
|---|---|---|
| STAR index build | 38 GB | 29 GB |
| Salmon decoy-aware index build | 20 GB | 18 GB |
| STAR alignment, per sample | 30 GB (the index, loaded) | ~1.5 GB of BAM per 10M pairs |
| Salmon quantification, per sample | 27 GB (the decoy-aware index holds the genome) | ~30 MB |
| RNA-SeQC, per sample | 0.4 GB | ~5 MB |

Peaks measured on the human tutorial (12M pairs per sample; `max_rss` in `benchmarks/`).
The defaults request 48 GB for the STAR index, 32 GB for Salmon and 40 GB for STAR
alignment, which leaves headroom. The `*_mem_mb` keys in the config set what each rule
requests. On a cluster, a job killed for memory is retried with twice, then three times,
as much (`-rt`, `scaled()`).

## Final checklist

Before the first real run:

- [ ] Genome and GTF come from **one** release directory (`MANIFEST.tsv` if you used the
      script).
- [ ] The genome is the **primary assembly** (no `_alt` contigs:
      `zcat genome.fa.gz | grep -c '_alt'` prints 0).
- [ ] The GTF is the **comprehensive** annotation, not `basic`.
- [ ] Every GTF sequence name is in the genome. The script checks this; by hand:
      `comm -13 <(zcat genome.fa.gz | grep '^>' | sed 's/^>//;s/ .*//' | sort -u) <(zcat genes.gtf.gz | grep -v '^#' | cut -f1 | sort -u)`
      prints nothing.
- [ ] `reference.dir` points at a shared location if other projects use the same genome.
- [ ] A dry-run is clean: `./run.py -w reference -c config_my_project.yaml -q 16 -n`.
- [ ] After the build, `logs/reference/prepare_annotation.log` reports few or no dropped
      records. Many dropped records mean the two files do not match.
