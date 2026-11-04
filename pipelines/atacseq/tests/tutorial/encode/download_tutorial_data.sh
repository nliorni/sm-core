#!/usr/bin/env bash
## ============================================================================
## sm-core/atacseq validation tutorial: ENCODE ATAC-seq of GM12878 and K562
## (Snyder lab, 2x100 bp, two biological replicates each, ~38 M pairs per
## library, ~19 GB), with ENCODE's own processed results as the truth set:
##
##   truth/<cell>_conservative_idr.bed.gz   reproducible (IDR) peaks, GRCh38
##   truth/<cell>_pseudoreplicated.bed.gz   ENCODE's default peak set, GRCh38
##   truth/encode_qc.tsv                     per-replicate QC from the ENCODE
##                                           pipeline: mito fraction, NRF,
##                                           PBC1/2, duplicates, TSS enrichment,
##                                           FRiP
##
##   tests/tutorial/encode/download_tutorial_data.sh [--genome G.fa.gz --gtf A.gtf.gz]
##
## --genome / --gtf: an existing GENCODE GRCh38 primary-assembly genome and
## GTF; fetched with rnaseq's get_reference.sh otherwise. Then validate with
## tests/tutorial/encode/validate.py (README.md).
## ============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE="$(cd "${HERE}/../../.." && pwd)"
GENOME=""
GTF=""
while (( $# )); do
  case "$1" in
    --genome) GENOME="$2"; shift 2 ;;
    --gtf)    GTF="$2"; shift 2 ;;
    *) echo "unknown option $1" >&2; exit 1 ;;
  esac
done
mkdir -p "${HERE}/reference" "${HERE}/fastq" "${HERE}/truth"

if [[ -z "${GENOME}" ]]; then
  "${PIPELINE}/../rnaseq/resources/get_reference.sh" --dest "${HERE}/reference/gencode"
  GENOME="${HERE}/reference/gencode/GRCh38.primary_assembly.genome.fa.gz"
  GTF="${HERE}/reference/gencode/gencode.v50.primary_assembly.annotation.gtf.gz"
fi
ln -sf "${GENOME}" "${HERE}/reference/genome.fa.gz"
ln -sf "${GTF}" "${HERE}/reference/genes.gtf.gz"
"${PIPELINE}/resources/get_atac_resources.sh" --build hg38 --dest "${HERE}/reference" > /dev/null

ENCODE="https://www.encodeproject.org"
fetch() {
  # ENCODE publishes each file's MD5; a download that does not match is deleted.
  local acc="$1" out="$2" md5
  md5="$(curl -sfL -H 'Accept: application/json' "${ENCODE}/files/${acc}/?format=json" \
         | python3 -c 'import json, sys; print(json.load(sys.stdin)["md5sum"])')"
  if [[ -s "${out}" ]] && [[ "$(md5sum < "${out}" | cut -d' ' -f1)" == "${md5}" ]]; then
    echo "  already present: ${out##*/}"; return
  fi
  echo "  ${acc} -> ${out##*/}"
  local ext="${out##*.}"; [[ "${out}" == *.fastq.gz ]] && ext=fastq.gz; [[ "${out}" == *.bed.gz ]] && ext=bed.gz
  curl -sfL --retry 3 -C - -o "${out}.part" "${ENCODE}/files/${acc}/@@download/${acc}.${ext}"
  [[ "$(md5sum < "${out}.part" | cut -d' ' -f1)" == "${md5}" ]] || { rm -f "${out}.part"; echo "${acc}: MD5 mismatch" >&2; exit 1; }
  mv "${out}.part" "${out}"
}

echo "== reads"
fetch ENCFF411MVZ "${HERE}/fastq/GM12878_1_R1.fastq.gz"
fetch ENCFF736ZWY "${HERE}/fastq/GM12878_1_R2.fastq.gz"
fetch ENCFF196VAZ "${HERE}/fastq/GM12878_2_R1.fastq.gz"
fetch ENCFF195FRB "${HERE}/fastq/GM12878_2_R2.fastq.gz"
fetch ENCFF391BFJ "${HERE}/fastq/K562_1_R1.fastq.gz"
fetch ENCFF186CQZ "${HERE}/fastq/K562_1_R2.fastq.gz"
fetch ENCFF440UAD "${HERE}/fastq/K562_2_R1.fastq.gz"
fetch ENCFF350ZZR "${HERE}/fastq/K562_2_R2.fastq.gz"

echo "== ENCODE peaks (GRCh38)"
fetch ENCFF346CZA "${HERE}/truth/GM12878_conservative_idr.bed.gz"
fetch ENCFF470YYO "${HERE}/truth/GM12878_pseudoreplicated.bed.gz"
fetch ENCFF695IGF "${HERE}/truth/K562_conservative_idr.bed.gz"
fetch ENCFF558BLC "${HERE}/truth/K562_pseudoreplicated.bed.gz"

echo "== ENCODE per-replicate QC"
python3 - "${HERE}/truth/encode_qc.tsv" <<'PY'
import json, sys, urllib.request
out = sys.argv[1]
experiments = {"GM12878": "ENCSR095QNB", "K562": "ENCSR483RKN"}
rows = []
for cell, acc in experiments.items():
    req = urllib.request.Request(f"https://www.encodeproject.org/experiments/{acc}/?frame=embedded",
                                 headers={"Accept": "application/json"})
    exp = json.load(urllib.request.urlopen(req, timeout=120))
    per_rep = {}
    for f in exp["files"]:
        reps = f.get("biological_replicates") or []
        if f.get("assembly") != "GRCh38" or len(reps) != 1 or f.get("status") != "released":
            continue
        rep = per_rep.setdefault(reps[0], {})
        for q in f.get("quality_metrics", []):
            kind = q["@type"][0]
            if kind == "AtacLibraryQualityMetric":
                rep.update(NRF=q["NRF"], PBC1=q["PBC1"], PBC2=q["PBC2"],
                           duplicate_fraction=q["pct_duplicate_reads"] / 100)
            elif kind == "AtacAlignmentQualityMetric" and q.get("processing_stage") == "unfiltered":
                rep.update(aligned_rate=q["pct_mapped_reads"] / 100, mito_fraction=q["frac_mito_reads"],
                           nucleosome_free_fraction=q["frac_reads_in_nfr"])
            elif kind == "AtacAlignmentEnrichmentQualityMetric":
                rep.update(tss_enrichment=q["tss_enrichment"])
            elif kind == "AtacPeakEnrichmentQualityMetric" and f.get("output_type") == "pseudoreplicated peaks":
                rep.update(FRiP=q["frip"])
    for r, m in sorted(per_rep.items()):
        rows.append({"sample": f"{cell}_{r}", **m})
cols = ["sample", "aligned_rate", "mito_fraction", "duplicate_fraction", "NRF", "PBC1", "PBC2",
        "FRiP", "tss_enrichment", "nucleosome_free_fraction"]
with open(out, "w") as h:
    h.write("\t".join(cols) + "\n")
    for row in rows:
        h.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")
print(open(out).read())
PY
echo "Done. Run: ./run.py -w all -c tests/tutorial/encode/config_tutorial.yaml -q 32"
