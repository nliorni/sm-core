#!/usr/bin/env bash
## ============================================================================
## Download the ATAC-specific reference files for sm-core/atacseq and print
## the config lines for them. The genome and GTF are the same as rnaseq's:
## fetch them with pipelines/rnaseq/resources/get_reference.sh.
##
##   get_atac_resources.sh --build hg38 --dest /refs/GRCh38_gencode50
##   get_atac_resources.sh --build mm10 --dest /refs/GRCm38 --ensembl
##
## Options
##   --build    hg38 (default) | hg19 | mm10
##   --dest     directory to download into (created; required)
##   --ensembl  rename sequences to Ensembl style (chr1 -> 1, chrM -> MT), for
##              an Ensembl genome; GENCODE genomes use the UCSC names as is
##
## What it fetches:
##   - the ENCODE blacklist v2 (Amemiya, Kundaje & Boyle 2019), regions with
##     anomalous signal in every functional-genomics experiment;
##   - hg38 only: somalier's ~17,000 common SNP sites, for the optional
##     sample identity check (qc.somalier_sites).
##
## Needs curl and gzip. Works behind an HTTP proxy (https_proxy) as is.
## ============================================================================
set -euo pipefail

BUILD=hg38
DEST=""
ENSEMBL=0
die() { echo "ERROR: $*" >&2; exit 1; }

while (( $# )); do
  case "$1" in
    --build)   BUILD="$2"; shift 2 ;;
    --dest)    DEST="$2"; shift 2 ;;
    --ensembl) ENSEMBL=1; shift ;;
    -h|--help) sed -n '3,/^## ====/{/^## ====/!p}' "$0" | sed 's/^## \{0,1\}//'; exit 0 ;;
    *) die "unknown option $1 (see --help)" ;;
  esac
done
[[ -n "${DEST}" ]] || die "--dest is required (see --help)"
case "${BUILD}" in hg38|hg19|mm10) ;; *) die "--build must be hg38, hg19 or mm10" ;; esac

BLACKLIST_URL="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/${BUILD}-blacklist.v2.bed.gz"
SITES_URL="https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz"

mkdir -p "${DEST}"
DEST="$(cd "${DEST}" && pwd)"
cd "${DEST}"

fetch() {
  local url="$1" out="$2"
  if [[ -s "${out}" ]] && gzip -t "${out}" 2>/dev/null; then
    echo "  already present: ${out}"
    return
  fi
  echo "  fetching ${out}"
  curl -sfL --retry 3 --connect-timeout 30 -o "${out}.part" "${url}" || { rm -f "${out}.part"; die "download failed: ${url}"; }
  gzip -t "${out}.part" || { rm -f "${out}.part"; die "${out}: not a valid gzip file"; }
  mv "${out}.part" "${out}"
}

blacklist="${BUILD}-blacklist.v2.bed.gz"
fetch "${BLACKLIST_URL}" "${blacklist}"
if (( ENSEMBL )); then
  out="${BUILD}-blacklist.v2.ensembl.bed.gz"
  gzip -cd "${blacklist}" | sed -e 's/^chrM\t/MT\t/' -e 's/^chr//' | gzip > "${out}"
  blacklist="${out}"
fi
echo "  blacklist: $(gzip -cd "${blacklist}" | wc -l) regions"

sites=""
if [[ "${BUILD}" == hg38 ]]; then
  sites="sites.hg38.vcf.gz"
  fetch "${SITES_URL}" "${sites}"
  if (( ENSEMBL )); then
    out="sites.hg38.ensembl.vcf.gz"
    gzip -cd "${sites}" | sed -e 's/^chr//' -e 's/^##contig=<ID=chr/##contig=<ID=/' | gzip > "${out}"
    sites="${out}"
  fi
fi

cat <<EOF

Done. Put these in the project config:

reference:
  blacklist: "${DEST}/${blacklist}"
EOF
if [[ -n "${sites}" ]]; then
  cat <<EOF
qc:
  somalier_sites: "${DEST}/${sites}"
EOF
fi
