#!/usr/bin/env bash
# ncbi_orthologs_group1.sh
# Run on the CSD3 LOGIN node (needs internet). Independent, synteny-aware read
# of the same 8 genes via NCBI Orthologs: does NCBI also see no chicken/duck
# ortholog? Writes a SEPARATE results file — never the curated table.
#
# Method: E-utilities to resolve human GeneIDs, then the NCBI Orthologs pairs
# file (gene_orthologs.gz) to read chicken/duck ortholog calls. No 'datasets'
# CLI needed (it isn't installed on CSD3).
set -uo pipefail

BASE="${BASE:-$HOME/comparativegenom/aviced/blast_disambiguation}"
DATA="$BASE/ncbi"; RES="$BASE/results"; LOGDIR="$BASE/logs"
mkdir -p "$DATA" "$RES" "$LOGDIR"
LOG="$LOGDIR/ncbi_orthologs_$(date +%Y%m%d_%H%M%S).log"
OUT="$RES/ncbi_orthologs_group1.tsv"

CHICKEN=9031; DUCK1=8839; DUCK2=8840           # G. gallus ; A. platyrhynchos (+ subspecies)
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
KEY=""; [ -n "${NCBI_API_KEY:-}" ] && KEY="&api_key=${NCBI_API_KEY}"   # optional, raises rate limit
GENES="BAX MAVS RBCK1 BID BAD BBC3 PMAIP1 MLKL"

# 1) one-time download of the NCBI Orthologs pairs file
ORTHO="$DATA/gene_orthologs.gz"
if [ ! -s "$ORTHO" ]; then
  echo "downloading gene_orthologs.gz (one-time, ~hundreds of MB)…" | tee -a "$LOG"
  curl -fSL --retry 3 --retry-delay 5 \
    "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz" -o "$ORTHO" 2>>"$LOG" \
    || { echo "download failed — see $LOG"; exit 1; }
fi

esearch_geneid () {  # $1=symbol -> human GeneID (or empty)
  curl -fsSL --retry 3 --retry-delay 2 \
    "${EUTILS}/esearch.fcgi?db=gene&term=${1}%5BPreferred+Symbol%5D+AND+9606%5Btaxid%5D&retmode=json${KEY}" \
    | grep -o '"idlist":\[[^]]*\]' | grep -oE '[0-9]+' | head -1
}
symbol_for_geneid () {  # $1=GeneID -> symbol (cosmetic; GeneID is the authoritative bit)
  [ -z "${1:-}" ] && { echo ""; return; }
  curl -fsSL --retry 3 --retry-delay 2 \
    "${EUTILS}/esummary.fcgi?db=gene&id=${1}&retmode=json${KEY}" \
    | tr ',' '\n' | grep -m1 '"name"' | sed 's/.*"name":"//; s/".*//'
}

printf "gene\thuman_geneid\tchicken_ortholog\tduck_ortholog\tduck_taxid\tsource\n" > "$OUT"

for g in $GENES; do
  echo "[$(date +%T)] $g" | tee -a "$LOG"
  hid=$(esearch_geneid "$g"); sleep 0.5
  if [ -z "${hid:-}" ]; then
    printf "%s\t-\t-\t-\t-\tHUMAN_GENEID_LOOKUP_FAILED\n" "$g" >> "$OUT"; continue
  fi
  # gene_orthologs cols: tax_id GeneID relationship Other_tax_id Other_GeneID
  # human may sit in (col1,col2) or (col4,col5) — check both orientations
  partners=$(zcat "$ORTHO" | awk -v G="$hid" \
      '($1==9606 && $2==G){print $4"\t"$5} ($4==9606 && $5==G){print $1"\t"$2}')
  ck=$(awk -v C="$CHICKEN" '$1==C{print $2}' <<<"$partners" | head -1)
  dline=$(awk -v D1="$DUCK1" -v D2="$DUCK2" '$1==D1||$1==D2{print $1"\t"$2}' <<<"$partners" | head -1)
  dtax=$(cut -f1 <<<"$dline"); dk=$(cut -f2 <<<"$dline")
  cks=$(symbol_for_geneid "$ck"); sleep 0.4
  dks=$(symbol_for_geneid "$dk"); sleep 0.4
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$g" "$hid" \
    "${ck:-none}${cks:+ ($cks)}" "${dk:-none}${dks:+ ($dks)}" "${dtax:-none}" \
    "NCBI Orthologs (gene_orthologs.gz)" >> "$OUT"
  echo "   human=$hid chicken=${ck:-none} duck=${dk:-none}" | tee -a "$LOG"
done

echo "----"; column -t -s$'\t' "$OUT" 2>/dev/null || cat "$OUT"
echo "written: $OUT  (separate file — curated table untouched)"
