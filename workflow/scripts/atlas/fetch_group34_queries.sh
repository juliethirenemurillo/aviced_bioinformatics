#!/usr/bin/env bash
# fetch_group34_queries.sh   — run on the CSD3 LOGIN node (needs internet)
set -uo pipefail
BASE="${BASE:-$HOME/comparativegenom/aviced/blast_disambiguation}"
QDIR="$BASE/queries"; LOGDIR="$BASE/logs"; mkdir -p "$QDIR" "$LOGDIR"
COMBINED="$QDIR/group34_human_queries.fasta"; MANIFEST="$QDIR/group34_manifest.tsv"
LOG="$LOGDIR/fetch_group34_$(date +%Y%m%d_%H%M%S).log"
GENES=$(cat <<'EOF'
CASP1	P29466
CASP4	P49662
CASP5	P51878
TRIF	Q8IUC6
NLRC4	Q9NPP4
TNFRSF10A	O00220
TNFRSF10B	O14763
EOF
)
: > "$COMBINED"; printf "gene\taccession\tfetched\theader\n" > "$MANIFEST"; fail=0
while IFS=$'\t' read -r gene acc; do
  [ -z "${gene:-}" ] && continue
  out="$QDIR/${gene}_${acc}.fasta"
  echo "[$(date +%T)] fetch $gene ($acc)" | tee -a "$LOG"
  if curl -fsSL --retry 3 --retry-delay 2 "https://rest.uniprot.org/uniprotkb/${acc}.fasta" -o "$out" 2>>"$LOG"; then
    hdr=$(head -1 "$out")
    if [[ "$hdr" == ">"* ]] && grep -q "$acc" <<<"$hdr"; then
      cat "$out" >> "$COMBINED"; printf "%s\t%s\tYES\t%s\n" "$gene" "$acc" "$hdr" >> "$MANIFEST"; echo "   ok: $hdr"|tee -a "$LOG"
    else printf "%s\t%s\tBAD_CONTENT\t%s\n" "$gene" "$acc" "${hdr:-<empty>}" >> "$MANIFEST"; echo "   !! bad content"|tee -a "$LOG"; fail=1; fi
  else printf "%s\t%s\tFAILED\t-\n" "$gene" "$acc" >> "$MANIFEST"; echo "   !! FAILED"|tee -a "$LOG"; fail=1; fi
  sleep 1
done <<< "$GENES"
n=$(grep -c '^>' "$COMBINED" || true)
echo "----"; echo "sequences: $n (expected 7)"; grep '^>' "$COMBINED" || true
if [ "$n" -eq 7 ] && [ "$fail" -eq 0 ]; then echo "RESULT: ALL 7 OK"; else echo "RESULT: CHECK NEEDED — $MANIFEST"; exit 1; fi
