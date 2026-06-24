#!/usr/bin/env bash
# fetch_group2_queries.sh   — run on the CSD3 LOGIN node (needs internet)
# Fetches the 16 Group-2 human reference proteins from UniProt with retry,
# builds a combined query FASTA, validates headers, reports pass/fail.
set -uo pipefail

BASE="${BASE:-$HOME/comparativegenom/aviced/blast_disambiguation}"
QDIR="$BASE/queries"; LOGDIR="$BASE/logs"; mkdir -p "$QDIR" "$LOGDIR"
COMBINED="$QDIR/group2_human_queries.fasta"
MANIFEST="$QDIR/group2_manifest.tsv"
LOG="$LOGDIR/fetch_group2_$(date +%Y%m%d_%H%M%S).log"

GENES=$(cat <<'EOF'
AIM2	O14862
ASC	Q9ULZ3
BIM	O43521
CARD8	Q9Y2G2
DDX58	O95786
GSDMA	Q96QA5
GSDMB	Q8TAX9
GSDMC	Q9BYG8
GSDMD	P57764
NAIP	Q13075
NLRP1	Q9C000
NLRP3	Q96P20
PYRIN	O15553
RIPK3	Q9Y572
RNF31	Q96EP0
ZBP1	Q9H171
EOF
)

: > "$COMBINED"
printf "gene\taccession\tfetched\theader\n" > "$MANIFEST"
fail=0

while IFS=$'\t' read -r gene acc; do
  [ -z "${gene:-}" ] && continue
  out="$QDIR/${gene}_${acc}.fasta"
  echo "[$(date +%T)] fetch $gene ($acc)" | tee -a "$LOG"
  if curl -fsSL --retry 3 --retry-delay 2 \
        "https://rest.uniprot.org/uniprotkb/${acc}.fasta" -o "$out" 2>>"$LOG"; then
    hdr=$(head -1 "$out")
    if [[ "$hdr" == ">"* ]] && grep -q "$acc" <<<"$hdr"; then
      cat "$out" >> "$COMBINED"
      printf "%s\t%s\tYES\t%s\n" "$gene" "$acc" "$hdr" >> "$MANIFEST"
      echo "   ok: $hdr" | tee -a "$LOG"
    else
      printf "%s\t%s\tBAD_CONTENT\t%s\n" "$gene" "$acc" "${hdr:-<empty>}" >> "$MANIFEST"
      echo "   !! unexpected content for $gene ($acc)" | tee -a "$LOG"; fail=1
    fi
  else
    printf "%s\t%s\tFAILED\t-\n" "$gene" "$acc" >> "$MANIFEST"
    echo "   !! fetch FAILED for $gene ($acc)" | tee -a "$LOG"; fail=1
  fi
  sleep 1
done <<< "$GENES"

n=$(grep -c '^>' "$COMBINED" || true)
echo "----"; echo "sequences in combined file: $n (expected 16)"
echo "headers:"; grep '^>' "$COMBINED" || true
echo "combined: $COMBINED"; echo "manifest: $MANIFEST"
if [ "$n" -eq 16 ] && [ "$fail" -eq 0 ]; then echo "RESULT: ALL 16 OK"; else echo "RESULT: CHECK NEEDED — see $MANIFEST / $LOG"; exit 1; fi
