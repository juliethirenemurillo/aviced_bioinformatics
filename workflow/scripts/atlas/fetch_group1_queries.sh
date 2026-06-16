#!/usr/bin/env bash
# fetch_group1_queries.sh
# Run on the CSD3 LOGIN node (needs internet). Fetches the 8 Group-1 human
# reference proteins from UniProt with retry, builds a combined query FASTA,
# validates each header, and reports pass/fail.
#
# Fix vs last session: curl --retry 3 --retry-delay 2 + sleep 1 (UniProt 503 rate-limit).
set -uo pipefail   # deliberately NOT -e: attempt all 8, then report

BASE="${BASE:-$HOME/comparativegenom/aviced/blast_disambiguation}"
QDIR="$BASE/queries"; LOGDIR="$BASE/logs"
mkdir -p "$QDIR" "$LOGDIR"
COMBINED="$QDIR/group1_human_queries.fasta"
MANIFEST="$QDIR/group1_manifest.tsv"
LOG="$LOGDIR/fetch_group1_$(date +%Y%m%d_%H%M%S).log"

# gene <TAB> human UniProt accession  (MLKL = human query for the duck-only search)
GENES=$(cat <<'EOF'
BAX	Q07812
MAVS	Q7Z434
RBCK1	Q9BYM8
BID	P55957
BAD	Q92934
BBC3	Q96PG8
PMAIP1	Q13794
MLKL	Q8NB16
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
      # soft check: does the gene symbol show up in the header? (warn only)
      grep -qiE "(GN=${gene}\b|\b${gene}_)" <<<"$hdr" || \
        echo "   ~ note: '$gene' not literally in header (check protein name)" | tee -a "$LOG"
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
echo "----"
echo "sequences in combined file: $n (expected 8)"
echo "headers:"; grep '^>' "$COMBINED" || true
echo "combined: $COMBINED"
echo "manifest: $MANIFEST"
if [ "$n" -eq 8 ] && [ "$fail" -eq 0 ]; then
  echo "RESULT: ALL 8 OK"
else
  echo "RESULT: CHECK NEEDED — see $MANIFEST and $LOG"; exit 1
fi
