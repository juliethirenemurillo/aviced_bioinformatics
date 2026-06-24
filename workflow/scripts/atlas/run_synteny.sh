#!/bin/bash
# run_synteny.sh -- AviCeD Paper 1 synteny step driver
#
# Loops 7 target genes x 3 species through synteny_call.py (module C).
# NLRP1 is excluded -- needs gene tree first (separate track).
#
# Outputs:
#   results/synteny_results.tsv  -- one row per (target, species), appended
#
# Run from:
#   ~/comparativegenom/aviced/blast_disambiguation/
#
# Requires: conda activate aviced; python in PATH.
# GFFs must exist in gff/ and reconciled TSV in flanker_inputs/.

set -euo pipefail

BASE="$HOME/comparativegenom/aviced/blast_disambiguation"
SCRIPTS="$BASE/scripts"
GFFDIR="$BASE/gff"
ORTHOLOGS="$BASE/flanker_inputs/reconciled_v2.tsv"
OUTDIR="$BASE/results"
mkdir -p "$OUTDIR"
OUT="$OUTDIR/synteny_results.tsv"

# Remove stale output so we get a clean header on first row
rm -f "$OUT"

# 7 targets (NLRP1 excluded -- gene tree first)
TARGETS=(BAX RBCK1 CARD8 AIM2 NAIP MEFV RIPK3)
# Note: PYRIN's GFF Name= is MEFV -- using MEFV here so the GFF lookup hits

# 3 species and their GFFs
declare -A GFF=(
    [chicken]="$GFFDIR/chicken_GRCg7b.gff"
    [duck_ZJU1.0]="$GFFDIR/duck_ZJU1.0.gff"
    [duck_T2T]="$GFFDIR/duck_T2T.gff"
)

HUMAN_GFF="$GFFDIR/human.gff"
LOG="$OUTDIR/synteny_run.log"
: > "$LOG"

echo "[$(date '+%H:%M:%S')] synteny run starting" | tee -a "$LOG"
echo "  targets: ${TARGETS[*]}" | tee -a "$LOG"
echo "  orthologs: $ORTHOLOGS" | tee -a "$LOG"
echo "  output: $OUT" | tee -a "$LOG"
echo "" | tee -a "$LOG"

for target in "${TARGETS[@]}"; do
    for species in chicken duck_ZJU1.0 duck_T2T; do
        echo "[$(date '+%H:%M:%S')] ${target} / ${species}" | tee -a "$LOG"
        python "$SCRIPTS/synteny_call.py" \
            --target "$target" \
            --human-gff "$HUMAN_GFF" \
            --bird-gff "${GFF[$species]}" \
            --orthologs "$ORTHOLOGS" \
            --species "$species" \
            --out-tsv "$OUT" \
            2>>"$LOG" || {
                echo "  ERROR: synteny_call failed for $target/$species -- check $LOG" | tee -a "$LOG"
            }
    done
    echo "" | tee -a "$LOG"
done

echo "[$(date '+%H:%M:%S')] done" | tee -a "$LOG"
echo ""
echo "=== results ==="
# print key columns only: target, species, tier, conserved, slot_status, slot_inner
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) col[$i]=i; print}
            NR>1 {printf "%-10s %-12s %-16s %-10s %-26s %s\n",
                  $col["target"],$col["species"],$col["tier"],
                  $col["conserved"],$col["slot_status"],$col["slot_inner"]}' "$OUT"