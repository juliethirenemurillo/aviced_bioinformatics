#!/bin/bash
# run_synteny.sh -- AviCeD Paper 1 synteny step driver (v2: with protein map)
#
# 1. Builds protein_accession -> GFF gene Name maps for each bird GFF
# 2. Loops 7 target genes x 3 species through synteny_call.py
#
# NLRP1 excluded -- needs gene tree first.
#
# Run from: ~/comparativegenom/aviced/blast_disambiguation/
# Requires: conda activate aviced; python in PATH.

set -euo pipefail

BASE="$HOME/comparativegenom/aviced/blast_disambiguation"
SCRIPTS="$BASE/scripts"
GFFDIR="$BASE/gff"
ORTHOLOGS="$BASE/flanker_inputs/reconciled_v2.tsv"
MAPDIR="$BASE/flanker_inputs/protein_maps"
OUTDIR="$BASE/results"
mkdir -p "$OUTDIR" "$MAPDIR"
OUT="$OUTDIR/synteny_results_v2.tsv"
rm -f "$OUT"

TARGETS=(BAX RBCK1 CARD8 AIM2 NAIP MEFV RIPK3)

declare -A GFF=(
    [chicken]="$GFFDIR/chicken_GRCg7b.gff"
    [duck_ZJU1.0]="$GFFDIR/duck_ZJU1.0.gff"
    [duck_T2T]="$GFFDIR/duck_T2T.gff"
)

declare -A PROT_MAP=(
    [chicken]="$MAPDIR/chicken_protein_map.tsv"
    [duck_ZJU1.0]="$MAPDIR/duck_ZJU_protein_map.tsv"
    [duck_T2T]="$MAPDIR/duck_T2T_protein_map.tsv"
)

HUMAN_GFF="$GFFDIR/human.gff"
LOG="$OUTDIR/synteny_run_v2.log"
: > "$LOG"

echo "[$(date '+%H:%M:%S')] === STEP 1: building protein maps ===" | tee -a "$LOG"
for species in chicken duck_ZJU1.0 duck_T2T; do
    map="${PROT_MAP[$species]}"
    if [[ -f "$map" ]]; then
        echo "  $species: map exists ($(wc -l < "$map") rows) -- skipping" | tee -a "$LOG"
    else
        echo "  $species: building protein map from ${GFF[$species]}..." | tee -a "$LOG"
        python "$SCRIPTS/build_protein_gene_map.py" "${GFF[$species]}" > "$map" 2>>"$LOG"
        echo "  $species: $(wc -l < "$map") rows" | tee -a "$LOG"
    fi
done

echo "" | tee -a "$LOG"
echo "[$(date '+%H:%M:%S')] === STEP 2: synteny calls ===" | tee -a "$LOG"
echo "  targets: ${TARGETS[*]}" | tee -a "$LOG"

for target in "${TARGETS[@]}"; do
    for species in chicken duck_ZJU1.0 duck_T2T; do
        echo "[$(date '+%H:%M:%S')] ${target} / ${species}" | tee -a "$LOG"
        python "$SCRIPTS/synteny_call.py" \
            --target "$target" \
            --human-gff "$HUMAN_GFF" \
            --bird-gff "${GFF[$species]}" \
            --orthologs "$ORTHOLOGS" \
            --species "$species" \
            --protein-map "${PROT_MAP[$species]}" \
            --out-tsv "$OUT" \
            2>>"$LOG" || {
                echo "  ERROR: $target/$species -- see $LOG" | tee -a "$LOG"
            }
    done
    echo "" | tee -a "$LOG"
done

echo "[$(date '+%H:%M:%S')] done" | tee -a "$LOG"
echo ""
echo "=== results ==="
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) col[$i]=i; print}
            NR>1 {printf "%-10s %-12s %-16s %-10s %-26s %-40s %s\n",
                  $col["target"],$col["species"],$col["tier"],
                  $col["conserved"],$col["slot_status"],$col["slot_inner"],$col["reason"]}' "$OUT"
