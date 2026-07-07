#!/bin/bash
# build_flanker_inputs.sh  --  AviCeD synteny step (v3)
#
# Three-tier flanker ortholog identification, matching the same logic used
# for the 31-gene BLAST disambiguation:
#
#   Tier 1: Compara (Ensembl tree-based orthology)
#   Tier 2: blastp reciprocal best, with top-N candidate search
#   Tier 3: tblastn against bird GENOME (catches unannotated genes)
#
# v3 design changes vs v2:
#   * blast.tsv captures the TOP FORWARD HIT regardless of reciprocity.
#     Reverse-hit symbol recorded separately. This lets downstream code
#     accept paralogs as positional flankers (with transparent labelling).
#   * tblastn added for cases where blastp finds nothing AND Compara is
#     also empty (i.e. tier 3 fallback).
#   * Richer blast.tsv schema (see below).
#
# Inputs:
#   FLANKERS_TSV  -- TSV with cols: human_gene, human_ensembl, human_protein_fasta
#
# Outputs:
#   compara.tsv  -- same schema as v2 (human_gene, species, bird_gene,
#                   ortholog_type, identity)
#   blast.tsv    -- NEW SCHEMA:
#                   human_gene, species,
#                   forward_top_hit, forward_pident, forward_evalue, forward_qcov,
#                   reverse_symbol, reciprocal_match,
#                   tblastn_seqid, tblastn_pident, tblastn_evalue, tblastn_qcov
#
# Tunables (env vars):
#   N_FORWARD_HITS  -- how many top forward hits to search for a reciprocal best
#                     (default 10). The TOP hit is always recorded as
#                     forward_top_hit; reciprocal_match=YES is set if any of
#                     the top-N reciprocates.

set -euo pipefail

FLANKERS_TSV="${1:-flankers.tsv}"
OUTDIR="${2:-flanker_inputs}"
BLASTDB="${BLASTDB:-$HOME/comparativegenom/aviced/blast_disambiguation/blastdb}"
N_FORWARD_HITS="${N_FORWARD_HITS:-10}"

mkdir -p "$OUTDIR"
COMPARA_OUT="$OUTDIR/compara.tsv"
BLAST_OUT="$OUTDIR/blast_v3.tsv"
LOG="$OUTDIR/build_v3.log"

# Species mappings
declare -A COMPARA_SPECIES=(
    [chicken]="gallus_gallus"
    [duck_ZJU1.0]="anas_platyrhynchos_platyrhynchos"
)

declare -A BIRD_PROT_DB=(
    [chicken]="chicken_prot"
    [duck_ZJU1.0]="duck_zju_prot"
    [duck_T2T]="duck_T2T_protein"
)

declare -A BIRD_GENOME_DB=(
    [chicken]="chicken_gen"
    [duck_ZJU1.0]="duck_zju_gen"
    [duck_T2T]="duck_t2t_gen"
)

HUMAN_DB="human_protein"

# Headers
echo -e "human_gene\tspecies\tbird_gene\tortholog_type\tidentity" > "$COMPARA_OUT"
echo -e "human_gene\tspecies\tforward_top_hit\tforward_pident\tforward_evalue\tforward_qcov\treverse_symbol\treciprocal_match\ttblastn_seqid\ttblastn_pident\ttblastn_evalue\ttblastn_qcov" > "$BLAST_OUT"
: > "$LOG"

echo "[$(date)] starting flanker input build v3 (3-tier: compara -> blastp -> tblastn)" | tee -a "$LOG"

# ---------- Compara fetch ----------
fetch_compara () {
    local human_gene="$1" ensembl_id="$2" species_label="$3" ensembl_species="$4"
    local url="https://rest.ensembl.org/homology/id/homo_sapiens/${ensembl_id}?target_species=${ensembl_species};type=orthologues;format=condensed"
    local resp
    if ! resp=$(curl -s -H 'Accept: application/json' "$url" 2>>"$LOG"); then
        echo "[compara] curl failed: $human_gene $species_label" >> "$LOG"
        return
    fi
    echo "$resp" | jq -r --arg sp "$ensembl_species" --arg g "$human_gene" --arg s "$species_label" '
        .data[0].homologies[]?
        | select(.species == $sp)
        | [$g, $s, (.id // ""), (.type // ""), ((.target.perc_id // 0)|tostring)]
        | @tsv
    ' >> "$COMPARA_OUT" 2>>"$LOG" || true
}

# ---------- BLAST helpers ----------
blast_topN () {
    local query_fasta="$1" db="$2" n="$3"
    blastp -query "$query_fasta" -db "$BLASTDB/$db" \
        -outfmt "6 sseqid pident evalue qcovs" \
        -max_target_seqs "$n" -max_hsps 1 -evalue 1e-5 2>>"$LOG"
}

reverse_gene_symbol () {
    local query_fasta="$1"
    blastp -query "$query_fasta" -db "$BLASTDB/$HUMAN_DB" \
        -outfmt "6 stitle" -max_target_seqs 1 -max_hsps 1 -evalue 1e-5 2>>"$LOG" \
        | head -n 1 \
        | grep -oE 'GN=[A-Za-z0-9_.-]+' \
        | head -n 1 \
        | cut -d= -f2 || true
}

# Returns: forward_top_hit, forward_pident, forward_evalue, forward_qcov,
#          reverse_symbol_of_top_hit, reciprocal_match (YES/NO)
# Empty if forward returned no hits at all.
do_blastp_capture () {
    local human_gene="$1" human_fasta="$2" bird_db="$3"
    local fwd_lines top_line top_hit top_pident top_evalue top_qcov
    local bird_id_clean bird_fa rev_symbol top_rev_symbol reciprocal

    fwd_lines=$(blast_topN "$human_fasta" "$bird_db" "$N_FORWARD_HITS")
    [[ -z "$fwd_lines" ]] && return

    # Capture the TOP forward hit (first line) as the definitive forward record
    top_line=$(echo "$fwd_lines" | head -n 1)
    top_hit=$(echo "$top_line" | cut -f1)
    top_pident=$(echo "$top_line" | cut -f2)
    top_evalue=$(echo "$top_line" | cut -f3)
    top_qcov=$(echo "$top_line" | cut -f4)

    # Get the reverse symbol of the TOP hit (for the reverse_symbol column)
    bird_id_clean=$(echo "$top_hit" | sed 's/^ref|//;s/|$//')
    bird_fa=$(mktemp)
    if blastdbcmd -db "$BLASTDB/$bird_db" -entry "$bird_id_clean" > "$bird_fa" 2>>"$LOG"; then
        top_rev_symbol=$(reverse_gene_symbol "$bird_fa")
    else
        top_rev_symbol=""
    fi
    rm -f "$bird_fa"

    # Check if any of top-N reciprocates (sets reciprocal=YES if found)
    reciprocal="NO"
    while IFS=$'\t' read -r bird_id ident evalue qcov; do
        [[ -z "$bird_id" ]] && continue
        bird_id_clean=$(echo "$bird_id" | sed 's/^ref|//;s/|$//')
        bird_fa=$(mktemp)
        if blastdbcmd -db "$BLASTDB/$bird_db" -entry "$bird_id_clean" > "$bird_fa" 2>>"$LOG"; then
            rev_symbol=$(reverse_gene_symbol "$bird_fa")
            if [[ "$rev_symbol" == "$human_gene" ]]; then
                reciprocal="YES"
                rm -f "$bird_fa"
                break
            fi
        fi
        rm -f "$bird_fa"
    done <<< "$fwd_lines"

    echo -e "$top_hit\t$top_pident\t$top_evalue\t$top_qcov\t$top_rev_symbol\t$reciprocal"
}

# tblastn against bird genome: returns top hit "sseqid\tpident\tevalue\tqcov"
# or empty.
do_tblastn () {
    local human_fasta="$1" genome_db="$2"
    tblastn -query "$human_fasta" -db "$BLASTDB/$genome_db" \
        -outfmt "6 sseqid pident evalue qcovs" \
        -max_target_seqs 1 -max_hsps 1 -evalue 1e-5 2>>"$LOG" \
        | head -n 1
}

# ---------- main loop ----------
tail -n +2 "$FLANKERS_TSV" | while IFS=$'\t' read -r human_gene ensembl_id human_fasta; do
    [[ -z "$human_gene" ]] && continue
    echo "[$(date)] $human_gene" | tee -a "$LOG"

    for species in chicken duck_ZJU1.0 duck_T2T; do
        # Tier 1: Compara (skips duck_T2T -- not in Ensembl)
        if [[ -n "${COMPARA_SPECIES[$species]:-}" ]]; then
            fetch_compara "$human_gene" "$ensembl_id" "$species" "${COMPARA_SPECIES[$species]}"
            sleep 0.1
        fi

        # Tier 2: blastp (always run; captures top hit regardless of reciprocity)
        blastp_row=$(do_blastp_capture "$human_gene" "$human_fasta" "${BIRD_PROT_DB[$species]}")

        # Tier 3: tblastn against genome (only if blastp empty)
        tblastn_cols=""
        if [[ -z "$blastp_row" ]]; then
            tblastn_hit=$(do_tblastn "$human_fasta" "${BIRD_GENOME_DB[$species]}")
            if [[ -n "$tblastn_hit" ]]; then
                tblastn_cols="\t$tblastn_hit"
                blastp_row=$(printf '\t\t\t\t\t')  # empty blastp columns
            fi
        fi

        # Write blast.tsv row only if we have any blast evidence (blastp OR tblastn)
        if [[ -n "$blastp_row" ]]; then
            if [[ -n "$tblastn_cols" ]]; then
                echo -e "$human_gene\t$species\t${blastp_row}${tblastn_cols}" >> "$BLAST_OUT"
            else
                # blastp succeeded; tblastn columns left empty
                echo -e "$human_gene\t$species\t${blastp_row}\t\t\t\t" >> "$BLAST_OUT"
            fi
        fi
        # If neither blastp nor tblastn found anything, no row written
        # (downstream will treat as NO_BLAST_EVIDENCE)
    done
done

echo "[$(date)] done" | tee -a "$LOG"
echo "  compara.tsv: $(wc -l < $COMPARA_OUT) lines"
echo "  blast_v3.tsv: $(wc -l < $BLAST_OUT) lines"
