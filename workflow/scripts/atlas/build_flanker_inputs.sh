#!/bin/bash
# build_flanker_inputs.sh  --  AviCeD synteny step (v2 refactor)
#
# Reciprocal-best-BLAST (RBB) refactor:
#   v1 took only the top forward hit and tested it for reciprocity.
#   That fails systematically when the top hit is a close paralog -- e.g.
#   human RUVBL2 BLASTed against chicken returns RUVBL1 as the top hit
#   (43% identity), reciprocity fails (RUVBL1 -> human RUVBL1), and the
#   real RUVBL2 ortholog (a few percent lower identity) was never checked.
#
# v2 strategy: pull the top N forward hits, then for each in identity
# order do the reverse BLAST. Report the FIRST forward hit whose reverse
# top hit maps back to the query gene. This is the canonical RBB definition
# and correctly recovers the true ortholog even when paralogs dominate.
#
# Inputs:
#   FLANKERS_TSV  -- TSV with cols: human_gene, human_ensembl, human_protein_fasta
# Outputs:
#   compara.tsv  -- human_gene, species, bird_gene, ortholog_type, identity
#   blast.tsv    -- human_gene, species, bird_gene, identity, evalue, qcov
#                   (now correctly captures non-top-hit reciprocal best)
#
# Reciprocity check: matches the reverse top hit's "GN=" gene symbol against
# the query gene name. Requires human_protein blastdb to be the SwissProt
# proteome (SwissProt FASTA headers contain GN=SYMBOL).
#
# Tunables:
#   N_FORWARD_HITS  -- how many top forward hits to test (default 10)
#                     Larger = catches more diverged orthologs but slower.
#
# Conda env: aviced (blastp, curl, jq).

set -euo pipefail

FLANKERS_TSV="${1:-flankers.tsv}"
OUTDIR="${2:-flanker_inputs}"
BLASTDB="${BLASTDB:-$HOME/comparativegenom/aviced/blast_disambiguation/blastdb}"
N_FORWARD_HITS="${N_FORWARD_HITS:-10}"

mkdir -p "$OUTDIR"
COMPARA_OUT="$OUTDIR/compara.tsv"
BLAST_OUT="$OUTDIR/blast.tsv"
LOG="$OUTDIR/build.log"

declare -A COMPARA_SPECIES=(
    [chicken]="gallus_gallus"
    [duck_ZJU1.0]="anas_platyrhynchos_platyrhynchos"
)

declare -A BIRD_DB=(
    [chicken]="chicken_prot"
    [duck_ZJU1.0]="duck_zju_prot"
    [duck_T2T]="duck_T2T_protein"
)
HUMAN_DB="human_protein"

echo -e "human_gene\tspecies\tbird_gene\tortholog_type\tidentity" > "$COMPARA_OUT"
echo -e "human_gene\tspecies\tbird_gene\tidentity\tevalue\tqcov"  > "$BLAST_OUT"
: > "$LOG"

echo "[$(date)] starting flanker input build (v2, top-${N_FORWARD_HITS} RBB)" | tee -a "$LOG"

# ---------- Compara fetch (unchanged) ----------
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
# Forward BLAST: returns top N hits as one "sseqid\tpident\tevalue\tqcov" per line.
blast_topN () {
    local query_fasta="$1" db="$2" n="$3"
    blastp -query "$query_fasta" -db "$BLASTDB/$db" \
        -outfmt "6 sseqid pident evalue qcovs" \
        -max_target_seqs "$n" -max_hsps 1 -evalue 1e-5 2>>"$LOG"
}

# Reverse BLAST: returns top 1 hit's "GN=" gene symbol from SwissProt title.
# Returns empty if no hit or no GN= field.
reverse_gene_symbol () {
    local query_fasta="$1"
    blastp -query "$query_fasta" -db "$BLASTDB/$HUMAN_DB" \
        -outfmt "6 stitle" -max_target_seqs 1 -max_hsps 1 -evalue 1e-5 2>>"$LOG" \
        | head -n 1 \
        | grep -oE 'GN=[A-Za-z0-9_.-]+' \
        | head -n 1 \
        | cut -d= -f2 || true
}

# Reciprocal best blast (v2): iterate top-N forward hits, return FIRST whose
# reverse top hit's gene symbol matches the query gene.
# Returns "bird_gene\tident\tevalue\tqcov" or empty.
do_rbb () {
    local human_gene="$1" human_fasta="$2" bird_db="$3"
    local fwd_lines bird_id bird_id_clean ident evalue qcov bird_fa rev_symbol

    fwd_lines=$(blast_topN "$human_fasta" "$bird_db" "$N_FORWARD_HITS")
    [[ -z "$fwd_lines" ]] && return

    # iterate forward hits in identity order
    while IFS=$'\t' read -r bird_id ident evalue qcov; do
        [[ -z "$bird_id" ]] && continue
        # strip "ref|...|" wrapper for blastdbcmd lookup
        bird_id_clean=$(echo "$bird_id" | sed 's/^ref|//;s/|$//')

        # extract this bird hit as fasta for the reverse query
        bird_fa=$(mktemp)
        if ! blastdbcmd -db "$BLASTDB/$bird_db" -entry "$bird_id_clean" > "$bird_fa" 2>>"$LOG"; then
            rm -f "$bird_fa"
            continue
        fi

        # reverse: get top human hit's GN= symbol
        rev_symbol=$(reverse_gene_symbol "$bird_fa")
        rm -f "$bird_fa"

        # reciprocity test
        if [[ "$rev_symbol" == "$human_gene" ]]; then
            echo -e "$bird_id\t$ident\t$evalue\t$qcov"
            return
        fi

        echo "[rbb] $human_gene -> $bird_id (rev=$rev_symbol, no match, trying next)" >> "$LOG"
    done <<< "$fwd_lines"

    # no forward hit reciprocated -- record the failure context
    echo "[rbb] $human_gene -> no reciprocal hit in top-$N_FORWARD_HITS for $bird_db" >> "$LOG"
}

# ---------- main loop ----------
tail -n +2 "$FLANKERS_TSV" | while IFS=$'\t' read -r human_gene ensembl_id human_fasta; do
    [[ -z "$human_gene" ]] && continue
    echo "[$(date)] $human_gene" | tee -a "$LOG"

    for species in chicken duck_ZJU1.0 duck_T2T; do
        if [[ -n "${COMPARA_SPECIES[$species]:-}" ]]; then
            fetch_compara "$human_gene" "$ensembl_id" "$species" "${COMPARA_SPECIES[$species]}"
            sleep 0.1
        fi

        rbb=$(do_rbb "$human_gene" "$human_fasta" "${BIRD_DB[$species]}")
        if [[ -n "$rbb" ]]; then
            echo -e "$human_gene\t$species\t$rbb" >> "$BLAST_OUT"
        fi
    done
done

echo "[$(date)] done" | tee -a "$LOG"
wc -l "$COMPARA_OUT" "$BLAST_OUT"
