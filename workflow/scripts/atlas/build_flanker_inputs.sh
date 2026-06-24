#!/bin/bash
# build_flanker_inputs.sh -- run on CSD3 to produce the inputs that
# flanker_orthologs.py (module B) consumes.
#
# Produces:
#   compara.tsv  -- human_gene, species, bird_gene, ortholog_type, identity
#   blast.tsv    -- human_gene, species, bird_gene, identity, evalue, qcov
#
# Inputs you provide:
#   FLANKERS_TSV  -- TSV with cols: human_gene, human_ensembl, human_protein_fasta
#                    one row per flanker gene to query (typically ~80 flankers:
#                    8 target genes x ~10 flankers each, deduplicated).
#
# Assumes the layout you already have on CSD3:
#   ~/comparativegenom/aviced/blast_disambiguation/blastdb/
#       chicken_protein, duck_ZJU_protein, duck_T2T_protein  (note: T2T protein
#       db must exist for RBB to T2T; if it's currently only a genome db, build
#       the protein db once -- see TODO at bottom)
#   ~/comparativegenom/aviced/blast_disambiguation/blastdb/human_protein
#       (for the reciprocal direction; build once if missing)
#
# Conda env: `conda activate aviced` (BLAST+, curl, jq).
#
# Run once. ~80 flankers x 3 species x 2 BLAST directions = ~480 blastp calls.
# On the login node, expect ~10-20 min total. Single squeue check at the end
# if you submit it as a job; loops over genes have no slurm calls.

set -euo pipefail

FLANKERS_TSV="${1:-flankers.tsv}"
OUTDIR="${2:-flanker_inputs}"
BLASTDB="${BLASTDB:-$HOME/comparativegenom/aviced/blast_disambiguation/blastdb}"

mkdir -p "$OUTDIR"
COMPARA_OUT="$OUTDIR/compara.tsv"
BLAST_OUT="$OUTDIR/blast.tsv"
LOG="$OUTDIR/build.log"

# Compara species labels (Ensembl REST naming)
declare -A COMPARA_SPECIES=(
    [chicken]="gallus_gallus"
    [duck_ZJU1.0]="anas_platyrhynchos_platyrhynchos"
    # duck_T2T: not in Ensembl -> Compara skipped, RBB will cover it
)

# Protein DB names for BLAST (these must exist under $BLASTDB)
declare -A BIRD_DB=(
    [chicken]="chicken_protein"
    [duck_ZJU1.0]="duck_ZJU_protein"
    [duck_T2T]="duck_T2T_protein"
)
HUMAN_DB="human_protein"

echo -e "human_gene\tspecies\tbird_gene\tortholog_type\tidentity" > "$COMPARA_OUT"
echo -e "human_gene\tspecies\tbird_gene\tidentity\tevalue\tqcov"  > "$BLAST_OUT"
: > "$LOG"

echo "[$(date)] starting flanker input build" | tee -a "$LOG"

# ---------- Compara fetch (Ensembl REST) ----------
fetch_compara () {
    local human_gene="$1" ensembl_id="$2" species_label="$3" ensembl_species="$4"
    # /homology/id/{species}/{id} returns all orthologs; filter to target species
    local url="https://rest.ensembl.org/homology/id/homo_sapiens/${ensembl_id}?target_species=${ensembl_species};type=orthologues;format=condensed"
    local resp
    if ! resp=$(curl -s -H 'Accept: application/json' "$url" 2>>"$LOG"); then
        echo "[compara] curl failed: $human_gene $species_label" >> "$LOG"
        return
    fi
    # Parse with jq: take all orthologs in target species
    # Fields: target.id (bird Ensembl ID), type (one2one/one2many/etc), target.perc_id
    echo "$resp" | jq -r --arg sp "$ensembl_species" --arg g "$human_gene" --arg s "$species_label" '
        .data[0].homologies[]?
        | select(.species == $sp)
        | [$g, $s, (.id // ""), (.type // ""), ((.target.perc_id // 0)|tostring)]
        | @tsv
    ' >> "$COMPARA_OUT" 2>>"$LOG" || true
}

# ---------- Reciprocal BLAST ----------
# blastp helper: returns top hit as "subject_id<TAB>pident<TAB>evalue<TAB>qcov"
blast_top () {
    local query_fasta="$1" db="$2"
    blastp -query "$query_fasta" -db "$BLASTDB/$db" \
        -outfmt "6 sseqid pident evalue qcovs" \
        -max_target_seqs 1 -max_hsps 1 -evalue 1e-5 2>>"$LOG" | head -n 1
}

# reciprocal best blast: returns "bird_gene\tident\tevalue\tqcov" if reciprocal,
# else empty.
do_rbb () {
    local human_gene="$1" human_fasta="$2" bird_db="$3"
    local fwd bird_id ident evalue qcov
    fwd=$(blast_top "$human_fasta" "$bird_db")
    [[ -z "$fwd" ]] && return
    bird_id=$(echo "$fwd" | cut -f1)
    ident=$(echo  "$fwd" | cut -f2)
    evalue=$(echo "$fwd" | cut -f3)
    qcov=$(echo   "$fwd" | cut -f4)

    # extract bird protein as fasta for reverse query
    local bird_fa
    bird_fa=$(mktemp)
    blastdbcmd -db "$BLASTDB/$bird_db" -entry "$bird_id" > "$bird_fa" 2>>"$LOG" || { rm -f "$bird_fa"; return; }

    local rev rev_id
    rev=$(blast_top "$bird_fa" "$HUMAN_DB")
    rm -f "$bird_fa"
    [[ -z "$rev" ]] && return
    rev_id=$(echo "$rev" | cut -f1)

    # reciprocity check: does the reverse top hit map back to the query gene?
    # We match by substring against the query header / Ensembl ID; you may want
    # to refine this with a human accession -> gene name map.
    if [[ "$rev_id" == *"$human_gene"* ]] || [[ "$human_fasta" == *"$rev_id"* ]]; then
        echo -e "$bird_id\t$ident\t$evalue\t$qcov"
    fi
}

# ---------- main loop ----------
tail -n +2 "$FLANKERS_TSV" | while IFS=$'\t' read -r human_gene ensembl_id human_fasta; do
    [[ -z "$human_gene" ]] && continue
    echo "[$(date)] $human_gene" | tee -a "$LOG"

    for species in chicken duck_ZJU1.0 duck_T2T; do
        # Compara (skips duck_T2T -- not in Ensembl)
        if [[ -n "${COMPARA_SPECIES[$species]:-}" ]]; then
            fetch_compara "$human_gene" "$ensembl_id" "$species" "${COMPARA_SPECIES[$species]}"
            sleep 0.1  # be polite to Ensembl REST
        fi

        # Reciprocal BLAST
        rbb=$(do_rbb "$human_gene" "$human_fasta" "${BIRD_DB[$species]}")
        if [[ -n "$rbb" ]]; then
            echo -e "$human_gene\t$species\t$rbb" >> "$BLAST_OUT"
        fi
    done
done

echo "[$(date)] done" | tee -a "$LOG"
wc -l "$COMPARA_OUT" "$BLAST_OUT"

# TODO -- one-off setup (do once, before first run):
#   1. Build human protein blastdb if not present:
#        makeblastdb -in human_proteome.fa -dbtype prot -out $BLASTDB/human_protein
#   2. Build duck T2T protein blastdb if you only have the genome db:
#        makeblastdb -in duck_T2T_proteins.fa -dbtype prot -out $BLASTDB/duck_T2T_protein
#   3. The reciprocity check uses substring match on rev_id ~ human_gene; this
#      is robust if your human protein FASTA headers contain the gene symbol.
#      Otherwise build a small accession->gene map and match against that.
