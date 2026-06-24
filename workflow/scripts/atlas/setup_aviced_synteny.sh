#!/bin/bash
# setup_aviced_synteny.sh -- one-time prerequisites for AviCeD Paper 1 synteny step
#
# Builds:
#   $BLASTDB/human_protein.{phr,pin,psq,...}    (UniProt SwissProt human, reviewed)
#   $BLASTDB/duck_T2T_protein.{phr,pin,psq,...}  (NCBI RefSeq T2T duck proteins)
#
# Fetches GFF3 annotations into $GFFDIR:
#   human.gff           (GRCh38.p14)
#   chicken_GRCg7b.gff  (GCF_016699485.2)
#   duck_ZJU1.0.gff     (GCF_015476345.1)
#   duck_T2T.gff        (GCF_047663525.1, annotation release RS_2025_03)
#
# Touches nothing else: the chicken/duck-ZJU protein+genome blastdbs that
# already exist in $BLASTDB are left alone.
#
# Idempotent: each step skips if its output already exists. Safe to re-run
# after a failure -- it picks up where it stopped.
#
# Requires (in the active conda env): ncbi-datasets-cli, blast, wget, unzip, jq.

set -euo pipefail

# ============ paths ============
BASE="$HOME/comparativegenom/aviced/blast_disambiguation"
BLASTDB="$BASE/blastdb"
GFFDIR="$BASE/gff"
WORKDIR="$BASE/setup_tmp"

mkdir -p "$BLASTDB" "$GFFDIR" "$WORKDIR"

# ============ assembly accessions ============
HUMAN_ASSEMBLY="GCF_000001405.40"   # GRCh38.p14
CHICKEN_ASSEMBLY="GCF_016699485.2"  # GRCg7b
DUCK_ZJU_ASSEMBLY="GCF_015476345.1" # ZJU1.0
DUCK_T2T_ASSEMBLY="GCF_047663525.1" # T2T, RS_2025_03 annotation

# ============ helpers ============
log() { echo "[$(date '+%H:%M:%S')] $*"; }

require_tool() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "ERROR: missing required tool '$1' on PATH"
        echo "       inside the aviced env, run:"
        echo "       mamba install -c conda-forge -c bioconda $1"
        exit 1
    fi
}

for t in datasets unzip makeblastdb blastdbcmd wget jq; do
    require_tool "$t"
done

log "setup starting; BASE=$BASE"

# ============================================================
# STEP 1 -- human SwissProt -> $BLASTDB/human_protein
# ============================================================
HUMAN_FA="$WORKDIR/human_swissprot.fasta"
HUMAN_DB="$BLASTDB/human_protein"

if [[ -f "${HUMAN_DB}.phr" && -f "${HUMAN_DB}.pin" && -f "${HUMAN_DB}.psq" ]]; then
    log "STEP 1: human_protein blastdb already exists -- skipping"
else
    log "STEP 1: downloading human SwissProt from UniProt..."
    # organism_id 9606 = Homo sapiens; reviewed:true = SwissProt (canonical, ~20k)
    URL='https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=organism_id%3A9606+AND+reviewed%3Atrue'
    wget --quiet --show-progress -O "${HUMAN_FA}.gz" "$URL"
    gunzip -f "${HUMAN_FA}.gz"

    NSEQ=$(grep -c '^>' "$HUMAN_FA")
    log "STEP 1: downloaded $NSEQ sequences"
    # sanity: SwissProt human is ~20,400; warn if wildly off
    if [[ "$NSEQ" -lt 15000 || "$NSEQ" -gt 30000 ]]; then
        echo "WARNING: $NSEQ sequences is outside the expected 15k-30k range"
        echo "         check the UniProt URL or query syntax before continuing"
    fi

    log "STEP 1: building human_protein blastdb..."
    makeblastdb -in "$HUMAN_FA" -dbtype prot -out "$HUMAN_DB" \
        -title "Human SwissProt (UniProt UP000005640 reviewed)" \
        -parse_seqids
    log "STEP 1: human_protein blastdb built"
fi

# ============================================================
# STEP 2 -- fetch GFF3 + protein set for each assembly
# ============================================================
# Helper: download an assembly bundle (gff3 + protein), copy gff to GFFDIR
# and protein.faa to WORKDIR (consumed by step 3 for T2T only).
fetch_assembly() {
    local acc="$1" label="$2"
    local zipfile="$WORKDIR/${label}.zip"
    local extractdir="$WORKDIR/${label}_data"
    local gff_target="$GFFDIR/${label}.gff"
    local prot_target="$WORKDIR/${label}_protein.faa"

    if [[ -f "$gff_target" ]]; then
        log "  ${label}: GFF already present (${gff_target}) -- skipping fetch"
        # but ensure protein.faa exists in WORKDIR if previous run cleaned tmp
        if [[ ! -f "$prot_target" && -d "$extractdir" ]]; then
            local prot_src
            prot_src=$(find "$extractdir" -name 'protein.faa' 2>/dev/null | head -1 || true)
            [[ -n "$prot_src" ]] && cp "$prot_src" "$prot_target"
        fi
        return
    fi

    log "  ${label}: downloading ${acc} (gff3 + protein)..."
    rm -rf "$extractdir"
    datasets download genome accession "$acc" \
        --include gff3,protein \
        --filename "$zipfile" 2>&1 | tail -2

    log "  ${label}: extracting..."
    unzip -q -o "$zipfile" -d "$extractdir"

    local gff_src
    gff_src=$(find "$extractdir" -name 'genomic.gff' | head -1)
    if [[ -z "$gff_src" ]]; then
        echo "ERROR: no genomic.gff found in ${extractdir}"
        return 1
    fi
    cp "$gff_src" "$gff_target"

    local prot_src
    prot_src=$(find "$extractdir" -name 'protein.faa' | head -1 || true)
    if [[ -n "$prot_src" ]]; then
        cp "$prot_src" "$prot_target"
    fi

    local ngenes
    ngenes=$(awk -F'\t' '$3=="gene"' "$gff_target" | wc -l)
    log "  ${label}: ${ngenes} gene features -> ${gff_target}"
    rm -f "$zipfile"
}

log "STEP 2: fetching GFFs and protein sets"
fetch_assembly "$HUMAN_ASSEMBLY"    "human"
fetch_assembly "$CHICKEN_ASSEMBLY"  "chicken_GRCg7b"
fetch_assembly "$DUCK_ZJU_ASSEMBLY" "duck_ZJU1.0"
fetch_assembly "$DUCK_T2T_ASSEMBLY" "duck_T2T"

# ============================================================
# STEP 3 -- duck_T2T proteins -> $BLASTDB/duck_T2T_protein
# ============================================================
T2T_PROT_FA="$WORKDIR/duck_T2T_protein.faa"
T2T_DB="$BLASTDB/duck_T2T_protein"

if [[ -f "${T2T_DB}.phr" && -f "${T2T_DB}.pin" && -f "${T2T_DB}.psq" ]]; then
    log "STEP 3: duck_T2T_protein blastdb already exists -- skipping"
elif [[ -f "$T2T_PROT_FA" ]]; then
    log "STEP 3: building duck_T2T_protein blastdb..."
    NSEQ=$(grep -c '^>' "$T2T_PROT_FA")
    log "STEP 3: $NSEQ T2T proteins to index"
    makeblastdb -in "$T2T_PROT_FA" -dbtype prot -out "$T2T_DB" \
        -title "Duck T2T (GCF_047663525.1) RefSeq proteins, RS_2025_03" \
        -parse_seqids
    log "STEP 3: duck_T2T_protein blastdb built"
else
    echo "ERROR: ${T2T_PROT_FA} not found -- step 2 should have produced it"
    exit 1
fi

# ============================================================
# VERIFY
# ============================================================
log "VERIFY:"
log "  blastdbs in $BLASTDB:"
blastdbcmd -list "$BLASTDB" -list_outfmt "    %f  %p  %t" || true
log "  GFFs in $GFFDIR:"
ls -lh "$GFFDIR" | awk 'NR>1 {printf "    %s  %s\n",$5,$9}'

log "setup complete"
log ""
log "next: run build_flanker_inputs.sh once the flankers.tsv is prepared"
