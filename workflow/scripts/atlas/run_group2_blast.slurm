#!/usr/bin/env bash
#SBATCH -J group2_blast
#SBATCH -A BRYANT-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH -o logs/group2_blast_%j.out
#SBATCH -e logs/group2_blast_%j.err
#
# run_group2_blast.slurm   — reuses the DBs built by run_group1_blast.slurm.
#   cd ~/comparativegenom/aviced/blast_disambiguation
#   sbatch scripts/run_group2_blast.slurm
#
# Phase A: 16 human queries -> blastp (chicken+duck proteins) + tblastn (3 genomes).
# Phase B: for single-species genes, derive the PRESENT species' protein (top hit
#          of the human query there) and search the ABSENT species with it
#          (bird-to-bird = more sensitive). The derived query is RECORDED, not
#          interpreted; verdicts stay pending for joint review.
set -uo pipefail
module load blast/2.17

BASE="$HOME/comparativegenom/aviced/blast_disambiguation"; cd "$BASE"
DBDIR="$BASE/blastdb"; G2="$BASE/results/group2"; QDIR="$BASE/queries"; LOG="${SLURM_JOB_ID:+$BASE/logs/group2_blast_$SLURM_JOB_ID.log}"
mkdir -p "$G2"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
QUERY="$QDIR/group2_human_queries.fasta"
[ -s "$QUERY" ] || { echo "MISSING $QUERY (run fetch_group2_queries.sh first)"; exit 1; }

OF="6 qseqid sseqid pident length qlen slen qcovs evalue bitscore stitle"
HDR="qseqid\tsseqid\tpident\tlength\tqlen\tslen\tqcovs\tevalue\tbitscore\tstitle"

search () {  # $1=prog $2=db $3=query $4=out $5=mask(default|nomask)
  local seg=""; [ "${5:-default}" = nomask ] && seg="-seg no"
  printf "%b\n" "$HDR" > "$4"
  "$1" -query "$3" -db "$2" -num_threads "$THREADS" -evalue 1 $seg -max_target_seqs 50 -outfmt "$OF" >> "$4"
  echo "    $(( $(wc -l < "$4") - 1 )) hits -> $(basename "$4")"
}

echo "[$(date +%T)] PHASE A: human queries (16) vs both species"
search blastp  "$DBDIR/chicken_prot"   "$QUERY" "$G2/blastp_chicken.tsv"
search blastp  "$DBDIR/duck_zju_prot"  "$QUERY" "$G2/blastp_duck_ZJU1.0.tsv"
search tblastn "$DBDIR/chicken_gen"    "$QUERY" "$G2/tblastn_chicken.tsv"
search tblastn "$DBDIR/duck_zju_gen"   "$QUERY" "$G2/tblastn_duck_ZJU1.0.tsv"
search tblastn "$DBDIR/duck_t2t_gen"   "$QUERY" "$G2/tblastn_duck_T2T.tsv"

echo "[$(date +%T)] PHASE B: bird-to-bird for single-species genes"
declare -A ACC=( [CARD8]=Q9Y2G2 [GSDMA]=Q96QA5 [GSDMB]=Q8TAX9 [GSDMC]=Q9BYG8 \
                 [GSDMD]=P57764 [RNF31]=Q96EP0 [BIM]=O43521 [DDX58]=O95786 )
# gene present_species  (present species supplies the query; absent species is searched)
B2B="CARD8 chicken
GSDMA chicken
GSDMB chicken
GSDMC chicken
GSDMD chicken
RNF31 chicken
BIM duck
DDX58 duck"

printf "gene\tpresent_species\tderived_query_acc\tderived_query_desc\ttarget_species\n" > "$G2/bird2bird_manifest.tsv"

while read -r gene pres; do
  [ -z "${gene:-}" ] && continue
  hacc="${ACC[$gene]}"
  if [ "$pres" = chicken ]; then
    srcbp="$G2/blastp_chicken.tsv"; srcdb="$DBDIR/chicken_prot"
    tgt_bp="$DBDIR/duck_zju_prot"; tgt_tn="$DBDIR/duck_zju_gen $DBDIR/duck_t2t_gen"; tlab="duck_ZJU1.0 duck_T2T"; tspec="duck"
  else
    srcbp="$G2/blastp_duck_ZJU1.0.tsv"; srcdb="$DBDIR/duck_zju_prot"
    tgt_bp="$DBDIR/chicken_prot"; tgt_tn="$DBDIR/chicken_gen"; tlab="chicken"; tspec="chicken"
  fi
  row=$(grep "$hacc" "$srcbp" | sort -t$'\t' -k9,9gr | head -1)
  if [ -z "$row" ]; then echo "  $gene: no $pres hit to derive a bird query — skipped"; continue; fi
  sseqid=$(cut -f2 <<<"$row"); desc=$(cut -f10 <<<"$row")
  qacc=$(grep -oE '[A-Z]P_[0-9]+\.[0-9]+' <<<"$sseqid" | head -1)
  if [ -z "$qacc" ]; then echo "  $gene: could not parse accession from '$sseqid' — skipped"; continue; fi
  q="$QDIR/bird_${gene}_${qacc}.faa"
  if ! blastdbcmd -db "$srcdb" -entry "$qacc" > "$q" 2>/dev/null; then
    echo "  $gene: blastdbcmd extract failed for $qacc — skipped"; continue; fi
  printf "%s\t%s\t%s\t%s\t%s\n" "$gene" "$pres" "$qacc" "$desc" "$tspec" >> "$G2/bird2bird_manifest.tsv"
  echo "  $gene: $pres query $qacc ($desc) -> search $tspec"
  search blastp "$tgt_bp" "$q" "$G2/bird2bird_${gene}_blastp_${tspec}.tsv" nomask
  set -- $tgt_tn; i=1
  for db in $tgt_tn; do
    lab=$(echo $tlab | cut -d' ' -f$i)
    search tblastn "$db" "$q" "$G2/bird2bird_${gene}_tblastn_${lab}.tsv" nomask
    i=$((i+1))
  done
done <<< "$B2B"

echo "[$(date +%T)] best-hit-per-query summary (Phase A)"
SUM="$G2/group2_besthit_summary.tsv"
printf "search\tqseqid\tpident\tqcov\tevalue\tstitle\n" > "$SUM"
for f in blastp_chicken blastp_duck_ZJU1.0 tblastn_chicken tblastn_duck_ZJU1.0 tblastn_duck_T2T; do
  tail -n +2 "$G2/$f.tsv" | sort -t$'\t' -k1,1 -k9,9gr | awk -F'\t' '!s[$1]++{print FN"\t"$1"\t"$3"\t"$7"\t"$8"\t"$10}' FN="$f" >> "$SUM"
done
echo "  -> $(basename "$SUM")"
echo "[$(date +%T)] done. results in $G2/"
