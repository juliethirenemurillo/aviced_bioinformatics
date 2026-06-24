#!/usr/bin/env bash
#SBATCH -J group34_blast
#SBATCH -A BRYANT-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=12G
#SBATCH -o logs/group34_blast_%j.out
#SBATCH -e logs/group34_blast_%j.err
#
# run_group34_blast.slurm   — reuses DBs from run_group1_blast.slurm.
#   cd ~/comparativegenom/aviced/blast_disambiguation
#   sbatch scripts/run_group34_blast.slurm     (or: bash scripts/... on login node)
#
# G4 (NLRC4, TNFRSF10A/B): standard human-query pass (confirm-and-dismiss).
# G3 (CASP1/4/5, TRIF): keep ALL hits, add ppos (positives %), dump the caspase
#   inventory per species, and run bird-to-bird with the annotated chicken
#   caspase-1 protein to test the ~14%-vs-published-~45% CASP1 contradiction.
set -uo pipefail
module load blast/2.17
BASE="$HOME/comparativegenom/aviced/blast_disambiguation"; cd "$BASE"
DBDIR="$BASE/blastdb"; G="$BASE/results/group34"; QDIR="$BASE/queries"; mkdir -p "$G"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
QUERY="$QDIR/group34_human_queries.fasta"
[ -s "$QUERY" ] || { echo "MISSING $QUERY (run fetch_group34_queries.sh first)"; exit 1; }
CHK_PROT="$BASE/proteins/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_protein.faa"
DUCK_PROT="$BASE/proteins/GCF_015476345.1_ZJU1.0_protein.faa"

OF="6 qseqid sseqid pident ppos length qlen slen qcovs evalue bitscore stitle"
HDR="qseqid\tsseqid\tpident\tppos\tlength\tqlen\tslen\tqcovs\tevalue\tbitscore\tstitle"

search () {  # $1=prog $2=db $3=query $4=out $5=mask(default|nomask)
  local seg=""; [ "${5:-default}" = nomask ] && seg="-seg no"
  printf "%b\n" "$HDR" > "$4"
  "$1" -query "$3" -db "$2" -num_threads "$THREADS" -evalue 1 $seg -max_target_seqs 100 -outfmt "$OF" >> "$4"
  echo "    $(( $(wc -l < "$4") - 1 )) hits -> $(basename "$4")"
}

echo "[$(date +%T)] PHASE A: 7 human queries vs both species (all hits, with ppos)"
search blastp  "$DBDIR/chicken_prot"   "$QUERY" "$G/blastp_chicken.tsv"
search blastp  "$DBDIR/duck_zju_prot"  "$QUERY" "$G/blastp_duck_ZJU1.0.tsv"
search tblastn "$DBDIR/chicken_gen"    "$QUERY" "$G/tblastn_chicken.tsv"
search tblastn "$DBDIR/duck_zju_gen"   "$QUERY" "$G/tblastn_duck_ZJU1.0.tsv"
search tblastn "$DBDIR/duck_t2t_gen"   "$QUERY" "$G/tblastn_duck_T2T.tsv"

echo "[$(date +%T)] caspase-cluster hit views (human caspase queries vs bird caspases)"
for f in blastp_chicken blastp_duck_ZJU1.0; do
  { head -1 "$G/$f.tsv"; awk -F'\t' 'NR>1 && tolower($11) ~ /caspase/' "$G/$f.tsv"; } > "$G/caspasehits_$f.tsv"
  echo "    $(( $(wc -l < "$G/caspasehits_$f.tsv") - 1 )) caspase hits -> caspasehits_$f.tsv"
done

echo "[$(date +%T)] annotated caspase inventory per species"
grep -iE '^>.*caspase'  "$CHK_PROT"  | sed 's/^>//' > "$G/caspase_inventory_chicken.txt"
grep -iE '^>.*caspase'  "$DUCK_PROT" | sed 's/^>//' > "$G/caspase_inventory_duck.txt"
echo "    chicken: $(wc -l < "$G/caspase_inventory_chicken.txt")  duck: $(wc -l < "$G/caspase_inventory_duck.txt") caspase proteins"

echo "[$(date +%T)] PHASE B: bird-to-bird with annotated chicken caspase-1"
c1hdr=$(grep -iE '^>.*caspase-1([^0-9]|$)' "$CHK_PROT" | head -1)
if [ -n "$c1hdr" ]; then
  c1acc=$(sed 's/^>//; s/ .*//' <<<"$c1hdr")
  echo "    chicken caspase-1 = $c1acc  ($(sed 's/^[^ ]* //' <<<"$c1hdr"))"
  q="$QDIR/chicken_CASP1_${c1acc}.faa"
  awk -v acc="$c1acc" '/^>/{p=($0 ~ acc)} p' "$CHK_PROT" > "$q"
  printf "%s\t%s\n" "chicken_CASP1_query" "$c1acc" > "$G/bird2bird_CASP1_query.txt"
  echo "$c1hdr" >> "$G/bird2bird_CASP1_query.txt"
  search blastp  "$DBDIR/duck_zju_prot" "$q" "$G/bird2bird_CASP1_blastp_duck.tsv" nomask
  search tblastn "$DBDIR/duck_zju_gen"  "$q" "$G/bird2bird_CASP1_tblastn_duck_ZJU1.0.tsv" nomask
  search tblastn "$DBDIR/duck_t2t_gen"  "$q" "$G/bird2bird_CASP1_tblastn_duck_T2T.tsv" nomask
else
  echo "    !! no protein annotated 'caspase-1' in chicken set — inspect caspase_inventory_chicken.txt by hand"
fi

echo "[$(date +%T)] best-hit-per-query summary (Phase A, with ppos)"
SUM="$G/group34_besthit_summary.tsv"
printf "search\tqseqid\tpident\tppos\tqcov\tevalue\tstitle\n" > "$SUM"
for f in blastp_chicken blastp_duck_ZJU1.0 tblastn_chicken tblastn_duck_ZJU1.0 tblastn_duck_T2T; do
  tail -n +2 "$G/$f.tsv" | sort -t$'\t' -k1,1 -k10,10gr | awk -F'\t' '!s[$1]++{print FN"\t"$1"\t"$3"\t"$4"\t"$8"\t"$9"\t"$11}' FN="$f" >> "$SUM"
done
echo "    -> $(basename "$SUM")"
echo "[$(date +%T)] done. results in $G/"
