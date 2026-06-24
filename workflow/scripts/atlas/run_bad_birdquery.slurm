#!/usr/bin/env bash
#SBATCH -J bad_birdquery
#SBATCH -A BRYANT-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:15:00
#SBATCH --mem=8G
#SBATCH -o logs/bad_birdquery_%j.out
#SBATCH -e logs/bad_birdquery_%j.err
#
# run_bad_birdquery.slurm
# Sensitivity check for BAD. The human BAD query found nothing in duck, but BAD
# is a short, fast-evolving BH3-only protein, so a distant human query is a weak
# test. Here we query the CHICKEN BAD protein (bird-to-bird is far more sensitive)
# against the duck genomes. Chicken genome = positive control (must hit).
# Reuses the nucleotide DBs built by run_group1_blast.slurm (no rebuild).
#   cd ~/comparativegenom/aviced/blast_disambiguation
#   sbatch scripts/run_bad_birdquery.slurm
set -euo pipefail
module load blast/2.17

BASE="$HOME/comparativegenom/aviced/blast_disambiguation"
cd "$BASE"
DBDIR="$BASE/blastdb"; RES="$BASE/results"; QDIR="$BASE/queries"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
CHK_PROT="$BASE/proteins/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_protein.faa"
ACC="XP_046762038.1"                         # chicken BAD, from the blastp result
QY="$QDIR/chicken_BAD_${ACC}.faa"

# extract the chicken BAD record from the local protein FASTA
awk -v acc="$ACC" '/^>/{p=($0 ~ acc)} p' "$CHK_PROT" > "$QY"
n=$(grep -c '^>' "$QY" || true)
[ "$n" -eq 1 ] || { echo "expected 1 sequence for $ACC, got $n"; exit 1; }
echo "query: $(head -1 "$QY")"

OF="6 qseqid sseqid pident length qlen slen qcovs evalue bitscore stitle"
HDR="qseqid\tsseqid\tpident\tlength\tqlen\tslen\tqcovs\tevalue\tbitscore\tstitle"

# deliberately maximal sensitivity for a short query:
#   -evalue 10  permissive ceiling (we WANT to see any trace)
#   -seg no     don't low-complexity-mask the short query
run () {  # $1=db  $2=out
  printf "%b\n" "$HDR" > "$2"
  tblastn -query "$QY" -db "$1" -num_threads "$THREADS" \
          -evalue 10 -seg no -max_target_seqs 50 -outfmt "$OF" >> "$2"
  echo "  $(( $(wc -l < "$2") - 1 )) hit rows -> $2"
}

echo "[$(date +%T)] chicken BAD -> chicken genome (positive control)"
run "$DBDIR/chicken_gen"  "$RES/tblastn_BAD_chickenq_chicken.tsv"
echo "[$(date +%T)] chicken BAD -> duck ZJU1.0"
run "$DBDIR/duck_zju_gen" "$RES/tblastn_BAD_chickenq_duck_ZJU1.0.tsv"
echo "[$(date +%T)] chicken BAD -> duck T2T"
run "$DBDIR/duck_t2t_gen" "$RES/tblastn_BAD_chickenq_duck_T2T.tsv"
echo "[$(date +%T)] done."
