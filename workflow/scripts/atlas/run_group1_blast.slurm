#!/usr/bin/env bash
#SBATCH -J group1_blast
#SBATCH -A BRYANT-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH -o logs/group1_blast_%j.out
#SBATCH -e logs/group1_blast_%j.err
#
# run_group1_blast.slurm
# Build BLAST databases and run blastp + tblastn for the 8 Group-1 human queries
# against chicken GRCg7b, duck ZJU1.0 (primary), and duck T2T (confirmation).
# Submit from the blast_disambiguation directory so logs/ resolves:
#   cd ~/comparativegenom/aviced/blast_disambiguation
#   sbatch scripts/run_group1_blast.slurm
#
# Raw, exploratory output: permissive -evalue is a REPORTING ceiling, NOT the
# present/absent decision threshold. Thresholds get set together after review.
set -euo pipefail

module load blast/2.17

BASE="$HOME/comparativegenom/aviced/blast_disambiguation"
cd "$BASE"
DBDIR="$BASE/blastdb"; RES="$BASE/results"; mkdir -p "$DBDIR" "$RES"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
QUERY="$BASE/queries/group1_human_queries.fasta"

CHK_PROT="$BASE/proteins/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_protein.faa"
DUCK_PROT="$BASE/proteins/GCF_015476345.1_ZJU1.0_protein.faa"
CHK_GEN="$BASE/genomes/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna"
DUCK_GEN="$BASE/genomes/GCF_015476345.1_ZJU1.0_genomic.fna"
T2T_GEN="$BASE/genomes/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.fna"

# sanity: every input present before we start
for f in "$QUERY" "$CHK_PROT" "$DUCK_PROT" "$CHK_GEN" "$DUCK_GEN" "$T2T_GEN"; do
  [ -s "$f" ] || { echo "MISSING INPUT: $f"; exit 1; }
done
echo "[$(date +%T)] inputs OK; threads=$THREADS"

# build a db only if its index is absent (idempotent: safe to resubmit)
build () {  # $1=fasta  $2=prot|nucl  $3=db-out  $4=title
  local probe; [ "$2" = prot ] && probe="$3.pin" || probe="$3.nin"
  if [ -s "$probe" ]; then echo "  db exists: $3"; return; fi
  echo "  makeblastdb ($2): $3"
  makeblastdb -in "$1" -dbtype "$2" -parse_seqids -out "$3" -title "$4"
}

echo "[$(date +%T)] building databases"
build "$CHK_PROT"  prot "$DBDIR/chicken_prot"   chicken_GRCg7b_prot
build "$DUCK_PROT" prot "$DBDIR/duck_zju_prot"  duck_ZJU1.0_prot
build "$CHK_GEN"   nucl "$DBDIR/chicken_gen"    chicken_GRCg7b_gen
build "$DUCK_GEN"  nucl "$DBDIR/duck_zju_gen"   duck_ZJU1.0_gen
build "$T2T_GEN"   nucl "$DBDIR/duck_t2t_gen"   duck_T2T_gen

OF="6 qseqid sseqid pident length qlen slen qcovs evalue bitscore stitle"
HDR="qseqid\tsseqid\tpident\tlength\tqlen\tslen\tqcovs\tevalue\tbitscore\tstitle"

run () {  # $1=blastp|tblastn  $2=db  $3=out.tsv
  printf "%b\n" "$HDR" > "$3"
  "$1" -query "$QUERY" -db "$2" -num_threads "$THREADS" \
       -evalue 1 -max_target_seqs 50 -outfmt "$OF" >> "$3"
  echo "  $(( $(wc -l < "$3") - 1 )) hit rows -> $3"
}

echo "[$(date +%T)] blastp (8 queries vs protein sets)"
run blastp  "$DBDIR/chicken_prot"  "$RES/blastp_chicken.tsv"
run blastp  "$DBDIR/duck_zju_prot" "$RES/blastp_duck_ZJU1.0.tsv"

echo "[$(date +%T)] tblastn (8 queries vs genomes)"
run tblastn "$DBDIR/chicken_gen"   "$RES/tblastn_chicken.tsv"
run tblastn "$DBDIR/duck_zju_gen"  "$RES/tblastn_duck_ZJU1.0.tsv"
run tblastn "$DBDIR/duck_t2t_gen"  "$RES/tblastn_duck_T2T.tsv"

echo "[$(date +%T)] done. results in $RES/"
