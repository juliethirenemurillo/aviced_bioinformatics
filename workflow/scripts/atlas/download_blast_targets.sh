#!/usr/bin/env bash
# download_blast_targets.sh
# Run on the CSD3 LOGIN node (needs internet; compute nodes have none).
# Downloads BLAST targets straight from the NCBI FTP site (no 'datasets' CLI,
# which isn't installed on CSD3), verifies md5, decompresses, and logs
# accession + assembly name + date for provenance.
#
#   blastp target  = RefSeq protein set  (_protein.faa)  -> proteins/
#   tblastn target = genome assembly     (_genomic.fna)  -> genomes/
#
# Targets:
#   chicken GRCg7b      GCF_016699485.2   protein + genome
#   duck    ZJU1.0      GCF_015476345.1   protein + genome   (primary)
#   duck    T2T         GCF_047663525.1   genome only        (no-hit confirmation)
set -uo pipefail

BASE="${BASE:-$HOME/comparativegenom/aviced/blast_disambiguation}"
PROT="$BASE/proteins"; GEN="$BASE/genomes"; LOGDIR="$BASE/logs"
mkdir -p "$PROT" "$GEN" "$LOGDIR"
LOG="$LOGDIR/download_targets_$(date +%Y%m%d_%H%M%S).log"
PROV="$BASE/genomes/DOWNLOAD_PROVENANCE.tsv"
FTP="https://ftp.ncbi.nlm.nih.gov/genomes/all"
fail=0
printf "target\taccession\tassembly_name\tassembly_date\tfile\tkind\tmd5\tdownloaded\n" > "$PROV"

# accession -> NCBI FTP parent path, e.g. GCF_016699485.2 -> GCF/016/699/485
acc_path () {
  local n="${1#GCF_}"; n="${n%.*}"
  echo "GCF/${n:0:3}/${n:3:3}/${n:6:3}"
}

verify_md5 () {  # $1=local file  $2=md5checksums file  $3=remote basename
  local exp got
  exp=$(grep -E "[[:space:]]\./?${3}\$" "$2" 2>/dev/null | awk '{print $1}' | head -1)
  [ -z "$exp" ] && exp=$(grep -F "$3" "$2" 2>/dev/null | awk '{print $1}' | head -1)
  got=$(md5sum "$1" | awk '{print $1}')
  [ -n "$exp" ] && [ "$exp" = "$got" ]
}

# $1=label  $2=accession  $3=want_protein(0/1)  $4=want_genome(0/1)
fetch_assembly () {
  local label="$1" acc="$2" wp="$3" wg="$4"
  local parent dir base report md5f asm date
  parent="$FTP/$(acc_path "$acc")"
  echo "[$(date +%T)] $label  $acc" | tee -a "$LOG"

  # resolve the assembly subdirectory name (contains the assembly name we don't hardcode)
  dir=$(curl -fsSL "$parent/" 2>>"$LOG" | grep -oE "${acc}_[._A-Za-z0-9-]+/" | head -1 | tr -d '/')
  if [ -z "$dir" ]; then
    echo "   !! could not resolve directory for $acc under $parent" | tee -a "$LOG"; fail=1; return
  fi
  base="$parent/$dir"
  echo "   dir: $dir" | tee -a "$LOG"

  # provenance: assembly name + date from the assembly report
  report="$GEN/${dir}_assembly_report.txt"
  curl -fsSL --retry 3 --retry-delay 3 "$base/${dir}_assembly_report.txt" -o "$report" 2>>"$LOG"
  asm=$(grep -m1 "Assembly name:"      "$report" 2>/dev/null | sed 's/.*: *//' | tr -d '\r')
  date=$(grep -m1 "Date:"              "$report" 2>/dev/null | sed 's/.*: *//' | tr -d '\r')
  : "${asm:=$dir}"; : "${date:=unknown}"

  # md5 manifest for integrity checking
  md5f="$GEN/${dir}_md5checksums.txt"
  curl -fsSL --retry 3 --retry-delay 3 "$base/md5checksums.txt" -o "$md5f" 2>>"$LOG"

  # download one file type, verify, decompress, record provenance
  grab () {  # $1=suffix (protein.faa / genomic.fna)  $2=destdir  $3=kind
    local rb="${dir}_$1.gz" dest="$2/${dir}_$1.gz"
    echo "   fetch $1.gz" | tee -a "$LOG"
    if curl -fSL --retry 3 --retry-delay 3 "$base/$rb" -o "$dest" 2>>"$LOG"; then
      local ok="md5_skip"
      if verify_md5 "$dest" "$md5f" "$rb"; then ok="md5_ok"; else ok="md5_FAIL"; fail=1; fi
      gunzip -f "$dest"
      local plain="$2/${dir}_$1"
      local n; n=$(grep -c '^>' "$plain" 2>/dev/null || echo "?")
      echo "      $ok  sequences=$n  -> $plain" | tee -a "$LOG"
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$label" "$acc" "$asm" "$date" "${dir}_$1" "$3" "$ok" "$(date +%F)" >> "$PROV"
    else
      echo "      !! download failed: $rb" | tee -a "$LOG"; fail=1
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$label" "$acc" "$asm" "$date" "${dir}_$1" "$3" "DL_FAIL" "$(date +%F)" >> "$PROV"
    fi
  }

  [ "$wp" = 1 ] && grab "protein.faa" "$PROT" "protein"
  [ "$wg" = 1 ] && grab "genomic.fna" "$GEN"  "genome"
}

fetch_assembly "chicken_GRCg7b" GCF_016699485.2 1 1
fetch_assembly "duck_ZJU1.0"    GCF_015476345.1 1 1
fetch_assembly "duck_T2T"       GCF_047663525.1 0 1   # genome only, for no-hit confirmation

echo "----" | tee -a "$LOG"
echo "provenance:" | tee -a "$LOG"; column -t -s$'\t' "$PROV" 2>/dev/null | tee -a "$LOG" || cat "$PROV"
echo "proteins/:"; ls -lh "$PROT"/*.faa 2>/dev/null
echo "genomes/:";  ls -lh "$GEN"/*.fna  2>/dev/null
if [ "$fail" -eq 0 ]; then echo "RESULT: ALL TARGETS OK"; else echo "RESULT: CHECK NEEDED — see $LOG and $PROV"; exit 1; fi
