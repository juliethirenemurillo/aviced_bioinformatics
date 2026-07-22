#!/usr/bin/env python
"""
00_fetch_sequences.py
=====================
Collect RIPK1 and MLKL protein + CDS sequences across vertebrate species
for the avian cell-death selection analysis (AviCeD Paper 1).

What this script does
---------------------
For each species in species_list.tsv, it uses the NCBI datasets CLI to
find the orthologue of the human gene (RIPK1 or MLKL). "Orthologue" here
means a gene in another species that descended from the same ancestral gene
— NCBI maintains a curated database of these relationships (Gene Orthologs).

The script downloads the protein sequence and the CDS (coding DNA sequence)
for each orthologue. When a gene has multiple transcript variants, it picks
the longest protein — this is a rough proxy for the "canonical" isoform,
which is usually the full-length one.

Output: multi-species FASTA files ready for alignment with MAFFT.

Requirements
------------
    conda activate aviced
    # ncbi-datasets-cli is already in the aviced env

Usage
-----
    python 00_fetch_sequences.py species_list.tsv data/

    First argument:  path to the species list (TSV, with columns
                     common_name, species, ncbi_taxid, group, notes)
    Second argument: output directory (created if it doesn't exist)

Output files
------------
    data/RIPK1_protein.fasta    multi-species protein alignment input
    data/RIPK1_cds.fasta        multi-species CDS (for codon back-translation)
    data/MLKL_protein.fasta
    data/MLKL_cds.fasta
    data/fetch_summary.tsv      one row per species x gene: what was found
"""

import subprocess
import json
import os
import sys
import csv
import zipfile
import tempfile
from pathlib import Path

# ---------------------------------------------------------------------------
# Human gene IDs — the anchors for the orthology lookup.
#
# The datasets CLI --ortholog flag says: "find me the orthologue of this
# human gene in taxon X." So we start from the human Gene ID and project
# outward to each target species.
# ---------------------------------------------------------------------------
GENES = {
    "RIPK1": 8737,    # NCBI Gene ID for human RIPK1
    "MLKL":  197259,  # NCBI Gene ID for human MLKL
}


# ---------------------------------------------------------------------------
# FASTA helpers
# ---------------------------------------------------------------------------

def parse_fasta(path):
    """
    Read a FASTA file. Returns a list of (header, sequence) tuples.

    A FASTA file is a text format where each sequence starts with a line
    beginning with '>', followed by the header/description, then one or more
    lines of sequence letters. Example:

        >NP_003792.2 RIPK1 [Homo sapiens]
        MSDAAQFPQP...
        KLLENDGDVI...
    """
    seqs = []
    header = ""
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    seqs.append((header, "".join(seq_lines)))
                header = line[1:]  # strip the '>'
                seq_lines = []
            elif line:
                seq_lines.append(line)
    if header:
        seqs.append((header, "".join(seq_lines)))
    return seqs


def write_fasta(path, records):
    """
    Write a list of (header, sequence) tuples as a FASTA file.
    Wraps sequences at 80 characters per line (standard convention).
    """
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


# ---------------------------------------------------------------------------
# Sequence fetching
# ---------------------------------------------------------------------------

def fetch_orthologue(human_gene_id, species_name, taxid, tmpdir):
    """
    Fetch the orthologue of a human gene in a target species.

    Uses the NCBI datasets CLI with the --ortholog flag, which queries
    NCBI's Gene Orthologs database. This is more reliable than searching
    by gene symbol (which can differ across species — e.g. a bird gene
    might be annotated as LOC12345 rather than RIPK1).

    Parameters
    ----------
    human_gene_id : int
        NCBI Gene ID for the human gene (e.g. 8737 for RIPK1)
    species_name : str
        Binomial species name (e.g. "Gallus gallus")
    taxid : str or int
        NCBI Taxonomy ID for the species (e.g. 9031 for chicken)
    tmpdir : str
        Temporary directory for downloads

    Returns
    -------
    dict with keys:
        found (bool), protein_id, protein_seq, cds_seq, gene_symbol, method
    """
    result = {
        "found": False,
        "protein_id": "",
        "protein_seq": "",
        "cds_seq": "",
        "gene_symbol": "",
        "method": "",
    }

    zippath = os.path.join(tmpdir, f"gene_{human_gene_id}_{taxid}.zip")
    extractdir = os.path.join(tmpdir, f"gene_{human_gene_id}_{taxid}")

    # ------------------------------------------------------------------
    # Strategy 1: --ortholog (preferred — uses curated orthology)
    # ------------------------------------------------------------------
    cmd = [
        "datasets", "download", "gene", "gene-id", str(human_gene_id),
        "--ortholog", str(taxid),
        "--include", "cds,protein",
        "--filename", zippath,
    ]

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=120
        )

        if proc.returncode != 0 or not os.path.exists(zippath):
            # Strategy 1 failed — no orthologue found in this taxon
            return result

        # Unzip the download
        os.makedirs(extractdir, exist_ok=True)
        with zipfile.ZipFile(zippath, "r") as zf:
            zf.extractall(extractdir)

        # The download contains a data report (JSONL) with gene metadata
        report_path = os.path.join(
            extractdir, "ncbi_dataset", "data", "data_report.jsonl"
        )
        if not os.path.exists(report_path):
            return result

        # Read gene metadata
        with open(report_path) as f:
            for line in f:
                entry = json.loads(line.strip())
                result["gene_symbol"] = entry.get("symbol", "unknown")
                break  # take the first entry

        result["method"] = "ortholog"

        # Extract protein sequences
        protein_path = os.path.join(
            extractdir, "ncbi_dataset", "data", "protein.faa"
        )
        if os.path.exists(protein_path):
            seqs = parse_fasta(protein_path)
            if seqs:
                # Pick the longest protein as the canonical isoform
                longest = max(seqs, key=lambda x: len(x[1]))
                result["protein_id"] = longest[0].split()[0]
                result["protein_seq"] = longest[1]

        # Extract CDS sequences
        cds_path = os.path.join(
            extractdir, "ncbi_dataset", "data", "cds.fna"
        )
        if os.path.exists(cds_path):
            seqs = parse_fasta(cds_path)
            if seqs:
                longest = max(seqs, key=lambda x: len(x[1]))
                result["cds_seq"] = longest[1]

        if result["protein_seq"]:
            result["found"] = True

    except subprocess.TimeoutExpired:
        print("  TIMEOUT", file=sys.stderr)
    except Exception as e:
        print(f"  WARNING: {e}", file=sys.stderr)

    return result


def fetch_human_gene(gene_id, tmpdir):
    """
    Fetch the human gene directly (not as orthologue of itself).
    Same logic, just without the --ortholog flag.
    """
    result = {
        "found": False,
        "protein_id": "",
        "protein_seq": "",
        "cds_seq": "",
        "gene_symbol": "",
        "method": "direct",
    }

    zippath = os.path.join(tmpdir, f"gene_{gene_id}_human.zip")
    extractdir = os.path.join(tmpdir, f"gene_{gene_id}_human")

    cmd = [
        "datasets", "download", "gene", "gene-id", str(gene_id),
        "--include", "cds,protein",
        "--filename", zippath,
    ]

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=120
        )
        if proc.returncode != 0 or not os.path.exists(zippath):
            return result

        os.makedirs(extractdir, exist_ok=True)
        with zipfile.ZipFile(zippath, "r") as zf:
            zf.extractall(extractdir)

        report_path = os.path.join(
            extractdir, "ncbi_dataset", "data", "data_report.jsonl"
        )
        if os.path.exists(report_path):
            with open(report_path) as f:
                for line in f:
                    entry = json.loads(line.strip())
                    result["gene_symbol"] = entry.get("symbol", "unknown")
                    break

        protein_path = os.path.join(
            extractdir, "ncbi_dataset", "data", "protein.faa"
        )
        if os.path.exists(protein_path):
            seqs = parse_fasta(protein_path)
            if seqs:
                longest = max(seqs, key=lambda x: len(x[1]))
                result["protein_id"] = longest[0].split()[0]
                result["protein_seq"] = longest[1]

        cds_path = os.path.join(
            extractdir, "ncbi_dataset", "data", "cds.fna"
        )
        if os.path.exists(cds_path):
            seqs = parse_fasta(cds_path)
            if seqs:
                longest = max(seqs, key=lambda x: len(x[1]))
                result["cds_seq"] = longest[1]

        if result["protein_seq"]:
            result["found"] = True

    except Exception as e:
        print(f"  WARNING: {e}", file=sys.stderr)

    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 3:
        print("Usage: python 00_fetch_sequences.py species_list.tsv outdir/")
        sys.exit(1)

    species_file = sys.argv[1]
    outdir = sys.argv[2]

    # Read species list
    species_list = []
    with open(species_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            species_list.append(row)

    print(f"Loaded {len(species_list)} species from {species_file}")
    os.makedirs(outdir, exist_ok=True)

    # Create a temp directory for raw downloads
    tmpdir = tempfile.mkdtemp(prefix="aviced_fetch_")
    print(f"Temp directory: {tmpdir}")

    summary_rows = []
    collected = {gene: {"protein": [], "cds": []} for gene in GENES}

    for gene_name, human_gene_id in GENES.items():
        print(f"\n{'='*60}")
        print(f"  {gene_name}  (Human Gene ID {human_gene_id})")
        print(f"{'='*60}\n")

        for sp in species_list:
            common   = sp["common_name"]
            binomial = sp["species"]
            taxid    = sp["ncbi_taxid"]
            group    = sp["group"]

            tag = f"  {common:<16s} ({binomial})"
            print(f"{tag}  ...", end=" ", flush=True)

            # Fetch
            if binomial == "Homo sapiens":
                res = fetch_human_gene(human_gene_id, tmpdir)
            else:
                res = fetch_orthologue(human_gene_id, binomial, taxid, tmpdir)

            # Collect
            if res["found"]:
                label = common  # clean label for FASTA header
                prot_len = len(res["protein_seq"])
                cds_len  = len(res["cds_seq"])

                collected[gene_name]["protein"].append(
                    (f"{label}", res["protein_seq"])
                )
                if res["cds_seq"]:
                    collected[gene_name]["cds"].append(
                        (f"{label}", res["cds_seq"])
                    )

                print(
                    f"OK  protein {prot_len:>5d} aa  "
                    f"CDS {cds_len:>6d} nt  "
                    f"symbol={res['gene_symbol']}"
                )
            else:
                prot_len = 0
                cds_len  = 0
                print("NOT FOUND")

            summary_rows.append({
                "gene":        gene_name,
                "common_name": common,
                "species":     binomial,
                "taxid":       taxid,
                "group":       group,
                "found":       res["found"],
                "protein_id":  res.get("protein_id", ""),
                "protein_len": prot_len,
                "cds_len":     cds_len,
                "gene_symbol": res.get("gene_symbol", ""),
                "method":      res.get("method", ""),
            })

    # ------------------------------------------------------------------
    # Write output files
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print("  OUTPUT FILES")
    print(f"{'='*60}\n")

    for gene_name in GENES:
        n_prot = len(collected[gene_name]["protein"])
        n_cds  = len(collected[gene_name]["cds"])

        if n_prot:
            p = os.path.join(outdir, f"{gene_name}_protein.fasta")
            write_fasta(p, collected[gene_name]["protein"])
            print(f"  {p}  ({n_prot} sequences)")

        if n_cds:
            p = os.path.join(outdir, f"{gene_name}_cds.fasta")
            write_fasta(p, collected[gene_name]["cds"])
            print(f"  {p}  ({n_cds} sequences)")

    # Summary table
    summary_path = os.path.join(outdir, "fetch_summary.tsv")
    fields = [
        "gene", "common_name", "species", "taxid", "group",
        "found", "protein_id", "protein_len", "cds_len",
        "gene_symbol", "method",
    ]
    with open(summary_path, "w") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"  {summary_path}")

    # ------------------------------------------------------------------
    # Summary to screen
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print("  SUMMARY")
    print(f"{'='*60}\n")

    for gene_name in GENES:
        gene_rows = [r for r in summary_rows if r["gene"] == gene_name]
        found   = sum(1 for r in gene_rows if r["found"])
        total   = len(gene_rows)
        missing = [r["common_name"] for r in gene_rows if not r["found"]]

        print(f"  {gene_name}: {found}/{total} species")
        if missing:
            print(f"    missing: {', '.join(missing)}")

    print(f"\nDone. Temp files in {tmpdir} (can be deleted).\n")


if __name__ == "__main__":
    main()
