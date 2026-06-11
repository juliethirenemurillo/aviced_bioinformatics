"""
02_uniprot_crosscheck.py
------------------------
Stage 2 — UniProt verification and canonical isoform selection.

For each (gene_symbol, organism) pair in the NCBI Gene search results,
queries UniProt to:
  1. Confirm the gene identity (catches aliasing errors like ZBP1/IGF2BP1)
  2. Retrieve the canonical isoform (SwissProt reviewed > TrEMBL)
  3. Cross-reference to the NCBI RefSeq accession
  4. Flag discrepancies between NCBI Gene search and UniProt

Input:  data/reference/ncbi_gene_search_results.tsv
Output: data/reference/uniprot_verification.tsv

Status codes:
  VERIFIED_SWISSPROT       UniProt reviewed entry found, RefSeq matches
  VERIFIED_TREMBL          UniProt unreviewed entry found
  MISMATCH_ALIAS           NCBI returned a different gene (e.g. ZBP1 → IGF2BP1)
  NO_UNIPROT_ENTRY         no UniProt entry exists
  CONFIRMED_ABSENT         no NCBI hit and no UniProt entry — true absence
"""

import requests
import csv
import os
import time

# ── configuration ──────────────────────────────────────────────────────────────

INPUT_FILE  = "data/reference/ncbi_gene_search_results.tsv"
OUTPUT_FILE = "data/reference/uniprot_verification.tsv"

# UniProt REST API endpoint
UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"

# organism to taxonomy ID mapping (UniProt uses these)
TAXON_IDS = {
    "human":   "9606",    # Homo sapiens
    "chicken": "9031",    # Gallus gallus
    "duck":    "8839",    # Anas platyrhynchos
}

# ── functions ──────────────────────────────────────────────────────────────────

def query_uniprot(gene_symbol, organism_key):
    """
    Query UniProt for a gene symbol in a specific organism.
    Prefers SwissProt (reviewed) over TrEMBL.
    Returns a dict with verification details.
    """
    taxon_id = TAXON_IDS[organism_key]
    
    # query 1: SwissProt (reviewed) first
    params = {
        "query":  f"gene_exact:{gene_symbol} AND organism_id:{taxon_id} AND reviewed:true",
        "format": "json",
        "fields": "accession,id,gene_names,protein_name,length,sequence,xref_refseq",
        "size":   5,
    }
    
    swissprot_hits = _do_query(params)
    
    if swissprot_hits:
        return _format_hit(swissprot_hits[0], "SWISSPROT", gene_symbol, organism_key)
    
    # query 2: TrEMBL (unreviewed) as fallback
    params["query"] = f"gene_exact:{gene_symbol} AND organism_id:{taxon_id} AND reviewed:false"
    trembl_hits = _do_query(params)
    
    if trembl_hits:
        return _format_hit(trembl_hits[0], "TREMBL", gene_symbol, organism_key)
    
    # no hits at all
    return {
        "uniprot_accession":  "",
        "uniprot_id":         "",
        "uniprot_protein":    "",
        "uniprot_gene_names": "",
        "uniprot_length":     "",
        "uniprot_refseq":     "",
        "uniprot_tier":       "NONE",
    }


def _do_query(params):
    """Perform a UniProt REST query and return the results list."""
    try:
        response = requests.get(UNIPROT_API, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()
        return data.get("results", [])
    except Exception as e:
        print(f"    UniProt query error: {e}")
        return []


def _format_hit(hit, tier, gene_symbol, organism_key):
    """Extract the fields we need from a UniProt hit."""
    
    # gene names (catches aliases)
    gene_names = ""
    if "genes" in hit:
        names = []
        for g in hit["genes"]:
            if "geneName" in g:
                names.append(g["geneName"].get("value", ""))
            if "synonyms" in g:
                for syn in g["synonyms"]:
                    names.append(syn.get("value", ""))
        gene_names = ";".join([n for n in names if n])
    
    # protein name
    protein_name = ""
    if "proteinDescription" in hit:
        rec = hit["proteinDescription"].get("recommendedName", {})
        if rec:
            protein_name = rec.get("fullName", {}).get("value", "")
    
    # RefSeq cross-references
    refseq_accs = []
    for xref in hit.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "RefSeq":
            refseq_accs.append(xref.get("id", ""))
    
    # sequence length
    length = hit.get("sequence", {}).get("length", "")
    
    return {
        "uniprot_accession":  hit.get("primaryAccession", ""),
        "uniprot_id":         hit.get("uniProtkbId", ""),
        "uniprot_protein":    protein_name,
        "uniprot_gene_names": gene_names,
        "uniprot_length":     str(length),
        "uniprot_refseq":     ";".join(refseq_accs),
        "uniprot_tier":       tier,
    }


def determine_status(ncbi_row, uniprot_result):
    """
    Determine verification status based on comparison between
    NCBI Gene search and UniProt query.
    """
    ncbi_symbol = ncbi_row["official_symbol"].upper() if ncbi_row["official_symbol"] else ""
    ncbi_status = ncbi_row["status"]
    target_gene = ncbi_row["gene_symbol"].upper()
    uniprot_tier = uniprot_result["uniprot_tier"]
    uniprot_genes = uniprot_result["uniprot_gene_names"].upper()
    
    # case 1: neither NCBI nor UniProt found the gene
    if ncbi_status == "NO_HIT" and uniprot_tier == "NONE":
        return "CONFIRMED_ABSENT"
    
    # case 2: NCBI found nothing but UniProt did — worth investigating
    if ncbi_status == "NO_HIT" and uniprot_tier != "NONE":
        return "UNIPROT_ONLY_CHECK"
    
    # case 3: UniProt has nothing but NCBI does
    if uniprot_tier == "NONE" and ncbi_status == "OK":
        return "NCBI_ONLY_CHECK"
    
    # case 4: NCBI returned a different gene symbol (aliasing error)
    if ncbi_symbol and ncbi_symbol != target_gene:
        # is the target gene in UniProt's gene name list?
        if target_gene in uniprot_genes.split(";"):
            # UniProt has the correct gene — the NCBI row for this hit is wrong
            return "MISMATCH_ALIAS"
    
    # case 5: good agreement
    if uniprot_tier == "SWISSPROT":
        return "VERIFIED_SWISSPROT"
    if uniprot_tier == "TREMBL":
        return "VERIFIED_TREMBL"
    
    return "UNCLEAR"


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("AviCeD — UniProt cross-check (Stage 2)")
    print("=" * 70)
    
    # read NCBI search results
    with open(INPUT_FILE, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        ncbi_rows = list(reader)
    
    print(f"Loaded {len(ncbi_rows)} rows from NCBI search")
    print()
    
    # deduplicate to unique (gene_symbol, organism) pairs
    # we only want to query UniProt once per pair, even if NCBI returned multiple hits
    seen_pairs = set()
    unique_queries = []
    for row in ncbi_rows:
        key = (row["gene_symbol"], row["organism"])
        if key not in seen_pairs:
            seen_pairs.add(key)
            unique_queries.append(row)
    
    # enrich each unique pair with UniProt data
    uniprot_cache = {}
    for i, row in enumerate(unique_queries, 1):
        gene = row["gene_symbol"]
        org  = row["organism"]
        
        print(f"[{i:3d}/{len(unique_queries)}] {gene} ({org})... ", end="", flush=True)
        
        result = query_uniprot(gene, org)
        uniprot_cache[(gene, org)] = result
        
        tier = result["uniprot_tier"]
        if tier == "SWISSPROT":
            print(f"✓ SwissProt {result['uniprot_accession']} ({result['uniprot_length']} aa)")
        elif tier == "TREMBL":
            print(f"~ TrEMBL {result['uniprot_accession']} ({result['uniprot_length']} aa)")
        else:
            print(f"✗ no UniProt entry")
        
        time.sleep(0.4)  # polite delay
    
    # now produce the output file combining NCBI and UniProt data
    output_rows = []
    for row in ncbi_rows:
        key = (row["gene_symbol"], row["organism"])
        uniprot_result = uniprot_cache.get(key, {})
        status = determine_status(row, uniprot_result)
        
        combined = {
            "gene_symbol":       row["gene_symbol"],
            "organism":          row["organism"],
            "ncbi_gene_id":      row["gene_id"],
            "ncbi_symbol":       row["official_symbol"],
            "ncbi_description":  row["description"],
            "ncbi_refseq_count": len(row["refseq_protein"].split(";")) if row["refseq_protein"] else 0,
            "uniprot_accession": uniprot_result.get("uniprot_accession", ""),
            "uniprot_id":        uniprot_result.get("uniprot_id", ""),
            "uniprot_protein":   uniprot_result.get("uniprot_protein", ""),
            "uniprot_length":    uniprot_result.get("uniprot_length", ""),
            "uniprot_refseq":    uniprot_result.get("uniprot_refseq", ""),
            "uniprot_tier":      uniprot_result.get("uniprot_tier", ""),
            "status":            status,
        }
        output_rows.append(combined)
    
    # write it
    fieldnames = [
        "gene_symbol", "organism", "ncbi_gene_id", "ncbi_symbol",
        "ncbi_description", "ncbi_refseq_count",
        "uniprot_accession", "uniprot_id", "uniprot_protein",
        "uniprot_length", "uniprot_refseq", "uniprot_tier", "status",
    ]
    
    with open(OUTPUT_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    
    print()
    print("=" * 70)
    print(f"Output: {OUTPUT_FILE}")
    print("=" * 70)
    
    # summary
    statuses = {}
    for r in output_rows:
        statuses[r["status"]] = statuses.get(r["status"], 0) + 1
    
    print("\nStatus summary:")
    for status, count in sorted(statuses.items()):
        print(f"  {status:25s} {count}")
    
    print("\nCheck especially:")
    print("  MISMATCH_ALIAS      — rows where NCBI returned the wrong gene")
    print("  NCBI_ONLY_CHECK     — NCBI has it, UniProt doesn't")
    print("  UNIPROT_ONLY_CHECK  — UniProt has it, NCBI didn't")


if __name__ == "__main__":
    main()