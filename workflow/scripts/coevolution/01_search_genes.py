"""
01_search_genes.py
------------------
Stage 1 of sequence curation pipeline.

Takes a list of gene symbols and searches NCBI Gene for each one across
three organisms (human, chicken, duck). Does NOT fetch sequences —
only returns identifying information for manual verification.

Output: data/reference/ncbi_gene_search_results.tsv

Columns:
  gene_symbol | organism | gene_id | official_symbol | description |
  chromosome | refseq_protein | status

After running, YOU inspect this table, verify each row, and mark
which entries to use for the sequence fetch in stage 2.
"""

from Bio import Entrez
import time
import csv
import os
import xml.etree.ElementTree as ET

# ── configuration ──────────────────────────────────────────────────────────────

Entrez.email = "your_email@cam.ac.uk"  # UPDATE THIS

OUTPUT_DIR = "data/reference"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "ncbi_gene_search_results.tsv")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# target organisms with their NCBI taxonomy names
ORGANISMS = {
    "human":   "Homo sapiens",
    "chicken": "Gallus gallus",
    "duck":    "Anas platyrhynchos",
}

# simple gene symbol list — only symbols, no accessions assumed
TARGET_GENES = [
    # pyroptosis
    "CASP1", "CASP3", "CASP4", "CASP5", "CASP7", "CASP8", "CASP9",
    "GSDMD", "GSDMA", "GSDME",
    "PYCARD",  # also known as ASC
    "NLRP3", "CARD8", "NAIP", "NLRC4", "AIM2",
    # necroptosis
    "RIPK1", "RIPK3", "ZBP1", "MLKL",
    # PRR adaptors
    "TICAM1",  # TRIF
    "TICAM2",  # TRAM
    "MYD88", "TRADD", "FADD",
    # viral sensing
    "DDX58",   # RIG-I
    "IFIH1",   # MDA5
    "STING1",  # TMEM173
    # lytic effector
    "NINJ1",
]

# ── functions ──────────────────────────────────────────────────────────────────

def search_gene_in_organism(gene_symbol, organism_full_name):
    """
    Search NCBI Gene database for a gene symbol in a specific organism.
    Returns the Gene ID(s) if found, or None.
    """
    query = f"{gene_symbol}[gene symbol] AND {organism_full_name}[orgn]"
    
    try:
        handle = Entrez.esearch(db="gene", term=query, retmax=5)
        record = Entrez.read(handle)
        handle.close()
        
        gene_ids = record.get("IdList", [])
        return gene_ids
    
    except Exception as e:
        print(f"  ✗ SEARCH ERROR for {gene_symbol} in {organism_full_name}: {e}")
        return []


def get_gene_details(gene_id):
    """
    Retrieve full details for a Gene ID.
    Returns a dict with: symbol, description, chromosome, refseq_proteins
    """
    try:
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="gene_table", retmode="xml")
        xml_data = handle.read()
        handle.close()
        
        root = ET.fromstring(xml_data)
        
        # extract the key fields from the XML
        details = {
            "gene_id": gene_id,
            "official_symbol": "",
            "description": "",
            "chromosome": "",
            "refseq_proteins": [],
        }
        
        # official symbol
        symbol_elem = root.find(".//Gene-ref_locus")
        if symbol_elem is not None:
            details["official_symbol"] = symbol_elem.text or ""
        
        # description
        desc_elem = root.find(".//Gene-ref_desc")
        if desc_elem is not None:
            details["description"] = desc_elem.text or ""
        
        # chromosome
        chrom_elem = root.find(".//Gene-source_subtype[@value='chromosome']/..")
        if chrom_elem is not None:
            name = chrom_elem.find("Gene-source_subtype_name")
            if name is not None:
                details["chromosome"] = name.text or ""
        
        # refseq proteins — look for products with NP_ or XP_ prefixes
        for product in root.findall(".//Gene-commentary"):
            accession_elem = product.find(".//Gene-commentary_accession")
            if accession_elem is not None and accession_elem.text:
                acc = accession_elem.text
                if acc.startswith(("NP_", "XP_")):
                    if acc not in details["refseq_proteins"]:
                        details["refseq_proteins"].append(acc)
        
        return details
    
    except Exception as e:
        print(f"  ✗ FETCH ERROR for gene ID {gene_id}: {e}")
        return None


def process_gene(gene_symbol, organism_key, organism_name):
    """Process a single gene-organism pair. Returns list of result rows."""
    print(f"\n→ {gene_symbol} in {organism_name}")
    
    gene_ids = search_gene_in_organism(gene_symbol, organism_name)
    time.sleep(0.4)  # be kind to NCBI
    
    if not gene_ids:
        print(f"  — no Gene record found")
        return [{
            "gene_symbol": gene_symbol,
            "organism": organism_key,
            "gene_id": "",
            "official_symbol": "",
            "description": "",
            "chromosome": "",
            "refseq_protein": "",
            "status": "NO_HIT",
        }]
    
    rows = []
    for gene_id in gene_ids:
        details = get_gene_details(gene_id)
        time.sleep(0.4)
        
        if details is None:
            continue
        
        refseq_list = details["refseq_proteins"]
        refseq_str = ";".join(refseq_list) if refseq_list else ""
        
        status = "OK"
        if len(gene_ids) > 1:
            status = "MULTIPLE_HITS"
        if not refseq_list:
            status = "NO_REFSEQ_PROTEIN"
        
        row = {
            "gene_symbol":     gene_symbol,
            "organism":        organism_key,
            "gene_id":         details["gene_id"],
            "official_symbol": details["official_symbol"],
            "description":     details["description"],
            "chromosome":      details["chromosome"],
            "refseq_protein":  refseq_str,
            "status":          status,
        }
        rows.append(row)
        
        print(f"  ✓ Gene ID {details['gene_id']} | {details['official_symbol']} | {details['description'][:60]}")
        if refseq_list:
            print(f"    RefSeq: {', '.join(refseq_list[:3])}{' ...' if len(refseq_list) > 3 else ''}")
    
    return rows


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("AviCeD — NCBI Gene search (Stage 1)")
    print("=" * 70)
    print(f"Genes:      {len(TARGET_GENES)}")
    print(f"Organisms:  {', '.join(ORGANISMS.values())}")
    print(f"Output:     {OUTPUT_FILE}")
    print("=" * 70)
    
    all_rows = []
    
    for gene_symbol in TARGET_GENES:
        for organism_key, organism_name in ORGANISMS.items():
            rows = process_gene(gene_symbol, organism_key, organism_name)
            all_rows.extend(rows)
    
    # write the TSV
    with open(OUTPUT_FILE, "w", newline="") as f:
        fieldnames = [
            "gene_symbol", "organism", "gene_id", "official_symbol",
            "description", "chromosome", "refseq_protein", "status",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)
    
    print("\n" + "=" * 70)
    print(f"Done. {len(all_rows)} rows written to {OUTPUT_FILE}")
    print("=" * 70)
    print("\nNEXT STEP: Open the TSV and manually verify each row.")
    print("Pay attention to:")
    print("  - NO_HIT       → confirms absence (cross-check with literature)")
    print("  - MULTIPLE_HITS → ambiguity to resolve manually")
    print("  - description  → confirms it is actually the right gene")


if __name__ == "__main__":
    main()