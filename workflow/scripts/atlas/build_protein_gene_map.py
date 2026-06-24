#!/usr/bin/env python
"""
build_protein_gene_map.py  --  AviCeD synteny step

Builds a map: protein_accession -> GFF_gene_Name  for a given GFF3 file.

Chain:
  CDS features have  Name=XP_xxx  and  Dbxref=...,GeneID:YYY
  gene features have Name=SYMBOL  and  Dbxref=...,GeneID:YYY
  Join on GeneID -> protein_accession -> GFF_gene_Name

Output TSV: protein_accession  gene_name  seqid

Usage:
  python build_protein_gene_map.py input.gff > protein_gene_map.tsv

The output is consumed by synteny_call.py to resolve BLAST hit accessions
(ref|XP_xxx|) to the actual Name= value in the GFF, handling LOC-annotated
genes that wouldn't be found by human gene symbol alone.
"""
import argparse, csv, sys


def extract_geneid(dbxref):
    """Find GeneID:XXXXX in a comma-separated Dbxref value."""
    for entry in dbxref.split(','):
        if entry.strip().startswith('GeneID:'):
            return entry.strip()[7:]
    return None


def extract_protein_id(attrs):
    """Extract protein accession from CDS attributes.
    Tries Name= first (standard for NCBI RefSeq GFFs), then protein_id=."""
    name = attrs.get('Name', '')
    if name and (name.startswith('XP_') or name.startswith('NP_')):
        return name
    pid = attrs.get('protein_id', '')
    if pid and (pid.startswith('XP_') or pid.startswith('NP_')):
        return pid
    return None


def parse_attrs(s):
    d = {}
    for kv in s.split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k] = v
    return d


def build_map(gff_path):
    """Build protein_accession -> (gene_name, seqid) via GeneID join."""
    # pass 1: gene features -> GeneID -> (gene_name, seqid)
    geneid_to_gene = {}
    # pass 2: CDS features -> protein_accession -> GeneID
    prot_to_geneid = {}

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 9:
                continue
            ftype = c[2]
            attrs = parse_attrs(c[8])
            dbxref = attrs.get('Dbxref', '')
            gid = extract_geneid(dbxref)

            if ftype == 'gene' and gid:
                name = attrs.get('Name', '')
                seqid = c[0]
                geneid_to_gene[gid] = (name, seqid)

            elif ftype == 'CDS' and gid:
                prot = extract_protein_id(attrs)
                if prot:
                    prot_to_geneid.setdefault(prot, gid)

    # join: protein -> geneid -> gene_name
    result = {}
    for prot, gid in prot_to_geneid.items():
        if gid in geneid_to_gene:
            name, seqid = geneid_to_gene[gid]
            result[prot] = (name, seqid)
    return result


def main():
    p = argparse.ArgumentParser(
        description='Build protein_accession -> GFF gene Name map from GFF3')
    p.add_argument('gff', help='GFF3 file path')
    a = p.parse_args()

    mapping = build_map(a.gff)
    w = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
    w.writerow(['protein_accession', 'gene_name', 'seqid'])
    for prot, (name, seqid) in sorted(mapping.items()):
        w.writerow([prot, name, seqid])

    print(f"# {len(mapping)} protein -> gene mappings", file=sys.stderr)


if __name__ == '__main__':
    main()
