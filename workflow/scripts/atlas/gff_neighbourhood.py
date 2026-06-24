#!/usr/bin/env python
"""
gff_neighbourhood.py  --  AviCeD synteny step (Layer B: GFF-based)

Reconstructs gene order from a GFF3 annotation file and extracts local
neighbourhoods, for synteny-based ortholog/absence calling.

WHY: BLAST left 8 genes ambiguous (paralog-only or no hit). Genomic POSITION
is conserved even when sequence diverges past BLAST recognition. By anchoring
on a gene's conserved flanking neighbours and reading what sits between them in
the bird genome, we distinguish:
  (1) conserved locus + diverged candidate in slot -> PRESENT (diverged)
  (2) conserved locus + EMPTY slot                 -> TRUE biological absence
  (3) locus not recoverable / neighbourhood broken  -> UNRESOLVED (assembly/annotation)

Two modes:
  neighbours  : list the k genes flanking a named anchor gene (build human reference)
  between     : list genes located between two named flanker genes (empty-slot test)

Reads plain or gzipped GFF3. No network, no dependencies beyond stdlib.

Terms:
  GFF3      tab-delimited genome annotation; one row per feature.
  feature   col3, e.g. 'gene','mRNA','exon'. We use 'gene' rows only.
  attributes col9, semicolon-separated key=value, holds the gene name.
"""
import argparse, gzip, sys
from dataclasses import dataclass


@dataclass
class Gene:
    name: str
    seqid: str        # chromosome/scaffold
    start: int
    end: int
    strand: str

    def __repr__(self):
        return f"{self.name}({self.seqid}:{self.start}-{self.end}{self.strand})"


def _open(path):
    """Open plain or gzipped text."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _name_from_attrs(attrs, name_keys):
    """Pull a gene name out of the GFF3 attributes field.
    Different annotators use different keys (Name=, gene=, gene_id=, locus_tag=),
    so we try a priority list and take the first that hits."""
    fields = {}
    for kv in attrs.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            fields[k.strip()] = v.strip()
    for k in name_keys:
        if k in fields and fields[k]:
            return fields[k]
    return None


def load_genes(path, name_keys=("Name", "gene", "gene_id", "locus_tag", "ID")):
    """Parse all gene-type rows into Gene objects."""
    genes = []
    with _open(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            name = _name_from_attrs(c[8], name_keys)
            if name is None:
                continue
            try:
                genes.append(Gene(name, c[0], int(c[3]), int(c[4]), c[6]))
            except ValueError:
                continue
    return genes


def order_by_position(genes):
    """Return {seqid: [genes sorted by start]} -- reconstructs gene order."""
    by_seq = {}
    for g in genes:
        by_seq.setdefault(g.seqid, []).append(g)
    for seqid in by_seq:
        by_seq[seqid].sort(key=lambda g: g.start)
    return by_seq


def _find(by_seq, target, case_insensitive=True):
    """Locate a gene by name; returns (seqid, index_in_ordered_list) or None."""
    t = target.lower() if case_insensitive else target
    for seqid, glist in by_seq.items():
        for i, g in enumerate(glist):
            gn = g.name.lower() if case_insensitive else g.name
            if gn == t:
                return seqid, i
    return None


def neighbours(by_seq, anchor, k=3):
    """k genes either side of `anchor`, in genomic order. Anchor marked with *."""
    hit = _find(by_seq, anchor)
    if hit is None:
        return None
    seqid, i = hit
    glist = by_seq[seqid]
    lo, hi = max(0, i - k), min(len(glist), i + k + 1)
    out = []
    for j in range(lo, hi):
        tag = "*" if j == i else " "
        out.append((tag, glist[j]))
    return seqid, out


def between(by_seq, flank_a, flank_b):
    """Genes strictly between two flankers on the same seqid (empty-slot test).
    Returns (seqid, [genes_between]) or a reason string if it can't be resolved."""
    ha, hb = _find(by_seq, flank_a), _find(by_seq, flank_b)
    if ha is None and hb is None:
        return "BOTH_FLANKERS_ABSENT (locus not recoverable -> UNRESOLVED)"
    if ha is None:
        return f"FLANKER_ABSENT:{flank_a} (-> UNRESOLVED / try other flank or T2T)"
    if hb is None:
        return f"FLANKER_ABSENT:{flank_b} (-> UNRESOLVED / try other flank or T2T)"
    (sa, ia), (sb, ib) = ha, hb
    if sa != sb:
        return f"FLANKERS_ON_DIFFERENT_SEQIDS ({sa} vs {sb} -> broken locus / UNRESOLVED)"
    glist = by_seq[sa]
    lo, hi = sorted((ia, ib))
    inner = glist[lo + 1:hi]
    return sa, inner


def main(argv=None):
    p = argparse.ArgumentParser(description="GFF3 neighbourhood extractor (AviCeD synteny)")
    p.add_argument("gff", help="GFF3 file (.gff/.gff3, optionally .gz)")
    sub = p.add_subparsers(dest="mode", required=True)

    pn = sub.add_parser("neighbours", help="k genes flanking an anchor gene")
    pn.add_argument("anchor")
    pn.add_argument("-k", type=int, default=3)

    pb = sub.add_parser("between", help="genes between two flanker genes")
    pb.add_argument("flank_a")
    pb.add_argument("flank_b")

    a = p.parse_args(argv)
    by_seq = order_by_position(load_genes(a.gff))

    if a.mode == "neighbours":
        res = neighbours(by_seq, a.anchor, a.k)
        if res is None:
            print(f"ANCHOR_NOT_FOUND: {a.anchor}")
            return 1
        seqid, rows = res
        print(f"# neighbourhood of {a.anchor} on {seqid} (k={a.k})")
        for tag, g in rows:
            print(f"{tag} {g}")
    elif a.mode == "between":
        res = between(by_seq, a.flank_a, a.flank_b)
        if isinstance(res, str):
            print(res)
            return 1
        seqid, inner = res
        print(f"# between {a.flank_a} and {a.flank_b} on {seqid}: {len(inner)} gene(s)")
        if not inner:
            print("  <EMPTY SLOT>  -> candidate TRUE ABSENCE (verify flankers are the real orthologs)")
        for g in inner:
            print(f"  {g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
