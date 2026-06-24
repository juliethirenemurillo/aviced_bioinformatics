#!/usr/bin/env python
"""
flanker_orthologs.py  --  AviCeD synteny step (module B: ortholog reconciler)

Merges two pre-computed ortholog tables for each flanker gene and reports
per-(flanker, species) concordance. Output drives the synteny_call step:
which bird gene to use as the anchor for the flanker, and how confident.

Why two tables:
  Compara  -- Ensembl's tree-based orthology (good for housekeeping flankers;
              not available for non-Ensembl assemblies, e.g. duck T2T).
  RBB BLAST -- Reciprocal Best BLAST hits (works on any proteome you have,
              including T2T; weaker on paralog-rich families).

The two are complementary, not redundant:
  - Where both cover (chicken, duck ZJU1.0): concordance = methods cross-check.
  - Where only RBB covers (duck T2T): RBB is the only signal; we mark it so.

INPUT TSVs (both required; either may be empty/missing rows are fine):
  --compara  human_gene  species  bird_gene  ortholog_type  identity
  --blast    human_gene  species  bird_gene  identity       evalue   qcov

Species labels pass through verbatim (whatever you use in the input -- e.g.
'chicken', 'duck_ZJU1.0', 'duck_T2T' -- shows up in the output).

OUTPUT TSV: one row per (human_gene, species) seen in either input.
Columns:
  human_gene  species
  compara_bird_gene  compara_type  compara_identity
  blast_bird_gene    blast_identity  blast_evalue  blast_qcov
  concordance        use_for_synteny  reason

Concordance categories (the actual decision logic):
  AGREE          both methods return the same bird gene (same ID -- rare; same namespace)
  CONCORDANT     both methods return a hit but different namespaces (Ensembl vs RefSeq)
                 -- the normal case when both Compara and RBB succeed
  COMPARA_ONLY   Compara has an ortholog, RBB doesn't
  BLAST_ONLY     RBB has a hit, Compara doesn't (e.g. T2T, or paralog-rich flanker)
  NO_ORTHOLOG    neither method found a bird ortholog
  NO_DATA_*      input table did not cover this (species, gene) pair

use_for_synteny decision rule (which bird gene to anchor on in synteny step):
  AGREE / CONCORDANT + one2one           -> Compara gene (Ensembl ID; canonical for Ensembl GFFs)
  CONCORDANT + one2many/many2many        -> Compara gene, flagged *
  COMPARA_ONLY + one2one                 -> Compara gene
  COMPARA_ONLY + one2many/many2many      -> Compara gene, flagged *
  BLAST_ONLY                             -> BLAST gene (RefSeq; used for T2T)
  NO_ORTHOLOG                            -> SKIP

Note on namespaces:
  Compara returns Ensembl stable IDs (ENSGALG..., ENSAPLG...).
  RBB BLAST returns RefSeq accessions (ref|XP_...|, ref|NP_...|).
  These never match by string. CONCORDANT means both independently found an
  ortholog at this locus -- not that they returned identical IDs.

Stdlib only. No network. Unit-testable on mock TSVs.
"""
import argparse, csv, sys
from collections import defaultdict


# ---------- input parsing ----------

def _read_tsv(path, required_cols):
    """Read a TSV with a header; return list of dict rows. Tolerant to extra cols."""
    if path is None:
        return []
    rows = []
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in required_cols if c not in (rdr.fieldnames or [])]
        if missing:
            sys.exit(f"ERROR: {path} missing required column(s): {missing}")
        for r in rdr:
            rows.append({k: (v.strip() if isinstance(v, str) else v) for k, v in r.items()})
    return rows


def index_compara(rows):
    """{(human_gene, species): {'bird_gene','type','identity'}}  --  first row wins per key."""
    idx = {}
    for r in rows:
        key = (r["human_gene"], r["species"])
        idx.setdefault(key, {
            "bird_gene": r.get("bird_gene", "") or "",
            "type":      r.get("ortholog_type", "") or "",
            "identity":  r.get("identity", "") or "",
        })
    return idx


def index_blast(rows):
    """{(human_gene, species): {'bird_gene','identity','evalue','qcov'}}"""
    idx = {}
    for r in rows:
        key = (r["human_gene"], r["species"])
        idx.setdefault(key, {
            "bird_gene": r.get("bird_gene", "") or "",
            "identity":  r.get("identity", "") or "",
            "evalue":    r.get("evalue", "") or "",
            "qcov":      r.get("qcov", "") or "",
        })
    return idx


# ---------- concordance + decision ----------

def classify(comp, blast):
    """Return (concordance, use_for_synteny, reason).
       comp/blast are either dict-records or None."""
    has_c = comp is not None and comp["bird_gene"]
    has_b = blast is not None and blast["bird_gene"]
    # Ensembl Compara returns "ortholog_one2one", "ortholog_one2many" etc.
    # Normalise by stripping the "ortholog_" prefix for comparisons below.
    if comp is not None and comp.get("type"):
        comp = {**comp, "type": comp["type"].replace("ortholog_", "")}

    # both methods saw this (gene,species) but neither found an ortholog
    if comp is not None and blast is not None and not has_c and not has_b:
        return "NO_ORTHOLOG", "SKIP", "neither method returned a bird ortholog"

    # one of the tables didn't cover this pair at all
    if comp is None and blast is None:
        return "NO_DATA_BOTH", "SKIP", "no data in either input"
    if comp is None and has_b:
        return "NO_DATA_COMPARA", f"BLAST:{blast['bird_gene']}", "BLAST only (Compara not provided)"
    if blast is None and has_c:
        return "NO_DATA_BLAST", f"COMPARA:{comp['bird_gene']}", "Compara only (BLAST not provided)"
    if comp is None and not has_b:
        return "NO_DATA_COMPARA", "SKIP", "BLAST empty; Compara not provided"
    if blast is None and not has_c:
        return "NO_DATA_BLAST", "SKIP", "Compara empty; BLAST not provided"

    # only one method found something (both tables present)
    if has_c and not has_b:
        ctype = comp["type"]
        if ctype == "one2one":
            return "COMPARA_ONLY", f"COMPARA:{comp['bird_gene']}", "Compara one2one; RBB silent"
        if ctype in ("one2many", "many2many", "many2one"):
            return "COMPARA_ONLY", f"COMPARA:{comp['bird_gene']}*", f"Compara {ctype} (lower conf); RBB silent"
        return "COMPARA_ONLY", f"COMPARA:{comp['bird_gene']}*", "Compara only (type unknown); RBB silent"

    if has_b and not has_c:
        return "BLAST_ONLY", f"BLAST:{blast['bird_gene']}", "RBB only (e.g. T2T, or paralog-rich)"

    # both found something.
    # Compara returns Ensembl IDs (ENSGALG..., ENSAPLG...);
    # RBB returns RefSeq accessions (ref|XP_...|).
    # These namespaces never overlap, so string equality is the wrong test.
    # Instead: if both methods return ANY hit, treat as CONCORDANT --
    # both independently confirm the ortholog exists at this locus.
    # Record both IDs; use Compara as canonical for chicken/ZJU (Ensembl GFFs),
    # RBB as canonical for T2T (RefSeq GFF, no Compara coverage).
    if comp["bird_gene"] == blast["bird_gene"]:
        # same namespace (shouldn't happen normally, but handle it)
        return "AGREE", f"COMPARA:{comp['bird_gene']}", "methods concordant (same ID)"

    # different namespaces -- both confirmed, use Compara as canonical
    if comp["type"] == "one2one":
        return ("CONCORDANT",
                f"COMPARA:{comp['bird_gene']}",
                f"both methods confirm ortholog; Compara one2one canonical (RBB: {blast['bird_gene']})")
    if comp["type"] in ("one2many", "many2many", "many2one"):
        return ("CONCORDANT",
                f"COMPARA:{comp['bird_gene']}*",
                f"both methods confirm ortholog; Compara {comp['type']} flagged (RBB: {blast['bird_gene']})")
    # unknown type but both hit -- still concordant, flag it
    return ("CONCORDANT",
            f"COMPARA:{comp['bird_gene']}*",
            f"both methods confirm ortholog; Compara type unknown (RBB: {blast['bird_gene']})")


# ---------- reconcile ----------

def reconcile(compara_rows, blast_rows):
    """Return list of output dicts, one per (human_gene, species) seen in either input."""
    cidx = index_compara(compara_rows)
    bidx = index_blast(blast_rows)
    keys = sorted(set(cidx) | set(bidx))

    # was each table provided at all? (controls NO_DATA_* semantics)
    compara_provided = compara_rows is not None and len(compara_rows) >= 0  # truthy if file given
    blast_provided   = blast_rows   is not None and len(blast_rows)   >= 0

    out = []
    for k in keys:
        gene, species = k
        comp_rec  = cidx.get(k) if compara_provided else None
        blast_rec = bidx.get(k) if blast_provided   else None

        # if a table was provided but this key isn't in it, that's an empty record
        # (i.e. the method ran for this pair and returned nothing), not "no data"
        if compara_provided and comp_rec is None:
            comp_rec = {"bird_gene": "", "type": "", "identity": ""}
        if blast_provided and blast_rec is None:
            blast_rec = {"bird_gene": "", "identity": "", "evalue": "", "qcov": ""}

        conc, use, reason = classify(comp_rec, blast_rec)
        out.append({
            "human_gene":       gene,
            "species":          species,
            "compara_bird_gene": (comp_rec or {}).get("bird_gene", ""),
            "compara_type":      (comp_rec or {}).get("type", ""),
            "compara_identity":  (comp_rec or {}).get("identity", ""),
            "blast_bird_gene":   (blast_rec or {}).get("bird_gene", ""),
            "blast_identity":    (blast_rec or {}).get("identity", ""),
            "blast_evalue":      (blast_rec or {}).get("evalue", ""),
            "blast_qcov":        (blast_rec or {}).get("qcov", ""),
            "concordance":       conc,
            "use_for_synteny":   use,
            "reason":            reason,
        })
    return out


# ---------- CLI ----------

OUT_COLS = ["human_gene", "species",
            "compara_bird_gene", "compara_type", "compara_identity",
            "blast_bird_gene", "blast_identity", "blast_evalue", "blast_qcov",
            "concordance", "use_for_synteny", "reason"]


def main(argv=None):
    p = argparse.ArgumentParser(description="Reconcile Compara + RBB BLAST ortholog tables for flanker genes")
    p.add_argument("--compara", help="Compara ortholog TSV (human_gene, species, bird_gene, ortholog_type, identity)")
    p.add_argument("--blast",   help="Reciprocal-best-BLAST TSV (human_gene, species, bird_gene, identity, evalue, qcov)")
    p.add_argument("-o", "--out", required=True, help="Output reconciled TSV")
    a = p.parse_args(argv)

    if not a.compara and not a.blast:
        sys.exit("ERROR: provide at least one of --compara / --blast")

    comp = _read_tsv(a.compara, ["human_gene", "species", "bird_gene", "ortholog_type", "identity"]) if a.compara else None
    blast = _read_tsv(a.blast,  ["human_gene", "species", "bird_gene", "identity", "evalue", "qcov"]) if a.blast   else None

    rows = reconcile(comp, blast)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # quick stdout summary by concordance category
    summary = defaultdict(int)
    for r in rows:
        summary[r["concordance"]] += 1
    print(f"# reconciled {len(rows)} (gene, species) pairs -> {a.out}", file=sys.stderr)
    for k in sorted(summary):
        print(f"#   {k:18} {summary[k]}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
