#!/usr/bin/env python
"""
synteny_call.py  --  AviCeD synteny step (module C: tier + slot)

For each target gene, decides:
  (1) Is the human flanker neighbourhood recovered in the bird genome,
      and at which CONFIDENCE TIER?
  (2) Inside the recovered block, is there a gene in the expected slot
      (PRESENT-DIVERGED), is the slot empty (TRUE BIOLOGICAL ABSENCE),
      or is the locus broken / not recoverable (UNRESOLVED)?

Confidence tiers (from agreed framework):
  HIGH     window +/-5  (10 flankers)  threshold >=3-4 conserved
  MEDIUM   window +/-10 (20 flankers)  threshold >=3 conserved
  RESCUE   window +/-20 (40 flankers)  used only when assembly is poor; weak signal

The call assigns the highest tier whose threshold is met. If even RESCUE fails
to meet threshold, the locus is UNRESOLVED (annotation/assembly issue, not loss).

Inputs:
  --target          target gene symbol (e.g. NAIP)
  --human-gff       human GFF3 (defines the reference window)
  --bird-gff        bird GFF3
  --orthologs       module-B output TSV (per-flanker concordance)
  --species         bird species label (must match column in orthologs TSV)
  --protein-map     TSV from build_protein_gene_map.py (protein_accession -> GFF Name)
                    Resolves LOC-annotated genes that wouldn't be found by symbol alone.
  --high-window     default 5
  --high-thresh     default 3
  --med-window      default 10
  --med-thresh      default 3
  --rescue-window   default 20

Output: a single TSV row per (target, species) with the call.

Stdlib only. Imports gff_neighbourhood for GFF parsing.
"""
import argparse, csv, os, sys

# allow running from the synteny dir
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gff_neighbourhood import load_genes, order_by_position, _find  # noqa: E402


# ---------- protein map loading ----------

def load_protein_map(path):
    """Load protein_accession -> gene_name map from build_protein_gene_map.py output.
    Returns {protein_accession: gene_name}."""
    if path is None:
        return {}
    m = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m[r["protein_accession"]] = r["gene_name"]
    return m


def strip_ref(acc):
    """Strip 'ref|...|' wrapper from RefSeq accessions.
    'ref|XP_015128553.2|' -> 'XP_015128553.2'
    'XP_015128553.2' -> 'XP_015128553.2' (no-op)"""
    if acc.startswith("ref|") and acc.endswith("|"):
        return acc[4:-1]
    return acc


# ---------- ortholog loading ----------

def load_orthologs(path, species, protein_map=None):
    """Read module-B TSV; return {human_gene: bird_gff_lookup_key} for one species.

    Lookup key resolution priority:
      1. If blast_bird_gene is available AND protein_map has it -> use the GFF Name
         from the map (resolves LOC-annotated genes).
      2. Fall back to the human gene symbol (works when bird GFF uses the same symbol).
      3. SKIP rows -> None.
    """
    if protein_map is None:
        protein_map = {}
    m = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["species"] != species:
                continue
            use = r.get("use_for_synteny", "") or ""
            if use == "SKIP" or use == "":
                m[r["human_gene"]] = None
                continue
            # try protein map first (resolves LOC genes)
            blast_acc = r.get("blast_bird_gene", "") or ""
            if blast_acc:
                clean_acc = strip_ref(blast_acc)
                if clean_acc in protein_map:
                    m[r["human_gene"]] = protein_map[clean_acc]
                    continue
            # fallback: use human gene symbol directly
            m[r["human_gene"]] = r["human_gene"]
    return m


# ---------- window extraction ----------

def human_window(by_seq_human, target, k):
    """Return the ordered list of gene names in a +/-k window around target
    in the human GFF, plus the index of the target inside that window."""
    hit = _find(by_seq_human, target)
    if hit is None:
        return None, None
    seqid, i = hit
    glist = by_seq_human[seqid]
    lo, hi = max(0, i - k), min(len(glist), i + k + 1)
    window = [g.name for g in glist[lo:hi]]
    target_idx = i - lo
    return window, target_idx


# ---------- cluster finding ----------

def find_bird_cluster(window_birds, by_seq_bird, max_span_genes=80):
    """Locate the largest cluster of `window_birds` in the bird genome.
       Returns (seqid, (lo_idx, hi_idx), conserved_count, names) or None."""
    positions = []
    for bird in window_birds:
        if not bird:
            continue
        hit = _find(by_seq_bird, bird)
        if hit is not None:
            positions.append(hit + (bird,))

    if not positions:
        return None

    by_chr = {}
    for seqid, idx, name in positions:
        by_chr.setdefault(seqid, []).append((idx, name))

    best = None
    for seqid, hits in by_chr.items():
        hits.sort()
        for i in range(len(hits)):
            j = i
            while j + 1 < len(hits) and hits[j + 1][0] - hits[i][0] <= max_span_genes:
                j += 1
            names_here = [h[1] for h in hits[i:j + 1]]
            unique = len(set(names_here))
            if best is None or unique > best[0]:
                best = (unique, seqid, hits[i][0], hits[j][0], names_here)
    if best is None:
        return None
    unique, seqid, lo, hi, names = best
    return seqid, (lo, hi), unique, names


# ---------- slot test ----------

def slot_status(by_seq_bird, seqid, lo_idx, hi_idx, target_window_idx,
                window_len, bird_flanker_positions):
    """Inside the recovered bird block, classify the slot for the target.

    Splits recovered flankers into 'left of target' and 'right of target'
    based on their position in the human window. The slot is between the
    innermost left-flanker and innermost right-flanker in bird genomic order.

    Handles inversions: if flanker order is reversed relative to human
    (left-flankers appear after right-flankers in bird coords), the slot
    boundaries are swapped. This is a local inversion -- common in immune
    gene neighbourhoods -- and doesn't invalidate the synteny call.
    """
    glist = by_seq_bird[seqid]
    left  = [bi for wi, bi in bird_flanker_positions if wi < target_window_idx]
    right = [bi for wi, bi in bird_flanker_positions if wi > target_window_idx]

    if not left or not right:
        # report which side is missing
        if not left and not right:
            return "UNRESOLVED_SLOT", [], "no flankers recovered in cluster"
        missing = "left" if not left else "right"
        return "UNRESOLVED_SLOT", [], f"missing flankers on {missing} side; slot undefined"

    # innermost flankers: rightmost of left-set, leftmost of right-set
    inner_left  = max(left)
    inner_right = min(right)

    # handle inversion: if bird order is reversed, swap boundaries
    inverted = False
    if inner_left >= inner_right:
        inner_left, inner_right = inner_right, inner_left
        inverted = True

    inner = glist[inner_left + 1:inner_right]
    inv_note = " (local inversion detected)" if inverted else ""

    if not inner:
        return ("EMPTY_SLOT", [],
                f"no annotated gene between recovered flankers{inv_note} (candidate true absence)")
    return ("PRESENT_DIVERGED_CANDIDATE",
            [g.name for g in inner],
            f"{len(inner)} gene(s) in slot{inv_note}")


# ---------- per-tier run ----------

def run_one_tier(target, human_by_seq, bird_by_seq, ortho_map, k):
    """Compute window, map to bird, find cluster, slot test. Returns dict."""
    window_names, target_idx = human_window(human_by_seq, target, k)
    if window_names is None:
        return {"error": f"target {target} not found in human GFF"}

    bird_names = [ortho_map.get(g) for g in window_names]
    flanker_pairs = [(i, b) for i, (g, b) in enumerate(zip(window_names, bird_names))
                     if i != target_idx and b]

    cluster = find_bird_cluster([b for _, b in flanker_pairs], bird_by_seq)
    if cluster is None:
        return {
            "window": 2 * k + 1, "conserved": 0,
            "locus_seqid": "", "locus_span": "",
            "slot_status": "UNRESOLVED_LOCUS",
            "slot_inner": "",
            "reason": "no flanker orthologs located in bird GFF",
            "window_names": ";".join(window_names),
        }

    seqid, (lo, hi), unique, recovered_names = cluster
    flanker_positions = []
    for wi, b in flanker_pairs:
        hit = _find(bird_by_seq, b)
        if hit is not None and hit[0] == seqid and lo <= hit[1] <= hi:
            flanker_positions.append((wi, hit[1]))

    slot, inner_names, slot_reason = slot_status(
        bird_by_seq, seqid, lo, hi, target_idx, 2 * k + 1, flanker_positions)
    return {
        "window": 2 * k + 1, "conserved": unique,
        "locus_seqid": seqid, "locus_span": f"{lo}-{hi}",
        "slot_status": slot,
        "slot_inner": ";".join(inner_names) if inner_names else "",
        "reason": slot_reason,
        "window_names": ";".join(window_names),
    }


# ---------- main ----------

def main(argv=None):
    p = argparse.ArgumentParser(description="Synteny call: tier + slot test for one target")
    p.add_argument("--target", required=True)
    p.add_argument("--human-gff", required=True)
    p.add_argument("--bird-gff", required=True)
    p.add_argument("--orthologs", required=True)
    p.add_argument("--species", required=True)
    p.add_argument("--protein-map", default=None,
                   help="TSV from build_protein_gene_map.py (resolves LOC genes)")
    p.add_argument("--high-window", type=int, default=5)
    p.add_argument("--high-thresh", type=int, default=3)
    p.add_argument("--med-window", type=int, default=10)
    p.add_argument("--med-thresh", type=int, default=3)
    p.add_argument("--rescue-window", type=int, default=20)
    p.add_argument("--rescue-thresh", type=int, default=3)
    p.add_argument("--out-tsv",
                   help="append result as TSV row (header written if file doesn't exist)")
    a = p.parse_args(argv)

    human_by_seq = order_by_position(load_genes(a.human_gff))
    bird_by_seq  = order_by_position(load_genes(a.bird_gff))
    prot_map     = load_protein_map(a.protein_map)
    ortho_map    = load_orthologs(a.orthologs, a.species, protein_map=prot_map)

    tiers = [
        ("HIGH",   a.high_window,   a.high_thresh),
        ("MEDIUM", a.med_window,    a.med_thresh),
        ("RESCUE", a.rescue_window, a.rescue_thresh),
    ]

    chosen = None
    per_tier = {}
    for name, w, t in tiers:
        res = run_one_tier(a.target, human_by_seq, bird_by_seq, ortho_map, w)
        per_tier[name] = res
        if "error" in res:
            chosen = (name, res)
            break
        if res["conserved"] >= t and chosen is None:
            chosen = (name, res)

    if chosen is None:
        chosen = ("UNRESOLVED_TIER", per_tier["RESCUE"])

    tier_name, res = chosen
    row = {
        "target":        a.target,
        "species":       a.species,
        "tier":          tier_name,
        "window":        res.get("window", ""),
        "conserved":     res.get("conserved", 0),
        "locus_seqid":   res.get("locus_seqid", ""),
        "locus_span":    res.get("locus_span", ""),
        "slot_status":   res.get("slot_status", ""),
        "slot_inner":    res.get("slot_inner", ""),
        "reason":        res.get("reason", res.get("error", "")),
        "window_genes":  res.get("window_names", ""),
    }

    cols = ["target", "species", "tier", "window", "conserved",
            "locus_seqid", "locus_span", "slot_status", "slot_inner",
            "reason", "window_genes"]

    if a.out_tsv:
        write_header = not os.path.exists(a.out_tsv)
        with open(a.out_tsv, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            if write_header:
                w.writeheader()
            w.writerow(row)
    w = csv.DictWriter(sys.stdout, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerow(row)
    return 0


if __name__ == "__main__":
    sys.exit(main())
