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
  --high-window     default 5
  --high-thresh     default 3
  --med-window      default 10
  --med-thresh      default 3
  --rescue-window   default 20

Output: a single TSV row per (target, species) with the call. Suitable for
appending to a synteny_results.tsv across all 8 targets x 3 species.

Stdlib only. Imports gff_neighbourhood for GFF parsing.
"""
import argparse, csv, os, sys
from collections import Counter

# allow running from the synteny dir
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gff_neighbourhood import load_genes, order_by_position, _find  # noqa: E402


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


def load_orthologs(path, species):
    """Read module-B TSV; return {human_gene: gff_lookup_key} for one species.

    The bird GFFs (RefSeq-annotated) use gene symbols as Name= values, not
    Ensembl IDs or RefSeq accessions. Since gene symbols are conserved between
    human and bird for well-annotated genes, we use the human_gene symbol
    directly as the GFF lookup key. The Ensembl/RefSeq IDs in use_for_synteny
    confirm the ortholog exists but are not used for GFF lookup.

    Returns {human_gene: human_gene} for usable rows, {human_gene: None} for SKIP.
    """
    m = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["species"] != species:
                continue
            use = r.get("use_for_synteny", "") or ""
            if use == "SKIP" or use == "":
                m[r["human_gene"]] = None
            else:
                # use the human gene symbol as the GFF Name= lookup key
                m[r["human_gene"]] = r["human_gene"]
    return m


def find_bird_cluster(window_birds, by_seq_bird, max_span_genes=80):
    """Locate the largest cluster of `window_birds` in the bird genome.
       Returns (seqid, indices_of_cluster_in_bird_order, conserved_count) or None.

       Strategy:
         1. Find every position in the bird genome where any window flanker sits.
         2. Group by seqid.
         3. On each seqid, find the densest window of size ~max_span_genes that
            contains the most flankers. Report it.
       max_span_genes guards against picking up scattered random hits across the
       whole chromosome; the syntenic block should be tight."""
    # all positions (seqid, idx) of every bird flanker
    positions = []
    for bird in window_birds:
        if not bird:
            continue
        hit = _find(by_seq_bird, bird)
        if hit is not None:
            positions.append(hit + (bird,))  # (seqid, idx, name)

    if not positions:
        return None

    # group by seqid
    by_chr = {}
    for seqid, idx, name in positions:
        by_chr.setdefault(seqid, []).append((idx, name))

    # for each chromosome, find tightest window with most flankers
    best = None  # (count, seqid, lo_idx, hi_idx, names)
    for seqid, hits in by_chr.items():
        hits.sort()
        # sliding window over hit-indices
        for i in range(len(hits)):
            j = i
            while j + 1 < len(hits) and hits[j + 1][0] - hits[i][0] <= max_span_genes:
                j += 1
            cnt = j - i + 1
            names_here = [h[1] for h in hits[i:j + 1]]
            unique = len(set(names_here))
            if best is None or unique > best[0]:
                best = (unique, seqid, hits[i][0], hits[j][0], names_here)
    if best is None:
        return None
    unique, seqid, lo, hi, names = best
    return seqid, (lo, hi), unique, names


def slot_status(by_seq_bird, seqid, lo_idx, hi_idx, target_window_idx, window_len, bird_flanker_positions):
    """Inside the recovered bird block, classify the slot for the target.

    Logic:
      - The human window has the target at index `target_window_idx`.
      - We split the recovered bird flankers into 'left of target' and 'right of
        target' based on their original window indices.
      - The expected bird slot is the region between the rightmost left-flanker
        and the leftmost right-flanker in the bird genome.
      - Genes in that slot region: candidate(s) for diverged ortholog.
      - No genes in that slot region: EMPTY SLOT.
      - One side has no recovered flanker: cannot define slot -> UNRESOLVED.
    """
    glist = by_seq_bird[seqid]
    # bird_flanker_positions: list of (window_idx_in_human, bird_idx) for recovered flankers
    left  = [bi for wi, bi in bird_flanker_positions if wi < target_window_idx]
    right = [bi for wi, bi in bird_flanker_positions if wi > target_window_idx]
    if not left or not right:
        return "UNRESOLVED_SLOT", [], "missing flankers on one side; slot undefined"
    # define slot in bird coords
    slot_lo = max(left)
    slot_hi = min(right)
    if slot_lo >= slot_hi:
        # rearranged: left flanker is to the right of right flanker
        return "UNRESOLVED_SLOT", [], "flankers rearranged (left > right in bird order)"
    inner = glist[slot_lo + 1:slot_hi]
    if not inner:
        return "EMPTY_SLOT", [], "no annotated gene between recovered flankers (candidate true absence)"
    return "PRESENT_DIVERGED_CANDIDATE", [g.name for g in inner], f"{len(inner)} gene(s) in slot"


def assign_tier(conserved, window_total, tiers):
    """Given conserved count and total window size, return highest tier that
    meets its threshold. tiers = [(name, window, threshold), ...]."""
    for name, w, t in tiers:
        if window_total == 2 * w + 1 and conserved >= t:
            return name
    return "UNRESOLVED_TIER"


def run_one_tier(target, human_by_seq, bird_by_seq, ortho_map, k):
    """Compute window, map to bird, find cluster, slot test. Returns dict or None."""
    window_names, target_idx = human_window(human_by_seq, target, k)
    if window_names is None:
        return {"error": f"target {target} not found in human GFF"}
    # build bird_names list aligned to window_names
    bird_names = [ortho_map.get(g) for g in window_names]
    # exclude the target itself from the flanker set
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
    # build (window_idx, bird_idx) list for the slot test
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


def main(argv=None):
    p = argparse.ArgumentParser(description="Synteny call: tier + slot test for one target")
    p.add_argument("--target", required=True)
    p.add_argument("--human-gff", required=True)
    p.add_argument("--bird-gff", required=True)
    p.add_argument("--orthologs", required=True)
    p.add_argument("--species", required=True, help="must match species label in orthologs TSV")
    p.add_argument("--high-window", type=int, default=5)
    p.add_argument("--high-thresh", type=int, default=3)
    p.add_argument("--med-window", type=int, default=10)
    p.add_argument("--med-thresh", type=int, default=3)
    p.add_argument("--rescue-window", type=int, default=20)
    p.add_argument("--rescue-thresh", type=int, default=3)
    p.add_argument("--out-tsv", help="if given, append result as TSV row (with header if file doesn't exist)")
    a = p.parse_args(argv)

    human_by_seq = order_by_position(load_genes(a.human_gff))
    bird_by_seq  = order_by_position(load_genes(a.bird_gff))
    ortho_map    = load_orthologs(a.orthologs, a.species)

    tiers = [
        ("HIGH",   a.high_window,   a.high_thresh),
        ("MEDIUM", a.med_window,    a.med_thresh),
        ("RESCUE", a.rescue_window, a.rescue_thresh),
    ]

    # try tiers in order of strictness; accept highest that meets threshold
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
            # don't break: still record higher-window tiers for transparency

    if chosen is None:
        # no tier met threshold -- report rescue-tier result as UNRESOLVED
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
    # also stdout
    w = csv.DictWriter(sys.stdout, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerow(row)
    return 0


if __name__ == "__main__":
    sys.exit(main())
