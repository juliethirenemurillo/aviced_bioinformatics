#!/usr/bin/env python
"""
flanker_orthologs.py  --  AviCeD synteny step (reconciler, v3)

Reconciles three sources of evidence for each (flanker, species) pair:
  1. Compara orthology (Ensembl tree-based)
  2. blastp top forward hit (with reciprocal-or-paralog flag)
  3. tblastn against bird genome (catches unannotated genes)

This v3 lets us accept paralogs as POSITIONAL flankers for synteny --
they're "the gene at the right place even if it's the wrong sister-gene
family member" -- which is what we actually need for the synteny step.
For target genes (the 8 we're investigating), strict reciprocity matters;
for flankers, position matters.

Input schemas:
  compara.tsv: human_gene, species, bird_gene, ortholog_type, identity
  blast_v3.tsv: human_gene, species,
                forward_top_hit, forward_pident, forward_evalue, forward_qcov,
                reverse_symbol, reciprocal_match,
                tblastn_seqid, tblastn_pident, tblastn_evalue, tblastn_qcov

Output schema (richer than v2):
  human_gene, species,
  compara_bird_gene, compara_type, compara_identity,
  blast_forward_top_hit, blast_forward_pident, blast_forward_qcov,
  blast_reverse_symbol, blast_reciprocal,
  tblastn_seqid, tblastn_pident,
  concordance, use_for_synteny, reason

Concordance categories:
  AGREE              compara + blastp_reciprocal both return same gene
  CONCORDANT         compara + blastp_reciprocal both return (different namespaces)
  COMPARA_ONLY       compara has ortholog, no blastp reciprocal
  RECIPROCAL_ONLY    blastp reciprocal, no compara
  PARALOG_FLANKER    blastp returns paralog (not reciprocal); usable as positional flanker
  GENOME_HIT         tblastn found a hit in genome (unannotated gene)
  NO_ORTHOLOG        no evidence anywhere

use_for_synteny decision rule:
  AGREE / CONCORDANT one2one / many2one     -> Compara bird_gene
  AGREE / CONCORDANT one2many / many2many   -> Compara bird_gene, flagged *
  COMPARA_ONLY one2one / many2one           -> Compara bird_gene
  COMPARA_ONLY one2many / many2many         -> Compara bird_gene, flagged *
  RECIPROCAL_ONLY                           -> blast_forward_top_hit
  PARALOG_FLANKER                           -> blast_forward_top_hit, flagged ~paralog
  GENOME_HIT                                -> SKIP (synteny module can't use coord-only yet)
  NO_ORTHOLOG                               -> SKIP

Reciprocity is checked by SwissProt GN= symbol match, requiring exact gene
symbol equality. Paralog hits are detected by reciprocal_match=NO with a
non-empty reverse_symbol.
"""

import argparse, csv, sys


def normalise_compara_type(t):
    """Strip 'ortholog_' prefix Ensembl REST uses (e.g. 'ortholog_one2one' -> 'one2one')."""
    return (t or '').replace('ortholog_', '')


COMPARA_HIGH_CONF = {'one2one', 'many2one'}
COMPARA_LOW_CONF = {'one2many', 'many2many'}


def load_compara(path):
    """{(human_gene, species): row}"""
    if not path:
        return {}
    out = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for r in reader:
            key = (r['human_gene'], r['species'])
            out[key] = {
                'bird_gene': r['bird_gene'],
                'ortholog_type': normalise_compara_type(r['ortholog_type']),
                'identity': r.get('identity', '0'),
            }
    return out


def load_blast(path):
    """{(human_gene, species): row}  -- new v3 schema"""
    if not path:
        return {}
    out = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for r in reader:
            key = (r['human_gene'], r['species'])
            out[key] = {
                'forward_top_hit':  r.get('forward_top_hit', ''),
                'forward_pident':   r.get('forward_pident', ''),
                'forward_qcov':     r.get('forward_qcov', ''),
                'reverse_symbol':   r.get('reverse_symbol', ''),
                'reciprocal_match': r.get('reciprocal_match', ''),
                'tblastn_seqid':    r.get('tblastn_seqid', ''),
                'tblastn_pident':   r.get('tblastn_pident', ''),
            }
    return out


def classify(compara_rec, blast_rec, human_gene):
    """Return (concordance, use_for_synteny, reason)."""
    has_compara = compara_rec and compara_rec['bird_gene']
    has_blast_reciprocal = (
        blast_rec
        and blast_rec.get('forward_top_hit')
        and blast_rec.get('reciprocal_match') == 'YES'
    )
    has_blast_paralog = (
        blast_rec
        and blast_rec.get('forward_top_hit')
        and blast_rec.get('reciprocal_match') == 'NO'
    )
    has_tblastn = blast_rec and blast_rec.get('tblastn_seqid')

    if has_compara and has_blast_reciprocal:
        ctype = compara_rec['ortholog_type']
        compara_gene = compara_rec['bird_gene']
        blast_gene = blast_rec['forward_top_hit']
        if compara_gene == blast_gene:
            concordance = 'AGREE'
        else:
            concordance = 'CONCORDANT'
        if ctype in COMPARA_HIGH_CONF:
            use = f'COMPARA:{compara_gene}'
            reason = f'compara {ctype}, blastp reciprocal confirms (top hit {blast_gene})'
        elif ctype in COMPARA_LOW_CONF:
            use = f'COMPARA:{compara_gene}*'
            reason = f'compara {ctype} (low conf), blastp reciprocal confirms'
        else:
            use = f'COMPARA:{compara_gene}*'
            reason = f'compara type unknown ({ctype}), blastp reciprocal confirms'
        return concordance, use, reason

    if has_compara and not has_blast_reciprocal:
        ctype = compara_rec['ortholog_type']
        compara_gene = compara_rec['bird_gene']
        if ctype in COMPARA_HIGH_CONF:
            use = f'COMPARA:{compara_gene}'
            reason = f'compara {ctype}, blastp silent or paralog'
        elif ctype in COMPARA_LOW_CONF:
            use = f'COMPARA:{compara_gene}*'
            reason = f'compara {ctype} only'
        else:
            use = f'COMPARA:{compara_gene}*'
            reason = f'compara type unknown ({ctype})'
        return 'COMPARA_ONLY', use, reason

    if has_blast_reciprocal:
        bg = blast_rec['forward_top_hit']
        return ('RECIPROCAL_ONLY',
                f'BLAST:{bg}',
                f'blastp reciprocal best hit')

    if has_blast_paralog:
        bg = blast_rec['forward_top_hit']
        rev = blast_rec.get('reverse_symbol', '')
        return ('PARALOG_FLANKER',
                f'BLAST:{bg}~{rev}',
                f'top blastp hit is paralog (reverses to {rev}); positional use only')

    if has_tblastn:
        seqid = blast_rec['tblastn_seqid']
        pident = blast_rec.get('tblastn_pident', '')
        return ('GENOME_HIT',
                'SKIP',
                f'tblastn hit at {seqid} ({pident}%); unannotated, synteny lookup-by-coord not implemented')

    return 'NO_ORTHOLOG', 'SKIP', 'no evidence in compara, blastp, or tblastn'


def reconcile(compara_path, blast_path, out_path):
    compara = load_compara(compara_path)
    blast = load_blast(blast_path)

    all_keys = sorted(set(compara) | set(blast))

    out_cols = [
        'human_gene', 'species',
        'compara_bird_gene', 'compara_type', 'compara_identity',
        'blast_forward_top_hit', 'blast_forward_pident', 'blast_forward_qcov',
        'blast_reverse_symbol', 'blast_reciprocal',
        'tblastn_seqid', 'tblastn_pident',
        'concordance', 'use_for_synteny', 'reason',
    ]

    counts = {}
    with open(out_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=out_cols, delimiter='\t')
        writer.writeheader()
        for human_gene, species in all_keys:
            compara_rec = compara.get((human_gene, species))
            blast_rec = blast.get((human_gene, species))
            concordance, use, reason = classify(compara_rec, blast_rec, human_gene)
            row = {
                'human_gene': human_gene,
                'species': species,
                'compara_bird_gene': compara_rec['bird_gene'] if compara_rec else '',
                'compara_type': compara_rec['ortholog_type'] if compara_rec else '',
                'compara_identity': compara_rec['identity'] if compara_rec else '',
                'blast_forward_top_hit': blast_rec['forward_top_hit'] if blast_rec else '',
                'blast_forward_pident': blast_rec.get('forward_pident', '') if blast_rec else '',
                'blast_forward_qcov': blast_rec.get('forward_qcov', '') if blast_rec else '',
                'blast_reverse_symbol': blast_rec.get('reverse_symbol', '') if blast_rec else '',
                'blast_reciprocal': blast_rec.get('reciprocal_match', '') if blast_rec else '',
                'tblastn_seqid': blast_rec.get('tblastn_seqid', '') if blast_rec else '',
                'tblastn_pident': blast_rec.get('tblastn_pident', '') if blast_rec else '',
                'concordance': concordance,
                'use_for_synteny': use,
                'reason': reason,
            }
            writer.writerow(row)
            counts[concordance] = counts.get(concordance, 0) + 1

    print(f'# reconciled {len(all_keys)} (gene, species) pairs -> {out_path}', file=sys.stderr)
    for cat in sorted(counts):
        print(f'#   {cat:18} {counts[cat]}', file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description='Reconcile flanker ortholog evidence v3')
    p.add_argument('--compara', help='compara.tsv from build_flanker_inputs.sh')
    p.add_argument('--blast', help='blast_v3.tsv from build_flanker_inputs.sh')
    p.add_argument('-o', '--out', required=True, help='reconciled output TSV')
    a = p.parse_args()
    reconcile(a.compara, a.blast, a.out)


if __name__ == '__main__':
    main()
