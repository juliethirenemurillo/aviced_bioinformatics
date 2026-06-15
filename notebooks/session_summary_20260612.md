# AviCeD Bioinformatics — Session Summary
## Date: 12 June 2026 (resumed 15 June)
## Purpose: Pick up seamlessly in a new chat

---

## ONE-LINE RESUME
Re-run the retry-enabled FASTA fetch for the 8 Group-1 queries, confirm 8 headers, then download BLAST targets (RefSeq protein + newest genome assemblies) and submit the blastp+tblastn disambiguation job to `BRYANT-SL3-CPU`.

---

## WHAT WE DID THIS SESSION

### 1. CSD3 infrastructure — all confirmed working
- **Compute:** `BRYANT-SL3-CPU`, ~200k core-hours free (used <1%). Submit with `#SBATCH -A BRYANT-SL3-CPU`. Service level SL3 (448 cores / 12 h max — far more than needed).
- **Storage:** bought RDS exists → `~/comparativegenom` (symlink to `/rds/project/rds-ofC3KNq1BSo`), 1 TB live tier, empty. NOT backed up — GitHub is the only real code backup. Free `~/rds/hpc-work` (rds-d7) holds ATACseq, ~474 GB / 1 TB used.
- **Working tree built:** `~/comparativegenom/aviced/blast_disambiguation/{genomes,proteins,queries,results,logs,scripts}`.
- **BLAST:** `module load blast/2.17` → blastp, tblastn, makeblastdb. NCBI `datasets` NOT installed (fallback needed). Modules don't persist across sessions → the `module load` line must also go INSIDE the Slurm script.
- **Login vs compute nodes:** login has internet (do downloads here), compute has none (BLAST runs here against pre-downloaded files).

### 2. Task — Option A: disambiguate "absent" calls (Group 1)
8 query genes, verified UniProt **protein** accessions (pulled from `human_cell_death_gene_universe_v1.xlsx`):
`BAX Q07812 · MAVS Q7Z434 · RBCK1 Q9BYM8 · BID P55957 · BAD Q92934 · BBC3 Q96PG8 · PMAIP1 Q13794 · MLKL Q8NB16` (MLKL duck-only).
- **Status:** FASTA fetch partially done — BAX/MAVS/RBCK1 fetched, last 5 failed on UniProt 503 (rate-limit). Fix = retry-enabled fetch with `--retry 3 --retry-delay 2` + `sleep 1`. Just needs re-running.
- **Plan:** blastp (vs RefSeq protein) + tblastn (vs genome) together → 3-way read per gene:
  - hit + hit = annotated protein exists, Compara missed it (divergence/ID issue)
  - tblastn hit + no blastp = sequence present, no gene model → **annotation gap**
  - no + no = **true biological absence** (modulo assembly/sensitivity)
- Newest RefSeq assemblies, accession + date documented. No cutoffs pre-applied (raw results, decide thresholds together). Synteny stays later. Results go to a SEPARATE file, never into the curated table without asking.

### 3. Methodology insights (the conceptual core)
- **Co-orthology / many-to-1:** two human paralogues → one bird gene is biology, not error. Resolve duplication-vs-loss with the gene tree; **outgroup-count shortcut** uses spotted gar (`lepisosteus_oculatus`) or *Xenopus* — NEVER zebrafish (teleost whole-genome duplication inflates copy counts).
- **Gasdermins resolved:** GSDME + PJVK deep-ancestral; GSDMA = amniote clade-founder (the bird inflammatory one, caspase-1-cleaved — YVAD-linker story); GSDMB/C/D = mammalian inventions → bird absence is **mammalian gain, not avian loss** → off the interesting list. Keep bird GSDME and GSDMA as distinct rows (different depth, different biology).
- **Caspases:** CASP4/5 = mammalian/primate expansion; bird has ancestral CASP1. **FLAG:** published chicken caspase-1 is ~45% identical to human (cloned 1998), contradicting the table's "DUBIOUS 10–20%" → likely mapping/metric artifact. Recheck: pull actual target Ensembl ID, recompute identity per-domain, read both id% AND pos%. Caspase-1's real story = substrate-site drift (variant active-site pentapeptide, IL-1β site), NOT detectability → it's a POOR poster-child for "structure beats sequence." Reassign that headline role to NLRP1/3, ASC.
- **id% vs pos% (correction):** both are pairwise alignment stats, NOT orthology scores. Orthology call = tree reconciliation (separate column). The id%/pos% gap signals conservative-substitution divergence; low on BOTH = probably aligned to wrong gene.
- **Engines:** Compara orthology = HMM family-clustering (TreeFam-style) + tree reconciliation. InterPro = HMM-based member methods (Pfam/SMART/PANTHER) or pre-computed ID lookup, depending on route taken in Step 0.

### 4. Cowork — parallel literature track
- 6-batch literature-table prompt drafted (gene-symbol-keyed; columns include DOI + "citation real? confidence" + "overlapping vs adjacent" novelty tags). Run as ONE block, desk-side; review Batch 1 (caspases) first. Companion dedupe prompt drafted for AFTER (cross-checks against Billman/Tanzer/Newton/Steinegger; reports which known anchors it failed to surface = blind-spot test).
- Cowork is desktop-only; phone (Dispatch) can monitor a running desktop session but isn't a standalone route — run literature work at the laptop.

---

## NEXT ACTIONS (in order)
1. Re-run retry-enabled FASTA fetch → confirm 8 headers naming the right proteins.
2. Download targets on login node: newest RefSeq protein sets + genome assemblies for chicken + duck (document accessions/date).
3. Write download script + Slurm script (commit/push to repo; genomes & results git-ignored).
4. Submit blastp+tblastn, charge `BRYANT-SL3-CPU`.
5. Review raw results together; THEN decide thresholds and what (if anything) updates the curated table.
6. Re-check the DUBIOUS caspase row as part of the BLAST pass (resolves which bird gene is the true CASP1 orthologue).

## STANDING RULES (carried)
- Exploratory work: label hypotheses as hypotheses, keep separate from data.
- Ask before adding any interpretive column to curated tables.
- Don't re-derive grant / three-paper strategy — settled.
- Explain HPC/Slurm/Snakemake from first principles; define jargon on first use.
- No squeue/slurmctld RPC calls in loops (use sacct or throttled polling).
