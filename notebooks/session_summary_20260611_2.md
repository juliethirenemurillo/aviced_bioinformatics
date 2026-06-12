# AviCeD Bioinformatics — Session Summary
## Date: 11 June 2026 (updated end-of-day)
## Purpose: Pick up seamlessly in a new chat

---

## 1. PAPER STRATEGY — THREE PAPERS DEFINED AND ORDERED

**Paper 1 — Structural Genomics Atlas (PRIMARY, start now)**
- Research question: What cell death pathway components do birds actually have, and what structural evidence exists for novel or diverged effectors invisible to sequence homology?
- Scope: Pan-avian (~50 species), full cell death landscape (pyroptosis, necroptosis, apoptosis, PANoptosis, sensors, regulators)
- Two-arm Foldseek screen: Arm 1 (reference-anchored, recovers diverged orthologues) + Arm 2 (unanchored, discovers novel structural candidates with three-condition filter: high TM-score + no sequence homology + no synteny)
- Target: MBE or eLife

**Paper 3 — Protease-substrate co-evolution (can start in parallel)**
- Research question: How do caspase active-site specificity and substrate cleavage motifs co-evolve across the avian phylogeny, and does co-evolutionary breakdown explain the rewiring of pyroptotic signalling?
- General aim: Map the co-evolutionary dynamics of protease pocket specificity and substrate tetrapeptide motifs across birds.
- Specific aim 1: Active site pocket profiling — AlphaFold structures of caspase-1 across birds, S1–S4 subsite classification. Tools: AlphaFold, PyMOL/ChimeraX, fpocket or SiteMap.
- Specific aim 2: Substrate motif scanning — scan GSDMA, GSDME, IL-1β, IL-18 across birds for tetrapeptide cleavage motifs. Tools: custom scripts, MAFFT.
- Specific aim 3: Co-evolutionary analysis — paired dN/dS, mirror tree / DCA. Tools: HyPhy, BayesTraits, plmDCA or GREMLIN.
- Start broad; conclusions drawn from what the data shows. Known published context (Billman/Miao 2024): the GSDMD→GSDMA substrate switch and songbird caspase-1 loss are already described — these are backdrop, not the paper's contribution. The co-evolutionary coupling quantification is the candidate novel angle, but let the data decide.
- Uses only public data + AlphaFold structures — no dependency on Paper 1
- May fold into Paper 1 or stand alone depending on signal strength
- Target: MBE standalone or chapter in Paper 1

**Paper 2 — Knowledge Graph / Systems Biology (starts after Paper 1 atlas)**
- Research question: Can interaction rules learned from well-characterised mammalian cell death networks predict how avian cell death pathways are wired?
- Architecture: KG link prediction (RotatE/ComplEx) trained on human + mouse + zebrafish, applied to bird proteins from Paper 1. Node features: ESM2 embeddings + InterPro domain architecture + expression. Edge types: cleaves, activates, scaffolds, inhibits, phosphorylates.
- Validation: Boucher cleavage assays, Stella KO phenotypes, Digital Fowl simulations.
- Computationally cheap (ESM2 embedding ~30 min GPU; RotatE training ~10–30 min).
- Target: Mol Sys Bio or PLOS Comp Bio

**Order: Paper 1 first → Paper 3 in parallel → Paper 2 after atlas.**

---

## 2. COMPETITIVE LITERATURE — WHAT IS ALREADY PUBLISHED

Critical reading before writing any Paper 1/3 claims. Both papers below are cited in the AviCeD grant (refs 15–16):

- **Billman et al. 2024, eLife** ("Caspase-1 activates gasdermin A in non-mammals"): Birds/reptiles/amphibians lack GSDMD; caspase-1 cleaves GSDMA instead. The substrate switch GSDME (fish) → GSDMA (amphibians/reptiles/birds) → GSDMD (mammals) is fully described. Bird caspase-1 uses tetrapeptide specificity (YVAD/FASD); mammalian caspase-1 uses exosite interactions with GSDMD C-terminus.
- **Billman et al. 2024, GBE** ("Unanticipated Loss of Inflammasomes in Birds"): Songbirds lost caspase-1 but retain IL-1β, IL-18, and YVAD-linker GSDMA. ASC lost in multiple independent events across bird phylogeny.

**What this means for your papers:**
- Do NOT frame findings around the GSDMD→GSDMA switch or songbird caspase-1 loss — already published.
- DO use Foldseek structural recovery to find things sequence-based surveys missed (Billman is largely sequence-based).
- The open question Billman explicitly flags: how do caspase-1-deficient birds mature and secrete IL-1β and IL-18? Your D80 (IL-1β) vs LESD (IL-18) asymmetry is an entry point.
- RIG-I/DDX58 chicken-absent/duck-present is also published (Campbell & Magor 2020, your grant ref 17) — use as positive control, not a result.

---

## 3. PAPER 1 PIPELINE — SEVEN STEPS

```
Step 0 (done) → Step 1 → Steps 2+3 (parallel) → Step 4 (checkpoint) → Step 5 → Step 6 → Step 7
```

### STEP 0 — Human reference scaffold ✅ DONE
67 genes compiled, UniProt verified, Ensembl Compara orthologues checked for chicken + duck, InterPro domain verification run, clean summary produced.

### STEP 1 — Foldseek domain reference panel → NEXT
Select best 3D structure per domain family: CARD, PYD, DD, DED, FIND, gasdermin-N, RHIM, BIR, BH3.
- PDB experimental > AFDB predicted; highest pLDDT over the domain region specifically
- Isolated domain structures for Arm 2 (novel discovery); full-length for Arm 1 (orthologue recovery)
- Output: directory of .pdb/.cif files + metadata table (domain family, source protein, resolution/pLDDT)
- Local work, ~2–3 days

### STEP 2 — BLAST baseline
- Build one taxid-tagged avian protein DB; blastp 67 queries against it → 1 Slurm job
- Per-species results come from staxids column — no per-species job arrays needed
- Defines what's findable by sequence; Foldseek recovers what BLAST misses

### BUILD afdb_aves (prerequisite for Step 3)
- Do NOT fetch structures per-accession (millions of API calls — wrong approach)
- Use Foldseek's prebuilt taxonomy-tagged database:
  ```
  foldseek databases Alphafold/UniProt resources/afdb_uniprot tmp
  foldseek filtertaxseqdb resources/afdb_uniprot/afdb resources/afdb_aves/afdb_aves --taxon-list 8782
  rm -rf resources/afdb_uniprot   # reclaim ~1TB transient disk
  ```
- Aves taxon ID = 8782 (includes all descendant species automatically)
- Use Alphafold/UniProt (not UniProt50) to preserve per-species resolution
- Watch transient disk peak on RDS (5TB, no backup)

### STEP 3 — Foldseek structural screen (2 arms)
- Arm 1: 67 full-length human structures vs afdb_aves → 1 Slurm job
- Arm 2: ~9 isolated domain structures vs afdb_aves → 1 Slurm job
- Output carries taxid/taxname columns; split by species in post-processing
- Three-condition novelty bar for Arm 2 hits: high TM-score + no sequence homology + no synteny

### STEP 4 — Synteny validation (Snakemake checkpoint)
- Validates "absent" calls and confirms Arm 2 candidates are genuinely novel (not syntenic to any human gene)
- Partly manual via Genomicus/Ensembl — triage programmatically, eyeball ambiguous cases
- Decision rule: scaffold <N kb → "inconclusive", not "absent"
- This is a Snakemake checkpoint: gates which genes proceed to Step 5

### STEP 5 — Phylogenetics and selection (fan-out across genes)
- Per gene: MAFFT alignment → trimAl → IQ-TREE 2 (ModelFinder + 1000 bootstrap) → HyPhy BUSTED/MEME
- Loss dating: map gene loss events onto Stiller 2024 species tree
- Requires ≥10 species per gene for selection analysis — transparency required for genes below threshold
- Slurm job array, one job per gene — the only truly parallel step

### STEP 6 — Expression overlay
- Stella's HD-11 RNAseq (chicken, IAV, in hand) + Elleder's libraries (pending)
- Supporting evidence only — never a discovery criterion
- A novel candidate absent from macrophage expression still appears in the atlas

### STEP 7 — Integration and per-gene tiering
- Tier 1: Conserved (BLAST + synteny + expressed)
- Tier 2: Diverged (Foldseek Arm 1, poor BLAST, synteny confirmed)
- Tier 3: Validated absent (no BLAST, no Foldseek, syntenic region confirmed empty)
- Tier 4: Novel candidate (Arm 2, three-condition filter passed)
- Tier 5: Inconclusive (bad assembly, no synteny resolution)
- This species × gene × evidence matrix is the main figure of Paper 1

---

## 4. SNAKEMAKE + CSD3 DESIGN (plain-language summary)

**What Snakemake does:** You write each pipeline step as a "rule" — what it needs and what it produces. Snakemake works out the running order, runs independent steps in parallel, and only reruns steps whose inputs changed. You never write "do step 1 then step 2" — dependencies are inferred automatically.

**DAG** = Directed Acyclic Graph. Just a dependency map: "to make this, I first need that." The pipeline diagram is a DAG. Snakemake builds it internally from your rules.

**Slurm** = CSD3's queue manager. You submit jobs; it decides when and where they run on the cluster. Never poll it in a loop. Snakemake's executor plugin handles all status-checking automatically using sacct (adaptive back-off, 40s–180s intervals) — as long as you let Snakemake manage jobs and never write your own polling wrapper, you're safe by design.

**Login node vs compute node:** When you log into CSD3 you land on the login node — fine for editing files, moving data, submitting jobs. Heavy analysis (BLAST, Foldseek, IQ-TREE) must go to compute nodes via Slurm submission.

**Conda env per tool group** (not one giant env): foldseek.yaml, blast.yaml, phylo.yaml (iqtree2+mafft+trimal), hyphy.yaml, pyanalysis.yaml. Snakemake activates the right one per rule.

**Profile** = saved settings file for CSD3 (`config/csd3/profile.yaml`). Includes partition, memory defaults, account string, max concurrent jobs, status-check rate. Written once; Snakemake reads it every run.

**Checkpoint** = a deliberate pause where Snakemake looks at results so far and decides what to do next. Step 4 (synteny) is a checkpoint: it writes the validated gene list, and only genes with enough species proceed to Step 5 phylogenetics. You don't manually maintain that list — the pipeline decides at runtime.

---

## 5. SCRIPTS WRITTEN AND RUN

All scripts in `workflow/scripts/atlas/` (renamed from track_a this session).

| Script | Purpose | Status |
|--------|---------|--------|
| 01_verify_uniprot.py | Verify UniProt accessions | ✅ Done |
| 02_bird_orthologue_check_v2.py | Ensembl Compara orthologues chicken + duck | ✅ Done |
| 02_diagnostic.py | Diagnosed Ensembl endpoint bug | ✅ Done (kept for reference) |
| 04_domain_verification.py | InterPro domain check for bird orthologues | ✅ Done |
| 06_merge_summary.py | Clean color-coded summary xlsx | ✅ Done |

Scripts 03 and 05 were superseded. Script 04 takes file path as argument — path-safe, won't break after reorganisation.

**Key bug fixes (already applied):**
- Ensembl REST: `/homology/id/homo_sapiens/{id}` not `/homology/id/{id}`
- Duck species: `anas_platyrhynchos_platyrhynchos` (CAU_duck1.0) not `anas_platyrhynchos`
- Use `python` not `python3` inside conda aviced env

---

## 6. CURRENT STATE OF THE ORTHOLOGUE TABLE

**Chicken (67 genes):** 44 PRESENT / 3 DUBIOUS (CASP1, CASP4, CASP5) / 3 NOISE (<10% id, treat as absent) / 17 ABSENT
**Duck (67 genes):** 37 PRESENT / 4 DUBIOUS (CASP1, CASP4, CASP5, TRIF) / 2 NOISE / 24 ABSENT

**Domain verdicts (chicken):** 22 FULL_MATCH, 15 GOOD_MATCH, 4 PARTIAL_MATCH, 2 NO_MATCH, 7 NO_UNIPROT, 17 NO_ORTHOLOGUE
**Domain verdicts (duck):** 19 FULL_MATCH, 7 GOOD_MATCH, 2 PARTIAL_MATCH, 2 WEAK_MATCH, 13 NO_UNIPROT, 24 NO_ORTHOLOGUE

**NO_UNIPROT cases (20 total):** Ensembl→UniProt xref failed. Orthologues are real. Fix: search UniProt directly (species + gene name). Notable: CASP7 chicken (71.6% id), CASP9 chicken (51%), MDA5/IFIH1 chicken (60.2%), BCL2L1 duck (63.5%).

**Suspicious absences needing BLAST on CSD3:** BAX, MAVS, ASC/PYCARD, BID, RBCK1/HOIL-1

**Chicken-duck discrepancies:** BIM (duck present, chicken absent), MLKL (chicken present, duck absent — annotation gap confirmed), GSDMA/B/C/D (all need manual gasdermin curation), CARD8, RNF31/HOIP.

**DDX58/RIG-I:** Confirmed from literature — chicken ABSENT, duck PRESENT. Add manually. Use as positive control, not a result.

---

## 7. REPO REORGANISATION — COMPLETED THIS SESSION

**Before:** folders organised by data type (track_a_rnaseq, track_b_genomics, track_c_alphafold) — didn't match the three-paper structure.

**After:** organised by purpose. Committed and pushed to GitHub (origin/main).

```
aviced_bioinformatics/
├── README.md
├── analysis/
│   ├── coevolution/         # Paper 3 (future)
│   └── kg/                  # Paper 2 (future)
├── config/                  # Snakemake settings + CSD3 profile (to build)
├── docs/
│   └── pipeline_decisions.md
├── notebooks/               # exploratory only
├── resources/
│   ├── reference/           # curated inputs tracked by Git
│   │   ├── human_cell_death_gene_universe_v1.xlsx  ← moved here (was in results/)
│   │   ├── orthologue_domain_summary.xlsx           ← moved here (was in results/)
│   │   ├── ncbi_gene_search_results.tsv
│   │   └── uniprot_verification.tsv
│   ├── structures/          # AFDB query structures + afdb_aves target (future)
│   ├── genomes/             # ignored by Git (large files)
│   └── rnaseq/              # ignored by Git (large files)
├── results/                 # pipeline outputs — ignored by Git, regeneratable
│   ├── atlas/  blast/  foldseek/  synteny/  phylo/  expression/
│   └── [intermediate xlsx files remain here, ignored]
└── workflow/
    ├── Snakefile
    ├── envs/
    │   └── aviced_env.yml
    ├── rules/               # Snakemake rules (to build, one per step)
    └── scripts/
        ├── atlas/           # Paper 1 scripts (was track_a)
        ├── coevolution/     # Paper 3 scripts (was track_b)
        └── kg/              # Paper 2 scripts (was track_c)
```

**Important:** `results/` is still Git-ignored (correct — pipeline outputs shouldn't be in Git). The two curated reference tables moved to `resources/reference/` where they ARE tracked by Git and protected.

**Backup:** Full `results/` folder copied to Google Drive before reorganisation.

---

## 8. WHAT TO DO NEXT (priority order)

**Next session — Step 1: build the Foldseek domain reference panel (local work)**
For each domain family, select the single best representative 3D structure:
- CARD, PYD, DD, DED, FIND, gasdermin-N, RHIM, BIR, BH3
- Rule: PDB experimental > AFDB predicted; best resolution / highest pLDDT over the domain region
- Isolated domain structures → Arm 2 queries
- Full-length structures → Arm 1 queries
- Output: `resources/structures/domain_panel/` + metadata table

**Before CSD3 work (can do locally):**
- Resolve 20 NO_UNIPROT cases by searching UniProt directly (species + gene name)
- Manually curate gasdermin assignments (which chicken gene(s) do GSDMA/B/C/D all map to?)
- Add DDX58/RIG-I manually to the orthologue table
- BLAST check suspicious absences in NCBI Gene: BAX, MAVS, ASC, BID, RBCK1

**CSD3 setup (when ready for Steps 2–3):**
- Build `config/csd3/profile.yaml` (Snakemake cluster profile)
- Write Snakemake rules for Steps 2 (BLAST) and the afdb_aves build
- Test with a dry-run (`snakemake -n`) before submitting real jobs
- Never run heavy analysis on the login node

---

## 9. KEY DECISIONS AND PRINCIPLES (permanent reference)

- **Three-condition novelty bar is non-negotiable** for Arm 2: high TM-score + no sequence homology + no synteny. All three must hold simultaneously.
- **Expression is supporting evidence only** — never a discovery criterion.
- **Ensembl "absent" ≠ biologically absent** — it means the orthologue pipeline couldn't assign. Each layer (Ensembl → BLAST → Foldseek) recovers what the previous missed.
- **Database annotation ≠ biological absence** — NCBI-only hits reflect annotation gaps, not true absence; synteny validation required.
- **Chickens ≠ ducks for pathogen tolerance** — chickens are susceptible to highly pathogenic AIV; ducks are the AIV-tolerant example. This distinction matters for framing.
- **DAG** = Directed Acyclic Graph — the pipeline's dependency map. Snakemake builds it from your rules automatically.
- **git mv** for all file moves inside the repo (not Finder/drag-and-drop) — keeps Git's rename history clean.
- **results/ is Git-ignored** (regeneratable outputs). **resources/reference/ is Git-tracked** (curated inputs that can't be regenerated by running a script).
