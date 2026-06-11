# AviCeD Bioinformatics — Session Summary
## Date: 11 June 2026
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
- General aim: Map the co-evolutionary dynamics of protease pocket specificity and substrate tetrapeptide motifs across birds to explain the GSDMD→GSDMA substrate switch and the IL-1β processing divergence.
- Specific aim 1: Active site pocket profiling — AlphaFold structures of caspase-1 across birds, S1–S4 subsite classification, compare pocket geometry for YVAD compatibility. Tools: AlphaFold structures, PyMOL/ChimeraX, fpocket or SiteMap.
- Specific aim 2: Substrate motif scanning — scan GSDMA, GSDME, IL-1β, IL-18 sequences across birds for tetrapeptide cleavage motifs, classify as YVAD-like, FASD-like, novel, or degraded. Tools: custom scripts, MAFFT alignments.
- Specific aim 3: Co-evolutionary analysis — paired dN/dS, mirror tree / DCA for protease-substrate co-evolution signal, test whether species where caspase-1 is lost show accelerated drift in substrate motifs. Tools: HyPhy, BayesTraits, DCA tools (plmDCA or GREMLIN).
- Uses only public data + AlphaFold structures — no dependency on Paper 1
- May fold into Paper 1 or stand alone depending on signal strength
- Target: MBE standalone or a chapter in Paper 1
- Strength: genuinely novel computational biology — nobody has done systematic protease-substrate co-evolution across avian phylogeny for cell death proteins
- Weakness: might be too narrow for standalone paper — decision made after seeing the data

**Paper 2 — Knowledge Graph / Systems Biology (starts after Paper 1 atlas)**
- Research question: Can interaction rules learned from well-characterised mammalian cell death networks predict how avian cell death pathways are wired, including functional substitutions for absent human orthologues?
- General aim: Build a cross-species knowledge graph of cell death protein interactions, train a link prediction model on human/mouse/zebrafish, and predict novel avian pathway wiring.
- Specific aim 1 — Knowledge graph construction: Nodes = all cell death + innate immunity proteins in human, mouse, zebrafish (from Reactome, STRING physical ≥0.7, IntAct, BioGRID). Node features: ESM2 embeddings (650M model), InterPro domain architecture, tissue/stimulus expression profiles. Edges typed: cleaves, activates, scaffolds, inhibits, phosphorylates. Sourced from Reactome (curated) + literature curation. Tools: ESM2 (fair-esm), PyKEEN, NetworkX, STRING/Reactome APIs.
- Specific aim 2 — Model training and validation: Train RotatE or ComplEx on mammalian + zebrafish KG. Validation: hold out 20% known edges, test link prediction recovery (hits@10, MRR). Biological validation: does model recover known interactions not in training set?
- Specific aim 3 — Avian wiring prediction: Add bird proteins from Paper 1 atlas as new nodes (same feature pipeline). Predict edges: which avian proteins interact, with what relation type? Rank by confidence.
- Specific aim 4 — Validation against experimental data: Cross-reference with Boucher's cleavage assays, Stella's KO phenotypes (RIPK1KO, CASP8KO), Digital Fowl simulations. Failed predictions flagged as biology, not model error.
- Target: Mol Sys Bio or PLOS Comp Bio
- Strength: genuinely novel application — nobody has applied KG link prediction to comparative innate immunity
- Weakness: depends entirely on Paper 1 being done first; curation bottleneck building the triple table; purely computational prediction paper faces sceptical immunology reviewers without experimental validation
- Critical question: whether KG is standalone paper or embedded in larger systems biology paper with graph topology analysis

**Order: Paper 1 first (empirical foundation) → Paper 3 in parallel → Paper 2 after atlas is done. KG model can be built and validated on human/mouse/zebrafish during Paper 1 work.**

---

## 2. CONCEPTUAL GROUNDWORK — EMBEDDINGS, TRANSFORMERS, KNOWLEDGE GRAPHS

We covered in depth:

**Embeddings:**
- A list of numbers that describes something. ESM2 produces 1,280 numbers per protein.
- Learned as side effect of fill-in-the-blank training: mask 15% of amino acids, predict what's hidden. To get good at this, the model develops internal representations that capture structural and functional properties.
- The embedding and the blank-filling emerge together. Day zero: random vectors. After millions of training steps: meaningful representations that cluster functionally related proteins together.
- ESM2 was trained on 250M protein sequences from UniRef50. Raw amino acid strings only. No 3D structure, no domain annotation, no GO terms, no species labels. Everything it knows was inferred from learning to predict masked residues.
- Difference from BLAST: BLAST gives one number (score) comparing two sequences. Embedding gives 1,280 features per protein capturing deep properties invisible to sequence alignment. Two proteins at 15% identity (BLAST says "unrelated") can sit next to each other in embedding space.

**How ESM2 relates to AlphaFold:**
- ESM2 asks: what is each residue's role in context? Output: vector per protein. It's a description.
- AlphaFold asks: what are the 3D coordinates? Output: structure with pLDDT/PAE. It's a prediction of physical shape.
- Both use transformer architecture. AlphaFold trained on ~170K experimental structures from PDB. ESM2 trained purely on sequences with no structural supervision, yet its representations encode structural information.
- For the project: AlphaFold gives 3D structures for Foldseek. ESM2 gives embedding vectors for KG node features. Different tools, complementary roles.

**How all relate to Claude/ChatGPT:**
- Same transformer architecture, different training data and objectives.
- Claude/ChatGPT: sequences of text tokens, predict next word, learns grammar/facts/reasoning.
- ESM2: sequences of amino acid tokens, predict masked residues, learns protein grammar/motifs/constraints.
- AlphaFold: MSA tokens + sequence, predict 3D coordinates of known structures.
- Key insight: sequences are sequences — attention mechanism doesn't care what tokens mean, learns patterns of co-occurrence and context-dependency.

**Transformers:**
- A model architecture (not an algorithm, not software — a blueprint for information flow).
- Core innovation: attention mechanism. Every position can "talk to" every other position simultaneously. Each residue computes Query ("what am I looking for?"), Key ("what do I offer?"), Value ("what info do I carry?"). High Query-Key match = high attention = relevant positions share information regardless of distance.
- Stack 33 attention layers (ESM2 t33): early layers capture local chemistry, middle layers secondary structure, late layers domain/whole-protein function.
- Hierarchy: Machine Learning > Deep Learning > Architecture (transformer is one; CNN, RNN are others) > Model (specific trained instance: ESM2, Claude, AlphaFold) > Algorithm (training procedure: SGD, Adam, backpropagation).

**RotatE/ComplEx:**
- NOT transformers. Much simpler embedding models for knowledge graphs.
- Thousands to low millions of parameters (not billions). Train in minutes.
- RotatE represents each relation as a rotation in complex vector space. If "CASP1 cleaves GSDMD" is real, rotating CASP1 vector by "cleaves" rotation lands near GSDMD vector.
- These learn interaction rules between protein descriptions. Transformers (ESM2) produce the descriptions; RotatE/ComplEx learn the relationships.

**Key insight for Paper 2:**
- Julieth's original dream: "train a model to learn cell death pathway wiring from scratch."
- Problem: n=1 (human) is insufficient. Can't learn rules of pathway assembly from one example.
- Reframed to: "learn general protein interaction grammar from all of protein biology across well-characterised species (human/mouse/zebrafish), then apply to predict bird wiring specifically."
- This works because interaction rules (complementary domain interfaces, CARD-CARD recruitment, zymogen activation) aren't unique to cell death — they're general protein interaction rules across thousands of pathways.
- The honest framing — transferring interaction knowledge across species — is actually stronger than "discovers rules from scratch" because it's defensible.

**Compute requirements (reassuringly small):**
- ESM2 embedding of all ~500-800 proteins: ~30 min on single GPU, or overnight on laptop CPU. Run once, save as numpy array.
- RotatE/ComplEx training on ~5K-10K nodes, ~50K-100K edges: 10-30 min on single GPU, 1-2 hours on CPU. <4GB RAM.
- The bottleneck is curation (weeks of human time), not compute.

---

## 3. PAPER 1 PIPELINE — SEVEN STEPS WITH FULL DETAIL

```
Step 0 → Step 1 (parallel with Step 2)
Step 0 → Step 2
Steps 1+2 → Step 3
Steps 2+3 → Step 4
Steps 2+3+4 → Step 5
Steps 3+4+5 → Step 6
Everything → Step 7
```

### STEP 0 — Human reference scaffold ✅ DONE

**What you do:** Compile the complete human cell death gene universe with UniProt accessions, one canonical isoform per gene.

**Input:** Literature (Newton 2024 compendium, Broz/Dixit reviews, Man & Kanneganti 2024). Your own knowledge of the field.

**Output:** A table. One row per gene. Columns: gene symbol, UniProt accession, AFDB structure ID, pathway membership (pyroptosis / necroptosis / apoptosis / PANoptosis hub / sensor), domain architecture (from InterPro).

**Feeds into:** Step 1, Step 2, Step 4, Step 5. Everything downstream depends on this being complete and correct. If you miss a gene here, it's invisible for the rest of the paper.

**Critical check:** Verify against Newton 2024 — if they include something you don't, you need a reason for excluding it.

**Status:** DONE. 67 genes compiled, UniProt accessions verified (script 01), Ensembl Compara orthologues checked for chicken and duck (script 02), InterPro domain verification run (script 04), clean summary produced (script 06). See section 5 below for full results.

### STEP 1 — Foldseek domain reference panel → NEXT

**What you do:** For each domain family in your scaffold (CARD, PYD, DD, DED, FIND, gasdermin-N, RHIM, BIR, BH3), select the single best representative 3D structure. This is your query set for all Foldseek searches.

**Input:** Step 0 table + PDB/AFDB. For each domain, pick the structure with best resolution (PDB experimental > AFDB predicted) and highest pLDDT for the domain region specifically.

**Output:** A directory of structure files (.pdb or .cif), one per domain. Plus a metadata table mapping each structure to its domain family, source protein, and confidence metrics.

**Feeds into:** Step 3 (Foldseek Arm 1 and Arm 2).

**Critical decision:** Do you use full-length protein structures or isolated domain structures as queries? Isolated domains give cleaner hits (no off-target matching to unrelated regions of multi-domain proteins). Full-length structures might catch domain arrangements you didn't anticipate. Decision made: isolated domains for Arm 2 (novel discovery) and full-length for Arm 1 (orthologue recovery). Different queries, same target databases.

**Time:** 2–3 days. Requires structural biology judgment, not just downloading.

**Depends on:** Step 0 only.

### STEP 2 — Sequence-based orthologue mapping (BLAST arm)

**What you do:** For every gene in Step 0, run BLAST (protein) against all ~50 avian proteomes. This is your classical homology search — fast, cheap, and catches everything that's still recognisable by sequence.

**Input:** Step 0 protein sequences + avian proteomes from NCBI/Ensembl/UniProt.

**Output:** Per gene × per species: best hit, e-value, percent identity, coverage. Many genes will have clear orthologues (CASP8, CASP3, TRIF). Some will return nothing (ZBP1, RIPK3 in chicken). Some will return ambiguous hits (CASP1, ASC).

**Feeds into:** Step 4 (synteny validation — the ambiguous and absent cases need synteny). Step 5 (the orthologue sequences feed the phylogenetics). Step 3 (tells you which genes Foldseek needs to recover that BLAST missed).

**Critical weakness of this step alone:** BLAST misses diverged orthologues. That's the whole justification for Foldseek. But BLAST is still essential because it defines the baseline: what's findable by sequence, what isn't, and therefore what Foldseek adds. If you skip BLAST, a reviewer will ask why you didn't try the simple approach first.

**Time:** Computationally trivial. 1 day for setup, minutes to run on CSD3. The real time is parsing and curating the output table.

**Depends on:** Step 0 only. Runs in parallel with Step 1.

### STEP 3 — Foldseek structural screen (both arms)

**What you do:** Two separate searches against the same target databases (~50 avian AFDB proteomes).

**Arm 1 (reference-anchored):** Query = full-length AFDB structures of every Step 0 gene. Target = avian proteomes. Purpose: recover diverged orthologues that BLAST missed. A hit here with high TM-score but low sequence identity is a diverged orthologue.

**Arm 2 (unanchored novel discovery):** Query = isolated domain structures from Step 1. Target = same avian proteomes. Purpose: find proteins with death-related domains that have no orthologue relationship to any human gene. The three-condition filter applies: high TM-score + no sequence homology (from Step 2) + no syntenic relationship (from Step 4).

**Input:** Step 1 (domain structures), Step 0 (full-length structures), avian AFDB proteomes.

**Output:** Two hit tables. Arm 1: gene × species × TM-score × sequence identity × target protein ID. Arm 2: domain × species × TM-score × target protein ID × domain architecture of the hit.

**Feeds into:** Step 4 (Arm 1 hits need synteny validation; Arm 2 hits need synteny exclusion). Step 5 (Arm 1 confirmed orthologues feed phylogenetics). Step 6 (Arm 2 validated novel candidates feed the final integration). Step 7 (both arms feed expression overlay).

**Critical risk:** Foldseek on ~50 proteomes is substantial compute. Each avian AFDB proteome is ~15,000–25,000 structures. This is a CSD3 job, not your laptop. Estimate 1–2 hours per proteome per query set. With ~50 proteomes and two query sets, that's ~100–200 CPU-hours. Submit as a Slurm job array, one job per species.

**Critical analytical risk:** False positives in Arm 2. Death-related domains (especially DD superfamily) are structurally similar to non-death domains. You will get hits to proteins that have death-fold-like topology but no biological role in cell death. The pLDDT filter helps but doesn't solve this. Domain architecture filtering (does the hit have other domains consistent with cell death function?) is essential.

**Depends on:** Steps 0, 1, and 2 (you need the BLAST results to classify Arm 1 vs Arm 2 hits).

### STEP 4 — Synteny validation

**What you do:** For every gene called "absent" by BLAST (Step 2) and every ambiguous case, check the syntenic region. Are the flanking genes conserved? Is there a degraded pseudogene remnant? Or is the locus genuinely empty?

For Arm 2 novel candidates: confirm they are NOT in a syntenic position relative to any human cell death gene. If they are, they're diverged orthologues, not novel genes. Reclassify them to Arm 1.

**Input:** Step 2 (absent/ambiguous gene list), Step 3 Arm 2 hits, Ensembl/UCSC/Genomicus genome browsers.

**Output:** Per gene × per species verdict: validated absent (syntenic region found, gene missing) / annotation gap (syntenic region has a predicted gene not in UniProt) / pseudogene remnant / reclassified to Arm 1.

**Feeds into:** Step 5 (validated absences are dated on the phylogeny — when was the gene lost?). Step 6 (final tiering). Step 3 feedback (Arm 2 hits reclassified to Arm 1).

**Critical bottleneck:** This is partly manual. Genomicus and Ensembl Compara help, but for poorly assembled genomes (many non-model birds), the syntenic region may be on a scaffold too short to confirm flanking genes. You need a decision rule: if the flanking gene scaffold is <N kb, the synteny call is "inconclusive" not "absent."

**Depends on:** Steps 2 and 3.

### STEP 5 — Phylogenetics and selection

**What you do:** For every gene confirmed present in ≥10 species (from Steps 2 + 3 Arm 1), build a gene tree and run selection analyses.

**Input:** Orthologue sequences from Steps 2 and 3 Arm 1. Species tree (published avian phylogeny, e.g. Jarvis 2014 or Stiller 2024).

**Sub-step 5a — Alignment:** MAFFT on the orthologue set per gene. Trim with trimAl or BMGE. Inspect in Jalview — automated trimming can destroy informative columns in fast-evolving genes.

**Sub-step 5b — Gene trees:** IQ-TREE 2, ModelFinder for substitution model selection, ultrafast bootstrap (1000 replicates). Compare gene tree to species tree — discordance indicates incomplete lineage sorting, duplication, or lateral transfer.

**Sub-step 5c — Selection:** HyPhy BUSTED (is there episodic positive selection anywhere on the tree?) and MEME (which specific sites are under selection?). Run on the Datamonkey server or locally on CSD3.

**Sub-step 5d — Loss dating:** For genes absent in some lineages, map the loss event onto the species tree. This tells you: did birds lose GSDMD once (ancestral loss) or multiple times (convergent loss)?

**Output:** Per gene: gene tree, selection results (which sites, which branches), loss dating. This is the evolutionary narrative of the paper.

**Feeds into:** Step 7 (evolutionary context for the final atlas). This is where the paper becomes more than a table — the evolutionary story is what makes it MBE rather than a database paper.

**Critical weakness:** Selection analysis requires adequate taxon sampling. Genes present in only 3–4 bird species won't have enough power for BUSTED/MEME. You'll know this only after Steps 2–4 produce the orthologue counts. For rare genes, report the phylogeny without selection analysis and be transparent about it.

**Depends on:** Steps 2, 3, 4. This is the latest step to start and the most computationally intensive (IQ-TREE on 50 taxa × 40 genes = hundreds of tree searches).

### STEP 6 — Expression overlay

**What you do:** Check whether Arm 2 novel candidates and diverged Arm 1 orthologues are actually expressed in macrophages, especially under infection.

**Input:** Candidate list from Steps 3–5. Stella's HD-11 RNAseq (chicken, IAV, in hand). Elleder's libraries (chicken + duck, Salmonella + IAV, pending).

**Output:** Per candidate: expressed yes/no, differentially expressed under infection yes/no, fold change, adjusted p-value. This is supporting evidence, not a discovery criterion. A novel candidate that isn't expressed in macrophages isn't necessarily wrong — it might function in a different cell type.

**Feeds into:** Step 7 (final integration).

**Critical rule we already set:** Expression does NOT gate inclusion in the atlas. A structurally validated novel candidate with no macrophage expression still appears in the paper — you just note the expression status. This prevents you from discarding real genes because of limited tissue sampling.

**Depends on:** Steps 3–5 for the candidate list. Stella's data for chicken. Elleder's data for duck (external dependency, not in your control).

### STEP 7 — Integration and per-gene tiering

**What you do:** Combine all evidence into the final species × gene × evidence matrix.

**Input:** Everything upstream.

**Output per gene per species:**
- **Tier 1 — Conserved:** BLAST hit + synteny confirmed + expressed
- **Tier 2 — Diverged:** Foldseek Arm 1 hit (poor BLAST) + synteny confirmed + expressed
- **Tier 3 — Validated absent:** No BLAST + no Foldseek + syntenic region confirmed empty/pseudogene
- **Tier 4 — Novel candidate:** Foldseek Arm 2 hit + three-condition filter passed + expression status noted
- **Tier 5 — Inconclusive:** Insufficient data (bad genome assembly, no synteny resolution)

This matrix is the main figure of the paper. Everything else (phylogenies, selection, expression) is supporting context.

**Depends on:** All steps.

---

## 4. SCRIPTS WRITTEN AND RUN

All scripts in the user's repo at `~/aviced_bioinformatics/`. Run with `python` (not `python3`) inside conda `aviced` environment.

| Script | Purpose | Status |
|--------|---------|--------|
| 01_verify_uniprot.py | Verify UniProt accessions against REST API | ✅ Run. 62/67 OK, 5 alias mismatches (expected: ASC→PYCARD, BIM→BCL2L11, DDX58→RIGI, PYRIN→MEFV, TRIF→TICAM1), 1 wrong accession fixed (DFNB59: Q8N8I3 → Q0ZLH3, the old accession was SPTY2D1) |
| 02_bird_orthologue_check_v2.py | Ensembl Compara orthologues for chicken + duck | ✅ Run. Two bugs fixed: (1) `anas_platyrhynchos` → `anas_platyrhynchos_platyrhynchos` (old BGI assembly frozen at release 80 vs CAU_duck1.0 current in Compara); (2) `homo_sapiens` added to homology endpoint path (`/homology/id/homo_sapiens/{id}` not `/homology/id/{id}`). Chicken: 50 present, 16 absent, 1 not in Ensembl. Duck: 43 present, 23 absent, 1 not in Ensembl. |
| 03_chicken_to_duck_orthologues.py | Check chicken→duck orthologues | Written but superseded by fix to script 02 |
| 04_domain_verification.py | InterPro domain check for bird orthologues | ✅ Run. Maps Ensembl ID → UniProt (via Ensembl xrefs, prefers Swiss-Prot over TrEMBL) → InterPro domains → compares with human. Added columns to the bird_orthologues xlsx. |
| 05_orthologue_domain_summary.py | All-in-one: reads xlsx, runs domain verification, outputs clean table | Written but not needed (04 already produced the data) |
| 06_merge_summary.py | Merges existing data into clean color-coded summary (no API calls) | ✅ Run. Output: orthologue_domain_summary.xlsx |

**Key bug fixes during session:**
- Ensembl REST homology endpoint needed species in path: `/homology/id/homo_sapiens/{id}` not `/homology/id/{id}` — diagnosed via 02_diagnostic.py which showed 404 on the homology call while lookup worked fine
- Duck species name: `anas_platyrhynchos_platyrhynchos` (CAU_duck1.0, current Compara) not `anas_platyrhynchos` (old BGI assembly, frozen at release 80, legacy, not in current Compara). Julieth found this by reading the Ensembl Compara documentation and TreeFam docs.
- `python` vs `python3` in conda: use `python` inside the aviced conda env — `python3` points to system Python without openpyxl
- DFNB59 accession Q8N8I3 was actually SPTY2D1 (a splicing factor). Correct accession: Q0ZLH3 (PJVK_HUMAN)

---

## 5. CURRENT STATE OF THE ORTHOLOGUE TABLE

**File: `results/orthologue_domain_summary.xlsx`** — the clean deliverable with color coding.
**Full data: `results/human_cell_death_gene_universe_v1_bird_orthologues_domains.xlsx`** — all 36 columns.

**Chicken (67 genes):**
- 44 PRESENT (≥20% identity)
- 3 DUBIOUS (10-20%): CASP1 (14.6%), CASP4 (14.6%), CASP5 (13.6%) — real biology, single diverged bird caspase mapping many2many to all three human inflammatory caspases
- 3 NOISE (<10%): NLRC4 (3.3%), TNFRSF10A (1.3%), TNFRSF10B (7.5%) — Ensembl misassignments confirmed by NO_MATCH domain verdict, treat as absent
- 17 ABSENT

**Duck (67 genes):**
- 37 PRESENT
- 4 DUBIOUS: CASP1 (15.1%), CASP4 (14.6%), CASP5 (13.1%), TRIF (19.0%)
- 2 NOISE: NLRC4 (2.5%), TNFRSF10B (7.3%)
- 24 ABSENT

**Domain verification verdicts (chicken):**
- 22 FULL_MATCH, 15 GOOD_MATCH, 4 PARTIAL_MATCH, 2 NO_MATCH, 7 NO_UNIPROT, 17 NO_ORTHOLOGUE

**Domain verification verdicts (duck):**
- 19 FULL_MATCH, 7 GOOD_MATCH, 2 PARTIAL_MATCH, 2 WEAK_MATCH, 13 NO_UNIPROT, 24 NO_ORTHOLOGUE

**NO_UNIPROT cases (20 total):** Ensembl→UniProt xref mapping failed. Orthologues are real (confirmed by Compara). Fixable by searching UniProt directly with species + gene name. Notable cases: CASP7 chicken (71.6% id), CASP9 chicken (51%), IFIH1/MDA5 chicken (60.2%), BCL2L1 duck (63.5% — Swiss-Prot Q07816 confirmed manually).

**Chicken-duck discrepancies (one present, one absent):**
- BIM: duck present (37.9% id), chicken absent — if real, chickens lost a BH3-only activator that ducks retained
- MLKL: chicken present (40.8%), duck absent — confirmed annotation gap, TrEMBL F6S337 exists for duck MLKL
- GSDMA/B/C/D: all "present" in chicken (one2many at 20-30% identity, all mapping to same chicken gasdermin(s)), all absent in duck — needs manual gasdermin curation for both species
- CARD8: chicken present (25.5%), duck absent — low identity, needs domain validation (FIND + CARD?)
- RNF31/HOIP: chicken present (24%), duck absent

**Suspicious absences in BOTH species (need BLAST on CSD3):**
- BAX — fundamental apoptosis effector, present in fish. Almost certainly present in birds but too diverged for Ensembl's orthologue pipeline.
- MAVS — required for MDA5 signalling, both species have MDA5. Without MAVS, MDA5 cannot signal. Either annotation gap or genuinely important rewiring.
- ASC/PYCARD — chicken candidate XP_015129143.1 exists with PYD and DD domains. Ensembl can't assign orthology.
- BID — bridges extrinsic→intrinsic apoptosis. Small, fast-evolving BH3-only protein.
- RBCK1/HOIL-1 — LUBAC component. The other two LUBAC subunits (RNF31/HOIP, SHARPIN) are both present. LUBAC without HOIL-1 is non-functional in mammals.

**BH3-only proteins pattern:** BAD, BID, BIM, BBC3/PUMA, PMAIP1/NOXA all absent in Ensembl Compara. These are small (~100-200 aa), intrinsically disordered, fast-evolving. Ensembl's pipeline struggles with them. Some may be genuinely absent (reduced BH3-only repertoire in birds) but don't trust Ensembl's "absent" call without BLAST + synteny verification.

**DDX58/RIG-I:** Not in Ensembl (HGNC updated symbol to RIGI, Ensembl hasn't caught up). Confirmed from literature: chicken ABSENT, duck PRESENT. The canonical chicken vs duck difference for AIV tolerance. Add manually.

---

## 6. KEY DECISIONS MADE

- Duck assembly: CAU_duck1.0 (in Ensembl, has Compara, assembly GCA_002743455.1) for orthologue queries. SKLA1.0 (better contiguity, complete MHC, built for influenza/MHC work) potentially for CSD3 pipeline. Can use both — discrepancies between assemblies are informative.
- Duck proteome: UniProt reference proteome UP000016666, taxon 8840, 27,053 entries, 16,599 genes. BUSCO 89.7% complete.
- AFDB duck proteome: needs verification. Check https://alphafold.ebi.ac.uk/download for UP000016666. The AlphaFold API script returned 400 — check download page manually. 
- Expression is supporting evidence only, never a discovery criterion.
- Three-condition novelty bar is non-negotiable for Arm 2 candidates (structural similarity + no sequence homology + no synteny).
- OrthoFinder not needed yet — save for pan-avian expansion in Step 5.
- Domain verification (InterPro) is essential complement to % identity for orthologue calls — % identity alone doesn't tell you functional equivalence.
- Ensembl Compara "absent" ≠ biologically absent. It means orthologue pipeline couldn't confidently assign. Each pipeline layer (Ensembl → BLAST → Foldseek) recovers what the previous missed.

---

## 7. WHAT TO DO NEXT

**Immediate (Step 1):**
1. Build the Foldseek domain reference panel — select best 3D structure per domain family
2. Structural biology judgment required: PDB experimental > AFDB predicted, highest pLDDT for domain region
3. Isolated domain structures for Arm 2, full-length for Arm 1

**Before CSD3 work:**
- Verify AFDB has duck proteome (download page check). update: Alphafold only has 48 species with their downloadable folded proteome. Shall we have fetch AlphaFold structures for each ID?
- Resolve NO_UNIPROT cases by direct UniProt search (species + gene name)
- Manually curate gasdermin assignments (which chicken gene(s) do GSDMA/B/C/D all map to?)
- Add DDX58/RIG-I manually to the table
- BLAST check in NCBI Gene: BAX, MAVS, ASC, BID, RBCK1

**On CSD3 (Steps 2-3):**
- BLAST all 67 human proteins against ~50 avian proteomes
- Foldseek screen (both arms) against avian AFDB proteomes
- Submit as Slurm job arrays, never loop on squeue

---

## 8. FILES AND LOCATIONS

- **Repository:** ~/aviced_bioinformatics (github.com/juliethirenemurillo/aviced_bioinformatics, private)
- **Scripts:** in aviced_bioinformatics/workflow/scripts/track_a/(01-06_*.py)
- **Results:** results/ directory
- **Key output:** results/orthologue_domain_summary.xlsx (clean color-coded summary, 17 columns)
- **Full data:** results/human_cell_death_gene_universe_v1_bird_orthologues_domains.xlsx (all 36 columns)
- **Reference table:** human_cell_death_gene_universe_v1.xlsx (67 genes, verified UniProt accessions)
