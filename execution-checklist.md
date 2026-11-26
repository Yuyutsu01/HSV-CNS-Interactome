# Phased Execution Checklist

Use this as a week-by-week tracking document.

---

## PHASE 1 — Dataset Construction and Disorder Mapping (Weeks 1–4)

### Week 1: Sequence Retrieval
- [ ] Download ICP34.5 (P36313), US3 (P04413), pUL37 (P10218) from UniProt
- [ ] Retrieve alphaherpesvirus orthologs for each (HSV-2, VZV, PRV, BoHV-1, EHV-1)
- [ ] Run quick IUPred3 scan on all 72 HSV-1 ORFs to identify control candidates
- [ ] Select 8–10 negative control proteins; download their sequences
- [ ] Download AF2 predicted structures from AlphaFold DB for all proteins

### Week 2: Disorder Prediction
- [ ] Run IUPred3 (long disorder + ANCHOR2) on all target and control proteins
- [ ] Run metapredict V2 on all proteins
- [ ] Extract pLDDT scores from AF2 models for all proteins
- [ ] Calculate consensus disorder (2/3 agreement) for all proteins
- [ ] Define IDR boundaries (≥ 20 contiguous consensus-disordered residues)
- [ ] Assign confidence tiers (Tier 1/2/3) to each protein's IDR map

### Week 3: Accessibility and Secondary Structure
- [ ] Run PSIPRED on all proteins (secondary structure prediction)
- [ ] Run s2D on IDR segments (disorder-specific secondary structure)
- [ ] Run NetSurfP-3.0 on all proteins (solvent accessibility)
- [ ] Run localCIDER on IDR segments (charge patterning, kappa)
- [ ] Compile per-residue annotation matrix for each protein

### Week 4: Phase 1 Quality Control and Documentation
- [ ] Generate disorder profile plots (3-method overlay) for all proteins
- [ ] Create Figure 2 draft (multi-method disorder maps)
- [ ] Write Phase 1 methods section draft
- [ ] Document all IDR boundaries in a master table
- [ ] Verify ICP34.5 IDRs match UniProt annotations (validation check)

**Phase 1 Milestone:** Complete IDR maps with confidence tiers for all target and control proteins.

---

## PHASE 2 — SLiM Discovery and Evolutionary Filtering (Weeks 5–10)

### Week 5: SLiM Scanning
- [ ] Download full ELM class definitions (regex patterns + annotations)
- [ ] Run ELM pattern scanning on all IDR segments (targets + controls)
- [ ] Run SLiMSearch on full-length proteins
- [ ] Compile raw hit list per protein
- [ ] Apply Filter 1: retain only hits within consensus IDR boundaries

### Week 6: SLiM Filtering
- [ ] Apply Filter 2: accessibility threshold (NetSurfP score > 0.4)
- [ ] Apply Filter 3: flanking sequence context check
- [ ] Build multiple sequence alignments (MAFFT L-INS-i) for each protein
- [ ] Calculate per-residue conservation scores (Shannon entropy, Valdar)
- [ ] Apply Filter 4: conservation — retain SLiMs conserved in ≥ 2 genera

### Week 7: Evolutionary Analysis
- [ ] Prepare codon alignments for PAML (PAL2NAL or manual)
- [ ] Run PAML codeml: M0, M1a vs M2a, M7 vs M8
- [ ] Calculate per-site dN/dS ratios
- [ ] Compare dN/dS in IDR segments vs. structured domains
- [ ] Identify sites under purifying selection within IDRs
- [ ] Check for sites under positive selection

### Week 8: Mimicry and ML Scoring
- [ ] Cross-reference all predicted SLiMs against human ELM instances
- [ ] Flag molecular mimicry candidates (viral SLiM matches human SLiM class)
- [ ] Retrieve known host binding partners for mimicked SLiMs
- [ ] Build ML training dataset (positive: known functional viral SLiMs; negative: random IDR peptides)
- [ ] Extract feature vectors for all SLiM hits
- [ ] Train XGBoost classifier with 5-fold stratified CV

### Week 9: ML Scoring and Ranking
- [ ] Score all predicted SLiMs with trained model
- [ ] Calculate AUC, precision-recall, and feature importances
- [ ] Generate SHAP plots for interpretability
- [ ] Rank all SLiMs per protein by ML score
- [ ] Create Figure 3 draft (SLiM landscape maps)
- [ ] Create Figure 4 draft (conservation on phylogenetic tree)

### Week 10: Phase 2 Quality Control
- [ ] Cross-validate top SLiM predictions against published mutagenesis data
- [ ] Check ICP34.5: do predicted SLiMs overlap with known Beclin-1/PP1α/TBK1 binding sites?
- [ ] Check US3: do predicted SLiMs overlap with known virulence variants (T52P, A153T)?
- [ ] Write Phase 2 results section draft
- [ ] Compile master SLiM catalog table (Supplementary Table)

**Phase 2 Milestone:** Ranked, filtered SLiM catalog with conservation, mimicry, and ML scores for all proteins.

---

## PHASE 3 — Conformational Sampling and Structural Modeling (Weeks 8–14)

*Note: Phase 3 starts in Week 8, overlapping with Phase 2.*

### Week 8–9: MD System Setup
- [ ] Install CALVADOS 2 / OpenMM environment
- [ ] Identify all IDR segments ≥ 30 residues for simulation
- [ ] Build CG systems for each IDR segment
- [ ] Run test simulations (100 ns) to verify setup
- [ ] Submit production runs (10 μs × 3 replicates per IDR)

### Week 10–11: MD Analysis
- [ ] Monitor simulation convergence (Rg time series, block averaging)
- [ ] Calculate Rg distributions for each IDR
- [ ] Generate contact maps from ensemble
- [ ] Calculate per-residue SASA from CG ensemble
- [ ] Back-map representative frames to all-atom for DSSP analysis

### Week 12: SLiM Exposure Analysis
- [ ] Calculate SLiM exposure frequency for each predicted SLiM
- [ ] Identify SLiMs with low exposure (< 20% of frames) — flag as likely non-functional
- [ ] For IDRs with multiple SLiMs: calculate joint exposure probabilities
- [ ] Test for anti-correlated exposure (competitive exclusion)
- [ ] Identify alpha-MoRFs (transient helical segments in IDRs)

### Week 13: Peptide Docking (Top Candidates)
- [ ] Select top 5 SLiM–host domain pairs by ML score
- [ ] Search PDB for co-crystal structures of the same SLiM class + domain
- [ ] For template-available cases: run FlexPepDock refinement
- [ ] For no-template cases: run HADDOCK peptide docking
- [ ] Calculate predicted binding energies
- [ ] Compare viral vs. human SLiM binding scores (mimicry test)

### Week 14: Phase 3 Integration
- [ ] Create Figure 5 draft (conformational ensembles with SLiM coloring)
- [ ] Update SLiM catalog with exposure frequencies and docking scores
- [ ] Write Phase 3 methods and results section draft

**Phase 3 Milestone:** Conformational ensemble properties, SLiM exposure frequencies, and top-candidate docking models.

---

## PHASE 4 — CNS-Weighted Interactome and Network Analysis (Weeks 12–18)

### Week 12–13: Network Construction
- [ ] Download HintHSV dataset
- [ ] Download IntAct and BioGRID HSV-1 interactions
- [ ] Merge databases, remove duplicates, assign confidence scores
- [ ] Build bipartite virus–host network
- [ ] Add SLiM-predicted edges (from Phase 2 catalog)
- [ ] Annotate each edge with evidence type (experimental PPI, SLiM prediction, both)

### Week 14–15: CNS Tissue Weighting
- [ ] Download GTEx v8 median TPM expression matrix
- [ ] Calculate tissue specificity (tau) for all human proteins in network
- [ ] Score each virus–host edge by host protein's CNS specificity
- [ ] Generate CNS-weighted network
- [ ] Identify top CNS-enriched interaction partners per viral protein

### Week 16: Pathway Analysis
- [ ] Run GO enrichment (clusterProfiler) per viral protein's interaction set
- [ ] Run KEGG/Reactome pathway enrichment
- [ ] Focus analysis: autophagy, innate immunity, apoptosis, axonal transport
- [ ] Run HotNet2 or NeST network perturbation analysis
- [ ] Identify significantly perturbed subnetworks

### Week 17: The Core Hypothesis Test
- [ ] Calculate SLiM density per IDR residue (all targets + controls)
- [ ] Calculate CNS interaction fraction (all targets + controls)
- [ ] Compute composite CNS-targeting score per protein
- [ ] Run Mann-Whitney U test: neurovirulence vs. controls
- [ ] Calculate effect size (rank-biserial r)
- [ ] Run permutation test (10,000 permutations) as robustness check
- [ ] Create Figure 7 (enrichment boxplot — the key figure)

### Week 18: Phase 4 Integration
- [ ] Create Figure 6 draft (CNS-weighted network, Cytoscape)
- [ ] Write Phase 4 methods and results section draft
- [ ] Compile pathway enrichment tables (Supplementary)

**Phase 4 Milestone:** CNS-weighted network, pathway enrichment, and statistical enrichment test result.

---

## PHASE 5 — Competitive-Binding Switch Identification and Validation (Weeks 16–22)

### Week 16–17: Switch Discovery
- [ ] For each IDR with ≥ 2 high-confidence SLiMs, run the 4-test decision tree
- [ ] Test 1: Distance check (≤ 30 residues apart)
- [ ] Test 2: Domain overlap check (from 3did/PDB)
- [ ] Test 3: Conformational exclusion (from MD data)
- [ ] Test 4: PTM toggle (NetPhos/GPS-Prot predictions between SLiMs)
- [ ] Classify each candidate: Strong / Moderate / Weak switch

### Week 18–19: Switch Annotation
- [ ] Annotate each switch with competing SLiMs and host targets
- [ ] Annotate toggle mechanism
- [ ] Annotate biological consequences of each switch state
- [ ] Compare to ICP34.5 Beclin-1/PP1α dual-SLiM architecture and E1A tethered-SLiM paradigm
- [ ] Create Figure 8 draft (switch models for top candidates)

### Week 20–21: Clinical Variant Validation
- [ ] Retrieve HSV-1 clinical genome sequences (encephalitis + non-neurological)
- [ ] Map nonsynonymous variants onto SLiM predictions
- [ ] Fisher's exact test: variant clustering in SLiMs (encephalitis vs. other)
- [ ] Specific check: US3 T52P and A153T vs. predicted SLiMs
- [ ] Literature cross-validation for all high-confidence predictions

### Week 22: Phase 5 Integration
- [ ] Compile switch catalog with all annotations
- [ ] Write Phase 5 results and discussion draft

**Phase 5 Milestone:** Annotated molecular switch catalog with clinical and literature validation.

---

## PHASE 6 — Integration, Visualization, and Manuscript (Weeks 20–26)

### Week 20–22: Figures and Tables
- [ ] Finalize all 8 main figures
- [ ] Create all supplementary figures
- [ ] Build supplementary tables (SLiM catalog, network edges, enrichment results)
- [ ] Create graphical abstract (Fig. 1)

### Week 23–24: Writing
- [ ] Write Introduction (1,500 words)
- [ ] Integrate all Results section drafts (4,000 words total)
- [ ] Write Discussion: biological implications, limitations, experimental roadmap (2,000 words)
- [ ] Write comprehensive Methods section (3,000 words)
- [ ] Write Abstract (250 words, structured)
- [ ] Compile references (target: 60–80 citations)

### Week 25: Internal Review
- [ ] Complete self-review pass for scientific accuracy
- [ ] Check all statistics for correctness
- [ ] Verify all figure legends match figure content
- [ ] Proofread and format for target journal (PLoS Comp Biol or Virus Evolution)
- [ ] Prepare cover letter

### Week 26: Submission
- [ ] Final formatting check against journal guidelines
- [ ] Upload to preprint server (bioRxiv) if desired
- [ ] Submit to primary target journal
- [ ] Archive all code and data to GitHub/Zenodo with DOI

**Phase 6 Milestone:** Submitted manuscript.

---

## Critical Decision Points

| Week | Decision | Impact |
|------|----------|--------|
| 2 | If pUL37 shows no consensus IDR overlap with effector region | Demote pUL37 to supplementary analysis only |
| 5 | If raw SLiM hits are > 500 per protein | Tighten filtering thresholds or add additional filter step |
| 10 | If ML model AUC < 0.7 | Revise feature set or switch to random forest; consider reporting without ML scoring |
| 14 | If MD shows most SLiMs have < 10% exposure | Revise accessibility thresholds; discuss as negative result |
| 17 | If enrichment test p > 0.05 | Report honestly; reframe paper around individual protein findings rather than enrichment claim |
| 19 | If no strong switch candidates found | Drop H2; focus paper on H1 (enrichment) and H3 (mimicry) |

---

## Data Management Plan

| Data Type | Storage | Sharing |
|-----------|---------|---------|
| Sequences and alignments | GitHub repository | Public at submission |
| Disorder predictions | GitHub repository | Public at submission |
| SLiM catalogs | GitHub + Supplementary Tables | Public at submission |
| MD trajectories | Local storage + institutional HPC | Available on request (large files) |
| Network files | GitHub (Cytoscape sessions) | Public at submission |
| ML models | GitHub (trained model + training data) | Public at submission |
| Analysis scripts | GitHub repository | Public at submission, MIT license |
