# Research Proposal

## IDR-Encoded Molecular Switches in HSV-1 Neurovirulence: A Computational Discovery of SLiM-Mediated Host Rewiring by ICP34.5, US3, and pUL37

---

## 1. Abstract

Herpes simplex virus type 1 (HSV-1) establishes lifelong latency in sensory neurons and can reactivate to cause severe neurological disease, including encephalitis. Three viral proteins — ICP34.5, US3 kinase, and pUL37 tegument protein — are established neurovirulence determinants, yet the molecular logic by which they simultaneously engage multiple host defense pathways in the central nervous system (CNS) remains poorly understood. Intrinsically disordered regions (IDRs) in viral proteins can encode short linear motifs (SLiMs) that mimic host signaling peptides, enabling compact multivalent hijacking of cellular machinery. ICP34.5 itself exemplifies this paradigm: its disordered N-terminus binds Beclin-1 to block autophagy, while its C-terminal PP1α-binding motif (an RVxF-class SLiM) drives eIF2α dephosphorylation to restore translation — two SLiM-mediated functions housed in a single predominantly disordered protein, separated by a variable-length repeat linker. The adenovirus E1A protein, where two tethered SLiMs in a disordered region achieve picomolar Rb binding through cooperative displacement, further validates the competitive-binding switch concept in viral IDRs (Glavina et al., 2022). We propose a systematic computational investigation of IDR-embedded SLiMs in ICP34.5, US3, and pUL37, integrating multi-method disorder prediction, cross-alphaherpesvirus evolutionary analysis, conformational sampling, and CNS-weighted interactome mapping. A matched negative-control analysis of non-neurovirulence HSV-1 proteins will test the central hypothesis: that neurovirulence determinants are specifically enriched for CNS-targeting SLiMs arranged in competitive-binding architectures. The outcome will be a ranked, experimentally testable catalog of predicted molecular switches and a generalizable framework linking IDR architecture to viral tissue tropism.

---

## 2. Background and Rationale

### 2.1 HSV-1 Neurovirulence: An Unsolved Multitarget Problem

HSV-1 infects over 3.7 billion people globally. While most infections are subclinical, the virus is the leading cause of sporadic viral encephalitis, with mortality exceeding 70% if untreated. Neurovirulence — the capacity to replicate in and damage neural tissue — depends on a small set of viral gene products that counteract host innate immunity, autophagy, and apoptosis specifically within the CNS. ICP34.5, US3, and pUL37 have each been shown through deletion and mutation studies to be required for full neurovirulence in animal models, yet each engages multiple host pathways simultaneously. How a single viral protein can multitask across several host defense systems is a central unsolved question.

### 2.2 The IDR–SLiM Paradigm in Viral Pathogenesis

Intrinsically disordered regions lack stable three-dimensional structure under physiological conditions but are far from functionless. IDRs serve as interaction platforms by encoding short linear motifs — typically 3–10 residue peptide segments that bind structured domains in host proteins. Because SLiMs are compact and degenerate, a single IDR can pack multiple, overlapping interaction interfaces into a short sequence, enabling:

- **Multivalency**: simultaneous or sequential engagement of different host targets
- **Competitive switching**: two SLiMs in the same IDR compete for binding, creating a binary regulatory switch controlled by local concentration or post-translational modification
- **Molecular mimicry**: viral SLiMs that match host SLiM consensus sequences can competitively displace endogenous interactions

The SLiM paradigm has been established for several viral families (e.g., papillomavirus E7, adenovirus E1A, Nipah V protein), but systematic analysis in herpesviruses — particularly for neurovirulence — is lacking.

### 2.3 The ICP34.5 Dual-SLiM Architecture and the Adenovirus E1A Precedent

ICP34.5 itself provides the clearest mechanistic template within HSV-1. The protein is predominantly disordered: its N-terminal domain (residues ~1–147) lacks stable structure and contains the Beclin-1 binding domain (BBD, residues 68–87), through which ICP34.5 inhibits autophagy in neurons — a function essential for neurovirulence in vivo (Orvedahl et al., 2007). The C-terminal domain harbors a canonical RVxF-class PP1α-binding SLiM that recruits protein phosphatase 1α to dephosphorylate eIF2α, restoring host translation shut down by PKR (He et al., 1997). These two SLiM-mediated functions — autophagy inhibition and translational rescue — are separated by a variable-length linker containing Ala-Thr-Pro (ATP) repeats, a hallmark of low-complexity disordered sequence. Whether additional SLiMs reside within ICP34.5's IDRs, and whether any of these motifs compete for overlapping host partners or are regulated by post-translational modifications, has never been systematically investigated.

Beyond HSV-1, the adenovirus E1A oncoprotein demonstrates that competitive SLiM-based switching within a viral IDR is a proven mechanism: two SLiMs tethered by a disordered linker cooperatively bind Rb with picomolar affinity, displacing host E2F transcription factors (Glavina et al., 2022). This architecture — multiple SLiMs within a single IDR, engaging host targets through competitive or cooperative binding — is exactly what we propose to search for systematically in ICP34.5, US3, and pUL37.

### 2.4 Knowledge Gap

Despite strong individual evidence for IDRs in ICP34.5 (UniProt-annotated disorder spanning residues 1–147, 162–191, 252–295), US3 (predicted disordered N- and C-terminal tails flanking the kinase domain), and pUL37 (structurally characterized disordered C-terminal tail overlapping the neuroinvasion effector region), no study has:

1. Systematically mapped SLiMs within these IDRs using modern computational tools
2. Assessed whether these SLiMs are evolutionarily conserved across alphaherpesviruses
3. Tested whether neurovirulence proteins are specifically enriched for CNS-relevant SLiMs compared to other IDR-containing HSV-1 proteins
4. Modeled competitive-binding architectures that could explain multitarget host manipulation

This proposal addresses all four gaps.

---

## 3. Hypotheses

**Primary Hypothesis (H1):** The IDRs of ICP34.5, US3, and pUL37 encode functionally constrained SLiMs that mediate interactions with CNS-expressed host proteins, and these SLiMs are enriched relative to IDRs in non-neurovirulence HSV-1 proteins.

**Secondary Hypothesis (H2):** A subset of these SLiMs are arranged in competitive-binding architectures (overlapping or adjacent motifs targeting mutually exclusive host partners), forming molecular switches analogous to the ICP34.5 Beclin-1/PP1α dual-function architecture and the adenovirus E1A cooperative-binding paradigm.

**Tertiary Hypothesis (H3):** Viral SLiMs in neurovirulence IDRs mimic endogenous human SLiMs, enabling competitive displacement of host-host interactions in autophagy, innate immunity, and axonal transport pathways.

---

## 4. Research Objectives

| # | Objective | Deliverable |
|---|-----------|-------------|
| O1 | Establish consensus IDR maps for ICP34.5, US3, and pUL37 using multi-method prediction | Per-protein disorder profile with confidence tiers |
| O2 | Discover and filter SLiMs within mapped IDRs | Ranked SLiM catalog with conservation, accessibility, and mimicry scores |
| O3 | Build a CNS-weighted host–virus interaction network | Tissue-specific interaction map with pathway enrichment |
| O4 | Test enrichment of CNS-targeting SLiMs in neurovirulence vs. control proteins | Statistical comparison (the core hypothesis test) |
| O5 | Identify competitive-binding switch architectures | Predicted molecular switches with structural and conformational evidence |
| O6 | Validate predictions against clinical and experimental data | Cross-referencing with encephalitis isolate variants and published mutagenesis |

---

## 5. Methodology — Phased Execution

### ═══════════════════════════════════════════
### PHASE 1: Dataset Construction and Disorder Mapping
### (Weeks 1–4)
### ═══════════════════════════════════════════

#### 5.1.1 Sequence Retrieval and Curation

**Target proteins:**
- ICP34.5 (gene RL1): UniProt P36313 (strain 17), plus orthologs from HSV-2, VZV (ORF47 analog), PRV (IE180 context), BoHV-1, EHV-1
- US3: UniProt P04413 (strain 17), plus alphaherpesvirus orthologs
- pUL37: UniProt P10218 (strain 17), plus alphaherpesvirus orthologs

**Negative control set** (critical for H1 testing):
- Select 8–10 HSV-1 proteins that: (a) contain predicted IDRs of comparable length, (b) are NOT classified as neurovirulence determinants, and (c) span functional categories (DNA replication, capsid, envelope)
- Candidates: UL29/ICP8 (ssDNA-binding, large IDRs), UL42 (processivity factor), gD/US6 (receptor binding, flexible regions), UL46/VP11/12 (tegument), UL48/VP16 (transactivator with known IDR)

**Data sources:** UniProt, NCBI RefSeq, AlphaFold Protein Structure Database, PDB (for any experimentally solved structures or fragments).

#### 5.1.2 Multi-Method Disorder Prediction

Run three independent disorder predictors on all target and control proteins:

| Method | Tool | Rationale |
|--------|------|-----------|
| Physics-based energy estimation | IUPred3 (long disorder + ANCHOR2 for binding regions) | Gold standard; detects both constitutive disorder and conditional disorder |
| Machine-learning consensus | metapredict V2 | Fast, accurate, trained on disprot, avoids IUPred's known biases |
| Structure-based proxy | AlphaFold2 pLDDT scores (pLDDT < 50 = disordered) | Independent of sequence-based predictors; structural validation |

**Consensus rule:** A residue is classified as disordered if ≥ 2 of 3 methods agree. IDR boundaries are defined as contiguous stretches of ≥ 20 consensus-disordered residues (standard threshold for functional IDRs).

#### 5.1.3 Confidence Tier Assignment

Each protein's IDR map receives a confidence tier:

- **Tier 1 (High):** UniProt-annotated disorder + multi-method consensus + experimental support (NMR, SAXS, or hydrogen-deuterium exchange). ICP34.5 qualifies.
- **Tier 2 (Medium-High):** Multi-method consensus + indirect experimental support (mutagenesis mapping functional residues to predicted IDR). US3 qualifies.
- **Tier 3 (Medium):** Multi-method consensus + structural inference (AF2 model, limited experimental data). pUL37 qualifies.

This tiering is reported transparently in the paper.

#### 5.1.4 Secondary Structure and Accessibility Profiling

For each IDR, calculate:
- **Secondary structure propensity:** PSIPRED and s2D (disorder-specific secondary structure predictor) to identify alpha-MoRFs (molecular recognition features that transiently fold upon binding)
- **Solvent accessibility:** NetSurfP-3.0 predictions on the full-length sequence
- **Charge distribution and sequence composition:** CIDER/localCIDER for charge patterning, proline content, and sequence complexity — these physicochemical features influence whether SLiMs within the IDR are accessible or buried in transient compact states

**Deliverable:** A per-residue annotation matrix for each protein combining disorder score, pLDDT, secondary structure propensity, accessibility, and charge pattern.

---

### ═══════════════════════════════════════════
### PHASE 2: SLiM Discovery and Evolutionary Filtering
### (Weeks 5–10)
### ═══════════════════════════════════════════

#### 5.2.1 Pattern-Based SLiM Scanning

Scan all IDRs against established motif databases:

| Resource | What It Provides |
|----------|-----------------|
| ELM (Eukaryotic Linear Motif) database | Curated regex patterns for ~320 SLiM classes with biological context |
| SLiMSearch / SLiMFinder | Proteome-wide SLiM annotation and de novo motif discovery |
| VirusMINT / PHISTO | Virus-specific interaction motifs |
| 3did (Domain-Domain Interaction) | Domain-motif interaction templates |

**Filtering cascade:**
1. Raw hits from ELM regex scanning (expect high false-positive rate)
2. Filter 1: Retain only hits within consensus IDR boundaries
3. Filter 2: Retain only hits with solvent accessibility score > 0.4 (exposed)
4. Filter 3: Retain only hits where flanking sequence context matches known SLiM grammar (e.g., correct charge/hydrophobicity flanking the core motif)
5. Filter 4: Conservation filter (see 5.2.2)

#### 5.2.2 Cross-Alphaherpesvirus Evolutionary Analysis

**Alignment construction:**
- For each target protein, build a multiple sequence alignment (MSA) using MAFFT (L-INS-i mode for accuracy) across:
  - HSV-1 strains (17, KOS, F, McKrae — representing lab and clinical isolates)
  - HSV-2 (HG52)
  - VZV (Dumas, Oka)
  - PRV (Becker, Kaplan)
  - BoHV-1 (Cooper)
  - EHV-1 (Ab4)

**Conservation metrics:**
- **Per-residue conservation:** Shannon entropy and Valdar conservation score across the MSA
- **SLiM-level conservation:** For each predicted SLiM, assess whether the core motif regex is satisfied in orthologs (allowing degenerate positions). A SLiM conserved across ≥ 2 alphaherpesvirus genera is flagged as high-confidence.
- **dN/dS analysis:** Using PAML (codeml, site models M1a vs. M2a and M7 vs. M8), calculate per-site dN/dS ratios across the alignment. Critically, compare dN/dS within IDR segments vs. structured domains. Purifying selection (dN/dS < 1) on IDR residues is strong evidence of functional constraint, since IDRs are otherwise expected to evolve rapidly. Positive selection on specific IDR sites may indicate adaptive SLiM diversification.
- **Convergent SLiM acquisition:** Screen for cases where the same ELM-class SLiM appears in orthologous IDR positions across distantly related alphaherpesviruses but in non-homologous sequence context — evidence of convergent evolution toward the same host-interaction strategy.

#### 5.2.3 Molecular Mimicry Analysis

For each high-confidence SLiM prediction:
1. Query the human ELM instance database: does any human protein contain the same SLiM class in a functionally characterized context?
2. If yes, retrieve the human SLiM and its binding partner
3. Assess whether the viral SLiM could competitively displace the human SLiM based on: (a) motif similarity, (b) flanking sequence compatibility, (c) predicted binding affinity (using structure-based scoring if a co-crystal structure exists for the human SLiM–domain interaction)

**Priority targets for mimicry:**
- LC3-interacting regions (LIR motifs) — ICP34.5 is known to target autophagy via Beclin-1; does it also mimic LIR motifs?
- TRAF-interacting motifs (TIMs) — innate immune signaling
- Nuclear export/import signals — relevant for US3 nuclear functions
- Dynein/kinesin-binding motifs — relevant for pUL37 axonal transport

#### 5.2.4 Machine Learning SLiM Scoring

Train a gradient-boosting classifier (XGBoost) to score each predicted SLiM:

**Training data:** Known functional viral SLiMs from ELM (positive set, ~200 instances) vs. random IDR peptides matching the same regex in non-interacting contexts (negative set, 5:1 ratio).

**Features (per SLiM instance):**
- Disorder score (IUPred3 average over motif)
- pLDDT score (AF2 average over motif)
- Conservation score (Valdar, across alphaherpesvirus MSA)
- Solvent accessibility (NetSurfP average)
- Flanking charge asymmetry (localCIDER kappa)
- Flanking disorder extent (length of contiguous IDR around motif)
- Secondary structure propensity (helix/strand/coil from PSIPRED)
- Sequence complexity (Shannon entropy of motif + 5-residue flanks)

**Output:** A probability score (0–1) for each SLiM instance, enabling ranking rather than binary classification. Report AUC, precision-recall, and feature importances.

**Deliverable:** A ranked SLiM catalog per protein, with conservation tier, mimicry flag, ML score, and accessibility annotation.

---

### ═══════════════════════════════════════════
### PHASE 3: Conformational Sampling and Structural Modeling
### (Weeks 8–14, overlapping with Phase 2)
### ═══════════════════════════════════════════

#### 5.3.1 Coarse-Grained Molecular Dynamics of IDRs

For each IDR segment (≥ 30 residues), run coarse-grained MD simulations using the CALVADOS 2 force field (specifically parameterized for intrinsically disordered proteins):

- **System setup:** Single-chain IDR in implicit solvent, 150 mM NaCl
- **Simulation length:** 10 μs effective time per IDR (achievable in ~24–48 hours on GPU for CG models)
- **Replicates:** 3 independent runs per IDR, different random seeds
- **Analysis:**
  - Radius of gyration (Rg) distribution — compact vs. extended ensemble
  - Contact maps — which residue pairs transiently interact?
  - Per-residue solvent-accessible surface area (SASA) from the ensemble
  - SLiM exposure frequency — for each predicted SLiM, what fraction of the ensemble has the motif fully solvent-exposed?
  - Transient secondary structure — DSSP analysis on CG-to-AA backmapped snapshots to identify alpha-MoRFs

#### 5.3.2 Competitive Accessibility Analysis

For IDRs containing multiple predicted SLiMs:
- Calculate the joint exposure probability: P(SLiM_A exposed AND SLiM_B exposed) vs. P(SLiM_A exposed) × P(SLiM_B exposed)
- If joint exposure is significantly lower than the product of marginals, the SLiMs are in a competitive-exclusion regime — the IDR cannot simultaneously present both motifs, creating a conformational switch
- Identify post-translational modification sites (predicted by NetPhos, GPS-Prot) between competing SLiMs that could toggle the conformational equilibrium

#### 5.3.3 Peptide–Domain Docking (Top Candidates Only)

For the top 5 highest-scoring SLiM–host domain pairs:
- Extract the SLiM peptide (motif + 5 flanking residues on each side)
- If a co-crystal structure exists for the same SLiM class bound to the host domain (from PDB), use it as a template for homology-based peptide docking (FlexPepDock refinement in Rosetta)
- If no template exists, use HADDOCK or ClusPro peptide-docking protocols
- Calculate predicted binding affinity (Rosetta interface energy, FoldX ΔΔG)
- Compare viral SLiM binding scores to the endogenous human SLiM for the same domain — if the viral SLiM scores comparably or better, the mimicry/displacement hypothesis is supported

**Deliverable:** Conformational ensemble properties per IDR, SLiM exposure frequencies, competitive-exclusion flags, and top-candidate docking models.

---

### ═══════════════════════════════════════════
### PHASE 4: CNS-Weighted Interactome and Network Analysis
### (Weeks 12–18, overlapping with Phase 3)
### ═══════════════════════════════════════════

#### 5.4.1 Host–Virus PPI Network Construction

**Data sources:**
- HintHSV database (10,432 human–HSV-1 interactions, 4,546 human proteins, 72 HSV-1 proteins)
- IntAct and BioGRID for additional experimentally validated interactions
- VirusMentha for curated virus–host PPIs
- NeuroViOme ORFeome for neurotropic virus interactome context

**Network construction:**
- Build a bipartite network: HSV-1 proteins ↔ human proteins
- For each predicted SLiM, add an edge from the viral protein to the predicted host binding partner (the domain-containing human protein)
- Edge confidence scoring: combine evidence from (a) SLiM ML score, (b) existing PPI evidence from databases, (c) domain-motif interaction evidence from 3did/ELM

#### 5.4.2 CNS Tissue Weighting

**Expression data:** GTEx v8 bulk RNA-seq, focusing on:
- Brain regions: cortex, hippocampus, cerebellum, substantia nigra, spinal cord (cervical C1)
- Peripheral tissues for comparison: liver, lung, skin, whole blood

**Weighting scheme:**
- For each human protein in the network, calculate a CNS specificity score: τ (tissue specificity index) using the method of Yanai et al., restricted to brain vs. non-brain tissues
- Weight each virus–host edge by the host protein's CNS specificity score
- A CNS-weighted network retains all edges but highlights interactions with CNS-enriched host partners

#### 5.4.3 Pathway Enrichment and Perturbation Analysis

**Pathway enrichment:**
- For each viral protein's predicted SLiM-mediated interactors, run Gene Ontology (GO) and KEGG/Reactome pathway enrichment (clusterProfiler, with BH-corrected p-values)
- Focus pathways: autophagy (KEGG hsa04140), innate immune signaling (RIG-I, cGAS-STING, NF-κB), apoptosis (hsa04210), axonal transport (GO:0008088), neurotrophin signaling (hsa04722)

**Network perturbation scoring:**
- Use HotNet2 or NeST to identify significantly perturbed subnetworks when the SLiM-mediated interactions are added to the basal human PPI network
- This tests whether the predicted viral interactions converge on specific network modules rather than scattering randomly across the interactome

#### 5.4.4 The Core Hypothesis Test: Neurovirulence Enrichment

**Statistical framework for H1:**

For each protein (3 neurovirulence targets + 8–10 negative controls):
1. Count the number of high-confidence SLiMs (ML score > 0.7) per IDR residue (SLiM density)
2. Calculate the fraction of SLiM-mediated interactions targeting CNS-enriched host proteins (CNS interaction fraction)
3. Compute a composite "CNS-targeting score" = SLiM density × CNS interaction fraction

**Test:** Mann-Whitney U test (or permutation test if sample size is limiting) comparing CNS-targeting scores of neurovirulence proteins vs. negative controls.

**Expected outcome:** If H1 is supported, neurovirulence proteins will show significantly higher CNS-targeting scores. If not, the paper honestly reports this and discusses alternative explanations (e.g., CNS tropism may depend on protein localization or expression timing rather than SLiM content).

**Deliverable:** CNS-weighted interaction network, pathway enrichment results, and statistical test of the enrichment hypothesis.

---

### ═══════════════════════════════════════════
### PHASE 5: Competitive-Binding Switch Identification and Clinical Validation
### (Weeks 16–22)
### ═══════════════════════════════════════════

#### 5.5.1 Switch Architecture Discovery

For each IDR containing ≥ 2 high-confidence SLiMs:

1. **Spatial proximity test:** Are the SLiMs within 30 residues of each other in primary sequence? (Threshold informed by the E1A tethered-SLiM architecture, where cooperative motifs are connected by a short disordered linker)
2. **Domain compatibility test:** Do the predicted host-binding domains for each SLiM have overlapping binding surfaces? (Check from co-crystal structures or domain family binding-mode annotations in 3did)
3. **Conformational exclusion test:** From Phase 3 MD data, do the SLiMs show anti-correlated exposure? (One exposed → other buried)
4. **PTM toggle test:** Is there a predicted phosphorylation, ubiquitination, or acetylation site between or within the competing SLiMs that could switch binding preference?

A candidate switch must pass tests 1 + 2 + at least one of (3, 4).

#### 5.5.2 Switch Annotation

Each predicted switch is annotated with:
- The two (or more) competing SLiMs and their host targets
- The predicted toggle mechanism (PTM site, conformational change, or concentration-dependent competition)
- The biological consequence of each switch state (e.g., State A: autophagy inhibition via Beclin-1 binding; State B: immune evasion via TBK1 sequestration)
- Confidence score based on the number of tests passed

#### 5.5.3 Clinical Variant Validation

**Data source:** Published HSV-1 genome sequences from:
- Encephalitis isolates (NCBI BioProject collections of clinical HSV-1 genomes)
- Keratitis isolates
- Genital/labial isolates (non-neurological disease)

**Analysis:**
- Map nonsynonymous variants from clinical isolates onto the SLiM predictions
- Test whether variants in encephalitis isolates cluster in or near predicted SLiMs more than variants in non-neurological isolates (Fisher's exact test)
- For US3 specifically: the known virulence-associated variants T52P and A153T — do these map to predicted SLiMs or switch regions?

#### 5.5.4 Literature Cross-Validation

Systematically cross-reference all high-confidence predictions against published mutagenesis studies:
- ICP34.5: Beclin-1 binding domain (residues 68–87 in strain 17), PP1α binding motif, TBK1 interaction region — do these correspond to predicted SLiMs?
- US3: substrate-binding residues, virulence-mapping mutations
- pUL37: neuroinvasion effector region, depalmitoylation site

Concordance between predicted SLiMs and experimentally mapped functional regions is the strongest validation available without new wet-lab data.

**Deliverable:** Annotated switch catalog, clinical variant overlay, and literature concordance analysis.

---

### ═══════════════════════════════════════════
### PHASE 6: Integration, Visualization, and Manuscript Preparation
### (Weeks 20–26)
### ═══════════════════════════════════════════

#### 5.6.1 Figure Design

| Figure | Content | Purpose |
|--------|---------|---------|
| Fig. 1 | Conceptual model: IDR–SLiM switches in neurovirulence (schematic) | Graphical abstract / introduction figure |
| Fig. 2 | Multi-method disorder maps for ICP34.5, US3, pUL37 with confidence tiers | Phase 1 results |
| Fig. 3 | SLiM landscape: per-protein motif maps with conservation heatmap overlay | Phase 2 results |
| Fig. 4 | Cross-alphaherpesvirus SLiM conservation mapped onto phylogenetic tree | Evolutionary evidence |
| Fig. 5 | Conformational ensemble snapshots with SLiM exposure coloring | Phase 3 results |
| Fig. 6 | CNS-weighted interaction network (Cytoscape, with pathway module coloring) | Phase 4 results |
| Fig. 7 | Neurovirulence enrichment: CNS-targeting score boxplot, targets vs. controls | The H1 test — the key figure |
| Fig. 8 | Competitive-binding switch models for top candidates | Phase 5 results — the payoff figure |

**Supplementary figures:** Full disorder prediction comparisons, ML model performance (ROC/PR curves), docking models, clinical variant maps.

#### 5.6.2 Manuscript Structure

**Target journal:** PLoS Computational Biology (primary) or Virus Evolution (secondary)

| Section | Content | Approximate Length |
|---------|---------|-------------------|
| Title | As above | — |
| Abstract | Structured: Background, Methods, Results, Conclusions | 250 words |
| Introduction | IDR–SLiM paradigm → HSV-1 neurovirulence → ICP34.5 dual-SLiM architecture and E1A precedent → knowledge gap → this study | 1,500 words |
| Results | Follows Phases 1–5 as results sections; leads with the enrichment test (Fig. 7) before detailing individual protein findings | 4,000 words |
| Discussion | Biological implications, comparison to other viral IDR studies, limitations (prediction-only, no wet-lab validation), experimental roadmap | 2,000 words |
| Methods | Detailed computational protocol (reproducibility-focused) | 3,000 words |
| Supplementary | Extended data tables, full MSAs, additional network analyses, ML training details | As needed |

#### 5.6.3 Proposed Experimental Validation Section (In Paper Discussion)

The paper is computational, but credibility requires proposing specific, minimal experiments:

1. **Fluorescence polarization binding assays:** Synthesize the top 3 predicted SLiM peptides, measure binding to recombinant host domains (e.g., Beclin-1 BARA domain, 14-3-3 proteins, dynein light chain)
2. **Motif mutation in viral context:** Introduce alanine substitutions at core SLiM residues in HSV-1 BAC constructs; assess replication and spread in primary neuronal cultures
3. **Competition assays:** If a molecular switch is predicted (SLiM_A vs. SLiM_B), add increasing concentrations of SLiM_B peptide to a SLiM_A–domain binding assay and measure displacement
4. **Clinical isolate phenotyping:** For any clinical variant mapping to a predicted SLiM, test the variant's effect on the SLiM–domain interaction by binding assay

---

## 6. Tools and Software

| Category | Tools |
|----------|-------|
| Disorder prediction | IUPred3, metapredict V2, AlphaFold2 (ColabFold) |
| SLiM discovery | ELM server, SLiMSearch 2.0, SLiMFinder |
| Sequence analysis | MAFFT, PAML (codeml), MEGA, Jalview |
| Structure prediction | AlphaFold2/ColabFold, PSIPRED, NetSurfP-3.0 |
| MD simulation | CALVADOS 2 (OpenMM), MDAnalysis, MDTraj |
| Peptide docking | Rosetta FlexPepDock, HADDOCK 2.4 |
| PTM prediction | NetPhos 3.1, GPS-Prot, PhosphoSitePlus |
| Network analysis | Cytoscape 3.10, HotNet2, clusterProfiler |
| ML scoring | XGBoost, scikit-learn, SHAP (for interpretability) |
| Expression data | GTEx Portal (v8), BioMart |
| PPI databases | HintHSV, IntAct, BioGRID, VirusMentha |
| Visualization | PyMOL, ChimeraX, matplotlib, seaborn |
| IDP analysis | CIDER/localCIDER, s2D |

---

## 7. Timeline

```
Week:  1----4----8----10----12----14----16----18----20----22----24----26
       |=========|                                                      Phase 1: Dataset + Disorder
              |==============|                                          Phase 2: SLiM Discovery
                    |===============|                                    Phase 3: MD + Structural
                              |===============|                          Phase 4: Network Analysis
                                        |===============|                Phase 5: Switches + Validation
                                                  |===============|      Phase 6: Writing + Figures
```

Phases overlap deliberately. Phase 3 (MD) is computationally intensive and starts as soon as Phase 1 IDR boundaries are established. Phase 4 (network) begins as Phase 2 SLiM results become available.

---

## 8. Expected Outcomes and Significance

### 8.1 Direct Outputs
1. A ranked, scored catalog of predicted SLiMs in three HSV-1 neurovirulence determinants
2. Statistical evidence for (or against) CNS-specific SLiM enrichment in neurovirulence proteins
3. Predicted molecular switches with conformational and evolutionary support
4. A CNS-weighted virus–host interaction network centered on SLiM-mediated contacts
5. Specific, testable experimental predictions

### 8.2 Broader Impact
- **Generalizable framework:** The pipeline (multi-method disorder mapping → SLiM discovery with ML scoring → tissue-weighted network → enrichment testing against controls) can be applied to any viral system where tissue tropism is linked to specific gene products
- **Therapeutic implications:** If SLiM-mediated interactions are confirmed as drivers of neurovirulence, short peptide mimetics or stapled peptides targeting these motifs become candidate antivirals
- **Evolutionary biology:** Cross-alphaherpesvirus SLiM conservation data will inform understanding of how viruses co-evolve with host signaling networks

### 8.3 Limitations (Acknowledged Proactively)
- All predictions are computational; no new experimental data are generated
- pUL37 IDR boundaries are inference-based (Tier 3 confidence)
- Coarse-grained MD captures ensemble properties but not atomic-resolution binding
- The negative control set size (8–10 proteins) limits statistical power; effect sizes must be large to reach significance
- SLiM prediction inherently has high false-positive rates; the multi-filter pipeline mitigates but does not eliminate this

---

## 9. Novelty Statement

This work is novel in four respects:

1. **First systematic SLiM analysis of HSV-1 neurovirulence determinants** — individual proteins have been studied for specific interactions, but no study has applied a unified SLiM-discovery pipeline across all three
2. **First negative-control enrichment test** — no prior study of viral IDR–SLiMs has included a matched control analysis to distinguish "viral IDRs have SLiMs" from "neurovirulence IDRs are specifically enriched for CNS-targeting SLiMs"
3. **First integration of conformational sampling with viral SLiM prediction** — prior viral SLiM studies are sequence-only; MD-based accessibility filtering adds a structural dimension
4. **First cross-alphaherpesvirus SLiM conservation analysis** — extends beyond single-species strain comparison to test evolutionary constraint across the subfamily

---

## 10. References (Key Citations)

> Note: Full reference list to be compiled during manuscript preparation. Key works to cite:

- Davison AJ (2010). Herpesvirus systematics. *Vet Microbiol*. — Alphaherpesvirus phylogeny
- He et al. (1997). The γ₁34.5 protein of herpes simplex virus 1 complexes with protein phosphatase 1α. *PNAS*. — PP1α-binding SLiM in ICP34.5
- Orvedahl et al. (2007). HSV-1 ICP34.5 confers neurovirulence by targeting the Beclin 1 autophagy protein. *Cell Host Microbe*. — Beclin-1 binding domain in ICP34.5; neurovirulence link
- Glavina et al. (2022). How viral proteins bind short linear motifs and hijack host protein function. *Essays Biochem*. — E1A tethered-SLiM competitive binding paradigm
- Chou et al. (1990). ICP34.5 neurovirulence function. — Foundational neurovirulence studies
- Mozzi et al. (2020). Intrinsically disordered regions are abundant in simplexvirus proteomes. *Virus Evol*. — IDR enrichment and positive selection in simplexviruses
- Tompa et al. (2014). A million peptide motifs for the molecular biologist. *Mol Cell*. — SLiM paradigm review
- Kumar et al. (2022). ELM — the Eukaryotic Linear Motif resource. *Nucleic Acids Res*. — SLiM database
- Tesei et al. (2024). CALVADOS 2 force field. *PNAS*. — IDP simulation methodology
- Lum & Bhatt (2018). HintHSV interactome. — Primary PPI dataset
- GTEx Consortium (2020). GTEx v8. *Science*. — Tissue expression data
- Davey et al. (2012). SLiMSearch. *Nucleic Acids Res*. — SLiM discovery tool
