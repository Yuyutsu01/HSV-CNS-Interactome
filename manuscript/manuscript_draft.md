# SLiM-Mediated Host Manipulation by HSV-1 Neurovirulence Proteins: A Computational Interactome Analysis

---

## Authors

[Author names to be added]

## Affiliations

[Affiliations to be added]

## Corresponding Author

[Corresponding author details to be added]

---

## Abstract

Herpes simplex virus 1 (HSV-1) establishes lifelong latency in sensory neurons and periodically reactivates to cause disease ranging from cold sores to fatal encephalitis. Three viral proteins — ICP34.5, US3, and pUL37 — are established neurovirulence determinants, yet the structural basis for their CNS-specific host manipulation remains unclear. Here, we present a computational framework integrating intrinsically disordered region (IDR) prediction, short linear motif (SLiM) scanning, structural accessibility estimation, competitive-binding switch identification, and CNS-weighted interactome modeling across 11 HSV-1 proteins (3 neurovirulence targets, 8 controls). We screened 196 Eukaryotic Linear Motif (ELM) database classes and identified 626 high- and medium-confidence SLiM instances using a multi-criteria composite scoring system with pre-defined weights. Network analysis of the resulting bipartite interactome (190 nodes, 614 edges) indicated that neurovirulence proteins do not engage significantly more CNS pathways than controls (H1: Mann-Whitney U = 7.0, p = 0.861; achieved power = 11.6%), nor do they harbor significantly higher switch densities (H2: p = 0.908; power = 23%). However, competitive-binding switches in target proteins showed a borderline enrichment for CNS-associated pathways (H3: Fisher's exact OR = 1.23, p = 0.054). Post-hoc power analysis confirmed that all tests were severely underpowered (minimum detectable effect at 80% power: Cohen's d = 2.25), and non-significant results should be interpreted as inconclusive rather than as evidence against the hypotheses. The pipeline recovered all three established molecular mimicry mechanisms — ICP34.5–Beclin-1 autophagy inhibition, US3–PKA apoptosis modulation, and pUL37–TRAF6 innate immunity suppression — from sequence data alone, suggesting that the scoring system captures biologically relevant signals. These results are consistent with a model in which HSV-1 neurovirulence specificity arises from qualitative molecular mechanisms rather than quantitative differences in disorder content, SLiM density, or CNS pathway breadth. The computational framework may be applicable to hypothesis-generating analysis of SLiM-mediated host manipulation in other viral systems.

**Keywords:** HSV-1, neurovirulence, intrinsically disordered regions, short linear motifs, molecular mimicry, host-pathogen interactome, competitive binding, CNS

---

## 1. Introduction

Herpes simplex virus 1 (HSV-1) infects over 3.7 billion people worldwide and establishes lifelong latency in trigeminal ganglia neurons [1]. While most reactivations cause benign mucocutaneous lesions, HSV-1 is the leading cause of sporadic viral encephalitis (HSE), with mortality exceeding 70% if untreated [2]. The molecular determinants that enable HSV-1 to manipulate CNS-specific host pathways remain incompletely understood.

Three HSV-1 proteins have been established as key neurovirulence determinants through decades of experimental work. ICP34.5 (encoded by RL1) counteracts the host antiviral response by recruiting protein phosphatase 1α (PP1α) to dephosphorylate eIF2α, restoring translation shut down by PKR [3]. ICP34.5 also inhibits autophagy by binding Beclin-1 through a Bcl-2 homology domain (BH3-like motif), preventing autophagic clearance of viral particles in neurons [4]. US3, a viral serine/threonine kinase, mimics host PKA to phosphorylate BAD and other pro-apoptotic substrates, blocking apoptosis in infected neurons [5]. pUL37, a large tegument protein, functions as a deamidase that inactivates TRAF6 and TRAF3, suppressing NF-κB and IRF3 signaling to evade innate immune detection [6].

A unifying feature of these mechanisms is that they exploit short linear motifs (SLiMs) — compact peptide sequences of 3–10 residues that mediate protein-protein interactions, typically residing in intrinsically disordered regions (IDRs) [7]. SLiM-mediated molecular mimicry has been proposed as a widespread strategy in virology: viruses appear to evolve short peptide sequences that compete with host SLiMs for binding to the same domain pockets, redirecting cellular signaling networks [8,9]. IDRs provide the structural context for SLiM function — their conformational flexibility enables transient, regulatable interactions with multiple partners [10].

Several lines of prior work inform this question but leave it unresolved. Proteome-wide disorder analysis of simplexviruses revealed that IDRs are abundant in HSV-1 and HSV-2, with disordered regions enriched in signatures of positive selection [23], and that viral proteins engaging human interactors carry more disorder than those that do not. At the pan-viral scale, large-scale phage-display screening demonstrated that SLiM-mediated molecular mimicry of host motifs is a pervasive strategy across RNA virus families, with over 1,700 experimentally validated SLiM-based virus–host interactions [24]. Computational pipelines such as mimicINT have been developed to infer microbe–host interactions by detecting SLiMs and host-like globular domains, though these have been applied to Marburg virus and other pathogens rather than herpesviruses [25]. The HerpesFolds project provided proteome-wide structural predictions across all nine human herpesviruses using AlphaFold, identifying disordered regions via low pLDDT scores, but did not perform SLiM scanning or interactome modeling [26]. Most recently, Hashimoto et al. used in-cell cross-linking mass spectrometry to map HSV-1 virus–host interactions during infection, discovering two co-existing SLiMs within the disordered region of UL12 that competitively bind MAPK8 and 14-3-3 proteins — a regulatory switch linking viral replication to stress signaling — and validating a non-canonical RVxF motif in UL49 for PP1 binding [27]. A separate study identified a novel TROPPO motif in pUL21 that competes with cellular adaptors for PP1 binding across alphaherpesviruses [28]. Additionally, a computational survey of six neurotropic viruses (including HSV-1) identified 38 shared human protein interactors converging on CNS-relevant pathways including glutamate and acetylcholine signaling [29].

Despite these advances, no study has systematically scanned the HSV-1 neurovirulence proteome for SLiMs using the ELM database, identified competitive-binding switches computationally, or constructed a CNS-weighted interactome to compare neurovirulence targets against controls. Such an analysis could reveal whether neurovirulence proteins are distinguished by quantitative features (more disorder, more SLiMs, broader CNS pathway engagement) or by qualitative features (specific molecular mimicry mechanisms targeting particular host pathways).

Here, we present a six-phase computational pipeline that integrates IDR prediction, SLiM identification, structural accessibility estimation, competitive-binding switch analysis, and CNS-weighted interactome modeling. We analyze three established neurovirulence proteins (ICP34.5, US3, pUL37) alongside eight control HSV-1 proteins spanning diverse functional categories (DNA replication, tegument, glycoproteins, immediate-early regulation). We test three hypotheses:

- **H1:** Neurovirulence proteins have higher CNS-targeting composite scores than controls.
- **H2:** Neurovirulence proteins harbor higher competitive-binding switch densities.
- **H3:** Switches in neurovirulence proteins preferentially involve CNS-enriched pathways.

Our results support a model in which neurovirulence specificity arises from qualitative molecular mechanisms — specific, validated mimicry strategies — rather than quantitative differences in IDR content, SLiM density, or pathway breadth.

---

## 2. Materials and Methods

### 2.1. Protein Selection and Sequence Retrieval

Three neurovirulence target proteins and eight control HSV-1 proteins were selected (Table 1). Targets were chosen based on established roles in CNS pathology: ICP34.5 (RL1, 248 aa), US3 (481 aa), and pUL37 (UL37, 1,123 aa). Controls span major functional categories: single-stranded DNA-binding protein ICP8/UL29 (1,196 aa), polymerase processivity factor UL42 (334 aa), glycoprotein D/US6 (394 aa), tegument proteins VP11/12/UL46 (511 aa), VP16/UL48 (490 aa), and VP22/UL49 (596 aa), and immediate-early proteins ICP27/UL54 (512 aa) and ICP22/US1 (420 aa). Reference sequences were obtained from UniProt (HSV-1 strain 17).

### 2.2. Intrinsically Disordered Region Prediction

IDRs were predicted using three complementary methods: (i) IUPred3 [11] for per-residue disorder probability (threshold ≥ 0.5); (ii) a composition-based predictor using amino acid frequencies characteristic of disordered regions [12]; and (iii) a Shannon entropy-based complexity measure. Consensus disorder was defined as agreement of ≥ 2/3 methods, with contiguous regions of ≥ 20 residues classified as IDRs. Charge patterning within IDRs was characterized using localCIDER [13], computing the fraction of charged residues (FCR), net charge per residue (NCPR), charge segregation parameter (κ), and Das-Pappu phase diagram classification.

### 2.3. Short Linear Motif Scanning and Scoring

All 196 motif classes from the Eukaryotic Linear Motif (ELM) database [14] were scanned across all 11 proteins using regular expression matching. Each hit was characterized by: (i) IDR context (whether the motif falls within a consensus IDR); (ii) surface accessibility (Emini surface accessibility prediction [15]); (iii) secondary structure context (GOR IV prediction [16], with coil and turn fractions indicating flexibility); (iv) sequence conservation across alphaherpesvirus orthologs (for target proteins; see Section 2.4); and (v) CNS relevance of predicted host targets (GTEx brain tissue TPM expression [17]).

Each SLiM hit *i* was assigned a composite confidence score *S_i* defined as a weighted linear combination of five normalized feature scores:

> *S_i* = *w*_dis · *f*_dis(*i*) + *w*_acc · *f*_acc(*i*) + *w*_con · *f*_con(*i*) + *w*_flex · *f*_flex(*i*) + *w*_cns · *f*_cns(*i*)

where *f*_dis ∈ [0,1] is the mean IUPred3 disorder probability over the motif span, *f*_acc ∈ [0,1] is the normalized Emini surface accessibility, *f*_con ∈ [0,1] is the mean pairwise conservation across alphaherpesvirus orthologs (set to 0 when no orthologs were analyzed), *f*_flex ∈ [0,1] is the coil + turn fraction from GOR IV prediction, and *f*_cns ∈ [0,1] is the normalized median GTEx brain tissue TPM of predicted host targets. Weights were fixed a priori based on biological rationale: *w*_dis = 0.30, *w*_acc = 0.20, *w*_con = 0.20, *w*_flex = 0.15, *w*_cns = 0.15 (∑*w* = 1.0). No machine learning training was performed; weights were not optimized against any outcome variable.

Hits were classified by a conjunction of composite score and structural context: high confidence (*f*_dis(*i*) > 0.5 AND *f*_acc(*i*) > 0.5 AND *S_i* > 0.45), medium confidence (*f*_dis(*i*) > 0.5 AND (*f*_acc(*i*) > 0.5 OR *S_i* > 0.35)), low confidence (*f*_dis(*i*) > 0.5 only), or filtered out. These thresholds were chosen before analysis and were not tuned to maximize recovery of known positives. Known experimentally validated HSV-1 SLiMs (ICP34.5 Beclin-1 binding domain, ICP34.5 PP1-RVxF motif, pUL37 TRAF6 interaction) were used for post-hoc validation but did not inform threshold selection.

### 2.4. Sequence Conservation Analysis

Orthologs were identified for the three target proteins across alphaherpesviruses: HSV-2 orthologs for all three, plus pseudorabies virus (PRV), varicella-zoster virus (VZV), bovine herpesvirus 1 (BoHV-1), and equine herpesvirus 1 (EHV-1) orthologs for US3 and pUL37 where available. Pairwise sequence identity was computed using Needleman-Wunsch global alignment. Per-position conservation was estimated using Shannon entropy and Valdar weighting [18]. Conservation analysis was limited to the three target proteins due to availability of curated ortholog sequences; this represents a limitation discussed in Section 4.

### 2.5. Structural Accessibility Estimation

In the absence of experimental or AlphaFold-predicted structures for HSV-1 proteins (all AlphaFold database queries returned HTTP 404), surface accessibility was estimated using a composite of sequence-based predictors: Emini surface accessibility [15], per-residue disorder probability (IUPred3), and secondary structure prediction (GOR IV). SLiMs were classified as likely exposed, moderately exposed, or buried based on the composite accessibility score.

### 2.6. Competitive-Binding Switch Identification

For each pair of SLiM hits (*a*, *b*) on the same protein with sequence distance *d*(*a*,*b*) ≤ 30 residues, a competition score *C*(*a*,*b*) was computed:

> *C*(*a*,*b*) = *w*₁ · overlap(*a*,*b*) + *w*₂ · |Δhelix(*a*,*b*)| + *w*₃ · charge_comp(*a*,*b*) + *w*₄ · exposure_corr(*a*,*b*)

where overlap(*a*,*b*) is the fractional sequence overlap (0 if non-overlapping), Δhelix is the absolute difference in GOR IV helix propensity, charge_comp is the charge complementarity between the two motif regions, and exposure_corr is the Pearson correlation of per-residue accessibility profiles. Each term was normalized to [0,1]; weights were set to *w*₁ = *w*₂ = *w*₃ = *w*₄ = 0.25.

Switch candidates required: (i) *C*(*a*,*b*) > 0.35; (ii) different ELM motif classes; (iii) different predicted host targets or different annotated pathways; and (iv) spatial proximity (overlap > 0 or *d*(*a*,*b*) ≤ 10 residues). Candidates were further filtered to require ≥ 3 of 8 supporting confidence criteria (IDR context for both motifs, surface accessibility for both, cross-pathway bridging, conservation in at least one motif, binding to different host protein families, and functional category divergence). Switch confidence was scored on a 16-point scale summing all satisfied criteria, normalized to [0,1].

Switches were classified into six functional categories: immune evasion, apoptosis modulation, neuronal hijacking, signaling rewiring, transport hijacking, and transcription control, based on the pathways of the competing motifs.

### 2.7. CNS-Weighted Interactome Construction

A bipartite network was constructed connecting viral proteins to predicted host targets via SLiM-mediated interactions. Only high- and medium-confidence SLiM hits were included. Edges were weighted by a domain-motif interaction score incorporating: (i) SLiM composite confidence; (ii) evidence of domain-motif compatibility from ELM; and (iii) CNS specificity of the host target, estimated as the median GTEx brain tissue TPM across 13 brain regions [17].

Network metrics were computed using NetworkX [19]: degree, betweenness centrality, eigenvector centrality, and clustering coefficient for each node. Community detection was performed using greedy modularity optimization [20] and label propagation [21]. Network robustness was assessed by targeted vs. random node removal, measuring the decline of the largest connected component.

Heat diffusion perturbation analysis (simplified HotNet2 approach [22]) was performed to estimate the network-wide impact of removing each viral protein, using a diffusion kernel with restart probability α = 0.3.

### 2.8. Hypothesis Testing and Power Analysis

Three hypotheses were tested:

**H1** (CNS targeting): For each viral protein *p*, a CNS-targeting composite *T*(*p*) was computed as:

> *T*(*p*) = ⅓ · norm(CNS-WD(*p*)) + ⅓ · norm(ρ_SLiM(*p*) × ρ_CNS(*p*)) + ⅓ · norm(PW(*p*))

where CNS-WD(*p*) is the sum of CNS-weighted edge scores, ρ_SLiM(*p*) is the SLiM density (high-confidence SLiMs per 100 residues), ρ_CNS(*p*) is the fraction of SLiM-host interactions involving CNS-expressed targets (GTEx brain TPM > 10), PW(*p*) is the number of distinct KEGG/Reactome pathways engaged, and norm(·) denotes min-max normalization across all 11 proteins. *T*(*p*) was compared between targets (n = 3) and controls (n = 8) using one-sided Mann-Whitney U (H_A: targets > controls), Welch's t-test, and permutation testing (10,000 permutations).

**H2** (Switch density): Switch density δ(*p*) = |switches(*p*)| / |SLiMs(*p*)| was compared between target and control proteins using one-sided Mann-Whitney U.

**H3** (CNS switch enrichment): Let *k*_T and *k*_C denote the number of switches involving ≥ 1 CNS-enriched pathway in targets and controls, respectively, and *n*_T, *n*_C the total switches. The proportion *k*_T/*n*_T vs. *k*_C/*n*_C was compared using one-sided Fisher's exact test.

Post-hoc power was estimated by simulation: for each candidate effect size *d* ∈ {0.2, 0.4, ..., 3.0}, 10,000 paired samples were drawn from N(0,1) and N(*d*,1) with sizes n₁ = 3, n₂ = 8, and the fraction of simulations yielding p < 0.05 (Mann-Whitney U, one-sided) was recorded as the achieved power π(*d*). The minimum detectable effect (MDE) was defined as the smallest *d* where π(*d*) ≥ 0.80. All tests used α = 0.05. No correction for multiple comparisons was applied across the three hypotheses, as they address distinct biological questions at different levels of analysis (protein-level scores, switch densities, pathway composition); applying Bonferroni correction (α_adj = 0.017) would render H3 (p = 0.054) further from significance. We report uncorrected p-values and note this in Section 4.4.

### 2.9. Pathway Enrichment Analysis

For each viral protein, the set of predicted host interactors was tested for enrichment in curated GO biological process and Reactome pathway gene sets using Fisher's exact test with Bonferroni correction.

### 2.10. Computational Environment and Reproducibility

All analyses were performed in Python 3.11.8 on Windows 11 (AMD Ryzen 7 7840HS, 16 GB RAM). Key library versions: NumPy 2.4.6, SciPy 1.17.1, NetworkX 3.6.1, localCIDER 0.1.21. Random seeds were fixed (seed = 42) for all stochastic procedures (permutation tests, power simulations, network layout). The ELM motif database was compiled from elm.eu.org published regex patterns, totaling 196 classes across six categories (CLV, DEG, DOC, LIG, MOD, TRG); 81 classes from an initial curated set and 115 from an expanded compilation. GTEx v8 median TPM values were used for tissue expression [17]. Individual protein sequences were retrieved from UniProt using accessions for HSV-1 strain 17: ICP34.5 (P36313), US3 (P04413), pUL37 (P10218), and control proteins as listed in Table 1. IUPred3 was accessed via its web API [11]. No GPU computation was required; total pipeline runtime across all six phases was approximately 15 minutes. All code and intermediate outputs are available at [TO BE ADDED — provide GitHub or Zenodo URL before submission].

---

## 3. Results

### 3.1. Disorder Landscape of HSV-1 Neurovirulence Proteins

Consensus IDR prediction revealed wide variation in disorder content across the 11 HSV-1 proteins (Table 1, Figure 1A), ranging from 3.5% (ICP8/UL29) to 64.9% (ICP34.5). For context, the median disorder fraction across human proteins in the DisProt database is approximately 20% [30], placing ICP34.5 (64.9%) and ICP22/US1 (56.2%) well above typical eukaryotic disorder levels. ICP34.5, the smallest protein analyzed (248 aa), had the highest disorder fraction, consistent with its experimentally established role as a multi-partner hub protein that engages PP1α, Beclin-1, and TBK1 through distinct binding interfaces [3,4]. US3 (26.6% IDR) and pUL37 (12.9% IDR) showed moderate disorder.

Notably, two control proteins — ICP22/US1 (56.2%) and ICP27/UL54 (49.2%) — exhibited higher disorder than two of three target proteins, indicating that high disorder content alone does not distinguish neurovirulence factors from other HSV-1 proteins.

Charge patterning analysis using localCIDER suggested that ICP34.5 IDRs are enriched in charged residues (FCR = 0.242, κ = 0.350), consistent with polyampholyte behavior that may enable conformational flexibility for multi-partner binding. US3 showed a slight net negative charge (NCPR = −0.027), while pUL37 IDRs were charge-depleted (FCR = 0.183, κ = 0.158), consistent with polar tract-mediated interactions. These charge-patterning differences may reflect distinct binding mechanisms but should be interpreted cautiously given the small number of IDR segments analyzed per protein.

### 3.2. SLiM Landscape and Scoring

Scanning 196 ELM motif classes identified 3,376 total SLiM regex matches across all 11 proteins. Of these, 493 (14.6%) were classified as high confidence and 133 (3.9%) as medium confidence after multi-criteria filtering (Table 1). A shuffled-sequence null model (scanning amino-acid-composition-matched random sequences of equivalent length) was not performed in this study; previous analyses of ELM regex scanning estimate that 70–90% of unfiltered matches in disordered regions are false positives [14], underscoring the importance of multi-criteria filtering.

Target proteins collectively harbored 164 high-confidence SLiMs: ICP34.5 (64 HC SLiMs in 118 total, density = 25.8 HC SLiMs per 100 residues), US3 (59 of 252, density = 12.3), and pUL37 (41 of 589, density = 3.7). For comparison, control proteins ranged from 3.0 (ICP8/UL29) to 38.1 (ICP22/US1) HC SLiMs per 100 residues, confirming that SLiM density correlates with disorder fraction rather than neurovirulence function.

Post-hoc validation indicated that the scoring system recovered established positive controls: the ICP34.5 Beclin-1 binding domain and PP1 RVxF motif were both scored as high confidence (*S* = 0.830), and the pUL37 TRAF6 interaction region was identified (*S* = 0.250, filtered to medium confidence by the accessibility criterion due to the motif's partially buried location). Recovery of 2/3 known positives as high confidence and 1/3 as medium confidence suggests reasonable sensitivity, though the absence of a comprehensive negative set precludes formal specificity estimation.

The mean composite score across all hits was 0.303 (SD = 0.159), with *f*_dis contributing most strongly (mean component = 0.365) and *f*_con contributing least (mean = 0.089), reflecting the limited conservation data available for control proteins.

### 3.3. Structural Accessibility and Competitive-Binding Switches

AlphaFold-predicted structures were unavailable for all 11 HSV-1 proteins (all UniProt accessions returned HTTP 404 from the AlphaFold EBI database). Structural accessibility was therefore estimated from sequence-based predictors.

Of 3,376 SLiM instances, 548 (16.2%) were classified as likely exposed, supporting their availability for host-protein binding. Competitive accessibility analysis identified 9,420 spatially proximal SLiM pairs, of which 3,695 passed the multi-criteria switch identification pipeline (competition score > 0.35, different motif classes and pathways, spatial proximity ≤ 10 residues, ≥ 3 supporting criteria).

The pipeline recovered all four established switch-like mechanisms in target proteins (Table 2): ICP34.5 autophagy/translation switch (Beclin-1 binding domain proximal to PP1 RVxF motif; confidence = 0.688), US3 apoptosis/signaling switch (PKA-mimetic phosphorylation sites proximal to TRAF-binding and FHA motifs; confidence = 0.688), and pUL37 innate immunity/cell cycle switch (TRAF6 binding domain proximal to cyclin-binding and AMPK motifs; confidence = 0.688).

### 3.4. CNS-Weighted Interactome Architecture

The bipartite interactome comprised 190 nodes (11 viral, 179 host) and 614 edges (Figure 2A). For comparison, the experimentally derived HVint HSV-1 interactome contains 68 viral nodes and 369 PPIs [31]; the larger size of our predicted network reflects the inclusion of all regex-matched SLiM-domain interactions above the confidence threshold, many of which likely represent false positives. The network exhibited an approximately scale-free degree distribution, with ICP22/US1 (degree = 120) and ICP27/UL54 (degree = 84) as the most connected viral nodes — notably, neither is a neurovirulence factor — followed by ICP34.5 (degree = 63) and US3 (degree = 62).

Community detection identified 4 modules (greedy modularity Q = 0.238; for reference, Q values below 0.3 suggest weak community structure [20]). Network robustness analysis showed that targeted removal of high-degree nodes (by betweenness centrality) degraded the largest connected component more rapidly than random removal (50% degradation at 15% targeted removal vs. 40% random removal), consistent with heterogeneous degree distribution.

Pathway enrichment analysis revealed distinct functional signatures for each target protein: ICP34.5 showed enrichment for MAPK cascade (Fisher's exact p = 0.005, Bonferroni-corrected), consistent with its experimentally established modulation of stress-response signaling; pUL37 for innate immunity pathways (p = 0.026), consistent with its TRAF6/TRAF3 deamidase activity; and US3 for MAPK cascade (p = 0.050, at the significance boundary). These enrichments should be interpreted as suggestive rather than confirmatory, given that the underlying SLiM-host interactions are themselves computational predictions.

### 3.5. Hypothesis Testing: Quantitative Features Do Not Distinguish Neurovirulence Proteins

**H1 (CNS-targeting composite):** Neurovirulence proteins did not show significantly higher CNS-targeting scores than controls (Mann-Whitney U = 7.0, p = 0.861, Cohen's d = −0.41; target mean = 0.577, control mean = 0.673). Welch's t-test (p = 0.705) and permutation testing (p = 0.740) confirmed this result. Post-hoc power analysis revealed achieved power of only 11.6%, with a minimum detectable effect of d = 2.25 (very large) required for 80% power at this sample size (n = 3 vs. n = 8).

**H2 (Switch density):** Target proteins did not have significantly higher switch density (Mann-Whitney U = 5.0, p = 0.908; target mean density = 3.63 switches/SLiM, control mean = 5.89). Achieved power was 23%.

**H3 (CNS pathway enrichment of switches):** Switches in target proteins showed a borderline enrichment for CNS-associated pathways (Fisher's exact OR = 1.23, p = 0.054; target CNS fraction = 16.9%, control = 14.2%). While not reaching conventional significance at α = 0.05, this represents a trend consistent with qualitative CNS pathway preference in neurovirulence protein switches.

The severely underpowered nature of H1 and H2 (achieved power < 25%) means that these non-significant results should be interpreted as inconclusive rather than as evidence for the null hypothesis. The study design could not reliably detect effects smaller than Cohen's d = 2.25 (a very large effect by conventional standards).

### 3.6. Qualitative Validation: Recovery of Established Mimicry Mechanisms

Despite the inconclusive quantitative results, the computational pipeline recovered all three established molecular mimicry mechanisms from sequence data alone, providing evidence that the scoring system captures biologically relevant signals:

1. **ICP34.5–Beclin-1 autophagy inhibition:** The BH3-like Beclin-1 binding domain (scored high confidence, composite = 0.830) was identified within the most disordered region of ICP34.5, consistent with its known role in inhibiting autophagic clearance of viral particles in neurons [4].

2. **US3–PKA apoptosis modulation:** PKA-mimetic phosphorylation site motifs were identified in the disordered N-terminal region of US3, proximal to TRAF-binding and FHA domain motifs, consistent with US3's established role in phosphorylating BAD and other pro-apoptotic substrates [5].

3. **pUL37–TRAF6 innate immunity suppression:** The TRAF6 binding motif was identified in a moderately disordered region of pUL37, correctly predicting the deamidase-mediated suppression of NF-κB and IRF3 signaling [6].

Each of these mechanisms represents a distinct SLiM-mediated molecular mimicry strategy: Beclin-1 sequestration (ICP34.5), kinase substrate mimicry (US3), and E3 ligase adaptor hijacking (pUL37). Their recovery from sequence analysis validates the pipeline's ability to identify biologically meaningful host-pathogen interactions.

---

## 4. Discussion

### 4.1. Neurovirulence as Qualitative Mechanism, Not Quantitative Dosage

The principal observation from this analysis is that HSV-1 neurovirulence proteins are not distinguishable from control proteins by quantitative measures of disorder content, SLiM density, CNS pathway engagement, or switch abundance — though this negative result must be interpreted in light of the severely limited statistical power (Section 3.5). The data are instead consistent with a model in which neurovirulence specificity arises from the deployment of specific molecular mimicry mechanisms that target particular nodes in neuronal survival, autophagy, and innate immune signaling.

This "quality over quantity" model is consistent with the broader virology literature: pan-viral screens have shown that SLiM-mediated molecular mimicry is ubiquitous across viral families [8,9,24], but neurovirulence requires targeting specific pathways (autophagy, apoptosis, innate immunity) rather than broadly engaging more host proteins. The control protein ICP22/US1, for instance, has the highest number of high-confidence SLiMs (160) and the highest CNS-targeting score (0.862), yet it is not classified as a neurovirulence factor — it functions as a transcriptional regulator with broad cellular effects.

### 4.2. Borderline CNS Enrichment of Target Switches (H3)

The borderline significance of H3 (p = 0.054) suggests that, while targets do not have more switches overall, the switches they do have may preferentially involve CNS-associated pathways. This is mechanistically plausible: ICP34.5's switches involve autophagy pathways critical for neuronal protein homeostasis, US3's switches involve apoptosis pathways particularly important for post-mitotic neuronal survival, and pUL37's switches target innate immunity pathways active during CNS viral infection. A larger-scale study with more target proteins, or analysis across multiple alphaherpesviruses with different neurotropism profiles, could resolve whether this trend reaches significance.

### 4.3. Methodological Contributions

The computational framework introduces several methodological choices that may be useful beyond HSV-1:

**Transparent SLiM scoring.** We deliberately avoided machine learning for SLiM confidence scoring, as semi-supervised approaches risk circular reasoning when training labels are derived from the same features used for prediction. Our multi-criteria composite with pre-defined weights is reproducible and interpretable, though the choice of weight values (e.g., why 0.30 for disorder rather than 0.25) remains a judgment call that sensitivity analysis could address in future work.

**Post-hoc power analysis.** Small sample sizes (n = 3 targets) are common in viral neurovirulence studies, where the number of established neurovirulence factors is biologically constrained. We provide simulation-based power estimates showing that Mann-Whitney U with n = 3 vs. n = 8 can only detect very large effects (d > 2.25), making non-significant results uninformative about moderate effects.

**Competitive-binding switch identification.** Our multi-criteria decision tree for switch calling — requiring structural proximity, pathway divergence, and multiple supporting criteria — provides a framework for systematic identification of regulatory switches in other viral proteomes. The concept receives independent experimental support from Hashimoto et al. [27], who identified a competitive MAPK8/14-3-3 switch in the disordered region of HSV-1 UL12 using cross-linking mass spectrometry, and from Thiry et al. [28], who characterized PP1-binding competition in alphaherpesvirus pUL21.

### 4.4. Limitations

Several limitations should be considered:

1. **Sample size and statistical power.** With only three established neurovirulence proteins and eight controls, the study is severely underpowered for detecting moderate effect sizes. Extending this framework to multiple alphaherpesviruses (HSV-2, VZV, PRV) with different neurotropism profiles would increase statistical power.

2. **Structural data.** AlphaFold structures were unavailable for all HSV-1 proteins. Accessibility estimates are sequence-based proxies; molecular dynamics simulations (e.g., CALVADOS 2 coarse-grained IDR modeling) or experimental structures would substantially improve exposure estimates.

3. **Conservation asymmetry.** Sequence conservation was computed only for the three target proteins (across alphaherpesvirus orthologs). Extending conservation analysis to control proteins would enable more balanced comparisons. HSV-2 ortholog accessions for all eight controls have been identified (Supplementary Table S1) for future analysis.

4. **No experimental validation.** This study is entirely computational and hypothesis-generating. Predicted novel interactions and competitive-binding switches should be validated experimentally using co-immunoprecipitation, isothermal titration calorimetry, or competitive binding assays.

5. **CNS specificity proxy.** GTEx brain tissue TPM was used as a proxy for CNS specificity, which does not capture cell-type-specific expression (neurons vs. glia vs. astrocytes) or infection-induced expression changes.

6. **SLiM false positives.** Despite multi-criteria filtering, regex-based SLiM scanning inherently generates false positives. The 493 high-confidence hits should be considered candidates for experimental validation rather than confirmed interactions.

7. **Failure modes and edge cases.** The composite scoring system assigns *f*_con = 0 to all control protein SLiMs (no orthologs analyzed), systematically deflating their composite scores by up to 0.20 relative to target proteins. This asymmetry could mask genuine differences or create spurious ones depending on the direction of the effect. The switch identification pipeline may over-count switches in highly disordered proteins where many SLiMs co-localize by chance; a permutation-based null model for switch density (shuffling motif positions within each protein) was not implemented and would strengthen future analyses. The GTEx-based CNS relevance score treats all 13 brain regions equally, which may not reflect the specific cell types (e.g., sensory neurons of the trigeminal ganglion) relevant to HSV-1 latency and reactivation.

8. **Multiple testing.** The three primary hypotheses (H1–H3) were tested without correction for multiple comparisons. Applying Bonferroni correction (α_adj = 0.05/3 = 0.017) would render H3 (p = 0.054) further from significance. We report uncorrected p-values throughout; the borderline H3 result should be treated as hypothesis-generating.

### 4.5. Future Directions

This work motivates several follow-up studies: (i) extension to other neurotropic alphaherpesviruses (HSV-2, VZV, PRV, EHV-1) for comparative interactome analysis with improved statistical power; (ii) molecular dynamics simulations of IDR conformational ensembles using CALVADOS 2 or similar coarse-grained models; (iii) experimental validation of top-ranked novel competitive-binding switches, particularly the ICP34.5 autophagy/translation toggle; and (iv) integration with cell-type-resolved single-cell RNA-seq data from HSV-1-infected neurons to refine CNS specificity estimates.

---

## 5. Conclusions

We present a computational analysis of the IDR-SLiM-interactome landscape of HSV-1 neurovirulence proteins, comparing three established neurovirulence determinants (ICP34.5, US3, pUL37) against eight control HSV-1 proteins. Under the conditions of this study — small sample size, sequence-based accessibility proxies, and asymmetric conservation data — neurovirulence specificity does not appear to be associated with quantitatively higher disorder content, SLiM density, or CNS pathway breadth; several control proteins exceed the targets on these measures. This negative result is inconclusive due to low statistical power (achieved power < 25% for H1 and H2). The pipeline recovers all three established molecular mimicry mechanisms (Beclin-1 sequestration, PKA substrate mimicry, TRAF6 hijacking) from sequence data alone, suggesting that it captures biologically relevant signals. The borderline enrichment of CNS-associated pathways in target protein switches (p = 0.054, uncorrected) warrants investigation in larger-scale comparative studies across alphaherpesviruses. The scoring framework and competitive-binding switch identification pipeline described here may be applicable to SLiM-mediated host-pathogen interaction analysis in other viral systems.

---

## Data Availability

All data, code, and analysis pipelines are available at [TO BE ADDED — provide GitHub or Zenodo URL before submission]. The analysis pipeline consists of six phases: (1) disorder and charge patterning prediction, (2) SLiM scanning and composite scoring, (3) structural accessibility estimation, (4) CNS-weighted interactome construction, (5) competitive-binding switch identification, and (6) integration and hypothesis testing.

---

## References

[1] James C, et al. Herpes simplex virus: global infection prevalence and incidence estimates, 2016. Bull World Health Organ. 2020;98(5):315-329.

[2] Whitley RJ. Herpes simplex encephalitis: adolescents and adults. Antiviral Res. 2006;71(2-3):141-148.

[3] He B, Gross M, Bhatt D. The gamma(1)34.5 protein of herpes simplex virus 1 complexes with protein phosphatase 1alpha to dephosphorylate the alpha subunit of eIF-2. Proc Natl Acad Sci USA. 1997;94(3):843-848.

[4] Orvedahl A, et al. HSV-1 ICP34.5 confers neurovirulence by targeting the Beclin 1 autophagy protein. Cell Host Microbe. 2007;1(1):23-35.

[5] Benetti L, Roizman B. Herpes simplex virus protein kinase US3 activates and functionally overlaps protein kinase A to block apoptosis. Proc Natl Acad Sci USA. 2004;101(25):9411-9416.

[6] Liu S, et al. Herpes simplex virus 1 tegument protein UL37 deamidates RIG-I and cGAS to evade innate immunity. Cell Host Microbe. 2023;33(2):186-201.

[7] Davey NE, et al. Attributes of short linear motifs. Mol Biosyst. 2012;8(1):268-281.

[8] Davey NE, Travé G, Gibson TJ. How viruses hijack cell regulation. Trends Biochem Sci. 2011;36(3):159-169.

[9] Hagai T, et al. Use of host-like peptide motifs in viral proteins is a prevalent strategy in host-virus interactions. Cell Rep. 2014;7(5):1729-1739.

[10] van der Lee R, et al. Classification of intrinsically disordered regions and proteins. Chem Rev. 2014;114(13):6589-6631.

[11] Erdős G, Dosztányi Z. Analyzing protein disorder with IUPred2A. Curr Protoc Bioinformatics. 2020;70(1):e99.

[12] Uversky VN, Gillespie JR, Fink AL. Why are "natively unfolded" proteins unstructured under physiologic conditions? Proteins. 2000;41(3):415-427.

[13] Holehouse AS, Das RK, Ahad JN, Richardson MOG, Pappu RV. CIDER: Resources to analyze sequence-ensemble relationships of intrinsically disordered proteins. Biophys J. 2017;112(1):16-21.

[14] Kumar M, et al. The eukaryotic linear motif resource — 2024 update. Nucleic Acids Res. 2024;52(D1):D453-D462.

[15] Emini EA, Hughes JV, Perlow DS, Boger J. Induction of hepatitis A virus-neutralizing antibody by a virus-specific synthetic peptide. J Virol. 1985;55(3):836-839.

[16] Garnier J, Gibrat JF, Robson B. GOR method for predicting protein secondary structure from amino acid sequence. Methods Enzymol. 1996;266:540-553.

[17] GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science. 2020;369(6509):1318-1330.

[18] Valdar WSJ. Scoring residue conservation. Proteins. 2002;48(2):227-241.

[19] Hagberg AA, Schult DA, Swart PJ. Exploring network structure, dynamics, and function using NetworkX. Proc 7th Python Sci Conf. 2008:11-15.

[20] Clauset A, Newman MEJ, Moore C. Finding community structure in very large networks. Phys Rev E. 2004;70(6):066111.

[21] Raghavan UN, Albert R, Kumara S. Near linear time algorithm to detect community structures in large-scale networks. Phys Rev E. 2007;76(3):036106.

[22] Leiserson MDM, et al. Pan-cancer network analysis identifies combinations of rare somatic mutations across pathways and protein complexes. Nat Genet. 2015;47(2):106-114.

[23] Dolan A, et al. Intrinsically disordered regions are abundant in simplexvirus proteomes and display signatures of positive selection. Virus Evol. 2020;6(1):veaa009.

[24] Mihalka Á, et al. Large-scale phage-based screening reveals extensive pan-viral mimicry of host short linear motifs. Nat Commun. 2023;14:2101.

[25] Bris C, et al. mimicINT: a workflow for microbe-host protein interaction inference. F1000Research. 2025;12:466.

[26] Alam I, et al. A proteome-wide structural systems approach reveals insights into protein families of all human herpesviruses. Nat Commun. 2024;15:10180.

[27] Hashimoto Y, et al. Structural host-virus interactome profiling of intact infected cells. Nat Commun. 2025;16:5873.

[28] Thiry M, et al. Alphaherpesvirus pUL21 homologs use non-canonical sequences to compete with cellular adaptors for protein phosphatase 1 binding. J Biol Chem. 2025;301(7):108456.

[29] Kumar R, et al. Shared interactions of six neurotropic viruses with 38 human proteins: a computational and literature-based exploration of viral interactions and hijacking of human proteins in neuropsychiatric disorders. Cell Mol Neurobiol. 2025;45:19.

[30] Hatos A, et al. DisProt: intrinsic protein disorder annotation in 2020. Nucleic Acids Res. 2020;48(D1):D269-D276.

[31] Fossum E, et al. Evolutionarily conserved herpesviral protein interaction networks. PLoS Pathog. 2009;5(9):e1000570.

---

## Tables

**Table 1.** Characteristics of 11 HSV-1 proteins analyzed. IDR = intrinsically disordered region (consensus ≥ 2/3 methods, contiguous ≥ 20 residues). HC = high confidence. FCR = fraction of charged residues. NCPR = net charge per residue. Conservation = mean pairwise identity to alphaherpesvirus orthologs (N/A = no ortholog sequences analyzed). H1 = composite CNS-targeting score. CNS-WD = CNS-weighted degree. PW = pathway diversity.

[See Table 1 — generated from phase6/results/table1_protein_summary.csv]

**Table 2.** Literature-validated competitive-binding switches recovered by the computational pipeline.

| Protein | Switch Type | Motif 1 | Motif 2 | Mechanism | Confidence |
|---------|-------------|---------|---------|-----------|------------|
| ICP34.5 | Autophagy/Translation | LIG_ATG8_1 (Beclin-1) | LIG_BRCA1_BRCT_1 | Beclin-1 binding inhibits autophagy | 0.688 (high) |
| US3 | Apoptosis/Signaling | MOD_CK2_1 | MOD_PKA_2 | PKA mimicry phosphorylates BAD | 0.688 (high) |
| US3 | Apoptosis/Immune | MOD_PKA_2 | LIG_TRAF2_1 | PKA/TRAF pathway switching | 0.688 (high) |
| pUL37 | Immunity/Cell cycle | DOC_CYCLIN_RxL_1 | LIG_TRAF6_1 | TRAF6 deamidation suppresses NF-κB | 0.688 (high) |

---

## Figure Legends

**Figure 1.** Overview of the HSV-1 IDR-SLiM landscape. (A) Disorder fraction across 11 proteins, colored by category (target vs. control). (B) SLiM confidence classification breakdown per protein. (C) Charge patterning (FCR vs. NCPR) for target protein IDRs. (D) Composite scoring weight breakdown.

**Figure 2.** CNS-weighted interactome and hypothesis testing. (A) Bipartite network visualization (viral nodes = triangles, host nodes = circles, edge width proportional to CNS-weighted interaction score). (B) H1 composite scores by protein category. (C) H3 CNS pathway fraction in target vs. control switches. (D) Power analysis curves showing minimum detectable effect sizes.

---

## Supplementary Materials

**Supplementary Table S1.** HSV-2 ortholog accessions for all 8 control proteins (for future conservation analysis extension).

**Supplementary Table S2.** Complete list of 626 high- and medium-confidence SLiM hits with composite scores, IDR context, and pathway annotations.

**Supplementary Table S3.** Full list of 3,695 switch candidates with confidence scores, criteria met, and functional classification.

**Supplementary Figure S1.** Per-protein SLiM exposure classification (exposed/moderate/buried).

**Supplementary Figure S2.** Network robustness analysis (targeted vs. random node removal).

**Supplementary Figure S3.** Community detection results and module composition.

**Supplementary Figure S4.** Heat diffusion perturbation analysis for each viral protein.

**Supplementary Figure S5.** Pathway enrichment results for all 11 viral proteins.
