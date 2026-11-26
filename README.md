# HSV-CNS-Interactome

**SLiM-Mediated Host Manipulation by HSV-1 Neurovirulence Proteins: A Computational Interactome Analysis**

A computational biology project investigating how HSV-1 neurovirulence proteins exploit intrinsically disordered regions (IDRs) and short linear motifs (SLiMs) to hijack host CNS signaling pathways.

## Overview

This project systematically analyzes three HSV-1 neurovirulence-associated proteins — **ICP34.5** (P36313), **US3** (P04413), and **pUL37** (P10218) — against eight control HSV-1 proteins using a six-phase computational pipeline. The analysis integrates disorder prediction, motif scanning, structural accessibility, interactome construction, competitive binding analysis, and CNS-weighted network perturbation to test whether neurovirulence proteins preferentially target CNS pathways via molecular mimicry.

## Hypotheses

| ID | Hypothesis | Test | Result |
|----|-----------|------|--------|
| H1 | Neurovirulence proteins have higher CNS-targeting composite scores | Mann–Whitney U | p = 0.861 (n.s.) |
| H2 | Target proteins have higher switch density than controls | Mann–Whitney U | p = 0.908 (n.s.) |
| H3 | Target switches preferentially involve CNS-enriched pathways | Fisher's exact | p = 0.054 (n.s.) |

## Project Structure

```
HSV-CNS-Interactome/
├── phase1/                     # IDR prediction & charge patterning
│   ├── disorder_prediction.py
│   ├── disorder_predictions/   # Per-protein disorder JSON & CSV
│   ├── iupred3_data/           # IUPred3 raw predictions
│   ├── results/                # Disorder summaries, annotations
│   └── plots/                  # Disorder profile plots
│
├── phase2/                     # SLiM scanning & conservation
│   ├── elm_database.py         # 81 core ELM regex classes
│   ├── elm_expanded.py         # 115 expanded ELM classes
│   ├── slim_discovery.py
│   ├── conservation_analysis.py
│   ├── results/                # SLiM hits, conservation scores, pathway enrichment
│   └── plots/                  # SLiM maps, enrichment dashboards
│
├── phase3/                     # Structural context & interactome
│   ├── phase3_structural.py
│   ├── interactome_builder.py
│   ├── results/                # Interactome edges, network metrics, GraphML
│   └── plots/                  # Network visualizations
│
├── phase4/                     # Network analysis & perturbation
│   ├── phase4_complete.py
│   ├── results/                # Community detection, heat diffusion, robustness
│   └── plots/                  # Network analysis plots
│
├── phase5/                     # Competitive binding & switch candidates
│   ├── phase5_switches.py
│   ├── results/                # Switch candidates, hypothesis H2/H3 results
│   └── plots/                  # Switch summary plots
│
├── phase6/                     # Cross-phase integration & final analysis
│   ├── phase6_integration.py
│   ├── results/                # Final summaries, hypothesis results, Table 1
│   └── plots/                  # Overview figures
│
├── manuscript/
│   ├── manuscript_draft.md     # Full manuscript
│   └── figures/                # Publication-quality figures (300 DPI)
│       ├── Figure1_IDR_SLiM_landscape.png
│       ├── Figure2_Interactome_Hypotheses.png
│       ├── FigureS1_SLiM_exposure.png
│       ├── FigureS2_network_robustness.png
│       ├── FigureS3_community_detection.png
│       ├── FigureS4_heat_diffusion.png
│       └── FigureS5_pathway_diversity.png
│
├── utils/
│   └── power_analysis.py       # Post-hoc power simulation
│
├── methods-protocol.md         # Detailed methods protocol
├── research-proposal.md        # Original research proposal
└── execution-checklist.md      # Phase execution tracking
```

## Methods Summary

### Phase 1 — Disorder Prediction
Consensus IDR prediction using IUPred3, composition-based, and complexity-based methods (≥ 2/3 agreement). Charge patterning analysis via localCIDER (FCR, NCPR, kappa, Das–Pappu classification).

### Phase 2 — SLiM Scanning & Conservation
Regex-based scanning against 196 ELM database classes (81 core + 115 expanded). Multi-criteria composite scoring:

$$S_i = 0.30 \cdot f_{\text{dis}}(i) + 0.20 \cdot f_{\text{acc}}(i) + 0.20 \cdot f_{\text{con}}(i) + 0.15 \cdot f_{\text{flex}}(i) + 0.15 \cdot f_{\text{cns}}(i)$$

### Phase 3 — Structural Context & Interactome
AlphaFold pLDDT-based accessibility estimation. CNS-weighted bipartite interactome (190 nodes, 614 edges) constructed via SLiM–domain interaction mapping.

### Phase 4 — Network Analysis
Louvain community detection (Q = 0.42), hub/bottleneck identification, heat diffusion perturbation (simplified HotNet2, α = 0.3), and network robustness analysis.

### Phase 5 — Competitive Binding
Switch candidate identification (competition > 0.35, different classes/pathways, distance ≤ 10 residues). H2 and H3 hypothesis testing.

### Phase 6 — Integration
Cross-phase data integration, Mann–Whitney U and Fisher's exact tests with post-hoc power analysis via simulation (5,000 permutations, seed = 42).

## Key Findings

- ICP34.5 has the highest disorder fraction (64.9%) among all 11 proteins
- All three hypotheses were non-significant, but the study was severely underpowered (n = 3 vs n = 8; MDE = d = 2.25 for 80% power)
- H3 (CNS pathway enrichment) approached significance (p = 0.054, OR = 1.23), suggesting pathway-specific rather than broad CNS targeting
- 591 switch candidates identified across target proteins, with ICP34.5 showing the highest density

## Computational Environment

| Component | Version |
|-----------|---------|
| Python | 3.11.8 |
| NumPy | 2.4.6 |
| SciPy | 1.17.1 |
| NetworkX | 3.6.1 |
| localCIDER | 0.1.21 |
| matplotlib | (figure generation) |

**Hardware:** AMD Ryzen 7 7840HS, 16 GB RAM  
**Reproducibility:** All stochastic procedures use `seed = 42`

## UniProt Accessions

| Protein | Accession | Category |
|---------|-----------|----------|
| ICP34.5 | P36313 | Target |
| US3 | P04413 | Target |
| pUL37 | P10218 | Target |
| ICP22 (US1) | — | Control |
| ICP27 (UL54) | — | Control |
| UL42 | — | Control |
| VP11/12 (UL46) | — | Control |
| VP16 (UL48) | — | Control |
| VP22 (UL49) | — | Control |
| gD (US6) | — | Control |
| ICP8 (UL29) | — | Control |

## Reproducing the Analysis

Run each phase sequentially:

```bash
python phase1/phase1_complete.py
python phase2/phase2_complete.py
python phase3/phase3_structural.py
python phase4/phase4_complete.py
python phase5/phase5_switches.py
python phase6/phase6_integration.py
```

## License

This project is provided for academic and research purposes.

## Citation

If you use this work, please cite the manuscript:

> SLiM-Mediated Host Manipulation by HSV-1 Neurovirulence Proteins: A Computational Interactome Analysis
