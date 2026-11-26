# Supplementary Methods Protocol
## Ready-to-Execute Computational Steps

This document provides concrete commands, scripts, and parameter settings for each phase of the research proposal.

---

## Phase 1: Sequence Retrieval and Disorder Prediction

### 1A. UniProt Sequence Retrieval

Target protein accessions (HSV-1 strain 17):

| Protein | UniProt ID | Gene | Length |
|---------|-----------|------|--------|
| ICP34.5 | P36313 | RL1 | 263 aa |
| US3 | P04413 | US3 | 481 aa |
| pUL37 | P10218 | UL37 | 1123 aa |

Retrieve sequences via UniProt REST API or manual download in FASTA format.

### 1B. Negative Control Selection Criteria

Select HSV-1 proteins meeting ALL of:
1. Predicted IDR length ≥ 30 residues (run IUPred3 quick scan on all 72 HSV-1 ORFs)
2. NOT listed as neurovirulence determinant in literature
3. Represent distinct functional classes

Recommended controls:
- UL29/ICP8 (P04296) — ssDNA binding, large N-terminal IDR
- UL42 (P10186) — polymerase processivity factor
- US6/gD (P57083) — envelope glycoprotein, flexible N-terminus
- UL46/VP11/12 (P10228) — tegument phosphoprotein
- UL48/VP16 (P06492) — transactivator, acidic activation domain (known IDR)
- UL49/VP22 (P10216) — tegument, intercellular spread
- UL54/ICP27 (P10238) — mRNA export, RGG-box IDR
- US1/ICP22 (P04485) — regulatory protein

### 1C. IUPred3 Protocol

Run via web server (https://iupred3.elte.hu/) or local installation:
- Mode: "IUPred3 long disorder" (for constitutive disorder)
- Also run: "ANCHOR2" (for disordered binding regions / MoRFs)
- Threshold: 0.5 (standard)
- Save per-residue scores as CSV

### 1D. metapredict V2 Protocol

```python
# Python environment setup
# pip install metapredict

import metapredict as meta

sequences = {
    "ICP34.5": "SEQUENCE_HERE",
    "US3": "SEQUENCE_HERE",
    "pUL37": "SEQUENCE_HERE",
    # ... add controls
}

for name, seq in sequences.items():
    scores = meta.predict_disorder(seq)
    # Consensus IDRs (contiguous regions with score > 0.5 for >= 20 residues)
    idrs = meta.predict_disorder_domains(seq)
    print(f"{name}: {len(idrs.disordered_domain_boundaries)} IDRs found")
    # Save per-residue scores
```

### 1E. AlphaFold2 pLDDT Extraction

Download pre-computed AF2 models from AlphaFold Protein Structure Database:
- HSV-1 proteome models are available via UniProt integration
- Extract pLDDT from B-factor column of the PDB file

```python
from Bio.PDB import PDBParser

parser = PDBParser()
structure = parser.get_structure("ICP34.5", "AF-P36313-F1-model_v4.pdb")

plddt_scores = []
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.get_name() == "CA":
                    plddt_scores.append(atom.get_bfactor())

# Residues with pLDDT < 50 are predicted disordered
disordered = [i+1 for i, s in enumerate(plddt_scores) if s < 50]
```

### 1F. Consensus Disorder Calculation

A residue is consensus-disordered if ≥ 2/3 methods agree:

```python
import numpy as np

def consensus_disorder(iupred_scores, metapred_scores, plddt_scores,
                       iupred_thresh=0.5, metapred_thresh=0.5, plddt_thresh=50):
    n = len(iupred_scores)
    consensus = np.zeros(n)
    for i in range(n):
        votes = 0
        if iupred_scores[i] > iupred_thresh: votes += 1
        if metapred_scores[i] > metapred_thresh: votes += 1
        if plddt_scores[i] < plddt_thresh: votes += 1  # Note: low pLDDT = disordered
        consensus[i] = 1 if votes >= 2 else 0

    # Extract contiguous IDRs of >= 20 residues
    idrs = []
    start = None
    for i in range(n):
        if consensus[i] == 1 and start is None:
            start = i
        elif consensus[i] == 0 and start is not None:
            if i - start >= 20:
                idrs.append((start + 1, i))  # 1-indexed
            start = None
    if start is not None and n - start >= 20:
        idrs.append((start + 1, n))
    return consensus, idrs
```

---

## Phase 2: SLiM Discovery

### 2A. ELM Scanning

Use the ELM web server (http://elm.eu.org/) or batch mode:
- Input: FASTA sequences of IDR segments only (not full proteins)
- Settings: filter by "cell compartment" relevant to viral infection (cytoplasm, nucleus)
- Download results as TSV

For programmatic scanning, use the ELM REST API or local regex matching:

```python
import re

# Example ELM patterns (subset — full list from ELM database)
elm_patterns = {
    "LIG_LIR_Gen_1": r"[EDST].{0,2}[WFY]...[ILV]",       # LC3-interacting region
    "LIG_14-3-3_CanoR_1": r"R[^E]([^E])[ST][^E]P",         # 14-3-3 binding (mode 1)
    "LIG_TRAF2_1": r"[PSAT].[QE]E",                          # TRAF2 binding
    "LIG_PP1_RVXF_1": r"[RK].{0,1}[VIL][^P][FW]",          # PP1 binding
    "MOD_CDK_SPxK_1": r"[ST]P.[KR]",                         # CDK phosphorylation
    "DOC_MAPK_DCC_7": r"[KR]{2,4}.{2,6}[LIVMP].[LIVMP]",   # MAPK D-motif
    # Add all ~320 ELM classes
}

def scan_elm(sequence, patterns, idr_boundaries):
    """Scan IDR regions only for ELM matches."""
    hits = []
    for name, pattern in patterns.items():
        for match in re.finditer(pattern, sequence):
            start, end = match.start() + 1, match.end()  # 1-indexed
            # Check if match falls within an IDR
            in_idr = any(idr_start <= start and end <= idr_end
                        for idr_start, idr_end in idr_boundaries)
            if in_idr:
                hits.append({
                    "motif_class": name,
                    "start": start,
                    "end": end,
                    "sequence": match.group(),
                })
    return hits
```

### 2B. SLiMSearch Protocol

Use SLiMSearch web server (http://slim.ucd.ie/slimsearch/) or command line:
- Input: Full protein sequence
- Motif source: ELM motifs
- IUPred filter: ON (only report motifs in disordered regions)
- Conservation filter: OFF at this stage (applied separately in 2C)

### 2C. MAFFT Alignment for Conservation

```bash
# For each protein, align orthologs across alphaherpesviruses
mafft --localpair --maxiterate 1000 ICP34.5_orthologs.fasta > ICP34.5_aligned.fasta
mafft --localpair --maxiterate 1000 US3_orthologs.fasta > US3_aligned.fasta
mafft --localpair --maxiterate 1000 pUL37_orthologs.fasta > pUL37_aligned.fasta
```

### 2D. dN/dS Analysis with PAML

```
# codeml control file (codeml.ctl)
seqfile = US3_codon_aligned.phy    * codon alignment in PHYLIP format
treefile = alphaherpes_tree.nwk     * species tree
outfile = US3_paml_output.txt

noisy = 3
verbose = 1
runmode = 0         * user tree
seqtype = 1         * codons
CodonFreq = 2       * F3x4
model = 0           * one omega for all branches
NSsites = 0 1 2 7 8 * M0, M1a, M2a, M7, M8

fix_omega = 0
omega = 0.4         * initial value
```

Compare M1a (neutral) vs M2a (selection) and M7 (beta) vs M8 (beta+selection) using likelihood ratio tests.

### 2E. ML SLiM Scoring — Feature Extraction

```python
import pandas as pd
from xgboost import XGBClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score

def extract_features(slim_hit, disorder_scores, plddt_scores,
                     accessibility_scores, conservation_scores, sequence):
    """Extract feature vector for one SLiM hit."""
    start, end = slim_hit["start"] - 1, slim_hit["end"]  # 0-indexed
    flank = 5  # residues on each side

    features = {
        "disorder_mean": np.mean(disorder_scores[start:end]),
        "plddt_mean": np.mean(plddt_scores[start:end]),
        "conservation_mean": np.mean(conservation_scores[start:end]),
        "accessibility_mean": np.mean(accessibility_scores[start:end]),
        "motif_length": end - start,
        "flanking_disorder": np.mean(disorder_scores[max(0,start-flank):end+flank]),
        "flanking_charge_asymmetry": calculate_kappa(
            sequence[max(0,start-flank):end+flank]),
        "helix_propensity": calculate_helix_prop(start, end),  # from PSIPRED
        "sequence_complexity": shannon_entropy(
            sequence[max(0,start-flank):end+flank]),
    }
    return features

# Train with 5-fold stratified CV
model = XGBClassifier(n_estimators=200, max_depth=4, learning_rate=0.1,
                       scale_pos_weight=5)  # handle class imbalance
```

---

## Phase 3: Molecular Dynamics

### 3A. CALVADOS 2 Setup

```python
# Requires: OpenMM, CALVADOS2 package
# See: https://github.com/KULL-Centre/CALVADOS

from calvados import build_system, run_simulation

# Extract IDR sequence (e.g., ICP34.5 residues 1-147)
idr_sequence = "MSTPNP..."  # IDR sequence here

# Build CG system
system = build_system(
    sequence=idr_sequence,
    forcefield="calvados2",
    ionic_strength=0.15,  # 150 mM NaCl
    temperature=300,       # Kelvin
    box_size=50.0,         # nm (large enough for extended chain)
)

# Run 10 μs simulation (3 replicates)
for replicate in range(3):
    run_simulation(
        system=system,
        n_steps=10_000_000,    # with 1 fs timestep in CG = ~10 μs effective
        output_prefix=f"ICP34.5_IDR1_rep{replicate}",
        save_interval=10000,
    )
```

### 3B. SLiM Exposure Analysis

```python
import mdtraj as md

def slim_exposure(trajectory_file, topology_file, slim_residues):
    """Calculate fraction of frames where SLiM is solvent-exposed."""
    traj = md.load(trajectory_file, top=topology_file)
    sasa = md.shrake_rupley(traj, mode='residue')

    # SLiM exposure: all motif residues have SASA > threshold
    threshold = 0.1  # nm^2, minimal exposure
    exposure_frames = 0
    for frame_idx in range(traj.n_frames):
        slim_sasa = sasa[frame_idx, slim_residues]
        if all(s > threshold for s in slim_sasa):
            exposure_frames += 1

    return exposure_frames / traj.n_frames
```

---

## Phase 4: Network Analysis

### 4A. GTEx CNS Specificity Score

```python
import pandas as pd
import numpy as np

def tissue_specificity_tau(expression_row, brain_tissues, other_tissues):
    """Calculate tau (tissue specificity index) for brain vs. non-brain."""
    all_tissues = brain_tissues + other_tissues
    values = expression_row[all_tissues].values.astype(float)
    max_val = values.max()
    if max_val == 0:
        return 0
    normalized = values / max_val
    tau = (len(normalized) - sum(normalized)) / (len(normalized) - 1)
    # Positive if brain-enriched
    brain_mean = expression_row[brain_tissues].mean()
    other_mean = expression_row[other_tissues].mean()
    sign = 1 if brain_mean > other_mean else -1
    return tau * sign

brain_tissues = [
    "Brain_Cortex", "Brain_Hippocampus", "Brain_Cerebellum",
    "Brain_Substantia_nigra", "Brain_Spinal_cord_cervical_c-1",
    "Brain_Frontal_Cortex_BA9", "Brain_Anterior_cingulate_cortex_BA24"
]

other_tissues = [
    "Liver", "Lung", "Skin_Sun_Exposed_Lower_leg",
    "Whole_Blood", "Heart_Left_Ventricle", "Kidney_Cortex"
]
```

### 4B. Enrichment Test

```python
from scipy.stats import mannwhitneyu

def test_cns_enrichment(neurovirulence_scores, control_scores):
    """Mann-Whitney U test for CNS-targeting score enrichment."""
    stat, pvalue = mannwhitneyu(
        neurovirulence_scores, control_scores,
        alternative='greater'  # one-sided: neurovirulence > controls
    )
    effect_size = stat / (len(neurovirulence_scores) * len(control_scores))
    return {
        "U_statistic": stat,
        "p_value": pvalue,
        "rank_biserial_r": 2 * effect_size - 1,  # effect size
        "n_neurovirulence": len(neurovirulence_scores),
        "n_control": len(control_scores),
    }
```

---

## Phase 5: Switch Identification Criteria

### Decision Tree for Switch Calling

```
For each IDR with >= 2 high-confidence SLiMs:

1. Distance check:
   |position(SLiM_A) - position(SLiM_B)| <= 30 residues?
   NO  → Not a switch candidate
   YES → Continue

2. Domain overlap check:
   Do host binding domains for SLiM_A and SLiM_B have overlapping
   binding pockets (from 3did/PDB structures)?
   NO  → Possible independent functions, not a switch
   YES → Continue

3. Conformational exclusion check (from MD):
   Pearson correlation of SLiM_A exposure vs SLiM_B exposure < -0.3?
   YES → STRONG switch candidate (anti-correlated exposure)
   NO  → Continue

4. PTM toggle check:
   Predicted phosphosite within 10 residues of either SLiM?
   YES → MODERATE switch candidate (PTM-regulated)
   NO  → WEAK switch candidate (concentration-dependent only)
```

---

## Computational Resource Estimates

| Task | Resource | Time |
|------|----------|------|
| Disorder prediction (all proteins) | Standard laptop | 1–2 hours |
| ELM scanning + filtering | Standard laptop | 2–4 hours |
| MAFFT alignments | Standard laptop | < 1 hour |
| PAML dN/dS analysis | Standard workstation | 4–8 hours |
| XGBoost training + CV | Standard laptop | 1–2 hours |
| AF2 structure prediction (if needed) | Google Colab (free GPU) | 2–4 hours per protein |
| CALVADOS MD (10 μs × 3 reps × ~6 IDRs) | GPU workstation or HPC | 3–7 days total |
| FlexPepDock (top 5 candidates) | CPU cluster (50 cores) | 1–2 days |
| Network analysis | Standard laptop | 2–4 hours |
| Visualization and figures | Standard laptop | 3–5 days |

**Total compute:** Achievable with one GPU workstation + Google Colab. No supercomputer access required.
