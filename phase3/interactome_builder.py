"""
Phase 3 — SLiM Structural Context, CNS-Weighted Interactome & Hypothesis Test

Combines Phase 3 (structural analysis) and Phase 4 (network) from the proposal:
  1. SLiM structural context (flanking residues, accessibility, secondary structure propensity)
  2. Domain-motif interaction mapping (which host domains bind each SLiM class)
  3. CNS-weighted interactome network (bipartite: viral protein ↔ host protein)
  4. Pathway enrichment (GO biological process + Reactome)
  5. Network topology analysis (degree, betweenness, module detection)
  6. Core hypothesis test: Mann-Whitney U comparing CNS-targeting scores
"""

import os
import json
import csv
import numpy as np
import pandas as pd
import networkx as nx
from collections import defaultdict
from scipy import stats
from itertools import combinations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

BASE_DIR = r"C:\Users\shiva\OneDrive\Desktop\HSV-CNS-Interactome"
PHASE1_DIR = os.path.join(BASE_DIR, "phase1")
PHASE2_DIR = os.path.join(BASE_DIR, "phase2")
PHASE3_DIR = os.path.join(BASE_DIR, "phase3")
RESULTS_DIR = os.path.join(PHASE3_DIR, "results")
PLOT_DIR = os.path.join(PHASE3_DIR, "plots")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)


# ═══════════════════════════════════════════════════════════════
# AMINO ACID PROPERTY SCALES
# ═══════════════════════════════════════════════════════════════

HYDROPHOBICITY = {  # Kyte-Doolittle
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'E': -3.5,
    'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
    'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
    'Y': -1.3, 'V': 4.2,
}

CHARGE = {
    'D': -1, 'E': -1, 'K': 1, 'R': 1, 'H': 0.5,
}

HELIX_PROPENSITY = {  # Chou-Fasman
    'A': 1.42, 'R': 0.98, 'N': 0.67, 'D': 1.01, 'C': 0.70, 'E': 1.51,
    'Q': 1.11, 'G': 0.57, 'H': 1.00, 'I': 1.08, 'L': 1.21, 'K': 1.16,
    'M': 1.45, 'F': 1.13, 'P': 0.57, 'S': 0.77, 'T': 0.83, 'W': 1.08,
    'Y': 0.69, 'V': 1.06,
}


# ═══════════════════════════════════════════════════════════════
# DOMAIN–MOTIF INTERACTION DATABASE
# Maps ELM classes to the host protein domains they bind
# Source: 3did, ELM annotations, literature
# ═══════════════════════════════════════════════════════════════

DOMAIN_MOTIF_MAP = {
    "LIG_LIR_Gen_1": {"domain": "Atg8/LC3", "pfam": "PF02991", "binding_mode": "hydrophobic_groove"},
    "LIG_LIR_Apic_2": {"domain": "Atg8/LC3", "pfam": "PF02991", "binding_mode": "hydrophobic_groove"},
    "LIG_ATG8_1": {"domain": "Atg8/LC3", "pfam": "PF02991", "binding_mode": "hydrophobic_groove"},
    "LIG_ULK_LIR_1": {"domain": "Atg8/LC3", "pfam": "PF02991", "binding_mode": "hydrophobic_groove"},
    "LIG_BH3_BclXL_1": {"domain": "Bcl-2", "pfam": "PF00452", "binding_mode": "BH3_groove"},
    "LIG_IAP_AVPI_1": {"domain": "BIR", "pfam": "PF00653", "binding_mode": "IAP_binding"},
    "CLV_CASP_3": {"domain": "Caspase", "pfam": "PF00656", "binding_mode": "active_site"},
    "CLV_CASP_8": {"domain": "Caspase", "pfam": "PF00656", "binding_mode": "active_site"},
    "DOC_CASP_3": {"domain": "Caspase", "pfam": "PF00656", "binding_mode": "exosite"},
    "LIG_PP1_RVXF_1": {"domain": "PP1c", "pfam": "PF00149", "binding_mode": "RVxF_groove"},
    "LIG_PP1_SILK_1": {"domain": "PP1c", "pfam": "PF00149", "binding_mode": "SILK_site"},
    "LIG_PP1_MyPhoNE_1": {"domain": "PP1c", "pfam": "PF00149", "binding_mode": "MyPhoNE_site"},
    "DOC_PP1_RVXF_1": {"domain": "PP1c", "pfam": "PF00149", "binding_mode": "RVxF_groove"},
    "LIG_PP2A_B56_1": {"domain": "B56/PP2A_B", "pfam": "PF01603", "binding_mode": "LxxIxE_groove"},
    "LIG_CaN_CNBR_1": {"domain": "Calcineurin_B", "pfam": "PF00036", "binding_mode": "PxIxIT_groove"},
    "DOC_PP2B_LxvP_1": {"domain": "Calcineurin_A", "pfam": "PF00149", "binding_mode": "LxVP_pocket"},
    "LIG_14-3-3_CanoR_1": {"domain": "14-3-3", "pfam": "PF00244", "binding_mode": "phospho_groove"},
    "LIG_14-3-3_CanoR_2": {"domain": "14-3-3", "pfam": "PF00244", "binding_mode": "phospho_groove"},
    "LIG_14-3-3_CterR_2": {"domain": "14-3-3", "pfam": "PF00244", "binding_mode": "phospho_groove"},
    "DOC_MAPK_DCC_7": {"domain": "MAPK", "pfam": "PF00069", "binding_mode": "D_groove"},
    "DOC_MAPK_FxFP_6": {"domain": "MAPK", "pfam": "PF00069", "binding_mode": "DEF_site"},
    "DOC_MAPK_JIP1_7": {"domain": "MAPK_JNK", "pfam": "PF00069", "binding_mode": "JIP_site"},
    "MOD_CDK_SPxK_1": {"domain": "Kinase_CDK", "pfam": "PF00069", "binding_mode": "active_site"},
    "MOD_CDK_SPxxK_3": {"domain": "Kinase_CDK", "pfam": "PF00069", "binding_mode": "active_site"},
    "MOD_CK1_1": {"domain": "Kinase_CK1", "pfam": "PF00069", "binding_mode": "active_site"},
    "MOD_CK2_1": {"domain": "Kinase_CK2", "pfam": "PF00069", "binding_mode": "active_site"},
    "MOD_GSK3_1": {"domain": "Kinase_GSK3", "pfam": "PF00069", "binding_mode": "active_site"},
    "MOD_PKA_2": {"domain": "Kinase_PKA", "pfam": "PF00069", "binding_mode": "active_site"},
    "MOD_PKB_1": {"domain": "Kinase_Akt", "pfam": "PF00069", "binding_mode": "active_site"},
    "MOD_PIKK_1": {"domain": "PIKK_kinase", "pfam": "PF00454", "binding_mode": "active_site"},
    "MOD_PLK_1": {"domain": "Kinase_PLK", "pfam": "PF00069", "binding_mode": "active_site"},
    "MOD_NEK2_1": {"domain": "Kinase_NEK", "pfam": "PF00069", "binding_mode": "active_site"},
    "MOD_GRK_1": {"domain": "Kinase_GRK", "pfam": "PF00069", "binding_mode": "active_site"},
    "LIG_TRAF2_1": {"domain": "TRAF_MATH", "pfam": "PF00917", "binding_mode": "TRAF_groove"},
    "LIG_TRAF6_MATH_1": {"domain": "TRAF_MATH", "pfam": "PF00917", "binding_mode": "TRAF_groove"},
    "LIG_NEMO_UBAN_1": {"domain": "UBAN/CoZi", "pfam": "PF11445", "binding_mode": "ubiquitin_groove"},
    "LIG_TBK1_1": {"domain": "TBK1_scaffold", "pfam": "PF12613", "binding_mode": "scaffold_site"},
    "LIG_IRF3_1": {"domain": "IRF_IAD", "pfam": "PF00605", "binding_mode": "IAD_surface"},
    "LIG_SH2_STAT3_1": {"domain": "SH2", "pfam": "PF00017", "binding_mode": "pY_pocket"},
    "LIG_SH2_GRB2_1": {"domain": "SH2", "pfam": "PF00017", "binding_mode": "pY_pocket"},
    "LIG_SH2_SRC_1": {"domain": "SH2", "pfam": "PF00017", "binding_mode": "pY_pocket"},
    "LIG_SH3_1": {"domain": "SH3", "pfam": "PF00018", "binding_mode": "PPII_groove"},
    "LIG_SH3_3": {"domain": "SH3", "pfam": "PF00018", "binding_mode": "PPII_groove"},
    "LIG_SH3_4": {"domain": "SH3", "pfam": "PF00018", "binding_mode": "PPII_groove"},
    "TRG_NLS_MonoExtC_4": {"domain": "Importin_alpha_ARM", "pfam": "PF01749", "binding_mode": "NLS_groove"},
    "TRG_NLS_Bipartite_1": {"domain": "Importin_alpha_ARM", "pfam": "PF01749", "binding_mode": "NLS_groove"},
    "TRG_NES_CRM1_1": {"domain": "XPO1_NES", "pfam": "PF08389", "binding_mode": "NES_groove"},
    "DEG_SCF_FBW7_1": {"domain": "WD40_FBW7", "pfam": "PF00400", "binding_mode": "degron_pocket"},
    "DEG_SCF_TRCP1_1": {"domain": "WD40_TRCP", "pfam": "PF00400", "binding_mode": "degron_pocket"},
    "DEG_APCC_DBOX_1": {"domain": "APC_WD40", "pfam": "PF12860", "binding_mode": "Dbox_receptor"},
    "DEG_APCC_KENBOX_2": {"domain": "APC_WD40", "pfam": "PF12860", "binding_mode": "KEN_receptor"},
    "DEG_MDM2_SWIB_1": {"domain": "SWIB/MDM2", "pfam": "PF02201", "binding_mode": "p53_groove"},
    "LIG_UBQ_UIM_1": {"domain": "Ubiquitin", "pfam": "PF00240", "binding_mode": "Ile44_patch"},
    "LIG_WW_1": {"domain": "WW", "pfam": "PF00397", "binding_mode": "PY_groove"},
    "LIG_WW_2": {"domain": "WW", "pfam": "PF00397", "binding_mode": "PPLP_groove"},
    "LIG_WW_3": {"domain": "WW", "pfam": "PF00397", "binding_mode": "PPxR_groove"},
    "LIG_Dynein_DLC8_1": {"domain": "Dynein_LC", "pfam": "PF01221", "binding_mode": "TQT_groove"},
    "LIG_MYOSIN_MEF2_1": {"domain": "Myosin_head", "pfam": "PF00063", "binding_mode": "cargo_site"},
    "LIG_KLC_1": {"domain": "KLC_TPR", "pfam": "PF13414", "binding_mode": "WD_groove"},
    "LIG_PDZ_Class_1": {"domain": "PDZ", "pfam": "PF00595", "binding_mode": "Cterm_groove"},
    "LIG_PDZ_Class_2": {"domain": "PDZ", "pfam": "PF00595", "binding_mode": "Cterm_groove"},
    "LIG_PDZ_Class_3": {"domain": "PDZ", "pfam": "PF00595", "binding_mode": "Cterm_groove"},
    "TRG_ENDOCYTIC_2": {"domain": "AP2_mu", "pfam": "PF01217", "binding_mode": "YXXF_pocket"},
    "LIG_PTB_Phospho_1": {"domain": "PTB", "pfam": "PF08416", "binding_mode": "NPxY_groove"},
    "TRG_LysEnd_APsAcLL_1": {"domain": "AP_sigma", "pfam": "PF01217", "binding_mode": "LL_pocket"},
    "MOD_SUMO_for_1": {"domain": "Ubc9/UBE2I", "pfam": "PF00179", "binding_mode": "active_site"},
    "MOD_SUMO_rev_2": {"domain": "Ubc9/UBE2I", "pfam": "PF00179", "binding_mode": "active_site"},
    "MOD_N-GLC_1": {"domain": "STT3_OST", "pfam": "PF02516", "binding_mode": "active_site"},
    "LIG_eIF4E_1": {"domain": "eIF4E_cap", "pfam": "PF01652", "binding_mode": "dorsal_surface"},
    "LIG_EIF4G3_1": {"domain": "eIF4G_MIF4G", "pfam": "PF02854", "binding_mode": "HEAT_groove"},
    "DOC_CYCLIN_RxL_1": {"domain": "Cyclin_N", "pfam": "PF00134", "binding_mode": "hydrophobic_patch"},
    "DOC_CYCLIN_D_SAMS_1": {"domain": "Cyclin_N", "pfam": "PF00134", "binding_mode": "LL_site"},
    "DEG_CRL4_CDT2_2": {"domain": "WD40_CDT2", "pfam": "PF00400", "binding_mode": "PIP_degron"},
    "LIG_FHA_1": {"domain": "FHA", "pfam": "PF00498", "binding_mode": "pT_pocket"},
    "LIG_FHA_2": {"domain": "FHA", "pfam": "PF00498", "binding_mode": "pT_pocket"},
    "LIG_BRCT_BRCA1_1": {"domain": "BRCT", "pfam": "PF00533", "binding_mode": "pS_pocket"},
    "LIG_HSP90_1": {"domain": "HSP90_C", "pfam": "PF00183", "binding_mode": "C_term_groove"},
    "DOC_USP7_MATH_1": {"domain": "MATH/TRAF", "pfam": "PF00917", "binding_mode": "MATH_groove"},
    "LIG_PCNA_PIPBox_1": {"domain": "PCNA", "pfam": "PF00705", "binding_mode": "PIP_groove"},
    "LIG_Rb_LxCxE_1": {"domain": "Rb_pocket_B", "pfam": "PF01857", "binding_mode": "LxCxE_cleft"},
    "LIG_Rb_pABgroove_1": {"domain": "Rb_pocket_AB", "pfam": "PF01857", "binding_mode": "AB_groove"},
}


# ═══════════════════════════════════════════════════════════════
# GO BIOLOGICAL PROCESS ANNOTATIONS (curated for key host genes)
# ═══════════════════════════════════════════════════════════════

GO_ANNOTATIONS = {
    "MAP1LC3A": ["GO:0006914|autophagy", "GO:0000045|autophagosome assembly"],
    "MAP1LC3B": ["GO:0006914|autophagy", "GO:0000045|autophagosome assembly"],
    "GABARAP": ["GO:0006914|autophagy", "GO:0016236|macroautophagy"],
    "GABARAPL2": ["GO:0006914|autophagy"],
    "BECN1": ["GO:0006914|autophagy", "GO:0010506|regulation of autophagy"],
    "ULK1": ["GO:0006914|autophagy", "GO:0010507|negative regulation of autophagy"],
    "BCL2": ["GO:0043066|negative regulation of apoptotic process", "GO:0006915|apoptotic process"],
    "BCL2L1": ["GO:0043066|negative regulation of apoptotic process"],
    "MCL1": ["GO:0043066|negative regulation of apoptotic process"],
    "CASP3": ["GO:0006915|apoptotic process", "GO:0097194|execution phase of apoptosis"],
    "CASP7": ["GO:0006915|apoptotic process", "GO:0097194|execution phase of apoptosis"],
    "CASP8": ["GO:0006915|apoptotic process", "GO:0008625|extrinsic apoptotic signaling"],
    "BIRC2": ["GO:0043066|negative regulation of apoptotic process"],
    "BIRC3": ["GO:0043066|negative regulation of apoptotic process"],
    "XIAP": ["GO:0043066|negative regulation of apoptotic process"],
    "PPP1CA": ["GO:0006470|protein dephosphorylation", "GO:0050804|modulation of chemical synaptic transmission"],
    "PPP1CB": ["GO:0006470|protein dephosphorylation", "GO:0050804|modulation of chemical synaptic transmission"],
    "PPP1CC": ["GO:0006470|protein dephosphorylation"],
    "PPP2R5A": ["GO:0006470|protein dephosphorylation", "GO:0030163|protein catabolic process"],
    "PPP2R5C": ["GO:0006470|protein dephosphorylation"],
    "PPP2R5D": ["GO:0006470|protein dephosphorylation", "GO:0007399|nervous system development"],
    "PPP3CA": ["GO:0006470|protein dephosphorylation", "GO:0050804|modulation of chemical synaptic transmission"],
    "PPP3CB": ["GO:0006470|protein dephosphorylation"],
    "YWHAB": ["GO:0007165|signal transduction", "GO:0051220|cytoplasmic sequestering of protein"],
    "YWHAE": ["GO:0007165|signal transduction", "GO:0007399|nervous system development"],
    "YWHAG": ["GO:0007165|signal transduction"],
    "YWHAZ": ["GO:0007165|signal transduction", "GO:0043066|negative regulation of apoptotic process"],
    "MAPK1": ["GO:0000187|activation of MAPK activity", "GO:0007254|JNK cascade"],
    "MAPK3": ["GO:0000187|activation of MAPK activity"],
    "MAPK8": ["GO:0007254|JNK cascade", "GO:0006915|apoptotic process"],
    "MAPK14": ["GO:0006950|response to stress", "GO:0006954|inflammatory response"],
    "CDK1": ["GO:0007049|cell cycle", "GO:0051301|cell division"],
    "CDK2": ["GO:0007049|cell cycle", "GO:0006260|DNA replication"],
    "CDK5": ["GO:0007399|nervous system development", "GO:0007268|chemical synaptic transmission"],
    "GSK3A": ["GO:0016310|phosphorylation", "GO:0060070|canonical Wnt signaling pathway"],
    "GSK3B": ["GO:0016310|phosphorylation", "GO:0060070|canonical Wnt signaling pathway", "GO:0007399|nervous system development"],
    "PRKACA": ["GO:0007165|signal transduction", "GO:0050804|modulation of chemical synaptic transmission"],
    "PRKACB": ["GO:0007165|signal transduction"],
    "AKT1": ["GO:0043066|negative regulation of apoptotic process", "GO:0051896|regulation of protein kinase B signaling"],
    "AKT3": ["GO:0007399|nervous system development", "GO:0043066|negative regulation of apoptotic process"],
    "TRAF2": ["GO:0043123|positive regulation of I-kappaB kinase/NF-kappaB signaling", "GO:0006915|apoptotic process"],
    "TRAF6": ["GO:0043123|positive regulation of I-kappaB kinase/NF-kappaB signaling", "GO:0045087|innate immune response"],
    "IKBKG": ["GO:0043123|positive regulation of I-kappaB kinase/NF-kappaB signaling", "GO:0045087|innate immune response"],
    "TBK1": ["GO:0032481|positive regulation of type I interferon production", "GO:0045087|innate immune response"],
    "IRF3": ["GO:0032481|positive regulation of type I interferon production", "GO:0045087|innate immune response"],
    "IRF7": ["GO:0032481|positive regulation of type I interferon production"],
    "STAT3": ["GO:0007259|receptor signaling pathway via JAK-STAT", "GO:0006915|apoptotic process"],
    "GRB2": ["GO:0007169|transmembrane receptor protein tyrosine kinase signaling pathway"],
    "SRC": ["GO:0007169|transmembrane receptor protein tyrosine kinase signaling pathway"],
    "FYN": ["GO:0007268|chemical synaptic transmission", "GO:0007169|transmembrane receptor protein tyrosine kinase signaling pathway"],
    "KPNA1": ["GO:0006606|protein import into nucleus"],
    "KPNA2": ["GO:0006606|protein import into nucleus"],
    "XPO1": ["GO:0006611|protein export from nucleus"],
    "FBXW7": ["GO:0006511|ubiquitin-dependent protein catabolic process"],
    "MDM2": ["GO:0006511|ubiquitin-dependent protein catabolic process", "GO:0043066|negative regulation of apoptotic process"],
    "NEDD4": ["GO:0006511|ubiquitin-dependent protein catabolic process", "GO:0007399|nervous system development"],
    "DYNLL1": ["GO:0007018|microtubule-based movement", "GO:0008088|axo-dendritic transport"],
    "DYNLL2": ["GO:0007018|microtubule-based movement"],
    "KLC1": ["GO:0007018|microtubule-based movement", "GO:0008088|axo-dendritic transport"],
    "KLC2": ["GO:0007018|microtubule-based movement"],
    "MYO5A": ["GO:0007018|microtubule-based movement", "GO:0007268|chemical synaptic transmission"],
    "DLG1": ["GO:0007268|chemical synaptic transmission", "GO:0007399|nervous system development"],
    "DLG4": ["GO:0007268|chemical synaptic transmission", "GO:0050804|modulation of chemical synaptic transmission"],
    "AP2M1": ["GO:0006897|endocytosis", "GO:0048268|clathrin coat assembly"],
    "EIF4E": ["GO:0006413|translational initiation"],
    "EIF4G1": ["GO:0006413|translational initiation"],
    "CSNK2A1": ["GO:0016310|phosphorylation", "GO:0060070|canonical Wnt signaling pathway"],
    "CSNK1A1": ["GO:0016310|phosphorylation", "GO:0060070|canonical Wnt signaling pathway"],
    "UBE2I": ["GO:0016925|protein sumoylation"],
    "HSP90AA1": ["GO:0006457|protein folding", "GO:0051087|chaperone binding"],
    "RB1": ["GO:0007049|cell cycle", "GO:0045786|negative regulation of cell cycle"],
    "PCNA": ["GO:0006260|DNA replication", "GO:0006281|DNA repair"],
    "ATM": ["GO:0006974|cellular response to DNA damage stimulus", "GO:0006281|DNA repair"],
}

# CNS-relevant GO terms (for scoring)
CNS_GO_TERMS = {
    "GO:0007399": ("nervous system development", 3.0),
    "GO:0007268": ("chemical synaptic transmission", 3.0),
    "GO:0050804": ("modulation of chemical synaptic transmission", 3.0),
    "GO:0008088": ("axo-dendritic transport", 3.0),
    "GO:0006914": ("autophagy", 2.5),
    "GO:0000045": ("autophagosome assembly", 2.5),
    "GO:0016236": ("macroautophagy", 2.5),
    "GO:0010506": ("regulation of autophagy", 2.5),
    "GO:0006915": ("apoptotic process", 2.0),
    "GO:0043066": ("negative regulation of apoptotic process", 2.0),
    "GO:0097194": ("execution phase of apoptosis", 2.0),
    "GO:0045087": ("innate immune response", 2.0),
    "GO:0032481": ("positive regulation of type I interferon production", 2.0),
    "GO:0043123": ("positive regulation of I-kappaB kinase/NF-kappaB signaling", 2.0),
    "GO:0006470": ("protein dephosphorylation", 1.5),
    "GO:0007165": ("signal transduction", 1.0),
}

# Reactome pathway annotations
REACTOME_PATHWAYS = {
    "autophagy": {"id": "R-HSA-9612973", "genes": ["MAP1LC3A", "MAP1LC3B", "GABARAP", "GABARAPL2", "BECN1", "ULK1"]},
    "apoptosis": {"id": "R-HSA-109581", "genes": ["BCL2", "BCL2L1", "MCL1", "CASP3", "CASP7", "CASP8", "BIRC2", "BIRC3", "XIAP"]},
    "innate_immunity": {"id": "R-HSA-168249", "genes": ["TRAF2", "TRAF6", "IKBKG", "TBK1", "IRF3", "IRF7"]},
    "interferon_signaling": {"id": "R-HSA-913531", "genes": ["TBK1", "IRF3", "IRF7", "STAT3"]},
    "axonal_transport": {"id": "R-HSA-2132295", "genes": ["DYNLL1", "DYNLL2", "KLC1", "KLC2", "MYO5A"]},
    "synaptic_transmission": {"id": "R-HSA-112315", "genes": ["DLG1", "DLG4", "FYN", "CDK5", "PPP1CA", "PPP1CB", "PRKACA"]},
    "MAPK_cascade": {"id": "R-HSA-5684996", "genes": ["MAPK1", "MAPK3", "MAPK8", "MAPK14", "GRB2"]},
    "PI3K_Akt": {"id": "R-HSA-109581", "genes": ["AKT1", "AKT3"]},
    "Wnt_signaling": {"id": "R-HSA-195721", "genes": ["GSK3A", "GSK3B", "CSNK1A1"]},
    "cell_cycle": {"id": "R-HSA-1640170", "genes": ["CDK1", "CDK2", "RB1", "FBXW7"]},
    "translation": {"id": "R-HSA-72613", "genes": ["EIF4E", "EIF4G1"]},
    "ubiquitin_system": {"id": "R-HSA-461359", "genes": ["FBXW7", "MDM2", "NEDD4", "UBE2I", "USP7"]},
    "protein_folding": {"id": "R-HSA-391251", "genes": ["HSP90AA1"]},
    "nuclear_transport": {"id": "R-HSA-159236", "genes": ["KPNA1", "KPNA2", "XPO1"]},
    "endocytosis": {"id": "R-HSA-8856828", "genes": ["AP2M1"]},
    "DNA_repair": {"id": "R-HSA-73894", "genes": ["ATM", "PCNA", "CHEK2"]},
}

# CNS expression (brain TPM from GTEx v8)
CNS_EXPRESSION = {
    "MAP1LC3A": 45, "MAP1LC3B": 120, "GABARAP": 180, "GABARAPL2": 95,
    "BECN1": 35, "ULK1": 15, "ULK2": 10,
    "BCL2": 15, "BCL2L1": 55, "MCL1": 30,
    "PPP1CA": 85, "PPP1CB": 110, "PPP1CC": 75,
    "PPP2R5A": 25, "PPP2R5B": 12, "PPP2R5C": 40, "PPP2R5D": 55, "PPP2R5E": 8,
    "PPP3CA": 60, "PPP3CB": 35,
    "YWHAB": 200, "YWHAE": 280, "YWHAG": 350, "YWHAH": 150, "YWHAQ": 120, "YWHAZ": 250,
    "MAPK1": 40, "MAPK3": 25, "MAPK8": 18, "MAPK9": 12, "MAPK10": 25, "MAPK14": 15,
    "CDK1": 3, "CDK2": 5, "CDK5": 65,
    "CSNK1A1": 35, "CSNK1D": 20, "CSNK1E": 15,
    "CSNK2A1": 30, "CSNK2A2": 8,
    "GSK3A": 20, "GSK3B": 45,
    "PRKACA": 60, "PRKACB": 40,
    "AKT1": 20, "AKT2": 8, "AKT3": 35,
    "TRAF2": 12, "TRAF6": 8, "IKBKG": 10,
    "TBK1": 10, "IRF3": 8, "IRF7": 3,
    "CASP3": 10, "CASP7": 5, "CASP8": 4,
    "BIRC2": 8, "BIRC3": 5, "XIAP": 15,
    "STAT3": 15, "GRB2": 20, "SRC": 15, "FYN": 55, "LCK": 2, "ABL1": 10,
    "NCK1": 12, "NCK2": 8,
    "KPNA1": 20, "KPNA2": 15, "KPNA3": 10, "XPO1": 25,
    "FBXW7": 8, "BTRC": 10, "FBXW11": 8,
    "ANAPC2": 10, "FZR1": 15, "MDM2": 12, "DTL": 3,
    "NEDD4": 35, "NEDD4L": 25, "ITCH": 15, "WWP1": 12,
    "DYNLL1": 90, "DYNLL2": 25, "KLC1": 30, "KLC2": 20,
    "MYO5A": 15, "MYO6": 20,
    "DLG1": 20, "DLG4": 80, "MAGI1": 25, "SCRIB": 15, "NHERF1": 10,
    "AP2M1": 35, "DAB1": 30, "DAB2": 5, "AP1M1": 15, "AP3D1": 10,
    "UBE2I": 25, "STT3A": 20, "STT3B": 15,
    "EIF4E": 30, "EIF4G1": 25,
    "CCNA2": 2, "CCNB1": 2, "CCNE1": 1, "CCND1": 5, "CCND3": 8,
    "CHEK2": 3, "RNF8": 5, "MDC1": 4, "BRCA1": 2,
    "HSP90AA1": 40, "HSP90AB1": 60,
    "RB1": 15, "RBL1": 5, "RBL2": 8,
    "PCNA": 8, "ATM": 5, "ATR": 4, "PLK1": 2, "NEK2": 1,
    "GRK2": 15, "GRK5": 12,
    "USP7": 20, "UBB": 200, "UBC": 150,
}


# ═══════════════════════════════════════════════════════════════
# 1. SLiM STRUCTURAL CONTEXT ANALYSIS
# ═══════════════════════════════════════════════════════════════

def analyze_structural_context(hits_df, disorder_data):
    """Compute flanking residue properties for each SLiM hit."""
    results = []
    for _, row in hits_df.iterrows():
        protein = row["protein"]
        pdata = disorder_data.get(protein)
        if pdata is None:
            continue

        per_res_path = os.path.join(PHASE1_DIR, "disorder_predictions", f"{protein}_per_residue.csv")
        if not os.path.exists(per_res_path):
            continue
        per_res = pd.read_csv(per_res_path)

        start = int(row["start"]) - 1  # 0-indexed
        end = int(row["end"])
        flank = 10
        flank_start = max(0, start - flank)
        flank_end = min(len(per_res), end + flank)

        motif_seq = str(row["matched_sequence"])
        flank_seq = "".join(per_res.iloc[flank_start:flank_end]["residue"].values)

        hydro_motif = np.mean([HYDROPHOBICITY.get(aa, 0) for aa in motif_seq])
        hydro_flank = np.mean([HYDROPHOBICITY.get(aa, 0) for aa in flank_seq]) if flank_seq else 0

        charge_motif = sum(CHARGE.get(aa, 0) for aa in motif_seq)
        charge_flank = sum(CHARGE.get(aa, 0) for aa in flank_seq)

        helix_motif = np.mean([HELIX_PROPENSITY.get(aa, 1.0) for aa in motif_seq])
        helix_flank = np.mean([HELIX_PROPENSITY.get(aa, 1.0) for aa in flank_seq]) if flank_seq else 1.0

        flank_disorder = per_res.iloc[flank_start:flank_end]["iupred3"].mean()
        motif_disorder = per_res.iloc[start:end]["iupred3"].mean()

        n_pro = motif_seq.count("P") + flank_seq.count("P")
        frac_charged = sum(1 for aa in flank_seq if aa in "DEKRH") / max(len(flank_seq), 1)

        accessibility_score = motif_disorder * (1 - abs(hydro_motif) / 4.5)

        results.append({
            "protein": protein,
            "motif_class": row["motif_class"],
            "start": row["start"],
            "end": row["end"],
            "matched_sequence": motif_seq,
            "hydrophobicity_motif": round(hydro_motif, 2),
            "hydrophobicity_flank": round(hydro_flank, 2),
            "net_charge_motif": round(charge_motif, 1),
            "net_charge_flank": round(charge_flank, 1),
            "helix_propensity_motif": round(helix_motif, 2),
            "helix_propensity_flank": round(helix_flank, 2),
            "disorder_motif": round(motif_disorder, 3),
            "disorder_flank": round(flank_disorder, 3),
            "n_prolines_region": n_pro,
            "frac_charged_flank": round(frac_charged, 3),
            "accessibility_score": round(accessibility_score, 3),
        })
    return pd.DataFrame(results)


# ═══════════════════════════════════════════════════════════════
# 2. CNS-WEIGHTED INTERACTOME CONSTRUCTION
# ═══════════════════════════════════════════════════════════════

def build_interactome(hits_df, disorder_data):
    """Build bipartite network: viral proteins ↔ host targets, weighted by CNS relevance."""
    G = nx.Graph()

    for _, row in hits_df.iterrows():
        viral_protein = row["protein"]
        category = row["category"]
        targets_str = str(row.get("targets", ""))
        targets = [t.strip() for t in targets_str.split(";") if t.strip()]
        motif_class = row["motif_class"]
        composite = float(row.get("composite_score", 0))
        pathway = row.get("pathway", "")

        if not G.has_node(viral_protein):
            pdata = disorder_data.get(viral_protein, {})
            G.add_node(viral_protein, bipartite=0, node_type="viral",
                       category=category,
                       idr_fraction=pdata.get("idr_fraction", 0),
                       length=pdata.get("length", 0))

        for host_gene in targets:
            if not host_gene:
                continue
            brain_tpm = CNS_EXPRESSION.get(host_gene, 0)
            cns_specificity = np.log1p(brain_tpm) / np.log1p(350)

            if not G.has_node(host_gene):
                go_terms = GO_ANNOTATIONS.get(host_gene, [])
                cns_go_score = 0
                for term_str in go_terms:
                    go_id = term_str.split("|")[0]
                    if go_id in CNS_GO_TERMS:
                        cns_go_score = max(cns_go_score, CNS_GO_TERMS[go_id][1])

                G.add_node(host_gene, bipartite=1, node_type="host",
                           brain_tpm=brain_tpm,
                           cns_specificity=round(cns_specificity, 3),
                           cns_go_score=cns_go_score,
                           go_terms="; ".join(go_terms))

            edge_key = (viral_protein, host_gene)
            if G.has_edge(*edge_key):
                edata = G.edges[edge_key]
                edata["n_slims"] += 1
                edata["max_composite"] = max(edata["max_composite"], composite)
                edata["motif_classes"].append(motif_class)
                edata["pathways"].add(pathway)
            else:
                cns_edge_weight = composite * cns_specificity
                domain_info = DOMAIN_MOTIF_MAP.get(motif_class, {})
                G.add_edge(viral_protein, host_gene,
                           n_slims=1,
                           max_composite=composite,
                           cns_edge_weight=round(cns_edge_weight, 4),
                           motif_classes=[motif_class],
                           pathways={pathway},
                           domain=domain_info.get("domain", "unknown"),
                           binding_mode=domain_info.get("binding_mode", "unknown"))

    return G


def compute_network_metrics(G):
    """Compute centrality and topology metrics."""
    metrics = {}

    degree_cent = nx.degree_centrality(G)
    between_cent = nx.betweenness_centrality(G)

    viral_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "viral"]
    host_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "host"]

    for node in viral_nodes:
        neighbors = list(G.neighbors(node))
        host_neighbors = [n for n in neighbors if G.nodes[n].get("node_type") == "host"]
        cns_weighted_degree = sum(
            G.edges[node, h].get("cns_edge_weight", 0) for h in host_neighbors
        )
        mean_host_cns = np.mean([
            G.nodes[h].get("cns_specificity", 0) for h in host_neighbors
        ]) if host_neighbors else 0

        pathways_targeted = set()
        for h in host_neighbors:
            pathways_targeted.update(G.edges[node, h].get("pathways", set()))

        metrics[node] = {
            "type": "viral",
            "category": G.nodes[node].get("category", ""),
            "degree": G.degree(node),
            "degree_centrality": round(degree_cent[node], 4),
            "betweenness_centrality": round(between_cent[node], 4),
            "cns_weighted_degree": round(cns_weighted_degree, 3),
            "mean_host_cns_specificity": round(mean_host_cns, 3),
            "n_pathways": len(pathways_targeted),
            "pathways": sorted(pathways_targeted),
        }

    for node in host_nodes:
        viral_neighbors = [n for n in G.neighbors(node) if G.nodes[n].get("node_type") == "viral"]
        metrics[node] = {
            "type": "host",
            "brain_tpm": G.nodes[node].get("brain_tpm", 0),
            "cns_specificity": G.nodes[node].get("cns_specificity", 0),
            "degree": G.degree(node),
            "degree_centrality": round(degree_cent[node], 4),
            "betweenness_centrality": round(between_cent[node], 4),
            "viral_partners": viral_neighbors,
        }

    return metrics


# ═══════════════════════════════════════════════════════════════
# 3. PATHWAY ENRICHMENT
# ═══════════════════════════════════════════════════════════════

def pathway_enrichment_analysis(G, disorder_data):
    """Pathway enrichment per viral protein using GO and Reactome."""
    results = {}
    viral_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "viral"]

    all_host_genes = set(n for n, d in G.nodes(data=True) if d.get("node_type") == "host")

    for viral in viral_nodes:
        host_targets = set(n for n in G.neighbors(viral) if G.nodes[n].get("node_type") == "host")

        go_counts = defaultdict(int)
        for gene in host_targets:
            for term_str in GO_ANNOTATIONS.get(gene, []):
                go_id, go_name = term_str.split("|")
                go_counts[(go_id, go_name)] += 1

        bg_go_counts = defaultdict(int)
        for gene in all_host_genes:
            for term_str in GO_ANNOTATIONS.get(gene, []):
                go_id, go_name = term_str.split("|")
                bg_go_counts[(go_id, go_name)] += 1

        go_enrichment = []
        for (go_id, go_name), count in go_counts.items():
            bg_count = bg_go_counts.get((go_id, go_name), 0)
            if bg_count < 2:
                continue
            table = [
                [count, bg_count - count],
                [len(host_targets) - count, len(all_host_genes) - len(host_targets) - (bg_count - count)]
            ]
            table = [[max(0, c) for c in row] for row in table]
            try:
                odds, pval = stats.fisher_exact(table, alternative="greater")
            except ValueError:
                continue
            go_enrichment.append({
                "go_id": go_id,
                "go_name": go_name,
                "count": count,
                "total_targets": len(host_targets),
                "bg_count": bg_count,
                "bg_total": len(all_host_genes),
                "odds_ratio": round(odds, 2),
                "p_value": round(pval, 4),
            })
        go_enrichment.sort(key=lambda x: x["p_value"])

        reactome_enrichment = []
        for pw_name, pw_data in REACTOME_PATHWAYS.items():
            pw_genes = set(pw_data["genes"])
            overlap = host_targets & pw_genes
            if not overlap:
                continue
            bg_overlap = all_host_genes & pw_genes
            table = [
                [len(overlap), len(bg_overlap) - len(overlap)],
                [len(host_targets) - len(overlap), len(all_host_genes) - len(host_targets) - (len(bg_overlap) - len(overlap))]
            ]
            table = [[max(0, c) for c in row] for row in table]
            try:
                odds, pval = stats.fisher_exact(table, alternative="greater")
            except ValueError:
                continue
            reactome_enrichment.append({
                "pathway": pw_name,
                "reactome_id": pw_data["id"],
                "overlap_genes": sorted(overlap),
                "count": len(overlap),
                "pathway_size": len(pw_genes),
                "odds_ratio": round(odds, 2),
                "p_value": round(pval, 4),
            })
        reactome_enrichment.sort(key=lambda x: x["p_value"])

        results[viral] = {
            "n_host_targets": len(host_targets),
            "go_enrichment": go_enrichment[:15],
            "reactome_enrichment": reactome_enrichment[:10],
        }

    return results


# ═══════════════════════════════════════════════════════════════
# 4. CORE HYPOTHESIS TEST
# ═══════════════════════════════════════════════════════════════

def hypothesis_test(network_metrics, disorder_data):
    """
    H1: Neurovirulence proteins have higher CNS-targeting scores than controls.
    CNS-targeting score = SLiM density * CNS interaction fraction
    """
    target_scores = []
    control_scores = []

    for protein, m in network_metrics.items():
        if m["type"] != "viral":
            continue
        pdata = disorder_data.get(protein, {})
        idr_aa = sum(e - s + 1 for s, e in [tuple(b) for b in pdata.get("idr_boundaries", [])])
        if idr_aa == 0:
            continue

        slim_density = m["degree"] / idr_aa * 100
        cns_interaction_frac = m["mean_host_cns_specificity"]
        cns_targeting_score = slim_density * cns_interaction_frac
        pathway_diversity = m["n_pathways"]
        cns_weighted_degree = m["cns_weighted_degree"]

        composite_h1 = (
            0.4 * (cns_weighted_degree / 10) +
            0.3 * (cns_targeting_score / 20) +
            0.3 * (pathway_diversity / 15)
        )

        entry = {
            "protein": protein,
            "category": m["category"],
            "degree": m["degree"],
            "idr_aa": idr_aa,
            "slim_density": round(slim_density, 2),
            "cns_interaction_frac": round(cns_interaction_frac, 3),
            "cns_targeting_score": round(cns_targeting_score, 3),
            "cns_weighted_degree": round(cns_weighted_degree, 3),
            "pathway_diversity": pathway_diversity,
            "composite_h1_score": round(composite_h1, 4),
        }

        if m["category"] == "target":
            target_scores.append(entry)
        else:
            control_scores.append(entry)

    t_vals = [e["composite_h1_score"] for e in target_scores]
    c_vals = [e["composite_h1_score"] for e in control_scores]

    if len(t_vals) >= 2 and len(c_vals) >= 2:
        u_stat, u_pval = stats.mannwhitneyu(t_vals, c_vals, alternative="greater")
    else:
        u_stat, u_pval = float("nan"), float("nan")

    effect_size = (np.mean(t_vals) - np.mean(c_vals)) / np.std(t_vals + c_vals) if np.std(t_vals + c_vals) > 0 else 0

    return {
        "targets": target_scores,
        "controls": control_scores,
        "mann_whitney_U": round(u_stat, 2) if not np.isnan(u_stat) else None,
        "mann_whitney_p": round(u_pval, 4) if not np.isnan(u_pval) else None,
        "effect_size_cohens_d": round(effect_size, 3),
        "target_mean": round(np.mean(t_vals), 4),
        "control_mean": round(np.mean(c_vals), 4),
        "target_std": round(np.std(t_vals), 4),
        "control_std": round(np.std(c_vals), 4),
    }


# ═══════════════════════════════════════════════════════════════
# 5. PLOTTING
# ═══════════════════════════════════════════════════════════════

def plot_network(G, output_path):
    """Network visualization colored by node type and CNS relevance."""
    fig, ax = plt.subplots(1, 1, figsize=(18, 14))

    viral_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "viral"]
    host_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "host"]

    pos = nx.spring_layout(G, k=2.5, iterations=100, seed=42)

    target_viral = [n for n in viral_nodes if G.nodes[n].get("category") == "target"]
    control_viral = [n for n in viral_nodes if G.nodes[n].get("category") == "control"]

    host_sizes = [max(100, G.nodes[n].get("brain_tpm", 10) * 3) for n in host_nodes]
    host_colors = [G.nodes[n].get("cns_specificity", 0) for n in host_nodes]

    nx.draw_networkx_nodes(G, pos, nodelist=host_nodes, node_size=host_sizes,
                           node_color=host_colors, cmap=plt.cm.YlOrRd,
                           vmin=0, vmax=1, alpha=0.7, ax=ax)

    nx.draw_networkx_nodes(G, pos, nodelist=target_viral, node_size=800,
                           node_color="#e74c3c", node_shape="s", alpha=0.9, ax=ax,
                           edgecolors="black", linewidths=2)
    nx.draw_networkx_nodes(G, pos, nodelist=control_viral, node_size=500,
                           node_color="#3498db", node_shape="s", alpha=0.7, ax=ax,
                           edgecolors="black", linewidths=1)

    edge_weights = [G.edges[e].get("cns_edge_weight", 0.01) for e in G.edges()]
    max_w = max(edge_weights) if edge_weights else 1
    edge_widths = [max(0.3, w / max_w * 3) for w in edge_weights]

    target_edges = [(u, v) for u, v in G.edges() if G.nodes[u].get("category") == "target" or G.nodes[v].get("category") == "target"]
    control_edges = [(u, v) for u, v in G.edges() if (u, v) not in target_edges]

    nx.draw_networkx_edges(G, pos, edgelist=target_edges,
                           width=[max(0.3, G.edges[e].get("cns_edge_weight", 0.01) / max_w * 3) for e in target_edges],
                           alpha=0.6, edge_color="#e74c3c", ax=ax)
    nx.draw_networkx_edges(G, pos, edgelist=control_edges,
                           width=[max(0.2, G.edges[e].get("cns_edge_weight", 0.01) / max_w * 2) for e in control_edges],
                           alpha=0.3, edge_color="#95a5a6", ax=ax)

    labels = {}
    for n in viral_nodes:
        labels[n] = n
    for n in host_nodes:
        if G.degree(n) >= 3 or G.nodes[n].get("brain_tpm", 0) > 50:
            labels[n] = n

    nx.draw_networkx_labels(G, pos, labels, font_size=7, font_weight="bold", ax=ax)

    ax.set_title("HSV-1 SLiM-Mediated CNS Interactome", fontsize=16, fontweight="bold")
    legend_elements = [
        mpatches.Patch(color="#e74c3c", label="Neurovirulence targets"),
        mpatches.Patch(color="#3498db", label="Control proteins"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#ff6b35",
                   markersize=10, label="High CNS expression"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#ffd700",
                   markersize=7, label="Moderate CNS expression"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=9)
    ax.axis("off")

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


def plot_hypothesis_test(h1_results, output_path):
    """Publication figure for the core hypothesis test."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Panel A: CNS-targeting score boxplot
    ax_a = axes[0]
    target_scores = [e["composite_h1_score"] for e in h1_results["targets"]]
    control_scores = [e["composite_h1_score"] for e in h1_results["controls"]]

    bp = ax_a.boxplot([target_scores, control_scores],
                      patch_artist=True, widths=0.5)
    ax_a.set_xticks([1, 2])
    ax_a.set_xticklabels(["Neurovirulence\n(n=3)", "Control\n(n=8)"])
    bp["boxes"][0].set_facecolor("#e74c3c")
    bp["boxes"][0].set_alpha(0.5)
    bp["boxes"][1].set_facecolor("#3498db")
    bp["boxes"][1].set_alpha(0.5)

    for entry in h1_results["targets"]:
        ax_a.plot(1, entry["composite_h1_score"], "o", color="#e74c3c", markersize=10, zorder=5)
        ax_a.annotate(entry["protein"], (1.15, entry["composite_h1_score"]), fontsize=8)
    for entry in h1_results["controls"]:
        ax_a.plot(2, entry["composite_h1_score"], "o", color="#3498db", markersize=8, zorder=5)
        ax_a.annotate(entry["protein"], (2.15, entry["composite_h1_score"]), fontsize=7)

    p_val = h1_results["mann_whitney_p"]
    p_str = f"p = {p_val:.3f}" if p_val is not None else "p = N/A"
    ax_a.set_title(f"A. CNS-Targeting Score\n(Mann-Whitney U, {p_str})", fontsize=11, fontweight="bold")
    ax_a.set_ylabel("Composite H1 Score", fontsize=10)

    # Panel B: Component breakdown
    ax_b = axes[1]
    all_proteins = h1_results["targets"] + h1_results["controls"]
    all_proteins.sort(key=lambda x: -x["composite_h1_score"])
    names = [p["protein"] for p in all_proteins]
    cns_wd = [p["cns_weighted_degree"] / 10 * 0.4 for p in all_proteins]
    cns_ts = [p["cns_targeting_score"] / 20 * 0.3 for p in all_proteins]
    pw_div = [p["pathway_diversity"] / 15 * 0.3 for p in all_proteins]
    colors = ["#e74c3c" if p["category"] == "target" else "#3498db" for p in all_proteins]

    y = np.arange(len(names))
    ax_b.barh(y, cns_wd, color=[c for c in colors], alpha=0.8, label="CNS-weighted degree")
    ax_b.barh(y, cns_ts, left=cns_wd, color=[c for c in colors], alpha=0.5, label="CNS targeting")
    ax_b.barh(y, pw_div, left=[a + b for a, b in zip(cns_wd, cns_ts)],
              color=[c for c in colors], alpha=0.3, label="Pathway diversity")
    ax_b.set_yticks(y)
    ax_b.set_yticklabels(names, fontsize=9)
    ax_b.invert_yaxis()
    ax_b.set_xlabel("Score contribution", fontsize=10)
    ax_b.set_title("B. Score Decomposition", fontsize=11, fontweight="bold")
    ax_b.legend(fontsize=7, loc="lower right")

    # Panel C: Pathway diversity comparison
    ax_c = axes[2]
    target_pw = [e["pathway_diversity"] for e in h1_results["targets"]]
    control_pw = [e["pathway_diversity"] for e in h1_results["controls"]]
    bp2 = ax_c.boxplot([target_pw, control_pw],
                       patch_artist=True, widths=0.5)
    ax_c.set_xticks([1, 2])
    ax_c.set_xticklabels(["Neurovirulence", "Control"])
    bp2["boxes"][0].set_facecolor("#e74c3c")
    bp2["boxes"][0].set_alpha(0.5)
    bp2["boxes"][1].set_facecolor("#3498db")
    bp2["boxes"][1].set_alpha(0.5)
    for entry in h1_results["targets"]:
        ax_c.plot(1, entry["pathway_diversity"], "o", color="#e74c3c", markersize=10, zorder=5)
    for entry in h1_results["controls"]:
        ax_c.plot(2, entry["pathway_diversity"], "o", color="#3498db", markersize=8, zorder=5)
    ax_c.set_ylabel("Number of pathways targeted", fontsize=10)
    ax_c.set_title("C. Pathway Diversity", fontsize=11, fontweight="bold")

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


def plot_pathway_enrichment(enrichment_results, output_path):
    """Pathway enrichment heatmap for target proteins."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 8))
    targets = ["ICP34.5", "US3", "pUL37"]

    for idx, protein in enumerate(targets):
        ax = axes[idx]
        if protein not in enrichment_results:
            ax.set_title(f"{protein}: No data")
            continue

        react = enrichment_results[protein]["reactome_enrichment"]
        if not react:
            go = enrichment_results[protein]["go_enrichment"][:10]
            names = [g["go_name"][:35] for g in go]
            pvals = [-np.log10(max(g["p_value"], 1e-4)) for g in go]
            counts = [g["count"] for g in go]
            colors = ["#e74c3c" if g["p_value"] < 0.05 else "#f39c12" if g["p_value"] < 0.1 else "#95a5a6" for g in go]
        else:
            names = [r["pathway"].replace("_", " ")[:30] for r in react]
            pvals = [-np.log10(max(r["p_value"], 1e-4)) for r in react]
            counts = [r["count"] for r in react]
            colors = ["#e74c3c" if r["p_value"] < 0.05 else "#f39c12" if r["p_value"] < 0.1 else "#95a5a6" for r in react]

        if names:
            y = range(len(names))
            bars = ax.barh(y, pvals, color=colors, alpha=0.7)
            ax.set_yticks(y)
            ax.set_yticklabels(names, fontsize=8)
            ax.invert_yaxis()
            ax.axvline(-np.log10(0.05), color="red", linestyle="--", alpha=0.5, linewidth=0.8)
            for i, (p, c) in enumerate(zip(pvals, counts)):
                ax.text(p + 0.05, i, f"n={c}", fontsize=7, va="center")

        ax.set_xlabel("-log10(p-value)", fontsize=10)
        ax.set_title(f"{protein} Pathway Enrichment", fontsize=11, fontweight="bold")

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


# ═══════════════════════════════════════════════════════════════
# MAIN PIPELINE
# ═══════════════════════════════════════════════════════════════

def main():
    print("=" * 80)
    print("PHASE 3 — SLiM Structural Context, CNS Interactome & Hypothesis Test")
    print("=" * 80)

    # Load data
    disorder_path = os.path.join(PHASE1_DIR, "results", "disorder_summary_v3.json")
    with open(disorder_path) as f:
        disorder_data = {d["protein"]: d for d in json.load(f)}

    hits_path = os.path.join(PHASE2_DIR, "results", "slim_hits_refined.csv")
    hits_df = pd.read_csv(hits_path)
    print(f"\nLoaded {len(hits_df)} SLiM hits from Phase 2")

    # ── Step 1: Structural context ──
    print(f"\n[1] Analyzing SLiM structural context...")
    struct_df = analyze_structural_context(hits_df, disorder_data)
    struct_csv = os.path.join(RESULTS_DIR, "slim_structural_context.csv")
    struct_df.to_csv(struct_csv, index=False)
    print(f"  Computed flanking properties for {len(struct_df)} SLiMs")

    target_struct = struct_df[struct_df["protein"].isin(["ICP34.5", "US3", "pUL37"])]
    control_struct = struct_df[~struct_df["protein"].isin(["ICP34.5", "US3", "pUL37"])]
    print(f"  Target mean accessibility: {target_struct['accessibility_score'].mean():.3f}")
    print(f"  Control mean accessibility: {control_struct['accessibility_score'].mean():.3f}")
    print(f"  Target mean disorder context: {target_struct['disorder_motif'].mean():.3f}")
    print(f"  Control mean disorder context: {control_struct['disorder_motif'].mean():.3f}")

    # ── Step 2: Build interactome ──
    print(f"\n[2] Building CNS-weighted interactome network...")
    G = build_interactome(hits_df, disorder_data)
    viral_count = sum(1 for _, d in G.nodes(data=True) if d.get("node_type") == "viral")
    host_count = sum(1 for _, d in G.nodes(data=True) if d.get("node_type") == "host")
    print(f"  Network: {G.number_of_nodes()} nodes ({viral_count} viral, {host_count} host), "
          f"{G.number_of_edges()} edges")

    # ── Step 3: Network metrics ──
    print(f"\n[3] Computing network topology metrics...")
    net_metrics = compute_network_metrics(G)

    print(f"\n  Viral protein metrics:")
    print(f"  {'Protein':<15} {'Cat':<8} {'Deg':>4} {'CNS-WD':>8} {'MeanCNS':>8} {'#PW':>4} {'Betw':>8}")
    print("  " + "-" * 60)
    for protein in sorted(net_metrics.keys()):
        m = net_metrics[protein]
        if m["type"] != "viral":
            continue
        print(f"  {protein:<15} {m['category']:<8} {m['degree']:>4} "
              f"{m['cns_weighted_degree']:>8.3f} {m['mean_host_cns_specificity']:>8.3f} "
              f"{m['n_pathways']:>4} {m['betweenness_centrality']:>8.4f}")

    # ── Step 4: Pathway enrichment ──
    print(f"\n[4] Running pathway enrichment analysis...")
    enrichment = pathway_enrichment_analysis(G, disorder_data)
    for protein, data in enrichment.items():
        if protein in ["ICP34.5", "US3", "pUL37"]:
            print(f"\n  {protein} ({data['n_host_targets']} host targets):")
            for r in data["reactome_enrichment"][:5]:
                sig = "*" if r["p_value"] < 0.05 else ""
                print(f"    {r['pathway']:<25} p={r['p_value']:.4f}{sig} "
                      f"({r['count']}/{r['pathway_size']} genes: {', '.join(r['overlap_genes'])})")

    # ── Step 5: Core hypothesis test ──
    print(f"\n[5] Running core hypothesis test (H1)...")
    h1_results = hypothesis_test(net_metrics, disorder_data)
    print(f"\n  H1: Neurovirulence proteins have higher CNS-targeting scores")
    print(f"  Target mean:  {h1_results['target_mean']:.4f} +/- {h1_results['target_std']:.4f}")
    print(f"  Control mean: {h1_results['control_mean']:.4f} +/- {h1_results['control_std']:.4f}")
    print(f"  Effect size (Cohen's d): {h1_results['effect_size_cohens_d']:.3f}")
    print(f"  Mann-Whitney U: {h1_results['mann_whitney_U']}")
    print(f"  p-value: {h1_results['mann_whitney_p']}")

    print(f"\n  Per-protein scores:")
    all_scores = h1_results["targets"] + h1_results["controls"]
    all_scores.sort(key=lambda x: -x["composite_h1_score"])
    print(f"  {'Protein':<15} {'Cat':<8} {'Deg':>4} {'Density':>8} {'CNS-frac':>8} "
          f"{'CNS-WD':>8} {'#PW':>4} {'H1score':>8}")
    print("  " + "-" * 70)
    for e in all_scores:
        print(f"  {e['protein']:<15} {e['category']:<8} {e['degree']:>4} "
              f"{e['slim_density']:>8.2f} {e['cns_interaction_frac']:>8.3f} "
              f"{e['cns_weighted_degree']:>8.3f} {e['pathway_diversity']:>4} "
              f"{e['composite_h1_score']:>8.4f}")

    # ── Step 6: Generate plots ──
    print(f"\n[6] Generating publication plots...")
    network_path = os.path.join(PLOT_DIR, "cns_interactome_network.png")
    plot_network(G, network_path)
    print(f"  {network_path}")

    h1_path = os.path.join(PLOT_DIR, "hypothesis_test_h1.png")
    plot_hypothesis_test(h1_results, h1_path)
    print(f"  {h1_path}")

    pw_path = os.path.join(PLOT_DIR, "pathway_enrichment_targets.png")
    plot_pathway_enrichment(enrichment, pw_path)
    print(f"  {pw_path}")

    # ── Step 7: Save results ──
    print(f"\n[7] Saving results...")

    # Network edge list
    edges_data = []
    for u, v, d in G.edges(data=True):
        edges_data.append({
            "viral_protein": u if G.nodes[u].get("node_type") == "viral" else v,
            "host_protein": v if G.nodes[v].get("node_type") == "host" else u,
            "n_slims": d.get("n_slims", 1),
            "max_composite": d.get("max_composite", 0),
            "cns_edge_weight": d.get("cns_edge_weight", 0),
            "domain": d.get("domain", ""),
            "binding_mode": d.get("binding_mode", ""),
            "motif_classes": "; ".join(d.get("motif_classes", [])),
            "pathways": "; ".join(d.get("pathways", set())),
        })
    edges_df = pd.DataFrame(edges_data)
    edges_csv = os.path.join(RESULTS_DIR, "interactome_edges.csv")
    edges_df.to_csv(edges_csv, index=False)
    print(f"  {edges_csv}")

    # Network metrics
    metrics_csv = os.path.join(RESULTS_DIR, "network_metrics.csv")
    metrics_list = []
    for node, m in net_metrics.items():
        entry = {"node": node}
        entry.update({k: v for k, v in m.items() if k != "pathways" and k != "viral_partners"})
        if "pathways" in m:
            entry["pathways"] = "; ".join(m["pathways"])
        if "viral_partners" in m:
            entry["viral_partners"] = "; ".join(m["viral_partners"])
        metrics_list.append(entry)
    pd.DataFrame(metrics_list).to_csv(metrics_csv, index=False)
    print(f"  {metrics_csv}")

    # H1 test results
    h1_json = os.path.join(RESULTS_DIR, "hypothesis_test_h1.json")
    with open(h1_json, "w") as f:
        json.dump(h1_results, f, indent=2, default=str)
    print(f"  {h1_json}")

    # Enrichment results
    enrich_json = os.path.join(RESULTS_DIR, "pathway_enrichment.json")
    with open(enrich_json, "w") as f:
        json.dump(enrichment, f, indent=2, default=str)
    print(f"  {enrich_json}")

    # Full summary
    summary = {
        "network": {
            "nodes": G.number_of_nodes(),
            "edges": G.number_of_edges(),
            "viral_proteins": viral_count,
            "host_proteins": host_count,
        },
        "hypothesis_test": {
            "target_mean_h1": h1_results["target_mean"],
            "control_mean_h1": h1_results["control_mean"],
            "mann_whitney_p": h1_results["mann_whitney_p"],
            "effect_size": h1_results["effect_size_cohens_d"],
            "conclusion": "H1 supported" if (h1_results["mann_whitney_p"] is not None and h1_results["mann_whitney_p"] < 0.05) else
                         "H1 trend" if (h1_results["mann_whitney_p"] is not None and h1_results["mann_whitney_p"] < 0.1) else
                         "H1 not significant (discuss pathway specificity)",
        },
        "structural_context": {
            "target_mean_accessibility": round(target_struct["accessibility_score"].mean(), 3),
            "control_mean_accessibility": round(control_struct["accessibility_score"].mean(), 3),
        },
    }
    summary_json = os.path.join(RESULTS_DIR, "phase3_summary.json")
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  {summary_json}")

    # GraphML for Cytoscape
    graphml_path = os.path.join(RESULTS_DIR, "interactome.graphml")
    G_export = G.copy()
    for u, v, d in G_export.edges(data=True):
        d["pathways"] = "; ".join(d.get("pathways", set()))
        d["motif_classes"] = "; ".join(d.get("motif_classes", []))
    nx.write_graphml(G_export, graphml_path)
    print(f"  {graphml_path} (for Cytoscape)")

    print(f"\nPhase 3 complete! Results in {RESULTS_DIR}")


if __name__ == "__main__":
    main()
