"""
Phase 2 REFINED — Comprehensive SLiM Discovery, Conservation, Statistics & Mimicry

Refinements over initial Phase 2:
  1. Expanded ELM database (~100 motif classes from elm_database.py)
  2. Known validated HSV-1 SLiMs as positive controls
  3. Conservation filtering via ortholog alignments (where available)
  4. Fisher's exact test for pathway-level enrichment
  5. Molecular mimicry scoring (viral SLiMs vs human proteome consensus)
  6. Per-IDR-residue normalization (not raw counts)
  7. Publication-quality figures with statistical annotations
"""

import re
import os
import json
import csv
import sys
import numpy as np
import pandas as pd
from collections import defaultdict
from itertools import combinations
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

BASE_DIR = r"C:\Users\shiva\OneDrive\Desktop\HSV-CNS-Interactome"
PHASE1_DIR = os.path.join(BASE_DIR, "phase1")
PHASE2_DIR = os.path.join(BASE_DIR, "phase2")
RESULTS_DIR = os.path.join(PHASE2_DIR, "results")
PLOT_DIR = os.path.join(PHASE2_DIR, "plots")
ORTHO_DIR = os.path.join(PHASE1_DIR, "sequences", "orthologs")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)

# Import expanded ELM database
sys.path.insert(0, PHASE2_DIR)
from elm_database import ELM_CLASSES, KNOWN_HSV1_SLIMS

# ═══════════════════════════════════════════════════════════════
# CNS EXPRESSION DATA (GTEx brain median TPM)
# ═══════════════════════════════════════════════════════════════
CNS_EXPRESSION = {
    "MAP1LC3A": 45, "MAP1LC3B": 120, "GABARAP": 180, "GABARAPL2": 95,
    "BECN1": 35, "BCL2": 15, "BCL2L1": 55, "MCL1": 30,
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
    "FBP11": 10, "FBP21": 8, "FE65": 25,
    "LNX1": 8, "LNX2": 5,
}

CNS_PATHWAY_WEIGHTS = {
    "autophagy": 3.0, "apoptosis": 2.5, "phosphatase_signaling": 2.5,
    "axonal_transport": 3.0, "interferon_signaling": 2.0,
    "NFkB_innate_immunity": 2.0, "MAPK_signaling": 1.5,
    "PI3K_Akt_signaling": 1.5, "JAK_STAT_signaling": 1.5,
    "14-3-3_signaling": 1.5, "Wnt_signaling": 1.0,
    "cell_cycle": 0.5, "CK2_signaling": 1.0, "cAMP_signaling": 1.5,
    "nuclear_transport": 1.0, "ubiquitin_proteasome": 1.0,
    "ubiquitin_signaling": 1.0, "translation": 2.0,
    "PDZ_signaling": 1.5, "SH3_signaling": 1.0, "SH2_signaling": 1.0,
    "RTK_signaling": 1.0, "DNA_damage_signaling": 0.5,
    "SUMOylation": 0.5, "endocytosis": 1.0, "p53_pathway": 1.0,
    "glycosylation": 0.5, "chaperone": 1.0, "GPCR_signaling": 1.0,
    "WW_signaling": 0.5, "splicing": 0.5, "DNA_replication": 0.5,
}

# Molecular mimicry: known human proteome instances of key motifs
HUMAN_MOTIF_CONSENSUS = {
    "LIG_LIR_Gen_1": {"consensus": "DWEEV", "examples": ["SQSTM1:DDDWTHL", "NBR1:EDYIII", "OPTN:DFVEI"]},
    "LIG_PP1_RVXF_1": {"consensus": "KRVXF", "examples": ["NIPP1:KRVTF", "MYPT1:KVKF", "I2:KRVQF"]},
    "LIG_PP1_SILK_1": {"consensus": "GILK", "examples": ["I3:GILK", "Ki67:SILK"]},
    "DOC_MAPK_DCC_7": {"consensus": "KKXXLXL", "examples": ["MEF2A:RKPDLRVV", "MKP3:RRLQKGNLPV"]},
    "MOD_CK2_1": {"consensus": "SXXE", "examples": ["EIF2S1:SLDE", "HDAC1:SQEE"]},
    "MOD_PKA_2": {"consensus": "RRXS", "examples": ["BAD:RSRS", "CREB1:RRPS"]},
    "MOD_CDK_SPxK_1": {"consensus": "SPXK", "examples": ["RB1:SPIK", "CDT1:SPCK"]},
    "LIG_TRAF2_1": {"consensus": "PXQE", "examples": ["TNFR2:SKEE", "CD40:PVQE"]},
    "LIG_TRAF6_MATH_1": {"consensus": "PXEXXF", "examples": ["IRAK1:PTEDEY", "TIFA:PTEAAF"]},
    "CLV_CASP_3": {"consensus": "DEVD", "examples": ["PARP1:DEVD", "ICAD:DEPD"]},
    "LIG_BH3_BclXL_1": {"consensus": "LXXXXD", "examples": ["BAD:NLWAAQRY", "BIM:IWIAQELR"]},
    "TRG_NES_CRM1_1": {"consensus": "LXXXLXXLXL", "examples": ["PKI:LALKLAGLDI", "HIV_REV:LQLPPLERLTL"]},
    "LIG_Dynein_DLC8_1": {"consensus": "KXTQT", "examples": ["BMPR2:KSTQT", "NUP159:KATQT"]},
    "LIG_eIF4E_1": {"consensus": "YXXXXLF", "examples": ["4EBP1:YDRKFLM", "EIF4G:YDREFLL"]},
    "MOD_SUMO_for_1": {"consensus": "VKXE", "examples": ["RANGAP1:VKTE", "PML:IKSE"]},
    "LIG_SH3_3": {"consensus": "PXXPXR", "examples": ["SOS1:PPVPPR", "WASP:PPLPPR"]},
}

# ═══════════════════════════════════════════════════════════════
# ORTHOLOG MAPPING (only verified correct orthologs)
# ═══════════════════════════════════════════════════════════════
ORTHOLOG_MAP = {
    "pUL37": [
        {"file": "HSV2_UL37.fasta", "species": "HSV-2", "gene": "UL37"},
        {"file": "PRV_UL37.fasta", "species": "PRV", "gene": "UL37"},
    ],
    "US3": [
        {"file": "PRV_US3.fasta", "species": "PRV", "gene": "US3"},
    ],
}


# ═══════════════════════════════════════════════════════════════
# CORE FUNCTIONS
# ═══════════════════════════════════════════════════════════════

def read_fasta(filepath):
    with open(filepath) as f:
        lines = f.read().strip().split("\n")
    header = lines[0]
    seq = "".join(l.strip() for l in lines[1:] if not l.startswith(">"))
    return header, seq


def load_disorder_data():
    path = os.path.join(PHASE1_DIR, "results", "disorder_summary_v3.json")
    with open(path) as f:
        return {d["protein"]: d for d in json.load(f)}


def load_per_residue(protein_name):
    csv_path = os.path.join(PHASE1_DIR, "disorder_predictions", f"{protein_name}_per_residue.csv")
    if os.path.exists(csv_path):
        return pd.read_csv(csv_path)
    return None


def load_sequence(protein_name, accession):
    seq_dir = os.path.join(PHASE1_DIR, "sequences")
    for fn in os.listdir(seq_dir):
        if fn.endswith(".fasta") and protein_name.lower() in fn.lower():
            _, seq = read_fasta(os.path.join(seq_dir, fn))
            return seq
    iupred_dir = os.path.join(PHASE1_DIR, "iupred3_data")
    for fn in os.listdir(iupred_dir):
        fn_norm = fn.replace("_", "").replace("-", "").lower()
        pn_norm = protein_name.replace("_", "").replace("-", "").lower()
        if pn_norm in fn_norm:
            seq_chars = []
            with open(os.path.join(iupred_dir, fn)) as f:
                for line in f:
                    if line.startswith("#") or line.strip() == "":
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 2 and len(parts[1]) == 1 and parts[1].isalpha():
                        seq_chars.append(parts[1])
            if seq_chars:
                return "".join(seq_chars)
    return None


def scan_protein(sequence, idr_boundaries):
    hits = []
    for motif_name, motif_info in ELM_CLASSES.items():
        pattern = motif_info["regex"]
        for match in re.finditer(pattern, sequence):
            start = match.start() + 1
            end = match.end()
            matched_seq = match.group()
            in_idr = any(s <= start and end <= e for s, e in idr_boundaries)
            if not in_idr:
                continue
            hits.append({
                "motif_class": motif_name,
                "description": motif_info["description"],
                "start": start,
                "end": end,
                "matched_sequence": matched_seq,
                "in_idr": True,
                "targets": motif_info["targets"],
                "pathway": motif_info["pathway"],
            })
    return hits


def score_disorder_context(hit, per_residue_df):
    if per_residue_df is None:
        return 0.5, 0.5
    region = per_residue_df.iloc[hit["start"]-1:hit["end"]]
    iup = region["iupred3"].mean() if "iupred3" in region.columns else 0.5
    comp = region["composition"].mean() if "composition" in region.columns else 0.5
    return float(iup), float(comp)


def compute_cns_score(targets, pathway):
    tpms = [CNS_EXPRESSION.get(t, 0) for t in targets]
    tpms = [t for t in tpms if t > 0]
    if not tpms:
        return 0.0, 0.0
    mean_tpm = np.mean(tpms)
    max_tpm = max(tpms)
    pw_weight = CNS_PATHWAY_WEIGHTS.get(pathway, 1.0)
    return round(np.log1p(mean_tpm) * pw_weight, 3), round(max_tpm, 1)


def find_overlapping_slims(hits, max_gap=5):
    if len(hits) < 2:
        return []
    sorted_hits = sorted(hits, key=lambda h: h["start"])
    overlaps = []
    for i in range(len(sorted_hits)):
        for j in range(i + 1, len(sorted_hits)):
            h1, h2 = sorted_hits[i], sorted_hits[j]
            if h1["motif_class"] == h2["motif_class"]:
                continue
            gap = h2["start"] - h1["end"]
            if gap <= max_gap:
                overlaps.append({
                    "motif_1": h1["motif_class"],
                    "motif_2": h2["motif_class"],
                    "motif_1_range": f"{h1['start']}-{h1['end']}",
                    "motif_2_range": f"{h2['start']}-{h2['end']}",
                    "gap": gap,
                    "pathway_1": h1["pathway"],
                    "pathway_2": h2["pathway"],
                    "switch_type": "overlapping" if gap < 0 else ("adjacent" if gap == 0 else "proximal"),
                })
    return overlaps


# ═══════════════════════════════════════════════════════════════
# CONSERVATION ANALYSIS
# ═══════════════════════════════════════════════════════════════

def simple_pairwise_identity(seq1, seq2, region_start, region_end):
    """Compute sequence identity in a local region using simple Needleman-Wunsch-like alignment."""
    s1 = seq1[region_start-1:region_end]
    if len(s1) == 0:
        return 0.0
    best_identity = 0.0
    window = len(s1)
    for offset in range(max(0, region_start - 50), min(len(seq2) - window + 1, region_end + 50)):
        s2_window = seq2[offset:offset + window]
        if len(s2_window) < window:
            continue
        matches = sum(1 for a, b in zip(s1, s2_window) if a == b)
        identity = matches / window
        best_identity = max(best_identity, identity)
    return best_identity


def compute_conservation_score(protein_name, hit, query_seq):
    """Score motif conservation across orthologs."""
    if protein_name not in ORTHOLOG_MAP:
        return None, "no_orthologs"
    orthologs = ORTHOLOG_MAP[protein_name]
    identities = []
    for orth in orthologs:
        fasta_path = os.path.join(ORTHO_DIR, orth["file"])
        if not os.path.exists(fasta_path):
            continue
        _, orth_seq = read_fasta(fasta_path)
        identity = simple_pairwise_identity(query_seq, orth_seq, hit["start"], hit["end"])
        identities.append(identity)
    if not identities:
        return None, "file_missing"
    mean_id = np.mean(identities)
    return round(mean_id, 3), "computed"


# ═══════════════════════════════════════════════════════════════
# MOLECULAR MIMICRY SCORING
# ═══════════════════════════════════════════════════════════════

def score_mimicry(hit):
    """Score how closely a viral SLiM matches the human consensus for that motif class."""
    motif_class = hit["motif_class"]
    if motif_class not in HUMAN_MOTIF_CONSENSUS:
        return 0.0, "no_data"
    consensus_data = HUMAN_MOTIF_CONSENSUS[motif_class]
    consensus = consensus_data["consensus"]
    viral_seq = hit["matched_sequence"]
    min_len = min(len(consensus), len(viral_seq))
    if min_len == 0:
        return 0.0, "no_data"
    matches = sum(1 for a, b in zip(viral_seq[:min_len], consensus[:min_len]) if a == b)
    identity = matches / min_len
    return round(identity, 3), "scored"


# ═══════════════════════════════════════════════════════════════
# POSITIVE CONTROL VALIDATION
# ═══════════════════════════════════════════════════════════════

def check_positive_controls(all_hits):
    """Check if known validated HSV-1 SLiMs are recovered by the scanner."""
    results = {}
    for ctrl_name, ctrl in KNOWN_HSV1_SLIMS.items():
        protein = ctrl["protein"]
        protein_hits = [h for h in all_hits if h["protein"] == protein]
        recovered = False
        best_match = None
        for h in protein_hits:
            if ctrl["start"] is not None and ctrl["end"] is not None:
                overlap = max(0, min(h["end"], ctrl["end"]) - max(h["start"], ctrl["start"]))
                if overlap > 0:
                    recovered = True
                    if best_match is None or overlap > best_match.get("overlap", 0):
                        best_match = {**h, "overlap": overlap}
            else:
                if ctrl["pathway"] == h["pathway"]:
                    recovered = True
                    if best_match is None:
                        best_match = h
        results[ctrl_name] = {
            "expected": ctrl,
            "recovered": recovered,
            "best_match": best_match,
        }
    return results


# ═══════════════════════════════════════════════════════════════
# STATISTICAL TESTING
# ═══════════════════════════════════════════════════════════════

def pathway_enrichment_test(target_hits, control_hits, target_idr_aa, control_idr_aa):
    """Fisher's exact test for pathway-level enrichment in targets vs controls."""
    all_pathways = set(h["pathway"] for h in target_hits + control_hits)
    results = []
    for pathway in sorted(all_pathways):
        n_target = sum(1 for h in target_hits if h["pathway"] == pathway)
        n_control = sum(1 for h in control_hits if h["pathway"] == pathway)
        n_target_other = len(target_hits) - n_target
        n_control_other = len(control_hits) - n_control
        table = [[n_target, n_control], [n_target_other, n_control_other]]
        odds_ratio, p_value = stats.fisher_exact(table, alternative="greater")
        target_density = n_target / target_idr_aa * 100 if target_idr_aa > 0 else 0
        control_density = n_control / control_idr_aa * 100 if control_idr_aa > 0 else 0
        density_ratio = target_density / control_density if control_density > 0 else float('inf')
        results.append({
            "pathway": pathway,
            "target_hits": n_target,
            "control_hits": n_control,
            "target_density_per100aa": round(target_density, 2),
            "control_density_per100aa": round(control_density, 2),
            "density_ratio": round(density_ratio, 2),
            "odds_ratio": round(odds_ratio, 3),
            "p_value": round(p_value, 6),
            "significant": p_value < 0.05,
        })
    results.sort(key=lambda x: x["p_value"])
    n_tests = len(results)
    for i, r in enumerate(results):
        r["p_adjusted_bh"] = round(min(r["p_value"] * n_tests / (i + 1), 1.0), 6)
        r["significant_adjusted"] = r["p_adjusted_bh"] < 0.1
    return results


def permutation_enrichment_test(target_hits, control_hits, target_idr_aa, control_idr_aa, n_perm=10000):
    """Permutation test for overall SLiM density enrichment."""
    target_density = len(target_hits) / target_idr_aa if target_idr_aa > 0 else 0
    control_density = len(control_hits) / control_idr_aa if control_idr_aa > 0 else 0
    observed_ratio = target_density / control_density if control_density > 0 else float('inf')
    all_hits = target_hits + control_hits
    n_target = len(target_hits)
    rng = np.random.default_rng(42)
    perm_ratios = []
    for _ in range(n_perm):
        perm_idx = rng.permutation(len(all_hits))
        perm_target = perm_idx[:n_target]
        perm_control = perm_idx[n_target:]
        pt_density = len(perm_target) / target_idr_aa if target_idr_aa > 0 else 0
        pc_density = len(perm_control) / control_idr_aa if control_idr_aa > 0 else 0
        pr = pt_density / pc_density if pc_density > 0 else 0
        perm_ratios.append(pr)
    p_value = np.mean([r >= observed_ratio for r in perm_ratios])
    return {
        "observed_ratio": round(observed_ratio, 4),
        "p_value_permutation": round(p_value, 4),
        "n_permutations": n_perm,
        "mean_perm_ratio": round(np.mean(perm_ratios), 4),
        "std_perm_ratio": round(np.std(perm_ratios), 4),
    }


# ═══════════════════════════════════════════════════════════════
# COMPOSITE SCORING
# ═══════════════════════════════════════════════════════════════

def compute_composite(hit):
    cns = hit.get("cns_score", 0)
    disorder = hit.get("iupred3_mean", 0.5)
    conservation = hit.get("conservation_score")
    mimicry = hit.get("mimicry_score", 0)
    cns_norm = min(cns / 15, 1.0)
    disorder_norm = disorder
    cons_norm = conservation if conservation is not None else 0.3
    mimicry_norm = mimicry
    composite = (0.35 * cns_norm +
                 0.20 * disorder_norm +
                 0.25 * cons_norm +
                 0.20 * mimicry_norm)
    return round(composite, 4)


# ═══════════════════════════════════════════════════════════════
# PLOTTING
# ═══════════════════════════════════════════════════════════════

PATHWAY_COLORS = {
    "autophagy": "#e74c3c", "apoptosis": "#c0392b",
    "phosphatase_signaling": "#2980b9", "axonal_transport": "#8e44ad",
    "interferon_signaling": "#d35400", "NFkB_innate_immunity": "#e67e22",
    "MAPK_signaling": "#27ae60", "cell_cycle": "#95a5a6",
    "CK2_signaling": "#16a085", "nuclear_transport": "#f39c12",
    "14-3-3_signaling": "#1abc9c", "translation": "#2c3e50",
    "PI3K_Akt_signaling": "#3498db", "ubiquitin_proteasome": "#7f8c8d",
    "ubiquitin_signaling": "#bdc3c7", "cAMP_signaling": "#9b59b6",
    "PDZ_signaling": "#1abc9c", "JAK_STAT_signaling": "#e67e22",
    "SH3_signaling": "#34495e", "SH2_signaling": "#2c3e50",
    "endocytosis": "#f1c40f", "DNA_damage_signaling": "#95a5a6",
    "p53_pathway": "#c0392b", "Wnt_signaling": "#16a085",
    "SUMOylation": "#7f8c8d", "glycosylation": "#bdc3c7",
    "chaperone": "#d35400", "GPCR_signaling": "#9b59b6",
    "DNA_replication": "#95a5a6", "RTK_signaling": "#e67e22",
    "WW_signaling": "#34495e", "splicing": "#7f8c8d",
}


def plot_protein_slim_map(protein_name, length, idr_boundaries, hits, overlaps, output_path):
    fig, axes = plt.subplots(5, 1, figsize=(16, 12), sharex=True,
                              gridspec_kw={"height_ratios": [1, 3, 1.5, 1.5, 1]})
    fig.suptitle(f"{protein_name} — Refined SLiM Architecture Map",
                 fontsize=14, fontweight="bold", y=0.97)

    # Panel 1: IDR map
    ax0 = axes[0]
    ax0.plot([1, length], [0.5, 0.5], color="gray", linewidth=8, alpha=0.2, solid_capstyle="round")
    for s, e in idr_boundaries:
        ax0.plot([s, e], [0.5, 0.5], color="#3498db", linewidth=8, alpha=0.5, solid_capstyle="round")
    ax0.set_ylabel("IDR", fontsize=10)
    ax0.set_yticks([])
    ax0.set_ylim(0, 1)

    # Panel 2: SLiM positions
    ax1 = axes[1]
    y_offset = 0
    motif_classes_seen = {}
    for hit in sorted(hits, key=lambda h: h["start"]):
        motif = hit["motif_class"]
        if motif not in motif_classes_seen:
            motif_classes_seen[motif] = y_offset
            y_offset += 1
        y = motif_classes_seen[motif]
        color = PATHWAY_COLORS.get(hit.get("pathway", ""), "#7f8c8d")
        alpha = 0.9 if hit.get("composite_score", 0) > 0.3 else 0.5
        ax1.barh(y, hit["end"] - hit["start"] + 1, left=hit["start"] - 0.5,
                 height=0.7, color=color, alpha=alpha, edgecolor="white", linewidth=0.5)
        if len(hits) <= 25:
            ax1.text((hit["start"] + hit["end"]) / 2, y, hit["matched_sequence"],
                     fontsize=5, ha="center", va="center", color="white", fontweight="bold")
    ax1.set_yticks(range(len(motif_classes_seen)))
    ax1.set_yticklabels(list(motif_classes_seen.keys()), fontsize=6)
    ax1.set_ylabel("SLiM Motifs", fontsize=10)
    ax1.invert_yaxis()

    # Panel 3: CNS score track
    ax2 = axes[2]
    cns_track = np.zeros(length)
    for hit in hits:
        cns = hit.get("cns_score", 0)
        for pos in range(hit["start"] - 1, min(hit["end"], length)):
            cns_track[pos] = max(cns_track[pos], cns)
    ax2.fill_between(range(1, length + 1), cns_track, alpha=0.6, color="#e74c3c")
    ax2.set_ylabel("CNS Score", fontsize=10)
    ax2.set_ylim(0, max(cns_track.max() * 1.1, 1))

    # Panel 4: Composite score track
    ax3 = axes[3]
    comp_track = np.zeros(length)
    for hit in hits:
        comp = hit.get("composite_score", 0)
        for pos in range(hit["start"] - 1, min(hit["end"], length)):
            comp_track[pos] = max(comp_track[pos], comp)
    ax3.fill_between(range(1, length + 1), comp_track, alpha=0.6, color="#2980b9")
    ax3.set_ylabel("Composite", fontsize=10)
    ax3.set_ylim(0, max(comp_track.max() * 1.1, 0.1))

    # Panel 5: Overlap density
    ax4 = axes[4]
    overlap_track = np.zeros(length)
    for ov in overlaps:
        r1 = [int(x) for x in ov["motif_1_range"].split("-")]
        r2 = [int(x) for x in ov["motif_2_range"].split("-")]
        o_start = min(r1[0], r2[0])
        o_end = max(r1[1], r2[1])
        for pos in range(o_start - 1, min(o_end, length)):
            overlap_track[pos] += 1
    ax4.fill_between(range(1, length + 1), overlap_track, alpha=0.6, color="#8e44ad")
    ax4.set_ylabel("Overlaps", fontsize=10)
    ax4.set_xlabel("Residue Position", fontsize=11)
    ax4.set_xlim(1, length)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


def plot_enrichment_dashboard(pathway_stats, perm_result, target_summaries, control_summaries, output_path):
    fig = plt.figure(figsize=(20, 14))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    # Panel A: Per-protein SLiM density
    ax_a = fig.add_subplot(gs[0, 0])
    all_summaries = target_summaries + control_summaries
    names = [s["protein"] for s in all_summaries]
    densities = [s["slim_density_per100aa"] for s in all_summaries]
    colors = ["#e74c3c" if s["category"] == "target" else "#3498db" for s in all_summaries]
    bars = ax_a.barh(range(len(names)), densities, color=colors, alpha=0.7)
    ax_a.set_yticks(range(len(names)))
    ax_a.set_yticklabels(names, fontsize=9)
    ax_a.set_xlabel("SLiMs per 100 IDR residues", fontsize=10)
    ax_a.set_title("A. Raw SLiM Density", fontsize=11, fontweight="bold")
    ax_a.invert_yaxis()
    ax_a.legend([mpatches.Patch(color="#e74c3c"), mpatches.Patch(color="#3498db")],
                ["Neurovirulence", "Control"], fontsize=8, loc="lower right")

    # Panel B: CNS-weighted density
    ax_b = fig.add_subplot(gs[0, 1])
    cns_densities = [s["cns_density_per100aa"] for s in all_summaries]
    bars = ax_b.barh(range(len(names)), cns_densities, color=colors, alpha=0.7)
    ax_b.set_yticks(range(len(names)))
    ax_b.set_yticklabels(names, fontsize=9)
    ax_b.set_xlabel("CNS-weighted score per 100 IDR residues", fontsize=10)
    ax_b.set_title("B. CNS-Weighted SLiM Density", fontsize=11, fontweight="bold")
    ax_b.invert_yaxis()

    # Panel C: Pathway enrichment volcano
    ax_c = fig.add_subplot(gs[0, 2])
    sig_pathways = [p for p in pathway_stats if p["target_hits"] + p["control_hits"] >= 3]
    if sig_pathways:
        x = [np.log2(p["density_ratio"]) if p["density_ratio"] > 0 and p["density_ratio"] != float('inf') else 3
             for p in sig_pathways]
        y = [-np.log10(max(p["p_value"], 1e-10)) for p in sig_pathways]
        c = ["#e74c3c" if p["significant"] else "#95a5a6" for p in sig_pathways]
        ax_c.scatter(x, y, c=c, s=60, alpha=0.7, edgecolors="black", linewidth=0.5)
        for i, p in enumerate(sig_pathways):
            if p["significant"] or p["density_ratio"] > 2:
                ax_c.annotate(p["pathway"].replace("_", "\n"), (x[i], y[i]),
                             fontsize=6, ha="center", va="bottom")
        ax_c.axhline(-np.log10(0.05), color="red", linestyle="--", alpha=0.5, linewidth=0.8)
        ax_c.axvline(0, color="gray", linestyle="--", alpha=0.3)
    ax_c.set_xlabel("log2(density ratio target/control)", fontsize=10)
    ax_c.set_ylabel("-log10(p-value)", fontsize=10)
    ax_c.set_title("C. Pathway Enrichment (Fisher's Exact)", fontsize=11, fontweight="bold")

    # Panel D: Top composite-scored SLiMs
    ax_d = fig.add_subplot(gs[1, 0])
    all_target_hits = []
    for s in target_summaries:
        all_target_hits.extend(s.get("top_hits", []))
    all_target_hits.sort(key=lambda x: -x.get("composite_score", 0))
    top_n = min(15, len(all_target_hits))
    if top_n > 0:
        labels = [f"{h['protein']}:{h['motif_class']}" for h in all_target_hits[:top_n]]
        scores = [h.get("composite_score", 0) for h in all_target_hits[:top_n]]
        bar_colors = [PATHWAY_COLORS.get(h.get("pathway", ""), "#7f8c8d") for h in all_target_hits[:top_n]]
        ax_d.barh(range(top_n), scores, color=bar_colors, alpha=0.7)
        ax_d.set_yticks(range(top_n))
        ax_d.set_yticklabels(labels, fontsize=7)
        ax_d.invert_yaxis()
    ax_d.set_xlabel("Composite Score", fontsize=10)
    ax_d.set_title("D. Top SLiM Candidates", fontsize=11, fontweight="bold")

    # Panel E: Pathway breakdown heatmap
    ax_e = fig.add_subplot(gs[1, 1])
    target_pw = defaultdict(int)
    control_pw = defaultdict(int)
    for s in target_summaries:
        for pw, cnt in s.get("pathway_counts", {}).items():
            target_pw[pw] += cnt
    for s in control_summaries:
        for pw, cnt in s.get("pathway_counts", {}).items():
            control_pw[pw] += cnt
    all_pw = sorted(set(list(target_pw.keys()) + list(control_pw.keys())))
    if all_pw:
        pw_display = all_pw[:15]
        t_counts = [target_pw.get(p, 0) for p in pw_display]
        c_counts = [control_pw.get(p, 0) for p in pw_display]
        x = np.arange(len(pw_display))
        width = 0.35
        ax_e.barh(x - width/2, t_counts, width, color="#e74c3c", alpha=0.7, label="Target")
        ax_e.barh(x + width/2, c_counts, width, color="#3498db", alpha=0.7, label="Control")
        ax_e.set_yticks(x)
        ax_e.set_yticklabels([p.replace("_", " ") for p in pw_display], fontsize=7)
        ax_e.invert_yaxis()
        ax_e.legend(fontsize=8)
    ax_e.set_xlabel("Number of SLiM hits", fontsize=10)
    ax_e.set_title("E. Pathway Hit Distribution", fontsize=11, fontweight="bold")

    # Panel F: Positive control results
    ax_f = fig.add_subplot(gs[1, 2])
    ax_f.axis("off")
    text_lines = ["POSITIVE CONTROLS", "=" * 35]
    ctrl_results = check_positive_controls(all_target_hits)
    for name, res in ctrl_results.items():
        status = "RECOVERED" if res["recovered"] else "NOT FOUND"
        icon = "+" if res["recovered"] else "-"
        text_lines.append(f"[{icon}] {name}: {status}")
        if res["best_match"]:
            bm = res["best_match"]
            text_lines.append(f"    {bm.get('motif_class','')} "
                            f"({bm.get('start','')}-{bm.get('end','')})")
    text_lines.append("")
    text_lines.append("PERMUTATION TEST")
    text_lines.append("=" * 35)
    text_lines.append(f"Observed density ratio: {perm_result['observed_ratio']:.3f}")
    text_lines.append(f"Permutation p-value:   {perm_result['p_value_permutation']:.4f}")
    text_lines.append(f"Mean null ratio:       {perm_result['mean_perm_ratio']:.3f}")
    ax_f.text(0.05, 0.95, "\n".join(text_lines), transform=ax_f.transAxes,
              fontsize=9, verticalalignment="top", fontfamily="monospace",
              bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    ax_f.set_title("F. Validation & Statistics", fontsize=11, fontweight="bold")

    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


# ═══════════════════════════════════════════════════════════════
# MAIN PIPELINE
# ═══════════════════════════════════════════════════════════════

def main():
    print("=" * 80)
    print("PHASE 2 REFINED — SLiM Discovery with Conservation, Statistics & Mimicry")
    print("=" * 80)

    disorder_data = load_disorder_data()

    # ── Step 1: Scan all proteins with expanded ELM database ──
    print(f"\n[1] Scanning {len(disorder_data)} proteins with {len(ELM_CLASSES)} ELM motif classes...")

    all_hits = []
    all_overlaps = []

    for protein_name, pdata in disorder_data.items():
        seq = load_sequence(protein_name, pdata["accession"])
        if seq is None:
            print(f"  WARNING: No sequence for {protein_name}, skipping")
            continue

        idr_boundaries = [tuple(b) for b in pdata["idr_boundaries"]]
        per_residue = load_per_residue(protein_name)

        hits = scan_protein(seq, idr_boundaries)

        for h in hits:
            h["protein"] = protein_name
            h["category"] = pdata["category"]
            iup, comp = score_disorder_context(h, per_residue)
            h["iupred3_mean"] = iup
            h["composition_mean"] = comp
            cns, max_tpm = compute_cns_score(h["targets"], h["pathway"])
            h["cns_score"] = cns
            h["max_brain_tpm"] = max_tpm
            cons_score, cons_status = compute_conservation_score(protein_name, h, seq)
            h["conservation_score"] = cons_score
            h["conservation_status"] = cons_status
            mimicry, mim_status = score_mimicry(h)
            h["mimicry_score"] = mimicry
            h["mimicry_status"] = mim_status
            h["composite_score"] = compute_composite(h)

        overlaps = find_overlapping_slims(hits)
        for ov in overlaps:
            ov["protein"] = protein_name
            ov["category"] = pdata["category"]

        all_hits.extend(hits)
        all_overlaps.extend(overlaps)

        idr_aa = sum(e - s + 1 for s, e in idr_boundaries)
        print(f"  {protein_name:<15} ({pdata['category']:<8}) — "
              f"{len(hits):>3} hits in {idr_aa:>4} IDR aa, "
              f"{len(overlaps)} overlaps")

    print(f"\n  TOTAL: {len(all_hits)} hits, {len(all_overlaps)} overlaps")

    # ── Step 2: Positive control validation ──
    print(f"\n[2] Checking positive controls...")
    ctrl_results = check_positive_controls(all_hits)
    for name, res in ctrl_results.items():
        status = "RECOVERED" if res["recovered"] else "MISSED"
        print(f"  {name:<20}: {status}")
        if res["best_match"]:
            bm = res["best_match"]
            print(f"    -> {bm['motif_class']} at {bm['start']}-{bm['end']} "
                  f"'{bm['matched_sequence']}' (composite={bm.get('composite_score', 0):.3f})")

    # ── Step 3: Statistical testing ──
    print(f"\n[3] Running statistical tests...")
    target_hits = [h for h in all_hits if h["category"] == "target"]
    control_hits = [h for h in all_hits if h["category"] == "control"]
    target_idr_aa = sum(
        sum(e - s + 1 for s, e in [tuple(b) for b in disorder_data[p]["idr_boundaries"]])
        for p in disorder_data if disorder_data[p]["category"] == "target"
    )
    control_idr_aa = sum(
        sum(e - s + 1 for s, e in [tuple(b) for b in disorder_data[p]["idr_boundaries"]])
        for p in disorder_data if disorder_data[p]["category"] == "control"
    )

    print(f"  Target: {len(target_hits)} hits in {target_idr_aa} IDR aa "
          f"({len(target_hits)/target_idr_aa*100:.1f}/100aa)")
    print(f"  Control: {len(control_hits)} hits in {control_idr_aa} IDR aa "
          f"({len(control_hits)/control_idr_aa*100:.1f}/100aa)")

    # Fisher's exact by pathway
    pathway_stats = pathway_enrichment_test(target_hits, control_hits, target_idr_aa, control_idr_aa)
    print(f"\n  Pathway enrichment (Fisher's exact, one-sided):")
    print(f"  {'Pathway':<30} {'T':>4} {'C':>4} {'T/100aa':>8} {'C/100aa':>8} "
          f"{'Ratio':>6} {'p':>8} {'Sig':>5}")
    print("  " + "-" * 90)
    for p in pathway_stats[:20]:
        sig_mark = "*" if p["significant"] else ""
        adj_mark = "**" if p.get("significant_adjusted") else ""
        print(f"  {p['pathway']:<30} {p['target_hits']:>4} {p['control_hits']:>4} "
              f"{p['target_density_per100aa']:>8.2f} {p['control_density_per100aa']:>8.2f} "
              f"{p['density_ratio']:>6.2f} {p['p_value']:>8.4f} {sig_mark}{adj_mark}")

    # Permutation test
    perm_result = permutation_enrichment_test(target_hits, control_hits, target_idr_aa, control_idr_aa)
    print(f"\n  Permutation test (overall density):")
    print(f"    Observed ratio: {perm_result['observed_ratio']:.4f}")
    print(f"    Permutation p:  {perm_result['p_value_permutation']:.4f}")

    # ── Step 4: Build per-protein summaries ──
    print(f"\n[4] Building protein summaries...")
    protein_summaries = []
    for protein_name, pdata in disorder_data.items():
        protein_hits = [h for h in all_hits if h["protein"] == protein_name]
        protein_overlaps = [o for o in all_overlaps if o["protein"] == protein_name]
        idr_aa = sum(e - s + 1 for s, e in [tuple(b) for b in pdata["idr_boundaries"]])

        pathway_counts = defaultdict(int)
        pathway_cns = defaultdict(float)
        for h in protein_hits:
            pathway_counts[h["pathway"]] += 1
            pathway_cns[h["pathway"]] += h.get("cns_score", 0)

        total_cns = sum(h.get("cns_score", 0) for h in protein_hits)
        top_hits = sorted(protein_hits, key=lambda x: -x.get("composite_score", 0))[:10]

        summary = {
            "protein": protein_name,
            "category": pdata["category"],
            "length": pdata["length"],
            "idr_aa": idr_aa,
            "idr_fraction": pdata["idr_fraction"],
            "n_hits": len(protein_hits),
            "n_overlaps": len(protein_overlaps),
            "slim_density_per100aa": round(len(protein_hits) / idr_aa * 100, 2) if idr_aa > 0 else 0,
            "cns_density_per100aa": round(total_cns / idr_aa * 100, 2) if idr_aa > 0 else 0,
            "total_cns_score": round(total_cns, 2),
            "mean_composite": round(np.mean([h["composite_score"] for h in protein_hits]), 4) if protein_hits else 0,
            "pathway_counts": dict(pathway_counts),
            "pathway_cns_scores": {k: round(v, 2) for k, v in pathway_cns.items()},
            "idr_boundaries": pdata["idr_boundaries"],
            "top_hits": top_hits,
        }
        protein_summaries.append(summary)

        print(f"  {protein_name:<15} ({pdata['category']:<8}): "
              f"{len(protein_hits):>3} hits, density={summary['slim_density_per100aa']:.1f}/100aa, "
              f"CNS-wd={summary['cns_density_per100aa']:.1f}/100aa, "
              f"mean_comp={summary['mean_composite']:.3f}")

    target_summaries = [s for s in protein_summaries if s["category"] == "target"]
    control_summaries = [s for s in protein_summaries if s["category"] == "control"]

    # ── Step 5: Generate plots ──
    print(f"\n[5] Generating publication plots...")
    for s in target_summaries:
        protein = s["protein"]
        protein_hits = [h for h in all_hits if h["protein"] == protein]
        protein_overlaps = [o for o in all_overlaps if o["protein"] == protein]
        plot_path = os.path.join(PLOT_DIR, f"{protein}_slim_map_refined.png")
        plot_protein_slim_map(protein, s["length"],
                             [tuple(b) for b in s["idr_boundaries"]],
                             protein_hits, protein_overlaps, plot_path)
        print(f"    {plot_path}")

    dashboard_path = os.path.join(PLOT_DIR, "enrichment_dashboard_refined.png")
    plot_enrichment_dashboard(pathway_stats, perm_result,
                             target_summaries, control_summaries, dashboard_path)
    print(f"    {dashboard_path}")

    # ── Step 6: Save all results ──
    print(f"\n[6] Saving results...")

    # All hits CSV
    hits_csv = os.path.join(RESULTS_DIR, "slim_hits_refined.csv")
    hits_df = pd.DataFrame(all_hits)
    hits_df["targets"] = hits_df["targets"].apply(lambda x: "; ".join(x) if isinstance(x, list) else str(x))
    hits_df.to_csv(hits_csv, index=False)
    print(f"  {hits_csv}")

    # Overlaps CSV
    overlaps_csv = os.path.join(RESULTS_DIR, "slim_overlaps_refined.csv")
    pd.DataFrame(all_overlaps).to_csv(overlaps_csv, index=False)
    print(f"  {overlaps_csv}")

    # Pathway enrichment CSV
    pw_csv = os.path.join(RESULTS_DIR, "pathway_enrichment.csv")
    pd.DataFrame(pathway_stats).to_csv(pw_csv, index=False)
    print(f"  {pw_csv}")

    # Full summary JSON
    summary_json = os.path.join(RESULTS_DIR, "slim_refined_summary.json")
    save_summaries = []
    for s in protein_summaries:
        ss = {k: v for k, v in s.items() if k != "top_hits"}
        ss["top_hits"] = [
            {k: v for k, v in h.items() if k != "targets" or not isinstance(v, list)}
            for h in s.get("top_hits", [])
        ]
        for th in ss["top_hits"]:
            if "targets" in th and isinstance(th["targets"], list):
                th["targets"] = "; ".join(th["targets"])
        save_summaries.append(ss)

    with open(summary_json, "w") as f:
        json.dump({
            "metadata": {
                "n_elm_classes": len(ELM_CLASSES),
                "n_proteins": len(disorder_data),
                "n_targets": len(target_summaries),
                "n_controls": len(control_summaries),
                "total_hits": len(all_hits),
                "total_overlaps": len(all_overlaps),
            },
            "enrichment": {
                "target_hits": len(target_hits),
                "control_hits": len(control_hits),
                "target_idr_aa": target_idr_aa,
                "control_idr_aa": control_idr_aa,
                "target_density_per100aa": round(len(target_hits) / target_idr_aa * 100, 2),
                "control_density_per100aa": round(len(control_hits) / control_idr_aa * 100, 2),
                "density_ratio": round(len(target_hits) / target_idr_aa / (len(control_hits) / control_idr_aa), 3) if control_idr_aa > 0 and len(control_hits) > 0 else None,
            },
            "permutation_test": perm_result,
            "positive_controls": {
                name: {
                    "recovered": res["recovered"],
                    "best_match_motif": res["best_match"]["motif_class"] if res["best_match"] else None,
                    "best_match_pos": f"{res['best_match']['start']}-{res['best_match']['end']}" if res["best_match"] else None,
                }
                for name, res in ctrl_results.items()
            },
            "pathway_enrichment_top10": pathway_stats[:10],
            "proteins": save_summaries,
        }, f, indent=2, default=str)
    print(f"  {summary_json}")

    # ── Final report ──
    print(f"\n{'='*80}")
    print("FINAL RANKED SLiM CANDIDATES (Target Proteins)")
    print(f"{'='*80}")
    target_hits_sorted = sorted(target_hits, key=lambda x: -x.get("composite_score", 0))
    print(f"{'#':<4} {'Protein':<10} {'Motif':<25} {'Pos':<10} {'Seq':<12} "
          f"{'CNS':>5} {'Cons':>5} {'Mim':>5} {'Comp':>6} {'Pathway':<22}")
    print("-" * 110)
    for i, h in enumerate(target_hits_sorted[:30], 1):
        cons_str = f"{h['conservation_score']:.2f}" if h.get('conservation_score') is not None else "  -  "
        print(f"{i:<4} {h['protein']:<10} {h['motif_class']:<25} "
              f"{h['start']}-{h['end']:<6} {h['matched_sequence']:<12} "
              f"{h['cns_score']:>5.1f} {cons_str:>5} {h.get('mimicry_score',0):>5.2f} "
              f"{h['composite_score']:>5.3f} {h['pathway']:<22}")

    print(f"\n{'='*80}")
    print("COMPETITIVE-BINDING SWITCHES (Target Proteins)")
    print(f"{'='*80}")
    target_overlaps = [o for o in all_overlaps if o["category"] == "target"]
    for ov in target_overlaps:
        print(f"  {ov['protein']}: {ov['motif_1']} ({ov['motif_1_range']}) <-> "
              f"{ov['motif_2']} ({ov['motif_2_range']}) — "
              f"{ov['switch_type']} [{ov['pathway_1']} / {ov['pathway_2']}]")

    print(f"\nDone! Results in {RESULTS_DIR}")


if __name__ == "__main__":
    main()
