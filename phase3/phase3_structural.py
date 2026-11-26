#!/usr/bin/env python3
"""
Phase 3 Structural Completion: AlphaFold pLDDT extraction, SLiM exposure
estimates, competitive accessibility analysis, and structural context profiling.

Note: Full MD simulations (CALVADOS 2) and peptide docking (FlexPepDock/HADDOCK)
require compute clusters and are documented as future work.
"""

import os, sys, json, csv, math, re, urllib.request
import numpy as np
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(BASE)
PHASE1 = os.path.join(PROJECT, "phase1")
PHASE2 = os.path.join(PROJECT, "phase2")
RESULTS = os.path.join(BASE, "results")
AF_DIR = os.path.join(BASE, "alphafold_data")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(AF_DIR, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════════════════

PROTEINS = [
    {"name": "ICP34.5", "accession": "P36313", "category": "target"},
    {"name": "US3", "accession": "P04413", "category": "target"},
    {"name": "pUL37", "accession": "P10221", "category": "target"},
    {"name": "ICP8_UL29", "accession": "P04296", "category": "control"},
    {"name": "UL42", "accession": "P10186", "category": "control"},
    {"name": "gD_US6", "accession": "P57083", "category": "control"},
    {"name": "VP11_12_UL46", "accession": "P10228", "category": "control"},
    {"name": "VP16_UL48", "accession": "P06492", "category": "control"},
    {"name": "VP22_UL49", "accession": "P10216", "category": "control"},
    {"name": "ICP27_UL54", "accession": "P10238", "category": "control"},
    {"name": "ICP22_US1", "accession": "P04485", "category": "control"},
]


def load_fasta(filepath):
    if not os.path.exists(filepath):
        return None
    seq = []
    with open(filepath) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)


def find_fasta(protein_name):
    for sub in ["targets", "controls"]:
        path = os.path.join(PHASE1, "sequences", sub, f"{protein_name}.fasta")
        if os.path.exists(path):
            return path
    pn = protein_name.replace("_", "").replace("-", "").lower()
    for sub in ["targets", "controls"]:
        d = os.path.join(PHASE1, "sequences", sub)
        if not os.path.isdir(d):
            continue
        for fn in os.listdir(d):
            if fn.endswith(".fasta"):
                fn_norm = fn.replace("_", "").replace("-", "").replace(".fasta", "").lower()
                if pn in fn_norm or fn_norm in pn:
                    return os.path.join(d, fn)
    return None


def load_iupred3(protein_name):
    iupred_dir = os.path.join(PHASE1, "iupred3_data")
    pn = protein_name.replace("_", "").replace("-", "").lower()
    for fn in os.listdir(iupred_dir):
        fn_norm = fn.replace("_", "").replace("-", "").lower()
        if pn in fn_norm:
            scores = []
            with open(os.path.join(iupred_dir, fn)) as f:
                for line in f:
                    if line.startswith("#") or line.strip() == "":
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        scores.append(float(parts[2]))
            return np.array(scores)
    return None


# ═══════════════════════════════════════════════════════════════════════
# ALPHAFOLD pLDDT EXTRACTION
# ═══════════════════════════════════════════════════════════════════════

def fetch_alphafold_plddt(accession, protein_name):
    """Fetch AlphaFold predicted structure and extract pLDDT scores."""
    pdb_path = os.path.join(AF_DIR, f"AF-{accession}-F1-model_v4.pdb")

    if os.path.exists(pdb_path):
        print(f"  Using cached AlphaFold PDB: {pdb_path}")
    else:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb"
        print(f"  Fetching AlphaFold: {url}")
        try:
            urllib.request.urlretrieve(url, pdb_path)
            print(f"  Downloaded: {pdb_path}")
        except Exception as e:
            print(f"  AlphaFold not available for {accession}: {e}")
            return None

    # Extract pLDDT from B-factor column of CA atoms
    plddt = []
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    bfactor = float(line[60:66].strip())
                    plddt.append(bfactor)
    except Exception as e:
        print(f"  Error parsing PDB: {e}")
        return None

    if plddt:
        print(f"  pLDDT scores: {len(plddt)} residues, mean={np.mean(plddt):.1f}")
    return np.array(plddt) if plddt else None


# ═══════════════════════════════════════════════════════════════════════
# STRUCTURAL CONTEXT ANALYSIS
# ═══════════════════════════════════════════════════════════════════════

# Amino acid property scales
HYDRO = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
}

CF_HELIX = {
    'A': 1.42, 'R': 0.98, 'N': 0.67, 'D': 1.01, 'C': 0.70,
    'Q': 1.11, 'E': 1.51, 'G': 0.57, 'H': 1.00, 'I': 1.08,
    'L': 1.21, 'K': 1.16, 'M': 1.45, 'F': 1.13, 'P': 0.57,
    'S': 0.77, 'T': 0.83, 'W': 1.08, 'Y': 0.69, 'V': 1.06,
}

EMINI = {
    'A': 0.815, 'R': 1.475, 'N': 1.296, 'D': 1.283, 'C': 0.394,
    'Q': 1.348, 'E': 1.445, 'G': 0.714, 'H': 1.180, 'I': 0.603,
    'L': 0.603, 'K': 1.596, 'M': 0.714, 'F': 0.603, 'P': 1.236,
    'S': 1.115, 'T': 1.018, 'W': 0.603, 'Y': 0.714, 'V': 0.603,
}

CHARGE = {'R': 1.0, 'K': 1.0, 'H': 0.1, 'D': -1.0, 'E': -1.0}

DISORDER_PROPENSITY = {
    'A': 0.06, 'R': -0.18, 'N': 0.01, 'D': 0.19, 'C': -0.34,
    'Q': 0.12, 'E': 0.30, 'G': 0.17, 'H': -0.01, 'I': -0.49,
    'L': -0.33, 'K': 0.09, 'M': -0.19, 'F': -0.42, 'P': 0.41,
    'S': 0.14, 'T': -0.05, 'W': -0.44, 'Y': -0.27, 'V': -0.46,
}


def compute_structural_context(sequence, iupred_scores, plddt_scores, start, end, flank=10):
    """Compute detailed structural context for a SLiM region."""
    n = len(sequence)
    s = max(0, start - 1)  # 0-indexed
    e = min(n, end)

    motif_seq = sequence[s:e]
    flank_start = max(0, s - flank)
    flank_end = min(n, e + flank)

    # Motif properties
    motif_hydro = np.mean([HYDRO.get(aa, 0) for aa in motif_seq])
    motif_charge = sum(CHARGE.get(aa, 0) for aa in motif_seq) / len(motif_seq)
    motif_helix = np.mean([CF_HELIX.get(aa, 1.0) for aa in motif_seq])
    motif_emini = np.mean([EMINI.get(aa, 0.8) for aa in motif_seq])
    motif_disorder_prop = np.mean([DISORDER_PROPENSITY.get(aa, 0) for aa in motif_seq])

    # Flanking properties
    flank_seq = sequence[flank_start:s] + sequence[e:flank_end]
    flank_hydro = np.mean([HYDRO.get(aa, 0) for aa in flank_seq]) if flank_seq else 0
    flank_charge = sum(CHARGE.get(aa, 0) for aa in flank_seq) / max(len(flank_seq), 1)

    # IUPred3 context
    motif_iupred = np.mean(iupred_scores[s:e]) if s < len(iupred_scores) else 0
    flank_iupred = 0
    if flank_seq:
        vals = list(iupred_scores[flank_start:s]) + list(iupred_scores[e:flank_end])
        flank_iupred = np.mean(vals) if vals else 0

    # pLDDT context
    motif_plddt = 0
    flank_plddt = 0
    if plddt_scores is not None and len(plddt_scores) > 0:
        if e <= len(plddt_scores):
            motif_plddt = np.mean(plddt_scores[s:e])
            plddt_flank_vals = list(plddt_scores[flank_start:s]) + list(plddt_scores[e:min(len(plddt_scores), flank_end)])
            flank_plddt = np.mean(plddt_flank_vals) if plddt_flank_vals else 0

    # Exposure estimate (combining Emini, disorder, and pLDDT)
    plddt_exposure = (100 - motif_plddt) / 100.0 if motif_plddt > 0 else 0.5
    exposure_score = 0.3 * motif_emini + 0.3 * motif_iupred + 0.2 * plddt_exposure + 0.2 * motif_disorder_prop

    return {
        "motif_hydrophobicity": round(motif_hydro, 3),
        "motif_charge": round(motif_charge, 3),
        "motif_helix_propensity": round(motif_helix, 3),
        "motif_surface_accessibility": round(motif_emini, 3),
        "motif_disorder_propensity": round(motif_disorder_prop, 3),
        "motif_iupred3": round(float(motif_iupred), 4),
        "motif_plddt": round(float(motif_plddt), 1),
        "flanking_hydrophobicity": round(flank_hydro, 3),
        "flanking_charge": round(flank_charge, 3),
        "flanking_iupred3": round(float(flank_iupred), 4),
        "flanking_plddt": round(float(flank_plddt), 1),
        "exposure_score": round(exposure_score, 4),
        "plddt_exposure": round(plddt_exposure, 3),
    }


# ═══════════════════════════════════════════════════════════════════════
# SLiM EXPOSURE ESTIMATES (without MD)
# ═══════════════════════════════════════════════════════════════════════

def estimate_slim_exposure(sequence, iupred_scores, plddt_scores, slim_hits):
    """
    Estimate SLiM exposure probabilities without MD simulation.
    Uses a composite of sequence-based accessibility, disorder prediction,
    and AlphaFold pLDDT confidence.

    This provides a static estimate; true exposure requires MD trajectories.
    """
    results = []
    for hit in slim_hits:
        ctx = compute_structural_context(sequence, iupred_scores, plddt_scores, hit["start"], hit["end"])

        # Exposure probability estimate
        # High disorder + low pLDDT + high surface accessibility → likely exposed
        p_exposed = ctx["exposure_score"]

        # Classify
        if p_exposed > 0.65:
            exposure_class = "likely_exposed"
        elif p_exposed > 0.45:
            exposure_class = "moderately_exposed"
        elif p_exposed > 0.3:
            exposure_class = "partially_buried"
        else:
            exposure_class = "likely_buried"

        result = {
            "protein": hit["protein"],
            "motif_class": hit["motif_class"],
            "start": hit["start"],
            "end": hit["end"],
            "matched_sequence": hit.get("matched_sequence", ""),
            "in_idr": hit.get("in_idr", False),
            "pathway": hit.get("pathway", ""),
            "p_exposed": round(p_exposed, 4),
            "exposure_class": exposure_class,
        }
        result.update(ctx)
        results.append(result)

    return results


# ═══════════════════════════════════════════════════════════════════════
# COMPETITIVE ACCESSIBILITY ANALYSIS
# ═══════════════════════════════════════════════════════════════════════

def analyze_competitive_accessibility(exposure_results, max_distance=30):
    """
    Identify SLiM pairs with potential competitive binding (overlapping or proximal SLiMs
    where exposure of one may occlude the other).

    Without MD, we estimate competition from:
    1. Spatial proximity (sequence distance < 30 residues)
    2. Overlapping binding surfaces (shared flanking residues)
    3. Opposite structural propensities (helix breaker vs helix former)
    """
    competitive_pairs = []

    # Group by protein
    by_protein = defaultdict(list)
    for exp in exposure_results:
        by_protein[exp["protein"]].append(exp)

    for protein, hits in by_protein.items():
        # Only consider IDR hits
        idr_hits = [h for h in hits if h["in_idr"]]

        for i in range(len(idr_hits)):
            for j in range(i + 1, len(idr_hits)):
                h1 = idr_hits[i]
                h2 = idr_hits[j]

                # Distance check
                dist = min(abs(h1["start"] - h2["end"]), abs(h2["start"] - h1["end"]))
                overlap = max(0, min(h1["end"], h2["end"]) - max(h1["start"], h2["start"]))

                if dist > max_distance and overlap == 0:
                    continue

                # Competition indicators
                # 1. Structural competition: different binding modes
                helix_diff = abs(h1["motif_helix_propensity"] - h2["motif_helix_propensity"])

                # 2. Charge complementarity (opposite charges → both attract same partner region)
                charge_product = h1["motif_charge"] * h2["motif_charge"]

                # 3. Exposure anti-correlation estimate
                exposure_diff = abs(h1["p_exposed"] - h2["p_exposed"])

                # Competition score
                competition_score = 0
                if overlap > 0:
                    competition_score += 0.4 * (overlap / max(h1["end"] - h1["start"] + 1, h2["end"] - h2["start"] + 1))
                else:
                    competition_score += 0.2 * max(0, 1 - dist / max_distance)

                competition_score += 0.2 * helix_diff
                competition_score += 0.2 * max(0, -charge_product)  # opposite charges compete
                competition_score += 0.2 * (1 - exposure_diff)  # similar exposure → steric clash

                # Classify
                if competition_score > 0.5:
                    comp_class = "strong_competition"
                elif competition_score > 0.3:
                    comp_class = "moderate_competition"
                elif competition_score > 0.15:
                    comp_class = "weak_competition"
                else:
                    comp_class = "independent"

                # Switch potential
                # Overlapping SLiMs with different targets → potential regulatory switch
                targets1 = set(h1.get("targets", "").split(",")) if h1.get("targets") else set()
                targets2 = set(h2.get("targets", "").split(",")) if h2.get("targets") else set()
                shared_targets = targets1 & targets2
                different_targets = (targets1 | targets2) - shared_targets

                different_pathways = h1["pathway"] != h2["pathway"]
                is_switch_candidate = (
                    competition_score > 0.35 and
                    h1["motif_class"] != h2["motif_class"] and
                    (len(different_targets) > 0 or different_pathways) and
                    (overlap > 0 or dist <= 10)
                )

                pair = {
                    "protein": protein,
                    "slim1_class": h1["motif_class"],
                    "slim1_pos": f"{h1['start']}-{h1['end']}",
                    "slim1_pathway": h1["pathway"],
                    "slim2_class": h2["motif_class"],
                    "slim2_pos": f"{h2['start']}-{h2['end']}",
                    "slim2_pathway": h2["pathway"],
                    "distance": dist,
                    "overlap": overlap,
                    "competition_score": round(competition_score, 4),
                    "competition_class": comp_class,
                    "is_switch_candidate": is_switch_candidate,
                    "helix_diff": round(helix_diff, 3),
                    "charge_product": round(charge_product, 3),
                    "exposure_diff": round(exposure_diff, 3),
                }
                competitive_pairs.append(pair)

    return competitive_pairs


# ═══════════════════════════════════════════════════════════════════════
# PLOTTING
# ═══════════════════════════════════════════════════════════════════════

def plot_structural_profiles(protein_name, sequence, iupred_scores, plddt_scores, idr_boundaries, exposure_results):
    """Plot structural context for one protein."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(sequence)
    positions = np.arange(1, n + 1)

    fig, axes = plt.subplots(4, 1, figsize=(14, 14), sharex=True)
    fig.suptitle(f"{protein_name} - Structural Context", fontsize=14, fontweight="bold")

    # Panel 1: IUPred3 + pLDDT
    ax = axes[0]
    ax.plot(positions[:len(iupred_scores)], iupred_scores, label="IUPred3", color="#E74C3C", linewidth=1)
    if plddt_scores is not None and len(plddt_scores) > 0:
        ax2 = ax.twinx()
        ax2.plot(positions[:len(plddt_scores)], plddt_scores, label="pLDDT", color="#3498DB", linewidth=1, alpha=0.7)
        ax2.set_ylabel("pLDDT", color="#3498DB")
        ax2.set_ylim(0, 100)
        ax2.tick_params(axis='y', labelcolor='#3498DB')
    for s, e in idr_boundaries:
        ax.axvspan(s, e, alpha=0.15, color="#F39C12")
    ax.set_ylabel("Disorder", color="#E74C3C")
    ax.set_ylim(0, 1.05)
    ax.set_title("Disorder + AlphaFold Confidence")

    # Panel 2: Hydrophobicity profile
    ax = axes[1]
    hydro_vals = [HYDRO.get(aa, 0) for aa in sequence]
    # Sliding window
    window = 9
    half_w = window // 2
    smoothed_hydro = []
    for i in range(n):
        s = max(0, i - half_w)
        e = min(n, i + half_w + 1)
        smoothed_hydro.append(np.mean(hydro_vals[s:e]))
    ax.fill_between(positions, smoothed_hydro, alpha=0.5, color="#2ECC71")
    ax.plot(positions, smoothed_hydro, color="#27AE60", linewidth=1)
    ax.axhline(0, color="gray", linewidth=0.5)
    for s, e in idr_boundaries:
        ax.axvspan(s, e, alpha=0.15, color="#F39C12")
    ax.set_ylabel("Hydrophobicity")
    ax.set_title("Kyte-Doolittle Hydropathy (window=9)")

    # Panel 3: Exposure estimates for SLiMs
    ax = axes[2]
    prot_exp = [e for e in exposure_results if e["protein"] == protein_name and e["in_idr"]]
    if prot_exp:
        for exp in prot_exp:
            mid = (exp["start"] + exp["end"]) / 2
            color = {"likely_exposed": "#2ECC71", "moderately_exposed": "#F39C12",
                     "partially_buried": "#E67E22", "likely_buried": "#E74C3C"}.get(exp["exposure_class"], "#BDC3C7")
            ax.bar(mid, exp["p_exposed"], width=exp["end"] - exp["start"] + 1, color=color, alpha=0.7, edgecolor="black", linewidth=0.3)
    for s, e in idr_boundaries:
        ax.axvspan(s, e, alpha=0.1, color="#F39C12")
    ax.set_ylabel("P(exposed)")
    ax.set_ylim(0, 1.05)
    ax.set_title("SLiM Exposure Probability (Green=exposed, Red=buried)")

    # Panel 4: Charge distribution
    ax = axes[3]
    charges = [CHARGE.get(aa, 0) for aa in sequence]
    colors_charge = ["#E74C3C" if c > 0 else "#3498DB" if c < 0 else "#BDC3C7" for c in charges]
    ax.bar(positions, charges, color=colors_charge, width=1.0, alpha=0.7)
    ax.axhline(0, color="gray", linewidth=0.5)
    for s, e in idr_boundaries:
        ax.axvspan(s, e, alpha=0.15, color="#F39C12")
    ax.set_ylabel("Charge")
    ax.set_xlabel("Residue Position")
    ax.set_title("Per-Residue Charge")

    plt.tight_layout()
    out_path = os.path.join(RESULTS, f"{protein_name}_structural_context.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Plot: {out_path}")
    return out_path


def plot_competition_summary(competitive_pairs):
    """Plot competitive accessibility summary."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Filter to target proteins
    target_pairs = [p for p in competitive_pairs if p["protein"] in ["ICP34.5", "US3", "pUL37"]]

    if not target_pairs:
        print("  No competitive pairs for target proteins")
        return None

    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    fig.suptitle("Competitive Accessibility Analysis (Target Proteins)", fontsize=14, fontweight="bold")

    # Panel 1: Competition score distribution
    ax = axes[0]
    for prot, color in [("ICP34.5", "#E74C3C"), ("US3", "#3498DB"), ("pUL37", "#2ECC71")]:
        scores = [p["competition_score"] for p in target_pairs if p["protein"] == prot]
        if scores:
            ax.hist(scores, bins=15, alpha=0.5, label=prot, color=color)
    ax.set_xlabel("Competition Score")
    ax.set_ylabel("Count")
    ax.set_title("Competition Score Distribution")
    ax.legend(fontsize=8)

    # Panel 2: Strong/moderate competition by protein
    ax = axes[1]
    from collections import Counter as Ctr
    strong_pairs = [p for p in target_pairs if p["competition_class"] in ["strong_competition", "moderate_competition"]]
    if strong_pairs:
        prots = [p["protein"] for p in strong_pairs]
        prot_counts = Ctr(prots)
        colors_map = {"ICP34.5": "#E74C3C", "US3": "#3498DB", "pUL37": "#2ECC71"}
        ax.bar(prot_counts.keys(), prot_counts.values(),
               color=[colors_map.get(k, "#BDC3C7") for k in prot_counts.keys()])
    ax.set_ylabel("Competitive Pairs")
    ax.set_title(f"Strong/Moderate Competition ({len(strong_pairs)} pairs)")

    # Panel 3: Top pathway pairs
    ax = axes[2]
    pathway_pairs = Ctr()
    for p in strong_pairs:
        pair = tuple(sorted([p["slim1_pathway"], p["slim2_pathway"]]))
        pathway_pairs[pair] += 1
    if pathway_pairs:
        top_pairs = pathway_pairs.most_common(10)
        labels = [f"{a}\nvs\n{b}" for (a, b), _ in top_pairs]
        counts = [c for _, c in top_pairs]
        ax.barh(range(len(labels)), counts, color="#9B59B6")
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel("Count")
    ax.set_title("Switch Pathway Pairs")

    plt.tight_layout()
    out_path = os.path.join(RESULTS, "competitive_accessibility_summary.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Competition plot: {out_path}")
    return out_path


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    # Load Phase 2 expanded hits
    slim_hits_path = os.path.join(PHASE2, "results", "slim_hits_expanded.csv")
    all_slim_hits = {}
    if os.path.exists(slim_hits_path):
        with open(slim_hits_path, encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                prot = row["protein"]
                if prot not in all_slim_hits:
                    all_slim_hits[prot] = []
                all_slim_hits[prot].append({
                    "protein": prot,
                    "motif_class": row["motif_class"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "matched_sequence": row.get("matched_sequence", ""),
                    "in_idr": row.get("in_idr", "False") == "True",
                    "pathway": row.get("pathway", ""),
                    "targets": row.get("targets", ""),
                    "category": row.get("category", ""),
                })
        print(f"Loaded {sum(len(v) for v in all_slim_hits.values())} SLiM hits from Phase 2")
    else:
        print("WARNING: No Phase 2 expanded hits found, using basic scan")

    # Load IDR boundaries
    with open(os.path.join(PHASE1, "results", "disorder_summary_v3.json")) as f:
        disorder_data = json.load(f)
    idr_map = {p["protein"]: p["idr_boundaries"] for p in disorder_data}

    all_exposure_results = []
    plddt_data = {}

    for pinfo in PROTEINS:
        name = pinfo["name"]
        accession = pinfo["accession"]
        category = pinfo["category"]
        print(f"\n{'='*60}")
        print(f"Structural analysis: {name} ({category})")
        print(f"{'='*60}")

        # Load sequence
        fasta_path = find_fasta(name)
        if not fasta_path:
            print(f"  ERROR: No FASTA found")
            continue
        sequence = load_fasta(fasta_path)

        # Load IUPred3
        iupred_scores = load_iupred3(name)
        if iupred_scores is None:
            print(f"  ERROR: No IUPred3 data")
            continue

        n = min(len(sequence), len(iupred_scores))
        sequence = sequence[:n]
        iupred_scores = iupred_scores[:n]

        # Fetch AlphaFold pLDDT
        plddt_scores = fetch_alphafold_plddt(accession, name)
        if plddt_scores is not None:
            plddt_data[name] = {
                "mean_plddt": round(float(np.mean(plddt_scores)), 1),
                "n_disordered_plddt": int(np.sum(plddt_scores < 50)),
                "plddt_disorder_frac": round(float(np.sum(plddt_scores < 50)) / len(plddt_scores), 3),
            }

        # SLiM exposure analysis
        slim_hits = all_slim_hits.get(name, [])
        if slim_hits:
            exposure_results = estimate_slim_exposure(sequence, iupred_scores, plddt_scores, slim_hits)
            all_exposure_results.extend(exposure_results)

            # Stats
            exposed = sum(1 for e in exposure_results if e["exposure_class"] == "likely_exposed")
            moderate = sum(1 for e in exposure_results if e["exposure_class"] == "moderately_exposed")
            buried = sum(1 for e in exposure_results if e["exposure_class"] in ["partially_buried", "likely_buried"])
            print(f"  SLiM exposure: {exposed} exposed, {moderate} moderate, {buried} buried (of {len(exposure_results)} total)")

        # Plot
        idr_bounds = idr_map.get(name, [])
        plot_structural_profiles(name, sequence, iupred_scores, plddt_scores, idr_bounds, all_exposure_results)

    # Competitive accessibility analysis
    print(f"\n{'='*60}")
    print("COMPETITIVE ACCESSIBILITY ANALYSIS")
    print(f"{'='*60}")
    competitive_pairs = analyze_competitive_accessibility(all_exposure_results)

    # Stats
    strong = sum(1 for p in competitive_pairs if p["competition_class"] == "strong_competition")
    moderate = sum(1 for p in competitive_pairs if p["competition_class"] == "moderate_competition")
    switches = sum(1 for p in competitive_pairs if p["is_switch_candidate"])
    print(f"  Total competitive pairs: {len(competitive_pairs)}")
    print(f"  Strong competition: {strong}")
    print(f"  Moderate competition: {moderate}")
    print(f"  Switch candidates: {switches}")

    # Target-specific switches
    for prot in ["ICP34.5", "US3", "pUL37"]:
        prot_switches = [p for p in competitive_pairs if p["protein"] == prot and p["is_switch_candidate"]]
        if prot_switches:
            print(f"\n  {prot} switch candidates:")
            for sw in sorted(prot_switches, key=lambda x: -x["competition_score"])[:5]:
                print(f"    {sw['slim1_class']} ({sw['slim1_pos']}, {sw['slim1_pathway']}) vs "
                      f"{sw['slim2_class']} ({sw['slim2_pos']}, {sw['slim2_pathway']}): "
                      f"score={sw['competition_score']:.3f}")

    # Plot competition
    plot_competition_summary(competitive_pairs)

    # Save outputs
    # 1. SLiM exposure results
    exp_path = os.path.join(RESULTS, "slim_exposure_estimates.csv")
    if all_exposure_results:
        keys = all_exposure_results[0].keys()
        with open(exp_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(all_exposure_results)
        print(f"\nExposure estimates: {exp_path} ({len(all_exposure_results)} entries)")

    # 2. Competitive pairs
    comp_path = os.path.join(RESULTS, "competitive_accessibility.csv")
    if competitive_pairs:
        keys = competitive_pairs[0].keys()
        with open(comp_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(competitive_pairs)
        print(f"Competitive pairs: {comp_path} ({len(competitive_pairs)} pairs)")

    # 3. pLDDT data
    plddt_path = os.path.join(RESULTS, "alphafold_plddt_summary.json")
    with open(plddt_path, "w") as f:
        json.dump(plddt_data, f, indent=2)
    print(f"AlphaFold pLDDT: {plddt_path}")

    # 4. Phase 3 summary
    summary = {
        "metadata": {
            "description": "Phase 3 Structural Analysis: AlphaFold pLDDT, SLiM exposure estimates, competitive accessibility",
            "methods": {
                "plddt": "AlphaFold2 predicted LDDT from EBI database",
                "exposure": "Composite of Emini accessibility, IUPred3 disorder, pLDDT confidence, and disorder propensity",
                "competition": "Sequence proximity (<30 residues), structural property differences, exposure correlation",
                "note": "Full MD (CALVADOS 2) and docking (FlexPepDock) require compute clusters"
            },
        },
        "alphafold_plddt": plddt_data,
        "exposure_summary": {
            "total_slims_scored": len(all_exposure_results),
            "likely_exposed": sum(1 for e in all_exposure_results if e["exposure_class"] == "likely_exposed"),
            "moderately_exposed": sum(1 for e in all_exposure_results if e["exposure_class"] == "moderately_exposed"),
            "partially_buried": sum(1 for e in all_exposure_results if e["exposure_class"] == "partially_buried"),
            "likely_buried": sum(1 for e in all_exposure_results if e["exposure_class"] == "likely_buried"),
        },
        "competition_summary": {
            "total_pairs": len(competitive_pairs),
            "strong_competition": strong,
            "moderate_competition": moderate,
            "switch_candidates": switches,
        },
    }

    summary_path = os.path.join(RESULTS, "phase3_structural_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
