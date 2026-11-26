#!/usr/bin/env python3
"""
Phase 1 Completion: Multi-method disorder, secondary structure propensities,
solvent accessibility estimates, charge patterning (CIDER), confidence tiers,
and per-residue annotation matrix for all 11 HSV-1 proteins.
"""

import os, sys, json, csv, math, re
import numpy as np
from collections import OrderedDict

BASE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(BASE)
RESULTS = os.path.join(BASE, "results")
IUPRED_DIR = os.path.join(BASE, "iupred3_data")
SEQ_DIR = os.path.join(BASE, "sequences")
os.makedirs(RESULTS, exist_ok=True)

# ── Protein metadata ──
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

# ═══════════════════════════════════════════════════════════════════════
# AMINO ACID PROPERTY SCALES
# ═══════════════════════════════════════════════════════════════════════

# Kyte-Doolittle hydrophobicity
HYDRO = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
}

# Charge at pH 7.0
CHARGE = {
    'R': 1.0, 'K': 1.0, 'H': 0.1,
    'D': -1.0, 'E': -1.0,
}

# Chou-Fasman helix propensity
CF_HELIX = {
    'A': 1.42, 'R': 0.98, 'N': 0.67, 'D': 1.01, 'C': 0.70,
    'Q': 1.11, 'E': 1.51, 'G': 0.57, 'H': 1.00, 'I': 1.08,
    'L': 1.21, 'K': 1.16, 'M': 1.45, 'F': 1.13, 'P': 0.57,
    'S': 0.77, 'T': 0.83, 'W': 1.08, 'Y': 0.69, 'V': 1.06,
}

# Chou-Fasman sheet propensity
CF_SHEET = {
    'A': 0.83, 'R': 0.93, 'N': 0.89, 'D': 0.54, 'C': 1.19,
    'Q': 1.10, 'E': 0.37, 'G': 0.75, 'H': 0.87, 'I': 1.60,
    'L': 1.30, 'K': 0.74, 'M': 1.05, 'F': 1.38, 'P': 0.55,
    'S': 0.75, 'T': 1.19, 'W': 1.37, 'Y': 1.47, 'V': 1.70,
}

# Chou-Fasman turn propensity
CF_TURN = {
    'A': 0.66, 'R': 0.95, 'N': 1.56, 'D': 1.46, 'C': 1.19,
    'Q': 0.98, 'E': 0.74, 'G': 1.56, 'H': 0.95, 'I': 0.47,
    'L': 0.59, 'K': 1.01, 'M': 0.60, 'F': 0.60, 'P': 1.52,
    'S': 1.43, 'T': 0.96, 'W': 0.96, 'Y': 1.14, 'V': 0.50,
}

# Janin interior-to-surface transfer (accessibility proxy; lower = more buried)
JANIN_ACC = {
    'A': 0.74, 'R': 0.64, 'N': 0.63, 'D': 0.62, 'C': 0.91,
    'Q': 0.62, 'E': 0.62, 'G': 0.72, 'H': 0.78, 'I': 0.88,
    'L': 0.85, 'K': 0.52, 'M': 0.85, 'F': 0.88, 'P': 0.64,
    'S': 0.66, 'T': 0.70, 'W': 0.85, 'Y': 0.76, 'V': 0.86,
}

# Emini surface accessibility scale (higher = more surface-exposed)
EMINI = {
    'A': 0.815, 'R': 1.475, 'N': 1.296, 'D': 1.283, 'C': 0.394,
    'Q': 1.348, 'E': 1.445, 'G': 0.714, 'H': 1.180, 'I': 0.603,
    'L': 0.603, 'K': 1.596, 'M': 0.714, 'F': 0.603, 'P': 1.236,
    'S': 1.115, 'T': 1.018, 'W': 0.603, 'Y': 0.714, 'V': 0.603,
}

# Disorder-promoting amino acids (Dunker et al. scale)
DISORDER_PROPENSITY = {
    'A': 0.06, 'R': -0.18, 'N': 0.01, 'D': 0.19, 'C': -0.34,
    'Q': 0.12, 'E': 0.30, 'G': 0.17, 'H': -0.01, 'I': -0.49,
    'L': -0.33, 'K': 0.09, 'M': -0.19, 'F': -0.42, 'P': 0.41,
    'S': 0.14, 'T': -0.05, 'W': -0.44, 'Y': -0.27, 'V': -0.46,
}


# ═══════════════════════════════════════════════════════════════════════
# SEQUENCE LOADING
# ═══════════════════════════════════════════════════════════════════════

def load_fasta(filepath):
    seq = []
    with open(filepath) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)


def find_fasta(protein_name):
    for sub in ["targets", "controls"]:
        path = os.path.join(SEQ_DIR, sub, f"{protein_name}.fasta")
        if os.path.exists(path):
            return path
    pn = protein_name.replace("_", "").replace("-", "").lower()
    for sub in ["targets", "controls"]:
        d = os.path.join(SEQ_DIR, sub)
        if not os.path.isdir(d):
            continue
        for fn in os.listdir(d):
            if fn.endswith(".fasta"):
                fn_norm = fn.replace("_", "").replace("-", "").replace(".fasta", "").lower()
                if pn == fn_norm or pn in fn_norm or fn_norm in pn:
                    return os.path.join(d, fn)
    return None


def load_iupred3(protein_name):
    pn = protein_name.replace("_", "").replace("-", "").lower()
    for fn in os.listdir(IUPRED_DIR):
        fn_norm = fn.replace("_", "").replace("-", "").lower()
        if pn in fn_norm:
            scores = []
            residues = []
            with open(os.path.join(IUPRED_DIR, fn)) as f:
                for line in f:
                    if line.startswith("#") or line.strip() == "":
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        residues.append(parts[1])
                        scores.append(float(parts[2]))
            return residues, scores
    return None, None


# ═══════════════════════════════════════════════════════════════════════
# METHOD 2: COMPOSITION-BASED DISORDER PREDICTION
# ═══════════════════════════════════════════════════════════════════════

def composition_disorder(sequence, window=21):
    """Predict disorder from amino acid composition using Dunker propensity scale."""
    n = len(sequence)
    scores = np.zeros(n)
    half_w = window // 2
    for i in range(n):
        start = max(0, i - half_w)
        end = min(n, i + half_w + 1)
        win_seq = sequence[start:end]
        props = [DISORDER_PROPENSITY.get(aa, 0.0) for aa in win_seq]
        scores[i] = np.mean(props)
    # Normalize to 0-1 range
    min_s, max_s = scores.min(), scores.max()
    if max_s > min_s:
        scores = (scores - min_s) / (max_s - min_s)
    return scores


# ═══════════════════════════════════════════════════════════════════════
# METHOD 3: SEQUENCE COMPLEXITY (LOW COMPLEXITY = DISORDER)
# ═══════════════════════════════════════════════════════════════════════

def complexity_disorder(sequence, window=21):
    """Shannon entropy-based disorder: low complexity windows → higher disorder."""
    n = len(sequence)
    scores = np.zeros(n)
    half_w = window // 2
    for i in range(n):
        start = max(0, i - half_w)
        end = min(n, i + half_w + 1)
        win_seq = sequence[start:end]
        # Shannon entropy
        aa_counts = {}
        for aa in win_seq:
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
        L = len(win_seq)
        entropy = 0
        for c in aa_counts.values():
            p = c / L
            if p > 0:
                entropy -= p * math.log2(p)
        # Max entropy for 20 AAs = log2(20) ≈ 4.32
        norm_entropy = entropy / 4.322
        # Low complexity → high disorder propensity
        # Combine with charge/hydro bias
        charges = sum(1 for aa in win_seq if aa in "RKDEH")
        charge_frac = charges / L
        hydro_mean = np.mean([HYDRO.get(aa, 0) for aa in win_seq])
        # Disorder score: low complexity + high charge + low hydrophobicity
        dis = (1 - norm_entropy) * 0.3 + charge_frac * 0.3 + max(0, (1 - (hydro_mean + 4.5) / 9.0)) * 0.4
        scores[i] = dis
    # Normalize
    min_s, max_s = scores.min(), scores.max()
    if max_s > min_s:
        scores = (scores - min_s) / (max_s - min_s)
    return scores


# ═══════════════════════════════════════════════════════════════════════
# CONSENSUS DISORDER
# ═══════════════════════════════════════════════════════════════════════

def consensus_disorder(iupred_scores, comp_scores, complex_scores, threshold=0.5):
    """Consensus: residue is disordered if ≥2/3 methods agree (score > threshold)."""
    n = len(iupred_scores)
    consensus = np.zeros(n)
    for i in range(n):
        votes = 0
        if iupred_scores[i] > threshold:
            votes += 1
        if comp_scores[i] > threshold:
            votes += 1
        if complex_scores[i] > threshold:
            votes += 1
        consensus[i] = 1.0 if votes >= 2 else 0.0
    # Extract contiguous IDRs of >= 20 residues
    idrs = []
    start = None
    for i in range(n):
        if consensus[i] == 1.0 and start is None:
            start = i
        elif consensus[i] == 0.0 and start is not None:
            if i - start >= 20:
                idrs.append((start + 1, i))  # 1-indexed
            start = None
    if start is not None and n - start >= 20:
        idrs.append((start + 1, n))
    return consensus, idrs


# ═══════════════════════════════════════════════════════════════════════
# SECONDARY STRUCTURE PROPENSITIES (Chou-Fasman windowed)
# ═══════════════════════════════════════════════════════════════════════

def secondary_structure_propensity(sequence, window=9):
    """Compute per-residue helix, sheet, turn propensities (Chou-Fasman, windowed)."""
    n = len(sequence)
    helix = np.zeros(n)
    sheet = np.zeros(n)
    turn = np.zeros(n)
    half_w = window // 2
    for i in range(n):
        start = max(0, i - half_w)
        end = min(n, i + half_w + 1)
        win = sequence[start:end]
        helix[i] = np.mean([CF_HELIX.get(aa, 1.0) for aa in win])
        sheet[i] = np.mean([CF_SHEET.get(aa, 1.0) for aa in win])
        turn[i] = np.mean([CF_TURN.get(aa, 1.0) for aa in win])
    # Assign SS state: highest propensity wins, with coil as default
    ss_state = []
    for i in range(n):
        h, s, t = helix[i], sheet[i], turn[i]
        if h > 1.03 and h >= s and h >= t:
            ss_state.append("H")
        elif s > 1.05 and s >= h and s >= t:
            ss_state.append("E")
        elif t > 1.0:
            ss_state.append("T")
        else:
            ss_state.append("C")
    return helix, sheet, turn, ss_state


# ═══════════════════════════════════════════════════════════════════════
# SOLVENT ACCESSIBILITY ESTIMATE
# ═══════════════════════════════════════════════════════════════════════

def accessibility_estimate(sequence, iupred_scores, window=9):
    """Estimate per-residue surface accessibility from Emini scale + disorder context."""
    n = len(sequence)
    emini_scores = np.zeros(n)
    janin_scores = np.zeros(n)
    half_w = window // 2
    for i in range(n):
        start = max(0, i - half_w)
        end = min(n, i + half_w + 1)
        win = sequence[start:end]
        emini_scores[i] = np.mean([EMINI.get(aa, 0.8) for aa in win])
        janin_scores[i] = np.mean([JANIN_ACC.get(aa, 0.7) for aa in win])
    # Composite accessibility: Emini (surface) weighted by disorder (IDR → more exposed)
    combined = np.zeros(n)
    for i in range(n):
        # In disordered regions, residues are more exposed
        disorder_boost = iupred_scores[i] * 0.4 if i < len(iupred_scores) else 0.2
        combined[i] = emini_scores[i] * 0.4 + (1 - janin_scores[i]) * 0.3 + disorder_boost * 0.3
    # Normalize to 0-1
    min_c, max_c = combined.min(), combined.max()
    if max_c > min_c:
        combined = (combined - min_c) / (max_c - min_c)
    return emini_scores, janin_scores, combined


# ═══════════════════════════════════════════════════════════════════════
# CHARGE PATTERNING (CIDER-like metrics)
# ═══════════════════════════════════════════════════════════════════════

def charge_patterning(sequence):
    """Compute CIDER-like charge patterning metrics."""
    n = len(sequence)
    # Per-residue charge
    charges = np.array([CHARGE.get(aa, 0.0) for aa in sequence])

    # Global metrics
    n_pos = sum(1 for aa in sequence if aa in "RK")
    n_neg = sum(1 for aa in sequence if aa in "DE")
    fcr = (n_pos + n_neg) / n  # fraction of charged residues
    ncpr = (n_pos - n_neg) / n  # net charge per residue

    # Kappa (charge patterning parameter, Das & Pappu 2013)
    # Simplified: measures how well-mixed positive and negative charges are
    # kappa = 0 → well-mixed, kappa = 1 → segregated
    if n_pos + n_neg < 4:
        kappa = 0.0
        sigma_max = 0.0
    else:
        # Sliding window sigma
        window = 5
        sigma_vals = []
        for i in range(n - window + 1):
            win_charges = charges[i:i+window]
            f_plus = sum(1 for c in win_charges if c > 0) / window
            f_minus = sum(1 for c in win_charges if c < 0) / window
            sigma = (f_plus - f_minus) ** 2
            sigma_vals.append(sigma)
        sigma_mean = np.mean(sigma_vals) if sigma_vals else 0
        sigma_max = ((n_pos/n) - (n_neg/n)) ** 2 if n > 0 else 0
        # For well-mixed: sigma ~ sigma_max; for segregated: sigma >> sigma_max
        # Approximate kappa
        if sigma_max > 0:
            kappa = min(1.0, sigma_mean / (sigma_max + 0.01))
        else:
            kappa = 0.0

    # Hydropathy (mean Kyte-Doolittle)
    mean_hydro = np.mean([HYDRO.get(aa, 0) for aa in sequence])

    # Uversky plot position (mean hydropathy vs |mean net charge|)
    # Boundary: <H> = 2.785 * |<R>| - 1.151
    abs_ncpr = abs(ncpr)
    uversky_boundary = 2.785 * abs_ncpr - 1.151
    # Normalize hydropathy to 0-1 for Uversky
    norm_hydro = (mean_hydro + 4.5) / 9.0  # scale from [-4.5, 4.5] to [0, 1]
    uversky_distance = norm_hydro - (2.785 * abs_ncpr - 1.151)
    # Negative distance → disordered side

    # Das-Pappu phase diagram classification
    if fcr < 0.25:
        if norm_hydro < 0.4:
            dp_class = "weak_polyampholyte"
        else:
            dp_class = "globule"
    elif fcr < 0.5:
        if abs(ncpr) < 0.1:
            dp_class = "strong_polyampholyte"
        else:
            dp_class = "polyelectrolyte"
    else:
        dp_class = "strong_polyelectrolyte"

    return {
        "fcr": round(fcr, 4),
        "ncpr": round(ncpr, 4),
        "kappa": round(kappa, 4),
        "mean_hydropathy": round(mean_hydro, 3),
        "norm_hydropathy": round(norm_hydro, 3),
        "n_positive": n_pos,
        "n_negative": n_neg,
        "uversky_distance": round(uversky_distance, 3),
        "das_pappu_class": dp_class,
        "per_residue_charge": charges.tolist(),
    }


def cider_analysis(sequence):
    """Full CIDER analysis using localCIDER if available, else manual."""
    try:
        from localcider.sequenceParameters import SequenceParameters
        sp = SequenceParameters(sequence)
        return {
            "kappa": round(sp.get_kappa(), 4),
            "fcr": round(sp.get_FCR(), 4),
            "ncpr": round(sp.get_NCPR(), 4),
            "mean_hydropathy": round(sp.get_mean_hydropathy(), 3),
            "delta": round(sp.get_delta(), 4),
            "deltaMax": round(sp.get_deltaMax(), 4),
            "source": "localCIDER",
        }
    except Exception as e:
        result = charge_patterning(sequence)
        result["source"] = "manual"
        return result


# ═══════════════════════════════════════════════════════════════════════
# CONFIDENCE TIERS
# ═══════════════════════════════════════════════════════════════════════

def assign_confidence_tiers(iupred_scores, comp_scores, complex_scores):
    """
    Assign confidence tiers per residue:
    Tier 1 (High): All 3 methods agree (all >0.5 or all ≤0.5)
    Tier 2 (Medium): 2/3 methods agree
    Tier 3 (Low): Only 1/3 methods predicts disorder
    """
    n = len(iupred_scores)
    tiers = []
    for i in range(n):
        votes_dis = 0
        if iupred_scores[i] > 0.5: votes_dis += 1
        if comp_scores[i] > 0.5: votes_dis += 1
        if complex_scores[i] > 0.5: votes_dis += 1

        if votes_dis == 3:
            tiers.append(("Tier1_High", "disordered"))
        elif votes_dis == 0:
            tiers.append(("Tier1_High", "ordered"))
        elif votes_dis == 2:
            tiers.append(("Tier2_Medium", "disordered"))
        else:
            tiers.append(("Tier3_Low", "disordered"))
    return tiers


# ═══════════════════════════════════════════════════════════════════════
# PER-RESIDUE ANNOTATION MATRIX
# ═══════════════════════════════════════════════════════════════════════

def build_annotation_matrix(protein_name, sequence, iupred_scores, comp_scores,
                            complex_scores, consensus, tiers, helix, sheet, turn,
                            ss_state, accessibility, charges):
    """Build the complete per-residue annotation matrix."""
    n = len(sequence)
    rows = []
    for i in range(n):
        row = {
            "protein": protein_name,
            "position": i + 1,
            "residue": sequence[i],
            "iupred3_score": round(iupred_scores[i], 4) if i < len(iupred_scores) else None,
            "composition_disorder": round(comp_scores[i], 4),
            "complexity_disorder": round(complex_scores[i], 4),
            "consensus_disorder": int(consensus[i]),
            "confidence_tier": tiers[i][0],
            "disorder_state": tiers[i][1],
            "helix_propensity": round(helix[i], 3),
            "sheet_propensity": round(sheet[i], 3),
            "turn_propensity": round(turn[i], 3),
            "ss_prediction": ss_state[i],
            "accessibility_score": round(accessibility[i], 4),
            "charge": round(charges[i], 2),
        }
        rows.append(row)
    return rows


# ═══════════════════════════════════════════════════════════════════════
# PLOTTING
# ═══════════════════════════════════════════════════════════════════════

def plot_protein_annotation(protein_name, sequence, iupred_scores, comp_scores,
                            complex_scores, consensus, helix, sheet, turn,
                            accessibility, charges, tiers, idr_boundaries):
    """Generate multi-panel annotation plot for one protein."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    n = len(sequence)
    positions = np.arange(1, n + 1)

    fig, axes = plt.subplots(5, 1, figsize=(14, 16), sharex=True,
                              gridspec_kw={"height_ratios": [2, 1.5, 1.5, 1.5, 1]})
    fig.suptitle(f"{protein_name} — Per-Residue Annotation", fontsize=14, fontweight="bold")

    # Panel 1: Multi-method disorder
    ax = axes[0]
    ax.plot(positions[:len(iupred_scores)], iupred_scores, label="IUPred3", color="#E74C3C", alpha=0.8, linewidth=1)
    ax.plot(positions, comp_scores, label="Composition", color="#3498DB", alpha=0.8, linewidth=1)
    ax.plot(positions, complex_scores, label="Complexity", color="#2ECC71", alpha=0.8, linewidth=1)
    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5)
    # Shade consensus IDRs
    for s, e in idr_boundaries:
        ax.axvspan(s, e, alpha=0.15, color="#F39C12")
    ax.set_ylabel("Disorder Score")
    ax.set_ylim(0, 1.05)
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("Multi-Method Disorder Prediction")

    # Panel 2: Secondary structure propensities
    ax = axes[1]
    ax.plot(positions, helix, label="Helix (Pα)", color="#E74C3C", alpha=0.8, linewidth=1)
    ax.plot(positions, sheet, label="Sheet (Pβ)", color="#3498DB", alpha=0.8, linewidth=1)
    ax.plot(positions, turn, label="Turn (Pt)", color="#2ECC71", alpha=0.8, linewidth=1)
    ax.axhline(1.0, color="gray", linestyle="--", alpha=0.5)
    for s, e in idr_boundaries:
        ax.axvspan(s, e, alpha=0.15, color="#F39C12")
    ax.set_ylabel("Propensity")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("Chou-Fasman Secondary Structure Propensities")

    # Panel 3: Accessibility
    ax = axes[2]
    ax.fill_between(positions, accessibility, alpha=0.5, color="#9B59B6")
    ax.plot(positions, accessibility, color="#8E44AD", linewidth=1)
    for s, e in idr_boundaries:
        ax.axvspan(s, e, alpha=0.15, color="#F39C12")
    ax.set_ylabel("Accessibility")
    ax.set_ylim(0, 1.05)
    ax.set_title("Estimated Solvent Accessibility")

    # Panel 4: Charge distribution
    ax = axes[3]
    colors_charge = ["#E74C3C" if c > 0 else "#3498DB" if c < 0 else "#BDC3C7" for c in charges]
    ax.bar(positions, charges, color=colors_charge, width=1.0, alpha=0.7)
    ax.axhline(0, color="gray", linewidth=0.5)
    for s, e in idr_boundaries:
        ax.axvspan(s, e, alpha=0.15, color="#F39C12")
    ax.set_ylabel("Charge")
    ax.set_title("Per-Residue Charge")

    # Panel 5: Confidence tier track
    ax = axes[4]
    tier_colors = {"Tier1_High": "#2ECC71", "Tier2_Medium": "#F39C12", "Tier3_Low": "#E74C3C"}
    for i, (tier, state) in enumerate(tiers):
        c = tier_colors.get(tier, "#BDC3C7")
        if state == "ordered":
            c = "#BDC3C7"
        ax.add_patch(Rectangle((i + 0.5, 0), 1, 1, facecolor=c, edgecolor="none"))
    ax.set_xlim(0.5, n + 0.5)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel("Residue Position")
    ax.set_title("Confidence Tier (Green=Tier1, Orange=Tier2, Red=Tier3, Gray=Ordered)")

    plt.tight_layout()
    out_path = os.path.join(RESULTS, f"{protein_name}_annotation.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Plot saved: {out_path}")
    return out_path


def plot_cider_summary(all_cider):
    """Plot CIDER charge patterning summary across all proteins."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    fig.suptitle("Charge Patterning Analysis (CIDER Metrics)", fontsize=14, fontweight="bold")

    # Normalize keys: localCIDER uses different names
    for name, d in all_cider.items():
        if "norm_hydropathy" not in d:
            mh = d.get("mean_hydropathy", 0)
            d["norm_hydropathy"] = (mh + 4.5) / 9.0  # normalize to 0-1
    targets = {k: v for k, v in all_cider.items() if v.get("category") == "target"}
    controls = {k: v for k, v in all_cider.items() if v.get("category") == "control"}

    # Panel 1: FCR vs NCPR (Das-Pappu diagram)
    ax = axes[0]
    for name, d in targets.items():
        ax.scatter(d["ncpr"], d["fcr"], c="#E74C3C", s=80, zorder=5, edgecolors="black", linewidth=0.5)
        ax.annotate(name, (d["ncpr"], d["fcr"]), fontsize=7, ha="left", va="bottom")
    for name, d in controls.items():
        ax.scatter(d["ncpr"], d["fcr"], c="#3498DB", s=60, zorder=4, edgecolors="black", linewidth=0.5)
        ax.annotate(name, (d["ncpr"], d["fcr"]), fontsize=6, ha="left", va="bottom")
    ax.set_xlabel("NCPR")
    ax.set_ylabel("FCR")
    ax.set_title("Das-Pappu Diagram")
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axvline(0, color="gray", linewidth=0.5)

    # Panel 2: Kappa comparison
    ax = axes[1]
    names_t = list(targets.keys())
    names_c = list(controls.keys())
    all_names = names_t + names_c
    kappas = [targets[n]["kappa"] for n in names_t] + [controls[n]["kappa"] for n in names_c]
    colors = ["#E74C3C"] * len(names_t) + ["#3498DB"] * len(names_c)
    bars = ax.barh(range(len(all_names)), kappas, color=colors)
    ax.set_yticks(range(len(all_names)))
    ax.set_yticklabels(all_names, fontsize=8)
    ax.set_xlabel("κ (charge patterning)")
    ax.set_title("Charge Segregation (κ)")

    # Panel 3: Hydropathy vs |NCPR| (Uversky plot)
    ax = axes[2]
    # Boundary line
    ncpr_range = np.linspace(0, 0.5, 100)
    boundary = 2.785 * ncpr_range - 1.151
    ax.plot(ncpr_range, boundary, "k--", alpha=0.5, label="Disorder boundary")
    ax.fill_between(ncpr_range, 0, boundary, alpha=0.1, color="orange")
    ax.text(0.05, 0.1, "Disordered", fontsize=9, color="orange", alpha=0.7)
    ax.text(0.05, 0.6, "Ordered", fontsize=9, color="blue", alpha=0.7)
    for name, d in targets.items():
        ax.scatter(abs(d["ncpr"]), d["norm_hydropathy"], c="#E74C3C", s=80, zorder=5,
                   edgecolors="black", linewidth=0.5)
        ax.annotate(name, (abs(d["ncpr"]), d["norm_hydropathy"]), fontsize=7)
    for name, d in controls.items():
        ax.scatter(abs(d["ncpr"]), d["norm_hydropathy"], c="#3498DB", s=60, zorder=4,
                   edgecolors="black", linewidth=0.5)
        ax.annotate(name, (abs(d["ncpr"]), d["norm_hydropathy"]), fontsize=6)
    ax.set_xlabel("|NCPR|")
    ax.set_ylabel("Mean Hydropathy (normalized)")
    ax.set_title("Uversky Plot")
    ax.set_xlim(-0.02, 0.5)
    ax.set_ylim(0, 1)

    plt.tight_layout()
    out_path = os.path.join(RESULTS, "cider_charge_patterning.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"CIDER summary plot saved: {out_path}")
    return out_path


# ═══════════════════════════════════════════════════════════════════════
# MAIN PIPELINE
# ═══════════════════════════════════════════════════════════════════════

def main():
    # Load existing Phase 1 summary for IDR boundaries
    with open(os.path.join(RESULTS, "disorder_summary_v3.json")) as f:
        existing = json.load(f)
    existing_map = {p["protein"]: p for p in existing}

    all_annotations = []
    all_cider = {}
    protein_summaries = []

    for pinfo in PROTEINS:
        name = pinfo["name"]
        accession = pinfo["accession"]
        category = pinfo["category"]
        print(f"\n{'='*60}")
        print(f"Processing: {name} ({category})")
        print(f"{'='*60}")

        # Load sequence
        fasta_path = find_fasta(name)
        if fasta_path is None:
            print(f"  WARNING: No FASTA found for {name}, skipping")
            continue
        sequence = load_fasta(fasta_path)
        print(f"  Sequence length: {len(sequence)}")

        # Load IUPred3 scores
        residues, iupred_scores = load_iupred3(name)
        if iupred_scores is None:
            print(f"  WARNING: No IUPred3 data for {name}, skipping")
            continue
        iupred_scores = np.array(iupred_scores)
        print(f"  IUPred3 scores: {len(iupred_scores)} residues, mean={np.mean(iupred_scores):.3f}")

        n = min(len(sequence), len(iupred_scores))
        sequence = sequence[:n]
        iupred_scores = iupred_scores[:n]

        # Method 2: Composition-based disorder
        comp_scores = composition_disorder(sequence)
        print(f"  Composition disorder: mean={np.mean(comp_scores):.3f}")

        # Method 3: Complexity-based disorder
        complex_scores = complexity_disorder(sequence)
        print(f"  Complexity disorder: mean={np.mean(complex_scores):.3f}")

        # Consensus disorder
        consensus, consensus_idrs = consensus_disorder(iupred_scores, comp_scores, complex_scores)
        n_dis = int(np.sum(consensus))
        print(f"  Consensus disorder: {n_dis}/{n} residues ({100*n_dis/n:.1f}%), {len(consensus_idrs)} IDRs")

        # Confidence tiers
        tiers = assign_confidence_tiers(iupred_scores, comp_scores, complex_scores)
        tier_counts = {}
        for t, s in tiers:
            key = f"{t}_{s}"
            tier_counts[key] = tier_counts.get(key, 0) + 1
        print(f"  Confidence tiers: {tier_counts}")

        # Secondary structure propensities
        helix, sheet, turn, ss_state = secondary_structure_propensity(sequence)
        ss_counts = {s: ss_state.count(s) for s in set(ss_state)}
        print(f"  SS prediction: {ss_counts}")

        # Solvent accessibility
        emini, janin, accessibility = accessibility_estimate(sequence, iupred_scores)
        print(f"  Accessibility: mean={np.mean(accessibility):.3f}")

        # Charge patterning (CIDER)
        charges = np.array([CHARGE.get(aa, 0.0) for aa in sequence])
        cider = cider_analysis(sequence)
        cider["category"] = category
        all_cider[name] = cider
        print(f"  CIDER: FCR={cider['fcr']}, NCPR={cider['ncpr']}, kappa={cider['kappa']}")

        # IDR boundaries from existing analysis
        idr_bounds = existing_map.get(name, {}).get("idr_boundaries", [])

        # Build annotation matrix
        matrix = build_annotation_matrix(
            name, sequence, iupred_scores, comp_scores, complex_scores,
            consensus, tiers, helix, sheet, turn, ss_state, accessibility, charges
        )
        all_annotations.extend(matrix)

        # Generate plot
        plot_protein_annotation(
            name, sequence, iupred_scores, comp_scores, complex_scores,
            consensus, helix, sheet, turn, accessibility, charges, tiers, idr_bounds
        )

        # IDR-level CIDER metrics
        idr_cider_data = []
        for s, e in idr_bounds:
            if s <= n and e <= n:
                idr_seq = sequence[s-1:e]
                if len(idr_seq) >= 10:
                    idr_c = charge_patterning(idr_seq)
                    idr_cider_data.append({
                        "region": f"{s}-{e}",
                        "length": len(idr_seq),
                        "fcr": idr_c["fcr"],
                        "ncpr": idr_c["ncpr"],
                        "mean_hydropathy": idr_c["mean_hydropathy"],
                        "das_pappu_class": idr_c["das_pappu_class"],
                    })

        # Protein summary
        summary = {
            "protein": name,
            "accession": accession,
            "category": category,
            "length": n,
            "iupred3_mean": round(float(np.mean(iupred_scores)), 4),
            "composition_mean": round(float(np.mean(comp_scores)), 4),
            "complexity_mean": round(float(np.mean(complex_scores)), 4),
            "consensus_idr_fraction": round(n_dis / n, 4),
            "consensus_idr_count": len(consensus_idrs),
            "consensus_idr_boundaries": consensus_idrs,
            "original_idr_boundaries": idr_bounds,
            "confidence_tiers": {
                "tier1_disordered": sum(1 for t, s in tiers if t == "Tier1_High" and s == "disordered"),
                "tier1_ordered": sum(1 for t, s in tiers if t == "Tier1_High" and s == "ordered"),
                "tier2_medium": sum(1 for t, s in tiers if t == "Tier2_Medium"),
                "tier3_low": sum(1 for t, s in tiers if t == "Tier3_Low"),
            },
            "secondary_structure": {
                "helix_residues": ss_state.count("H"),
                "sheet_residues": ss_state.count("E"),
                "turn_residues": ss_state.count("T"),
                "coil_residues": ss_state.count("C"),
                "helix_fraction": round(ss_state.count("H") / n, 3),
                "sheet_fraction": round(ss_state.count("E") / n, 3),
                "turn_fraction": round(ss_state.count("T") / n, 3),
                "coil_fraction": round(ss_state.count("C") / n, 3),
            },
            "accessibility_mean": round(float(np.mean(accessibility)), 4),
            "cider_metrics": {k: v for k, v in cider.items() if k not in ["per_residue_charge", "category"]},
            "idr_cider": idr_cider_data,
        }
        protein_summaries.append(summary)

    # Save all outputs
    # 1. Per-residue annotation matrix (CSV)
    csv_path = os.path.join(RESULTS, "per_residue_annotation_matrix.csv")
    if all_annotations:
        keys = all_annotations[0].keys()
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(all_annotations)
        print(f"\nAnnotation matrix: {csv_path} ({len(all_annotations)} rows)")

    # 2. CIDER summary
    cider_path = os.path.join(RESULTS, "cider_charge_patterning.json")
    cider_out = {k: {kk: vv for kk, vv in v.items() if kk != "per_residue_charge"}
                 for k, v in all_cider.items()}
    with open(cider_path, "w") as f:
        json.dump(cider_out, f, indent=2)
    print(f"CIDER metrics: {cider_path}")

    # 3. Full Phase 1 summary
    summary_path = os.path.join(RESULTS, "phase1_complete_summary.json")
    full_summary = {
        "metadata": {
            "description": "Phase 1 Complete: Multi-method disorder, SS propensities, accessibility, charge patterning, confidence tiers",
            "methods": {
                "disorder_method_1": "IUPred3 long disorder (threshold 0.5)",
                "disorder_method_2": "Composition-based (Dunker propensity scale, window=21)",
                "disorder_method_3": "Complexity-based (Shannon entropy + charge + hydropathy, window=21)",
                "consensus": "≥2/3 methods agree, contiguous ≥20 residues",
                "secondary_structure": "Chou-Fasman propensities (helix/sheet/turn), window=9",
                "accessibility": "Emini surface + Janin interior scales + disorder context",
                "charge_patterning": "localCIDER (FCR, NCPR, kappa, hydropathy, Das-Pappu class)",
            },
            "confidence_tiers": {
                "Tier1_High": "All 3 methods agree",
                "Tier2_Medium": "2/3 methods agree",
                "Tier3_Low": "Only 1/3 methods predicts disorder",
            },
            "n_proteins": len(protein_summaries),
        },
        "proteins": protein_summaries,
    }
    with open(summary_path, "w") as f:
        json.dump(full_summary, f, indent=2)
    print(f"Full summary: {summary_path}")

    # 4. CIDER summary plot
    plot_cider_summary(all_cider)

    # Print summary table
    print(f"\n{'='*80}")
    print("PHASE 1 COMPLETION SUMMARY")
    print(f"{'='*80}")
    print(f"{'Protein':<18} {'Cat':<8} {'Len':>5} {'IUPred':>7} {'Comp':>7} {'Cmplx':>7} {'Cons%':>7} {'H%':>5} {'E%':>5} {'T%':>5} {'FCR':>6} {'NCPR':>7}")
    print("-" * 100)
    for s in protein_summaries:
        ss = s["secondary_structure"]
        ci = s["cider_metrics"]
        print(f"{s['protein']:<18} {s['category']:<8} {s['length']:>5} "
              f"{s['iupred3_mean']:>7.3f} {s['composition_mean']:>7.3f} {s['complexity_mean']:>7.3f} "
              f"{s['consensus_idr_fraction']*100:>6.1f}% "
              f"{ss['helix_fraction']*100:>4.1f} {ss['sheet_fraction']*100:>4.1f} {ss['turn_fraction']*100:>4.1f} "
              f"{ci['fcr']:>6.3f} {ci['ncpr']:>7.4f}")


if __name__ == "__main__":
    main()
