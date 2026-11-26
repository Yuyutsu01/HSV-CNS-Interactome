"""
Phase 5 — Competitive-Binding Switch Identification
Proposal gaps filled:
  1. Decision tree for switch calling (multi-criteria)
  2. Switch classification (immune evasion, neuronal hijacking, apoptosis modulation)
  3. Host-viral competitive binding analysis
  4. Switch confidence scoring
  5. Biological interpretation and mechanism annotation
  6. Summary statistics and publication plots
"""

import os, sys, json, math
import numpy as np
import pandas as pd
from collections import defaultdict, Counter
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

BASE_DIR = r"C:\Users\shiva\OneDrive\Desktop\HSV-CNS-Interactome"
PHASE2_DIR = os.path.join(BASE_DIR, "phase2")
PHASE3_DIR = os.path.join(BASE_DIR, "phase3")
PHASE4_DIR = os.path.join(BASE_DIR, "phase4")
PHASE5_DIR = os.path.join(BASE_DIR, "phase5")
RESULTS_DIR = os.path.join(PHASE5_DIR, "results")
PLOT_DIR = os.path.join(PHASE5_DIR, "plots")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)

SWITCH_CATEGORIES = {
    "immune_evasion": {
        "pathways": {"innate_immunity", "interferon_signaling", "NFkB_innate_immunity",
                     "TLR_signaling", "antigen_presentation", "complement"},
        "description": "Viral SLiM competes with host immune signaling motif",
    },
    "apoptosis_modulation": {
        "pathways": {"apoptosis", "ubiquitin_proteasome", "caspase", "autophagy"},
        "description": "Viral SLiM modulates cell death pathways",
    },
    "neuronal_hijacking": {
        "pathways": {"synaptic_adhesion", "calcium_signaling", "axonal_transport",
                     "synaptic_transmission", "neuronal_signaling"},
        "description": "Viral SLiM targets neuron-specific signaling",
    },
    "signaling_rewiring": {
        "pathways": {"MAPK_signaling", "Wnt_signaling", "PI3K_Akt", "CK2_signaling",
                     "cell_cycle", "cell_signaling"},
        "description": "Viral SLiM redirects host signaling cascades",
    },
    "transport_hijacking": {
        "pathways": {"nuclear_transport", "Golgi_transport", "endocytosis",
                     "retrograde_transport", "vesicular_transport"},
        "description": "Viral SLiM subverts intracellular transport",
    },
    "transcription_control": {
        "pathways": {"transcription", "chromatin_remodeling", "RNA_processing",
                     "translation", "mRNA_regulation"},
        "description": "Viral SLiM modulates gene expression",
    },
}

KNOWN_SWITCH_MECHANISMS = {
    ("ICP34.5", "LIG_PP1_RVXF_1"): {
        "mechanism": "ICP34.5 recruits PP1 to dephosphorylate eIF2alpha, blocking PKR-mediated translation shutoff",
        "reference": "He et al., PNAS 1997",
        "validated": True,
    },
    ("ICP34.5", "LIG_ATG8_1"): {
        "mechanism": "ICP34.5 binds Beclin-1 to inhibit autophagy, preventing viral clearance",
        "reference": "Orvedahl et al., Cell Host Microbe 2007",
        "validated": True,
    },
    ("US3", "MOD_PKA_2"): {
        "mechanism": "US3 kinase mimics PKA to phosphorylate BAD, blocking apoptosis",
        "reference": "Cartier et al., J Virol 2003",
        "validated": True,
    },
    ("pUL37", "LIG_TRAF6_1"): {
        "mechanism": "pUL37 deamidase targets TRAF6/TRAF3 to suppress NF-kB and IRF3 activation",
        "reference": "Zhang et al., Cell Host Microbe 2013",
        "validated": True,
    },
}


def classify_switch(pw1, pw2):
    """Classify a competitive pair into switch category based on pathway annotations."""
    pw1_lower = pw1.lower() if pw1 else ""
    pw2_lower = pw2.lower() if pw2 else ""

    for cat, info in SWITCH_CATEGORIES.items():
        cat_pws = {p.lower() for p in info["pathways"]}
        if pw1_lower in cat_pws or pw2_lower in cat_pws:
            return cat
    return "other"


def compute_switch_confidence(row, exposure_data, conservation_data):
    """Multi-criteria decision tree for switch confidence scoring."""
    score = 0
    criteria_met = []

    competition = row.get("competition_score", 0)
    if competition >= 0.5:
        score += 3
        criteria_met.append("strong_competition")
    elif competition >= 0.35:
        score += 2
        criteria_met.append("moderate_competition")
    elif competition >= 0.25:
        score += 1
        criteria_met.append("weak_competition")

    if row.get("overlap", 0) > 0:
        score += 2
        criteria_met.append("overlapping_motifs")

    distance = row.get("distance", 999)
    if distance <= 5:
        score += 2
        criteria_met.append("very_proximal")
    elif distance <= 15:
        score += 1
        criteria_met.append("proximal")

    slim1_class = row.get("slim1_class", "")
    slim2_class = row.get("slim2_class", "")
    if slim1_class != slim2_class:
        score += 1
        criteria_met.append("different_motif_classes")

    pw1 = str(row.get("slim1_pathway", ""))
    pw2 = str(row.get("slim2_pathway", ""))
    if pw1 != pw2 and pw1 and pw2:
        score += 2
        criteria_met.append("cross_pathway")

    protein = row.get("protein", "")
    for key, info in KNOWN_SWITCH_MECHANISMS.items():
        if protein == key[0] and (slim1_class == key[1] or slim2_class == key[1]):
            score += 3
            criteria_met.append("literature_validated")
            break

    exp1 = exposure_data.get((protein, slim1_class, str(row.get("slim1_pos", ""))), {})
    exp2 = exposure_data.get((protein, slim2_class, str(row.get("slim2_pos", ""))), {})
    exp1_score = exp1.get("exposure_score", 0.5)
    exp2_score = exp2.get("exposure_score", 0.5)
    if exp1_score > 0.6 and exp2_score > 0.6:
        score += 2
        criteria_met.append("both_exposed")
    elif exp1_score > 0.5 or exp2_score > 0.5:
        score += 1
        criteria_met.append("one_exposed")

    cons1 = conservation_data.get((protein, slim1_class), {})
    cons2 = conservation_data.get((protein, slim2_class), {})
    cons1_score = cons1.get("conservation_score", 0)
    cons2_score = cons2.get("conservation_score", 0)
    if cons1_score > 0.5 or cons2_score > 0.5:
        score += 1
        criteria_met.append("conserved_motif")

    max_score = 16
    confidence = score / max_score

    if confidence >= 0.5:
        tier = "high"
    elif confidence >= 0.3:
        tier = "medium"
    else:
        tier = "low"

    return {
        "confidence_score": round(confidence, 3),
        "confidence_tier": tier,
        "raw_score": score,
        "max_score": max_score,
        "criteria_met": criteria_met,
    }


def build_exposure_lookup(exposure_df):
    """Build lookup dict from exposure estimates."""
    lookup = {}
    for _, row in exposure_df.iterrows():
        key = (row["protein"], row["motif_class"], f"{int(row['start'])}-{int(row['end'])}")
        lookup[key] = {
            "exposure_score": row.get("exposure_score", 0.5),
            "exposure_class": row.get("exposure_class", "unknown"),
        }
    return lookup


def build_conservation_lookup(conservation_df):
    """Build lookup dict from conservation data."""
    lookup = {}
    if conservation_df is not None:
        for _, row in conservation_df.iterrows():
            key = (row.get("protein", ""), row.get("motif_class", ""))
            lookup[key] = {
                "conservation_score": row.get("mean_conservation", 0),
            }
    return lookup


def identify_switches(competition_df, exposure_lookup, conservation_lookup):
    """Full switch identification pipeline with decision tree."""
    switches = []

    for _, row in competition_df.iterrows():
        if not row.get("is_switch_candidate", False):
            continue

        conf = compute_switch_confidence(row, exposure_lookup, conservation_lookup)

        # Require at least 3 supporting criteria for inclusion
        if len(conf["criteria_met"]) < 3:
            continue

        category = classify_switch(str(row.get("slim1_pathway", "")), str(row.get("slim2_pathway", "")))

        protein = row["protein"]
        mechanism = ""
        validated = False
        for key, info in KNOWN_SWITCH_MECHANISMS.items():
            if protein == key[0] and (row["slim1_class"] == key[1] or row["slim2_class"] == key[1]):
                mechanism = info["mechanism"]
                validated = info["validated"]
                break

        switches.append({
            "protein": protein,
            "slim1_class": row["slim1_class"],
            "slim1_pos": row["slim1_pos"],
            "slim1_pathway": row["slim1_pathway"],
            "slim2_class": row["slim2_class"],
            "slim2_pos": row["slim2_pos"],
            "slim2_pathway": row["slim2_pathway"],
            "distance": row.get("distance", 0),
            "overlap": row.get("overlap", 0),
            "competition_score": row.get("competition_score", 0),
            "switch_category": category,
            "confidence_score": conf["confidence_score"],
            "confidence_tier": conf["confidence_tier"],
            "criteria_met": "; ".join(conf["criteria_met"]),
            "n_criteria": len(conf["criteria_met"]),
            "known_mechanism": mechanism,
            "literature_validated": validated,
        })

    switches.sort(key=lambda x: -x["confidence_score"])
    return switches


def hypothesis_test_h2(switches, all_hits_df):
    """
    H2: Target proteins have more cross-pathway switches than controls.
    """
    target_proteins = {"ICP34.5", "US3", "pUL37"}

    target_switches = [s for s in switches if s["protein"] in target_proteins]
    control_switches = [s for s in switches if s["protein"] not in target_proteins]

    per_protein = defaultdict(lambda: {"total": 0, "high": 0, "cross_pathway": 0, "categories": set()})
    for s in switches:
        p = s["protein"]
        per_protein[p]["total"] += 1
        if s["confidence_tier"] == "high":
            per_protein[p]["high"] += 1
        if "cross_pathway" in s["criteria_met"]:
            per_protein[p]["cross_pathway"] += 1
        per_protein[p]["categories"].add(s["switch_category"])

    target_counts = []
    control_counts = []
    for p, data in per_protein.items():
        pdata = all_hits_df[all_hits_df["protein"] == p]
        n_slims = len(pdata)
        density = data["total"] / max(n_slims, 1)
        entry = {
            "protein": p,
            "n_switches": data["total"],
            "n_high_conf": data["high"],
            "n_cross_pathway": data["cross_pathway"],
            "n_categories": len(data["categories"]),
            "switch_density": round(density, 3),
            "categories": sorted(data["categories"]),
        }
        if p in target_proteins:
            target_counts.append(entry)
        else:
            control_counts.append(entry)

    t_vals = [e["switch_density"] for e in target_counts]
    c_vals = [e["switch_density"] for e in control_counts]

    if len(t_vals) >= 2 and len(c_vals) >= 2:
        u_stat, u_pval = stats.mannwhitneyu(t_vals, c_vals, alternative="greater")
    else:
        u_stat, u_pval = float("nan"), float("nan")

    # Post-hoc power analysis
    pooled_std = np.std(t_vals + c_vals) if len(t_vals + c_vals) > 1 else 1.0
    observed_d = (np.mean(t_vals) - np.mean(c_vals)) / pooled_std if pooled_std > 0 else 0

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from utils.power_analysis import full_power_report
    power = full_power_report(len(t_vals), len(c_vals), observed_d)
    print(f"  H2 power analysis: achieved power={power['achieved_power']:.1%}, "
          f"MDE={power['min_detectable_effect_d_at_80pct']}")

    return {
        "targets": sorted(target_counts, key=lambda x: -x["switch_density"]),
        "controls": sorted(control_counts, key=lambda x: -x["switch_density"]),
        "mann_whitney_U": round(float(u_stat), 2) if not np.isnan(u_stat) else None,
        "mann_whitney_p": round(float(u_pval), 4) if not np.isnan(u_pval) else None,
        "target_mean_density": round(np.mean(t_vals), 3) if t_vals else 0,
        "control_mean_density": round(np.mean(c_vals), 3) if c_vals else 0,
        "power_analysis": power,
    }


def hypothesis_test_h3(switches, exposure_lookup):
    """
    H3: Target protein switches preferentially target CNS-enriched pathways.
    """
    target_proteins = {"ICP34.5", "US3", "pUL37"}
    cns_pathways = {"synaptic_adhesion", "calcium_signaling", "axonal_transport",
                    "neuronal_signaling", "synaptic_transmission"}

    target_cns = sum(1 for s in switches
                     if s["protein"] in target_proteins
                     and (s["slim1_pathway"] in cns_pathways or s["slim2_pathway"] in cns_pathways))
    target_total = sum(1 for s in switches if s["protein"] in target_proteins)

    control_cns = sum(1 for s in switches
                      if s["protein"] not in target_proteins
                      and (s["slim1_pathway"] in cns_pathways or s["slim2_pathway"] in cns_pathways))
    control_total = sum(1 for s in switches if s["protein"] not in target_proteins)

    table = [[target_cns, target_total - target_cns],
             [control_cns, control_total - control_cns]]
    try:
        odds, pval = stats.fisher_exact(table, alternative="greater")
    except ValueError:
        odds, pval = float("nan"), float("nan")

    # Post-hoc power analysis for Fisher's exact
    p1 = target_cns / max(target_total, 1)
    p2 = control_cns / max(control_total, 1)

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from utils.power_analysis import fisher_power
    achieved_power = fisher_power(max(p1, p2 + 0.1), p2, target_total, control_total)
    print(f"  H3 power analysis: achieved power={achieved_power:.1%} "
          f"(for detecting 10pp difference in CNS fraction)")

    return {
        "target_cns_switches": target_cns,
        "target_total_switches": target_total,
        "target_cns_fraction": round(target_cns / max(target_total, 1), 3),
        "control_cns_switches": control_cns,
        "control_total_switches": control_total,
        "control_cns_fraction": round(control_cns / max(control_total, 1), 3),
        "fisher_odds_ratio": round(float(odds), 2) if not np.isnan(odds) else None,
        "fisher_p_value": round(float(pval), 4) if not np.isnan(pval) else None,
        "power_analysis": {
            "achieved_power_10pp_diff": round(achieved_power, 3),
            "target_cns_proportion": round(p1, 3),
            "control_cns_proportion": round(p2, 3),
            "interpretation": (
                f"With {target_total} target vs {control_total} control switches, "
                f"this test has {achieved_power:.0%} power to detect a 10 percentage-point "
                f"difference in CNS pathway fraction. Non-significant results are inconclusive."
            ),
        },
    }


# ═══════════════════════════════════════════════════════════════
# PLOTTING
# ═══════════════════════════════════════════════════════════════

def plot_switch_summary(switches, output_path):
    """Overview of switch candidates."""
    fig = plt.figure(figsize=(20, 12))
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

    # A: Switches per protein
    ax = fig.add_subplot(gs[0, 0])
    per_protein = Counter(s["protein"] for s in switches)
    proteins = sorted(per_protein.keys(), key=lambda x: -per_protein[x])
    counts = [per_protein[p] for p in proteins]
    target_set = {"ICP34.5", "US3", "pUL37"}
    colors = ["#e74c3c" if p in target_set else "#3498db" for p in proteins]
    ax.barh(range(len(proteins)), counts, color=colors, alpha=0.7)
    ax.set_yticks(range(len(proteins)))
    ax.set_yticklabels(proteins, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Number of switch candidates")
    ax.set_title("A. Switch Candidates per Protein", fontweight="bold")

    # B: Confidence distribution
    ax = fig.add_subplot(gs[0, 1])
    tiers = Counter(s["confidence_tier"] for s in switches)
    tier_names = ["high", "medium", "low"]
    tier_colors = ["#27ae60", "#f39c12", "#95a5a6"]
    tier_counts = [tiers.get(t, 0) for t in tier_names]
    ax.bar(tier_names, tier_counts, color=tier_colors, alpha=0.7)
    ax.set_ylabel("Count")
    ax.set_title("B. Confidence Distribution", fontweight="bold")
    for i, c in enumerate(tier_counts):
        ax.text(i, c + 5, str(c), ha="center", fontsize=10)

    # C: Category distribution
    ax = fig.add_subplot(gs[0, 2])
    cats = Counter(s["switch_category"] for s in switches)
    cat_names = sorted(cats.keys(), key=lambda x: -cats[x])
    cat_counts = [cats[c] for c in cat_names]
    cat_colors = plt.cm.Set2(np.linspace(0, 1, len(cat_names)))
    ax.barh(range(len(cat_names)), cat_counts, color=cat_colors, alpha=0.7)
    ax.set_yticks(range(len(cat_names)))
    ax.set_yticklabels([c.replace("_", " ") for c in cat_names], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Count")
    ax.set_title("C. Switch Categories", fontweight="bold")

    # D: Competition score distribution
    ax = fig.add_subplot(gs[1, 0])
    target_scores = [s["competition_score"] for s in switches if s["protein"] in target_set]
    control_scores = [s["competition_score"] for s in switches if s["protein"] not in target_set]
    ax.hist(target_scores, bins=20, alpha=0.6, color="#e74c3c", label="Targets", density=True)
    ax.hist(control_scores, bins=20, alpha=0.6, color="#3498db", label="Controls", density=True)
    ax.set_xlabel("Competition score")
    ax.set_ylabel("Density")
    ax.set_title("D. Competition Score Distribution", fontweight="bold")
    ax.legend()

    # E: High-confidence switches per protein
    ax = fig.add_subplot(gs[1, 1])
    high_conf = [s for s in switches if s["confidence_tier"] == "high"]
    high_per_protein = Counter(s["protein"] for s in high_conf)
    hp_names = sorted(high_per_protein.keys(), key=lambda x: -high_per_protein[x])
    hp_counts = [high_per_protein[p] for p in hp_names]
    hp_colors = ["#e74c3c" if p in target_set else "#3498db" for p in hp_names]
    ax.barh(range(len(hp_names)), hp_counts, color=hp_colors, alpha=0.7)
    ax.set_yticks(range(len(hp_names)))
    ax.set_yticklabels(hp_names, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("High-confidence switches")
    ax.set_title("E. High-Confidence Switches", fontweight="bold")

    # F: Cross-pathway vs same-pathway
    ax = fig.add_subplot(gs[1, 2])
    for protein in target_set:
        p_switches = [s for s in switches if s["protein"] == protein]
        cross = sum(1 for s in p_switches if "cross_pathway" in s["criteria_met"])
        same = len(p_switches) - cross
        ax.bar(protein, cross, color="#e74c3c", alpha=0.7, label="Cross-pathway" if protein == "ICP34.5" else "")
        ax.bar(protein, same, bottom=cross, color="#f39c12", alpha=0.5, label="Same-pathway" if protein == "ICP34.5" else "")
    ax.set_ylabel("Switch count")
    ax.set_title("F. Target Cross-Pathway Switches", fontweight="bold")
    ax.legend()

    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


def plot_top_switches(switches, output_path):
    """Detailed view of top 20 switches."""
    top = switches[:20]
    fig, ax = plt.subplots(figsize=(16, 10))

    target_set = {"ICP34.5", "US3", "pUL37"}
    y = range(len(top))
    labels = []
    for s in top:
        label = f"{s['protein']}: {s['slim1_class']}({s['slim1_pathway']}) vs {s['slim2_class']}({s['slim2_pathway']})"
        labels.append(label[:80])

    colors = ["#e74c3c" if s["protein"] in target_set else "#3498db" for s in top]
    bars = ax.barh(y, [s["confidence_score"] for s in top], color=colors, alpha=0.7)

    for i, s in enumerate(top):
        marker = " *" if s["literature_validated"] else ""
        ax.text(s["confidence_score"] + 0.01, i, f"{s['confidence_tier']}{marker}", fontsize=8, va="center")

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Confidence score", fontsize=11)
    ax.set_title("Top 20 Competitive-Binding Switch Candidates (* = literature validated)",
                 fontsize=12, fontweight="bold")
    ax.axvline(0.5, color="green", linestyle="--", alpha=0.5, label="High confidence threshold")
    ax.axvline(0.3, color="orange", linestyle="--", alpha=0.5, label="Medium confidence threshold")
    ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

def main():
    print("=" * 80)
    print("PHASE 5 — Competitive-Binding Switch Identification")
    print("=" * 80)

    # Load data
    competition_df = pd.read_csv(os.path.join(PHASE3_DIR, "results", "competitive_accessibility.csv"))
    print(f"Loaded {len(competition_df)} competitive pairs")

    exposure_df = pd.read_csv(os.path.join(PHASE3_DIR, "results", "slim_exposure_estimates.csv"))
    exposure_lookup = build_exposure_lookup(exposure_df)
    print(f"Loaded {len(exposure_lookup)} exposure estimates")

    cons_path = os.path.join(PHASE2_DIR, "results", "slim_conservation.csv")
    if os.path.exists(cons_path):
        cons_df = pd.read_csv(cons_path)
        conservation_lookup = build_conservation_lookup(cons_df)
        print(f"Loaded {len(conservation_lookup)} conservation scores")
    else:
        conservation_lookup = {}
        print("No conservation data found, proceeding without")

    hits_df = pd.read_csv(os.path.join(PHASE2_DIR, "results", "slim_hits_high_confidence.csv"), encoding="utf-8")
    print(f"Loaded {len(hits_df)} high-confidence SLiM hits")

    # 1. Switch identification
    print("\n[1] Identifying switches with decision tree...")
    switches = identify_switches(competition_df, exposure_lookup, conservation_lookup)
    print(f"  Total switch candidates: {len(switches)}")

    high = [s for s in switches if s["confidence_tier"] == "high"]
    medium = [s for s in switches if s["confidence_tier"] == "medium"]
    low = [s for s in switches if s["confidence_tier"] == "low"]
    print(f"  High confidence: {len(high)}")
    print(f"  Medium confidence: {len(medium)}")
    print(f"  Low confidence: {len(low)}")

    # 2. Category breakdown
    print("\n[2] Switch categories:")
    cat_counts = Counter(s["switch_category"] for s in switches)
    for cat, count in cat_counts.most_common():
        desc = SWITCH_CATEGORIES.get(cat, {}).get("description", "")
        print(f"  {cat:<25} {count:>5}  {desc}")

    # 3. Per-protein breakdown
    print("\n[3] Per-protein switch summary:")
    target_set = {"ICP34.5", "US3", "pUL37"}
    print(f"  {'Protein':<15} {'Cat':<8} {'Total':>6} {'High':>6} {'Med':>6} {'Low':>6} {'Categories'}")
    print("  " + "-" * 80)
    per_protein = defaultdict(lambda: {"total": 0, "high": 0, "medium": 0, "low": 0, "cats": set()})
    for s in switches:
        p = s["protein"]
        per_protein[p]["total"] += 1
        per_protein[p][s["confidence_tier"]] += 1
        per_protein[p]["cats"].add(s["switch_category"])

    for p in sorted(per_protein.keys(), key=lambda x: -per_protein[x]["total"]):
        d = per_protein[p]
        cat_label = "target" if p in target_set else "control"
        cats = ", ".join(sorted(d["cats"]))
        print(f"  {p:<15} {cat_label:<8} {d['total']:>6} {d['high']:>6} {d['medium']:>6} "
              f"{d['low']:>6} {cats}")

    # 4. Top switches
    print("\n[4] Top 15 switches:")
    print(f"  {'Protein':<12} {'SLiM1':<20} {'SLiM2':<20} {'Category':<18} {'Conf':>5} {'Tier':<6} {'Valid'}")
    print("  " + "-" * 95)
    for s in switches[:15]:
        val = "YES" if s["literature_validated"] else ""
        print(f"  {s['protein']:<12} {s['slim1_class']:<20} {s['slim2_class']:<20} "
              f"{s['switch_category']:<18} {s['confidence_score']:>5.3f} {s['confidence_tier']:<6} {val}")

    # 5. Literature-validated switches
    print("\n[5] Literature-validated switches:")
    validated = [s for s in switches if s["literature_validated"]]
    for s in validated:
        print(f"  {s['protein']}: {s['slim1_class']} vs {s['slim2_class']}")
        print(f"    Mechanism: {s['known_mechanism']}")
        print(f"    Confidence: {s['confidence_score']:.3f} ({s['confidence_tier']})")

    # 6. Hypothesis tests
    print("\n[6] Hypothesis test H2 (switch density)...")
    h2 = hypothesis_test_h2(switches, hits_df)
    print(f"  Target mean switch density: {h2['target_mean_density']}")
    print(f"  Control mean switch density: {h2['control_mean_density']}")
    print(f"  Mann-Whitney p: {h2['mann_whitney_p']}")

    print("\n[7] Hypothesis test H3 (CNS pathway switches)...")
    h3 = hypothesis_test_h3(switches, exposure_lookup)
    print(f"  Target CNS fraction: {h3['target_cns_fraction']}")
    print(f"  Control CNS fraction: {h3['control_cns_fraction']}")
    print(f"  Fisher's exact p: {h3['fisher_p_value']}")

    # 7. Plots
    print("\n[8] Generating plots...")
    plot_switch_summary(switches, os.path.join(PLOT_DIR, "switch_summary.png"))
    print("  Switch summary")
    plot_top_switches(switches, os.path.join(PLOT_DIR, "top_switches.png"))
    print("  Top switches")

    # 8. Save results
    print("\n[9] Saving results...")
    switches_df = pd.DataFrame(switches)
    switches_df.to_csv(os.path.join(RESULTS_DIR, "switch_candidates.csv"), index=False, encoding="utf-8")
    print(f"  {len(switches)} switch candidates saved")

    high_df = switches_df[switches_df["confidence_tier"].isin(["high", "medium"])]
    high_df.to_csv(os.path.join(RESULTS_DIR, "switch_candidates_confident.csv"), index=False, encoding="utf-8")
    print(f"  {len(high_df)} high/medium confidence switches")

    summary = {
        "total_switches": len(switches),
        "high_confidence": len(high),
        "medium_confidence": len(medium),
        "low_confidence": len(low),
        "categories": {cat: count for cat, count in cat_counts.most_common()},
        "per_protein": {p: {"total": d["total"], "high": d["high"], "categories": sorted(d["cats"])}
                        for p, d in per_protein.items()},
        "hypothesis_h2": h2,
        "hypothesis_h3": h3,
        "n_validated": len(validated),
    }
    with open(os.path.join(RESULTS_DIR, "phase5_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)

    with open(os.path.join(RESULTS_DIR, "hypothesis_h2_switches.json"), "w") as f:
        json.dump(h2, f, indent=2, default=str)

    with open(os.path.join(RESULTS_DIR, "hypothesis_h3_cns_switches.json"), "w") as f:
        json.dump(h3, f, indent=2, default=str)

    print(f"\nPhase 5 complete! Results in {RESULTS_DIR}")


if __name__ == "__main__":
    main()
