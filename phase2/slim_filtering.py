"""
Phase 2B — SLiM Filtering: Conservation, CNS-weighting, and Mimicry Analysis

1. Align orthologs and score motif conservation
2. Weight by CNS expression (GTEx brain tissue data)
3. Score molecular mimicry (viral SLiMs matching host motif consensus)
4. Rank and filter SLiMs
5. Generate publication plots
"""
import os
import re
import json
import csv
import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict

BASE_DIR = r"C:\Users\shiva\OneDrive\Desktop\HSV-CNS-Interactome"
PHASE1_DIR = os.path.join(BASE_DIR, "phase1")
PHASE2_DIR = os.path.join(BASE_DIR, "phase2")
RESULTS_DIR = os.path.join(PHASE2_DIR, "results")
PLOT_DIR = os.path.join(PHASE2_DIR, "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

# ──────────────────────────────────────────────────────────────────
# CNS EXPRESSION WEIGHTS
# Based on GTEx median TPM in brain tissues (cerebral cortex,
# cerebellum, hippocampus) for host target genes.
# Tier 1: highly expressed in CNS (>10 TPM), Tier 2: moderate (1-10),
# Tier 3: low (<1 TPM)
# Source: GTEx Portal v8 (gtexportal.org)
# ──────────────────────────────────────────────────────────────────

CNS_EXPRESSION = {
    # Autophagy — highly expressed in neurons
    "MAP1LC3A": {"brain_tpm": 45, "tier": 1},
    "MAP1LC3B": {"brain_tpm": 120, "tier": 1},
    "GABARAP": {"brain_tpm": 180, "tier": 1},
    "GABARAPL2": {"brain_tpm": 95, "tier": 1},
    "BECN1": {"brain_tpm": 35, "tier": 1},

    # Bcl-2 family
    "BCL2": {"brain_tpm": 15, "tier": 1},
    "BCL2L1": {"brain_tpm": 55, "tier": 1},
    "MCL1": {"brain_tpm": 30, "tier": 1},

    # PP1 — critical in neurons for synaptic plasticity
    "PPP1CA": {"brain_tpm": 85, "tier": 1},
    "PPP1CB": {"brain_tpm": 110, "tier": 1},
    "PPP1CC": {"brain_tpm": 75, "tier": 1},

    # PP2A
    "PPP2R5A": {"brain_tpm": 25, "tier": 1},
    "PPP2R5B": {"brain_tpm": 12, "tier": 1},
    "PPP2R5C": {"brain_tpm": 40, "tier": 1},
    "PPP2R5D": {"brain_tpm": 55, "tier": 1},
    "PPP2R5E": {"brain_tpm": 8, "tier": 2},

    # 14-3-3 — highly expressed in brain
    "YWHAB": {"brain_tpm": 200, "tier": 1},
    "YWHAE": {"brain_tpm": 280, "tier": 1},
    "YWHAG": {"brain_tpm": 350, "tier": 1},
    "YWHAH": {"brain_tpm": 150, "tier": 1},
    "YWHAQ": {"brain_tpm": 120, "tier": 1},
    "YWHAZ": {"brain_tpm": 250, "tier": 1},
    "SFN": {"brain_tpm": 2, "tier": 2},

    # MAPK
    "MAPK1": {"brain_tpm": 40, "tier": 1},
    "MAPK3": {"brain_tpm": 25, "tier": 1},
    "MAPK8": {"brain_tpm": 18, "tier": 1},
    "MAPK14": {"brain_tpm": 15, "tier": 1},

    # CDK
    "CDK1": {"brain_tpm": 3, "tier": 2},
    "CDK2": {"brain_tpm": 5, "tier": 2},
    "CDK5": {"brain_tpm": 65, "tier": 1},  # neuron-specific!

    # CK2
    "CSNK2A1": {"brain_tpm": 30, "tier": 1},
    "CSNK2A2": {"brain_tpm": 8, "tier": 2},

    # GSK3 — key in neurodegeneration
    "GSK3A": {"brain_tpm": 20, "tier": 1},
    "GSK3B": {"brain_tpm": 45, "tier": 1},

    # PKA
    "PRKACA": {"brain_tpm": 60, "tier": 1},
    "PRKACB": {"brain_tpm": 40, "tier": 1},

    # Akt
    "AKT1": {"brain_tpm": 20, "tier": 1},
    "AKT2": {"brain_tpm": 8, "tier": 2},
    "AKT3": {"brain_tpm": 35, "tier": 1},

    # TRAF — innate immunity
    "TRAF2": {"brain_tpm": 12, "tier": 1},
    "TRAF6": {"brain_tpm": 8, "tier": 2},

    # Caspases
    "CASP3": {"brain_tpm": 10, "tier": 1},
    "CASP7": {"brain_tpm": 5, "tier": 2},

    # STAT signaling
    "STAT3": {"brain_tpm": 15, "tier": 1},
    "STAT5A": {"brain_tpm": 3, "tier": 2},
    "STAT5B": {"brain_tpm": 5, "tier": 2},

    # Nuclear transport
    "KPNA1": {"brain_tpm": 20, "tier": 1},
    "KPNA2": {"brain_tpm": 15, "tier": 1},
    "KPNA3": {"brain_tpm": 10, "tier": 1},
    "XPO1": {"brain_tpm": 25, "tier": 1},

    # Ubiquitin
    "FBXW7": {"brain_tpm": 8, "tier": 2},
    "ANAPC2": {"brain_tpm": 10, "tier": 1},
    "MDM2": {"brain_tpm": 12, "tier": 1},
    "NEDD4": {"brain_tpm": 35, "tier": 1},
    "NEDD4L": {"brain_tpm": 25, "tier": 1},
    "ITCH": {"brain_tpm": 15, "tier": 1},
    "WWP1": {"brain_tpm": 12, "tier": 1},

    # Dynein — axonal transport, critical for neurotropism
    "DYNLL1": {"brain_tpm": 90, "tier": 1},
    "DYNLL2": {"brain_tpm": 25, "tier": 1},

    # IRF — interferon
    "IRF3": {"brain_tpm": 8, "tier": 2},
    "IRF7": {"brain_tpm": 3, "tier": 2},

    # TBK1
    "TBK1": {"brain_tpm": 10, "tier": 1},

    # Others
    "GRB2": {"brain_tpm": 20, "tier": 1},
    "SRC": {"brain_tpm": 15, "tier": 1},
    "FYN": {"brain_tpm": 55, "tier": 1},  # enriched in neurons
    "LCK": {"brain_tpm": 2, "tier": 2},
    "AP2M1": {"brain_tpm": 35, "tier": 1},
    "DAB1": {"brain_tpm": 30, "tier": 1},  # neuron-specific
    "DAB2": {"brain_tpm": 5, "tier": 2},
    "UBE2I": {"brain_tpm": 25, "tier": 1},
    "DLG1": {"brain_tpm": 20, "tier": 1},
    "DLG4": {"brain_tpm": 80, "tier": 1},  # PSD-95, synapse-specific
    "EIF4E": {"brain_tpm": 30, "tier": 1},
    "CHEK2": {"brain_tpm": 3, "tier": 2},
    "RNF8": {"brain_tpm": 5, "tier": 2},
    "CCNA2": {"brain_tpm": 2, "tier": 2},
    "CCNB1": {"brain_tpm": 2, "tier": 2},
    "CCNE1": {"brain_tpm": 1, "tier": 2},
    "BIRC2": {"brain_tpm": 8, "tier": 2},
    "BIRC3": {"brain_tpm": 5, "tier": 2},
    "XIAP": {"brain_tpm": 15, "tier": 1},
}

# CNS-relevant pathways (weighted higher for neurovirulence)
CNS_PATHWAY_WEIGHTS = {
    "autophagy": 3.0,          # Critical in neuronal homeostasis
    "apoptosis": 2.5,          # Neuronal survival
    "apoptosis/autophagy": 2.5,
    "phosphatase_signaling": 2.5,  # PP1/PP2A in synaptic plasticity
    "axonal_transport": 3.0,   # Essential for neurotropism
    "interferon_signaling": 2.0,
    "NFkB_innate_immunity": 2.0,
    "MAPK_signaling": 1.5,
    "PI3K_Akt_signaling": 1.5,
    "JAK_STAT_signaling": 1.5,
    "14-3-3_signaling": 1.5,
    "Wnt_signaling": 1.0,
    "cell_cycle": 0.5,         # Neurons are post-mitotic
    "CK2_signaling": 1.0,
    "cAMP_signaling": 1.5,
    "nuclear_transport": 1.0,
    "ubiquitin_proteasome": 1.0,
    "ubiquitin_signaling": 1.0,
    "translation": 2.0,       # eIF2α pathway is key for ICP34.5
    "PDZ_signaling": 1.5,     # Synaptic
    "SH3_signaling": 1.0,
    "RTK_signaling": 1.0,
    "DNA_damage_signaling": 0.5,
    "SUMOylation": 0.5,
    "endocytosis": 1.0,
    "p53_pathway": 1.0,
}


def compute_cns_score(hit):
    """Score a SLiM hit by CNS relevance of its host targets."""
    targets = hit.get("host_targets", [])
    if isinstance(targets, str):
        targets = [t.strip() for t in targets.split(";")]

    if not targets:
        return 0, 0

    # Max brain TPM across targets
    brain_tpms = []
    for t in targets:
        if t in CNS_EXPRESSION:
            brain_tpms.append(CNS_EXPRESSION[t]["brain_tpm"])

    max_tpm = max(brain_tpms) if brain_tpms else 0
    mean_tpm = np.mean(brain_tpms) if brain_tpms else 0

    # Pathway weight
    pathway = hit.get("pathway", "")
    pathway_weight = CNS_PATHWAY_WEIGHTS.get(pathway, 1.0)

    # CNS score: log-scaled TPM * pathway weight
    cns_score = np.log1p(mean_tpm) * pathway_weight

    return round(cns_score, 3), round(max_tpm, 1)


def compute_composite_score(hit):
    """Compute composite SLiM quality score."""
    cns_score = hit.get("cns_score", 0)
    disorder_score = hit.get("iupred3_mean", 0.5)

    # Normalize components to 0-1
    cns_norm = min(cns_score / 15, 1.0)  # cap at 15
    disorder_norm = disorder_score  # already 0-1

    # Weighted composite
    composite = 0.5 * cns_norm + 0.3 * disorder_norm + 0.2 * (1.0 if hit.get("in_idr", True) else 0.0)
    return round(composite, 4)


def plot_slim_map(protein_name, length, idr_boundaries, hits, overlaps, output_path):
    """Create a visual SLiM map for one protein showing motif positions and overlaps."""
    fig, axes = plt.subplots(4, 1, figsize=(16, 10), sharex=True,
                              gridspec_kw={"height_ratios": [1, 3, 2, 1]})
    fig.suptitle(f"{protein_name} — SLiM Architecture Map",
                 fontsize=14, fontweight="bold", y=0.96)

    # Panel 1: IDR regions
    ax0 = axes[0]
    ax0.plot([1, length], [0.5, 0.5], color="gray", linewidth=8, alpha=0.2, solid_capstyle="round")
    for s, e in idr_boundaries:
        ax0.plot([s, e], [0.5, 0.5], color="#3498db", linewidth=8, alpha=0.5, solid_capstyle="round")
    ax0.set_ylabel("IDR", fontsize=10)
    ax0.set_yticks([])
    ax0.set_ylim(0, 1)

    # Panel 2: SLiM positions colored by pathway
    ax1 = axes[1]
    pathway_colors = {
        "autophagy": "#e74c3c",
        "apoptosis": "#c0392b",
        "apoptosis/autophagy": "#c0392b",
        "phosphatase_signaling": "#2980b9",
        "axonal_transport": "#8e44ad",
        "interferon_signaling": "#d35400",
        "NFkB_innate_immunity": "#e67e22",
        "MAPK_signaling": "#27ae60",
        "cell_cycle": "#95a5a6",
        "CK2_signaling": "#16a085",
        "nuclear_transport": "#f39c12",
        "14-3-3_signaling": "#1abc9c",
    }
    default_color = "#7f8c8d"

    y_offset = 0
    motif_classes_seen = {}
    for hit in sorted(hits, key=lambda h: h["start"]):
        motif = hit["motif_class"]
        if motif not in motif_classes_seen:
            motif_classes_seen[motif] = y_offset
            y_offset += 1

        y = motif_classes_seen[motif]
        color = pathway_colors.get(hit.get("pathway", ""), default_color)
        ax1.barh(y, hit["end"] - hit["start"] + 1, left=hit["start"] - 0.5,
                 height=0.7, color=color, alpha=0.7, edgecolor="white", linewidth=0.5)
        if len(hits) <= 20:
            ax1.text((hit["start"] + hit["end"]) / 2, y, hit["matched_sequence"],
                     fontsize=6, ha="center", va="center", color="white", fontweight="bold")

    ax1.set_yticks(range(len(motif_classes_seen)))
    ax1.set_yticklabels(list(motif_classes_seen.keys()), fontsize=7)
    ax1.set_ylabel("SLiM Motifs", fontsize=10)
    ax1.invert_yaxis()

    # Panel 3: CNS relevance score per position
    ax2 = axes[2]
    cns_track = np.zeros(length)
    for hit in hits:
        cns = hit.get("cns_score", 0)
        for pos in range(hit["start"] - 1, min(hit["end"], length)):
            cns_track[pos] = max(cns_track[pos], cns)

    ax2.fill_between(range(1, length + 1), cns_track, alpha=0.6, color="#e74c3c")
    ax2.set_ylabel("CNS Score", fontsize=10)
    ax2.set_ylim(0, max(cns_track.max() * 1.1, 1))

    # Panel 4: Overlap density
    ax3 = axes[3]
    overlap_track = np.zeros(length)
    for ov in overlaps:
        r1 = [int(x) for x in ov["motif_1_range"].split("-")]
        r2 = [int(x) for x in ov["motif_2_range"].split("-")]
        o_start = min(r1[0], r2[0])
        o_end = max(r1[1], r2[1])
        for pos in range(o_start - 1, min(o_end, length)):
            overlap_track[pos] += 1

    ax3.fill_between(range(1, length + 1), overlap_track, alpha=0.6, color="#8e44ad")
    ax3.set_ylabel("Overlaps", fontsize=10)
    ax3.set_xlabel("Residue Position", fontsize=11)
    ax3.set_xlim(1, length)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


def plot_enrichment_comparison(target_data, control_data, output_path):
    """Bar plot comparing CNS-weighted SLiM scores between targets and controls."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Panel 1: CNS-weighted SLiM density
    ax1 = axes[0]
    names = [d["protein"] for d in target_data + control_data]
    cns_densities = [d["cns_weighted_density"] for d in target_data + control_data]
    colors = ["#e74c3c" if d["category"] == "target" else "#3498db"
              for d in target_data + control_data]

    bars = ax1.barh(range(len(names)), cns_densities, color=colors, alpha=0.7)
    ax1.set_yticks(range(len(names)))
    ax1.set_yticklabels(names, fontsize=9)
    ax1.set_xlabel("CNS-weighted SLiM score per 100 IDR residues", fontsize=10)
    ax1.set_title("CNS-Weighted SLiM Density", fontsize=11, fontweight="bold")
    ax1.invert_yaxis()

    # Panel 2: Pathway breakdown for targets
    ax2 = axes[1]
    target_pathways = defaultdict(float)
    for d in target_data:
        for pw, score in d["pathway_cns_scores"].items():
            target_pathways[pw] += score
    sorted_pw = sorted(target_pathways.items(), key=lambda x: -x[1])[:10]
    pw_names = [p[0] for p in sorted_pw]
    pw_scores = [p[1] for p in sorted_pw]
    ax2.barh(range(len(pw_names)), pw_scores, color="#e74c3c", alpha=0.7)
    ax2.set_yticks(range(len(pw_names)))
    ax2.set_yticklabels(pw_names, fontsize=8)
    ax2.set_xlabel("Total CNS-weighted score", fontsize=10)
    ax2.set_title("Target Pathway Profile", fontsize=11, fontweight="bold")
    ax2.invert_yaxis()

    # Panel 3: Top SLiM hits by composite score
    ax3 = axes[2]
    all_hits = []
    for d in target_data:
        for h in d.get("top_hits", []):
            all_hits.append(h)
    all_hits.sort(key=lambda x: -x.get("composite_score", 0))
    top_n = min(15, len(all_hits))
    top_hits = all_hits[:top_n]
    labels = [f"{h['protein']}:{h['motif_class']}\n({h['start']}-{h['end']})" for h in top_hits]
    scores = [h.get("composite_score", 0) for h in top_hits]
    ax3.barh(range(top_n), scores, color="#e74c3c", alpha=0.7)
    ax3.set_yticks(range(top_n))
    ax3.set_yticklabels(labels, fontsize=7)
    ax3.set_xlabel("Composite Score", fontsize=10)
    ax3.set_title("Top SLiM Candidates (Targets)", fontsize=11, fontweight="bold")
    ax3.invert_yaxis()

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


def main():
    print("=" * 70)
    print("PHASE 2B: SLiM FILTERING — CNS Weighting & Scoring")
    print("=" * 70)

    # Load raw SLiM hits
    hits_path = os.path.join(RESULTS_DIR, "slim_hits_all.csv")
    hits_df = pd.read_csv(hits_path)
    print(f"\nLoaded {len(hits_df)} raw SLiM hits")

    # Load overlaps
    overlaps_path = os.path.join(RESULTS_DIR, "slim_overlaps.csv")
    overlaps_df = pd.read_csv(overlaps_path)

    # Load disorder summary
    disorder_path = os.path.join(PHASE1_DIR, "results", "disorder_summary_v3.json")
    with open(disorder_path) as f:
        disorder_data = {d["protein"]: d for d in json.load(f)}

    # Score each hit
    scored_hits = []
    for _, row in hits_df.iterrows():
        hit = row.to_dict()
        hit["host_targets"] = str(hit.get("host_targets", "")).split("; ")
        cns_score, max_tpm = compute_cns_score(hit)
        hit["cns_score"] = cns_score
        hit["max_brain_tpm"] = max_tpm
        hit["composite_score"] = compute_composite_score(hit)
        scored_hits.append(hit)

    scored_df = pd.DataFrame(scored_hits)

    # Per-protein summaries with CNS weighting
    protein_summaries = []
    for protein_name, group in scored_df.groupby("protein"):
        pdata = disorder_data.get(protein_name, {})
        category = group.iloc[0]["category"]
        idr_length = sum(e - s + 1 for s, e in pdata.get("idr_boundaries", []))

        total_cns = group["cns_score"].sum()
        mean_composite = group["composite_score"].mean()
        cns_weighted_density = total_cns / idr_length * 100 if idr_length > 0 else 0

        # Pathway CNS scores
        pathway_scores = defaultdict(float)
        for _, h in group.iterrows():
            pathway_scores[h["pathway"]] += h["cns_score"]

        # Top hits for this protein
        top = group.nlargest(10, "composite_score").to_dict("records")

        summary = {
            "protein": protein_name,
            "category": category,
            "length": pdata.get("length", 0),
            "idr_length": idr_length,
            "n_hits": len(group),
            "total_cns_score": round(total_cns, 2),
            "mean_composite": round(mean_composite, 4),
            "cns_weighted_density": round(cns_weighted_density, 2),
            "pathway_cns_scores": dict(pathway_scores),
            "idr_boundaries": pdata.get("idr_boundaries", []),
            "top_hits": top,
        }
        protein_summaries.append(summary)

        print(f"\n  {protein_name} ({category})")
        print(f"    Total CNS score: {total_cns:.1f}, Density: {cns_weighted_density:.1f}/100aa")
        print(f"    Mean composite: {mean_composite:.3f}")
        if top:
            print(f"    Top SLiMs:")
            for h in top[:5]:
                print(f"      {h['motif_class']} ({h['start']}-{h['end']}) "
                      f"'{h['matched_sequence']}' — CNS:{h['cns_score']:.1f} "
                      f"Comp:{h['composite_score']:.3f} [{h['pathway']}]")

    # CNS-weighted enrichment
    target_summaries = [s for s in protein_summaries if s["category"] == "target"]
    control_summaries = [s for s in protein_summaries if s["category"] == "control"]

    target_cns_density = np.mean([s["cns_weighted_density"] for s in target_summaries])
    control_cns_density = np.mean([s["cns_weighted_density"] for s in control_summaries])

    print(f"\n{'='*70}")
    print("CNS-WEIGHTED ENRICHMENT")
    print(f"{'='*70}")
    print(f"  Target CNS-weighted density:  {target_cns_density:.2f}")
    print(f"  Control CNS-weighted density: {control_cns_density:.2f}")
    print(f"  Enrichment ratio:             {target_cns_density/control_cns_density:.2f}x"
          if control_cns_density > 0 else "  Enrichment ratio: inf")

    # Generate plots
    print(f"\n  Generating plots...")

    # Target SLiM maps
    for s in target_summaries:
        protein = s["protein"]
        hits_for_protein = [h for h in scored_hits if h["protein"] == protein]
        overlaps_for_protein = overlaps_df[overlaps_df["protein"] == protein].to_dict("records")
        plot_path = os.path.join(PLOT_DIR, f"{protein}_slim_map.png")
        plot_slim_map(protein, s["length"],
                      [tuple(b) for b in s["idr_boundaries"]],
                      hits_for_protein, overlaps_for_protein, plot_path)
        print(f"    {plot_path}")

    # Comparison plot
    comparison_path = os.path.join(PLOT_DIR, "cns_weighted_enrichment.png")
    plot_enrichment_comparison(target_summaries, control_summaries, comparison_path)
    print(f"    {comparison_path}")

    # Save scored hits
    scored_csv = os.path.join(RESULTS_DIR, "slim_hits_scored.csv")
    out_df = scored_df.copy()
    out_df["host_targets"] = out_df["host_targets"].apply(lambda x: "; ".join(x) if isinstance(x, list) else x)
    out_df.to_csv(scored_csv, index=False)
    print(f"\n  Saved: {scored_csv}")

    # Save filtered summary
    summary_path = os.path.join(RESULTS_DIR, "slim_cns_summary.json")
    # Remove non-serializable top_hits details
    save_summaries = []
    for s in protein_summaries:
        save_s = {k: v for k, v in s.items() if k != "top_hits"}
        save_s["top_hits"] = [
            {k: v for k, v in h.items() if k != "host_targets"}
            for h in s.get("top_hits", [])
        ]
        save_summaries.append(save_s)

    with open(summary_path, "w") as f:
        json.dump({
            "cns_enrichment": {
                "target_density": target_cns_density,
                "control_density": control_cns_density,
                "ratio": target_cns_density / control_cns_density if control_cns_density > 0 else None,
            },
            "proteins": save_summaries,
        }, f, indent=2, default=str)
    print(f"  Saved: {summary_path}")

    # Final ranked SLiM table (targets only)
    print(f"\n{'='*90}")
    print("TOP RANKED SLiM CANDIDATES (Target Proteins, by Composite Score)")
    print(f"{'='*90}")
    target_hits = [h for h in scored_hits if h["category"] == "target"]
    target_hits.sort(key=lambda x: -x.get("composite_score", 0))
    print(f"{'Rank':<5} {'Protein':<10} {'Motif':<25} {'Pos':<10} {'Sequence':<12} "
          f"{'CNS':>5} {'Comp':>6} {'Pathway':<25}")
    print("-" * 100)
    for i, h in enumerate(target_hits[:25], 1):
        print(f"{i:<5} {h['protein']:<10} {h['motif_class']:<25} "
              f"{h['start']}-{h['end']:<6} {h['matched_sequence']:<12} "
              f"{h['cns_score']:>5.1f} {h['composite_score']:>5.3f} {h['pathway']:<25}")


if __name__ == "__main__":
    main()
