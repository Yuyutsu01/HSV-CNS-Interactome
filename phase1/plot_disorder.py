"""
Phase 1 — Generate disorder profile plots for all target proteins.
Creates publication-quality multi-panel figures.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import json

BASE_DIR = r"C:\Users\shiva\OneDrive\Desktop\HSV-CNS-Interactome\phase1"
PRED_DIR = os.path.join(BASE_DIR, "disorder_predictions")
PLOT_DIR = os.path.join(BASE_DIR, "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

# Known functional regions for annotation
FUNCTIONAL_REGIONS = {
    "ICP34.5": [
        (68, 87, "Beclin-1\nbinding", "#e74c3c"),
        (233, 248, "PP1\u03b1\nbinding", "#3498db"),
        (160, 190, "ATP\nrepeats", "#95a5a6"),
    ],
    "US3": [
        (12, 63, "Disordered\nN-tail", "#e67e22"),
        (52, 52, "T52P\nvirulence", "#e74c3c"),
        (153, 153, "A153T\nvirulence", "#e74c3c"),
        (100, 400, "Kinase\ndomain", "#2ecc71"),
    ],
    "pUL37": [
        (1, 27, "N-terminal\nIDR", "#e67e22"),
        (1047, 1123, "C-terminal\ntail/effector", "#e74c3c"),
    ],
}

# UniProt disorder annotations for comparison
UNIPROT_DISORDER = {
    "ICP34.5": [(1, 129), (149, 174)],
    "US3": [(12, 63)],
    "pUL37": [(1, 27), (1047, 1123)],
}


def plot_target_disorder(name, csv_path, output_path):
    """Create a detailed disorder profile plot for one target protein."""
    df = pd.read_csv(csv_path)
    n = len(df)
    positions = df["position"].values

    fig, axes = plt.subplots(3, 1, figsize=(14, 8), sharex=True,
                              gridspec_kw={"height_ratios": [3, 1, 1]})
    fig.suptitle(f"{name} — Disorder Prediction Profile", fontsize=14, fontweight="bold", y=0.95)

    # Panel 1: Prediction scores
    ax1 = axes[0]
    ax1.plot(positions, df["composition"], color="#2980b9", linewidth=1.2,
             label="Composition (TOP-IDP)", alpha=0.9)
    ax1.plot(positions, df["iupred_like"], color="#e74c3c", linewidth=1.0,
             label="Energy estimation", alpha=0.7)
    if "plddt" in df.columns and df["plddt"].notna().any():
        plddt_norm = 1.0 - df["plddt"].values / 100.0  # Invert: low pLDDT = high disorder
        ax1.plot(positions, plddt_norm, color="#27ae60", linewidth=1.0,
                 label="ESMFold (1 - pLDDT/100)", alpha=0.7)

    ax1.axhline(y=0.5, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    ax1.set_ylabel("Disorder Score", fontsize=11)
    ax1.set_ylim(-0.05, 1.05)
    ax1.legend(loc="upper right", fontsize=9)
    ax1.text(0.02, 0.95, "Disordered", transform=ax1.transAxes, fontsize=9,
             color="gray", va="top")
    ax1.text(0.02, 0.05, "Ordered", transform=ax1.transAxes, fontsize=9,
             color="gray", va="bottom")

    # Panel 2: Consensus IDR map
    ax2 = axes[1]
    consensus = df["consensus_disordered"].values

    # Draw consensus as colored bar
    for i in range(n):
        if consensus[i] == 1:
            ax2.axvspan(positions[i] - 0.5, positions[i] + 0.5,
                       color="#3498db", alpha=0.6)

    # Overlay UniProt annotations
    if name in UNIPROT_DISORDER:
        for start, end in UNIPROT_DISORDER[name]:
            ax2.axvspan(start - 0.5, end + 0.5, ymin=0.0, ymax=0.4,
                       color="#e74c3c", alpha=0.4)

    ax2.set_ylabel("IDR Map", fontsize=11)
    ax2.set_yticks([])
    ax2.set_ylim(0, 1)

    # Legend for IDR map
    predicted_patch = mpatches.Patch(color="#3498db", alpha=0.6, label="Predicted IDR (consensus)")
    uniprot_patch = mpatches.Patch(color="#e74c3c", alpha=0.4, label="UniProt annotated disorder")
    ax2.legend(handles=[predicted_patch, uniprot_patch], loc="upper right", fontsize=8)

    # Panel 3: Functional regions
    ax3 = axes[2]
    if name in FUNCTIONAL_REGIONS:
        for start, end, label, color in FUNCTIONAL_REGIONS[name]:
            if start == end:
                # Point mutation
                ax3.axvline(x=start, color=color, linewidth=2, alpha=0.8)
                ax3.text(start, 0.7, label, fontsize=7, ha="center", va="bottom",
                        color=color, fontweight="bold")
            else:
                ax3.axvspan(start - 0.5, end + 0.5, color=color, alpha=0.3)
                mid = (start + end) / 2
                ax3.text(mid, 0.5, label, fontsize=7, ha="center", va="center",
                        fontweight="bold")

    ax3.set_ylabel("Functional\nRegions", fontsize=10)
    ax3.set_yticks([])
    ax3.set_ylim(0, 1)
    ax3.set_xlabel("Residue Position", fontsize=11)
    ax3.set_xlim(0, n + 1)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


def plot_comparison_summary(results_path, output_path):
    """Create a summary comparison plot of IDR fractions across all proteins."""
    with open(results_path) as f:
        data = json.load(f)

    targets = [d for d in data if d["category"] == "target"]
    controls = [d for d in data if d["category"] == "control"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6),
                                     gridspec_kw={"width_ratios": [1, 2]})

    # Left panel: IDR fraction bar chart
    all_names = [d["protein"] for d in targets + controls]
    all_fracs = [d["idr_fraction"] for d in targets + controls]
    colors = ["#e74c3c" if d["category"] == "target" else "#3498db"
              for d in targets + controls]

    bars = ax1.barh(range(len(all_names)), all_fracs, color=colors, alpha=0.7, edgecolor="white")
    ax1.set_yticks(range(len(all_names)))
    ax1.set_yticklabels(all_names, fontsize=9)
    ax1.set_xlabel("IDR Fraction", fontsize=11)
    ax1.set_title("Predicted IDR Content", fontsize=12, fontweight="bold")
    ax1.invert_yaxis()
    ax1.axvline(x=0.5, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)

    target_patch = mpatches.Patch(color="#e74c3c", alpha=0.7, label="Neurovirulence targets")
    control_patch = mpatches.Patch(color="#3498db", alpha=0.7, label="Negative controls")
    ax1.legend(handles=[target_patch, control_patch], loc="lower right", fontsize=9)

    # Right panel: IDR boundary map (linear protein diagrams)
    ax2.set_title("IDR Boundary Maps", fontsize=12, fontweight="bold")
    y_pos = 0
    for entry in targets + controls:
        name = entry["protein"]
        length = entry["length"]
        idrs = entry["idr_boundaries"]
        is_target = entry["category"] == "target"
        color = "#e74c3c" if is_target else "#3498db"

        # Draw protein backbone
        ax2.plot([0, length], [y_pos, y_pos], color="gray", linewidth=3, alpha=0.3, solid_capstyle="round")

        # Draw IDRs
        for start, end in idrs:
            ax2.plot([start, end], [y_pos, y_pos], color=color, linewidth=6, alpha=0.6, solid_capstyle="round")

        ax2.text(-30, y_pos, name, fontsize=8, ha="right", va="center",
                fontweight="bold" if is_target else "normal")
        ax2.text(length + 10, y_pos, f"{length} aa", fontsize=7, ha="left", va="center", color="gray")

        y_pos -= 1

    ax2.set_xlim(-150, max(d["length"] for d in data) + 80)
    ax2.set_ylim(y_pos - 0.5, 0.5)
    ax2.set_xlabel("Residue Position", fontsize=11)
    ax2.set_yticks([])
    ax2.spines["left"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    print("Generating disorder profile plots...")

    # Individual target plots
    for name in ["ICP34.5", "US3", "pUL37"]:
        csv_path = os.path.join(PRED_DIR, f"{name}_per_residue.csv")
        if os.path.exists(csv_path):
            plot_path = os.path.join(PLOT_DIR, f"{name}_disorder_profile.png")
            plot_target_disorder(name, csv_path, plot_path)

    # Summary comparison plot
    results_path = os.path.join(BASE_DIR, "results", "disorder_summary_v2.json")
    if os.path.exists(results_path):
        summary_path = os.path.join(PLOT_DIR, "disorder_comparison_summary.png")
        plot_comparison_summary(results_path, summary_path)

    print("\nAll plots generated.")


if __name__ == "__main__":
    main()
