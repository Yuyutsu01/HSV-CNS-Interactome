"""
Integrate real IUPred3 scores from downloaded text files into the disorder prediction pipeline.
Rebuilds consensus maps, per-residue CSVs, summary JSON, and regenerates all plots.
"""
import os
import csv
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BASE_DIR = r"C:\Users\shiva\OneDrive\Desktop\HSV-CNS-Interactome\phase1"
IUPRED_DIR = os.path.join(BASE_DIR, "iupred3_data")
SEQ_DIR = os.path.join(BASE_DIR, "sequences")
PRED_DIR = os.path.join(BASE_DIR, "disorder_predictions")
PLOT_DIR = os.path.join(BASE_DIR, "plots")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

for d in [PRED_DIR, PLOT_DIR, RESULTS_DIR]:
    os.makedirs(d, exist_ok=True)

# Map filenames to protein names
FILE_MAP = {
    "ICP34.5_iupred3.txt.txt": ("ICP34.5", "P36313", "target"),
    "US3_iupred3.txt.txt": ("US3", "P04413", "target"),
    "pUL37_iupred3.txt.txt": ("pUL37", "P10221", "target"),
    "ICP8_UL29_iupred.txt.txt": ("ICP8_UL29", "P04296", "control"),
    "UL42_iupred3.txt.txt": ("UL42", "P10186", "control"),
    "gD_US6_iupred3.txt.txt": ("gD_US6", "P57083", "control"),
    "VP11-12_UL46_iupred3.txt.txt": ("VP11_12_UL46", "P10228", "control"),
    "VP16_UL48_iupred3.txt.txt": ("VP16_UL48", "P06492", "control"),
    "VP22_UL49_iupred3.txt.txt": ("VP22_UL49", "P10216", "control"),
    "ICP27_UL54_iupred3.txt.txt": ("ICP27_UL54", "P10238", "control"),
    "ICP22_US1_iupred3.txt.txt": ("ICP22_US1", "P04485", "control"),
}

# TOP-IDP propensity scale (Campen et al. 2008)
PROPENSITY = {
    'A': 0.06, 'R': 0.18, 'N': -0.04, 'D': 0.10, 'C': -0.20,
    'Q': 0.18, 'E': 0.24, 'G': 0.17, 'H': -0.01, 'I': -0.49,
    'L': -0.34, 'K': 0.18, 'M': -0.19, 'F': -0.42, 'P': 0.34,
    'S': 0.15, 'T': -0.03, 'W': -0.44, 'Y': -0.27, 'V': -0.38,
}

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

UNIPROT_DISORDER = {
    "ICP34.5": [(1, 129), (149, 174)],
    "US3": [(12, 63)],
    "pUL37": [(1, 27), (1047, 1123)],
}


def parse_iupred3_file(filepath):
    """Parse IUPred3 text output file."""
    positions = []
    residues = []
    scores = []
    with open(filepath, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                try:
                    pos = int(parts[0])
                    aa = parts[1]
                    score = float(parts[2])
                    positions.append(pos)
                    residues.append(aa)
                    scores.append(score)
                except ValueError:
                    continue
    return positions, residues, scores


def composition_disorder(sequence, window=21):
    """Composition-based disorder prediction using TOP-IDP scale."""
    seq = sequence.upper()
    n = len(seq)
    raw = np.array([PROPENSITY.get(aa, 0.0) for aa in seq])
    scores = np.zeros(n)
    half_w = window // 2
    for i in range(n):
        start = max(0, i - half_w)
        end = min(n, i + half_w + 1)
        scores[i] = np.mean(raw[start:end])
    scores = 1.0 / (1.0 + np.exp(-10 * scores))
    return scores.tolist()


def extract_idr_boundaries(consensus, min_length=20):
    """Extract contiguous IDR boundaries (>= min_length residues)."""
    idrs = []
    start = None
    for i, val in enumerate(consensus):
        if val == 1 and start is None:
            start = i
        elif val == 0 and start is not None:
            if i - start >= min_length:
                idrs.append((start + 1, i))
            start = None
    if start is not None and len(consensus) - start >= min_length:
        idrs.append((start + 1, len(consensus)))
    return idrs


def process_protein(name, accession, category, iupred_positions, iupred_residues, iupred_scores):
    """Build full disorder prediction using real IUPred3 + composition."""
    n = len(iupred_scores)
    sequence = "".join(iupred_residues)

    # Method 1: Real IUPred3 scores
    iupred = np.array(iupred_scores)

    # Method 2: Composition-based (TOP-IDP)
    comp = np.array(composition_disorder(sequence))

    # Consensus: disordered if BOTH methods agree (2/2)
    consensus = np.zeros(n, dtype=int)
    for i in range(n):
        if iupred[i] > 0.5 and comp[i] > 0.5:
            consensus[i] = 1

    idrs = extract_idr_boundaries(consensus.tolist())
    total_idr = sum(end - start + 1 for start, end in idrs)
    idr_frac = round(total_idr / n, 3)

    # Save per-residue CSV
    csv_path = os.path.join(PRED_DIR, f"{name}_per_residue.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["position", "residue", "iupred3", "composition", "consensus_disordered"])
        for i in range(n):
            writer.writerow([
                iupred_positions[i],
                iupred_residues[i],
                f"{iupred[i]:.6f}",
                f"{comp[i]:.6f}",
                consensus[i],
            ])

    return {
        "protein": name,
        "accession": accession,
        "category": category,
        "length": n,
        "idr_boundaries": idrs,
        "idr_fraction": idr_frac,
        "n_idrs": len(idrs),
        "iupred3_mean": round(float(np.mean(iupred)), 3),
        "composition_mean": round(float(np.mean(comp)), 3),
        "iupred3_disordered_frac": round(float(np.mean(iupred > 0.5)), 3),
    }


def plot_target_disorder(name, csv_path, output_path):
    """Create detailed disorder profile plot for one target protein."""
    df = pd.read_csv(csv_path)
    n = len(df)
    positions = df["position"].values

    fig, axes = plt.subplots(3, 1, figsize=(14, 8), sharex=True,
                              gridspec_kw={"height_ratios": [3, 1, 1]})
    fig.suptitle(f"{name} — Disorder Prediction Profile (IUPred3 + Composition)",
                 fontsize=14, fontweight="bold", y=0.95)

    ax1 = axes[0]
    ax1.plot(positions, df["iupred3"], color="#e74c3c", linewidth=1.5,
             label="IUPred3 (long)", alpha=0.9)
    ax1.plot(positions, df["composition"], color="#2980b9", linewidth=1.2,
             label="Composition (TOP-IDP)", alpha=0.8)
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
    for i in range(n):
        if consensus[i] == 1:
            ax2.axvspan(positions[i] - 0.5, positions[i] + 0.5,
                       color="#3498db", alpha=0.6)

    if name in UNIPROT_DISORDER:
        for start, end in UNIPROT_DISORDER[name]:
            ax2.axvspan(start - 0.5, end + 0.5, ymin=0.0, ymax=0.4,
                       color="#e74c3c", alpha=0.4)

    ax2.set_ylabel("IDR Map", fontsize=11)
    ax2.set_yticks([])
    ax2.set_ylim(0, 1)
    predicted_patch = mpatches.Patch(color="#3498db", alpha=0.6, label="Predicted IDR (consensus)")
    uniprot_patch = mpatches.Patch(color="#e74c3c", alpha=0.4, label="UniProt annotated disorder")
    ax2.legend(handles=[predicted_patch, uniprot_patch], loc="upper right", fontsize=8)

    # Panel 3: Functional regions
    ax3 = axes[2]
    if name in FUNCTIONAL_REGIONS:
        for start, end, label, color in FUNCTIONAL_REGIONS[name]:
            if start == end:
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


def plot_comparison_summary(all_results, output_path):
    """Summary comparison plot of IDR fractions across all proteins."""
    targets = [d for d in all_results if d["category"] == "target"]
    controls = [d for d in all_results if d["category"] == "control"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7),
                                     gridspec_kw={"width_ratios": [1, 2.5]})

    all_entries = targets + controls
    all_names = [d["protein"] for d in all_entries]
    all_fracs = [d["idr_fraction"] for d in all_entries]
    colors = ["#e74c3c" if d["category"] == "target" else "#3498db" for d in all_entries]

    bars = ax1.barh(range(len(all_names)), all_fracs, color=colors, alpha=0.7, edgecolor="white")
    ax1.set_yticks(range(len(all_names)))
    ax1.set_yticklabels(all_names, fontsize=9)
    ax1.set_xlabel("IDR Fraction (consensus)", fontsize=11)
    ax1.set_title("Predicted IDR Content\n(IUPred3 + Composition consensus)", fontsize=11, fontweight="bold")
    ax1.invert_yaxis()
    ax1.axvline(x=0.5, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)

    target_patch = mpatches.Patch(color="#e74c3c", alpha=0.7, label="Neurovirulence targets")
    control_patch = mpatches.Patch(color="#3498db", alpha=0.7, label="Negative controls")
    ax1.legend(handles=[target_patch, control_patch], loc="lower right", fontsize=9)

    # Right panel: IDR boundary map
    ax2.set_title("IDR Boundary Maps", fontsize=12, fontweight="bold")
    y_pos = 0
    for entry in all_entries:
        name = entry["protein"]
        length = entry["length"]
        idrs = entry["idr_boundaries"]
        is_target = entry["category"] == "target"
        color = "#e74c3c" if is_target else "#3498db"

        ax2.plot([0, length], [y_pos, y_pos], color="gray", linewidth=3,
                 alpha=0.3, solid_capstyle="round")
        for start, end in idrs:
            ax2.plot([start, end], [y_pos, y_pos], color=color, linewidth=6,
                     alpha=0.6, solid_capstyle="round")

        ax2.text(-30, y_pos, name, fontsize=8, ha="right", va="center",
                fontweight="bold" if is_target else "normal")
        ax2.text(length + 10, y_pos, f"{length} aa", fontsize=7, ha="left",
                 va="center", color="gray")
        y_pos -= 1

    ax2.set_xlim(-180, max(d["length"] for d in all_entries) + 80)
    ax2.set_ylim(y_pos - 0.5, 0.5)
    ax2.set_xlabel("Residue Position", fontsize=11)
    ax2.set_yticks([])
    ax2.spines["left"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


def main():
    print("=" * 70)
    print("INTEGRATING REAL IUPred3 SCORES INTO DISORDER PREDICTIONS")
    print("=" * 70)

    all_results = []

    for filename, (name, accession, category) in FILE_MAP.items():
        filepath = os.path.join(IUPRED_DIR, filename)
        if not os.path.exists(filepath):
            print(f"\n  WARNING: {filename} not found, skipping {name}")
            continue

        positions, residues, scores = parse_iupred3_file(filepath)
        print(f"\n  {name} ({accession}): {len(scores)} residues from IUPred3")

        result = process_protein(name, accession, category, positions, residues, scores)
        all_results.append(result)

        disordered_pct = result["iupred3_disordered_frac"] * 100
        consensus_pct = result["idr_fraction"] * 100
        print(f"    IUPred3 > 0.5: {disordered_pct:.1f}%")
        print(f"    Consensus IDR: {consensus_pct:.1f}% ({result['n_idrs']} regions)")
        for s, e in result["idr_boundaries"]:
            print(f"      {s}-{e} ({e - s + 1} residues)")

    # Save summary JSON
    summary_path = os.path.join(RESULTS_DIR, "disorder_summary_v3.json")
    with open(summary_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n  Summary saved: {summary_path}")

    # Save summary CSV
    csv_summary_path = os.path.join(RESULTS_DIR, "disorder_summary_v3.csv")
    with open(csv_summary_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["protein", "accession", "category", "length", "n_idrs",
                         "idr_fraction", "iupred3_mean", "iupred3_disordered_frac",
                         "idr_boundaries"])
        for r in all_results:
            boundaries = "; ".join(f"{s}-{e}" for s, e in r["idr_boundaries"])
            writer.writerow([r["protein"], r["accession"], r["category"], r["length"],
                            r["n_idrs"], r["idr_fraction"], r["iupred3_mean"],
                            r["iupred3_disordered_frac"], boundaries])

    # Generate plots
    print("\n  Generating plots...")
    for name in ["ICP34.5", "US3", "pUL37"]:
        csv_path = os.path.join(PRED_DIR, f"{name}_per_residue.csv")
        if os.path.exists(csv_path):
            plot_path = os.path.join(PLOT_DIR, f"{name}_disorder_profile.png")
            plot_target_disorder(name, csv_path, plot_path)
            print(f"    {plot_path}")

    if all_results:
        summary_plot = os.path.join(PLOT_DIR, "disorder_comparison_summary.png")
        plot_comparison_summary(all_results, summary_plot)
        print(f"    {summary_plot}")

    # Final summary table
    print(f"\n{'='*80}")
    print("DISORDER PREDICTION SUMMARY (IUPred3 + Composition consensus)")
    print(f"{'='*80}")
    print(f"{'Protein':<18} {'Cat':<8} {'Len':>5} {'IUPred3':>8} {'Consens':>8} {'#IDRs':>6} {'IDR Boundaries'}")
    print("-" * 80)
    for r in all_results:
        boundaries = ", ".join(f"{s}-{e}" for s, e in r["idr_boundaries"]) or "None"
        print(f"{r['protein']:<18} {r['category']:<8} {r['length']:>5} "
              f"{r['iupred3_disordered_frac']:>7.1%} {r['idr_fraction']:>7.1%} "
              f"{r['n_idrs']:>6} {boundaries}")

    # Stats comparison
    target_fracs = [r["idr_fraction"] for r in all_results if r["category"] == "target"]
    control_fracs = [r["idr_fraction"] for r in all_results if r["category"] == "control"]
    if target_fracs and control_fracs:
        print(f"\n  Target mean IDR fraction:  {np.mean(target_fracs):.3f}")
        print(f"  Control mean IDR fraction: {np.mean(control_fracs):.3f}")
        print(f"  Ratio (target/control):    {np.mean(target_fracs)/np.mean(control_fracs):.2f}x")


if __name__ == "__main__":
    main()
