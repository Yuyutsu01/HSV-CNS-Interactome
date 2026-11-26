#!/usr/bin/env python3
"""
Phase 2 Conservation Analysis: MSA-based conservation scoring using
Shannon entropy, Valdar weighting, and dN/dS approximation for
SLiM-containing regions across alphaherpesvirus orthologs.
"""

import os, sys, json, csv, re, math
import numpy as np
from collections import Counter

BASE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(BASE)
PHASE1 = os.path.join(PROJECT, "phase1")
ORTHO_DIR = os.path.join(PHASE1, "sequences", "orthologs")
RESULTS = os.path.join(BASE, "results")
os.makedirs(RESULTS, exist_ok=True)

# Ortholog groupings (protein → list of ortholog FASTA files)
ORTHOLOG_MAP = {
    "ICP34.5": {
        "reference_fasta": os.path.join(PHASE1, "sequences", "targets", "ICP34.5.fasta"),
        "orthologs": [
            ("HSV-2_RL1", os.path.join(ORTHO_DIR, "HSV2_RL1.fasta")),
        ],
    },
    "US3": {
        "reference_fasta": os.path.join(PHASE1, "sequences", "targets", "US3.fasta"),
        "orthologs": [
            ("HSV-2_US3", os.path.join(ORTHO_DIR, "HSV2_US3.fasta")),
            ("PRV_US3", os.path.join(ORTHO_DIR, "PRV_US3.fasta")),
            ("VZV_ORF66", os.path.join(ORTHO_DIR, "VZV_ORF66.fasta")),
            ("BoHV1_US3", os.path.join(ORTHO_DIR, "BoHV1_US3.fasta")),
            ("EHV1_US3", os.path.join(ORTHO_DIR, "EHV1_US3.fasta")),
        ],
    },
    "pUL37": {
        "reference_fasta": os.path.join(PHASE1, "sequences", "targets", "pUL37.fasta"),
        "orthologs": [
            ("HSV-2_UL37", os.path.join(ORTHO_DIR, "HSV2_UL37.fasta")),
            ("PRV_UL37", os.path.join(ORTHO_DIR, "PRV_UL37.fasta")),
            ("VZV_ORF21", os.path.join(ORTHO_DIR, "VZV_ORF21.fasta")),
        ],
    },
}

# BLOSUM62 substitution matrix (subset for conservation weighting)
BLOSUM62 = {
    ('A','A'): 4, ('A','R'):-1, ('A','N'):-2, ('A','D'):-2, ('A','C'): 0,
    ('A','Q'):-1, ('A','E'):-1, ('A','G'): 0, ('A','H'):-2, ('A','I'):-1,
    ('A','L'):-1, ('A','K'):-1, ('A','M'):-1, ('A','F'):-2, ('A','P'):-1,
    ('A','S'): 1, ('A','T'): 0, ('A','W'):-3, ('A','Y'):-2, ('A','V'): 0,
    ('R','R'): 5, ('R','N'): 0, ('R','D'):-2, ('R','C'):-3, ('R','Q'): 1,
    ('R','E'): 0, ('R','G'):-2, ('R','H'): 0, ('R','I'):-3, ('R','L'):-2,
    ('R','K'): 2, ('R','M'):-1, ('R','F'):-3, ('R','P'):-2, ('R','S'):-1,
    ('R','T'):-1, ('R','W'):-3, ('R','Y'):-2, ('R','V'):-3,
    ('N','N'): 6, ('N','D'): 1, ('N','C'):-3, ('N','Q'): 0, ('N','E'): 0,
    ('N','G'): 0, ('N','H'): 1, ('N','I'):-3, ('N','L'):-3, ('N','K'): 0,
    ('N','M'):-2, ('N','F'):-3, ('N','P'):-2, ('N','S'): 1, ('N','T'): 0,
    ('N','W'):-4, ('N','Y'):-2, ('N','V'):-3,
    ('D','D'): 6, ('D','C'):-3, ('D','Q'): 0, ('D','E'): 2, ('D','G'):-1,
    ('D','H'):-1, ('D','I'):-3, ('D','L'):-4, ('D','K'):-1, ('D','M'):-3,
    ('D','F'):-3, ('D','P'):-1, ('D','S'): 0, ('D','T'):-1, ('D','W'):-4,
    ('D','Y'):-3, ('D','V'):-3,
    ('C','C'): 9, ('C','Q'):-3, ('C','E'):-4, ('C','G'):-3, ('C','H'):-3,
    ('C','I'):-1, ('C','L'):-1, ('C','K'):-3, ('C','M'):-1, ('C','F'):-2,
    ('C','P'):-3, ('C','S'):-1, ('C','T'):-1, ('C','W'):-2, ('C','Y'):-2,
    ('C','V'):-1,
    ('Q','Q'): 5, ('Q','E'): 2, ('Q','G'):-2, ('Q','H'): 0, ('Q','I'):-3,
    ('Q','L'):-2, ('Q','K'): 1, ('Q','M'): 0, ('Q','F'):-3, ('Q','P'):-1,
    ('Q','S'): 0, ('Q','T'):-1, ('Q','W'):-2, ('Q','Y'):-1, ('Q','V'):-2,
    ('E','E'): 5, ('E','G'):-2, ('E','H'): 0, ('E','I'):-3, ('E','L'):-3,
    ('E','K'): 1, ('E','M'):-2, ('E','F'):-3, ('E','P'):-1, ('E','S'): 0,
    ('E','T'):-1, ('E','W'):-3, ('E','Y'):-2, ('E','V'):-2,
    ('G','G'): 6, ('G','H'):-2, ('G','I'):-4, ('G','L'):-4, ('G','K'):-2,
    ('G','M'):-3, ('G','F'):-3, ('G','P'):-2, ('G','S'): 0, ('G','T'):-2,
    ('G','W'):-2, ('G','Y'):-3, ('G','V'):-3,
    ('H','H'): 8, ('H','I'):-3, ('H','L'):-3, ('H','K'):-1, ('H','M'):-2,
    ('H','F'):-1, ('H','P'):-2, ('H','S'):-1, ('H','T'):-2, ('H','W'):-2,
    ('H','Y'): 2, ('H','V'):-3,
    ('I','I'): 4, ('I','L'): 2, ('I','K'):-3, ('I','M'): 1, ('I','F'): 0,
    ('I','P'):-3, ('I','S'):-2, ('I','T'):-1, ('I','W'):-3, ('I','Y'):-1,
    ('I','V'): 3,
    ('L','L'): 4, ('L','K'):-2, ('L','M'): 2, ('L','F'): 0, ('L','P'):-3,
    ('L','S'):-2, ('L','T'):-1, ('L','W'):-2, ('L','Y'):-1, ('L','V'): 1,
    ('K','K'): 5, ('K','M'):-1, ('K','F'):-3, ('K','P'):-1, ('K','S'): 0,
    ('K','T'):-1, ('K','W'):-3, ('K','Y'):-2, ('K','V'):-2,
    ('M','M'): 5, ('M','F'): 0, ('M','P'):-2, ('M','S'):-1, ('M','T'):-1,
    ('M','W'):-1, ('M','Y'):-1, ('M','V'): 1,
    ('F','F'): 6, ('F','P'):-4, ('F','S'):-2, ('F','T'):-2, ('F','W'): 1,
    ('F','Y'): 3, ('F','V'):-1,
    ('P','P'): 7, ('P','S'):-1, ('P','T'):-1, ('P','W'):-4, ('P','Y'):-3,
    ('P','V'):-2,
    ('S','S'): 4, ('S','T'): 1, ('S','W'):-3, ('S','Y'):-2, ('S','V'):-2,
    ('T','T'): 5, ('T','W'):-2, ('T','Y'):-2, ('T','V'): 0,
    ('W','W'):11, ('W','Y'): 2, ('W','V'):-3,
    ('Y','Y'): 7, ('Y','V'):-1,
    ('V','V'): 4,
}

def blosum_score(a, b):
    if a == '-' or b == '-':
        return -4
    pair = (a.upper(), b.upper())
    if pair in BLOSUM62:
        return BLOSUM62[pair]
    pair = (b.upper(), a.upper())
    if pair in BLOSUM62:
        return BLOSUM62[pair]
    return -1


def load_fasta(filepath):
    if not os.path.exists(filepath):
        return None
    seq = []
    with open(filepath) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)


# ═══════════════════════════════════════════════════════════════════════
# PAIRWISE ALIGNMENT (Needleman-Wunsch)
# ═══════════════════════════════════════════════════════════════════════

def needleman_wunsch(seq1, seq2, gap_open=-10, gap_extend=-1):
    """Simple Needleman-Wunsch global alignment with affine gap penalties."""
    m, n = len(seq1), len(seq2)
    # Use simpler linear gap for speed with long sequences
    gap_pen = -2

    # For very long sequences, use Biopython if available
    if m > 500 or n > 500:
        try:
            from Bio import pairwise2
            alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -5, -0.5,
                                                   one_alignment_only=True)
            if alignments:
                return alignments[0].seqA, alignments[0].seqB, alignments[0].score
        except:
            pass

    # Simple NW for shorter sequences
    score = np.zeros((m + 1, n + 1))
    for i in range(1, m + 1):
        score[i][0] = i * gap_pen
    for j in range(1, n + 1):
        score[0][j] = j * gap_pen

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i-1][j-1] + blosum_score(seq1[i-1], seq2[j-1])
            delete = score[i-1][j] + gap_pen
            insert = score[i][j-1] + gap_pen
            score[i][j] = max(match, delete, insert)

    # Traceback
    aln1, aln2 = [], []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score[i][j] == score[i-1][j-1] + blosum_score(seq1[i-1], seq2[j-1]):
            aln1.append(seq1[i-1])
            aln2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and score[i][j] == score[i-1][j] + gap_pen:
            aln1.append(seq1[i-1])
            aln2.append('-')
            i -= 1
        else:
            aln1.append('-')
            aln2.append(seq2[j-1])
            j -= 1

    return ''.join(reversed(aln1)), ''.join(reversed(aln2)), score[m][n]


# ═══════════════════════════════════════════════════════════════════════
# MSA FROM PAIRWISE ALIGNMENTS (star alignment)
# ═══════════════════════════════════════════════════════════════════════

def star_msa(ref_seq, ortholog_seqs):
    """Build MSA using star alignment with reference as center."""
    if not ortholog_seqs:
        return [ref_seq], [ref_seq]

    pairwise_alns = []
    for name, seq in ortholog_seqs:
        aln_ref, aln_orth, score = needleman_wunsch(ref_seq, seq)
        pairwise_alns.append((name, aln_ref, aln_orth, score))
        print(f"    Aligned {name}: score={score:.0f}, identity={sum(1 for a,b in zip(aln_ref, aln_orth) if a==b and a!='-')/len(ref_seq)*100:.1f}%")

    # For star MSA: merge gap positions from all pairwise alignments
    # Simple approach: use the reference alignment positions as guide
    # This gives per-position conservation relative to the reference
    return pairwise_alns


# ═══════════════════════════════════════════════════════════════════════
# CONSERVATION METRICS
# ═══════════════════════════════════════════════════════════════════════

def shannon_entropy(column):
    """Shannon entropy for an alignment column (excluding gaps)."""
    residues = [r for r in column if r != '-']
    if len(residues) < 2:
        return 0.0
    counts = Counter(residues)
    n = len(residues)
    H = 0
    for count in counts.values():
        p = count / n
        if p > 0:
            H -= p * math.log2(p)
    return H


def valdar_conservation(column, weights=None):
    """Valdar conservation score using BLOSUM62-weighted identity."""
    residues = [(i, r) for i, r in enumerate(column) if r != '-']
    if len(residues) < 2:
        return 1.0
    n = len(residues)
    total = 0
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            w_i = weights[residues[i][0]] if weights else 1.0
            w_j = weights[residues[j][0]] if weights else 1.0
            s = blosum_score(residues[i][1], residues[j][1])
            # Normalize by self-score
            max_s = max(blosum_score(residues[i][1], residues[i][1]),
                       blosum_score(residues[j][1], residues[j][1]))
            if max_s > 0:
                norm_s = s / max_s
            else:
                norm_s = 0
            total += w_i * w_j * norm_s
            count += w_i * w_j
    return total / count if count > 0 else 0.0


def compute_per_position_conservation(ref_seq, pairwise_alns):
    """Compute per-position conservation metrics from pairwise alignments."""
    n = len(ref_seq)
    conservation = np.zeros(n)
    identity = np.zeros(n)
    n_aligned = np.zeros(n)

    for name, aln_ref, aln_orth, score in pairwise_alns:
        # Map alignment positions back to reference positions
        ref_pos = -1
        for aln_pos in range(len(aln_ref)):
            if aln_ref[aln_pos] != '-':
                ref_pos += 1
                if ref_pos < n:
                    n_aligned[ref_pos] += 1
                    if aln_ref[aln_pos] == aln_orth[aln_pos]:
                        identity[ref_pos] += 1
                    # BLOSUM62-weighted similarity
                    s = blosum_score(aln_ref[aln_pos], aln_orth[aln_pos])
                    max_s = blosum_score(aln_ref[aln_pos], aln_ref[aln_pos])
                    if max_s > 0:
                        conservation[ref_pos] += s / max_s

    # Normalize
    for i in range(n):
        if n_aligned[i] > 0:
            conservation[i] /= n_aligned[i]
            identity[i] /= n_aligned[i]
        else:
            conservation[i] = 0.0
            identity[i] = 0.0

    return conservation, identity, n_aligned


def compute_slim_conservation(ref_seq, conservation_scores, identity_scores, slim_hits):
    """Compute conservation metrics specifically for SLiM regions."""
    n = len(ref_seq)
    results = []
    for hit in slim_hits:
        s = hit["start"] - 1  # 0-indexed
        e = hit["end"]
        if e > n:
            e = n
        if s < 0:
            s = 0

        region_cons = conservation_scores[s:e]
        region_id = identity_scores[s:e]

        # Flanking conservation (10 residues each side)
        flank = 10
        flank_s = max(0, s - flank)
        flank_e = min(n, e + flank)
        flank_cons = np.concatenate([conservation_scores[flank_s:s], conservation_scores[e:flank_e]])

        result = {
            "motif_class": hit["motif_class"],
            "start": hit["start"],
            "end": hit["end"],
            "mean_conservation": round(float(np.mean(region_cons)), 4) if len(region_cons) > 0 else 0,
            "mean_identity": round(float(np.mean(region_id)), 4) if len(region_id) > 0 else 0,
            "min_conservation": round(float(np.min(region_cons)), 4) if len(region_cons) > 0 else 0,
            "flanking_conservation": round(float(np.mean(flank_cons)), 4) if len(flank_cons) > 0 else 0,
            "conservation_ratio": round(float(np.mean(region_cons)) / (float(np.mean(flank_cons)) + 0.001), 3) if len(region_cons) > 0 and len(flank_cons) > 0 else 0,
        }

        # Conservation classification
        mc = result["mean_conservation"]
        mi = result["mean_identity"]
        if mi >= 0.8:
            result["conservation_class"] = "highly_conserved"
        elif mi >= 0.5:
            result["conservation_class"] = "moderately_conserved"
        elif mi >= 0.2:
            result["conservation_class"] = "weakly_conserved"
        else:
            result["conservation_class"] = "not_conserved"

        results.append(result)
    return results


# ═══════════════════════════════════════════════════════════════════════
# dN/dS APPROXIMATION
# ═══════════════════════════════════════════════════════════════════════

def approximate_dnds(ref_seq, pairwise_alns, window=30):
    """
    Approximate dN/dS using protein-level conservation as proxy.
    True dN/dS requires codon alignments; this uses the ratio of
    non-synonymous-like changes to synonymous-like changes based on
    amino acid property conservation.
    """
    n = len(ref_seq)
    # Classify substitutions by biochemical similarity
    SIMILAR_GROUPS = [
        set("AILMV"),  # nonpolar aliphatic
        set("FWY"),    # aromatic
        set("DE"),     # acidic
        set("RKH"),    # basic
        set("STNQ"),   # polar uncharged
        set("CGP"),    # special
    ]

    per_position_ratio = np.zeros(n)
    n_subs = np.zeros(n)

    for name, aln_ref, aln_orth, score in pairwise_alns:
        ref_pos = -1
        for aln_pos in range(len(aln_ref)):
            if aln_ref[aln_pos] != '-':
                ref_pos += 1
                if ref_pos >= n:
                    break
                if aln_orth[aln_pos] == '-':
                    continue
                if aln_ref[aln_pos] != aln_orth[aln_pos]:
                    n_subs[ref_pos] += 1
                    # Check if substitution is conservative (synonymous-like)
                    aa1, aa2 = aln_ref[aln_pos], aln_orth[aln_pos]
                    conservative = any(aa1 in g and aa2 in g for g in SIMILAR_GROUPS)
                    if conservative:
                        per_position_ratio[ref_pos] += 0.3  # synonymous-like
                    else:
                        per_position_ratio[ref_pos] += 1.0  # non-synonymous-like

    # Smooth over windows
    smoothed = np.zeros(n)
    half_w = window // 2
    for i in range(n):
        s = max(0, i - half_w)
        e = min(n, i + half_w + 1)
        if np.sum(n_subs[s:e]) > 0:
            smoothed[i] = np.sum(per_position_ratio[s:e]) / np.sum(n_subs[s:e])
        else:
            smoothed[i] = 0.0  # No substitutions = perfectly conserved

    return smoothed, per_position_ratio, n_subs


# ═══════════════════════════════════════════════════════════════════════
# PLOTTING
# ═══════════════════════════════════════════════════════════════════════

def plot_conservation(protein_name, ref_seq, conservation, identity, dnds_ratio, idr_boundaries):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(ref_seq)
    positions = np.arange(1, n + 1)

    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    fig.suptitle(f"{protein_name} - Conservation Analysis", fontsize=14, fontweight="bold")

    # Panel 1: Conservation + Identity
    ax = axes[0]
    ax.plot(positions, conservation, label="BLOSUM62 conservation", color="#E74C3C", linewidth=1, alpha=0.8)
    ax.plot(positions, identity, label="Sequence identity", color="#3498DB", linewidth=1, alpha=0.8)
    for s, e in idr_boundaries:
        ax.axvspan(s, e, alpha=0.15, color="#F39C12")
    ax.set_ylabel("Score")
    ax.set_ylim(-0.1, 1.1)
    ax.legend(fontsize=8)
    ax.set_title("Per-Residue Conservation")

    # Panel 2: dN/dS approximation
    ax = axes[1]
    ax.fill_between(positions, dnds_ratio, alpha=0.5, color="#9B59B6")
    ax.plot(positions, dnds_ratio, color="#8E44AD", linewidth=1)
    ax.axhline(1.0, color="gray", linestyle="--", alpha=0.5, label="dN/dS = 1 (neutral)")
    for s, e in idr_boundaries:
        ax.axvspan(s, e, alpha=0.15, color="#F39C12")
    ax.set_ylabel("dN/dS proxy")
    ax.legend(fontsize=8)
    ax.set_title("Evolutionary Rate (dN/dS Approximation)")

    # Panel 3: Conservation heatmap
    ax = axes[2]
    # Color by conservation: green (high) → yellow → red (low)
    colors = []
    for i in range(n):
        c = conservation[i]
        if c >= 0.7:
            colors.append("#2ECC71")
        elif c >= 0.4:
            colors.append("#F39C12")
        else:
            colors.append("#E74C3C")
    ax.bar(positions, np.ones(n) * 0.8, color=colors, width=1.0, alpha=0.8)
    for s, e in idr_boundaries:
        ax.axvspan(s, e, alpha=0.15, color="#3498DB")
    ax.set_yticks([])
    ax.set_xlabel("Residue Position")
    ax.set_title("Conservation Heatmap (Green=High, Orange=Medium, Red=Low; Blue=IDR)")

    plt.tight_layout()
    out_path = os.path.join(RESULTS, f"{protein_name}_conservation.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Plot saved: {out_path}")
    return out_path


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    # Load existing slim hits
    slim_hits_path = os.path.join(RESULTS, "slim_hits_refined.csv")
    all_slim_hits = {}
    if os.path.exists(slim_hits_path):
        with open(slim_hits_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                prot = row["protein"]
                if prot not in all_slim_hits:
                    all_slim_hits[prot] = []
                all_slim_hits[prot].append({
                    "motif_class": row["motif_class"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                })

    # Load IDR boundaries
    with open(os.path.join(PHASE1, "results", "disorder_summary_v3.json")) as f:
        disorder_data = json.load(f)
    idr_map = {p["protein"]: p["idr_boundaries"] for p in disorder_data}

    all_conservation = {}
    all_slim_conservation = []

    for protein_name, ortho_info in ORTHOLOG_MAP.items():
        print(f"\n{'='*60}")
        print(f"Conservation analysis: {protein_name}")
        print(f"{'='*60}")

        ref_seq = load_fasta(ortho_info["reference_fasta"])
        if not ref_seq:
            print(f"  ERROR: Cannot load reference sequence")
            continue
        print(f"  Reference: {len(ref_seq)} residues")

        # Load orthologs
        ortholog_seqs = []
        for name, path in ortho_info["orthologs"]:
            seq = load_fasta(path)
            if seq:
                ortholog_seqs.append((name, seq))
                print(f"  Ortholog {name}: {len(seq)} residues")
            else:
                print(f"  WARNING: Cannot load {name}")

        if not ortholog_seqs:
            print(f"  No orthologs available, skipping")
            continue

        # Pairwise alignments
        print(f"  Running pairwise alignments...")
        pairwise_alns = star_msa(ref_seq, ortholog_seqs)

        # Per-position conservation
        conservation, identity, n_aligned = compute_per_position_conservation(ref_seq, pairwise_alns)
        print(f"  Mean conservation: {np.mean(conservation):.3f}")
        print(f"  Mean identity: {np.mean(identity):.3f}")

        # dN/dS approximation
        dnds_ratio, _, _ = approximate_dnds(ref_seq, pairwise_alns)
        print(f"  Mean dN/dS proxy: {np.mean(dnds_ratio):.3f}")

        # SLiM-specific conservation
        slim_hits = all_slim_hits.get(protein_name, [])
        if slim_hits:
            slim_cons = compute_slim_conservation(ref_seq, conservation, identity, slim_hits)
            for sc in slim_cons:
                sc["protein"] = protein_name
                all_slim_conservation.append(sc)
            print(f"  SLiM conservation: {len(slim_cons)} motifs scored")
            # Print top conserved
            sorted_sc = sorted(slim_cons, key=lambda x: x["mean_conservation"], reverse=True)
            for sc in sorted_sc[:5]:
                print(f"    {sc['motif_class']} ({sc['start']}-{sc['end']}): "
                      f"cons={sc['mean_conservation']:.3f}, id={sc['mean_identity']:.3f}, "
                      f"class={sc['conservation_class']}")

        # Store results
        all_conservation[protein_name] = {
            "length": len(ref_seq),
            "n_orthologs": len(ortholog_seqs),
            "mean_conservation": round(float(np.mean(conservation)), 4),
            "mean_identity": round(float(np.mean(identity)), 4),
            "mean_dnds_proxy": round(float(np.mean(dnds_ratio)), 4),
            "per_position_conservation": conservation.tolist(),
            "per_position_identity": identity.tolist(),
            "per_position_dnds": dnds_ratio.tolist(),
        }

        # IDR-specific conservation
        idr_bounds = idr_map.get(protein_name, [])
        idr_cons = []
        for s, e in idr_bounds:
            if s <= len(ref_seq) and e <= len(ref_seq):
                region_cons = conservation[s-1:e]
                region_id = identity[s-1:e]
                idr_cons.append({
                    "region": f"{s}-{e}",
                    "mean_conservation": round(float(np.mean(region_cons)), 4),
                    "mean_identity": round(float(np.mean(region_id)), 4),
                })
        all_conservation[protein_name]["idr_conservation"] = idr_cons

        # Plot
        plot_conservation(protein_name, ref_seq, conservation, identity, dnds_ratio, idr_bounds)

    # Save outputs
    # 1. Full conservation data
    cons_path = os.path.join(RESULTS, "conservation_analysis.json")
    cons_out = {k: {kk: vv for kk, vv in v.items()
                     if kk not in ["per_position_conservation", "per_position_identity", "per_position_dnds"]}
                for k, v in all_conservation.items()}
    with open(cons_path, "w") as f:
        json.dump(cons_out, f, indent=2)
    print(f"\nConservation summary: {cons_path}")

    # 2. SLiM conservation CSV
    if all_slim_conservation:
        slim_cons_path = os.path.join(RESULTS, "slim_conservation.csv")
        keys = all_slim_conservation[0].keys()
        with open(slim_cons_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(all_slim_conservation)
        print(f"SLiM conservation: {slim_cons_path} ({len(all_slim_conservation)} entries)")

    # 3. Per-position conservation CSV (for annotation matrix integration)
    pos_path = os.path.join(RESULTS, "per_position_conservation.csv")
    rows = []
    for prot, data in all_conservation.items():
        for i, (c, ident, d) in enumerate(zip(
            data.get("per_position_conservation", []),
            data.get("per_position_identity", []),
            data.get("per_position_dnds", [])
        )):
            rows.append({
                "protein": prot,
                "position": i + 1,
                "blosum62_conservation": round(c, 4),
                "sequence_identity": round(ident, 4),
                "dnds_proxy": round(d, 4),
            })
    if rows:
        with open(pos_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        print(f"Per-position conservation: {pos_path} ({len(rows)} rows)")

    # Summary
    print(f"\n{'='*60}")
    print("CONSERVATION SUMMARY")
    print(f"{'='*60}")
    for prot, data in cons_out.items():
        print(f"  {prot}: {data['n_orthologs']} orthologs, "
              f"mean_cons={data['mean_conservation']:.3f}, "
              f"mean_id={data['mean_identity']:.3f}, "
              f"dN/dS={data['mean_dnds_proxy']:.3f}")


if __name__ == "__main__":
    main()
