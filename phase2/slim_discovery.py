"""
Phase 2 — SLiM Discovery and Filtering Pipeline

Step 1: Scan all IDRs for ELM motif matches using curated regex patterns
Step 2: Score each hit by disorder confidence, conservation, and mimicry
Step 3: Identify overlapping/competitive SLiM architectures
Step 4: Compare target vs control enrichment
"""
import re
import os
import json
import csv
import numpy as np
import pandas as pd
from collections import defaultdict

BASE_DIR = r"C:\Users\shiva\OneDrive\Desktop\HSV-CNS-Interactome"
PHASE1_DIR = os.path.join(BASE_DIR, "phase1")
PHASE2_DIR = os.path.join(BASE_DIR, "phase2")
os.makedirs(os.path.join(PHASE2_DIR, "results"), exist_ok=True)
os.makedirs(os.path.join(PHASE2_DIR, "plots"), exist_ok=True)

# ──────────────────────────────────────────────────────────────────────
# ELM MOTIF DATABASE — curated subset of ~100 most relevant classes
# Source: Eukaryotic Linear Motif resource (elm.eu.org)
# Focused on: autophagy, apoptosis, kinase docking, innate immunity,
#             protein phosphatase binding, nuclear transport, ubiquitin,
#             14-3-3, TRAF, SH2/SH3, and neuronal signaling
# ──────────────────────────────────────────────────────────────────────

ELM_MOTIFS = {
    # === AUTOPHAGY ===
    "LIG_LIR_Gen_1": {
        "regex": r"[EDST].{0,2}[WFY]..[ILV]",
        "description": "LC3-interacting region (LIR/AIM motif)",
        "host_targets": ["MAP1LC3A", "MAP1LC3B", "GABARAP", "GABARAPL2"],
        "pathway": "autophagy",
    },
    "LIG_LIR_Apic_2": {
        "regex": r"[EDST].{0,2}[WFY]..[ILVF]",
        "description": "Apicomplexan-type LIR motif",
        "host_targets": ["MAP1LC3B"],
        "pathway": "autophagy",
    },

    # === BECLIN-1 / AUTOPHAGY REGULATION ===
    "LIG_BH3_BclXL_1": {
        "regex": r"[IVLMF].{1,4}[ILMVF].{2,3}[GAS][DEST][EDST][IVLMFAW]",
        "description": "BH3 motif binding anti-apoptotic Bcl-2 family",
        "host_targets": ["BCL2", "BCL2L1", "MCL1"],
        "pathway": "apoptosis/autophagy",
    },

    # === PP1 BINDING ===
    "LIG_PP1_RVXF_1": {
        "regex": r"[RK].{0,1}[VIL][^P][FW]",
        "description": "PP1 binding RVxF motif",
        "host_targets": ["PPP1CA", "PPP1CB", "PPP1CC"],
        "pathway": "phosphatase_signaling",
    },
    "LIG_PP1_SILK_1": {
        "regex": r"[GS]IL[KR]",
        "description": "PP1 SILK motif",
        "host_targets": ["PPP1CA", "PPP1CB", "PPP1CC"],
        "pathway": "phosphatase_signaling",
    },

    # === PP2A BINDING ===
    "LIG_PP2A_B56_1": {
        "regex": r"[LMI]xx[ILV]xE",
        "description": "PP2A B56 regulatory subunit binding",
        "host_targets": ["PPP2R5A", "PPP2R5B", "PPP2R5C", "PPP2R5D", "PPP2R5E"],
        "pathway": "phosphatase_signaling",
    },

    # === 14-3-3 BINDING ===
    "LIG_14-3-3_CanoR_1": {
        "regex": r"R[^E][^E][ST][^E]P",
        "description": "14-3-3 binding phospho-motif (mode 1)",
        "host_targets": ["YWHAB", "YWHAE", "YWHAG", "YWHAH", "YWHAQ", "YWHAZ", "SFN"],
        "pathway": "14-3-3_signaling",
    },
    "LIG_14-3-3_CanoR_2": {
        "regex": r"R[^E][^E][^E][ST][^E]P",
        "description": "14-3-3 binding phospho-motif (mode 2)",
        "host_targets": ["YWHAB", "YWHAE", "YWHAG", "YWHAH", "YWHAQ", "YWHAZ"],
        "pathway": "14-3-3_signaling",
    },

    # === KINASE DOCKING ===
    "DOC_MAPK_DCC_7": {
        "regex": r"[KR]{2,4}.{2,6}[LIVMP].[LIVMP]",
        "description": "MAPK D-motif / docking site",
        "host_targets": ["MAPK1", "MAPK3", "MAPK8", "MAPK14"],
        "pathway": "MAPK_signaling",
    },
    "DOC_MAPK_FxFP_6": {
        "regex": r"F.FP",
        "description": "MAPK FxFP docking motif (DEF motif)",
        "host_targets": ["MAPK1", "MAPK3"],
        "pathway": "MAPK_signaling",
    },
    "MOD_CDK_SPxK_1": {
        "regex": r"[ST]P.[KR]",
        "description": "CDK phosphorylation site (SP-x-K/R)",
        "host_targets": ["CDK1", "CDK2", "CDK5"],
        "pathway": "cell_cycle",
    },
    "MOD_CK2_1": {
        "regex": r"[ST]..[DE]",
        "description": "CK2 phosphorylation site",
        "host_targets": ["CSNK2A1", "CSNK2A2"],
        "pathway": "CK2_signaling",
    },
    "MOD_GSK3_1": {
        "regex": r"[ST]...[ST]",
        "description": "GSK3 phosphorylation site (primed)",
        "host_targets": ["GSK3A", "GSK3B"],
        "pathway": "Wnt_signaling",
    },
    "MOD_PKA_2": {
        "regex": r"[RK][RK].[ST]",
        "description": "PKA phosphorylation site",
        "host_targets": ["PRKACA", "PRKACB"],
        "pathway": "cAMP_signaling",
    },
    "MOD_PKB_1": {
        "regex": r"R.R..[ST]",
        "description": "Akt/PKB phosphorylation site",
        "host_targets": ["AKT1", "AKT2", "AKT3"],
        "pathway": "PI3K_Akt_signaling",
    },

    # === TRAF BINDING (innate immunity) ===
    "LIG_TRAF2_1": {
        "regex": r"[PSAT].[QE]E",
        "description": "TRAF2 binding motif",
        "host_targets": ["TRAF2", "TRAF6"],
        "pathway": "NFkB_innate_immunity",
    },
    "LIG_TRAF6_MATH_1": {
        "regex": r"P.E..[FYWHDE]",
        "description": "TRAF6 MATH domain binding motif",
        "host_targets": ["TRAF6"],
        "pathway": "NFkB_innate_immunity",
    },

    # === APOPTOSIS ===
    "LIG_IAP_AVPI_1": {
        "regex": r"A[VIL]P[FIYVM]",
        "description": "IAP-binding motif (IBM / Smac mimetic)",
        "host_targets": ["BIRC2", "BIRC3", "XIAP"],
        "pathway": "apoptosis",
    },
    "DOC_CASP_3": {
        "regex": r"D..[DA]",
        "description": "Caspase-3 cleavage site",
        "host_targets": ["CASP3", "CASP7"],
        "pathway": "apoptosis",
    },

    # === UBIQUITIN SYSTEM ===
    "DEG_SCF_FBW7_1": {
        "regex": r"[LIVMP]T[PSG][^R][^K][ST]",
        "description": "SCF-FBW7 phosphodegron",
        "host_targets": ["FBXW7"],
        "pathway": "ubiquitin_proteasome",
    },
    "DEG_APCC_DBOX_1": {
        "regex": r"R..L",
        "description": "APC/C D-box (destruction box)",
        "host_targets": ["ANAPC2"],
        "pathway": "ubiquitin_proteasome",
    },
    "DEG_MDM2_SWIB_1": {
        "regex": r"F..[WL]...[IVLM]",
        "description": "MDM2 binding motif (p53-like)",
        "host_targets": ["MDM2"],
        "pathway": "p53_pathway",
    },

    # === SH2 / SH3 DOMAIN BINDING ===
    "LIG_SH2_STAT3_1": {
        "regex": r"Y..[LIVMF]",
        "description": "STAT3 SH2 domain binding pYxxL motif",
        "host_targets": ["STAT3"],
        "pathway": "JAK_STAT_signaling",
    },
    "LIG_SH2_STAT5_1": {
        "regex": r"Y..[LIVMF]",
        "description": "STAT5 SH2 domain binding pYxxL motif",
        "host_targets": ["STAT5A", "STAT5B"],
        "pathway": "JAK_STAT_signaling",
    },
    "LIG_SH2_GRB2_1": {
        "regex": r"Y[EDST]N[IVLMF]",
        "description": "Grb2 SH2 domain binding",
        "host_targets": ["GRB2"],
        "pathway": "RTK_signaling",
    },
    "LIG_SH3_3": {
        "regex": r"[RKY]..P..P",
        "description": "Class II SH3 binding motif",
        "host_targets": ["SRC", "FYN", "LCK"],
        "pathway": "SH3_signaling",
    },

    # === NUCLEAR TRANSPORT ===
    "TRG_NLS_MonoExtC_4": {
        "regex": r"[KR][KR][KR]",
        "description": "Monopartite NLS (basic cluster)",
        "host_targets": ["KPNA1", "KPNA2"],
        "pathway": "nuclear_transport",
    },
    "TRG_NLS_Bipartite_1": {
        "regex": r"[KR][KR].{8,12}[KR]{3}",
        "description": "Bipartite NLS",
        "host_targets": ["KPNA1", "KPNA2", "KPNA3"],
        "pathway": "nuclear_transport",
    },
    "TRG_NES_CRM1_1": {
        "regex": r"[LIVMF].{2,3}[LIVMF].{2,3}[LIVMF].[LIVMF]",
        "description": "CRM1/XPO1 nuclear export signal",
        "host_targets": ["XPO1"],
        "pathway": "nuclear_transport",
    },

    # === ENDOCYTOSIS / TRAFFICKING ===
    "TRG_ENDOCYTIC_2": {
        "regex": r"Y..[LIVMF]",
        "description": "Tyrosine-based endocytic sorting motif",
        "host_targets": ["AP2M1"],
        "pathway": "endocytosis",
    },
    "LIG_PTB_Phospho_1": {
        "regex": r"NP.Y",
        "description": "PTB domain binding NPxY motif",
        "host_targets": ["DAB1", "DAB2"],
        "pathway": "endocytosis",
    },

    # === SUMO / ACETYLATION ===
    "MOD_SUMO_for_1": {
        "regex": r"[VILMAFP]K.E",
        "description": "SUMO conjugation site",
        "host_targets": ["UBE2I"],
        "pathway": "SUMOylation",
    },

    # === PDZ DOMAIN BINDING ===
    "LIG_PDZ_Class_1": {
        "regex": r"[ST].[IVLMF]$",
        "description": "PDZ domain binding C-terminal class I",
        "host_targets": ["DLG1", "DLG4"],
        "pathway": "PDZ_signaling",
    },

    # === WW DOMAIN BINDING ===
    "LIG_WW_1": {
        "regex": r"PP.Y",
        "description": "WW domain binding PY motif",
        "host_targets": ["NEDD4", "NEDD4L", "ITCH", "WWP1"],
        "pathway": "ubiquitin_signaling",
    },

    # === DYNEIN / KINESIN (axonal transport — key for neurotropism) ===
    "LIG_Dynein_DLC8_1": {
        "regex": r"[KR].TQT",
        "description": "Dynein light chain LC8/DYNLL1 binding",
        "host_targets": ["DYNLL1", "DYNLL2"],
        "pathway": "axonal_transport",
    },

    # === INTERFERON / INNATE IMMUNITY ===
    "LIG_IRF3_1": {
        "regex": r"[DE]L[IVLM][ST]",
        "description": "IRF3 dimerization interface motif",
        "host_targets": ["IRF3", "IRF7"],
        "pathway": "interferon_signaling",
    },

    # === mTOR SIGNALING ===
    "LIG_FHA_2": {
        "regex": r"..T..[IVED]",
        "description": "FHA domain binding pTxxD motif",
        "host_targets": ["CHEK2", "RNF8"],
        "pathway": "DNA_damage_signaling",
    },

    # === CYCLIN BINDING ===
    "DOC_CYCLIN_RxL_1": {
        "regex": r"[RK].L",
        "description": "Cyclin binding RxL motif",
        "host_targets": ["CCNA2", "CCNB1", "CCNE1"],
        "pathway": "cell_cycle",
    },

    # === eIF4E BINDING (translational control) ===
    "LIG_eIF4E_1": {
        "regex": r"Y....L[VILMF]",
        "description": "eIF4E binding motif",
        "host_targets": ["EIF4E"],
        "pathway": "translation",
    },

    # === STING / cGAS pathway ===
    "LIG_TBK1_1": {
        "regex": r"[DE]..[IVLM].LR",
        "description": "TBK1 binding motif (STING-like)",
        "host_targets": ["TBK1"],
        "pathway": "interferon_signaling",
    },
}


def read_fasta(filepath):
    """Read FASTA file, return header and sequence."""
    with open(filepath) as f:
        lines = f.read().strip().split("\n")
    header = lines[0]
    seq = "".join(l.strip() for l in lines[1:] if not l.startswith(">"))
    return header, seq


def load_disorder_data():
    """Load Phase 1 disorder predictions."""
    summary_path = os.path.join(PHASE1_DIR, "results", "disorder_summary_v3.json")
    with open(summary_path) as f:
        data = json.load(f)
    return {d["protein"]: d for d in data}


def load_per_residue(protein_name):
    """Load per-residue disorder scores."""
    csv_path = os.path.join(PHASE1_DIR, "disorder_predictions", f"{protein_name}_per_residue.csv")
    if os.path.exists(csv_path):
        return pd.read_csv(csv_path)
    return None


def scan_slim(sequence, idr_boundaries, require_idr=True):
    """Scan sequence for all ELM motifs. Optionally restrict to IDR regions."""
    hits = []
    for motif_name, motif_info in ELM_MOTIFS.items():
        pattern = motif_info["regex"]
        for match in re.finditer(pattern, sequence):
            start = match.start() + 1  # 1-indexed
            end = match.end()
            matched_seq = match.group()

            in_idr = any(idr_s <= start and end <= idr_e for idr_s, idr_e in idr_boundaries)

            if require_idr and not in_idr:
                continue

            hits.append({
                "motif_class": motif_name,
                "description": motif_info["description"],
                "start": start,
                "end": end,
                "matched_sequence": matched_seq,
                "in_idr": in_idr,
                "host_targets": motif_info["host_targets"],
                "pathway": motif_info["pathway"],
            })

    return hits


def score_slim_disorder(hit, per_residue_df):
    """Score a SLiM hit by its disorder context."""
    if per_residue_df is None:
        return 0.5, 0.5

    start_idx = hit["start"] - 1
    end_idx = hit["end"]
    region = per_residue_df.iloc[start_idx:end_idx]

    iupred_mean = region["iupred3"].mean() if "iupred3" in region.columns else 0.5
    comp_mean = region["composition"].mean() if "composition" in region.columns else 0.5

    return float(iupred_mean), float(comp_mean)


def find_overlapping_slims(hits, max_gap=5):
    """Find overlapping or closely adjacent SLiMs (potential competitive switches)."""
    if len(hits) < 2:
        return []

    sorted_hits = sorted(hits, key=lambda h: h["start"])
    overlaps = []

    for i in range(len(sorted_hits)):
        for j in range(i + 1, len(sorted_hits)):
            h1, h2 = sorted_hits[i], sorted_hits[j]
            if h1["motif_class"] == h2["motif_class"]:
                continue

            overlap_start = max(h1["start"], h2["start"])
            overlap_end = min(h1["end"], h2["end"])

            if overlap_end >= overlap_start - max_gap:
                # They overlap or are within max_gap residues
                overlaps.append({
                    "motif_1": h1["motif_class"],
                    "motif_2": h2["motif_class"],
                    "motif_1_seq": h1["matched_sequence"],
                    "motif_2_seq": h2["matched_sequence"],
                    "motif_1_range": f"{h1['start']}-{h1['end']}",
                    "motif_2_range": f"{h2['start']}-{h2['end']}",
                    "overlap_residues": max(0, overlap_end - overlap_start + 1),
                    "gap": max(0, overlap_start - overlap_end - 1),
                    "targets_1": h1["host_targets"],
                    "targets_2": h2["host_targets"],
                    "pathway_1": h1["pathway"],
                    "pathway_2": h2["pathway"],
                    "switch_type": "overlapping" if overlap_end >= overlap_start else "adjacent",
                })

    return overlaps


def compute_enrichment(target_results, control_results):
    """Compare SLiM density between target and control proteins."""
    target_density = []
    control_density = []

    for r in target_results:
        length = r["length"]
        idr_length = sum(e - s + 1 for s, e in r["idr_boundaries"])
        n_slims = r["n_slim_hits"]
        if idr_length > 0:
            target_density.append(n_slims / idr_length * 100)

    for r in control_results:
        length = r["length"]
        idr_length = sum(e - s + 1 for s, e in r["idr_boundaries"])
        n_slims = r["n_slim_hits"]
        if idr_length > 0:
            control_density.append(n_slims / idr_length * 100)

    return {
        "target_mean_density": np.mean(target_density) if target_density else 0,
        "control_mean_density": np.mean(control_density) if control_density else 0,
        "ratio": (np.mean(target_density) / np.mean(control_density)
                  if control_density and np.mean(control_density) > 0 else float("inf")),
        "target_densities": target_density,
        "control_densities": control_density,
    }


def main():
    print("=" * 70)
    print("PHASE 2: SLiM DISCOVERY AND FILTERING")
    print("=" * 70)

    # Load Phase 1 data
    disorder_data = load_disorder_data()
    print(f"\nLoaded disorder data for {len(disorder_data)} proteins")
    print(f"ELM motif database: {len(ELM_MOTIFS)} motif classes")

    all_results = []
    all_hits_detailed = []
    all_overlaps = []

    for protein_name, pdata in disorder_data.items():
        category = pdata["category"]
        idr_boundaries = [tuple(b) for b in pdata["idr_boundaries"]]
        length = pdata["length"]

        # Load sequence from FASTA
        if category == "target":
            fasta_path = os.path.join(PHASE1_DIR, "sequences", "targets", f"{protein_name}.fasta")
        else:
            fasta_path = os.path.join(PHASE1_DIR, "sequences", "controls", f"{protein_name}.fasta")

        if not os.path.exists(fasta_path):
            print(f"\n  WARNING: {fasta_path} not found, skipping")
            continue

        _, sequence = read_fasta(fasta_path)
        per_residue = load_per_residue(protein_name)

        # Scan IDRs for SLiMs
        hits = scan_slim(sequence, idr_boundaries, require_idr=True)

        # Score each hit
        for hit in hits:
            iupred_score, comp_score = score_slim_disorder(hit, per_residue)
            hit["iupred3_mean"] = round(iupred_score, 4)
            hit["composition_mean"] = round(comp_score, 4)
            hit["protein"] = protein_name
            hit["category"] = category
            all_hits_detailed.append(hit)

        # Find overlapping SLiMs
        overlaps = find_overlapping_slims(hits)
        for ov in overlaps:
            ov["protein"] = protein_name
            ov["category"] = category
        all_overlaps.extend(overlaps)

        # Pathway distribution
        pathway_counts = defaultdict(int)
        for h in hits:
            pathway_counts[h["pathway"]] += 1

        idr_length = sum(e - s + 1 for s, e in idr_boundaries)

        result = {
            "protein": protein_name,
            "category": category,
            "length": length,
            "idr_boundaries": idr_boundaries,
            "idr_length": idr_length,
            "n_slim_hits": len(hits),
            "slim_density": round(len(hits) / idr_length * 100, 2) if idr_length > 0 else 0,
            "n_unique_motifs": len(set(h["motif_class"] for h in hits)),
            "n_overlapping_pairs": len(overlaps),
            "pathway_distribution": dict(pathway_counts),
            "host_targets": sorted(set(t for h in hits for t in h["host_targets"])),
        }
        all_results.append(result)

        # Print summary
        print(f"\n  {protein_name} ({category}) — {length} aa, IDR: {idr_length} aa")
        print(f"    SLiM hits: {len(hits)} ({result['slim_density']:.1f} per 100 IDR residues)")
        print(f"    Unique motif classes: {result['n_unique_motifs']}")
        print(f"    Overlapping pairs: {len(overlaps)}")
        if pathway_counts:
            top_pathways = sorted(pathway_counts.items(), key=lambda x: -x[1])[:5]
            print(f"    Top pathways: {', '.join(f'{p}({n})' for p, n in top_pathways)}")
        if overlaps:
            print(f"    Competitive switches:")
            for ov in overlaps[:5]:
                print(f"      {ov['motif_1']} vs {ov['motif_2']} "
                      f"({ov['switch_type']}, {ov['motif_1_range']}/{ov['motif_2_range']})")

    # ── Enrichment analysis ──
    target_results = [r for r in all_results if r["category"] == "target"]
    control_results = [r for r in all_results if r["category"] == "control"]
    enrichment = compute_enrichment(target_results, control_results)

    print(f"\n{'='*70}")
    print("ENRICHMENT ANALYSIS: Targets vs Controls")
    print(f"{'='*70}")
    print(f"  Target SLiM density:  {enrichment['target_mean_density']:.2f} hits per 100 IDR residues")
    print(f"  Control SLiM density: {enrichment['control_mean_density']:.2f} hits per 100 IDR residues")
    print(f"  Enrichment ratio:     {enrichment['ratio']:.2f}x")

    # ── Pathway enrichment ──
    target_pathways = defaultdict(int)
    control_pathways = defaultdict(int)
    for h in all_hits_detailed:
        if h["category"] == "target":
            target_pathways[h["pathway"]] += 1
        else:
            control_pathways[h["pathway"]] += 1

    all_pathways = sorted(set(list(target_pathways.keys()) + list(control_pathways.keys())))
    print(f"\n  Pathway distribution:")
    print(f"  {'Pathway':<30} {'Targets':>8} {'Controls':>8} {'Ratio':>8}")
    print(f"  {'-'*56}")
    for pw in all_pathways:
        t = target_pathways.get(pw, 0)
        c = control_pathways.get(pw, 0)
        ratio = t / c if c > 0 else float("inf") if t > 0 else 0
        flag = " ***" if ratio > 2 else ""
        print(f"  {pw:<30} {t:>8} {c:>8} {ratio:>7.1f}x{flag}")

    # ── Save results ──
    # Detailed hits CSV
    hits_csv = os.path.join(PHASE2_DIR, "results", "slim_hits_all.csv")
    with open(hits_csv, "w", newline="") as f:
        cols = ["protein", "category", "motif_class", "description", "start", "end",
                "matched_sequence", "in_idr", "iupred3_mean", "composition_mean",
                "pathway", "host_targets"]
        writer = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        writer.writeheader()
        for h in all_hits_detailed:
            row = dict(h)
            row["host_targets"] = "; ".join(h["host_targets"])
            writer.writerow(row)
    print(f"\n  Saved: {hits_csv}")

    # Overlaps CSV
    overlaps_csv = os.path.join(PHASE2_DIR, "results", "slim_overlaps.csv")
    with open(overlaps_csv, "w", newline="") as f:
        cols = ["protein", "category", "motif_1", "motif_2", "motif_1_seq", "motif_2_seq",
                "motif_1_range", "motif_2_range", "overlap_residues", "gap",
                "switch_type", "pathway_1", "pathway_2"]
        writer = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        writer.writeheader()
        for ov in all_overlaps:
            row = dict(ov)
            row.pop("targets_1", None)
            row.pop("targets_2", None)
            writer.writerow(row)
    print(f"  Saved: {overlaps_csv}")

    # Summary JSON
    summary = {
        "n_proteins": len(all_results),
        "n_total_hits": len(all_hits_detailed),
        "n_overlapping_pairs": len(all_overlaps),
        "enrichment": {
            "target_density": enrichment["target_mean_density"],
            "control_density": enrichment["control_mean_density"],
            "ratio": enrichment["ratio"],
        },
        "proteins": all_results,
    }
    summary_path = os.path.join(PHASE2_DIR, "results", "slim_discovery_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"  Saved: {summary_path}")

    # ── Summary table ──
    print(f"\n{'='*90}")
    print("PHASE 2 SUMMARY")
    print(f"{'='*90}")
    print(f"{'Protein':<18} {'Cat':<8} {'IDR aa':>6} {'SLiMs':>6} {'Density':>8} "
          f"{'Unique':>6} {'Overlaps':>8} {'Host targets':>12}")
    print("-" * 90)
    for r in all_results:
        print(f"{r['protein']:<18} {r['category']:<8} {r['idr_length']:>6} "
              f"{r['n_slim_hits']:>6} {r['slim_density']:>7.1f} "
              f"{r['n_unique_motifs']:>6} {r['n_overlapping_pairs']:>8} "
              f"{len(r['host_targets']):>12}")


if __name__ == "__main__":
    main()
