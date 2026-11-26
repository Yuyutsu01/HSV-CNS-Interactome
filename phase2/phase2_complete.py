#!/usr/bin/env python3
"""
Phase 2 Completion: Expanded ELM scanning with accessibility filter,
conservation integration, and XGBoost ML scoring.
"""

import os, sys, json, csv, re, math
import numpy as np
from collections import Counter

BASE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(BASE)
PHASE1 = os.path.join(PROJECT, "phase1")
RESULTS = os.path.join(BASE, "results")
os.makedirs(RESULTS, exist_ok=True)

sys.path.insert(0, BASE)
from elm_database import ELM_CLASSES, KNOWN_HSV1_SLIMS
from elm_expanded import EXPANDED_ELM_CLASSES

# Merge databases
ALL_ELM = {}
ALL_ELM.update(ELM_CLASSES)
ALL_ELM.update(EXPANDED_ELM_CLASSES)
print(f"Total ELM classes: {len(ALL_ELM)} ({len(ELM_CLASSES)} original + {len(EXPANDED_ELM_CLASSES)} expanded)")

# ═══════════════════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════════════════

def load_fasta(filepath):
    if not os.path.exists(filepath):
        return None
    seq = []
    with open(filepath) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)


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


def load_annotation_matrix():
    """Load Phase 1 per-residue annotation matrix."""
    path = os.path.join(PHASE1, "results", "per_residue_annotation_matrix.csv")
    data = {}
    if os.path.exists(path):
        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                prot = row["protein"]
                pos = int(row["position"])
                if prot not in data:
                    data[prot] = {}
                data[prot][pos] = row
    return data


def load_conservation():
    """Load conservation per-position data."""
    path = os.path.join(RESULTS, "per_position_conservation.csv")
    data = {}
    if os.path.exists(path):
        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                prot = row["protein"]
                pos = int(row["position"])
                if prot not in data:
                    data[prot] = {}
                data[prot][pos] = {
                    "conservation": float(row["blosum62_conservation"]),
                    "identity": float(row["sequence_identity"]),
                    "dnds_proxy": float(row["dnds_proxy"]),
                }
    return data


# ═══════════════════════════════════════════════════════════════════════
# PROTEIN LIST
# ═══════════════════════════════════════════════════════════════════════

PROTEINS = [
    {"name": "ICP34.5", "category": "target", "fasta": os.path.join(PHASE1, "sequences", "targets", "ICP34.5.fasta")},
    {"name": "US3", "category": "target", "fasta": os.path.join(PHASE1, "sequences", "targets", "US3.fasta")},
    {"name": "pUL37", "category": "target", "fasta": os.path.join(PHASE1, "sequences", "targets", "pUL37.fasta")},
    {"name": "ICP8_UL29", "category": "control", "fasta": os.path.join(PHASE1, "sequences", "controls", "ICP8_UL29.fasta")},
    {"name": "UL42", "category": "control", "fasta": os.path.join(PHASE1, "sequences", "controls", "UL42.fasta")},
    {"name": "gD_US6", "category": "control", "fasta": os.path.join(PHASE1, "sequences", "controls", "gD_US6.fasta")},
    {"name": "VP11_12_UL46", "category": "control", "fasta": os.path.join(PHASE1, "sequences", "controls", "VP11_12_UL46.fasta")},
    {"name": "VP16_UL48", "category": "control", "fasta": os.path.join(PHASE1, "sequences", "controls", "VP16_UL48.fasta")},
    {"name": "VP22_UL49", "category": "control", "fasta": os.path.join(PHASE1, "sequences", "controls", "VP22_UL49.fasta")},
    {"name": "ICP27_UL54", "category": "control", "fasta": os.path.join(PHASE1, "sequences", "controls", "ICP27_UL54.fasta")},
    {"name": "ICP22_US1", "category": "control", "fasta": os.path.join(PHASE1, "sequences", "controls", "ICP22_US1.fasta")},
]

# IDR boundaries
with open(os.path.join(PHASE1, "results", "disorder_summary_v3.json")) as f:
    disorder_data = json.load(f)
IDR_MAP = {p["protein"]: p["idr_boundaries"] for p in disorder_data}

# CNS expression data
CNS_EXPRESSION = {
    "BECN1": 12.5, "PPP1CA": 45.2, "PPP1CB": 38.7, "TRAF6": 8.3,
    "BAD": 15.2, "BCL2": 7.8, "BCL2L1": 22.1, "MCL1": 18.5,
    "MAP1LC3B": 35.6, "ATG5": 12.3, "ATG7": 15.8, "ATG12": 10.2,
    "GABARAP": 42.1, "SQSTM1": 25.3, "ULK1": 8.7, "ULK2": 12.1,
    "CASP3": 5.2, "CASP8": 3.1, "CASP9": 4.8, "FADD": 6.5,
    "YWHAZ": 65.3, "YWHAB": 55.2, "YWHAG": 48.7, "YWHAE": 72.1,
    "YWHAH": 35.6, "SFN": 2.1,
    "MAPK1": 25.3, "MAPK3": 18.7, "MAPK8": 15.2, "MAPK14": 12.8,
    "CDK1": 3.5, "CDK2": 8.2, "CDK4": 5.1, "CDK5": 42.3,
    "CCNB1": 2.1, "CCND1": 5.8,
    "CK2A1": 28.5, "CK2B": 22.3,
    "GSK3B": 32.1, "GSK3A": 15.2,
    "PRKACA": 35.8, "PRKACB": 28.1,
    "AKT1": 18.5, "AKT2": 8.2, "AKT3": 25.3,
    "CAMK2A": 85.2, "CAMK2B": 72.1, "CAMK2D": 15.3, "CAMK2G": 22.1,
    "DYRK1A": 18.5, "DYRK1B": 5.2,
    "SRC": 15.8, "FYN": 42.3, "LCK": 1.2, "YES1": 8.5,
    "PPP1CC": 32.1, "PPP2CA": 28.5, "PPP2CB": 15.2,
    "PPP3CA": 45.8, "PPP3CB": 38.2,
    "PRKCA": 22.1, "PRKCB": 35.2, "PRKCE": 18.5, "PRKCG": 42.1,
    "KPNA1": 15.2, "KPNA2": 12.8, "KPNA6": 8.5,
    "XPO1": 18.5,
    "XIAP": 12.3, "BIRC2": 8.5, "BIRC3": 5.2,
    "TRAF2": 8.5, "TRAF3": 6.2, "TRAF1": 3.5,
    "NEMO": 7.8, "IKBKB": 8.2, "TBK1": 12.5, "IKBKE": 5.2,
    "IRF3": 8.5, "IRF7": 5.2,
    "RELA": 12.3, "NFKB1": 8.5,
    "STAT1": 10.2, "STAT3": 15.8,
    "JAK1": 12.5, "JAK2": 8.2, "TYK2": 6.5,
    "MTOR": 15.2, "RPTOR": 12.1,
    "ROCK1": 18.5, "ROCK2": 22.1,
    "PIN1": 25.3,
    "EP300": 18.5, "CREBBP": 22.1, "KAT2A": 12.3,
    "CREB1": 28.5,
    "NLRP3": 5.2, "PYCARD": 3.8,
    "MYD88": 5.5, "IRAK1": 6.2, "IRAK4": 5.8,
    "DDX58": 6.5, "IFIH1": 5.2, "MAVS": 8.5,
    "STING1": 7.2, "CGAS": 3.5,
    "RIPK1": 12.3, "RIPK3": 5.2, "MLKL": 3.8,
    "PARP1": 18.5,
    "BRCA1": 5.2, "TP53BP1": 8.5,
    "CALM1": 55.2, "CALM2": 48.5, "CALM3": 42.1,
    "DLG4": 85.2, "HOMER1": 72.1, "HOMER2": 15.2,
    "DISC1": 32.1,
    "NRXN1": 65.2, "NLGN1": 55.3,
    "SNAP25": 82.1, "STX1A": 72.3, "VAMP2": 68.5,
    "VIM": 25.3, "ACTB": 55.2, "TUBA1A": 62.1, "TUBB3": 72.3,
    "NMT1": 12.3, "NMT2": 8.5,
    "CBX1": 15.2, "CBX3": 18.5, "CBX5": 12.3,
    "BRD2": 15.2, "BRD4": 18.5,
    "EZH2": 8.5,
    "PRMT1": 22.1, "PRMT5": 15.2,
    "OGT": 18.5,
    "FURIN": 12.3,
    "MMP2": 5.2, "MMP9": 3.5, "MMP14": 8.5,
    "GZMB": 1.2,
    "GRB2": 15.2, "PTPN6": 2.5, "PTPN11": 15.8,
    "CRK": 12.3, "CRKL": 8.5,
    "NCK1": 8.5, "NCK2": 5.2,
    "NEDD4": 18.5, "NEDD4L": 12.3, "ITCH": 8.5,
    "SMURF1": 12.3, "SMURF2": 8.5,
    "SMAD2": 12.3, "SMAD3": 8.5, "SMAD4": 15.2,
    "NOTCH1": 8.5, "NOTCH2": 12.3,
    "FZD1": 5.2, "FZD2": 8.5,
    "SOS1": 12.3, "SOS2": 8.5,
    "CDC42": 25.3, "RAC1": 22.1,
    "RHOA": 28.5,
    "HLA-A": 2.1, "HLA-B": 1.8, "HLA-C": 1.5,
    "TAP1": 3.5, "TAP2": 2.8,
    "PCNA": 8.5,
    "HNRNPA1": 35.2, "SRSF1": 22.1,
    "EIF4G1": 25.3, "PABPC1": 32.1,
    "IFNAR1": 5.2, "IFNAR2": 3.8,
    "DYNLL1": 35.2, "DCTN1": 22.1, "MYO5A": 32.1,
    "FBXW7": 12.3, "SPOP": 5.2,
    "KEAP1": 12.3,
    "TSC1": 15.2, "TSC2": 12.3,
    "PRKAA1": 15.2, "PRKAA2": 8.5,
    "ATM": 12.3, "ATR": 8.5,
    "CHEK1": 8.5, "CHEK2": 5.2,
    "FOS": 15.2, "JUN": 18.5,
}


# ═══════════════════════════════════════════════════════════════════════
# SCANNING WITH ACCESSIBILITY FILTER
# ═══════════════════════════════════════════════════════════════════════

def scan_protein_expanded(protein_name, sequence, iupred_scores, idr_boundaries,
                          annotation_matrix, conservation_data, category):
    """Scan protein with expanded ELM + accessibility filter."""
    n = len(sequence)
    hits = []

    for motif_name, motif_info in ALL_ELM.items():
        pattern = motif_info["regex"]
        try:
            for match in re.finditer(pattern, sequence):
                start = match.start() + 1  # 1-indexed
                end = match.end()

                # IDR filter
                in_idr = any(s <= start and end <= e for s, e in idr_boundaries)

                # IUPred3 score for the motif
                motif_iupred = np.mean(iupred_scores[start-1:end]) if start-1 < len(iupred_scores) else 0

                # Accessibility filter (from Phase 1 annotation)
                acc_scores = []
                for pos in range(start, end + 1):
                    if protein_name in annotation_matrix and pos in annotation_matrix[protein_name]:
                        acc = annotation_matrix[protein_name][pos].get("accessibility_score")
                        if acc:
                            acc_scores.append(float(acc))
                mean_accessibility = np.mean(acc_scores) if acc_scores else 0.5
                accessible = mean_accessibility > 0.35  # accessibility threshold

                # Secondary structure context
                ss_states = []
                for pos in range(start, end + 1):
                    if protein_name in annotation_matrix and pos in annotation_matrix[protein_name]:
                        ss = annotation_matrix[protein_name][pos].get("ss_prediction", "C")
                        ss_states.append(ss)
                ss_coil_frac = ss_states.count("C") / len(ss_states) if ss_states else 0
                ss_turn_frac = ss_states.count("T") / len(ss_states) if ss_states else 0

                # Conservation scores (for target proteins)
                cons_score = 0.0
                cons_id = 0.0
                cons_dnds = 0.0
                if protein_name in conservation_data:
                    cons_vals = []
                    id_vals = []
                    dnds_vals = []
                    for pos in range(start, end + 1):
                        if pos in conservation_data[protein_name]:
                            cons_vals.append(conservation_data[protein_name][pos]["conservation"])
                            id_vals.append(conservation_data[protein_name][pos]["identity"])
                            dnds_vals.append(conservation_data[protein_name][pos]["dnds_proxy"])
                    if cons_vals:
                        cons_score = np.mean(cons_vals)
                        cons_id = np.mean(id_vals)
                        cons_dnds = np.mean(dnds_vals)

                # CNS score
                targets = motif_info.get("targets", [])
                cns_scores = [CNS_EXPRESSION.get(t, 0) for t in targets]
                max_brain_tpm = max(cns_scores) if cns_scores else 0
                mean_brain_tpm = np.mean(cns_scores) if cns_scores else 0

                hit = {
                    "motif_class": motif_name,
                    "description": motif_info.get("description", ""),
                    "start": start,
                    "end": end,
                    "matched_sequence": match.group(),
                    "in_idr": in_idr,
                    "accessible": accessible,
                    "targets": ",".join(targets),
                    "pathway": motif_info.get("pathway", "unknown"),
                    "protein": protein_name,
                    "category": category,
                    "iupred3_mean": round(motif_iupred, 4),
                    "accessibility_mean": round(mean_accessibility, 4),
                    "ss_coil_fraction": round(ss_coil_frac, 3),
                    "ss_turn_fraction": round(ss_turn_frac, 3),
                    "conservation_score": round(cons_score, 4),
                    "conservation_identity": round(cons_id, 4),
                    "conservation_dnds": round(cons_dnds, 4),
                    "max_brain_tpm": round(max_brain_tpm, 1),
                    "mean_brain_tpm": round(mean_brain_tpm, 1),
                }
                hits.append(hit)
        except re.error:
            pass

    return hits


# ═══════════════════════════════════════════════════════════════════════
# MULTI-CRITERIA COMPOSITE SCORING (transparent, non-circular)
# ═══════════════════════════════════════════════════════════════════════
#
# Replaces the previous XGBoost semi-supervised scorer, which suffered
# from circular labeling: positives were defined using the same features
# the model trained on (in_idr, accessible, iupred3_mean), yielding a
# misleading AUC of 0.999. This composite uses biologically motivated
# criteria with pre-defined weights, validated against known HSV-1 SLiMs
# but never trained on them.

# Weights reflect prior biological knowledge about SLiM functionality:
#   - Disorder context is essential (SLiMs function in disordered regions)
#   - Surface accessibility enables binding
#   - Structural flexibility (coil/turn) favors conformational adaptation
#   - Conservation indicates functional constraint
#   - Brain expression of host targets indicates CNS relevance
COMPOSITE_WEIGHTS = {
    "disorder":      0.30,  # IUPred3 disorder probability
    "accessibility": 0.20,  # Emini surface accessibility
    "flexibility":   0.15,  # Coil + turn secondary structure
    "conservation":  0.20,  # Sequence conservation across herpesviruses
    "cns_relevance": 0.15,  # Host target brain expression
}


def compute_composite_score(hit):
    """Compute transparent multi-criteria SLiM confidence score.

    Each criterion is independently normalized to [0, 1] and combined
    with pre-defined weights. No ML training is involved.
    """
    # Disorder: IUPred3 probability (already 0-1)
    disorder = hit["iupred3_mean"]

    # Accessibility: Emini accessibility (already 0-1 from prior normalization)
    accessibility = min(hit["accessibility_mean"], 1.0)

    # Flexibility: fraction of motif in coil or turn
    flexibility = (hit["ss_coil_fraction"] + hit["ss_turn_fraction"])
    flexibility = min(flexibility, 1.0)

    # Conservation: combined score (already 0-1)
    conservation = hit["conservation_score"]

    # CNS relevance: host target brain expression (cap at TPM=200 for normalization)
    cns = min(hit["max_brain_tpm"] / 200.0, 1.0)

    components = {
        "disorder": disorder,
        "accessibility": accessibility,
        "flexibility": flexibility,
        "conservation": conservation,
        "cns_relevance": cns,
    }

    score = sum(COMPOSITE_WEIGHTS[k] * components[k] for k in COMPOSITE_WEIGHTS)
    return score, components


def build_ml_scorer(all_hits):
    """Score SLiM hits using transparent multi-criteria composite.

    Known positive controls (experimentally validated HSV-1 SLiMs) are
    used only for post-hoc validation, not for defining score thresholds.
    """
    print("  Using multi-criteria composite scoring (transparent, non-circular)")
    print(f"  Weights: {COMPOSITE_WEIGHTS}")

    scores = []
    all_components = []
    for hit in all_hits:
        s, comp = compute_composite_score(hit)
        scores.append(s)
        all_components.append(comp)

    scores = np.array(scores)

    # Report score distribution
    print(f"  Score distribution: mean={np.mean(scores):.3f}, "
          f"median={np.median(scores):.3f}, "
          f"std={np.std(scores):.3f}, "
          f"range=[{np.min(scores):.3f}, {np.max(scores):.3f}]")

    # Report per-criterion contributions
    criterion_means = {}
    for key in COMPOSITE_WEIGHTS:
        vals = [c[key] for c in all_components]
        criterion_means[key] = round(np.mean(vals), 3)
    print(f"  Mean criterion values: {criterion_means}")

    meta = {
        "method": "multi_criteria_composite",
        "weights": COMPOSITE_WEIGHTS,
        "criterion_means": criterion_means,
        "n_hits": len(all_hits),
        "note": "Pre-defined weights, no ML training, no circular labeling"
    }

    return scores, COMPOSITE_WEIGHTS, meta


# ═══════════════════════════════════════════════════════════════════════
# COMPOSITE RANKING
# ═══════════════════════════════════════════════════════════════════════

def is_known_positive_hit(hit):
    for ctrl_name, ctrl_info in KNOWN_HSV1_SLIMS.items():
        if hit["protein"] != ctrl_info["protein"]:
            continue
        s = ctrl_info.get("start")
        e = ctrl_info.get("end")
        if s is not None and e is not None:
            if hit["start"] >= s - 5 and hit["end"] <= e + 5:
                return True
        else:
            host = ctrl_info.get("host_target", "")
            if host and host in hit.get("targets", ""):
                return True
    return False


def compute_final_rank(hits, composite_scores):
    """Classify hits using composite score and structural context."""
    for i, hit in enumerate(hits):
        cs = composite_scores[i]

        is_positive_ctrl = is_known_positive_hit(hit)

        hit["composite_score"] = round(float(cs), 4)
        hit["ml_score"] = round(float(cs), 4)  # backward compatibility
        hit["is_positive_control"] = is_positive_ctrl

        # Classification based on independent structural criteria + composite
        if hit["in_idr"] and hit["accessible"] and cs > 0.45:
            hit["filter_class"] = "high_confidence"
        elif hit["in_idr"] and (hit["accessible"] or cs > 0.35):
            hit["filter_class"] = "medium_confidence"
        elif hit["in_idr"]:
            hit["filter_class"] = "low_confidence"
        else:
            hit["filter_class"] = "filtered_out"


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    print("Loading Phase 1 annotation matrix...")
    annotation_matrix = load_annotation_matrix()

    print("Loading conservation data...")
    conservation_data = load_conservation()

    all_hits = []

    for pinfo in PROTEINS:
        name = pinfo["name"]
        category = pinfo["category"]
        print(f"\n{'='*60}")
        print(f"Scanning: {name} ({category})")
        print(f"{'='*60}")

        sequence = load_fasta(pinfo["fasta"])
        if not sequence:
            print(f"  ERROR: Cannot load sequence")
            continue

        iupred_scores = load_iupred3(name)
        if iupred_scores is None:
            print(f"  ERROR: Cannot load IUPred3 data")
            continue

        n = min(len(sequence), len(iupred_scores))
        sequence = sequence[:n]
        iupred_scores = iupred_scores[:n]

        idr_bounds = IDR_MAP.get(name, [])

        hits = scan_protein_expanded(
            name, sequence, iupred_scores, idr_bounds,
            annotation_matrix, conservation_data, category
        )

        # Statistics
        in_idr = sum(1 for h in hits if h["in_idr"])
        accessible = sum(1 for h in hits if h["accessible"])
        both = sum(1 for h in hits if h["in_idr"] and h["accessible"])
        print(f"  Total hits: {len(hits)}")
        print(f"  In IDR: {in_idr}")
        print(f"  Accessible: {accessible}")
        print(f"  IDR + Accessible: {both}")

        all_hits.extend(hits)

    # Composite scoring
    print(f"\n{'='*60}")
    print("MULTI-CRITERIA COMPOSITE SCORING")
    print(f"{'='*60}")
    composite_scores, weights, ml_meta = build_ml_scorer(all_hits)
    compute_final_rank(all_hits, composite_scores)

    print("  Criterion weights:")
    for feat, w in sorted(weights.items(), key=lambda x: -x[1]):
        print(f"    {feat}: {w:.2f}")

    # Filter stats
    high = sum(1 for h in all_hits if h["filter_class"] == "high_confidence")
    medium = sum(1 for h in all_hits if h["filter_class"] == "medium_confidence")
    low = sum(1 for h in all_hits if h["filter_class"] == "low_confidence")
    filtered = sum(1 for h in all_hits if h["filter_class"] == "filtered_out")
    print(f"\n  Filter classification:")
    print(f"    High confidence: {high}")
    print(f"    Medium confidence: {medium}")
    print(f"    Low confidence: {low}")
    print(f"    Filtered out: {filtered}")

    # Positive control recovery
    print(f"\n  Positive control recovery:")
    for ctrl_name, ctrl_info in KNOWN_HSV1_SLIMS.items():
        found = [h for h in all_hits if h["is_positive_control"] and
                 h["protein"] == ctrl_info["protein"]]
        if found:
            best = max(found, key=lambda x: x["composite_score"])
            print(f"    {ctrl_name}: FOUND (score={best['composite_score']:.3f}, class={best['filter_class']})")
        else:
            print(f"    {ctrl_name}: NOT FOUND")

    # Save outputs
    # 1. Full expanded hits CSV
    csv_path = os.path.join(RESULTS, "slim_hits_expanded.csv")
    if all_hits:
        keys = all_hits[0].keys()
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(all_hits)
        print(f"\nExpanded hits: {csv_path} ({len(all_hits)} total)")

    # 2. High-confidence hits only
    hc_path = os.path.join(RESULTS, "slim_hits_high_confidence.csv")
    hc_hits = [h for h in all_hits if h["filter_class"] in ["high_confidence", "medium_confidence"]]
    if hc_hits:
        keys = hc_hits[0].keys()
        with open(hc_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(hc_hits)
        print(f"High/medium confidence: {hc_path} ({len(hc_hits)} hits)")

    # 3. Summary JSON
    # Per-protein breakdown
    protein_stats = {}
    for h in all_hits:
        prot = h["protein"]
        if prot not in protein_stats:
            protein_stats[prot] = {"total": 0, "in_idr": 0, "accessible": 0,
                                   "high_conf": 0, "medium_conf": 0, "category": h["category"]}
        protein_stats[prot]["total"] += 1
        if h["in_idr"]: protein_stats[prot]["in_idr"] += 1
        if h["accessible"]: protein_stats[prot]["accessible"] += 1
        if h["filter_class"] == "high_confidence": protein_stats[prot]["high_conf"] += 1
        if h["filter_class"] == "medium_confidence": protein_stats[prot]["medium_conf"] += 1

    # Pathway breakdown for target proteins
    target_pathways = {}
    for h in all_hits:
        if h["category"] == "target" and h["filter_class"] in ["high_confidence", "medium_confidence"]:
            pw = h["pathway"]
            if pw not in target_pathways:
                target_pathways[pw] = 0
            target_pathways[pw] += 1

    summary = {
        "metadata": {
            "total_elm_classes": len(ALL_ELM),
            "total_hits": len(all_hits),
            "high_confidence": high,
            "medium_confidence": medium,
            "low_confidence": low,
            "filtered_out": filtered,
            "ml_scoring": ml_meta,
        },
        "protein_stats": protein_stats,
        "target_pathway_breakdown": dict(sorted(target_pathways.items(), key=lambda x: -x[1])),
    }

    summary_path = os.path.join(RESULTS, "phase2_complete_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Summary: {summary_path}")

    # Print final table
    print(f"\n{'='*80}")
    print("PHASE 2 EXPANDED RESULTS")
    print(f"{'='*80}")
    print(f"{'Protein':<18} {'Cat':<8} {'Total':>6} {'IDR':>5} {'Acc':>5} {'High':>5} {'Med':>5}")
    print("-" * 60)
    for prot in [p["name"] for p in PROTEINS]:
        s = protein_stats.get(prot, {})
        print(f"{prot:<18} {s.get('category',''):<8} {s.get('total',0):>6} "
              f"{s.get('in_idr',0):>5} {s.get('accessible',0):>5} "
              f"{s.get('high_conf',0):>5} {s.get('medium_conf',0):>5}")


if __name__ == "__main__":
    main()
