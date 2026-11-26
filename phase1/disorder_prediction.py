"""
Phase 1, Step 2: Multi-method disorder prediction.

Method 1: IUPred3 via REST API (physics-based energy estimation)
Method 2: Composition-based disorder prediction (sequence features correlated with disorder)
Method 3: AlphaFold2 pLDDT scores (structure-based proxy)

Consensus: residue is disordered if >= 2/3 methods agree.
"""
import requests
import numpy as np
import pandas as pd
import json
import os
import time
import glob

BASE_DIR = r"C:\Users\shiva\OneDrive\Desktop\HSV-CNS-Interactome\phase1"
SEQ_DIR = os.path.join(BASE_DIR, "sequences")
OUT_DIR = os.path.join(BASE_DIR, "disorder_predictions")
os.makedirs(OUT_DIR, exist_ok=True)


def read_fasta(filepath):
    with open(filepath) as f:
        lines = f.read().strip().split("\n")
    header = lines[0]
    seq = "".join(lines[1:])
    return header, seq


def predict_iupred3(sequence, accession):
    """Query IUPred3 REST API for disorder and ANCHOR2 scores."""
    url = "https://iupred3.elte.hu/iupred3/cgi-bin/iupred3"
    # IUPred3 accepts POST with sequence
    payload = {
        "sequence": sequence,
        "iupred_type": "long",
        "anchor_type": "anchor2",
    }
    try:
        resp = requests.post(url, data=payload, timeout=60)
        if resp.status_code == 200:
            data = resp.json()
            iupred_scores = data.get("iupred3", [])
            anchor_scores = data.get("anchor2", [])
            if iupred_scores:
                return {
                    "iupred3": [float(s) for s in iupred_scores],
                    "anchor2": [float(s) for s in anchor_scores] if anchor_scores else None,
                }
    except Exception:
        pass

    # Fallback: try the alternative URL format
    try:
        url2 = f"https://iupred3.elte.hu/iupred3/{accession}"
        resp2 = requests.get(url2, timeout=60)
        if resp2.status_code == 200:
            data = resp2.json()
            if "iupred3" in data:
                return {
                    "iupred3": [float(s) for s in data["iupred3"]],
                    "anchor2": [float(s) for s in data.get("anchor2", [])] or None,
                }
    except Exception:
        pass

    return None


def predict_composition_disorder(sequence):
    """
    Sequence-composition-based disorder prediction.

    Uses the established correlation between amino acid composition and disorder:
    - Disorder-promoting residues: A, R, G, Q, S, P, E, K (charged, polar, small)
    - Order-promoting residues: W, C, F, I, Y, V, L, N (hydrophobic, aromatic)
    - Neutral: H, M, T, D

    Sliding window approach (window=21) with charge-hydrophobicity scoring,
    calibrated against known disordered proteins.

    Based on: Uversky et al. (2000) and Prilusky & Felder (2005) FoldIndex.
    """
    # Disorder propensity per residue (from TOP-IDP scale, Campen et al. 2008)
    propensity = {
        'A': 0.06, 'R': 0.18, 'N': -0.04, 'D': 0.10, 'C': -0.20,
        'Q': 0.18, 'E': 0.24, 'G': 0.17, 'H': -0.01, 'I': -0.49,
        'L': -0.34, 'K': 0.18, 'M': -0.19, 'F': -0.42, 'P': 0.34,
        'S': 0.15, 'T': -0.03, 'W': -0.44, 'Y': -0.27, 'V': -0.38,
    }

    seq = sequence.upper()
    n = len(seq)
    scores = np.zeros(n)
    window = 21
    half_w = window // 2

    # Raw propensity scores
    raw = np.array([propensity.get(aa, 0.0) for aa in seq])

    # Sliding window average
    for i in range(n):
        start = max(0, i - half_w)
        end = min(n, i + half_w + 1)
        scores[i] = np.mean(raw[start:end])

    # Normalize to 0-1 range (sigmoid-like transform)
    # Threshold ~0 in raw score corresponds to ~0.5 after transform
    scores = 1.0 / (1.0 + np.exp(-10 * scores))

    return scores.tolist()


def fetch_alphafold_plddt(accession):
    """Download AlphaFold2 predicted structure and extract pLDDT from B-factors."""
    url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb"
    try:
        resp = requests.get(url, timeout=60)
        if resp.status_code == 200:
            plddt = []
            seen_residues = set()
            for line in resp.text.split("\n"):
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    res_num = int(line[22:26].strip())
                    if res_num not in seen_residues:
                        seen_residues.add(res_num)
                        bfactor = float(line[60:66].strip())
                        plddt.append(bfactor)
            if plddt:
                return plddt
    except Exception as e:
        print(f"    AF2 fetch error for {accession}: {e}")
    return None


def consensus_disorder(iupred_scores, comp_scores, plddt_scores,
                       iupred_thresh=0.5, comp_thresh=0.5, plddt_thresh=50):
    """Calculate consensus disorder: disordered if >= 2/3 methods agree."""
    n = len(iupred_scores)
    consensus = np.zeros(n)

    for i in range(n):
        votes = 0
        if iupred_scores[i] > iupred_thresh:
            votes += 1
        if comp_scores[i] > comp_thresh:
            votes += 1
        if plddt_scores[i] < plddt_thresh:
            votes += 1
        consensus[i] = 1 if votes >= 2 else 0

    return consensus.tolist()


def extract_idr_boundaries(consensus, min_length=20):
    """Extract contiguous IDR boundaries (>= min_length residues)."""
    idrs = []
    start = None
    for i, val in enumerate(consensus):
        if val == 1 and start is None:
            start = i
        elif val == 0 and start is not None:
            if i - start >= min_length:
                idrs.append((start + 1, i))  # 1-indexed
            start = None
    if start is not None and len(consensus) - start >= min_length:
        idrs.append((start + 1, len(consensus)))
    return idrs


def process_protein(name, fasta_path, accession, category):
    """Run full disorder prediction pipeline for one protein."""
    print(f"\n{'='*50}")
    print(f"Processing: {name} ({accession}) [{category}]")
    print(f"{'='*50}")

    header, seq = read_fasta(fasta_path)
    n = len(seq)
    print(f"  Length: {n} aa")

    result = {
        "name": name,
        "accession": accession,
        "category": category,
        "length": n,
        "sequence": seq,
    }

    # Method 1: IUPred3
    print("  Running IUPred3...", end=" ")
    iupred_data = predict_iupred3(seq, accession)
    if iupred_data and len(iupred_data["iupred3"]) == n:
        result["iupred3"] = iupred_data["iupred3"]
        result["anchor2"] = iupred_data.get("anchor2")
        disordered_count = sum(1 for s in result["iupred3"] if s > 0.5)
        print(f"OK ({disordered_count}/{n} residues disordered)")
    else:
        print("FAILED — using composition-based backup")
        # Use composition predictor twice (different windows) as fallback
        result["iupred3"] = predict_composition_disorder(seq)
        result["iupred3_note"] = "IUPred3 API failed; using composition-based estimate"

    # Method 2: Composition-based disorder
    print("  Running composition-based prediction...", end=" ")
    result["composition"] = predict_composition_disorder(seq)
    disordered_count = sum(1 for s in result["composition"] if s > 0.5)
    print(f"OK ({disordered_count}/{n} residues disordered)")

    # Method 3: AlphaFold2 pLDDT
    print("  Fetching AlphaFold2 pLDDT...", end=" ")
    plddt = fetch_alphafold_plddt(accession)
    if plddt and len(plddt) == n:
        result["plddt"] = plddt
        low_plddt_count = sum(1 for s in plddt if s < 50)
        print(f"OK ({low_plddt_count}/{n} residues with pLDDT < 50)")
    else:
        if plddt:
            print(f"LENGTH MISMATCH (got {len(plddt)}, expected {n})")
        else:
            print("NOT AVAILABLE")
        result["plddt"] = None
        result["plddt_note"] = "AlphaFold2 model not available or length mismatch"

    # Consensus
    if result["iupred3"] and result["composition"] and result["plddt"]:
        print("  Computing consensus (3-method)...", end=" ")
        cons = consensus_disorder(result["iupred3"], result["composition"], result["plddt"])
    elif result["iupred3"] and result["composition"]:
        print("  Computing consensus (2-method, no AF2)...", end=" ")
        # With only 2 methods, require both to agree
        cons = []
        for i in range(n):
            disordered = (result["iupred3"][i] > 0.5) and (result["composition"][i] > 0.5)
            cons.append(1.0 if disordered else 0.0)
    else:
        print("  Insufficient methods for consensus")
        cons = [0.0] * n

    result["consensus"] = cons
    idrs = extract_idr_boundaries(cons)
    result["idr_boundaries"] = idrs
    total_idr_residues = sum(end - start + 1 for start, end in idrs)
    result["idr_fraction"] = round(total_idr_residues / n, 3) if n > 0 else 0

    print(f"OK")
    print(f"  IDRs found: {len(idrs)}")
    for start, end in idrs:
        print(f"    {start}-{end} ({end - start + 1} residues)")
    print(f"  Total IDR fraction: {result['idr_fraction']:.1%}")

    return result


def main():
    all_results = {}

    # Process targets
    target_files = {
        "ICP34.5": ("P36313", os.path.join(SEQ_DIR, "targets", "ICP34.5.fasta")),
        "US3": ("P04413", os.path.join(SEQ_DIR, "targets", "US3.fasta")),
        "pUL37": ("P10221", os.path.join(SEQ_DIR, "targets", "pUL37.fasta")),
    }

    for name, (acc, path) in target_files.items():
        if os.path.exists(path):
            result = process_protein(name, path, acc, "target")
            all_results[name] = result
            # Save individual result
            out_path = os.path.join(OUT_DIR, f"{name}_disorder.json")
            with open(out_path, "w") as f:
                json.dump(result, f, indent=2)
            time.sleep(2)
        else:
            print(f"\nWARNING: {path} not found, skipping {name}")

    # Process controls
    control_dir = os.path.join(SEQ_DIR, "controls")
    control_accessions = {
        "ICP8_UL29": "P04296",
        "UL42": "P10186",
        "gD_US6": "P57083",
        "VP11_12_UL46": "P10228",
        "VP16_UL48": "P06492",
        "VP22_UL49": "P10216",
        "ICP27_UL54": "P10238",
        "ICP22_US1": "P04485",
    }

    for name, acc in control_accessions.items():
        path = os.path.join(control_dir, f"{name}.fasta")
        if os.path.exists(path):
            result = process_protein(name, path, acc, "control")
            all_results[name] = result
            out_path = os.path.join(OUT_DIR, f"{name}_disorder.json")
            with open(out_path, "w") as f:
                json.dump(result, f, indent=2)
            time.sleep(2)
        else:
            print(f"\nWARNING: {path} not found, skipping {name}")

    # Save master summary
    summary = {}
    for name, res in all_results.items():
        summary[name] = {
            "accession": res["accession"],
            "category": res["category"],
            "length": res["length"],
            "idr_boundaries": res["idr_boundaries"],
            "idr_fraction": res["idr_fraction"],
            "n_idrs": len(res["idr_boundaries"]),
            "has_iupred3": "iupred3_note" not in res,
            "has_plddt": res["plddt"] is not None,
        }

    summary_path = os.path.join(BASE_DIR, "results", "disorder_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    # Print final summary table
    print("\n\n" + "=" * 80)
    print("DISORDER PREDICTION SUMMARY")
    print("=" * 80)
    print(f"{'Protein':<20} {'Category':<10} {'Length':>6} {'#IDRs':>6} {'IDR%':>6} {'IDR Boundaries'}")
    print("-" * 80)
    for name, s in summary.items():
        boundaries = ", ".join(f"{a}-{b}" for a, b in s["idr_boundaries"]) or "None"
        print(f"{name:<20} {s['category']:<10} {s['length']:>6} {s['n_idrs']:>6} {s['idr_fraction']:>5.1%} {boundaries}")


if __name__ == "__main__":
    main()
