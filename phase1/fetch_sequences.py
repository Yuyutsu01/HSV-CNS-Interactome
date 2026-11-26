"""
Phase 1, Step 1: Fetch target, ortholog, and control protein sequences from UniProt.
"""
import requests
import json
import os
import time

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# ─── Target proteins (HSV-1 strain 17) ───
TARGETS = {
    "ICP34.5": {"uniprot": "P36313", "gene": "RL1", "description": "Neurovirulence factor ICP34.5"},
    "US3":     {"uniprot": "P04413", "gene": "US3", "description": "Serine/threonine-protein kinase US3"},
    "pUL37":   {"uniprot": "P10218", "gene": "UL37", "description": "Tegument protein pUL37"},
}

# ─── Negative control proteins (HSV-1, NOT neurovirulence determinants, with IDRs) ───
CONTROLS = {
    "ICP8_UL29":   {"uniprot": "P04296", "gene": "UL29", "description": "ssDNA-binding protein ICP8"},
    "UL42":        {"uniprot": "P10186", "gene": "UL42", "description": "DNA polymerase processivity factor"},
    "gD_US6":      {"uniprot": "P57083", "gene": "US6", "description": "Envelope glycoprotein D"},
    "VP11_12_UL46":{"uniprot": "P10228", "gene": "UL46", "description": "Tegument protein VP11/12"},
    "VP16_UL48":   {"uniprot": "P06492", "gene": "UL48", "description": "Transactivator VP16"},
    "VP22_UL49":   {"uniprot": "P10216", "gene": "UL49", "description": "Tegument protein VP22"},
    "ICP27_UL54":  {"uniprot": "P10238", "gene": "UL54", "description": "mRNA export factor ICP27"},
    "ICP22_US1":   {"uniprot": "P04485", "gene": "US1", "description": "Regulatory protein ICP22"},
}

# ─── Alphaherpesvirus orthologs for evolutionary analysis ───
# These are approximate orthologs — exact accessions verified from UniProt
ORTHOLOGS = {
    "ICP34.5": {
        "HSV2_RL1":    {"uniprot": "P89449", "species": "HSV-2 (HG52)"},
    },
    "US3": {
        "HSV2_US3":    {"uniprot": "P06485", "species": "HSV-2 (HG52)"},
        "VZV_ORF66":   {"uniprot": "P09259", "species": "VZV (Dumas)"},
        "PRV_US3":     {"uniprot": "Q01139", "species": "PRV (Becker)"},
        "BoHV1_US3":   {"uniprot": "Q65563", "species": "BoHV-1"},
        "EHV1_US3":    {"uniprot": "P28969", "species": "EHV-1 (Ab4)"},
    },
    "pUL37": {
        "HSV2_UL37":   {"uniprot": "G9I242", "species": "HSV-2 (HG52)"},
        "VZV_ORF21":   {"uniprot": "P09252", "species": "VZV (Dumas)"},
        "PRV_UL37":    {"uniprot": "Q04532", "species": "PRV"},
        "EHV1_UL37":   {"uniprot": "P28952", "species": "EHV-1 (Ab4)"},
    },
}


def fetch_uniprot_fasta(accession, max_retries=3):
    """Fetch FASTA sequence from UniProt REST API."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200 and resp.text.startswith(">"):
                return resp.text.strip()
            elif resp.status_code == 300:
                # Merged/demerged entry — follow redirect
                lines = resp.text.strip().split("\n")
                for line in lines:
                    if line.startswith("UniProtKB:"):
                        new_acc = line.split(":")[1].strip()
                        return fetch_uniprot_fasta(new_acc, max_retries=1)
            print(f"  WARNING: {accession} returned status {resp.status_code}, attempt {attempt+1}")
        except requests.exceptions.RequestException as e:
            print(f"  ERROR fetching {accession}: {e}, attempt {attempt+1}")
        time.sleep(2)
    return None


def fetch_uniprot_metadata(accession):
    """Fetch protein metadata (length, features) from UniProt JSON API."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            data = resp.json()
            length = data.get("sequence", {}).get("length", "?")
            organism = data.get("organism", {}).get("scientificName", "?")

            # Extract disorder annotations if present
            disorder_regions = []
            for feature in data.get("features", []):
                if feature.get("type") == "Region" and "disordered" in feature.get("description", "").lower():
                    start = feature.get("location", {}).get("start", {}).get("value", "?")
                    end = feature.get("location", {}).get("end", {}).get("value", "?")
                    disorder_regions.append(f"{start}-{end}")

            return {
                "length": length,
                "organism": organism,
                "disorder_annotations": disorder_regions,
            }
    except Exception as e:
        print(f"  Metadata fetch error for {accession}: {e}")
    return None


def parse_fasta_sequence(fasta_text):
    """Extract raw sequence from FASTA text."""
    lines = fasta_text.strip().split("\n")
    header = lines[0]
    sequence = "".join(line.strip() for line in lines[1:])
    return header, sequence


def save_fasta(fasta_text, filepath):
    """Save FASTA text to file."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        f.write(fasta_text + "\n")


def main():
    results = {"targets": {}, "controls": {}, "orthologs": {}}

    # ── Fetch target proteins ──
    print("=" * 60)
    print("FETCHING TARGET PROTEINS")
    print("=" * 60)
    for name, info in TARGETS.items():
        print(f"\n--- {name} ({info['uniprot']}) ---")
        fasta = fetch_uniprot_fasta(info["uniprot"])
        meta = fetch_uniprot_metadata(info["uniprot"])

        if fasta:
            filepath = os.path.join(BASE_DIR, "sequences", "targets", f"{name}.fasta")
            save_fasta(fasta, filepath)
            header, seq = parse_fasta_sequence(fasta)
            print(f"  Length: {len(seq)} aa")
            print(f"  Header: {header[:80]}...")
            if meta and meta["disorder_annotations"]:
                print(f"  UniProt disorder annotations: {', '.join(meta['disorder_annotations'])}")
            else:
                print(f"  No UniProt disorder annotations found")
            results["targets"][name] = {
                "accession": info["uniprot"],
                "length": len(seq),
                "disorder_annotations": meta["disorder_annotations"] if meta else [],
            }
        else:
            print(f"  FAILED to fetch sequence")
            results["targets"][name] = {"accession": info["uniprot"], "error": "fetch failed"}
        time.sleep(1)

    # ── Fetch control proteins ──
    print("\n" + "=" * 60)
    print("FETCHING NEGATIVE CONTROL PROTEINS")
    print("=" * 60)
    for name, info in CONTROLS.items():
        print(f"\n--- {name} ({info['uniprot']}) ---")
        fasta = fetch_uniprot_fasta(info["uniprot"])
        meta = fetch_uniprot_metadata(info["uniprot"])

        if fasta:
            filepath = os.path.join(BASE_DIR, "sequences", "controls", f"{name}.fasta")
            save_fasta(fasta, filepath)
            header, seq = parse_fasta_sequence(fasta)
            print(f"  Length: {len(seq)} aa")
            if meta and meta["disorder_annotations"]:
                print(f"  UniProt disorder annotations: {', '.join(meta['disorder_annotations'])}")
            results["controls"][name] = {
                "accession": info["uniprot"],
                "length": len(seq),
                "disorder_annotations": meta["disorder_annotations"] if meta else [],
            }
        else:
            print(f"  FAILED to fetch sequence")
        time.sleep(1)

    # ── Fetch orthologs ──
    print("\n" + "=" * 60)
    print("FETCHING ALPHAHERPESVIRUS ORTHOLOGS")
    print("=" * 60)
    for target_name, ortho_dict in ORTHOLOGS.items():
        print(f"\n--- Orthologs of {target_name} ---")
        results["orthologs"][target_name] = {}
        for ortho_name, ortho_info in ortho_dict.items():
            print(f"  {ortho_name} ({ortho_info['uniprot']}, {ortho_info['species']})")
            fasta = fetch_uniprot_fasta(ortho_info["uniprot"])
            if fasta:
                filepath = os.path.join(BASE_DIR, "sequences", "orthologs", f"{ortho_name}.fasta")
                save_fasta(fasta, filepath)
                _, seq = parse_fasta_sequence(fasta)
                print(f"    Length: {len(seq)} aa")
                results["orthologs"][target_name][ortho_name] = {
                    "accession": ortho_info["uniprot"],
                    "species": ortho_info["species"],
                    "length": len(seq),
                }
            else:
                print(f"    FAILED to fetch")
            time.sleep(1)

    # ── Save summary ──
    summary_path = os.path.join(BASE_DIR, "results", "sequence_retrieval_summary.json")
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n\nSummary saved to: {summary_path}")

    # ── Print summary table ──
    print("\n" + "=" * 60)
    print("SEQUENCE RETRIEVAL SUMMARY")
    print("=" * 60)
    print(f"\nTargets fetched: {sum(1 for v in results['targets'].values() if 'length' in v)}/{len(TARGETS)}")
    print(f"Controls fetched: {sum(1 for v in results['controls'].values() if 'length' in v)}/{len(CONTROLS)}")
    total_ortho = sum(len(v) for v in ORTHOLOGS.values())
    fetched_ortho = sum(1 for t in results['orthologs'].values() for v in t.values() if 'length' in v)
    print(f"Orthologs fetched: {fetched_ortho}/{total_ortho}")


if __name__ == "__main__":
    main()
