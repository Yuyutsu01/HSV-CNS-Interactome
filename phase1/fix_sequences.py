"""Fix incorrect accessions and re-fetch failed orthologs."""
import requests
import os
import json
import time

BASE = r"C:\Users\shiva\OneDrive\Desktop\HSV-CNS-Interactome\phase1\sequences"
RESULTS = r"C:\Users\shiva\OneDrive\Desktop\HSV-CNS-Interactome\phase1\results"

fixes = {
    ("targets", "pUL37"): ("P10221", "HSV-1 strain 17"),
    ("orthologs", "HSV2_UL37"): ("G9I259", "HSV-2 (HG52)"),
    ("orthologs", "PRV_US3"): ("A0A0S2MM84", "PRV (Becker)"),
    ("orthologs", "PRV_UL37"): ("Q911W0", "PRV"),
}

results = {}

for (subdir, name), (acc, species) in fixes.items():
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    r = requests.get(url, timeout=30)
    if r.status_code == 200 and r.text.startswith(">"):
        path = os.path.join(BASE, subdir, f"{name}.fasta")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            f.write(r.text.strip() + "\n")
        seq = "".join(r.text.strip().split("\n")[1:])
        print(f"{name} ({acc}, {species}): {len(seq)} aa - OK")
        results[name] = {"accession": acc, "species": species, "length": len(seq)}
    else:
        print(f"{name} ({acc}): FAILED (status {r.status_code})")
        results[name] = {"accession": acc, "error": "fetch failed"}
    time.sleep(1)

# Check pUL37 disorder annotations
print("\n--- pUL37 (P10221) disorder annotations ---")
r = requests.get("https://rest.uniprot.org/uniprotkb/P10221.json", timeout=30)
if r.status_code == 200:
    data = r.json()
    disorder = []
    for feat in data.get("features", []):
        if feat.get("type") == "Region" and "disordered" in feat.get("description", "").lower():
            s = feat["location"]["start"]["value"]
            e = feat["location"]["end"]["value"]
            disorder.append(f"{s}-{e}")
    if disorder:
        print(f"Disorder regions: {', '.join(disorder)}")
    else:
        print("No UniProt disorder annotations (expected — Tier 3 confidence)")

# Save corrections
with open(os.path.join(RESULTS, "sequence_corrections.json"), "w") as f:
    json.dump(results, f, indent=2)

# Verify all target sequences
print("\n--- Final target sequence verification ---")
for name in ["ICP34.5", "US3", "pUL37"]:
    path = os.path.join(BASE, "targets", f"{name}.fasta")
    if os.path.exists(path):
        with open(path) as f:
            lines = f.read().strip().split("\n")
            header = lines[0]
            seq = "".join(lines[1:])
            print(f"{name}: {len(seq)} aa | {header[:80]}")
    else:
        print(f"{name}: FILE NOT FOUND")
