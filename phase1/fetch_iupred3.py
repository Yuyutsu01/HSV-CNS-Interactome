"""
Fetch IUPred3 disorder predictions via the web interface.
Replicates the browser flow: POST accession -> GET plot_api -> parse Bokeh JSON.
"""
import requests
import json
import os
import time
import csv
import sys

BASE_DIR = r"C:\Users\shiva\OneDrive\Desktop\HSV-CNS-Interactome\phase1"
OUT_DIR = os.path.join(BASE_DIR, "iupred3_data")
os.makedirs(OUT_DIR, exist_ok=True)

PROTEINS = {
    # Targets
    "ICP34.5": "P36313",
    "US3": "P04413",
    "pUL37": "P10221",
    # Controls
    "ICP8_UL29": "P04296",
    "UL42": "P10186",
    "gD_US6": "P57083",
    "VP11_12_UL46": "P10228",
    "VP16_UL48": "P06492",
    "VP22_UL49": "P10216",
    "ICP27_UL54": "P10238",
    "ICP22_US1": "P04485",
}


def fetch_iupred3(name, accession, session):
    """Fetch IUPred3 predictions for one protein."""
    print(f"\n  {name} ({accession})...", end=" ", flush=True)

    # Step 1: POST the accession to /plot
    resp = session.post(
        "https://iupred3.elte.hu/plot",
        data={"accession": accession, "iupred_type": "long", "anchor_type": "anchor2"},
        timeout=60,
    )
    if resp.status_code != 200:
        print(f"POST failed: {resp.status_code}")
        return None

    # Step 2: GET plot_api (same session, server remembers the query)
    time.sleep(1)
    api_resp = session.get("https://iupred3.elte.hu/plot_api", timeout=60)
    if api_resp.status_code != 200:
        print(f"plot_api failed: {api_resp.status_code}")
        return None

    try:
        api_data = api_resp.json()
    except json.JSONDecodeError:
        print(f"JSON parse failed (got {len(api_resp.text)} bytes)")
        return None

    # Step 3: Parse Bokeh JSON to extract x,y data
    if "grid" not in api_data:
        print(f"No 'grid' key (keys: {list(api_data.keys())})")
        return None

    grid = json.loads(api_data["grid"])
    refs = grid["doc"]["roots"]["references"]

    # Find ColumnDataSource with x,y arrays (the IUPred line data)
    iupred_data = None
    anchor_data = None
    for ref in refs:
        attrs = ref.get("attributes", {})
        data = attrs.get("data", {})
        if "x" in data and "y" in data:
            x = data["x"]
            y = data["y"]
            if isinstance(x, list) and len(x) > 10:
                if iupred_data is None:
                    iupred_data = (x, y)
                elif anchor_data is None:
                    anchor_data = (x, y)

    if iupred_data is None:
        print("No IUPred data found in Bokeh JSON")
        return None

    x_vals, y_vals = iupred_data
    print(f"OK ({len(x_vals)} residues)")

    # Save to CSV
    out_path = os.path.join(OUT_DIR, f"{name}_iupred3.csv")
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["position", "iupred3_score"])
        for pos, score in zip(x_vals, y_vals):
            writer.writerow([int(pos), f"{score:.6f}"])

    # Also save ANCHOR2 if available
    if anchor_data:
        ax, ay = anchor_data
        anchor_path = os.path.join(OUT_DIR, f"{name}_anchor2.csv")
        with open(anchor_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["position", "anchor2_score"])
            for pos, score in zip(ax, ay):
                writer.writerow([int(pos), f"{score:.6f}"])
        print(f"    + ANCHOR2 saved ({len(ax)} residues)")

    return {"name": name, "accession": accession, "length": len(x_vals), "file": out_path}


def main():
    # Check which proteins we already have
    existing = set()
    for name in PROTEINS:
        path = os.path.join(OUT_DIR, f"{name}_iupred3.csv")
        if os.path.exists(path):
            with open(path) as f:
                lines = f.readlines()
            if len(lines) > 10:
                existing.add(name)
                print(f"  {name}: already have {len(lines)-1} residues, skipping")

    to_fetch = {k: v for k, v in PROTEINS.items() if k not in existing}
    if not to_fetch:
        print("\nAll proteins already fetched!")
        return

    print(f"\nFetching {len(to_fetch)} proteins from IUPred3...")

    session = requests.Session()
    session.headers.update({
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    })

    # Initialize session by visiting the main page
    session.get("https://iupred3.elte.hu/", timeout=30)
    time.sleep(1)

    results = []
    for name, acc in to_fetch.items():
        result = fetch_iupred3(name, acc, session)
        if result:
            results.append(result)
        else:
            print(f"    FAILED - will need manual entry")
        time.sleep(3)

    # Summary
    print(f"\n{'='*60}")
    print(f"Successfully fetched: {len(results)}/{len(to_fetch)}")
    for r in results:
        print(f"  {r['name']}: {r['length']} residues -> {r['file']}")

    failed = [n for n in to_fetch if n not in {r['name'] for r in results}]
    if failed:
        print(f"\nFailed: {', '.join(failed)}")


if __name__ == "__main__":
    main()
