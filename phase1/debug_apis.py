"""Debug API connectivity for IUPred3 and AlphaFold2."""
import requests

# Test AlphaFold DB API
print("=== AlphaFold2 API ===")
for acc in ["P04413", "P36313", "P10221"]:
    try:
        r = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=15)
        print(f"{acc}: status={r.status_code}")
        if r.status_code == 200:
            data = r.json()
            if isinstance(data, list) and data:
                pdb_url = data[0].get("pdbUrl", "none")
                print(f"  pdbUrl: {pdb_url}")
                # Try downloading
                r2 = requests.head(pdb_url, timeout=10)
                print(f"  pdb HEAD: status={r2.status_code}")
    except Exception as e:
        print(f"{acc}: {e}")

# Try IUPred3 with proper web scraping approach
print("\n=== IUPred3 Alternative ===")
# IUPred3 may need browser headers
headers = {
    "User-Agent": "Mozilla/5.0",
    "Accept": "application/json",
}
try:
    r = requests.get("https://iupred3.elte.hu", headers=headers, timeout=15)
    print(f"Main page: status={r.status_code}")
except Exception as e:
    print(f"Main page: {e}")

# Try pip installable iupred alternative
print("\n=== Checking for local IUPred alternatives ===")
try:
    import importlib
    for mod in ["iupred", "iupred3", "localcider"]:
        spec = importlib.util.find_spec(mod)
        print(f"  {mod}: {'installed' if spec else 'not installed'}")
except Exception as e:
    print(f"  Error: {e}")
