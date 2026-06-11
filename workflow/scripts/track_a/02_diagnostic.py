#!/usr/bin/env python3
"""
Quick diagnostic: query Ensembl REST API for CASP8 orthologues in chicken.
Run this and paste the output back to Claude.
"""
import requests
import json

ENSEMBL_REST = "https://rest.ensembl.org"

# Step 1: Look up CASP8 Ensembl ID
print("=== Step 1: Lookup CASP8 ===")
r1 = requests.get(
    f"{ENSEMBL_REST}/lookup/symbol/homo_sapiens/CASP8",
    headers={"Accept": "application/json"},
    timeout=20
)
print(f"Status: {r1.status_code}")
print(f"Content-Type returned: {r1.headers.get('Content-Type', 'NONE')}")

if r1.status_code != 200:
    print(f"Body (first 500 chars): {r1.text[:500]}")
    exit()

lookup = r1.json()
ensembl_id = lookup.get("id", "NOT_FOUND")
print(f"Ensembl ID: {ensembl_id}")

# Step 2: Query orthologues in chicken
print(f"\n=== Step 2: Orthologues for {ensembl_id} in chicken ===")

# Try with Accept header (correct)
r2 = requests.get(
    f"{ENSEMBL_REST}/homology/id/{ensembl_id}",
    params={"target_species": "gallus_gallus", "type": "orthologues"},
    headers={"Accept": "application/json"},
    timeout=20
)
print(f"Status: {r2.status_code}")
print(f"Content-Type returned: {r2.headers.get('Content-Type', 'NONE')}")

if r2.status_code != 200:
    print(f"Body (first 500 chars): {r2.text[:500]}")
else:
    data = r2.json()
    print(f"\nTop-level keys: {list(data.keys())}")
    data_list = data.get("data", [])
    print(f"data[] length: {len(data_list)}")
    if data_list:
        first = data_list[0]
        print(f"data[0] keys: {list(first.keys())}")
        homologies = first.get("homologies", [])
        print(f"homologies count: {len(homologies)}")
        if homologies:
            print(f"\nFirst homology (pretty):")
            print(json.dumps(homologies[0], indent=2)[:1500])
        else:
            print("\nNO HOMOLOGIES FOUND — printing full response:")
            print(json.dumps(data, indent=2)[:3000])
    else:
        print("\ndata[] is empty — printing full response:")
        print(json.dumps(data, indent=2)[:3000])

# Step 3: Try without target_species filter
print(f"\n=== Step 3: All orthologues (no species filter) ===")
r3 = requests.get(
    f"{ENSEMBL_REST}/homology/id/{ensembl_id}",
    params={"type": "orthologues"},
    headers={"Accept": "application/json"},
    timeout=20
)
print(f"Status: {r3.status_code}")
if r3.status_code == 200:
    data3 = r3.json()
    homologies3 = data3.get("data", [{}])[0].get("homologies", [])
    print(f"Total orthologues (all species): {len(homologies3)}")
    chicken_hits = [h for h in homologies3
                    if "gallus" in h.get("target", {}).get("species", "").lower()]
    print(f"Chicken hits: {len(chicken_hits)}")
    if chicken_hits:
        print(json.dumps(chicken_hits[0], indent=2)[:1000])
