"""Run STRING v12 enrichment; write one table, without assigning module annotations.

Uses frozen input identifiers, species 10029, and the default STRING proteome
background. STRING's returned FDR is retained without cross-module correction.
"""
from pathlib import Path
import argparse
import json
import time
import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
API = "https://version-12-0.string-db.org/api/json/enrichment"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-inputs", action="store_true")
    parser.add_argument("--output", type=Path, default=HERE / "discovery_string_functional_enrichment_results.csv")
    args = parser.parse_args()
    mapping = pd.read_csv(HERE.parent / "data/discovery_string_protein_mapping.csv")
    column = "STRING_identifier" if "STRING_identifier" in mapping else "String_name"
    assert not mapping.duplicated(["iModulon", "Gene"]).any()
    groups = list(mapping.groupby("iModulon", sort=True))
    if args.check_inputs:
        print(f"{len(groups)} modules; {mapping[column].notna().sum()} mapped membership rows")
        return
    rows = []
    for module, group in groups:
        identifiers = sorted(set(group[column].dropna().astype(str).str.strip()) - {""})
        count = len(identifiers)
        terms = []
        if count >= 2:
            for attempt in range(3):
                try:
                    response = requests.post(API, data={
                        "identifiers": "\r".join(identifiers), "species": 10029,
                        "caller_identity": "CHO_iModulon_reproducibility",
                    }, timeout=90)
                    response.raise_for_status()
                    terms = response.json()
                    if not isinstance(terms, list) or any("Error" in r for r in terms):
                        raise ValueError(str(terms)[:300])
                    break
                except (requests.RequestException, ValueError):
                    if attempt == 2:
                        raise
                    time.sleep(3 * (attempt + 1))
            time.sleep(1)
        status = "completed" if count >= 2 else "skipped_fewer_than_two_identifiers"
        for term in terms or [{}]:
            record = dict(term)
            for key, value in record.items():
                if isinstance(value, list):
                    record[key] = json.dumps(value)
            record.update(iModulon=int(module), n_identifiers=count, status=status)
            record["FDR05"] = record.get("fdr", 1) <= .05
            record["nonliterature_FDR05"] = record["FDR05"] and record.get("category") != "PMID"
            rows.append(record)
        print(f"iM{module}: {count} identifiers, {len(terms)} terms", flush=True)
    # Do not replace an existing result if any API request failed.
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
