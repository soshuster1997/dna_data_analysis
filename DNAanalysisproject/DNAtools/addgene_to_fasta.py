#!/usr/bin/env python3
"""
addgene_to_fasta.py
====================
Converts Addgene plasmid/sequence records into individual .fasta files, using
data obtained through Addgene's official Developers Portal (bulk JSON
download or live API) -- see https://developers.addgene.org/access-options/.
This does NOT scrape addgene.org; it requires an approved Developers Portal
account and access token, since Addgene gates programmatic/bulk sequence
access behind that process.

IMPORTANT -- field names are placeholders
------------------------------------------
The API endpoint itself is confirmed (Addgene's bulk plasmids_with_sequences
download, authenticated with a "Token" header). The JSON *field names* inside
each record (ID_FIELD, NAME_FIELD, SEQUENCE_FIELD_CANDIDATES below) are still
my best guess, since I haven't seen a real response. Run with --inspect
first -- it prints the raw keys of the first record so you can fix the
constants in one place before converting everything.

Usage
-----
    # 1) See what fields your data actually has (also works via --api)
    python addgene_to_fasta.py --api --token YOUR_TOKEN --inspect

    # 2) Convert straight from the live API -- randomly samples 1000 of
    #    however many are returned (Addgene's full dataset can be 100,000+)
    python addgene_to_fasta.py --api --token YOUR_TOKEN --out-dir fasta_out --limit 1000

    # Reproducible random sample (same 1000 every time you rerun with this seed)
    python addgene_to_fasta.py --api --token YOUR_TOKEN --out-dir fasta_out --limit 1000 --seed 42

    # Take the first 1000 in Addgene's own order instead of a random sample
    python addgene_to_fasta.py --api --token YOUR_TOKEN --out-dir fasta_out --limit 1000 --sample first

    # 3) Or from a bulk JSON file you've already downloaded/saved locally
    python addgene_to_fasta.py --bulk-json addgene_bulk.json --out-dir fasta_out --limit 1000
"""

from __future__ import annotations

import argparse
import json
import os
import random
import re
import sys
from typing import Any, Dict, Iterable, List, Optional, Tuple

# --------------------------------------------------------------------------
# CONFIG -- adjust these after running --inspect against your real data
# --------------------------------------------------------------------------

ID_FIELD = "id"                 # unique plasmid record identifier
NAME_FIELD = "name"              # human-readable plasmid name

# Confirmed from Addgene's actual response schema: sequence data lives under
# record["sequences"][<one of these 4 lists>][i]["sequence"]. Preference order:
# Addgene-confirmed sequences first, then depositor/user-submitted, and full
# sequences before partial ones.
SEQUENCE_LIST_PRIORITY = [
    "public_addgene_full_sequences",
    "public_user_full_sequences",
    "public_addgene_partial_sequences",
    "public_user_partial_sequences",
]

# Addgene's bulk-download API endpoint (confirmed working endpoint, not paginated --
# it returns the full plasmids-with-sequences dataset in one response).
ADDGENE_DOWNLOAD_URL = "https://api.developers.addgene.org/download/plasmids_with_sequences/"


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def safe_filename(name: str, fallback: str) -> str:
    name = (name or fallback).strip()
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    return name or fallback


def extract_sequence(record: Dict[str, Any]) -> Tuple[Optional[str], Optional[str]]:
    """Pull a nucleotide sequence out of record["sequences"], trying each of
    the 4 possible lists in priority order (Addgene-confirmed full sequence
    first, down to user-submitted partial sequence last). Returns
    (sequence, which_list_it_came_from) so callers can record provenance;
    (None, None) if no usable sequence was found anywhere."""
    seqs = record.get("sequences") or {}
    for list_name in SEQUENCE_LIST_PRIORITY:
        entries = seqs.get(list_name) or []
        for entry in entries:
            seq = entry.get("sequence") if isinstance(entry, dict) else None
            if seq:
                return seq, list_name
    return None, None


def clean_dna(seq: str) -> str:
    return re.sub(r"[^ACGTacgt]", "", seq).upper()


def write_fasta(seq: str, path: str, header: str, wrap: int = 70) -> None:
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), wrap):
            fh.write(seq[i:i + wrap] + "\n")


# --------------------------------------------------------------------------
# Input sources
# --------------------------------------------------------------------------

def parse_records_payload(data: Any) -> List[Dict[str, Any]]:
    """Handle the common shapes a bulk/API JSON payload might come in --
    either a top-level list, or a dict with a 'data'/'plasmids'/'results'/
    'items' list inside."""
    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        for key in ("data", "plasmids", "results", "items"):
            if key in data and isinstance(data[key], list):
                return data[key]
        raise ValueError(
            f"Could not find a records list in the JSON (top-level keys: {list(data.keys())}). "
            "Adjust parse_records_payload() to match your response's structure."
        )
    raise ValueError("Unrecognized JSON structure.")


def records_from_bulk_json(path: str) -> Iterable[Dict[str, Any]]:
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    yield from parse_records_payload(data)


def records_from_api(token: str) -> List[Dict[str, Any]]:
    """Addgene's plasmids_with_sequences endpoint returns the whole dataset
    in a single authenticated GET -- no pagination needed. Returns the full
    list; sampling/limiting happens in main() so --sample random can draw
    from the entire result set rather than just however the API ordered it."""
    try:
        import requests
    except ImportError:
        sys.exit("Live API mode requires the 'requests' package: pip install requests")

    headers = {"Authorization": f"Token {token}", "Accept": "application/json"}
    print(f"Requesting {ADDGENE_DOWNLOAD_URL} ...")
    resp = requests.get(ADDGENE_DOWNLOAD_URL, headers=headers, timeout=120)
    if resp.status_code != 200:
        raise RuntimeError(
            f"Addgene API request failed ({resp.status_code}): {resp.text[:300]}"
        )
    records = parse_records_payload(resp.json())
    print(f"Received {len(records)} record(s) from the API.")
    return records


# --------------------------------------------------------------------------
# Main conversion
# --------------------------------------------------------------------------

def convert(records: Iterable[Dict[str, Any]], out_dir: str, limit: int, inspect: bool) -> None:
    os.makedirs(out_dir, exist_ok=True)
    written = 0
    skipped_no_seq = 0
    manifest_rows = []

    for i, record in enumerate(records):
        if inspect and i == 0:
            print("First record's top-level keys:")
            print(json.dumps(list(record.keys()), indent=2))
            print("\nFull first record (for reference):")
            print(json.dumps(record, indent=2)[:3000])
            print("\n--inspect only shows the first record; re-run without --inspect to convert. Exiting.")
            return

        if written >= limit:
            break

        raw_seq, seq_source = extract_sequence(record)
        if not raw_seq:
            skipped_no_seq += 1
            continue

        seq = clean_dna(raw_seq)
        if not seq:
            skipped_no_seq += 1
            continue

        record_id = str(record.get(ID_FIELD, f"record_{i}"))
        name = str(record.get(NAME_FIELD, record_id))
        header = f"{record_id} {name}".strip()
        filename = safe_filename(f"{record_id}_{name}", f"record_{i}") + ".fasta"
        out_path = os.path.join(out_dir, filename)

        resistance = record.get("bacterial_resistance") or ";".join(record.get("resistance_markers") or [])
        vector_types = ";".join((record.get("cloning") or {}).get("vector_types") or [])

        write_fasta(seq, out_path, header=header)
        manifest_rows.append((filename, record_id, name, len(seq), seq_source, resistance, vector_types))
        written += 1

        if written % 50 == 0:
            print(f"  ...{written} FASTA files written so far")

    manifest_path = os.path.join(out_dir, "manifest.csv")
    with open(manifest_path, "w") as fh:
        fh.write("filename,addgene_id,name,length_bp,sequence_source,bacterial_resistance,vector_types\n")
        for row in manifest_rows:
            fh.write(",".join(str(x).replace(",", " ") for x in row) + "\n")

    print(f"\nDone. Wrote {written} FASTA file(s) to '{out_dir}/'.")
    print(f"Skipped {skipped_no_seq} record(s) with no usable sequence field.")
    print(f"Manifest written to {manifest_path}")


def main():
    parser = argparse.ArgumentParser(description="Convert Addgene records (bulk JSON or live API) into individual FASTA files.")
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--bulk-json", help="Path to a bulk JSON file downloaded from the Addgene Developers Portal.")
    src.add_argument("--api", action="store_true", help="Fetch records live from the Addgene API instead of a bulk file.")

    parser.add_argument("--token", help="Addgene API access token (required with --api).")
    parser.add_argument("--out-dir", default="addgene_fasta", help="Directory to write .fasta files into (default: addgene_fasta).")
    parser.add_argument("--limit", type=int, default=1000, help="Max number of records to convert (default: 1000).")
    parser.add_argument("--sample", choices=["random", "first"], default="random",
                         help="How to pick which records to convert when there are more than --limit available: "
                              "'random' draws a random subset from the whole result set (default), 'first' just "
                              "takes them in the order Addgene returned them.")
    parser.add_argument("--seed", type=int, default=None,
                         help="Random seed for reproducible --sample random selection (default: not fixed).")
    parser.add_argument("--inspect", action="store_true",
                         help="Print the structure of the first record only, then exit -- use this first "
                              "to confirm/fix the field name constants at the top of this script.")
    args = parser.parse_args()

    if args.api and not args.token:
        parser.error("--api requires --token")

    if args.bulk_json:
        records = list(records_from_bulk_json(args.bulk_json))
    else:
        records = records_from_api(args.token)

    if args.inspect:
        convert(records[:1], out_dir=args.out_dir, limit=1, inspect=True)
        return

    total = len(records)
    if total > args.limit:
        if args.sample == "random":
            rng = random.Random(args.seed)
            records = rng.sample(records, args.limit)
            print(f"Randomly sampled {args.limit} of {total} total record(s)"
                  + (f" (seed={args.seed})" if args.seed is not None else " (no fixed seed -- rerun with --seed N to reproduce this exact sample)") + ".")
        else:
            records = records[:args.limit]
            print(f"Took the first {args.limit} of {total} total record(s) (--sample first).")

    convert(records, out_dir=args.out_dir, limit=len(records), inspect=False)


if __name__ == "__main__":
    main()
