#!/usr/bin/env python3
"""
wikipedia_to_dna.py
======================
Fetches random Wikipedia articles, converts each article's text to binary
(UTF-8, 8 bits/byte), then converts that binary to DNA using the same 2-bit
encoding as dna_binary_converter.py (A=00 C=01 G=10 T=11), and writes out one
.fasta file per article. Each sequence is cut off at a random length between
--min-bp and --max-bp base pairs (default 300-5000).

Requires dna_binary_converter.py in the same directory (imports binary_to_dna
from it, so the encoding is guaranteed identical to that tool).

Usage
-----
    python wikipedia_to_dna.py --n 500 --out-dir wiki_dna_fasta
    python wikipedia_to_dna.py --n 500 --min-bp 300 --max-bp 5000 --seed 42
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Optional, Tuple

try:
    from dna_binary_converter import binary_to_dna
except ImportError:
    sys.exit("Could not import dna_binary_converter.py -- make sure it is in the same directory as this script.")

WIKI_API_URL = "https://en.wikipedia.org/w/api.php"
USER_AGENT = "orf-pipeline-wiki-to-dna/1.0 (educational/demo use; contact: none provided)"


# ---------------------------------------------------------------------------
# Wikipedia fetching
# ---------------------------------------------------------------------------

def fetch_random_articles(batch_size: int, max_retries: int = 6) -> List[Tuple[str, str, str]]:
    """Fetch a BATCH of random Wikipedia articles' plain-text extracts in a
    single API call (much friendlier to Wikipedia's rate limits than one
    request per article). Returns a list of (title, extract_text, page_url)
    tuples, skipping any pages with no extract. On 429/503, respects the
    Retry-After header (falling back to exponential backoff starting at 5s,
    per Wikimedia's API etiquette) and retries."""
    params = {
        "action": "query",
        "format": "json",
        "generator": "random",
        "grnnamespace": 0,   # main namespace only (real articles, not talk/category/etc.)
        "grnlimit": batch_size,
        "prop": "extracts",
        "explaintext": 1,
    }
    url = WIKI_API_URL + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})

    backoff = 5.0
    for attempt in range(max_retries):
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode("utf-8"))
            break
        except urllib.error.HTTPError as e:
            if e.code in (429, 503):
                retry_after = e.headers.get("Retry-After")
                wait = float(retry_after) if retry_after else backoff
                print(f"  rate limited ({e.code}), waiting {wait:.0f}s before retrying "
                      f"(attempt {attempt + 1}/{max_retries})...", file=sys.stderr)
                time.sleep(wait)
                backoff = min(backoff * 2, 60.0)
                continue
            print(f"  warning: request failed ({e}), skipping this batch", file=sys.stderr)
            return []
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            print(f"  warning: request failed ({e}), skipping this batch", file=sys.stderr)
            return []
    else:
        print("  giving up on this batch after repeated rate-limit errors", file=sys.stderr)
        return []

    results = []
    pages = data.get("query", {}).get("pages", {})
    for page in pages.values():
        title = page.get("title", "unknown")
        extract = (page.get("extract") or "").strip()
        page_id = page.get("pageid")
        page_url = f"https://en.wikipedia.org/?curid={page_id}" if page_id else ""
        if extract:
            results.append((title, extract, page_url))
    return results


# ---------------------------------------------------------------------------
# Text -> binary -> DNA
# ---------------------------------------------------------------------------

def text_to_binary(text: str) -> str:
    return "".join(f"{byte:08b}" for byte in text.encode("utf-8"))


def safe_filename(name: str) -> str:
    import re
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    return name or "article"


def write_fasta(seq: str, path: str, header: str, wrap: int = 70) -> None:
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), wrap):
            fh.write(seq[i:i + wrap] + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Fetch random Wikipedia articles, convert text -> binary -> DNA, write N FASTA files."
    )
    parser.add_argument("--n", type=int, default=500, help="Number of articles/sequences to generate (default: 500).")
    parser.add_argument("--min-bp", type=int, default=300, help="Minimum sequence length in bp (default: 300).")
    parser.add_argument("--max-bp", type=int, default=5000, help="Maximum sequence length in bp (default: 5000).")
    parser.add_argument("--out-dir", default="wiki_dna_fasta", help="Directory to write .fasta files into.")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducible length cutoffs (default: not fixed).")
    parser.add_argument("--batch-size", type=int, default=20,
                         help="Random articles requested per API call (default: 20). Wikipedia's anonymous-user "
                              "limit is typically 20-50 per request; fewer, larger batches means far fewer total "
                              "requests and much less chance of hitting rate limits than one request per article.")
    parser.add_argument("--delay", type=float, default=1.0,
                         help="Seconds to wait between Wikipedia API batch requests (default: 1.0).")
    parser.add_argument("--max-retries-per-batch", type=int, default=6,
                         help="Max retries (with backoff) per batch request if rate-limited (default: 6).")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    rng = random.Random(args.seed)

    manifest_rows = []
    written = 0
    used_titles = set()
    queue: List[Tuple[str, str, str]] = []
    empty_batch_streak = 0

    while written < args.n:
        if not queue:
            remaining = args.n - written
            batch_size = min(args.batch_size, max(remaining, 1))
            fetched = fetch_random_articles(batch_size, max_retries=args.max_retries_per_batch)
            time.sleep(args.delay)
            new_items = [item for item in fetched if item[0] not in used_titles]
            if not new_items:
                empty_batch_streak += 1
                if empty_batch_streak >= 10:
                    print("  10 consecutive empty/failed batches -- stopping early.", file=sys.stderr)
                    break
                continue
            empty_batch_streak = 0
            queue.extend(new_items)
            continue

        title, extract, page_url = queue.pop()
        used_titles.add(title)

        target_bp = rng.randint(args.min_bp, args.max_bp)
        needed_bits = target_bp * 2
        needed_bytes = -(-needed_bits // 8)  # ceil division: enough bytes to cover needed_bits

        text_for_encoding = extract
        # Wikipedia extracts are usually plenty long, but pad by repeating the
        # text if a short article can't supply enough bytes to hit target_bp.
        while len(text_for_encoding.encode("utf-8")) < needed_bytes:
            text_for_encoding += " " + extract

        binary = text_to_binary(text_for_encoding)
        binary = binary[:needed_bits]  # trim to an exact even-length multiple matching target_bp
        dna_seq = binary_to_dna(binary)
        actual_bp = len(dna_seq)

        filename = f"{written:04d}_{safe_filename(title)}.fasta"
        out_path = os.path.join(args.out_dir, filename)
        header = f"{title} | target_bp={target_bp} actual_bp={actual_bp} source={page_url}"
        write_fasta(dna_seq, out_path, header=header)

        manifest_rows.append([filename, title, target_bp, actual_bp, page_url])
        written += 1

        if written % 25 == 0:
            print(f"  ...{written}/{args.n} sequences written")

    manifest_path = os.path.join(args.out_dir, "manifest.csv")
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["filename", "wikipedia_title", "target_bp", "actual_bp", "source_url"])
        writer.writerows(manifest_rows)

    print(f"\nDone. Wrote {written} FASTA file(s) to '{args.out_dir}/'.")
    print(f"Manifest written to {manifest_path}")


if __name__ == "__main__":
    main()
