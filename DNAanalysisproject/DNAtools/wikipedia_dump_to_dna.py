#!/usr/bin/env python3
"""
wikipedia_dump_to_dna.py
===========================
Alternative to wikipedia_to_dna.py that avoids Wikipedia's live-API rate
limits entirely by reading a Wikipedia XML database dump you've already
downloaded, instead of calling the API per-article. Converts text -> binary
-> DNA using the same 2-bit encoding as dna_binary_converter.py, cutting each
sequence off at a random length between --min-bp and --max-bp.

Why a dump instead of the API
------------------------------
Wikipedia dumps are plain static file downloads from dumps.wikimedia.org --
there's no per-request rate limiting because there are no "requests" beyond
the single file download itself. The tradeoff is you download a file first.

Which dump to get
------------------
Full English Wikipedia's article dump is 25+ GB compressed (100+ GB
decompressed) -- overkill for generating mock/test data. I'd strongly
recommend Simple English Wikipedia instead: same official format, same
official source, but only a few hundred MB compressed, which is plenty to
draw 500 random articles from. Download one of these (right-click -> Save
As, or wget/curl -- these are plain HTTPS file downloads, not an API, so
there's no rate limit to hit):

    Simple English Wikipedia (recommended, much smaller):
    https://dumps.wikimedia.org/simplewiki/latest/simplewiki-latest-pages-articles-multistream.xml.bz2

    Full English Wikipedia (much larger, only if you specifically want it):
    https://dumps.wikimedia.org/enwiki/latest/enwiki-latest-pages-articles-multistream.xml.bz2

Either a `.bz2` compressed dump or an already-decompressed `.xml` file works
directly with this script -- no need to unzip it yourself.

Usage
-----
    python wikipedia_dump_to_dna.py --dump-file simplewiki-latest-pages-articles-multistream.xml.bz2 \
        --n 500 --min-bp 300 --max-bp 5000 --out-dir wiki_dna_fasta --seed 42
"""

from __future__ import annotations

import argparse
import bz2
import os
import random
import re
import sys
import xml.etree.ElementTree as ET
from typing import List, Optional, Tuple

try:
    from dna_binary_converter import binary_to_dna
except ImportError:
    sys.exit("Could not import dna_binary_converter.py -- make sure it is in the same directory as this script.")

MIN_ARTICLE_CHARS = 200  # skip very short stubs/disambiguation pages


def localname(tag: str) -> str:
    """Strip the MediaWiki export namespace prefix from an ElementTree tag,
    e.g. '{http://www.mediawiki.org/xml/export-0.10/}page' -> 'page'."""
    return tag.split("}", 1)[-1] if "}" in tag else tag


def open_dump(path: str):
    if path.endswith(".bz2"):
        return bz2.open(path, "rb")
    return open(path, "rb")


def iter_articles(dump_path: str, max_scan: Optional[int]):
    """Stream-parse a MediaWiki XML dump, yielding (title, text) for
    main-namespace, non-redirect pages with at least MIN_ARTICLE_CHARS of
    text. Uses iterparse + element.clear() so memory stays flat regardless
    of dump size."""
    scanned = 0
    with open_dump(dump_path) as fh:
        context = ET.iterparse(fh, events=("end",))
        page_ns = None
        page_title = None
        page_redirect = False
        page_text = None

        for event, elem in context:
            tag = localname(elem.tag)

            if tag == "ns":
                page_ns = elem.text
            elif tag == "title":
                page_title = elem.text
            elif tag == "redirect":
                page_redirect = True
            elif tag == "text":
                page_text = elem.text
            elif tag == "page":
                scanned += 1
                if page_ns == "0" and not page_redirect and page_text and len(page_text) >= MIN_ARTICLE_CHARS:
                    yield page_title or "untitled", page_text
                # reset for next page, and free memory
                page_ns, page_title, page_redirect, page_text = None, None, False, None
                elem.clear()
                if max_scan and scanned >= max_scan:
                    return


def reservoir_sample(iterable, k: int, rng: random.Random) -> List:
    """Algorithm R: unbiased random sample of size k from a stream of
    unknown length, without holding the whole stream in memory."""
    reservoir = []
    for i, item in enumerate(iterable):
        if i < k:
            reservoir.append(item)
        else:
            j = rng.randint(0, i)
            if j < k:
                reservoir[j] = item
    return reservoir


# ---------------------------------------------------------------------------
# Text -> binary -> DNA (identical logic to wikipedia_to_dna.py)
# ---------------------------------------------------------------------------

def text_to_binary(text: str) -> str:
    return "".join(f"{byte:08b}" for byte in text.encode("utf-8", errors="ignore"))


def safe_filename(name: str) -> str:
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    return name or "article"


def write_fasta(seq: str, path: str, header: str, wrap: int = 70) -> None:
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), wrap):
            fh.write(seq[i:i + wrap] + "\n")


def build_sequence(primary_text: str, pool: List[Tuple[str, str]], target_bp: int, rng: random.Random) -> str:
    """Encode primary_text to DNA, truncated/padded to exactly target_bp.
    If primary_text alone isn't long enough, pad with OTHER random articles
    from the pool (not repeats of the same text) to avoid periodic
    repetition artifacts in the resulting DNA."""
    needed_bits = target_bp * 2
    needed_bytes = -(-needed_bits // 8)

    text = primary_text
    guard = 0
    while len(text.encode("utf-8")) < needed_bytes and guard < 50:
        _, extra_text = pool[rng.randrange(len(pool))]
        text += " " + extra_text
        guard += 1

    binary = text_to_binary(text)[:needed_bits]
    return binary_to_dna(binary)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Convert random articles from a downloaded Wikipedia XML dump into DNA FASTA files (no API rate limits)."
    )
    parser.add_argument("--dump-file", required=True,
                         help="Path to a downloaded MediaWiki XML dump (.xml or .xml.bz2), e.g. "
                              "simplewiki-latest-pages-articles-multistream.xml.bz2")
    parser.add_argument("--n", type=int, default=500, help="Number of sequences to generate (default: 500).")
    parser.add_argument("--min-bp", type=int, default=300, help="Minimum sequence length in bp (default: 300).")
    parser.add_argument("--max-bp", type=int, default=5000, help="Maximum sequence length in bp (default: 5000).")
    parser.add_argument("--out-dir", default="wiki_dna_fasta", help="Directory to write .fasta files into.")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducible sampling/length cutoffs.")
    parser.add_argument("--pool-size", type=int, default=None,
                         help="Size of the random article pool to draw from (default: 3x --n). Larger pools give "
                              "more variety, especially for padding donor text, at the cost of more memory.")
    parser.add_argument("--max-scan", type=int, default=None,
                         help="Stop after scanning this many <page> elements in the dump (default: no limit, scan "
                              "the whole file). Use this to cap runtime on very large dumps like full enwiki; not "
                              "needed for Simple English Wikipedia, which is small enough to scan fully quickly.")
    args = parser.parse_args()

    if not os.path.exists(args.dump_file):
        sys.exit(f"Dump file not found: {args.dump_file}\n\n"
                  f"Download one first, e.g.:\n"
                  f"  wget https://dumps.wikimedia.org/simplewiki/latest/simplewiki-latest-pages-articles-multistream.xml.bz2\n")

    os.makedirs(args.out_dir, exist_ok=True)
    rng = random.Random(args.seed)
    pool_size = args.pool_size or max(args.n * 3, 100)

    print(f"Scanning {args.dump_file} for qualifying articles (this streams the file; may take a while for large dumps)...")
    pool = reservoir_sample(iter_articles(args.dump_file, args.max_scan), pool_size, rng)
    if len(pool) < args.n:
        sys.exit(f"Only found {len(pool)} qualifying articles in the dump, need at least --n={args.n}. "
                  f"Try a bigger dump, remove --max-scan, or lower --n.")
    print(f"Collected a pool of {len(pool)} articles. Generating {args.n} sequences...")

    primary_indices = rng.sample(range(len(pool)), args.n)
    manifest_rows = []

    for i, idx in enumerate(primary_indices):
        title, text = pool[idx]
        target_bp = rng.randint(args.min_bp, args.max_bp)
        dna_seq = build_sequence(text, pool, target_bp, rng)
        actual_bp = len(dna_seq)

        filename = f"{i:04d}_{safe_filename(title)}.fasta"
        out_path = os.path.join(args.out_dir, filename)
        header = f"{title} | target_bp={target_bp} actual_bp={actual_bp} source=wikipedia_dump"
        write_fasta(dna_seq, out_path, header=header)
        manifest_rows.append([filename, title, target_bp, actual_bp])

        if (i + 1) % 50 == 0:
            print(f"  ...{i + 1}/{args.n} sequences written")

    manifest_path = os.path.join(args.out_dir, "manifest.csv")
    with open(manifest_path, "w") as fh:
        fh.write("filename,wikipedia_title,target_bp,actual_bp\n")
        for row in manifest_rows:
            fh.write(",".join(str(x).replace(",", " ") for x in row) + "\n")

    print(f"\nDone. Wrote {len(manifest_rows)} FASTA file(s) to '{args.out_dir}/'.")
    print(f"Manifest written to {manifest_path}")


if __name__ == "__main__":
    main()
