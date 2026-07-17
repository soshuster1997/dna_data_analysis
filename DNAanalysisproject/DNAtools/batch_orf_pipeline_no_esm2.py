#!/usr/bin/env python3
"""
batch_orf_pipeline_no_esm2.py
================================
A variant of batch_orf_pipeline.py that:

  - Does NOT use ESM-2 at all (no torch/transformers dependency, no model
    download, much faster -- good for large batches or environments without
    GPU/internet access to Hugging Face).
  - Finds ORFs starting with ATG *and* common alternative bacterial start
    codons (GTG, TTG, CTG by default), flagging which ones used a
    non-traditional start so you can see how much that setting matters.
  - Adds two sequence-complexity checks per ORF and per whole file:
      * Shannon entropy (order-0, over A/C/G/T frequencies, 0-2 bits)
      * A k-mer repeat measure: how much of the sequence is made of
        k-mers (default k=6) that also occur elsewhere in the same
        sequence -- a simple, fast proxy for repetitiveness/low complexity.

Non-traditional start codons -- what this means
--------------------------------------------------
Biologically, the ribosome loads an initiator Met-tRNA at the start codon
regardless of whether that codon is ATG, GTG, or TTG -- so the resulting
protein still begins with methionine even when the start codon itself would
normally code for Val (GTG) or Leu (TTG/CTG). This script follows that
convention: for alternative-start ORFs, the first residue is reported as M,
and the rest of the ORF is translated normally. This mirrors how NCBI
ORFfinder's "ATG and alternative initiation codons" option behaves.

Usage
-----
    python batch_orf_pipeline_no_esm2.py --fasta-dir addgene_fasta --out-dir batch_results
    python batch_orf_pipeline_no_esm2.py --fasta-dir addgene_fasta --out-dir batch_results --atg-only
    python batch_orf_pipeline_no_esm2.py --fasta-dir addgene_fasta --out-dir batch_results \
        --start-codons ATG,GTG,TTG,CTG,ATT --kmer-size 8
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import statistics
import sys
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from orf_analyzer import (
    STANDARD_CODON_TABLE,
    STOP_CODONS,
    ECOLI_CODON_FREQ,
    reverse_complement,
    codon_adaptation_index,
    read_sequence,
)

FASTA_EXTENSIONS = (".fasta", ".fa", ".fna", ".txt")
TRADITIONAL_START = "ATG"
DEFAULT_ALT_STARTS = ["GTG", "TTG", "CTG"]  # common bacterial alternative initiation codons


# ---------------------------------------------------------------------------
# ORF data structure
# ---------------------------------------------------------------------------

@dataclass
class ORF:
    frame: int
    strand: str
    start_nt: int
    end_nt: int
    nt_seq: str
    protein: str
    start_codon: str
    is_nonstandard_start: bool
    length_aa: int = field(init=False)
    cai_percent: Optional[float] = None
    shannon_entropy: Optional[float] = None
    repeat_count: Optional[int] = None
    repeat_fraction: Optional[float] = None

    def __post_init__(self):
        self.length_aa = len(self.protein)


def translate_with_alt_start(nt_seq: str, start_codon: str) -> str:
    """Translate nt_seq to protein, stopping before the first stop codon.
    If start_codon is a non-ATG initiator, the first residue is forced to M
    (initiator Met), matching real translation biology -- see module docstring."""
    protein = []
    for i in range(0, len(nt_seq) - 2, 3):
        codon = nt_seq[i:i + 3]
        aa = STANDARD_CODON_TABLE.get(codon, "X")
        if aa == "*":
            break
        protein.append(aa)
    if protein and start_codon != TRADITIONAL_START:
        protein[0] = "M"
    return "".join(protein)


# ---------------------------------------------------------------------------
# ORF finding (parameterized start-codon set)
# ---------------------------------------------------------------------------

def _find_orfs_on_strand(seq: str, strand: str, min_aa: int, full_len: int,
                          start_codons: Tuple[str, ...]) -> List[ORF]:
    orfs: List[ORF] = []
    working_seq = seq if strand == "+" else reverse_complement(seq)

    for frame in range(3):
        frame_num = frame + 1 if strand == "+" else -(frame + 1)
        i = frame
        while i <= len(working_seq) - 3:
            codon = working_seq[i:i + 3]
            if codon in start_codons:
                j = i
                found_stop = False
                while j <= len(working_seq) - 3:
                    c = working_seq[j:j + 3]
                    if c in STOP_CODONS:
                        found_stop = True
                        break
                    j += 3
                if found_stop:
                    nt_seq = working_seq[i:j + 3]
                    protein = translate_with_alt_start(nt_seq, codon)
                    if len(protein) > min_aa:
                        if strand == "+":
                            start_nt = i + 1
                            end_nt = j + 3
                        else:
                            start_nt = full_len - (j + 3) + 1
                            end_nt = full_len - i
                        orfs.append(ORF(
                            frame=frame_num, strand=strand, start_nt=start_nt, end_nt=end_nt,
                            nt_seq=nt_seq, protein=protein, start_codon=codon,
                            is_nonstandard_start=(codon != TRADITIONAL_START),
                        ))
                    i = j + 3
                    continue
                else:
                    break
            i += 3
    return orfs


def find_all_orfs(seq: str, min_aa: int, start_codons: Tuple[str, ...]) -> List[ORF]:
    full_len = len(seq)
    orfs = _find_orfs_on_strand(seq, "+", min_aa, full_len, start_codons)
    orfs += _find_orfs_on_strand(seq, "-", min_aa, full_len, start_codons)
    orfs.sort(key=lambda o: o.length_aa, reverse=True)
    return orfs


# ---------------------------------------------------------------------------
# Sequence complexity checks: Shannon entropy + k-mer repeat content
# ---------------------------------------------------------------------------

def shannon_entropy(seq: str) -> Optional[float]:
    """Order-0 Shannon entropy over A/C/G/T frequency, in bits (0-2).
    2.0 = perfectly uniform base composition; lower = more skewed/predictable."""
    if not seq:
        return None
    counts = Counter(seq)
    n = len(seq)
    h = 0.0
    for base in "ACGT":
        p = counts.get(base, 0) / n
        if p > 0:
            h -= p * math.log2(p)
    return h


def kmer_repeat_stats(seq: str, k: int = 6) -> Tuple[int, Optional[float]]:
    """Count how many k-mers in seq are repeats of an earlier-seen k-mer
    (a simple, fast repetitiveness proxy -- not a full tandem-repeat finder).
    Returns (repeat_count, repeat_fraction), where repeat_fraction is
    repeat_count / total_kmers. None fraction if the sequence is shorter
    than k."""
    if len(seq) < k:
        return 0, None
    seen: Dict[str, int] = {}
    repeat_count = 0
    total = 0
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        total += 1
        if seen.get(kmer, 0) > 0:
            repeat_count += 1
        seen[kmer] = seen.get(kmer, 0) + 1
    fraction = repeat_count / total if total else None
    return repeat_count, fraction


# ---------------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------------

def find_fasta_files(fasta_dir: str) -> List[Path]:
    return sorted(p for p in Path(fasta_dir).iterdir() if p.is_file() and p.suffix.lower() in FASTA_EXTENSIONS)


def describe(values: List[float]) -> Dict[str, Optional[float]]:
    if not values:
        return {"n": 0, "mean": None, "median": None, "stdev": None, "min": None, "max": None}
    return {
        "n": len(values), "mean": statistics.mean(values), "median": statistics.median(values),
        "stdev": statistics.stdev(values) if len(values) > 1 else 0.0,
        "min": min(values), "max": max(values),
    }


def process_batch(fasta_dir: str, min_aa: int, start_codons: Tuple[str, ...], kmer_size: int, progress_every: int):
    files = find_fasta_files(fasta_dir)
    if not files:
        sys.exit(f"No FASTA files (.fasta/.fa/.fna/.txt) found in '{fasta_dir}'.")

    per_file_stats: List[Dict] = []
    all_orf_rows: List[List] = []
    failed_files: List[str] = []

    start = time.time()
    for i, path in enumerate(files, start=1):
        try:
            seq = read_sequence(str(path))
            orfs = find_all_orfs(seq, min_aa=min_aa, start_codons=start_codons)
            for orf in orfs:
                orf.cai_percent = codon_adaptation_index(orf.nt_seq)
                orf.shannon_entropy = shannon_entropy(orf.nt_seq)
                orf.repeat_count, orf.repeat_fraction = kmer_repeat_stats(orf.nt_seq, k=kmer_size)
        except Exception as e:
            failed_files.append(f"{path.name}: {e}")
            continue

        seq_len = len(seq)
        seq_entropy = shannon_entropy(seq)
        seq_repeat_count, seq_repeat_fraction = kmer_repeat_stats(seq, k=kmer_size)

        cai_vals = [o.cai_percent for o in orfs if o.cai_percent is not None]
        n_nonstandard = sum(1 for o in orfs if o.is_nonstandard_start)

        per_file_stats.append({
            "file": path.name,
            "seq_len_bp": seq_len,
            "n_orfs": len(orfs),
            "orfs_per_kb": (len(orfs) / seq_len) * 1000 if seq_len else 0.0,
            "mean_cai_percent": statistics.mean(cai_vals) if cai_vals else None,
            "n_nonstandard_start_orfs": n_nonstandard,
            "pct_nonstandard_start_orfs": (n_nonstandard / len(orfs) * 100) if orfs else None,
            "seq_shannon_entropy": seq_entropy,
            "seq_repeat_count": seq_repeat_count,
            "seq_repeat_fraction": seq_repeat_fraction,
        })

        for idx, orf in enumerate(orfs, start=1):
            all_orf_rows.append([
                path.name, idx, orf.strand, orf.frame, orf.start_nt, orf.end_nt, orf.length_aa,
                orf.start_codon, orf.is_nonstandard_start,
                f"{orf.cai_percent:.2f}" if orf.cai_percent is not None else "",
                f"{orf.shannon_entropy:.4f}" if orf.shannon_entropy is not None else "",
                orf.repeat_count if orf.repeat_count is not None else "",
                f"{orf.repeat_fraction:.4f}" if orf.repeat_fraction is not None else "",
            ])

        if i % progress_every == 0 or i == len(files):
            print(f"  processed {i}/{len(files)} files ({time.time() - start:.1f}s elapsed)")

    return per_file_stats, all_orf_rows, failed_files


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_csv(path: str, header: List[str], rows: List[List]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)


def write_per_file_csv(path: str, per_file_stats: List[Dict]) -> None:
    header = ["file", "seq_len_bp", "n_orfs", "orfs_per_kb", "mean_cai_percent",
              "n_nonstandard_start_orfs", "pct_nonstandard_start_orfs",
              "seq_shannon_entropy", "seq_repeat_count", "seq_repeat_fraction"]
    rows = []
    for r in per_file_stats:
        rows.append([
            r["file"], r["seq_len_bp"], r["n_orfs"], f'{r["orfs_per_kb"]:.3f}',
            f'{r["mean_cai_percent"]:.2f}' if r["mean_cai_percent"] is not None else "",
            r["n_nonstandard_start_orfs"],
            f'{r["pct_nonstandard_start_orfs"]:.1f}' if r["pct_nonstandard_start_orfs"] is not None else "",
            f'{r["seq_shannon_entropy"]:.4f}' if r["seq_shannon_entropy"] is not None else "",
            r["seq_repeat_count"],
            f'{r["seq_repeat_fraction"]:.4f}' if r["seq_repeat_fraction"] is not None else "",
        ])
    write_csv(path, header, rows)


def build_summary_report(fasta_dir: str, per_file_stats: List[Dict], all_orf_rows: List[List],
                          failed_files: List[str], min_aa: int, start_codons: Tuple[str, ...],
                          kmer_size: int) -> str:
    n_files = len(per_file_stats)
    total_bp = sum(r["seq_len_bp"] for r in per_file_stats)
    total_orfs = sum(r["n_orfs"] for r in per_file_stats)
    total_nonstandard = sum(r["n_nonstandard_start_orfs"] for r in per_file_stats)

    orfs_per_file_stats = describe([r["n_orfs"] for r in per_file_stats])
    density_stats = describe([r["orfs_per_kb"] for r in per_file_stats])
    seq_entropy_stats = describe([r["seq_shannon_entropy"] for r in per_file_stats if r["seq_shannon_entropy"] is not None])
    seq_repeat_stats = describe([r["seq_repeat_fraction"] for r in per_file_stats if r["seq_repeat_fraction"] is not None])

    all_cai = [float(row[9]) for row in all_orf_rows if row[9] != ""]
    all_orf_entropy = [float(row[10]) for row in all_orf_rows if row[10] != ""]
    all_orf_repeat_frac = [float(row[12]) for row in all_orf_rows if row[12] != ""]
    cai_stats = describe(all_cai)
    orf_entropy_stats = describe(all_orf_entropy)
    orf_repeat_stats = describe(all_orf_repeat_frac)

    files_with_zero_orfs = sum(1 for r in per_file_stats if r["n_orfs"] == 0)

    def fmt(s: Dict, d: int = 2) -> str:
        if s["n"] == 0:
            return "n=0 (no data)"
        return (f"n={s['n']}, mean={s['mean']:.{d}f}, median={s['median']:.{d}f}, "
                f"stdev={s['stdev']:.{d}f}, min={s['min']:.{d}f}, max={s['max']:.{d}f}")

    lines = []
    lines.append("# Batch ORF Analysis Summary (no ESM-2)")
    lines.append("")
    lines.append(f"- **Input directory**: {fasta_dir}")
    lines.append(f"- **FASTA files processed**: {n_files}" + (f" ({len(failed_files)} failed, see below)" if failed_files else ""))
    lines.append(f"- **Total bases across all files**: {total_bp:,} bp")
    lines.append(f"- **Start codons searched**: {', '.join(start_codons)} "
                 f"(non-ATG codons treated as alternative/non-traditional starts)")
    lines.append(f"- **Minimum ORF length filter**: > {min_aa} amino acids")
    lines.append(f"- **k-mer size for repeat detection**: {kmer_size} bp")
    lines.append(f"- **Total ORFs found**: {total_orfs}")
    lines.append(f"- **ORFs with a non-traditional (non-ATG) start codon**: {total_nonstandard} "
                 f"({total_nonstandard/total_orfs*100:.1f}% of all ORFs)" if total_orfs else "")
    lines.append(f"- **Files with zero qualifying ORFs**: {files_with_zero_orfs}")
    lines.append("")

    lines.append("## Descriptive Statistics")
    lines.append("")
    lines.append(f"- **ORFs per file**: {fmt(orfs_per_file_stats)}")
    lines.append(f"- **ORF density (per kb)**: {fmt(density_stats, 3)}")
    lines.append(f"- **CAI codon match (%), across all {cai_stats['n']} ORFs**: {fmt(cai_stats, 1)}")
    lines.append(f"- **Whole-sequence Shannon entropy (bits, 0-2)**: {fmt(seq_entropy_stats, 3)}")
    lines.append(f"- **Whole-sequence repeat fraction (k={kmer_size})**: {fmt(seq_repeat_stats, 3)}")
    lines.append(f"- **Per-ORF Shannon entropy (bits, 0-2)**: {fmt(orf_entropy_stats, 3)}")
    lines.append(f"- **Per-ORF repeat fraction (k={kmer_size})**: {fmt(orf_repeat_stats, 3)}")
    lines.append("")

    if failed_files:
        lines.append("## Files That Failed to Process")
        lines.append("")
        for msg in failed_files:
            lines.append(f"- {msg}")
        lines.append("")

    lines.append("## Per-File Results (first 50)")
    lines.append("")
    lines.append("(Full detail in `per_file_summary.csv`; individual ORFs in `all_orfs.csv`.)")
    lines.append("")
    lines.append("| File | Length (bp) | # ORFs | ORFs/kb | Mean CAI (%) | Non-ATG starts | Entropy | Repeat frac |")
    lines.append("|------|-------------|--------|---------|---------------|-----------------|---------|--------------|")
    for r in per_file_stats[:50]:
        cai_str = f'{r["mean_cai_percent"]:.1f}' if r["mean_cai_percent"] is not None else "N/A"
        ent_str = f'{r["seq_shannon_entropy"]:.3f}' if r["seq_shannon_entropy"] is not None else "N/A"
        rep_str = f'{r["seq_repeat_fraction"]:.3f}' if r["seq_repeat_fraction"] is not None else "N/A"
        lines.append(f'| {r["file"]} | {r["seq_len_bp"]} | {r["n_orfs"]} | {r["orfs_per_kb"]:.2f} | '
                     f'{cai_str} | {r["n_nonstandard_start_orfs"]} | {ent_str} | {rep_str} |')
    if len(per_file_stats) > 50:
        lines.append(f"\n*(showing first 50 of {len(per_file_stats)} files -- see per_file_summary.csv for all)*")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Batch ORF/CAI/entropy/repeat analysis across a directory of FASTA files (no ESM-2)."
    )
    parser.add_argument("--fasta-dir", required=True, help="Directory containing .fasta/.fa/.fna files.")
    parser.add_argument("--out-dir", default="batch_results_no_esm2", help="Output directory.")
    parser.add_argument("--min-aa", type=int, default=80, help="Minimum ORF length in amino acids (default: 80).")
    parser.add_argument("--start-codons", default=None,
                         help="Comma-separated list of start codons to search for (default: "
                              f"ATG,{','.join(DEFAULT_ALT_STARTS)}).")
    parser.add_argument("--atg-only", action="store_true",
                         help="Only search for traditional ATG starts (equivalent to --start-codons ATG).")
    parser.add_argument("--kmer-size", type=int, default=6,
                         help="k-mer size used for the repeat-content check (default: 6).")
    parser.add_argument("--progress-every", type=int, default=25, help="Print progress every N files (default: 25).")
    args = parser.parse_args()

    if args.atg_only:
        start_codons = (TRADITIONAL_START,)
    elif args.start_codons:
        start_codons = tuple(c.strip().upper() for c in args.start_codons.split(","))
    else:
        start_codons = (TRADITIONAL_START, *DEFAULT_ALT_STARTS)

    os.makedirs(args.out_dir, exist_ok=True)

    per_file_stats, all_orf_rows, failed_files = process_batch(
        args.fasta_dir, args.min_aa, start_codons, args.kmer_size, args.progress_every
    )

    write_per_file_csv(os.path.join(args.out_dir, "per_file_summary.csv"), per_file_stats)
    write_csv(
        os.path.join(args.out_dir, "all_orfs.csv"),
        ["file", "orf_index", "strand", "frame", "start_nt", "end_nt", "length_aa",
         "start_codon", "is_nonstandard_start", "cai_percent", "shannon_entropy",
         "repeat_count", "repeat_fraction"],
        all_orf_rows,
    )

    report = build_summary_report(args.fasta_dir, per_file_stats, all_orf_rows, failed_files,
                                   args.min_aa, start_codons, args.kmer_size)
    report_path = os.path.join(args.out_dir, "batch_summary_report.md")
    with open(report_path, "w") as fh:
        fh.write(report)

    print(f"\nWrote {report_path}")
    print(f"Wrote {os.path.join(args.out_dir, 'per_file_summary.csv')}")
    print(f"Wrote {os.path.join(args.out_dir, 'all_orfs.csv')}")


if __name__ == "__main__":
    main()
