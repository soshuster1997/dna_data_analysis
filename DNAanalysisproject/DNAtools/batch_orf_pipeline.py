#!/usr/bin/env python3
"""
batch_orf_pipeline.py
=======================
Runs orf_analyzer.py's ORF-finding + CAI + ESM-2 pipeline across every FASTA
file in a directory (e.g. the output of addgene_to_fasta.py), and produces:

  1. per_file_summary.csv   -- one row per FASTA file: seq length, # ORFs,
                                ORF-to-bp ratio, mean CAI, mean perplexity
  2. all_orfs.csv           -- one row per individual ORF found, across all files
  3. batch_summary_report.md -- descriptive statistics across the whole batch
                                (mean/median/stdev of ORFs per file, CAI, and
                                perplexity, plus totals)

Unlike calling orf_analyzer.py once per file from the shell, this script
loads the ESM-2 model ONCE and reuses it for every file -- important when
you're processing hundreds to thousands of files.

Usage
-----
    python batch_orf_pipeline.py --fasta-dir addgene_fasta --out-dir batch_results
    python batch_orf_pipeline.py --fasta-dir addgene_fasta --out-dir batch_results --no-esm2
    python batch_orf_pipeline.py --fasta-dir addgene_fasta --out-dir batch_results \
        --esm-model facebook/esm2_t12_35M_UR50D --min-aa 80
"""

from __future__ import annotations

import argparse
import csv
import os
import statistics
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import orf_analyzer as oa

FASTA_EXTENSIONS = (".fasta", ".fa", ".fna", ".txt")


def find_fasta_files(fasta_dir: str) -> List[Path]:
    return sorted(
        p for p in Path(fasta_dir).iterdir()
        if p.is_file() and p.suffix.lower() in FASTA_EXTENSIONS
    )


def describe(values: List[float]) -> Dict[str, Optional[float]]:
    """Basic descriptive statistics, tolerant of empty/singleton lists."""
    if not values:
        return {"n": 0, "mean": None, "median": None, "stdev": None, "min": None, "max": None}
    return {
        "n": len(values),
        "mean": statistics.mean(values),
        "median": statistics.median(values),
        "stdev": statistics.stdev(values) if len(values) > 1 else 0.0,
        "min": min(values),
        "max": max(values),
    }


def process_batch(
    fasta_dir: str,
    min_aa: int,
    esm_model: str,
    no_esm2: bool,
    accurate_perplexity: bool,
    progress_every: int = 25,
):
    files = find_fasta_files(fasta_dir)
    if not files:
        sys.exit(f"No FASTA files (.fasta/.fa/.fna/.txt) found in '{fasta_dir}'.")

    scorer = None
    esm2_error = None
    if not no_esm2:
        print(f"Loading ESM-2 model '{esm_model}' once for the whole batch...")
        try:
            scorer = oa.ESM2Scorer(model_name=esm_model)
            print("Model loaded.")
        except Exception as e:
            esm2_error = str(e)
            print(f"Warning: ESM-2 unavailable ({esm2_error}). "
                  f"Continuing without perplexity scoring.", file=sys.stderr)

    per_file_stats: List[Dict] = []
    all_orf_rows: List[List] = []
    failed_files: List[str] = []

    start = time.time()
    for i, path in enumerate(files, start=1):
        try:
            seq = oa.read_sequence(str(path))
            orfs = oa.find_all_orfs(seq, min_aa=min_aa)
            for orf in orfs:
                orf.cai_percent = oa.codon_adaptation_index(orf.nt_seq)
                if scorer is not None:
                    try:
                        orf.esm2_perplexity = scorer.perplexity(orf.protein, fast=not accurate_perplexity)
                    except Exception:
                        orf.esm2_perplexity = None
        except Exception as e:
            failed_files.append(f"{path.name}: {e}")
            continue

        cai_vals = [o.cai_percent for o in orfs if o.cai_percent is not None]
        ppl_vals = [o.esm2_perplexity for o in orfs if o.esm2_perplexity is not None]
        seq_len = len(seq)

        per_file_stats.append({
            "file": path.name,
            "seq_len_bp": seq_len,
            "n_orfs": len(orfs),
            "orf_to_bp_ratio": (len(orfs) / seq_len) if seq_len else 0.0,
            "orfs_per_kb": (len(orfs) / seq_len) * 1000 if seq_len else 0.0,
            "mean_cai_percent": statistics.mean(cai_vals) if cai_vals else None,
            "mean_esm2_perplexity": statistics.mean(ppl_vals) if ppl_vals else None,
        })

        for idx, orf in enumerate(orfs, start=1):
            all_orf_rows.append([
                path.name, idx, orf.strand, orf.frame, orf.start_nt, orf.end_nt,
                orf.length_aa,
                f"{orf.cai_percent:.2f}" if orf.cai_percent is not None else "",
                f"{orf.esm2_perplexity:.4f}" if orf.esm2_perplexity is not None else "",
            ])

        if i % progress_every == 0 or i == len(files):
            elapsed = time.time() - start
            print(f"  processed {i}/{len(files)} files ({elapsed:.1f}s elapsed)")

    return per_file_stats, all_orf_rows, failed_files, esm2_error


def write_csv(path: str, header: List[str], rows: List[List]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)


def write_per_file_csv(path: str, per_file_stats: List[Dict]) -> None:
    header = ["file", "seq_len_bp", "n_orfs", "orf_to_bp_ratio", "orfs_per_kb",
              "mean_cai_percent", "mean_esm2_perplexity"]
    rows = []
    for row in per_file_stats:
        rows.append([
            row["file"], row["seq_len_bp"], row["n_orfs"],
            f'{row["orf_to_bp_ratio"]:.6f}', f'{row["orfs_per_kb"]:.3f}',
            f'{row["mean_cai_percent"]:.2f}' if row["mean_cai_percent"] is not None else "",
            f'{row["mean_esm2_perplexity"]:.4f}' if row["mean_esm2_perplexity"] is not None else "",
        ])
    write_csv(path, header, rows)


def build_summary_report(
    fasta_dir: str,
    per_file_stats: List[Dict],
    all_orf_rows: List[List],
    failed_files: List[str],
    esm2_error: Optional[str],
    esm_model: str,
    min_aa: int,
) -> str:
    n_files = len(per_file_stats)
    total_bp = sum(r["seq_len_bp"] for r in per_file_stats)
    total_orfs = sum(r["n_orfs"] for r in per_file_stats)

    orfs_per_file_stats = describe([r["n_orfs"] for r in per_file_stats])
    density_stats = describe([r["orfs_per_kb"] for r in per_file_stats])

    all_cai = [float(row[7]) for row in all_orf_rows if row[7] != ""]
    all_ppl = [float(row[8]) for row in all_orf_rows if row[8] != ""]
    cai_stats = describe(all_cai)
    ppl_stats = describe(all_ppl)

    files_with_zero_orfs = sum(1 for r in per_file_stats if r["n_orfs"] == 0)

    def fmt_stats(s: Dict, decimals: int = 2) -> str:
        if s["n"] == 0:
            return "n=0 (no data)"
        return (f"n={s['n']}, mean={s['mean']:.{decimals}f}, median={s['median']:.{decimals}f}, "
                f"stdev={s['stdev']:.{decimals}f}, min={s['min']:.{decimals}f}, max={s['max']:.{decimals}f}")

    lines = []
    lines.append("# Batch ORF Analysis Summary")
    lines.append("")
    lines.append(f"- **Input directory**: {fasta_dir}")
    lines.append(f"- **FASTA files processed**: {n_files}" + (f" ({len(failed_files)} failed, see below)" if failed_files else ""))
    lines.append(f"- **Total bases across all files**: {total_bp:,} bp")
    lines.append(f"- **Total ORFs found (> {min_aa} aa)**: {total_orfs}")
    lines.append(f"- **Files with zero qualifying ORFs**: {files_with_zero_orfs}")
    lines.append(f"- **ESM-2 model**: {esm_model}" + (" (unavailable, see note)" if esm2_error else ""))
    lines.append("")

    if esm2_error:
        lines.append(f"> **Note:** ESM-2 scoring was unavailable ({esm2_error}). "
                      f"Perplexity statistics below reflect only files where scoring succeeded, if any.")
        lines.append("")

    lines.append("## Descriptive Statistics")
    lines.append("")
    lines.append(f"- **ORFs per file**: {fmt_stats(orfs_per_file_stats, decimals=2)}")
    lines.append(f"- **ORF density (ORFs per kb)**: {fmt_stats(density_stats, decimals=3)}")
    lines.append(f"- **CAI codon match (%), across all {cai_stats['n']} ORFs**: {fmt_stats(cai_stats, decimals=1)}")
    lines.append(f"- **ESM-2 perplexity, across all {ppl_stats['n']} ORFs**: {fmt_stats(ppl_stats, decimals=2)}")
    lines.append("")

    if failed_files:
        lines.append("## Files That Failed to Process")
        lines.append("")
        for msg in failed_files:
            lines.append(f"- {msg}")
        lines.append("")

    lines.append("## Per-File Results")
    lines.append("")
    lines.append("(Full detail in `per_file_summary.csv`; individual ORFs in `all_orfs.csv`.)")
    lines.append("")
    lines.append("| File | Length (bp) | # ORFs | ORFs/kb | Mean CAI (%) | Mean Perplexity |")
    lines.append("|------|-------------|--------|---------|---------------|-------------------|")
    for row in per_file_stats[:50]:
        cai_str = f'{row["mean_cai_percent"]:.1f}' if row["mean_cai_percent"] is not None else "N/A"
        ppl_str = f'{row["mean_esm2_perplexity"]:.2f}' if row["mean_esm2_perplexity"] is not None else "N/A"
        lines.append(f'| {row["file"]} | {row["seq_len_bp"]} | {row["n_orfs"]} | '
                     f'{row["orfs_per_kb"]:.2f} | {cai_str} | {ppl_str} |')
    if len(per_file_stats) > 50:
        lines.append(f"\n*(showing first 50 of {len(per_file_stats)} files -- see per_file_summary.csv for all)*")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Run orf_analyzer's pipeline across a directory of FASTA files and summarize results.")
    parser.add_argument("--fasta-dir", required=True, help="Directory containing .fasta/.fa/.fna files.")
    parser.add_argument("--out-dir", default="batch_results", help="Directory to write CSVs and the summary report (default: batch_results).")
    parser.add_argument("--min-aa", type=int, default=80, help="Minimum ORF length in amino acids (default: 80).")
    parser.add_argument("--esm-model", default="facebook/esm2_t6_8M_UR50D", help="ESM-2 model to use (loaded once for the whole batch).")
    parser.add_argument("--no-esm2", action="store_true", help="Skip ESM-2 scoring entirely (much faster).")
    parser.add_argument("--accurate-perplexity", action="store_true", help="Use slower true masked-residue perplexity instead of the fast approximation.")
    parser.add_argument("--progress-every", type=int, default=25, help="Print progress every N files (default: 25).")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    per_file_stats, all_orf_rows, failed_files, esm2_error = process_batch(
        args.fasta_dir, args.min_aa, args.esm_model, args.no_esm2,
        args.accurate_perplexity, args.progress_every,
    )

    write_per_file_csv(os.path.join(args.out_dir, "per_file_summary.csv"), per_file_stats)
    write_csv(
        os.path.join(args.out_dir, "all_orfs.csv"),
        ["file", "orf_index", "strand", "frame", "start_nt", "end_nt", "length_aa", "cai_percent", "esm2_perplexity"],
        all_orf_rows,
    )

    report = build_summary_report(
        args.fasta_dir, per_file_stats, all_orf_rows, failed_files, esm2_error, args.esm_model, args.min_aa
    )
    report_path = os.path.join(args.out_dir, "batch_summary_report.md")
    with open(report_path, "w") as fh:
        fh.write(report)

    print(f"\nWrote {report_path}")
    print(f"Wrote {os.path.join(args.out_dir, 'per_file_summary.csv')}")
    print(f"Wrote {os.path.join(args.out_dir, 'all_orfs.csv')}")


if __name__ == "__main__":
    main()
