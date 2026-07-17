#!/usr/bin/env python3
"""
visualize_orfs.py
==================
Companion script to orf_analyzer.py. Reads the Markdown report produced by
orf_analyzer.py (or re-runs the analysis directly) and plots histograms of:

  1. ESM-2 pseudo-perplexity across all detected ORFs
  2. E. coli codon usage match (CAI %) across all detected ORFs

Usage
-----
    # Parse an existing report (fast, no re-computation, no ESM-2 needed)
    python visualize_orfs.py --report report.md --out-prefix orf_plots

    # Or re-run the full analysis directly from a FASTA/sequence and plot
    python visualize_orfs.py --input my_sequence.fasta --out-prefix orf_plots
    python visualize_orfs.py --seq ATGCGT... --out-prefix orf_plots

Outputs (written next to --out-prefix):
    <prefix>_perplexity_hist.png
    <prefix>_cai_hist.png
    <prefix>_combined.png   (both side by side)
"""

from __future__ import annotations

import argparse
import re
import sys
from typing import List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")  # headless-safe; works in scripts, servers, and notebooks alike
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def parse_report(report_path: str) -> Tuple[List[Optional[float]], List[Optional[float]]]:
    """Extract ESM-2 perplexity and CAI% values from the summary table of a
    Markdown report produced by orf_analyzer.py. Rows are of the form:

        | 1 | + | 3 | 51 | 341 | 96 | 12.34 | 87.5% |

    'N/A' entries are kept as None so the caller can decide whether to warn
    or silently skip them.
    """
    perplexities: List[Optional[float]] = []
    cai_percents: List[Optional[float]] = []

    with open(report_path, "r") as fh:
        lines = fh.readlines()

    in_table = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("| #"):
            in_table = True
            continue
        if in_table:
            if not stripped.startswith("|"):
                break  # table ended
            if set(stripped.replace("|", "").replace("-", "").strip()) == set():
                continue  # the "|---|---|" separator row
            cells = [c.strip() for c in stripped.strip("|").split("|")]
            if len(cells) < 8:
                continue
            ppl_cell, cai_cell = cells[6], cells[7]

            ppl_val = None if ppl_cell == "N/A" else _to_float(ppl_cell)
            cai_val = None if cai_cell == "N/A" else _to_float(cai_cell.rstrip("%"))

            perplexities.append(ppl_val)
            cai_percents.append(cai_val)

    return perplexities, cai_percents


def _to_float(s: str) -> Optional[float]:
    try:
        return float(s)
    except ValueError:
        return None


def run_analysis_directly(source: str, min_aa: int, esm_model: str, no_esm2: bool,
                           accurate_perplexity: bool) -> Tuple[List[Optional[float]], List[Optional[float]]]:
    """Re-run orf_analyzer's pipeline in-process (used when the user passes
    --input/--seq instead of an existing --report)."""
    try:
        import orf_analyzer as oa
    except ImportError:
        sys.exit("Could not import orf_analyzer.py — make sure it is in the same directory "
                 "as visualize_orfs.py, or on your PYTHONPATH.")

    seq = oa.read_sequence(source)
    orfs = oa.find_all_orfs(seq, min_aa=min_aa)
    for orf in orfs:
        orf.cai_percent = oa.codon_adaptation_index(orf.nt_seq)

    if not no_esm2:
        err = oa.score_orfs_with_esm2(orfs, model_name=esm_model, fast=not accurate_perplexity)
        if err:
            print(f"Warning: ESM-2 scoring unavailable ({err}). "
                  f"Perplexity histogram will be skipped.", file=sys.stderr)

    perplexities = [orf.esm2_perplexity for orf in orfs]
    cai_percents = [orf.cai_percent for orf in orfs]
    return perplexities, cai_percents


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _clean(values: List[Optional[float]]) -> List[float]:
    return [v for v in values if v is not None]


def plot_histogram(values: List[float], title: str, xlabel: str, color: str,
                    out_path: str, bins: int = 15) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(values, bins=min(bins, max(len(values), 1)), color=color, edgecolor="white", alpha=0.85)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel("Number of ORFs", fontsize=11)
    ax.axvline(sum(values) / len(values), color="black", linestyle="--", linewidth=1,
               label=f"mean = {sum(values)/len(values):.2f}")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_combined(perplexities: List[float], cai_percents: List[float], out_path: str, bins: int = 15) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    if perplexities:
        axes[0].hist(perplexities, bins=min(bins, max(len(perplexities), 1)),
                      color="#4C72B0", edgecolor="white", alpha=0.85)
        axes[0].axvline(sum(perplexities) / len(perplexities), color="black", linestyle="--", linewidth=1,
                          label=f"mean = {sum(perplexities)/len(perplexities):.2f}")
        axes[0].legend(frameon=False)
        axes[0].set_title("ESM-2 Pseudo-Perplexity", fontsize=13, fontweight="bold")
        axes[0].set_xlabel("Perplexity (lower = more protein-like)")
    else:
        axes[0].text(0.5, 0.5, "No ESM-2 perplexity\ndata available", ha="center", va="center",
                      transform=axes[0].transAxes, fontsize=11, color="gray")
        axes[0].set_title("ESM-2 Pseudo-Perplexity", fontsize=13, fontweight="bold")
    axes[0].set_ylabel("Number of ORFs")
    axes[0].spines["top"].set_visible(False)
    axes[0].spines["right"].set_visible(False)

    if cai_percents:
        axes[1].hist(cai_percents, bins=min(bins, max(len(cai_percents), 1)),
                      color="#DD8452", edgecolor="white", alpha=0.85)
        axes[1].axvline(sum(cai_percents) / len(cai_percents), color="black", linestyle="--", linewidth=1,
                          label=f"mean = {sum(cai_percents)/len(cai_percents):.1f}%")
        axes[1].legend(frameon=False)
        axes[1].set_title("E. coli Codon Usage Match (CAI %)", fontsize=13, fontweight="bold")
        axes[1].set_xlabel("CAI match (%)")
    else:
        axes[1].text(0.5, 0.5, "No CAI data available", ha="center", va="center",
                      transform=axes[1].transAxes, fontsize=11, color="gray")
        axes[1].set_title("E. coli Codon Usage Match (CAI %)", fontsize=13, fontweight="bold")
    axes[1].set_ylabel("Number of ORFs")
    axes[1].spines["top"].set_visible(False)
    axes[1].spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Plot histograms of ESM-2 perplexity and E. coli CAI codon-match "
                    "percentages for ORFs found by orf_analyzer.py."
    )
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--report", help="Path to an existing Markdown report from orf_analyzer.py.")
    src.add_argument("--input", help="Path to a FASTA/text file (re-runs the full analysis).")
    src.add_argument("--seq", help="Raw DNA sequence (re-runs the full analysis).")

    parser.add_argument("--min-aa", type=int, default=80, help="Only used with --input/--seq (default: 80).")
    parser.add_argument("--esm-model", default="facebook/esm2_t6_8M_UR50D",
                         help="Only used with --input/--seq.")
    parser.add_argument("--no-esm2", action="store_true", help="Only used with --input/--seq.")
    parser.add_argument("--accurate-perplexity", action="store_true", help="Only used with --input/--seq.")

    parser.add_argument("--out-prefix", default="orf_plots",
                         help="Filename prefix for saved plots (default: orf_plots).")
    parser.add_argument("--bins", type=int, default=15, help="Number of histogram bins (default: 15).")
    args = parser.parse_args()

    if args.report:
        perplexities, cai_percents = parse_report(args.report)
    else:
        source = args.input if args.input else args.seq
        perplexities, cai_percents = run_analysis_directly(
            source, args.min_aa, args.esm_model, args.no_esm2, args.accurate_perplexity
        )

    n_total = max(len(perplexities), len(cai_percents))
    clean_ppl = _clean(perplexities)
    clean_cai = _clean(cai_percents)

    print(f"Loaded {n_total} ORF(s): {len(clean_ppl)} with perplexity values, "
          f"{len(clean_cai)} with CAI values.")

    if not clean_ppl and not clean_cai:
        sys.exit("No usable perplexity or CAI values found — nothing to plot.")

    if clean_ppl:
        plot_histogram(clean_ppl, "ESM-2 Pseudo-Perplexity Across ORFs",
                        "Perplexity (lower = more protein-like)", "#4C72B0",
                        f"{args.out_prefix}_perplexity_hist.png", bins=args.bins)
        print(f"Wrote {args.out_prefix}_perplexity_hist.png")
    else:
        print("Skipping perplexity histogram: no ESM-2 values available.")

    if clean_cai:
        plot_histogram(clean_cai, "E. coli Codon Usage Match (CAI %) Across ORFs",
                        "CAI match (%)", "#DD8452",
                        f"{args.out_prefix}_cai_hist.png", bins=args.bins)
        print(f"Wrote {args.out_prefix}_cai_hist.png")
    else:
        print("Skipping CAI histogram: no CAI values available.")

    plot_combined(clean_ppl, clean_cai, f"{args.out_prefix}_combined.png", bins=args.bins)
    print(f"Wrote {args.out_prefix}_combined.png")


if __name__ == "__main__":
    main()
