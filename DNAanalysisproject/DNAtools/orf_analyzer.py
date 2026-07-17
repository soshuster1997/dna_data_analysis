#!/usr/bin/env python3
"""
orf_analyzer.py
================
A self-contained ORFfinder-style tool that:

  1. Scans a DNA sequence in all 6 reading frames (3 forward, 3 reverse
     complement) and reports every ORF that runs from a start codon (ATG)
     to an in-frame stop codon and translates to > 80 amino acids
     (this mirrors NCBI ORFfinder's "ATG only" + minimum-length behavior,
     see https://www.ncbi.nlm.nih.gov/orffinder/).
  2. Scores each candidate protein with an ESM-2 protein language model
     (https://github.com/facebookresearch/esm), reporting a pseudo-perplexity
     -- lower is "more natural / more likely to be a real protein" in ESM-2's
     view.
  3. Compares the codon usage of each ORF's nucleotide sequence against a
     reference E. coli (K-12 / W3110) codon usage table and reports a
     Codon Adaptation Index (CAI)-based percentage match.
  4. Writes a Markdown report summarizing all of the above.

Usage
-----
    python orf_analyzer.py --input sequence.fasta --out report.md
    python orf_analyzer.py --seq ATGCGT...  --out report.md
    python orf_analyzer.py --input seq.fasta --min-aa 80 --esm-model facebook/esm2_t12_35M_UR50D

If ESM-2 (via the `transformers` + `torch` packages) is not installed or a
model cannot be downloaded (no internet access), the tool still runs the
ORF-finding and codon-usage analysis, and reports the perplexity column as
"N/A" with an explanation, rather than crashing.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# 1. Genetic code / translation
# ---------------------------------------------------------------------------

STANDARD_CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def translate(nt_seq: str) -> str:
    """Translate a nucleotide sequence (length multiple of 3) to protein,
    stopping at (and excluding) the first stop codon encountered."""
    protein = []
    for i in range(0, len(nt_seq) - 2, 3):
        codon = nt_seq[i:i + 3]
        aa = STANDARD_CODON_TABLE.get(codon, "X")
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)


# ---------------------------------------------------------------------------
# 2. Reference E. coli codon usage table
#    Source: Kazusa Codon Usage Database, E. coli W3110 [gbbct],
#    4332 CDS's / 1,372,057 codons (http://www.kazusa.or.jp/codon/)
#    Values are frequency per 1000 codons.
# ---------------------------------------------------------------------------

ECOLI_CODON_FREQ: Dict[str, float] = {
    "TTT": 22.2, "TTC": 16.5, "TTA": 13.8, "TTG": 13.6,
    "CTT": 11.0, "CTC": 11.1, "CTA": 3.8, "CTG": 53.1,
    "ATT": 30.4, "ATC": 25.2, "ATA": 4.2, "ATG": 27.8,
    "GTT": 18.2, "GTC": 15.3, "GTA": 10.9, "GTG": 26.3,
    "TCT": 8.4, "TCC": 8.6, "TCA": 7.0, "TCG": 8.9,
    "CCT": 7.0, "CCC": 5.5, "CCA": 8.4, "CCG": 23.4,
    "ACT": 8.8, "ACC": 23.5, "ACA": 6.9, "ACG": 14.4,
    "GCT": 15.2, "GCC": 25.7, "GCA": 20.1, "GCG": 33.9,
    "TAT": 16.1, "TAC": 12.2, "TAA": 2.0, "TAG": 0.2,
    "CAT": 13.0, "CAC": 9.8, "CAA": 15.4, "CAG": 29.0,
    "AAT": 17.6, "AAC": 21.6, "AAA": 33.6, "AAG": 10.3,
    "GAT": 32.2, "GAC": 19.1, "GAA": 39.7, "GAG": 18.0,
    "TGT": 5.1, "TGC": 6.4, "TGA": 0.9, "TGG": 15.2,
    "CGT": 21.0, "CGC": 22.3, "CGA": 3.5, "CGG": 5.4,
    "AGT": 8.7, "AGC": 16.1, "AGA": 2.0, "AGG": 1.1,
    "GGT": 24.7, "GGC": 29.8, "GGA": 7.9, "GGG": 11.0,
}

# Group codons by the amino acid (or stop) they encode, for CAI weighting.
_AA_TO_CODONS: Dict[str, List[str]] = {}
for _codon, _aa in STANDARD_CODON_TABLE.items():
    _AA_TO_CODONS.setdefault(_aa, []).append(_codon)

# Relative adaptiveness w(codon) = freq(codon) / max(freq of synonyms),
# per amino acid family (stop codons excluded from CAI, as is standard).
CODON_RELATIVE_ADAPTIVENESS: Dict[str, float] = {}
for _aa, _codons in _AA_TO_CODONS.items():
    if _aa == "*":
        continue
    max_freq = max(ECOLI_CODON_FREQ[c] for c in _codons)
    for c in _codons:
        CODON_RELATIVE_ADAPTIVENESS[c] = ECOLI_CODON_FREQ[c] / max_freq if max_freq > 0 else 0.0


def codon_adaptation_index(nt_seq: str) -> Optional[float]:
    """Compute the Codon Adaptation Index (Sharp & Li, 1987) of an ORF's
    nucleotide sequence against the E. coli reference table, expressed as
    a percentage match (0-100). Stop codon and any incomplete/ambiguous
    trailing codon are ignored. Met/Trp (single-codon families) contribute
    w=1.0 and are informative for length but not discriminative."""
    log_w_sum = 0.0
    n = 0
    for i in range(0, len(nt_seq) - 2, 3):
        codon = nt_seq[i:i + 3].upper()
        w = CODON_RELATIVE_ADAPTIVENESS.get(codon)
        if w is None or codon in STOP_CODONS:
            continue
        if w <= 0:
            w = 1e-6  # avoid log(0); codon essentially unused in E. coli
        log_w_sum += math.log(w)
        n += 1
    if n == 0:
        return None
    cai = math.exp(log_w_sum / n)
    return cai * 100.0


# ---------------------------------------------------------------------------
# 3. Sequence input handling
# ---------------------------------------------------------------------------

def clean_sequence(raw: str) -> str:
    """Uppercase and strip anything that isn't A/C/G/T (whitespace, digits,
    FASTA headers already removed)."""
    raw = raw.upper()
    return re.sub(r"[^ACGT]", "", raw)


def read_sequence(source: str) -> str:
    """`source` may be a path to a FASTA (or plain-text) file, or a raw
    DNA sequence string typed directly on the command line."""
    looks_like_path = ("\n" not in source) and (
        source.lower().endswith((".fa", ".fasta", ".fna", ".txt")) or "/" in source or "\\" in source
    )
    text = source
    if looks_like_path:
        try:
            with open(source, "r") as fh:
                text = fh.read()
        except FileNotFoundError:
            pass  # fall back to treating `source` itself as sequence text

    lines = text.splitlines()
    seq_lines = [ln for ln in lines if not ln.startswith(">")]
    seq = clean_sequence("".join(seq_lines)) if seq_lines else clean_sequence(text)
    if not seq:
        raise ValueError("No valid A/C/G/T sequence found in input.")
    return seq


# ---------------------------------------------------------------------------
# 4. ORF finding across all 6 frames
# ---------------------------------------------------------------------------

@dataclass
class ORF:
    frame: int              # 1,2,3 (forward) or -1,-2,-3 (reverse)
    strand: str              # "+" or "-"
    start_nt: int             # 1-based, in the ORIGINAL input sequence coordinates
    end_nt: int               # 1-based, inclusive, in ORIGINAL input sequence coordinates
    nt_seq: str
    protein: str
    length_aa: int = field(init=False)
    cai_percent: Optional[float] = None
    esm2_perplexity: Optional[float] = None

    def __post_init__(self):
        self.length_aa = len(self.protein)


def _find_orfs_on_strand(seq: str, strand: str, min_aa: int, full_len: int) -> List[ORF]:
    orfs: List[ORF] = []
    working_seq = seq if strand == "+" else reverse_complement(seq)

    for frame in range(3):
        frame_num = frame + 1 if strand == "+" else -(frame + 1)
        i = frame
        while i <= len(working_seq) - 3:
            codon = working_seq[i:i + 3]
            if codon == START_CODON:
                # walk forward from this ATG to the next in-frame stop
                j = i
                found_stop = False
                while j <= len(working_seq) - 3:
                    c = working_seq[j:j + 3]
                    if c in STOP_CODONS:
                        found_stop = True
                        break
                    j += 3
                if found_stop:
                    nt_seq = working_seq[i:j + 3]  # includes stop codon
                    protein = translate(nt_seq)
                    if len(protein) > min_aa:
                        if strand == "+":
                            start_nt = i + 1
                            end_nt = j + 3
                        else:
                            # convert reverse-complement coordinates back to
                            # 1-based positions on the ORIGINAL forward sequence
                            start_nt = full_len - (j + 3) + 1
                            end_nt = full_len - i
                        orfs.append(
                            ORF(
                                frame=frame_num,
                                strand=strand,
                                start_nt=start_nt,
                                end_nt=end_nt,
                                nt_seq=nt_seq,
                                protein=protein,
                            )
                        )
                    # continue scanning after this ORF's stop codon so
                    # overlapping/nested ATGs in the same frame are still found
                    i = j + 3
                    continue
                else:
                    # no stop codon found before the end of the sequence
                    break
            i += 3
    return orfs


def find_all_orfs(seq: str, min_aa: int = 80) -> List[ORF]:
    full_len = len(seq)
    orfs = _find_orfs_on_strand(seq, "+", min_aa, full_len)
    orfs += _find_orfs_on_strand(seq, "-", min_aa, full_len)
    orfs.sort(key=lambda o: o.length_aa, reverse=True)
    return orfs


# ---------------------------------------------------------------------------
# 5. ESM-2 perplexity scoring (optional; requires `transformers` + `torch`)
# ---------------------------------------------------------------------------

class ESM2Scorer:
    """Wraps an ESM-2 masked-language model to compute a pseudo-perplexity
    for a protein sequence. Lower perplexity = the model finds the sequence
    more "natural" (higher likelihood under the training distribution of
    real proteins)."""

    def __init__(self, model_name: str = "facebook/esm2_t6_8M_UR50D", device: Optional[str] = None):
        try:
            import torch
            from transformers import AutoTokenizer, AutoModelForMaskedLM
        except ImportError as e:
            raise RuntimeError(
                "ESM-2 scoring requires the 'transformers' and 'torch' packages. "
                "Install them with: pip install torch transformers"
            ) from e

        self.torch = torch
        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = AutoModelForMaskedLM.from_pretrained(model_name)
        self.model.to(self.device)
        self.model.eval()

    def perplexity(self, protein_seq: str, fast: bool = True) -> float:
        """
        fast=True  : single forward pass, log-likelihood of each residue read
                     off the unmasked logits at its own position. This is an
                     approximation (the model can "see" the true residue) but
                     is O(1) forward passes and fine for relative comparison
                     across many ORFs.
        fast=False : true masked pseudo-perplexity (Salazar et al., 2020) --
                     mask each position one at a time and re-run the model.
                     More accurate, but O(L) forward passes per sequence.
        """
        torch = self.torch
        with torch.no_grad():
            enc = self.tokenizer(protein_seq, return_tensors="pt").to(self.device)
            input_ids = enc["input_ids"]
            seq_len = input_ids.shape[1]
            # positions 0 and -1 are BOS/EOS special tokens for ESM tokenizers
            residue_positions = list(range(1, seq_len - 1))

            if fast:
                logits = self.model(**enc).logits[0]  # (seq_len, vocab)
                log_probs = torch.log_softmax(logits, dim=-1)
                nll_total = 0.0
                for pos in residue_positions:
                    true_id = input_ids[0, pos].item()
                    nll_total += -log_probs[pos, true_id].item()
                mean_nll = nll_total / len(residue_positions)
            else:
                mask_id = self.tokenizer.mask_token_id
                nll_total = 0.0
                for pos in residue_positions:
                    masked = input_ids.clone()
                    true_id = masked[0, pos].item()
                    masked[0, pos] = mask_id
                    logits = self.model(input_ids=masked, attention_mask=enc["attention_mask"]).logits[0]
                    log_probs = torch.log_softmax(logits, dim=-1)
                    nll_total += -log_probs[pos, true_id].item()
                mean_nll = nll_total / len(residue_positions)

        return math.exp(mean_nll)


def score_orfs_with_esm2(
    orfs: List[ORF], model_name: str = "facebook/esm2_t6_8M_UR50D", fast: bool = True
) -> Optional[str]:
    """Populate orf.esm2_perplexity in place. Returns an error message string
    if scoring could not be performed at all, else None."""
    if not orfs:
        return None
    try:
        scorer = ESM2Scorer(model_name=model_name)
    except RuntimeError as e:
        return str(e)
    except Exception as e:  # model download / network failure, etc.
        return f"Could not load ESM-2 model '{model_name}': {e}"

    for orf in orfs:
        try:
            orf.esm2_perplexity = scorer.perplexity(orf.protein, fast=fast)
        except Exception as e:
            orf.esm2_perplexity = None
    return None


# ---------------------------------------------------------------------------
# 6. Report generation
# ---------------------------------------------------------------------------

def build_report(
    input_name: str,
    seq_len: int,
    min_aa: int,
    orfs: List[ORF],
    esm2_error: Optional[str],
    esm2_model: str,
) -> str:
    orf_density_per_kb = (len(orfs) / seq_len) * 1000 if seq_len else 0.0
    orf_to_bp_ratio = len(orfs) / seq_len if seq_len else 0.0

    lines: List[str] = []
    lines.append(f"# ORF Analysis Report")
    lines.append("")
    lines.append(f"- **Input**: {input_name}")
    lines.append(f"- **Number of bases**: {seq_len} bp")
    lines.append(f"- **Minimum ORF length filter**: > {min_aa} amino acids")
    lines.append(f"- **ORFs found (ATG...stop, all 6 frames)**: {len(orfs)}")
    lines.append(f"- **ORF-to-bp ratio**: {orf_to_bp_ratio:.6f} ORFs/bp "
                 f"({orf_density_per_kb:.3f} ORFs per kb)")
    lines.append(f"- **ESM-2 model**: {esm2_model}" + (" (scoring failed, see note below)" if esm2_error else ""))
    lines.append(f"- **Codon usage reference**: E. coli W3110, Kazusa Codon Usage Database "
                 f"(4332 CDS's, 1,372,057 codons)")
    lines.append("")

    if esm2_error:
        lines.append("> **Note on ESM-2 scoring:** " + esm2_error)
        lines.append(">")
        lines.append("> Perplexity columns below are reported as N/A. Install "
                      "`torch` and `transformers` (with internet access to "
                      "download model weights from Hugging Face) to enable this step.")
        lines.append("")

    if not orfs:
        lines.append(f"No ORFs longer than {min_aa} amino acids were found.")
        return "\n".join(lines)

    lines.append("## Summary Table")
    lines.append("")
    lines.append("| # | Strand | Frame | Start | End | Length (aa) | ESM-2 Perplexity | E. coli Codon Match (%) |")
    lines.append("|---|--------|-------|-------|-----|--------------|-------------------|---------------------------|")
    for idx, orf in enumerate(orfs, start=1):
        ppl_str = f"{orf.esm2_perplexity:.2f}" if orf.esm2_perplexity is not None else "N/A"
        cai_str = f"{orf.cai_percent:.1f}%" if orf.cai_percent is not None else "N/A"
        lines.append(
            f"| {idx} | {orf.strand} | {orf.frame} | {orf.start_nt} | {orf.end_nt} | "
            f"{orf.length_aa} | {ppl_str} | {cai_str} |"
        )
    lines.append("")

    lines.append("## ORF Details")
    lines.append("")
    for idx, orf in enumerate(orfs, start=1):
        lines.append(f"### ORF {idx} (strand {orf.strand}, frame {orf.frame}, "
                      f"{orf.start_nt}-{orf.end_nt}, {orf.length_aa} aa)")
        lines.append("")
        lines.append("**Protein sequence:**")
        lines.append("```")
        lines.append(orf.protein)
        lines.append("```")
        ppl_str = f"{orf.esm2_perplexity:.2f}" if orf.esm2_perplexity is not None else "N/A"
        cai_str = f"{orf.cai_percent:.1f}%" if orf.cai_percent is not None else "N/A"
        lines.append(f"- ESM-2 pseudo-perplexity: **{ppl_str}** (lower = more protein-like per ESM-2)")
        lines.append(f"- E. coli codon usage match (CAI-based): **{cai_str}**")
        lines.append("")

    lines.append("---")
    lines.append("### Methodology notes")
    lines.append("")
    lines.append("- ORFs are defined as ATG-to-stop-codon in one of the 6 reading frames, "
                  "matching NCBI ORFfinder's default 'ATG only' start-codon setting.")
    lines.append("- 'ESM-2 pseudo-perplexity' is derived from a masked protein-language model "
                  "(Lin et al. 2022, *Language models of protein sequences at the scale of "
                  "evolution enable accurate structure prediction*); it reflects how likely "
                  "the model considers each residue given its context, not experimental stability or function.")
    lines.append("- 'E. coli codon usage match' is the Codon Adaptation Index (Sharp & Li, 1987) "
                  "of the ORF's codons relative to the E. coli W3110 reference table, expressed "
                  "as a percentage (100% = every codon is E. coli's most-preferred synonym for that amino acid).")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 7. CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Find ORFs (>80 aa) in all 6 reading frames, score them with ESM-2, "
                    "and compare codon usage against E. coli."
    )
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--input", help="Path to a FASTA or plain-text file containing the DNA sequence.")
    src.add_argument("--seq", help="Raw DNA sequence provided directly on the command line.")
    parser.add_argument("--min-aa", type=int, default=80, help="Minimum ORF length in amino acids (default: 80).")
    parser.add_argument("--out", default="orf_report.md", help="Path to write the Markdown report (default: orf_report.md).")
    parser.add_argument("--esm-model", default="facebook/esm2_t6_8M_UR50D",
                         help="Hugging Face model name for ESM-2 (default: facebook/esm2_t6_8M_UR50D, the smallest/fastest). "
                              "Larger options: facebook/esm2_t12_35M_UR50D, facebook/esm2_t30_150M_UR50D, "
                              "facebook/esm2_t33_650M_UR50D.")
    parser.add_argument("--no-esm2", action="store_true", help="Skip ESM-2 scoring entirely.")
    parser.add_argument("--accurate-perplexity", action="store_true",
                         help="Use true masked-residue perplexity (slower: O(length) forward passes per ORF) "
                              "instead of the fast single-pass approximation.")
    args = parser.parse_args()

    source = args.input if args.input else args.seq
    seq = read_sequence(source)

    orfs = find_all_orfs(seq, min_aa=args.min_aa)
    for orf in orfs:
        orf.cai_percent = codon_adaptation_index(orf.nt_seq)

    esm2_error = None
    if not args.no_esm2:
        esm2_error = score_orfs_with_esm2(orfs, model_name=args.esm_model, fast=not args.accurate_perplexity)

    report = build_report(
        input_name=args.input if args.input else "(raw sequence input)",
        seq_len=len(seq),
        min_aa=args.min_aa,
        orfs=orfs,
        esm2_error=esm2_error,
        esm2_model=args.esm_model,
    )

    with open(args.out, "w") as fh:
        fh.write(report)

    density_per_kb = (len(orfs) / len(seq)) * 1000 if len(seq) else 0.0
    print(f"Sequence length: {len(seq)} bp")
    print(f"Found {len(orfs)} ORF(s) > {args.min_aa} aa.")
    print(f"ORF-to-bp ratio: {len(orfs)/len(seq):.6f} ORFs/bp ({density_per_kb:.3f} ORFs per kb)")
    print(f"Report written to {args.out}")


if __name__ == "__main__":
    main()
