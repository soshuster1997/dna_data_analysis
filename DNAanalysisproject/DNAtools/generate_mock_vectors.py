#!/usr/bin/env python3
"""
generate_mock_vectors.py
==========================
Generates mock plasmid-like DNA sequences (300-3000 bp) that mimic real vector
architecture by mixing in common backbone elements -- a resistance marker
ORF, an origin-of-replication-like region, and a promoter/operator/MCS
cloning cassette -- interspersed with random filler DNA, and writes them out
in the same {"data": [...]} bulk JSON format used by addgene_to_fasta.py.

Element provenance (please read before treating this as "real" biology):

  - AmpR (bla / TEM-1 beta-lactamase): the amino acid sequence is the REAL,
    published TEM-1 protein (UniProt P62593, 286 aa) -- the DNA is a
    back-translation of that real protein using E. coli codon-usage-weighted
    codon choices, so it is a genuinely gene-like, high-CAI ORF, but it is a
    synthetic nucleotide sequence, not the literal natural bla gene sequence.
  - "KanR-like" marker: a 264-aa placeholder protein (chosen to match the
    real APH(3')-II / NPTII kanamycin-resistance protein's length) with a
    RANDOM amino acid sequence, back-translated the same way. This gives a
    second gene-like/high-CAI ORF for variety, but it is NOT the real KanR
    sequence -- treat it as illustrative only.
  - T7 promoter (TAATACGACTCACTATAGGG) and the classic lac operator
    (AATTGTGAGCGGATAACAATT) are real, well-established short regulatory
    sequences.
  - The multiple cloning site (MCS) is a real concatenation of common unique
    restriction sites (EcoRI, BamHI, HindIII, XhoI, NotI, SalI).
  - The "ori-like" region is a SYNTHETIC random-sequence placeholder sized
    like a typical ColE1/pUC origin (~600 bp) -- it does NOT reproduce the
    actual ColE1 origin sequence, it just occupies a realistic amount of
    space so sequences have plausible overall composition/length.

In short: this is mock/test data for exercising the ORF-finding, CAI, and
ESM-2 pipeline -- not a substitute for real plasmid sequences.

Usage
-----
    python generate_mock_vectors.py --n 500 --min-len 300 --max-len 3000 --out mock_vectors_backbone.json
"""

from __future__ import annotations

import argparse
import json
import random
from typing import Dict, List, Tuple

from orf_analyzer import ECOLI_CODON_FREQ, STANDARD_CODON_TABLE

# ---------------------------------------------------------------------------
# Backbone element definitions
# ---------------------------------------------------------------------------

# Real TEM-1 beta-lactamase (AmpR/bla) protein sequence, UniProt P62593 (286 aa)
AMPR_PROTEIN = (
    "MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRP"
    "EERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVREL"
    "CSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTM"
    "PAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGS"
    "RGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW"
)

T7_PROMOTER = "TAATACGACTCACTATAGGG"
LAC_OPERATOR = "AATTGTGAGCGGATAACAATT"
MCS = "GAATTCGGATCCAAGCTTCTCGAGGCGGCCGCGTCGAC"  # EcoRI, BamHI, HindIII, XhoI, NotI, SalI
ORI_LIKE_LENGTH = 600  # typical ColE1/pUC-ori-scale footprint (synthetic filler, see docstring)
KANR_LIKE_LENGTH_AA = 264  # matches real APH(3')-II length, but sequence itself is random/synthetic


def _aa_to_codons_table() -> Dict[str, List[str]]:
    table: Dict[str, List[str]] = {}
    for codon, aa in STANDARD_CODON_TABLE.items():
        table.setdefault(aa, []).append(codon)
    return table


AA_TO_CODONS = _aa_to_codons_table()
STOP_CODONS = [c for c, aa in STANDARD_CODON_TABLE.items() if aa == "*"]


def backtranslate(protein: str, rng: random.Random) -> str:
    """Back-translate a protein to DNA using E. coli codon-usage-weighted
    random choices among synonymous codons (biases toward realistic, high-CAI
    codon usage without being deterministic)."""
    codons = []
    for aa in protein:
        choices = AA_TO_CODONS[aa]
        weights = [ECOLI_CODON_FREQ[c] for c in choices]
        codons.append(rng.choices(choices, weights=weights, k=1)[0])
    stop_weights = [ECOLI_CODON_FREQ[c] for c in STOP_CODONS]
    codons.append(rng.choices(STOP_CODONS, weights=stop_weights, k=1)[0])
    return "".join(codons)


def random_protein(length: int, rng: random.Random) -> str:
    amino_acids = [aa for aa in AA_TO_CODONS if aa != "*"]
    return "".join(rng.choice(amino_acids) for _ in range(length))


def random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


# ---------------------------------------------------------------------------
# Sequence assembly
# ---------------------------------------------------------------------------

def build_elements(rng: random.Random) -> Dict[str, str]:
    """Build one fresh set of backbone element sequences. AmpR is fixed (real
    protein); the KanR-like decoy and ori-like filler are freshly randomized
    per call so not every 'KanR' or 'ori' region in the dataset is identical."""
    return {
        "AmpR_bla": backtranslate(AMPR_PROTEIN, rng),
        "KanR_like": backtranslate(random_protein(KANR_LIKE_LENGTH_AA, rng), rng),
        "T7_promoter": T7_PROMOTER,
        "lac_operator": LAC_OPERATOR,
        "MCS": MCS,
        "ori_like": random_dna(ORI_LIKE_LENGTH, rng),
    }


def choose_layout(target_length: int, rng: random.Random) -> List[str]:
    """Decide which named elements to include for a plasmid of a given
    target length, prioritizing a resistance marker, then ori, then the
    promoter/operator/MCS cassette, dropping anything that won't fit."""
    budget = int(target_length * 0.9)  # leave room for filler DNA
    layout: List[str] = []

    marker_roll = rng.random()
    if marker_roll < 0.45:
        candidate = "AmpR_bla"
    elif marker_roll < 0.85:
        candidate = "KanR_like"
    else:
        candidate = None  # ~15% of sequences have no resistance marker at all

    element_lengths = {
        "AmpR_bla": len(AMPR_PROTEIN) * 3 + 3,
        "KanR_like": KANR_LIKE_LENGTH_AA * 3 + 3,
        "ori_like": ORI_LIKE_LENGTH,
        "promoter_cassette": len(T7_PROMOTER) + len(LAC_OPERATOR) + len(MCS) + 6,  # +spacers
    }

    used = 0
    if candidate and used + element_lengths[candidate] <= budget:
        layout.append(candidate)
        used += element_lengths[candidate]

    if rng.random() < 0.5 and used + element_lengths["ori_like"] <= budget:
        layout.append("ori_like")
        used += element_lengths["ori_like"]

    if rng.random() < 0.5 and used + element_lengths["promoter_cassette"] <= budget:
        layout.append("promoter_cassette")
        used += element_lengths["promoter_cassette"]

    rng.shuffle(layout)
    return layout


def assemble_sequence(target_length: int, rng: random.Random) -> Tuple[str, List[Dict]]:
    elements = build_elements(rng)
    layout = choose_layout(target_length, rng)

    # expand "promoter_cassette" into its 3 constituent pieces with small random spacers
    pieces: List[Tuple[str, str]] = []  # (label, sequence)
    for name in layout:
        if name == "promoter_cassette":
            pieces.append(("T7_promoter", elements["T7_promoter"]))
            pieces.append(("lac_operator", elements["lac_operator"]))
            pieces.append(("MCS", elements["MCS"]))
        else:
            pieces.append((name, elements[name]))

    fixed_len = sum(len(seq) for _, seq in pieces)
    filler_total = max(target_length - fixed_len, 0)

    # split filler into len(pieces)+1 random-sized chunks (some can be 0)
    n_chunks = len(pieces) + 1
    if n_chunks == 1:
        chunk_sizes = [filler_total]
    else:
        cut_points = sorted(rng.randint(0, filler_total) for _ in range(n_chunks - 1))
        chunk_sizes = [b - a for a, b in zip([0] + cut_points, cut_points + [filler_total])]

    sequence_parts: List[str] = []
    features: List[Dict] = []
    pos = 0

    def add_filler(size: int):
        nonlocal pos
        if size > 0:
            seq = random_dna(size, rng)
            sequence_parts.append(seq)
            pos += size

    add_filler(chunk_sizes[0])
    for (label, seq), chunk_size in zip(pieces, chunk_sizes[1:]):
        start = pos + 1  # 1-based
        sequence_parts.append(seq)
        pos += len(seq)
        features.append({"label": label, "start": start, "end": pos})
        add_filler(chunk_size)

    full_seq = "".join(sequence_parts)
    # target_length is approximate due to integer rounding in filler split; that's fine for mock data
    return full_seq, features


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate mock plasmid-like DNA sequences with common vector backbone elements mixed in.")
    parser.add_argument("--n", type=int, default=500, help="Number of sequences to generate (default: 500).")
    parser.add_argument("--min-len", type=int, default=300, help="Minimum sequence length in bp (default: 300).")
    parser.add_argument("--max-len", type=int, default=3000, help="Maximum sequence length in bp (default: 3000).")
    parser.add_argument("--seed", type=int, default=2026, help="Random seed for reproducibility (default: 2026).")
    parser.add_argument("--out", default="mock_vectors_backbone.json", help="Output JSON path.")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    records = []
    vector_types = ["Mammalian Expression", "Bacterial Expression", "Lentiviral", "AAV", "CRISPR"]

    for i in range(args.n):
        target_length = rng.randint(args.min_len, args.max_len)
        seq, features = assemble_sequence(target_length, rng)
        records.append({
            "id": 200000 + i,
            "name": f"pMockBackbone{i:04d}",
            "full_sequence": seq,
            "vector_type": rng.choice(vector_types),
            "features": features,  # ground-truth annotation of which elements were inserted where
        })

    with open(args.out, "w") as fh:
        json.dump({"data": records}, fh)

    lengths = [len(r["full_sequence"]) for r in records]
    with_marker = sum(1 for r in records if any(f["label"] in ("AmpR_bla", "KanR_like") for f in r["features"]))
    print(f"Wrote {len(records)} records to {args.out}")
    print(f"Length range: {min(lengths)}-{max(lengths)} bp")
    print(f"Records containing a resistance marker ORF: {with_marker} ({with_marker/len(records)*100:.0f}%)")


if __name__ == "__main__":
    main()
