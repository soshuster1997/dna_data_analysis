# ORF Analyzer

Finds ORFs across all 6 reading frames (like NCBI's [ORFfinder](https://www.ncbi.nlm.nih.gov/orffinder/)
with the default ATG-only start setting), scores each candidate protein with
an ESM-2 protein language model, compares codon usage to *E. coli*, and
writes it all to a Markdown report.

## Install

```bash
pip install -r requirements.txt   # torch + transformers, only needed for ESM-2 scoring
```

The ORF-finding and codon-usage steps have **no dependencies** beyond the
Python standard library. `torch`/`transformers` are only needed for the
ESM-2 perplexity step, and the first run will download model weights from
Hugging Face (requires internet access).

## Usage

```bash
# From a FASTA file
python orf_analyzer.py --input my_sequence.fasta --out report.md

# From a raw sequence typed/pasted directly
python orf_analyzer.py --seq ATGCGTACG... --out report.md

# Skip ESM-2 (fast, no model download)
python orf_analyzer.py --input my_sequence.fasta --no-esm2

# Use a bigger, more accurate ESM-2 model (slower, bigger download)
python orf_analyzer.py --input my_sequence.fasta --esm-model facebook/esm2_t30_150M_UR50D

# Use true masked pseudo-perplexity instead of the fast approximation
# (O(protein length) forward passes per ORF instead of O(1) -- much slower)
python orf_analyzer.py --input my_sequence.fasta --accurate-perplexity
```

Other options:
- `--min-aa N` — minimum ORF length in amino acids (default: 80, i.e. ORFs
  reported are strictly longer than this).
- `--out path.md` — where to write the report (default: `orf_report.md`).

## What each metric means

**ESM-2 pseudo-perplexity.** ESM-2 (Lin et al., 2022) is a protein language
model trained on hundreds of millions of natural protein sequences. Lower
perplexity means the model finds a sequence more consistent with real
protein sequences it was trained on — a rough, model-based signal of
"protein-likeness," not a measurement of stability, folding, or function.
By default the tool uses the fast single-forward-pass approximation; pass
`--accurate-perplexity` for the slower, more rigorous masked-residue version
(Salazar et al., 2020) if you want tighter numbers, e.g. for a final
shortlist of a few ORFs.

**E. coli codon usage match.** This is the Codon Adaptation Index (CAI;
Sharp & Li, 1987) of each ORF's own codons against a reference *E. coli*
W3110 codon usage table (Kazusa Codon Usage Database, 4,332 CDS's). It's
expressed as a percentage, where 100% means every codon in the ORF is
*E. coli*'s single most-preferred synonymous codon for that amino acid.
This is a standard way to gauge whether a sequence's codon choices already
look "E. coli-like" (useful context if you're thinking about expressing the
protein in *E. coli*), not a measure of protein quality.

## Notes on ORF-finding logic

- ORFs must start with `ATG` and end at the first in-frame stop codon
  (`TAA`, `TAG`, `TGA`) — matching ORFfinder's default settings.
- All 6 frames are searched: 3 forward frames and 3 frames on the reverse
  complement.
- Only ORFs translating to **more than the `--min-aa` amino acids** (default
  80) are reported.
- Nested/overlapping ATGs within the same frame are each considered as
  separate candidate ORFs (the same behavior as ORFfinder).
- Reported start/end coordinates are always given in the coordinates of the
  original input sequence you provided (1-based, inclusive), even for
  reverse-strand ORFs.

## Files

- `orf_analyzer.py` — the tool itself (also importable as a module).
- `requirements.txt` — pip dependencies for the ESM-2 step.
