# DNA ORF Analysis Toolkit

A set of standalone Python scripts for finding open reading frames (ORFs) in DNA, scoring them against real biological references (E. coli codon usage, ESM-2 protein language model), and comparing real plasmid DNA against synthetic/random DNA. Built up incrementally to answer one core question: **can a handful of cheap, well-understood metrics tell real biological DNA apart from DNA that only looks like DNA?**

No script requires any other script to run standalone, except where noted (`orf_analyzer.py` is a shared dependency for several tools, and `dna_binary_converter.py` is a dependency of the Wikipedia-to-DNA scripts).

---

## Contents

| Script | What it does |
|---|---|
| `orf_analyzer.py` | Core library + standalone tool: finds ORFs in a single DNA sequence, scores CAI and (optionally) ESM-2 perplexity, writes a Markdown report. |
| `batch_orf_pipeline.py` | Runs `orf_analyzer`'s logic (ORFs + CAI + ESM-2) across a whole directory of FASTA files at once, loading the ESM-2 model only once for the whole batch. |
| `batch_orf_pipeline_no_esm2.py` | Same idea, but drops ESM-2 entirely (no GPU/model download needed) and adds alternative start-codon detection, Shannon entropy, and k-mer repeat-content checks. |
| `addgene_to_fasta.py` | Pulls real plasmid sequences from Addgene's Developers Portal API (or a local bulk JSON download) and writes them out as individual FASTA files. |
| `wikipedia_dump_to_dna.py` | Converts real Wikipedia article text into DNA (via binary) using a downloaded Wikipedia XML dump. No API rate limits since it's just parsing a local file. |
| `wikipedia_to_dna.py` | Same idea, but pulls articles live from Wikipedia's API instead of a downloaded dump. Simpler to run, but subject to Wikipedia's rate limits. |
| `dna_binary_converter.py` | Small utility: converts between binary and DNA using a fixed 2-bit code (A=00, C=01, G=10, T=11). Used internally by the Wikipedia-to-DNA scripts. |
| `visualize_orfs.py` | Plots histograms (ESM-2 perplexity and CAI%) from an `orf_analyzer.py` report, or by re-running the analysis directly. |
| `generate_mock_vectors.py` | Generates synthetic plasmid-like test sequences with real/realistic backbone elements (AmpR, promoters, MCS, etc.) mixed into random filler DNA. Useful for testing the pipeline without needing real data. |

---

## Setup

Only `batch_orf_pipeline.py` and `orf_analyzer.py`'s ESM-2 scoring need extra packages:

```bash
pip install torch transformers          # only needed for ESM-2 scoring
pip install requests                    # only needed for addgene_to_fasta.py --api mode
```

Everything else (`batch_orf_pipeline_no_esm2.py`, `addgene_to_fasta.py` with `--bulk-json`, `dna_binary_converter.py`, `wikipedia_dump_to_dna.py`) runs on the Python standard library alone.

Keep all scripts in the same folder. Several of them import from each other:
- `batch_orf_pipeline.py`, `batch_orf_pipeline_no_esm2.py`, and `visualize_orfs.py` import from `orf_analyzer.py`.
- `wikipedia_to_dna.py` and `wikipedia_dump_to_dna.py` import from `dna_binary_converter.py`.

---

## Script-by-script reference

### `orf_analyzer.py`
Finds every ORF (ATG to in-frame stop, all 6 reading frames, >80 aa by default) in a single DNA sequence, scores each one's Codon Adaptation Index (CAI) against a real E. coli codon usage table, optionally scores each one with ESM-2, and writes a Markdown report.

```bash
python orf_analyzer.py --input sequence.fasta --out report.md
python orf_analyzer.py --seq ATGCGTACG... --out report.md
python orf_analyzer.py --input seq.fasta --no-esm2                        # skip ESM-2, much faster
python orf_analyzer.py --input seq.fasta --accurate-perplexity            # slower, correct masked scoring (see note below)
python orf_analyzer.py --input seq.fasta --esm-model facebook/esm2_t30_150M_UR50D
```

Key flags: `--min-aa` (default 80), `--esm-model` (default `facebook/esm2_t6_8M_UR50D`, the smallest/fastest), `--no-esm2`, `--accurate-perplexity`.

**Important:** always use `--accurate-perplexity` for real analysis. The default fast mode reads each residue's probability off a single unmasked forward pass; because ESM-2 is bidirectional, a residue's own embedding leaks into its own prediction, which artificially compresses perplexity into a narrow, uninformative ~1.0-2.0 range regardless of input. `--accurate-perplexity` does true masked-residue scoring (one mask + forward pass per residue, so slower, but correct).

### `batch_orf_pipeline.py`
Same analysis as `orf_analyzer.py`, but across every FASTA file in a directory, loading the ESM-2 model once and reusing it (much faster than calling `orf_analyzer.py` once per file, which would reload the model every time).

```bash
python batch_orf_pipeline.py --fasta-dir my_fasta_folder --out-dir batch_results
python batch_orf_pipeline.py --fasta-dir my_fasta_folder --out-dir batch_results --accurate-perplexity
python batch_orf_pipeline.py --fasta-dir my_fasta_folder --out-dir batch_results --no-esm2   # skip ESM-2 entirely
```

Outputs in `--out-dir`:
- `batch_summary_report.md` — descriptive stats (mean/median/stdev/min/max) across the whole batch for ORF count, ORF density, CAI%, and perplexity.
- `per_file_summary.csv` — one row per FASTA file.
- `all_orfs.csv` — one row per individual ORF found, across all files.

Use `--progress-every N` to control how often progress prints during a large run (default: every 25 files).

If you skip `--accurate-perplexity`, both the console output and the report will print an explicit warning about the fast-mode bias described above.

### `batch_orf_pipeline_no_esm2.py`
A variant of the batch pipeline that drops ESM-2 completely and instead:
- Searches for ATG **and** the common bacterial alternative start codons GTG, TTG, and CTG by default (flagging which ORFs used a non-ATG start). Use `--atg-only` to restrict to ATG only, or `--start-codons ATG,GTG,TTG,CTG,ATT` to customize.
- Computes Shannon entropy (0-2 bits, based on A/C/G/T frequency) for each ORF and each whole sequence.
- Computes a k-mer repeat fraction (default k=6): what proportion of a sequence's k-mers are repeats of one seen earlier in the same sequence, a fast repetitiveness proxy.

```bash
python batch_orf_pipeline_no_esm2.py --fasta-dir my_fasta_folder --out-dir batch_results
python batch_orf_pipeline_no_esm2.py --fasta-dir my_fasta_folder --out-dir batch_results --atg-only
python batch_orf_pipeline_no_esm2.py --fasta-dir my_fasta_folder --out-dir batch_results --kmer-size 8
```

Same output structure as `batch_orf_pipeline.py` (`batch_summary_report.md`, `per_file_summary.csv`, `all_orfs.csv`, `--progress-every` works the same way), with extra columns for start codon, non-standard-start flag, entropy, and repeat fraction/count. No `torch`/`transformers` needed at all, so this one runs almost instantly even on a laptop CPU.

**Caveat on Shannon entropy:** as implemented, it's an order-0 measure of base composition only. It's blind to sequence order, so a repetitive sequence like `AAAACCCCGGGGTTTT` scores a perfect 2.0 (even base counts) despite being obviously structured. It answers "is the base composition balanced?", not "does this look random?" The k-mer repeat fraction is the metric that actually captures repetitiveness.

### `addgene_to_fasta.py`
Converts Addgene plasmid records into individual FASTA files, either from a local bulk JSON download or by calling Addgene's Developers Portal API directly (requires an approved API token).

```bash
# From a bulk JSON file downloaded from the Developers Portal
python addgene_to_fasta.py --bulk-json addgene_bulk.json --out-dir addgene_fasta --limit 1000

# From the live API (fetches everything in one authenticated request, then samples)
python addgene_to_fasta.py --api --token YOUR_TOKEN --out-dir addgene_fasta --limit 1000

# Inspect the raw record structure first (recommended before a big run)
python addgene_to_fasta.py --bulk-json addgene_bulk.json --inspect
```

By default, `--limit` draws a **random** sample from however many records are available (Addgene's full dataset can be 100,000+ plasmids) rather than just taking the first N. Use `--seed 42` for a reproducible sample, or `--sample first` to revert to taking them in whatever order Addgene returned them. Writes a `manifest.csv` alongside the FASTA files recording each plasmid's ID, name, which sequence field was used (Addgene's own vs. depositor-submitted, full vs. partial), resistance markers, and vector type.

### `wikipedia_dump_to_dna.py`
Generates synthetic "DNA" by taking real Wikipedia article text, converting it to binary, then to nucleotides. This version reads a Wikipedia XML database dump you've downloaded yourself, so there's no API rate limit to worry about (it's just parsing a static file).

```bash
# 1. Download a dump first (Simple English Wikipedia recommended: much smaller than full English Wikipedia)
wget https://dumps.wikimedia.org/simplewiki/latest/simplewiki-latest-pages-articles-multistream.xml.bz2

# 2. Run the script directly on the downloaded .bz2 (no need to unzip it yourself)
python wikipedia_dump_to_dna.py --dump-file simplewiki-latest-pages-articles-multistream.xml.bz2 \
    --n 500 --min-bp 300 --max-bp 5000 --out-dir wiki_dna_fasta --seed 42
```

Streams the dump (memory stays flat regardless of file size), filters out redirects/talk pages/stubs, uses reservoir sampling (`--pool-size`, default 3x `--n`) to build an unbiased random pool, and pads short articles with *other* random articles (not repeats of the same text) if needed to hit the target length, avoiding artificial periodicity in the output DNA. Use `--max-scan N` to cap how many pages get scanned before stopping, handy if you point this at the much larger full English Wikipedia dump instead of Simple English Wikipedia and don't want to wait for a full scan.

### `wikipedia_to_dna.py`
Same text-to-DNA conversion, but fetches articles live from Wikipedia's API instead of a downloaded dump. Simpler to get started (no multi-hundred-MB download first), but subject to Wikipedia's rate limits on large runs.

```bash
python wikipedia_to_dna.py --n 500 --min-bp 300 --max-bp 5000 --out-dir wiki_dna_fasta --seed 42
```

Fetches articles in batches (`--batch-size`, default 20 per API call, not one call per article) and automatically backs off and retries on rate-limit (429) responses (`--max-retries-per-batch`, default 6), respecting the `Retry-After` header. If you still hit rate limits frequently, lower `--batch-size` and raise `--delay`, or switch to `wikipedia_dump_to_dna.py`.

### `dna_binary_converter.py`
A small standalone utility. Converts between binary and DNA using a fixed 2-bit code: **A=00, C=01, G=10, T=11**.

```bash
python dna_binary_converter.py --to-dna 0001101100
python dna_binary_converter.py --to-binary ACGT
python dna_binary_converter.py --to-dna 0001101100 --out result.fasta --header my_sequence
```

### `visualize_orfs.py`
Plots histograms of ESM-2 perplexity and CAI% for a set of ORFs, either by parsing an existing `orf_analyzer.py` Markdown report or by re-running the analysis directly on a FASTA file.

```bash
python visualize_orfs.py --report report.md --out-prefix orf_plots
python visualize_orfs.py --input sequence.fasta --out-prefix orf_plots
```

Note: this parses the single-sequence report format from `orf_analyzer.py`, not the batch summary format from `batch_orf_pipeline.py`. For plotting batch results, load `all_orfs.csv` or `per_file_summary.csv` directly with pandas/matplotlib instead. Use `--bins N` to control histogram bin count (default 15).

### `generate_mock_vectors.py`
Generates synthetic plasmid-like sequences for testing the pipeline without needing real data or downloads. Mixes random filler DNA with a real backbone element (the actual TEM-1 beta-lactamase/AmpR protein, back-translated with E. coli-weighted codons), a same-length randomized decoy marker, real regulatory motifs (T7 promoter, lac operator, a multiple cloning site), and a synthetic origin-of-replication-sized filler region.

```bash
python generate_mock_vectors.py --n 500 --min-len 300 --max-len 3000 --out mock_vectors.json
```

Output is bulk JSON in the same `{"plasmids": [...]}`-style shape Addgene's API returns (specifically, a `{"data": [...]}` list with `id`/`name`/`full_sequence` fields), so it can be fed into `addgene_to_fasta.py` for testing before you have real API access.

---

## Example: full pipeline, start to finish

This example builds two comparable datasets, a real one and a synthetic one, and runs the same analysis on both.

```bash
# --- Real DNA: pull plasmids from Addgene ---
python addgene_to_fasta.py --api --token YOUR_TOKEN --out-dir addgene_fasta --limit 500 --seed 42

# --- Synthetic DNA: convert Wikipedia text into DNA ---
wget https://dumps.wikimedia.org/simplewiki/latest/simplewiki-latest-pages-articles-multistream.xml.bz2
python wikipedia_dump_to_dna.py --dump-file simplewiki-latest-pages-articles-multistream.xml.bz2 \
    --n 500 --min-bp 300 --max-bp 5000 --out-dir wiki_dna_fasta --seed 42

# --- Run the full pipeline (ORFs + CAI + ESM-2) on both ---
python batch_orf_pipeline.py --fasta-dir addgene_fasta   --out-dir addgene_results   --accurate-perplexity
python batch_orf_pipeline.py --fasta-dir wiki_dna_fasta  --out-dir wiki_results      --accurate-perplexity

# --- Or, run the faster no-ESM2 pipeline instead (adds entropy/repeat checks, alt start codons) ---
python batch_orf_pipeline_no_esm2.py --fasta-dir addgene_fasta  --out-dir addgene_results_v2
python batch_orf_pipeline_no_esm2.py --fasta-dir wiki_dna_fasta --out-dir wiki_results_v2

# --- Read the results ---
cat addgene_results/batch_summary_report.md
cat wiki_results/batch_summary_report.md
```

From there, `all_orfs.csv` and `per_file_summary.csv` in each output folder can be loaded directly (pandas, R, Excel, whatever you prefer) for further comparison, plotting, or statistical testing between the two datasets.


