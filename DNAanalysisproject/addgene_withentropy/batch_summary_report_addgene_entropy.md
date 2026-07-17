# Batch ORF Analysis Summary (no ESM-2)

- **Input directory**: addgene_fasta
- **FASTA files processed**: 483
- **Total bases across all files**: 3,151,658 bp
- **Start codons searched**: ATG, GTG, TTG, CTG (non-ATG codons treated as alternative/non-traditional starts)
- **Minimum ORF length filter**: > 80 amino acids
- **k-mer size for repeat detection**: 6 bp
- **Total ORFs found**: 8207
- **ORFs with a non-traditional (non-ATG) start codon**: 5892 (71.8% of all ORFs)
- **Files with zero qualifying ORFs**: 11

## Descriptive Statistics

- **ORFs per file**: n=483, mean=16.99, median=16.00, stdev=10.76, min=0.00, max=120.00
- **ORF density (per kb)**: n=483, mean=2.572, median=2.585, stdev=0.851, min=0.000, max=9.685
- **CAI codon match (%), across all 8207 ORFs**: n=8207, mean=55.6, median=55.6, stdev=6.1, min=33.5, max=96.6
- **Whole-sequence Shannon entropy (bits, 0-2)**: n=483, mean=1.993, median=1.997, stdev=0.012, min=1.859, max=2.000
- **Whole-sequence repeat fraction (k=6)**: n=483, mean=0.499, median=0.532, stdev=0.163, min=0.070, max=0.886
- **Per-ORF Shannon entropy (bits, 0-2)**: n=8207, mean=1.961, median=1.973, stdev=0.046, min=1.502, max=2.000
- **Per-ORF repeat fraction (k=6)**: n=8207, mean=0.117, median=0.084, stdev=0.106, min=0.004, max=0.900

## Per-File Results (first 50)

(Full detail in `per_file_summary.csv`; individual ORFs in `all_orfs.csv`.)

| File | Length (bp) | # ORFs | ORFs/kb | Mean CAI (%) | Non-ATG starts | Entropy | Repeat frac |
|------|-------------|--------|---------|---------------|-----------------|---------|--------------|
| 100643_YCe1917_HC_Kan_conn_D-F.fasta | 1889 | 4 | 2.12 | 58.8 | 2 | 1.998 | 0.225 |
| 100749_FRB-mCherry-FTH1.fasta | 5565 | 12 | 2.16 | 56.4 | 9 | 1.996 | 0.490 |
| 102618_p2lox-cAMPr_R739A.fasta | 6043 | 14 | 2.32 | 53.0 | 10 | 1.998 | 0.512 |
| 102747_A435S_EPHA3_pcDNA3.1.fasta | 8539 | 18 | 2.11 | 54.2 | 12 | 2.000 | 0.601 |
| 103001_pET22b_PfSSB_W166C.fasta | 5990 | 18 | 3.01 | 61.8 | 12 | 1.999 | 0.509 |
| 103596_LSB-hsa-miR-513c-5p.fasta | 11043 | 40 | 3.62 | 54.6 | 27 | 1.995 | 0.699 |
| 103615_LSB-hsa-miR-525-5p.fasta | 11039 | 42 | 3.80 | 54.2 | 28 | 1.995 | 0.699 |
| 103663_LSB-hsa-miR-582-5p.fasta | 11047 | 40 | 3.62 | 54.6 | 27 | 1.995 | 0.699 |
| 103708_LSB-hsa-miR-652-5p.fasta | 11055 | 39 | 3.53 | 54.6 | 27 | 1.994 | 0.699 |
| 104388_PB-eDIO-EYFP.fasta | 9600 | 30 | 3.12 | 55.1 | 21 | 1.994 | 0.639 |
| 105118_pFRETgc-2in1-CN.fasta | 13356 | 34 | 2.55 | 58.0 | 24 | 1.999 | 0.721 |
| 105964_hsp70_Arl13b-mKate2.fasta | 7776 | 19 | 2.44 | 54.9 | 13 | 1.988 | 0.590 |
| 106285_pJL1_Aquamarine.fasta | 2639 | 8 | 3.03 | 57.6 | 7 | 1.996 | 0.304 |
| 107232_pLT3REVIR_MARK3_1122.fasta | 10284 | 29 | 2.82 | 53.5 | 21 | 1.999 | 0.656 |
| 108100_LentiV_Cas9_puro.fasta | 11203 | 29 | 2.59 | 53.8 | 23 | 1.993 | 0.691 |
| 108743_pBW2603_pCAG-PV1-iCre-N256-L1-ABI-NLS-BGHpA.fasta | 6988 | 25 | 3.58 | 54.0 | 15 | 1.993 | 0.565 |
| 108759_pBW2280_pCAG-PV1-FlpO-N374-L1-FRB-NLS-BGHpA.fasta | 6703 | 22 | 3.28 | 54.7 | 13 | 1.992 | 0.562 |
| 109584_p3E-CDC42-DN.fasta | 3214 | 7 | 2.18 | 56.4 | 7 | 1.998 | 0.341 |
| 109636_actc1b-mKate2-rab3ab.fasta | 9604 | 16 | 1.67 | 55.2 | 10 | 1.966 | 0.656 |
| 109886_pHH0103_DOK6_IRS.fasta | 7044 | 18 | 2.56 | 60.2 | 12 | 2.000 | 0.550 |
| 109965_pHH0103_CBX1_Chromo-shadow.fasta | 6921 | 18 | 2.60 | 60.0 | 12 | 2.000 | 0.549 |
| 110402_psiCHECK2_3_UTR_KGA.fasta | 8803 | 16 | 1.82 | 54.3 | 13 | 1.999 | 0.618 |
| 110666_pCHKU33.1-2.2.fasta | 14222 | 37 | 2.60 | 54.9 | 28 | 2.000 | 0.729 |
| 110823_pB-CAGGS-dCas9.fasta | 11681 | 40 | 3.42 | 55.2 | 33 | 1.996 | 0.686 |
| 110837_pLenti-Cas9-P2A-Puro.fasta | 11203 | 27 | 2.41 | 53.9 | 23 | 1.993 | 0.691 |
| 111123_pCR1280.fasta | 8914 | 25 | 2.80 | 53.8 | 18 | 1.998 | 0.614 |
| 111730_pBacMam2-DiEx-LIC-C-flag_huntingtin_full-length_Q109.fasta | 16038 | 37 | 2.31 | 54.8 | 28 | 1.999 | 0.761 |
| 111975_Hsh155_N4A_3xFLAG.fasta | 8211 | 11 | 1.34 | 54.1 | 7 | 1.982 | 0.601 |
| 112984_mei-S332_5.6Kb_genomic_DNA_in_Casper4.fasta | 13515 | 30 | 2.22 | 58.5 | 22 | 1.998 | 0.717 |
| 113125_EcRfaHKOW.fasta | 6679 | 22 | 3.29 | 60.1 | 18 | 1.998 | 0.539 |
| 114942_AT1G51810_M2_pECIA2.fasta | 5411 | 14 | 2.59 | 54.0 | 12 | 1.998 | 0.462 |
| 115855_pAPM-D4_miR30-AGO2_ts3.fasta | 6319 | 15 | 2.37 | 54.2 | 12 | 1.997 | 0.532 |
| 116035_pVL90_-_VND7.fasta | 3663 | 8 | 2.18 | 55.2 | 5 | 1.999 | 0.356 |
| 116283_pHAGE-EGFR-P589L.fasta | 12122 | 35 | 2.89 | 54.2 | 25 | 1.998 | 0.697 |
| 116429_pHAGE-MAP2K1-S222A.fasta | 9694 | 21 | 2.17 | 53.8 | 14 | 1.998 | 0.642 |
| 116660_pHAGE-RAF1-S257L.fasta | 10459 | 24 | 2.29 | 53.5 | 15 | 1.999 | 0.659 |
| 117104_SLAH2.fasta | 5564 | 11 | 1.98 | 56.4 | 9 | 1.997 | 0.482 |
| 11752_pCMV5B-Flag-Smurf1_wt.fasta | 6847 | 17 | 2.48 | 55.0 | 9 | 2.000 | 0.543 |
| 118168_lentiSAMv2_EGFP-guide3.fasta | 13501 | 32 | 2.37 | 53.8 | 28 | 1.996 | 0.729 |
| 118678_pDONR_221-AIP2_NON-STOP_.fasta | 3477 | 9 | 2.59 | 55.2 | 9 | 1.998 | 0.359 |
| 118963_puc19-pCas13a.fasta | 9089 | 21 | 2.31 | 55.8 | 12 | 1.991 | 0.632 |
| 119775_pAtU6-sgRNA.fasta | 4629 | 13 | 2.81 | 55.7 | 7 | 1.999 | 0.423 |
| 121054_pJW1583.fasta | 10661 | 32 | 3.00 | 56.3 | 22 | 1.998 | 0.662 |
| 121066_pUHD-AID.fasta | 5039 | 16 | 3.18 | 54.7 | 8 | 1.995 | 0.463 |
| 122168_p8xEGFP-N1.fasta | 9848 | 16 | 1.62 | 54.0 | 14 | 1.975 | 0.737 |
| 122441_pHR-FUSN-miRFP670-Cry2WT.fasta | 12184 | 29 | 2.38 | 56.9 | 17 | 1.999 | 0.695 |
| 122589_Human_Inositol_1_4_5-Trisphosphate_3-Kinase_B.fasta | 8457 | 27 | 3.19 | 53.0 | 20 | 1.992 | 0.607 |
| 123152_pXWP11-gfp2.fasta | 3150 | 6 | 1.90 | 56.5 | 5 | 1.990 | 0.345 |
| 123197_pSB161_-_pL2_pSB90_2x35S_GUS_tMAS.fasta | 9010 | 28 | 3.11 | 58.7 | 20 | 1.997 | 0.627 |
| 123558_pMK-aadA.fasta | 3122 | 13 | 4.16 | 60.3 | 10 | 1.997 | 0.360 |

*(showing first 50 of 483 files -- see per_file_summary.csv for all)*