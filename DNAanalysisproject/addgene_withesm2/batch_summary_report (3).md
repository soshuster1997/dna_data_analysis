# Batch ORF Analysis Summary

- **Input directory**: addgene_fasta
- **FASTA files processed**: 484
- **Total bases across all files**: 3,203,992 bp
- **Total ORFs found (> 80 aa)**: 4064
- **Files with zero qualifying ORFs**: 20
- **ESM-2 model**: facebook/esm2_t6_8M_UR50D

## Descriptive Statistics

- **ORFs per file**: n=484, mean=8.40, median=8.00, stdev=5.60, min=0.00, max=61.00
- **ORF density (ORFs per kb)**: n=484, mean=1.247, median=1.249, stdev=0.467, min=0.000, max=3.501
- **CAI codon match (%), across all 4064 ORFs**: n=4064, mean=56.9, median=56.6, stdev=6.6, min=22.2, max=100.0
- **ESM-2 perplexity, across all 4064 ORFs**: n=4064, mean=13.99, median=14.45, stdev=3.05, min=1.17, max=19.91

## Per-File Results

(Full detail in `per_file_summary.csv`; individual ORFs in `all_orfs.csv`.)

| File | Length (bp) | # ORFs | ORFs/kb | Mean CAI (%) | Mean Perplexity |
|------|-------------|--------|---------|---------------|-------------------|
| 100032_pMT3-RhoAC.fasta | 7053 | 7 | 0.99 | 59.8 | 14.30 |
| 100635_YCe1721_HC_Kan_RFP-p25.fasta | 2935 | 3 | 1.02 | 67.3 | 15.06 |
| 100786_DsRed-CPI-17-T38A.fasta | 6500 | 9 | 1.38 | 56.4 | 13.45 |
| 102864_pLEX307-mFzd1.fasta | 10656 | 13 | 1.22 | 55.2 | 14.40 |
| 103013_S-HBsAg.fasta | 6118 | 8 | 1.31 | 56.1 | 13.90 |
| 103032_Lifeact-7-iRFP670.fasta | 4976 | 6 | 1.21 | 62.0 | 15.33 |
| 103400_LSB-hsa-miR-302a-3p.fasta | 11047 | 17 | 1.54 | 55.9 | 13.91 |
| 103900_PTPN13_PDZ_3.fasta | 5316 | 8 | 1.50 | 59.6 | 13.70 |
| 104466_MBP-hnRNPA2_LC_R191K_R254K.fasta | 6912 | 9 | 1.30 | 60.1 | 14.65 |
| 104517_LI_C-D_GUS.fasta | 4255 | 6 | 1.41 | 64.9 | 15.50 |
| 105168_pTRE-T2-_mIR-IRES-Puromycin.fasta | 4859 | 5 | 1.03 | 55.3 | 9.66 |
| 105414_pTEI041.fasta | 3211 | 8 | 2.49 | 56.3 | 16.06 |
| 105946_HS-aPKC-CAAX-GFP.fasta | 7532 | 7 | 0.93 | 55.4 | 13.82 |
| 106159_pCfB6627.fasta | 6327 | 9 | 1.42 | 57.1 | 14.62 |
| 106706_gh25.fasta | 3180 | 3 | 0.94 | 52.5 | 12.84 |
| 107182_pAAV-CAG-EGFPKir7.1.fasta | 7306 | 10 | 1.37 | 52.4 | 12.80 |
| 10792_1436_pcDNA3_Flag_HA.fasta | 5505 | 6 | 1.09 | 57.0 | 13.79 |
| 108513_pAJM.712.fasta | 2927 | 4 | 1.37 | 57.2 | 13.64 |
| 108785_pBW2585_pCAG-PV1-VCre-N269-L1-ABI-NLS-BGHpA.fasta | 7027 | 13 | 1.85 | 54.9 | 13.12 |
| 109018_pLenti-X1-Neo-GFP-ATL1.fasta | 11423 | 16 | 1.40 | 54.6 | 14.99 |
| 109432_Non-specific_sgRNA.fasta | 2279 | 3 | 1.32 | 55.0 | 13.55 |
| 10955_pUC8_mouse_IFN-gamma.fasta | 13146 | 8 | 0.61 | 53.6 | 15.81 |
| 109620_actc1b-mKate2-rab30.fasta | 9534 | 7 | 0.73 | 54.1 | 14.53 |
| 109670_p3E-rab15.fasta | 3277 | 2 | 0.61 | 58.7 | 10.29 |
| 109862_pHH0103_EEF1A2_GTP-EFTU.fasta | 7449 | 13 | 1.75 | 59.8 | 15.03 |
| 110866_pLenti-HF1RA-P2A-GFP-PGK-Puro.fasta | 12448 | 12 | 0.96 | 55.9 | 14.03 |
| 110963_MAEBL-COMP-blac-flag-his.fasta | 13559 | 11 | 0.81 | 53.7 | 9.44 |
| 111049_pCHKU28.1-2.2.fasta | 13384 | 27 | 2.02 | 55.4 | 14.33 |
| 111111_pCR1165.fasta | 8937 | 10 | 1.12 | 55.3 | 13.95 |
| 111175_LT3REVIN.fasta | 10319 | 14 | 1.36 | 55.6 | 14.46 |
| 111748_pBacMam2-DiEx-LIC-C-flag_huntingtin_full-length_Q52.fasta | 15867 | 20 | 1.26 | 55.4 | 14.77 |
| 111844_pCEP4-BG505.SOSIP-TM.fasta | 12541 | 17 | 1.36 | 56.2 | 13.54 |
| 112467_pELF4.1.0-gDNA.fasta | 8545 | 8 | 0.94 | 56.3 | 12.80 |
| 112965_pCI-Neo-JAM-A.fasta | 6374 | 7 | 1.10 | 56.9 | 13.35 |
| 113346_pClodAcytCh-GFP-P4M.fasta | 4640 | 4 | 0.86 | 55.4 | 11.12 |
| 113929_pCDNA3.1-H-Luc_cyto.fasta | 6936 | 9 | 1.30 | 58.6 | 14.19 |
| 114001_pOGG083.fasta | 4062 | 7 | 1.72 | 63.2 | 15.45 |
| 114261_pMXBP_pFa.fasta | 8065 | 10 | 1.24 | 58.8 | 15.35 |
| 114708_pTK711.fasta | 15949 | 22 | 1.38 | 55.6 | 12.77 |
| 115245_NODALvar-MYC-DYK.fasta | 5895 | 6 | 1.02 | 56.6 | 15.48 |
| 116291_pHAGE-EGFR-R108G.fasta | 12122 | 16 | 1.32 | 55.7 | 14.30 |
| 116332_pHAGE-ERBB2-G776S.fasta | 12256 | 15 | 1.22 | 54.2 | 14.34 |
| 116635_pHAGE-PTEN-S170N.fasta | 9724 | 13 | 1.34 | 54.5 | 15.08 |
| 116946_pMRXIP-GFP-YKT6.fasta | 7480 | 7 | 0.94 | 57.8 | 12.42 |
| 117157_pALPS_puro_miR30-IRF5.fasta | 6430 | 9 | 1.40 | 54.3 | 11.57 |
| 118085_p5-71.fasta | 3093 | 6 | 1.94 | 58.8 | 15.54 |
| 118273_pAAV-EF1a-FLEX_frt_-GCaMP6f-WPRE.fasta | 6954 | 10 | 1.44 | 55.0 | 14.06 |
| 118401_pBHDuet64.fasta | 4560 | 7 | 1.54 | 59.0 | 12.55 |
| 118577_pMRE-Tn7-165.fasta | 16609 | 20 | 1.20 | 58.4 | 13.92 |
| 119288_pcDNA3.1_CIP2A_1-905_L533E_V5_His.fasta | 8252 | 8 | 0.97 | 55.6 | 14.07 |

*(showing first 50 of 484 files -- see per_file_summary.csv for all)*