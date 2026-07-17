# Batch ORF Analysis Summary (no ESM-2)

- **Input directory**: wiki_dna_fasta
- **FASTA files processed**: 500
- **Total bases across all files**: 1,294,958 bp
- **Start codons searched**: ATG, GTG, TTG, CTG (non-ATG codons treated as alternative/non-traditional starts)
- **Minimum ORF length filter**: > 80 amino acids
- **k-mer size for repeat detection**: 6 bp
- **Total ORFs found**: 1221
- **ORFs with a non-traditional (non-ATG) start codon**: 1069 (87.6% of all ORFs)
- **Files with zero qualifying ORFs**: 113

## Descriptive Statistics

- **ORFs per file**: n=500, mean=2.44, median=2.00, stdev=2.45, min=0.00, max=15.00
- **ORF density (per kb)**: n=500, mean=0.878, median=0.791, stdev=0.779, min=0.000, max=4.730
- **CAI codon match (%), across all 1221 ORFs**: n=1221, mean=51.3, median=51.1, stdev=3.9, min=38.3, max=68.3
- **Whole-sequence Shannon entropy (bits, 0-2)**: n=500, mean=1.954, median=1.958, stdev=0.018, min=1.839, max=1.986
- **Whole-sequence repeat fraction (k=6)**: n=500, mean=0.685, median=0.716, stdev=0.117, min=0.258, max=0.849
- **Per-ORF Shannon entropy (bits, 0-2)**: n=1221, mean=1.925, median=1.938, stdev=0.054, min=1.523, max=1.996
- **Per-ORF repeat fraction (k=6)**: n=1221, mean=0.374, median=0.357, stdev=0.118, min=0.107, max=0.833

## Per-File Results (first 50)

(Full detail in `per_file_summary.csv`; individual ORFs in `all_orfs.csv`.)

| File | Length (bp) | # ORFs | ORFs/kb | Mean CAI (%) | Non-ATG starts | Entropy | Repeat frac |
|------|-------------|--------|---------|---------------|-----------------|---------|--------------|
| 0000_Sireethorn_Leearamwat.fasta | 1037 | 0 | 0.00 | N/A | 0 | 1.939 | 0.602 |
| 0001_Pawling_village_New_York.fasta | 1401 | 1 | 0.71 | 50.3 | 1 | 1.945 | 0.643 |
| 0002_Miguel_ngel_Asturias.fasta | 1740 | 0 | 0.00 | N/A | 0 | 1.949 | 0.683 |
| 0003_Singleton_New_South_Wales.fasta | 3640 | 0 | 0.00 | N/A | 0 | 1.963 | 0.788 |
| 0004_Municipalities_of_Chihuahua.fasta | 4693 | 6 | 1.28 | 51.6 | 6 | 1.958 | 0.777 |
| 0005_New_London_County_Connecticut.fasta | 1716 | 4 | 2.33 | 49.2 | 3 | 1.922 | 0.714 |
| 0006_Ruja_Ignatova.fasta | 1493 | 5 | 3.35 | 51.0 | 5 | 1.930 | 0.689 |
| 0007_Helping_hand.fasta | 307 | 0 | 0.00 | N/A | 0 | 1.928 | 0.318 |
| 0008_The_Viking_Queen.fasta | 2289 | 2 | 0.87 | 53.1 | 2 | 1.936 | 0.697 |
| 0009_Child_poverty.fasta | 3254 | 3 | 0.92 | 52.3 | 3 | 1.965 | 0.754 |
| 0010_Marina_Mohnen.fasta | 2976 | 2 | 0.67 | 51.7 | 2 | 1.969 | 0.782 |
| 0011_Baie-Comeau_Drakkar.fasta | 974 | 1 | 1.03 | 52.9 | 1 | 1.951 | 0.558 |
| 0012_University_of_Nebraska_Press.fasta | 654 | 1 | 1.53 | 52.6 | 1 | 1.980 | 0.424 |
| 0013_1926_Women_s_World_Games_standing_long_jump.fasta | 591 | 0 | 0.00 | N/A | 0 | 1.973 | 0.527 |
| 0014_John_McColl_British_Army_officer_.fasta | 1091 | 2 | 1.83 | 54.0 | 2 | 1.975 | 0.508 |
| 0015_First_Ladyship_of_Laura_Bush.fasta | 1640 | 4 | 2.44 | 48.9 | 4 | 1.956 | 0.637 |
| 0016_Meridian_arc.fasta | 4337 | 5 | 1.15 | 49.8 | 4 | 1.927 | 0.847 |
| 0017_James_Prescott_Joule.fasta | 2085 | 1 | 0.48 | 54.8 | 1 | 1.963 | 0.694 |
| 0018_Common_sage.fasta | 1549 | 2 | 1.29 | 52.4 | 2 | 1.963 | 0.685 |
| 0019_List_of_Medal_of_Honor_recipients_for_the_Second_Battle_of_Fort_Fisher.fasta | 1272 | 1 | 0.79 | 50.8 | 1 | 1.935 | 0.648 |
| 0020_Landslide.fasta | 3370 | 0 | 0.00 | N/A | 0 | 1.955 | 0.788 |
| 0021_Acoustic_Kitty.fasta | 4919 | 8 | 1.63 | 51.3 | 8 | 1.965 | 0.789 |
| 0022_Kristi_Yamaguchi.fasta | 4290 | 1 | 0.23 | 52.6 | 1 | 1.962 | 0.777 |
| 0023_National_Living_Wage.fasta | 4748 | 3 | 0.63 | 51.4 | 3 | 1.967 | 0.799 |
| 0024_2008_09_Calgary_Flames_season.fasta | 491 | 0 | 0.00 | N/A | 0 | 1.971 | 0.459 |
| 0025_Gord_Fashoway.fasta | 4002 | 3 | 0.75 | 48.2 | 3 | 1.960 | 0.781 |
| 0026_Kallithea.fasta | 3110 | 2 | 0.64 | 53.9 | 2 | 1.979 | 0.742 |
| 0027_Scarface_rapper_.fasta | 3694 | 4 | 1.08 | 47.6 | 4 | 1.974 | 0.756 |
| 0028_Vitrey-sur-Mance.fasta | 3639 | 1 | 0.27 | 48.4 | 1 | 1.949 | 0.789 |
| 0029_Dongguan.fasta | 468 | 0 | 0.00 | N/A | 0 | 1.954 | 0.495 |
| 0030_Electoral_district_of_Eureka.fasta | 2676 | 3 | 1.12 | 50.8 | 3 | 1.960 | 0.763 |
| 0031_Dan_Feuerriegel.fasta | 1528 | 0 | 0.00 | N/A | 0 | 1.961 | 0.638 |
| 0032_Stupa.fasta | 4731 | 4 | 0.85 | 51.2 | 4 | 1.963 | 0.802 |
| 0033_Car_battery.fasta | 4585 | 6 | 1.31 | 50.1 | 6 | 1.949 | 0.828 |
| 0034_Imperial_Japanese_Armed_Forces.fasta | 3794 | 4 | 1.05 | 51.7 | 4 | 1.947 | 0.775 |
| 0035_Newport_County_A.F.C..fasta | 2027 | 2 | 0.99 | 54.4 | 2 | 1.964 | 0.690 |
| 0036_Peter_Osgood.fasta | 917 | 0 | 0.00 | N/A | 0 | 1.946 | 0.580 |
| 0037_Espeja_de_San_Marcelino.fasta | 4924 | 5 | 1.02 | 51.8 | 5 | 1.963 | 0.833 |
| 0038_Epsom_New_Zealand_electorate_.fasta | 1699 | 2 | 1.18 | 49.6 | 1 | 1.946 | 0.689 |
| 0039_Vr_ly.fasta | 1537 | 4 | 2.60 | 51.6 | 4 | 1.873 | 0.680 |
| 0040_Paulo_Nagamura.fasta | 4939 | 4 | 0.81 | 44.9 | 3 | 1.880 | 0.836 |
| 0041_1862.fasta | 2792 | 1 | 0.36 | 48.8 | 1 | 1.944 | 0.725 |
| 0042_Martin_Schulz.fasta | 4596 | 9 | 1.96 | 50.6 | 9 | 1.970 | 0.793 |
| 0043_Brovello-Carpugnino.fasta | 4596 | 3 | 0.65 | 50.5 | 3 | 1.964 | 0.816 |
| 0044_V_lye.fasta | 850 | 0 | 0.00 | N/A | 0 | 1.936 | 0.609 |
| 0045_Kingman_Reef.fasta | 4754 | 7 | 1.47 | 50.3 | 6 | 1.839 | 0.841 |
| 0046_Bernhard_Lippert.fasta | 3343 | 3 | 0.90 | 51.5 | 3 | 1.950 | 0.751 |
| 0047_Eddie_Filgate.fasta | 2984 | 0 | 0.00 | N/A | 0 | 1.955 | 0.704 |
| 0048_Ruud_Gullit.fasta | 1822 | 2 | 1.10 | 53.8 | 2 | 1.961 | 0.631 |
| 0049_Wright_brothers.fasta | 3785 | 3 | 0.79 | 48.8 | 3 | 1.969 | 0.796 |

*(showing first 50 of 500 files -- see per_file_summary.csv for all)*