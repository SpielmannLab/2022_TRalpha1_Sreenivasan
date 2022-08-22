# 2022_TRalpha1_Sreenivasan
Code used to analyse data, as presented in the manuscript

1. Use CellRanger count pipeline (with include-introns options), to generate feature barcode matrices.
2. Run Main_Analysis/vs_sc_p_sample_main.sh to perform a standard analysis of the data. This wrapper calls a Nextflow script (src/sc_p_sample.nf), which includes multiple steps including importing the data, QC, filtering, clustering, dimensionality reduction, and DE gene analysis.
3. Run Main_Analysis/sc_multi_sample_main.sh to combine the samples. This wrapper calls a Nextflow script (src/sc_multi_sample.nf), which combines the data, performs dimensionality reduction, clustering and DE gene analysis.
4. Run Figure_Generation/make_figures.nf (nextflow run make_figures.nf) to perform custom analysis used in this manuscript as well as generate figures presented in the manuscript.


Note. Absolute paths are used at certain instances, which will need to be adapted, as needed
