# 2022_TRalpha1_Sreenivasan
Following are the scripts used to analyse data, as presented in the manuscript

1. Run Preprocess_fastq/script_count_wrap.sh to generate gene expression count matrices from sequencing fastq files
2. Run Main_Analysis/vs_sc_p_sample_main.sh to perform a standard analysis of the data. This wrapper calls a Nextflow script (Main_Analysis/src/sc_p_sample.nf), which includes multiple steps including importing the data, QC, filtering, clustering, dimensionality reduction, and DE gene analysis. This is run in the conda environment (Main_Analysis/ScRNA.yml)
3. Run Main_Analysis/sc_multi_sample_main.sh to combine the samples. This wrapper calls a Nextflow script (Main_Analysis/src/sc_multi_sample.nf), which combines the data, performs dimensionality reduction, clustering and DE gene analysis. This is run in the conda environment (Main_Analysis/ScRNA.yml)
4. Run Figure_Generation/make_figures.nf (nextflow run Figure_Generation/make_figures.nf) to perform custom analysis used in this manuscript as well as generate figures presented in the manuscript. This is run in the conda environment (Figure_Generation/HypthmsMittag.yml)


[1] Absolute paths are used at certain instances, which will need to be adapted, as needed.  
[2] Anaconda was used to set up the necessary environments
