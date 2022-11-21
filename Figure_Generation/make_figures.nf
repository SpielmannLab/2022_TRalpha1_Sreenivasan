// Usage : nextflow run /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/make_figures.nf -with-timeline
// do this either after launching an srun or sbatch, and going into $SCRATCH
// do this in the conda environment HypthmsMittag

// define input parameters
params.file_sc_obj = "/data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/multi_sample_set7_merge/cluster/hypthmsMittag-merge_clustering.rds"
params.file_markers = "/data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/multi_sample_set7_merge/cluster/hypthmsMittag-merge_markers.tsv"
file_functions="/data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/functions.R"
params.res = "0.008" //resolution to be used for making the figures. If changing this, be careful. Lots of dependencies
params.file_linnarson = "/data/humangen_mouse/hypthmsMittag/analysis/RefMapping/linnarson_data/l5_all.rds" // linnarson data
file_hypomap = "/data/humangen_mouse/hypthmsMittag/hypoMap/hypoMap.rds"

//define outpuut folder
params.outfolder = "/data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/figures/"

// make qc on non-neuronal cells based on markers
// make qc and characterisation figure for the paper.
code="/data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_nonneuronal.R"
process figure_nonneuronal {
  publishDir params.outfolder+'figure_nonneuronal', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  path file_markers from params.file_markers
  path "figure_nonneuronal.R" from code
  output:
  file "*.txt" into figure_nonneuronal
  """
  Rscript figure_nonneuronal.R ${file_markers}
  """
}

// make qc and characterisation figure for the paper.
code="/data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_fulldata.R"
process figure_fulldata {
  publishDir params.outfolder+'figure_fulldata', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  path file_sc_obj from params.file_sc_obj
  path file_markers from params.file_markers
  path file_functions
  val res from params.res
  path "figure_fulldata.R" from code
  output:
  file "*.pdf" into figure_qc_pdfs
  file "*.png" into figure_qc_pngs
  file "*x.rds" into file_sc_obj_ch, file_sc_obj_ch_copy1, file_sc_obj_ch_copy2, file_sc_obj_ch_copy3, file_sc_obj_ch_copy4, file_sc_obj_ch_copy5
  file "*.tsv" into cell_annotations_geo
  """
  Rscript figure_fulldata.R --file_sc_obj=${file_sc_obj} --file_markers=${file_markers} --file_functions=${file_functions} --res=${res}
  tree
  """
}

// integrate the entire dataset using harmony
process integrate_fulldata_seurat{
  publishDir params.outfolder+'figure_fulldata_int_seurat', mode: 'copy', pattern: "*_int.rds", overwrite: true
  input:
  file file_sc_obj from file_sc_obj_ch_copy4
  path file_functions
  output:
  file "*_int.rds" into rds_fulldata_int_seurat, rds_fulldata_int_seurat_copy1
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/integrate_clusters.R --file_sc_obj=${file_sc_obj} --file_functions=${file_functions} --task="integrate-seurat"
  """
}

// make figures from the full data(seurat integrated)
process figure_fulldata_int_seurat {
  publishDir params.outfolder+'figure_fulldata_int_seurat', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  path file_sc_obj from rds_fulldata_int_seurat_copy1
  path file_markers from params.file_markers
  path file_functions
  val res from params.res
  output:
  file "*.pdf" into figure_fulldata_int_seurat_pdfs
  file "*.png" into figure_fulldata_int_seurat_pngs
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_fulldata_int_seurat.R --file_sc_obj=${file_sc_obj} --file_markers=${file_markers} --file_functions=${file_functions} --res=${res}
  """
}

// integrate the entire dataset using harmony
process integrate_fulldata_harmony{
  publishDir params.outfolder+'figure_fulldata_int_harmony', mode: 'copy', pattern: "*_int.rds", overwrite: true
  input:
  file file_sc_obj from file_sc_obj_ch_copy3
  path file_functions
  output:
  file "*_int.rds" into rds_fulldata_int_harmony
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/integrate_clusters.R --file_sc_obj=${file_sc_obj} --file_functions=${file_functions} --task="integrate-harmony"
  """
}

// cluster the entire data set (integrated)
process cluster_fulldata_int_harmony{
  publishDir params.outfolder+'figure_fulldata_int_harmony', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  file rds_fulldata_int_harmony
  path file_functions
  output:
  file "*_clus.rds" into rds_fulldata_int_harmony_clusters
  file "*.pdf" into figure_fulldata_int_harmony_clusters
  file "Markers_resolution*.tsv" into markers_fulldata_int_harmony_clusters
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/subcluster_degenes.R --file_sc_obj=${rds_fulldata_int_harmony} --res="0.008" --task="integrate-harmony" --file_functions=${file_functions}
  """
}

// make figures from the full data(harmony integrated) as well as comparison
process figure_fulldata_int_harmony {
  publishDir params.outfolder+'figure_fulldata_int_harmony', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  path file_sc_obj from rds_fulldata_int_harmony_clusters
  path file_markers from params.file_markers
  path file_functions
  val res from params.res
  output:
  file "*.pdf" into figure_fulldata_int_harmony_pdfs
  file "*.png" into figure_fulldata_int_harmony_pngs
  file "*.txt" into code_output_int_harmony
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_fulldata_int_harmony.R --file_sc_obj=${file_sc_obj} --file_markers=${file_markers} --file_functions=${file_functions} --res=${res}
  """
}

// now separate the main rds into each cluster
process separate_clusters{
  publishDir params.outfolder+'rds_of_clusters', mode: 'copy', pattern: "*.rds", overwrite: true
  input:
  file file_sc_obj from file_sc_obj_ch
  val res from params.res
  output:
  file "*_cluster*.rds" into rdss_of_clusters, rdss_of_clusters_copy, rdss_of_clusters_copy2 // all clusters
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/separate_clusters.R --file_sc_obj=${file_sc_obj} --res=${res}
  """
}

rds_of_clusters = rdss_of_clusters.flatten() //split to individual channels
rds_of_clusters_copy = rdss_of_clusters_copy.flatten() //split to individual channels

// make UMAPs stratified by repeats and genotype per cluster before integration
process figure_umap_clusters{
  publishDir params.outfolder+'figure_umap_separate_clusters_merged', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  file file_sc_obj from rds_of_clusters
  output:
  file "*.pdf" into figure_umap_clusters_pdfs
  file "*.png" into figure_umap_clusters_pngs
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_umap_separate_clusters.R --file_sc_obj=${file_sc_obj}
  """
}

// integrate individual clusters using harmony
process integrate_clusters{
  publishDir params.outfolder+'rds_of_clusters', mode: 'copy', pattern: "*_int.rds", overwrite: true
  input:
  file file_sc_obj from rds_of_clusters_copy
  path file_functions
  output:
  file "*_int.rds" into rds_of_clusters_int, rds_of_clusters_int_copy, rds_of_clusters_int_copy2, rds_of_clusters_int_copy3, rds_of_clusters_int_copy4
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/integrate_clusters.R --file_sc_obj=${file_sc_obj} --file_functions=${file_functions} --task="integrate-harmony"
  """
}

// make UMAPs stratified by repeats and genotype per cluster after integration
process figure_umap_clusters_int{
  publishDir params.outfolder+'figure_umap_separate_clusters_integrated', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  file file_sc_obj from rds_of_clusters_int_copy
  output:
  file "*.pdf" into figure_umap_separate_clusters_integrated_pdfs
  file "*.png" into figure_umap_separate_clusters_integrated_pngs
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_umap_separate_clusters.R --file_sc_obj=${file_sc_obj}
  """
}

// only Neurons_1
rds_of_neurons_int = rds_of_clusters_int_copy3
  .filter(~/.*cluster0.*/)

// subcluster the Neurons_1
process subcluster_neurons{
  publishDir params.outfolder+'figure_neurons', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  file rds_of_neurons_int
  path file_functions
  output:
  file "*_clus.rds" into rds_of_neurons_clusters, rds_of_neurons_clusters_copy1
  file "*.pdf" into figure_neurons_clusters
  file "Markers_resolution*.tsv" into markers_neurons_clusters
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/subcluster_degenes.R --file_sc_obj=${rds_of_neurons_int} --res="0.1" --task="integrate-harmony" --file_functions=${file_functions}
  """
}

// make neuron figure. The resolution and cluster annotations are hard-coded in the R-script
process figure_neurons {
  publishDir params.outfolder+'figure_neurons', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  path "sc_obj.rds" from rds_of_neurons_clusters
  path file_functions
  file markers_neurons_clusters
  output:
  file "*.pdf" into figure_neurons_pdfs
  file "*.png" into figure_neurons_pngs
  file "annotations*.tsv" into annotations_neurons, annotations_neurons_copy1
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_neurons.R --file_sc_obj="sc_obj.rds" --file_functions=${file_functions} --file_markers=${markers_neurons_clusters} --res="0.1"
  """
}

// Make Parvalbumin figure
process figure_pvalb_neurons {
  publishDir params.outfolder+'figure_pvalb_neurons', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  path "sc_obj.rds" from rds_of_neurons_clusters_copy1
  output:
  file "*.pdf" into figure_pvalb_neurons_pdfs
  file "*.tsv" into highly_expressed_genes_pvalb
  file "*.txt" into stats_pvalb
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_pvalb_neurons.R "sc_obj.rds"
  """
}

// only Astrocytes and Tanycytes
rds_of_astrotany_int = rds_of_clusters_int_copy2
  .filter(~/.*cluster1.*/)


// subcluster the Astrocytes_Tanycytes
process subcluster_astrotany{
  publishDir params.outfolder+'figure_astrotany', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  file rds_of_astrotany_int
  path file_functions
  output:
  file "*_clus.rds" into rds_of_astrotany_clusters
  file "*.pdf" into figure_astrotany_clusters
  file "Markers_resolution*.tsv" into markers_astrotany_clusters
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/subcluster_degenes.R --file_sc_obj=${rds_of_astrotany_int} --res="0.01" --task="integrate-harmony" --file_functions=${file_functions}
  """
}

// make Astrocytes_Tanycytes figure. The resolution and cluster annotations are hard-coded in the R-script
process figure_astrotany {
  publishDir params.outfolder+'figure_astrotany', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  path "sc_obj.rds" from rds_of_astrotany_clusters
  path file_functions
  file markers_astrotany_clusters
  output:
  file "*.pdf" into figure_astrotany_pdfs
  file "*.png" into figure_astrotany_pngs
  file "annotations*.tsv" into annotations_astrotany, annotations_astrotany_copy1
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_astrotany.R --file_sc_obj="sc_obj.rds" --file_functions=${file_functions} --file_markers=${markers_astrotany_clusters} --res="0.01"
  """
}

// only oligodendrocytes
rds_of_olig_int = rds_of_clusters_int
  .filter(~/.*cluster2.*/)

rds_of_olig_int.into { rds_of_olig_int; rds_of_olig_int_copy1; rds_of_olig_int_copy2 }

// make plots for oligodendrocytes
code="/data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_oligo.R"
process figure_oligo{
  publishDir params.outfolder+'figure_oligo', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  file file_sc_obj_oligo from rds_of_olig_int  // odc data, harmony integrated pca is "harmony"
  path "figure_oligo.R" from code
  output:
  file "*.pdf" into figure_oligo_pdfs
  file "*.png" into figure_oligo_pngs
  file "*.tsv" into oligo_markers_wtvstra1
  """
  Rscript figure_oligo.R --file_sc_obj_oligo=${file_sc_obj_oligo}
  """
}

// only OPCs
rds_of_opc_int = rds_of_clusters_int_copy4
  .filter(~/.*cluster3.*/)
// make plots for opcs
process figure_opc{
  publishDir params.outfolder+'figure_opc', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  file file_sc_obj_opc from rds_of_opc_int  // opc data, harmony integrated pca is "harmony"
  output:
  file "*.pdf" into figure_opc_pdf
  file "*.png" into figure_opc_png
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_opc.R --file_sc_obj_opc=${file_sc_obj_opc}
  """
}

// add raw and seurat-integrated datasets into a channel for main_sub_cluster DE gene calculation vs genotypes and vs repeats
file_sc_objs_multiple = file_sc_obj_ch_copy1
  .mix( rds_fulldata_int_seurat )

file_sc_objs_multiple
  .into { file_sc_objs_multiple; file_sc_objs_multiple_copy1 }

//de genes vs genotypes
code="/data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_main_sub_cluster_genotype_degenes.R"
process figure_genotype_de_genes_main_sub_clusters{
  publishDir params.outfolder+'figure_genotype_de_genes_main_sub_clusters', mode: 'copy', overwrite: true // The output figures will be saved
  input:
  path "figure_main_sub_cluster_genotype_degenes.R" from code
  each file_sc_obj from file_sc_objs_multiple
  file annotations_neurons
  file annotations_astrotany
  output:
  file "*.tsv" into file_all_markers_genotype
  file "*.p??" into figure_nDEgenes_genotype
  file "*.tsv" into markers_all_genotype
  file "*.txt" into code_output_genotype
  """
  Rscript figure_main_sub_cluster_genotype_degenes.R --file_sc_obj=${file_sc_obj} --file_annotations_neurons=${annotations_neurons} --file_annotations_astrotany=${annotations_astrotany}
  """
}

//de genes vs repeats
code="/data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_main_sub_cluster_rep_degenes.R"
process figure_rep_de_genes_main_sub_clusters{
  publishDir params.outfolder+'figure_rep_de_genes_main_sub_clusters', mode: 'copy', overwrite: true // The output figures will be saved
  input:
  path "figure_main_sub_cluster_rep_degenes.R" from code
  each file_sc_obj from file_sc_objs_multiple_copy1
  file annotations_neurons_copy1
  file annotations_astrotany_copy1
  output:
  file "*.tsv" into file_all_markers_rep
  file "*.p??" into figure_nDEgenes_rep
  file "*.tsv" into markers_all_rep
  file "*.txt" into code_output_rep
  """
  Rscript figure_main_sub_cluster_rep_degenes.R --file_sc_obj=${file_sc_obj} --file_annotations_neurons=${annotations_neurons_copy1} --file_annotations_astrotany=${annotations_astrotany_copy1}
  """
}

// Integrate HypoMap ODC and OPCS to show that intermediate NFOL etc are pretty much absent.
process figure_hypoMap_integration{
  publishDir params.outfolder+'figure_hypomap', mode: 'copy', overwrite: true // The output figures will be saved
  input:
  file file_sc_obj from file_sc_obj_ch_copy4
  path file_hypomap
  path file_functions
  output:
  file "*.p??" into figure_hypomap
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_hypoMap_integration.R ${file_sc_obj} ${file_hypomap} ${file_functions}
  """
}
