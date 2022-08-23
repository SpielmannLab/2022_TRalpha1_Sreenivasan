// Usage : nextflow run /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/make_figures.nf
// do this either after launching an srun or sbatch, and going into $SCRATCH
// do this in the conda environment HypthmsMittag

// define input parameters
params.file_sc_obj = "/data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/multi_sample_set7_merge/cluster/hypthmsMittag-merge_clustering.rds"
params.file_markers = "/data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/multi_sample_set7_merge/cluster/hypthmsMittag-merge_markers.tsv"
params.file_functions="/data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/functions.R"
params.res = "0.008" //resolution to be used for making the figures. If changing this, be careful. Lots of dependencies
params.file_linnarson = "/data/humangen_mouse/hypthmsMittag/analysis/RefMapping/linnarson_data/l5_all.rds" // linnarson data
params.file_markers_neurons="/data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/figures/figure_neurons/Markers_resolution_0.1.tsv"
params.file_markers_astrotany="/data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/figures/figure_astrotany/Markers_resolution_0.01.tsv"

//define outpuut folder
params.outfolder = "/data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/figures/"

// make qc and characterisation figure for the paper.
process figure_fulldata {
  publishDir params.outfolder+'figure_fulldata', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  path file_sc_obj from params.file_sc_obj
  path file_markers from params.file_markers
  path file_functions from params.file_functions
  val res from params.res
  output:
  file "*.pdf" into figure_qc_pdfs
  file "*.png" into figure_qc_pngs
  file "*x.rds" into file_sc_obj_ch, file_sc_obj_ch_copy1, file_sc_obj_ch_copy2
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_fulldata.R --file_sc_obj=${file_sc_obj} --file_markers=${file_markers} --file_functions=${file_functions} --res=${res}
  tree
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
  path file_functions from params.file_functions
  output:
  file "*_int.rds" into rds_of_clusters_int, rds_of_clusters_int_copy, rds_of_clusters_int_copy2, rds_of_clusters_int_copy3, rds_of_clusters_int_copy4
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/integrate_clusters.R --file_sc_obj=${file_sc_obj} --file_functions=${file_functions}
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
  path file_functions from params.file_functions
  output:
  file "*_clus.rds" into rds_of_neurons_clusters
  file "*.pdf" into figure_neurons_clusters
  file "*.tsv" into markers_neurons_clusters
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/subcluster_degenes.R --file_sc_obj=${rds_of_neurons_int} --res="0.1" --task="integrate-harmony" --file_functions=${file_functions}
  """
}

// make neuron figure. The resolution and cluster annotations are hard-coded in the R-script
process figure_neurons {
  publishDir params.outfolder+'figure_neurons', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  file rds_of_neurons_clusters
  path file_functions from params.file_functions
  path file_markers_neurons from params.file_markers_neurons
  output:
  file "*.pdf" into figure_neurons_pdfs
  file "*.png" into figure_neurons_pngs
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_neurons.R --file_sc_obj=${rds_of_neurons_clusters} --file_functions=${file_functions} --file_markers=${file_markers_neurons}
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
  path file_functions from params.file_functions
  output:
  file "*_clus.rds" into rds_of_astrotany_clusters, rds_of_astrotany_clusters_copy1
  file "*.pdf" into figure_astrotany_clusters
  file "*.tsv" into markers_astrotany_clusters
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/subcluster_degenes.R --file_sc_obj=${rds_of_astrotany_int} --res="0.01" --task="integrate-harmony" --file_functions=${file_functions}
  """
}

// make Astrocytes_Tanycytes figure. The resolution and cluster annotations are hard-coded in the R-script
process figure_astrotany {
  publishDir params.outfolder+'figure_astrotany', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  file rds_of_astrotany_clusters
  path file_functions from params.file_functions
  path file_markers_astrotany from params.file_markers_astrotany
  output:
  file "*.pdf" into figure_astrotany_pdfs
  file "*.png" into figure_astrotany_pngs
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_astrotany.R --file_sc_obj=${rds_of_astrotany_clusters} --file_functions=${file_functions} --file_markers=${file_markers_astrotany} --res="0.01"
  """
}

// now separate astrocytes and tancytes into individual clusters
process separate_astrotany_subclusters{
  publishDir params.outfolder+'rds_of_astrotany_subclusters', mode: 'copy', pattern: "*.rds", overwrite: true
  input:
  file file_sc_obj from rds_of_astrotany_clusters_copy1
  output:
  file "*_clus_cluster*.rds" into rdss_of_astrotany_subclusters, rdss_of_astrotany_subclusters_copy1 // all clusters from astrotany cluster
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/separate_clusters.R --file_sc_obj=${file_sc_obj} --res="0.01"
  """
}

rds_of_astrotany_subclusters = rdss_of_astrotany_subclusters.flatten() //split to individual channels

// now integrate the astrocyte and tanycyte subclusters using harmony
process integrate_astrotany_subclusters{
  publishDir params.outfolder+'rds_of_astrotany_subclusters', mode: 'copy', pattern: "*_int.rds", overwrite: true
  input:
  file file_sc_obj from rds_of_astrotany_subclusters
  path file_functions from params.file_functions
  output:
  file "*_int.rds" into rds_of_astrotany_clusters_int, rds_of_astrotany_clusters_int_copy
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/integrate_clusters.R --file_sc_obj=${file_sc_obj} --file_functions=${file_functions}
  """
}

//only Astrocytes
rds_of_astrocyte_int = rds_of_astrotany_clusters_int
  .filter(~/.*clus_cluster0.*/)
rds_of_astrocyte_int.into { rds_of_astrocyte_int; rds_of_astrocyte_int_copy1 }
// make astrocyte figure. No subclustering. Just DE genes and stuff
process figure_astrocyte {
  publishDir params.outfolder+'figure_astrocyte', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  file rds_of_astrocyte_int
  path file_functions from params.file_functions
  output:
  file "*.pdf" into figure_astrocyte_pdfs
  file "*.png" into figure_astrocyte_pngs
  file "*.tsv" into markers_astrocyte
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_astrocyte.R --file_sc_obj=${rds_of_astrocyte_int} --file_functions=${file_functions}
  """
}

//only tanycytes
rds_of_tanycyte_int = rds_of_astrotany_clusters_int_copy
  .filter(~/.*clus_cluster1.*/)
// make tanycyte figure. No subclustering. Just DE genes and stuff
process figure_tanycyte {
  publishDir params.outfolder+'figure_tanycyte', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  file rds_of_tanycyte_int
  path file_functions from params.file_functions
  output:
  file "*.pdf" into figure_tanycyte_pdfs
  file "*.png" into figure_tanycyte_pngs
  file "*.tsv" into markers_tanycyte
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_tanycyte.R --file_sc_obj=${rds_of_tanycyte_int} --file_functions=${file_functions}
  """
}

// only oligodendrocytes
rds_of_olig_int = rds_of_clusters_int
  .filter(~/.*cluster2.*/)

rds_of_olig_int.into { rds_of_olig_int; rds_of_olig_int_copy1; rds_of_olig_int_copy2 }

// make plots for oligodendrocytes
process figure_oligo{
  cache true
  publishDir params.outfolder+'figure_oligo', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  file file_sc_obj_oligo from rds_of_olig_int  // odc data, harmony integrated pca is "harmony"
  file file_sc_obj from file_sc_obj_ch_copy2
  output:
  file "*.pdf" into figure_oligo_pdfs
  file "*.png" into figure_oligo_pngs
  file "*.tsv" into oligo_markers_wtvstra1
  """
  Rscript /data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/code/figure_oligo.R --file_sc_obj=${file_sc_obj} --file_sc_obj_oligo=${file_sc_obj_oligo}
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
