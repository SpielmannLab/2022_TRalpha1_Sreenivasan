// Usage : nextflow run /data/humangen_mouse/scpipeline/src/sc_p_sample.nf -params-file /data/humangen_mouse/scpipeline/sc_p_sample.yaml --id ${SCRATCH/"/scratch/"/}

in = Channel
	.fromPath(params.inputfolder)
	.view()
sample_names = Channel
	.from(params.samplename)
	.view()

process params{
publishDir params.outfolder, mode: 'copy', overwrite: true
output:
file "*.txt" into parameter
"""
echo ${params} | tr , '\n' > parameters_${params.id}.txt
"""
}

process read_10x{

input:

file input from in
val sample from sample_names

output:
file "*.rds" into scobj_read10x

"""
printf '${input}'
Rscript /data/humangen_mouse/scpipeline/src/read_10x.R --infolder='${input}' --samplename='${sample}'
"""

}
process QC_RAW {
publishDir params.outfolder+'read_10x', mode: 'copy', overwrite: true
errorStrategy 'ignore'
input:
file scobj from scobj_read10x

output:
file "*.rds" into scobj_qc_raw
file "*.pdf" into qc_raw_report
"""
Rscript /data/humangen_mouse/scpipeline/src/qc_sc.R --scobject=${scobj} --dataset='RAW'
"""

}

process filter_qc{
publishDir params.outfolder+'filtered',pattern:'*.pdf', mode: 'copy', overwrite: true
input:
file scobj from scobj_qc_raw

output:
file "*.pdf" into filter_report
file "*.mtx" into sc_matrix_file
file "*.rds" into scobj_filter_qc

"""
Rscript /data/humangen_mouse/scpipeline/src/filter_qc.R --sc_file=${scobj} --mincount_p_gene='${params.mincount_p_gene}' --maxcount_p_gene='${params.maxcount_p_gene}' --mincell_p_gene='${params.mincell_p_gene}' --maxcell_p_gene='${params.maxcell_p_gene}' --mincount_p_cell='${params.mincount_p_cell}' --maxcount_p_cell='${params.maxcount_p_cell}' --mingene_p_cell='${params.mingene_p_cell}' --maxgene_p_cell='${params.maxgene_p_cell}' --maxpct_mt='${params.maxpct_mt}' --maxpct_rb='${params.maxpct_rb}' --rm_mt='${params.rm_mt}' --rm_rb='${params.rm_rb}'
"""
}
process scrublet_filter_doublet{
publishDir params.outfolder+'filtered', pattern:"*.pdf", mode: 'copy', overwrite: true
input:
file sc_mtx from sc_matrix_file
path scobj from scobj_filter_qc
output:
file "*.csv" into scrublet_report
file "*.rds" into scobj_filter_doublet
file "*.pdf" into doublet_report
"""
python /data/humangen_mouse/scpipeline/src/run_scrublet.py ${sc_mtx} '${params.npcs}' '${params.exp_db_rate}'
Rscript /data/humangen_mouse/scpipeline/src/filter_doublets.R --sc_file=${scobj} --threshold='${params.threshold}'
"""
}

process QC_filtered{
publishDir params.outfolder+'filtered', mode: 'copy', overwrite: true
errorStrategy 'ignore'
input:
file scobj from scobj_filter_doublet

output:
file "*.pdf" into qc_filtered_report
file "*.rds" into scobj_qc_filtered
"""
Rscript /data/humangen_mouse/scpipeline/src/qc_sc.R --scobject=${scobj} --dataset='Filtered'
"""
}

process normalize{
errorStrategy 'ignore'
publishDir params.outfolder+'normalized', mode: 'copy', overwrite: true
input:
file scobj from scobj_qc_filtered
output:
file "*.rds" into scobj_normalize
"""
Rscript /data/humangen_mouse/scpipeline/src/normalize.R --scobject=${scobj} --method='${params.method}' --nhvg='${params.nhvg}' --covars='${params.covars}' --genes_regress='${params.genes_regress}' --ncores='${params.ncores}'
"""

}


process dim_reduction{
publishDir params.outfolder+'dim_reduc', mode: 'copy', overwrite: true
input:
file scobj from scobj_normalize
output:
file "*.rds" into scobj_dim_reduction
file "*.pdf" into dim_reduc_report
"""
Rscript /data/humangen_mouse/scpipeline/src/dimensionality_reduction.R --scobject=${scobj} --npcs='${params.npcs}' --ncores='${params.ncores}'
"""

}


process cluster{
publishDir params.outfolder+'cluster', mode: 'copy', overwrite: true
input:
file scobj from scobj_dim_reduction
output:
file "*.rds" into scobj_cluster
file "*.pdf" into cluster_report

"""
Rscript /data/humangen_mouse/scpipeline/src/cluster.R --scobject=${scobj} --res='${params.res}' --ncores='${params.ncores}'
"""
}


process de_genes{
errorStrategy 'ignore'
publishDir params.outfolder+'cluster', mode: 'copy', overwrite: true

input:
file scobj from scobj_cluster
output:
file '*.tsv' into de_genes_tsv
file '*.pdf' into de_gene_report
"""
Rscript /data/humangen_mouse/scpipeline/src/de_genes.R --sc_file=${scobj} --test.use='${params.test_use}' --min.cells.group='${params.min_cell_group}' --min.pct='${params.min_pct}' --logfc.threshold='${params.logfc_threshold}' --resolution='${params.res_de_genes}' --features='${params.features}' --no.of.pages='${params.no_of_pages}' --ncores='${params.ncores}'
"""
}
