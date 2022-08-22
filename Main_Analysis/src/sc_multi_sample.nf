// Usage : nextflow run /data/humangen_mouse/scpipeline/src/sc_multi_sample.nf -params-file /data/humangen_mouse/scpipeline/sc_multi_sample.yaml --id ${SCRATCH/"/scratch/"/}

inputfiles = Channel
    .fromPath(params.filenames)
    .collect()
	.view()


process params{
publishDir params.outfolder, mode: 'copy', overwrite: true
output:
file "*.txt" into parameter
"""
echo ${params} | tr , '\n' > parameters${params.id}.txt
"""
}

process combine{
publishDir params.outfolder+'combine', mode: 'copy', overwrite: true
input:
file sc_file from inputfiles
output:
file "*.rds" into scobj_combined
"""
Rscript /data/humangen_mouse/scpipeline/src/merge_or_integrate.R --task=${params.task} --method=${params.method} --project=${params.project} --nhvg=${params.nhvg} --ncores=${params.ncores} --npcs=${params.npcs} --reference=${params.reference} --integrateBy=${params.integrateBy} --sample.tree=${params.sampleTree} --nclust=${params.nclust} --covars=${params.covars} ${sc_file}
"""
}

process dim_reduction{
publishDir params.outfolder+'dim_reduc', mode: 'copy', overwrite: true
input:
file scobj from scobj_combined
output:
file "*.rds" into scobj_dim_reduction
file "*.pdf" into dim_reduc_report
"""
Rscript /data/humangen_mouse/scpipeline/src/dimensionality_reduction.R --scobject=${scobj} --npcs='${params.npcs}' --task=${params.task} --ncores='${params.ncores}'
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
Rscript /data/humangen_mouse/scpipeline/src/cluster.R --scobject=${scobj} --res='${params.res}' --task=${params.task} --ncores='${params.ncores}'
"""
}

process de_genes{

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
