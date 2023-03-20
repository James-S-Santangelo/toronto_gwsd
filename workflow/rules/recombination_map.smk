# Rules for fitting and interpolating recombination maps
# Relies on genetic positions from Olsen et al. 2021 and GBS marker sequences sent by Olsen

rule blast_markers:
    input:
        markers = ancient('{0}/{{map_pop}}_tags.fa'.format(GENMAP_RESOURCE_DIR)),
        db_flag = rules.makeblastdb_fromRef.output,
        ref = REFERENCE_GENOME 
    output:
        '{0}/{{map_pop}}_marker_blast.txt'.format(GENMAP_RESULTS_DIR)
    conda: '../envs/recombination_map.yaml'
    log: LOG_DIR + '/blast_markers/{map_pop}_marker_blast.log'
    params: 
        outfmt = "'6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovs qcovhsp'"
    shell:
        """
        blastn -db {input.ref} \
            -query {input.markers} \
            -out {output} \
            -outfmt {params.outfmt} \
            -num_threads {threads} \
            -evalue 1e-10 \
            -max_hsps 5 \
            -max_target_seqs 5 2> {log}  
        """

rule sites_toInterpolate_byChrom:
    input:
        lambda w: expand(rules.remove_duplicate_sites.output, chrom=w.chrom, site_type='snps', miss='0')
    output:
        '{0}/genMap_interpolation/{{chrom}}_forGenMapInterpolation.sites'.format(PROGRAM_RESOURCE_DIR)
    shell:
        """
        grep -v '^#' {input} | cut -f1,2 > {output}
        """

rule interpolate_genetic_map:
    input:
        DG_markers = expand(rules.blast_markers.output, map_pop='DG'),
        SG_markers = expand(rules.blast_markers.output, map_pop='SG'),
        DG_genMap = ancient('{0}/DG_genMap.csv'.format(GENMAP_RESOURCE_DIR)),
        SG_genMap = ancient('{0}/SG_genMap.csv'.format(GENMAP_RESOURCE_DIR)),
        DG_names = ancient('{0}/DG_marker_key.csv'.format(GENMAP_RESOURCE_DIR)),
        SG_names = ancient('{0}/SG_marker_key.csv'.format(GENMAP_RESOURCE_DIR)),
        sites = expand(rules.sites_toInterpolate_byChrom.output, chrom=CHROMOSOMES)
    output:
        markers_byPop_plot = '{0}/genMap/markers_bothPops_withOutliers_allChroms.pdf'.format(FIGURES_DIR),
        scamFits_plot = '{0}/genMap/scamFits_allChroms.pdf'.format(FIGURES_DIR),
        genMap_interp = '{0}/genMap_interpolated_allChroms.txt'.format(GENMAP_RESULTS_DIR)
    conda: '../envs/recombination_map.yaml'
    notebook:
        "../notebooks/genMap_interpolation.r.ipynb"

rule split_genMap:
    input:
        rules.interpolate_genetic_map.output.genMap_interp
    output:
        '{0}/{{chrom}}_genMap.txt'.format(GENMAP_RESULTS_DIR)
    shell:
        """
        head -n 1 {input} > {output} && grep {wildcards.chrom} {input} >> {output}
        """

rule recombination_map_done:
    input:
        expand(rules.blast_markers.output, map_pop = ['SG', 'DG']),
        expand(rules.sites_toInterpolate_byChrom.output, chrom = CHROMOSOMES),
        expand(rules.split_genMap.output, chrom = CHROMOSOMES)
    output:
        '{0}/recombination_map_done'.format(GENMAP_RESULTS_DIR)
    shell:
        """
        touch {output}
		"""
