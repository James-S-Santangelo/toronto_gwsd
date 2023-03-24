# Rules for fitting and interpolating recombination maps
# Relies on genetic positions from Olsen et al. 2021 and GBS marker sequences sent by Olsen

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
        linkage_map = config['genmap'],
        sites = expand(rules.sites_toInterpolate_byChrom.output, chrom=CHROMOSOMES)
    output:
        scamFits_plot = '{0}/genMap/scamFits_allChroms.pdf'.format(GENMAP_RESULTS_DIR),
        genMap_interp = '{0}/genMap_interpolated_allChroms.txt'.format(GENMAP_RESULTS_DIR)
    conda: '../envs/recombination_map.yaml'
    notebook:
        "../scripts/r/snakemake/genMap_interpolation.r"

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
        expand(rules.sites_toInterpolate_byChrom.output, chrom = CHROMOSOMES),
        expand(rules.split_genMap.output, chrom = CHROMOSOMES)
    output:
        '{0}/recombination_map_done'.format(GENMAP_RESULTS_DIR)
    shell:
        """
        touch {output}
		"""
