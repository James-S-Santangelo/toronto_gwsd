# Rules for fitting and interpolating recombination maps
# Relies on genetic positions from Olsen et al. 2021 and GBS marker sequences sent by Olsen

rule sites_toInterpolate_byChrom:
    """
    Extract sites from VCF to be used for genetic map interpolation
    """
    input:
        lambda w: expand(rules.remove_duplicate_sites.output.vcf, chrom=w.chrom, site_type='snps', miss='0')
    output:
        '{0}/genMap_interpolation/{{chrom}}_forGenMapInterpolation.sites'.format(PROGRAM_RESOURCE_DIR)
    shell:
        """
        zgrep -v '^#' {input} | cut -f1,2 > {output}
        """

rule interpolate_genetic_map:
    """
    Run script to interpolate genetic map
    """
    input:
        linkage_map = config['genmap'],
        sites = expand(rules.sites_toInterpolate_byChrom.output, chrom=CHROMOSOMES)
    output:
        scamFits_plot = '{0}/gen_map/scamFits_allChroms.pdf'.format(FIGURES_DIR),
        genMap_interp = '{0}/genMap_interpolated_allChroms.txt'.format(GENMAP_RESULTS_DIR)
    conda: '../envs/recombination_map.yaml'
    script:
        "../scripts/r/genMap_interpolation.R"

rule split_genMap:
    """
    Split genetic map by chromosome
    """
    input:
        rules.interpolate_genetic_map.output.genMap_interp
    output:
        '{0}/{{chrom}}_genMap.txt'.format(GENMAP_RESULTS_DIR)
    shell:
        """
        head -n 1 {input} > {output} && grep {wildcards.chrom} {input} >> {output}
        """

rule recombination_map_done:
    """
    Create empty flag file signaling completion of recombination map interpolation
    """
    input:
        expand(rules.sites_toInterpolate_byChrom.output, chrom = CHROMOSOMES),
        expand(rules.split_genMap.output, chrom = CHROMOSOMES)
    output:
        '{0}/recombination_map.done'.format(GENMAP_RESULTS_DIR)
    shell:
        """
        touch {output}
		"""
