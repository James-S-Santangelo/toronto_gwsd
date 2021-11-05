# Snakemake rules to filter freebayes VCFs

rule prop_sites_missing_bySample:
    """
    Estimates the proportion of sites with missing GT calls for each sample
    """
    input:
        vcf = rules.bcftools_split_variants.output,
        idx = rules.tabix_vcf.output
    output:
        '{0}/filtering/{{chrom}}_{{site_type}}_prop_sites_missing_bySamples.imiss'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/prop_sites_missing_bySample/{chrom}_{site_type}_prop_sites_missing_bySample.log'
    conda: '../envs/filtering.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '03:00:00'
    threads: 4
    params:
        out = '{0}/filtering/{{chrom}}_{{site_type}}_prop_sites_missing_bySamples'.format(PROGRAM_RESOURCE_DIR)
    shell:
        """
        vcftools --gzvcf {input.vcf} --missing-indv --out {params.out}
        """
        
