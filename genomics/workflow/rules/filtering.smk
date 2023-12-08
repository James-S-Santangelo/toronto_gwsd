# Snakemake rules to filter freebayes VCFs

rule bcftools_filter_vcfs:
    input:
        vcf = rules.bcftools_split_variants.output
    output:
        vcf = temp('{0}/vcf/{{chrom}}/{{chrom}}_allFinalSamples_{{site_type}}_miss{{miss}}_filtered.vcf.gz'.format(FREEBAYES_DIR)),
        tbi = temp('{0}/vcf/{{chrom}}/{{chrom}}_allFinalSamples_{{site_type}}_miss{{miss}}_filtered.vcf.gz.tbi'.format(FREEBAYES_DIR))
    log: LOG_DIR + '/bcftools_filter_vcfs/{chrom}_{site_type}_miss{miss}_filter.log'
    conda: '../envs/filtering.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '03:00:00'
    shell:
        """
        ( bcftools filter -i 'F_MISSING <= {wildcards.miss}' {input.vcf} |
            bcftools filter -i 'QUAL >= 30' |
            bcftools filter -i 'AB >= 0.25 & AB <= 0.75 | AB <= 0.01' |
            bcftools filter -i 'SAF > 0 & SAR > 0' |
            bcftools filter -i 'MQM >=30 & MQMR >= 30' |
            bcftools filter -i '((PAIRED > 0.05) & (PAIREDR > 0.05) & (PAIREDR / PAIRED < 1.75 ) & (PAIREDR / PAIRED > 0.25)) | ((PAIRED < 0.05) & (PAIREDR < 0.05))' |
            bcftools filter -O z -i '((AF > 0) & (AF < 1))' \
            > {output.vcf} && tabix {output.vcf} ) 2> {log}
        """

rule remove_duplicate_sites:
    input:
        rules.bcftools_filter_vcfs.output.vcf
    output:
        temp('{0}/vcf/{{chrom}}/{{chrom}}_allFinalSamples_{{site_type}}_miss{{miss}}_filtered_noDups.vcf'.format(FREEBAYES_DIR))
    log: LOG_DIR + '/remove_duplicate_sites/{chrom}_{site_type}_miss{miss}_removeDups.log'
    conda: '../envs/filtering.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '01:00:00'
    shell:
        """
        vcfuniq {input} > {output} 2> {log}
        """

rule vcf_filtering_done:
    input:
        expand(rules.remove_duplicate_sites.output, chrom=CHROMOSOMES, site_type='snps', miss=['0'])
    output:
        '{0}/filtering.done'.format(FREEBAYES_DIR)
    shell:
        """
        touch {output}
        """
