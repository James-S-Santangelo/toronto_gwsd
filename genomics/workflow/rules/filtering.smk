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
        if [ {wildcards.site_type} = 'invariant' ]; then
            bcftools filter -i 'F_MISSING <= {wildcards.miss}' {input.vcf} -Oz > {output.vcf} 2> {log} &&
            sleep 5
            tabix {output.vcf} 
        elif [ {wildcards.site_type} = 'snps' ]; then
            ( bcftools filter -i 'F_MISSING <= {wildcards.miss}' {input.vcf} |
                bcftools filter -i 'QUAL >= 30' |
                bcftools filter -i 'AB >= 0.25 & AB <= 0.75 | AB <= 0.01' |
                bcftools filter -i 'SAF > 0 & SAR > 0' |
                bcftools filter -i 'MQM >=30 & MQMR >= 30' |
                bcftools filter -i '((PAIRED > 0.05) & (PAIREDR > 0.05) & (PAIREDR / PAIRED < 1.75 ) & (PAIREDR / PAIRED > 0.25)) | ((PAIRED < 0.05) & (PAIREDR < 0.05))' |
                bcftools filter -o {output.vcf} -O z -i '((AF > 0) & (AF < 1))' && sleep 5 && tabix {output.vcf} ) 2> {log}
        fi
        """

rule remove_duplicate_sites:
    input:
        vcf = rules.bcftools_filter_vcfs.output.vcf,
        tbi = rules.bcftools_filter_vcfs.output.tbi
    output:
        vcf = temp('{0}/vcf/{{chrom}}/{{chrom}}_allFinalSamples_{{site_type}}_miss{{miss}}_filtered_noDups.vcf.gz'.format(FREEBAYES_DIR)),
        tbi = temp('{0}/vcf/{{chrom}}/{{chrom}}_allFinalSamples_{{site_type}}_miss{{miss}}_filtered_noDups.vcf.gz.tbi'.format(FREEBAYES_DIR))
    log: LOG_DIR + '/remove_duplicate_sites/{chrom}_{site_type}_miss{miss}_removeDups.log'
    conda: '../envs/filtering.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '01:00:00'
    shell:
        """
        vcfuniq {input.vcf} | bcftools view -Oz > {output.vcf} 2> {log}
        sleep 5
        tabix {output.vcf}
        """

rule concat_variant_invariant_sites:
    input:
        vcfs = lambda w: expand(rules.remove_duplicate_sites.output.vcf, chrom=w.chrom, site_type=['snps', 'invariant'], miss=w.miss),
        tbis= lambda w: expand(rules.remove_duplicate_sites.output.tbi, chrom=w.chrom, site_type=['snps', 'invariant'], miss=w.miss)
    output:
        vcf = temp(f"{FREEBAYES_DIR}/vcf/{{chrom}}/{{chrom}}_miss{{miss}}_allFinalSamples_allSites_filtered.vcf.gz"),
        tbi = temp(f"{FREEBAYES_DIR}/vcf/{{chrom}}/{{chrom}}_miss{{miss}}_allFinalSamples_allSites_filtered.vcf.gz.tbi")
    log: f"{LOG_DIR}/concat_variant_invariant_sites/{{chrom}}_miss{{miss}}_concat.log"
    conda: '../envs/filtering.yaml'
    shell:
        """
        bcftools concat -a {input.vcfs} | bcftools sort -Oz -o {output.vcf} 2> {log}
        sleep 5
        tabix {output.vcf}
        """
