# Rules to statistically phase VCFs using WhatsHap and SHAPEIT4

###############
#### SETUP ####
###############

rule bcftools_split_samples:
    input:
        lambda wildcards: expand(rules.remove_duplicate_sites.output, chrom=wildcards.chrom, sample=wildcards.sample, site_type=['snps'], miss=['0'])
    output:
        temp('{0}/vcf/{{chrom}}/by_sample/{{sample}}.vcf'.format(FREEBAYES_DIR))
    log: LOG_DIR + '/bcftools_split_samples/{chrom}/{chrom}_{sample}_split.log'
    conda: '../envs/phasing.yaml'
    params:
        out = '{0}/vcf/{{chrom}}/by_sample/'.format(FREEBAYES_DIR)
    shell:
        """
        bcftools +split {input} \
            --samples-file <(echo {wildcards.sample}) \
            --output-type v \
            --output {params.out} 2> {log}
        """

rule whatshap_phase:
    input:
        unpack(get_whatshap_phase_input)
    output:
        '{0}/vcf/{{chrom}}/by_sample/{{chrom}}_{{sample}}_whatshapPhased.vcf'.format(FREEBAYES_DIR)
    log: LOG_DIR + '/whatshap_phase/{chrom}/{chrom}_{sample}_whatshap_phase.log'
    conda: '../envs/phasing.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '03:00:00'
    shell:
        """
        whatshap phase \
            --output {output} \
            --reference {input.ref} \
            --chromosome {wildcards.chrom} \
            {input.vcf} {input.bam} 2> {log}
        """

rule bgzip_index_whatshap_vcf:
    input:
        rules.whatshap_phase.output
    output:
        vcf = temp('{0}/vcf/{{chrom}}/by_sample/{{chrom}}_{{sample}}_whatshapPhased.vcf.gz'.format(FREEBAYES_DIR)),
        idx = temp('{0}/vcf/{{chrom}}/by_sample/{{chrom}}_{{sample}}_whatshapPhased.vcf.gz.tbi'.format(FREEBAYES_DIR))
    log: LOG_DIR + '/bgzip_index_whatshap_vcf/{chrom}_{sample}_bgzip_index.log'
    conda: '../envs/phasing.yaml'
    shell:
        """
        ( bgzip {input} && tabix {output.vcf} 2> {log} ) 
        """

rule bcftools_remove_format_tags:
    input:
        rules.bgzip_index_whatshap_vcf.output.vcf
    output:
        vcf = temp('{0}/vcf/{{chrom}}/by_sample/{{chrom}}_{{sample}}_whatshapPhased_remTag.vcf.gz'.format(FREEBAYES_DIR)),
        idx = temp('{0}/vcf/{{chrom}}/by_sample/{{chrom}}_{{sample}}_whatshapPhased_remTag.vcf.gz.tbi'.format(FREEBAYES_DIR))
    log: LOG_DIR + '/bgzip_index_whatshap_vcf/{chrom}_{sample}_remAD.log'
    conda: '../envs/phasing.yaml'
    shell:
        """
        ( bcftools annotate \
            -x FORMAT/AD,FORMAT/AO,FORMAT/QA,FORMAT/GL \
            -O z {input} > {output.vcf} && tabix {output.vcf}) 2> {log}
        """

rule bcftools_merge_phased:
    input:
        lambda wildcards: expand(rules.bcftools_remove_format_tags.output.vcf, chrom=wildcards.chrom, sample=FINAL_SAMPLES)
    output:
        vcf = '{0}/vcf/{{chrom}}/{{chrom}}_allFinalSamples_whatshapPhased.vcf.gz'.format(FREEBAYES_DIR),
        idx = '{0}/vcf/{{chrom}}/{{chrom}}_allFinalSamples_whatshapPhased.vcf.gz.tbi'.format(FREEBAYES_DIR)
    log: LOG_DIR + '/bcftools_merge_phased/{chrom}_merge.log'
    conda: '../envs/phasing.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        ( bcftools merge \
            --info-rules - \
            -O z \
            {input} > {output.vcf} && tabix {output.vcf} )2> {log}
        """

rule phasing_done:
    input:
        expand(rules.bcftools_merge_phased.output, chrom=CHROMOSOMES)
    output:
        "{0}/phasing.done".format(FREEBAYES_DIR)
    shell:
        """
        touch {output}
        """

