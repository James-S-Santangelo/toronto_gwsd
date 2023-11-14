# Rules to statistically phase VCFs using WhatsHap and SHAPEIT4

###############
#### SETUP ####
###############

rule bcftools_concat_filtered_vcfs:
    input:
        lambda wildcards: expand(rules.remove_duplicate_sites.output, 
            chrom=CHROMOSOMES, 
            site_type=['snps'], 
            miss=['0'])
    output:
        vcf = temp('{0}/vcf/allChroms_allFinalSamples_filtered_noDups.vcf.gz'.format(FREEBAYES_DIR)),
        idx = temp('{0}/vcf/allChroms_allFinalSamples_filtered_noDups.vcf.gz.tbi'.format(FREEBAYES_DIR))
    log: LOG_DIR + '/bcftools_concat_filtered_vcfs/bcftools_concat_filtered_vcfs.log'
    conda: '../envs/phasing.yaml'
    threads: 4
    shell:
        """
        ( bcftools concat \
            --output-type z \
            --threads {threads} \
            --output {output.vcf} \
            {input} && tabix {output.vcf} ) 2> {log}
        """
    
rule bcftools_split_samples:
    input:
        rules.bcftools_concat_filtered_vcfs.output.vcf
    output:
        temp('{0}/vcf/by_sample/{{sample}}.vcf'.format(FREEBAYES_DIR))
    log: LOG_DIR + '/bcftools_split_samples/{sample}_split.log'
    conda: '../envs/phasing.yaml'
    params:
        out = '{0}/vcf/by_sample/'.format(FREEBAYES_DIR)
    shell:
        """
        bcftools +split {input} \
            --samples-file <(echo {wildcards.sample}) \
            --output-type v \
            --output {params.out} 2> {log}
        """

#############################
#### READ-BACKED PHASING ####
#############################

rule whatshap_phase:
    input:
        unpack(get_whatshap_phase_input)
    output:
        vcf = temp('{0}/vcf/by_sample_phased/{{sample}}_whatshapPhased.vcf.gz'.format(FREEBAYES_DIR)),
        idx = temp('{0}/vcf/by_sample_phased/{{sample}}_whatshapPhased.vcf.gz.tbi'.format(FREEBAYES_DIR))
    log: LOG_DIR + '/whatshap_phase/{sample}_whatshap_phase.log'
    conda: '../envs/phasing.yaml'
    params:
        vcf_out = '{0}/vcf/by_sample_phased/{{sample}}_whatshapPhased.vcf'.format(FREEBAYES_DIR),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '03:00:00'
    shell:
        """
        ( whatshap phase \
            --output {params.vcf_out} \
            --reference {input.ref} \
            {input.vcf} {input.bam} \
                && bgzip {params.vcf_out} \
                && tabix {output.vcf} ) 2> {log}
        """

rule bcftools_remove_format_tags:
    input:
        rules.whatshap_phase.output.vcf
    output:
        vcf = temp('{0}/vcf/by_sample_phased/{{sample}}_whatshapPhased_remTag.vcf.gz'.format(FREEBAYES_DIR)),
        idx = temp('{0}/vcf/by_sample_phased/{{sample}}_whatshapPhased_remTag.vcf.gz.tbi'.format(FREEBAYES_DIR))
    log: LOG_DIR + '/bgzip_index_whatshap_vcf/{sample}_remTags.log'
    conda: '../envs/phasing.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        ( bcftools annotate \
            -x FORMAT/AD,FORMAT/AO,FORMAT/QA,FORMAT/GL \
            -O z {input} > {output.vcf} && tabix {output.vcf} ) 2> {log}
        """

rule bcftools_merge_phased:
    input:
        lambda wildcards: expand(rules.bcftools_remove_format_tags.output.vcf, sample=FINAL_SAMPLES)
    output:
        vcf = temp('{0}/vcf/allChroms_allFinalSamples_whatshapPhased.vcf.gz'.format(FREEBAYES_DIR)),
        idx = temp('{0}/vcf/allChroms_allFinalSamples_whatshapPhased.vcf.gz.tbi'.format(FREEBAYES_DIR))
    log: LOG_DIR + '/bcftools_merge_phased/merge_phased.log'
    conda: '../envs/phasing.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '03:00:00'
    threads: 8
    shell:
        """
        ( bcftools merge \
            --info-rules - \
            -O z \
            --threads {threads} \
            {input} > {output.vcf} && tabix {output.vcf} ) 2> {log}
        """

rule split_phased_vcf_byChrom:
    input:
        rules.bcftools_merge_phased.output.vcf
    output:
        vcf = '{0}/vcf/{{chrom}}/{{chrom}}_allFinalSamples_whatshapPhased.vcf.gz'.format(FREEBAYES_DIR),
        idx = '{0}/vcf/{{chrom}}/{{chrom}}_allFinalSamples_whatshapPhased.vcf.gz.tbi'.format(FREEBAYES_DIR)
    log: LOG_DIR + '/split_phased_vcf_byChrom/{chrom}_phased_vcf_split.log'
    conda: '../envs/phasing.yaml'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
       ( bcftools filter \
            -r {wildcards.chrom} \
            -O z \
            -o {output.vcf} \
            --threads {threads} \
            {input} && tabix {output.vcf} ) 2> {log}
        """

############################
#### POPULATION PHASING ####
############################

rule shapeit_phase:
    input:
        vcf = rules.split_phased_vcf_byChrom.output.vcf,
        genMap = rules.split_genMap.output
    output:
        vcf = '{0}/vcf/{{chrom}}/{{chrom}}_allFinalSamples_whatshapPhased_shapeitPhased.vcf.gz'.format(FREEBAYES_DIR),
        idx = '{0}/vcf/{{chrom}}/{{chrom}}_allFinalSamples_whatshapPhased_shapeitPhased.vcf.gz.tbi'.format(FREEBAYES_DIR)
    log: LOG_DIR + '/shapeit_phase/{chrom}_shapeit.log'
    conda: '../envs/phasing.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    threads: 8
    shell:
        """
        shapeit4 --input {input.vcf} \
            --map {input.genMap} \
            --region {wildcards.chrom} \
            --use-PS 0.0001 \
            --thread {threads} \
            --log {log} \
            --output {output.vcf} && tabix {output.vcf}
        """

rule bcftools_concat_phased_vcfs:
    input:
        expand(rules.shapeit_phase.output.vcf, chrom=CHROMOSOMES)
    output:
        vcf = '{0}/vcf/allChroms_allFinalSamples_whatsHapPhased_shapeitPhased.vcf.gz'.format(FREEBAYES_DIR),
        idx = '{0}/vcf/allChroms_allFinalSamples_whatsHapPhased_shapeitPhased.vcf.gz.tbi'.format(FREEBAYES_DIR)
    log: LOG_DIR + '/bcftools_concat_filtered_vcfs/bcftools_concat_phased_vcfs.log'
    conda: '../envs/phasing.yaml'
    threads: 4 
    shell:
        """
        ( bcftools concat \
            --output-type z \
            --threads {threads} \
            --output {output.vcf} \
            {input} && tabix {output.vcf} ) 2> {log}
        """

rule phasing_done:
    input:
        rules.bcftools_concat_phased_vcfs.output
    output:
        "{0}/phasing.done".format(FREEBAYES_DIR)
    shell:
        """
        touch {output}
        """

