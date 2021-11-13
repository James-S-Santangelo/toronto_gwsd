# Rules to statistically phase VCFs using WhatsHap and SHAPEIT4

rule whatshap_phase:
    input:
        unpack(get_whatshap_phase_input)
    output:
        '{0}/vcf/{{chrom}}/{{chrom}}_allFinalSamples_{{site_type}}_whatshapPhased.vcf.gz'.format(FREEBAYES_DIR)
    log: 'logs/whatshap_phase/{chrom}_{site_type}_whatshap_phase.log'
    conda: '../envs/phasing.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 30000,
        time = '06:00:00'
    shell:
        """
        whatshap phase \
            --output {output} \
            --reference {input.ref} \
            --chromosome {wildcards.chrom} \
            {input.vcf} {input.bams} 2> {log}
        """

rule phasing_done:
    input:
        expand(rules.whatshap_phase.output, chrom=CHROMOSOMES, site_type='snps')
    output:
        "{0}/phasing.done".format(FREEBAYES_DIR)
    shell:
        """
        touch {output}
        """

