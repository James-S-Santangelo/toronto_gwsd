rule depth_hcn_loci:
    input:
        bams = rules.create_bam_list_highQualSamples.output,
        ref = REFERENCE_GENOME
    output:
        '{0}/depth/{{gene}}/{{gene}}_samtools.depth'.format(HCN_LOCI_DIR)
    log: 'logs/depth_hcn_loci/{gene}_samtools_depth.log'
    conda: '../envs/hcn_loci.yaml'
    params:
        region = lambda w: 'CM019103.1:18570521-20572232' if w.gene == 'ac' else 'CM019108.1:29227327-31233709'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '03:00:00'
    shell:
        """
        samtools depth -a \
            -f {input.bams} \
            -d 0 \
            -o {output} \
            -q 20 \
            -Q 30 \
            -r {params.region} 2> {log}
        """

rule hcn_loci_done:
    input:
        expand(rules.depth_hcn_loci.output, gene=['ac','li'])
    output:
        '{0}/hnc_loci.done'.format(HCN_LOCI_DIR)
    shell:
        """
        touch {output}
        """
