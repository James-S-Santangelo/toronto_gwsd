rule test_angsd_baq_GL:
    input:
        bams = rules.create_bam_list_allFinalSamples.output,
        ref = REFERENCE_GENOME
    output:
        saf = '{0}/test_params/GL{{GL}}_baq{{baq}}/{{chrom}}_allFinalSamples_{{site}}_GL{{GL}}_baq{{baq}}.saf.gz'.format(ANGSD_DIR),
        saf_idx = '{0}/test_params/GL{{GL}}_baq{{baq}}/{{chrom}}_allFinalSamples_{{site}}_GL{{GL}}_baq{{baq}}.saf.idx'.format(ANGSD_DIR),
        saf_pos = '{0}/test_params/GL{{GL}}_baq{{baq}}/{{chrom}}_allFinalSamples_{{site}}_GL{{GL}}_baq{{baq}}.saf.pos.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/test_angsd_baq_GL/{chrom}_allFinalSamples_{site}_GL{GL}_baq{baq}.log'
    conda: '../envs/angsd.yaml'
    params:
        out = '{0}/test_params/GL{{GL}}_baq{{baq}}/{{chrom}}_allFinalSamples_{{site}}_GL{{GL}}_baq{{baq}}'.format(ANGSD_DIR),
        max_dp = ANGSD_MAX_DP,
        min_dp_ind = ANGSD_MIN_DP_IND_SFS
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE,
        time = '12:00:00'
    wildcard_constraints:
        chrom='CM019101.1',
        site='allSites'
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND*80/100 ));
        angsd -GL {wildcards.GL} \
            -out {params.out} \
            -nThreads {resources.ntasks} \
            -doCounts 1 \
            -dumpCounts 2 \
            -setMinDepthInd {params.min_dp_ind} \
            -setMaxDepth {params.max_dp} \
            -baq {wildcards.baq} \
            -ref {input.ref} \
            -minInd $MIN_IND \
            -minQ 20 \
            -minMapQ 30 \
            -doSaf 1 \
            -anc {input.ref} \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """

rule test_angsd_sfs_baq_GL:
    input:
        rules.test_angsd_baq_GL.output.saf_idx 
    output:
        '{0}/test_params/GL{{GL}}_baq{{baq}}/{{chrom}}_allFinalSamples_{{site}}_GL{{GL}}_baq{{baq}}.sfs'.format(ANGSD_DIR)
    log: LOG_DIR + '/test_angsd_sfs_baq_GL/{chrom}_allFinalSamples_{site}_GL{GL}_baq{baq}.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 50000,
        time = '03:00:00'
    wildcard_constraints:
        chrom='CM019101.1',
        sample_set='finalSamples_relatedRemoved',
        site='allSites'
    shell:
        """
        realSFS {input} -P {threads} -fold 1 > {output} 2> {log}
        """

rule param_tests_done:
    input:
        expand(rules.test_angsd_sfs_baq_GL.output, chrom='CM019101.1', sample_set='finalSamples_relatedRemoved', site='allSites', GL=['1','2'], baq=['0','1','2'])
    output:
        '{0}/test_params/param_tests.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """
