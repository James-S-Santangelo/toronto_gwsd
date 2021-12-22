rule create_bam_list_highQualSamples:
    input:
        expand(rules.samtools_markdup.output.bam, sample=SAMPLES)
    output:
        '{0}/bam_lists/highQualSamples_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/create_bam_list/highQualSamples_bam_list.log'
    run:
        import os
        with open(output[0], 'w') as f:
            for bam in input:
                sample = os.path.basename(bam).split('_merged')[0]
                if sample not in LOWQUAL_SAMPLES_TO_EXCLUDE:
                    f.write('{0}\n'.format(bam))

rule angsd_gl_forNGSrelate:
    input:
        bams = rules.create_bam_list_highQualSamples.output,
        ref = rules.unzip_reference.output,
        sites = rules.split_angsd_sites_byChrom.output,
        idx = rules.angsd_index_sites.output
    output:
        gls = '{0}/gls/ngsrelate/{{chrom}}_ngsRelateSNPs_{{site}}_maf{{maf}}.glf.gz'.format(ANGSD_DIR),
        mafs = '{0}/gls/ngsrelate/{{chrom}}_ngsRelateSNPs_{{site}}_maf{{maf}}.mafs.gz'.format(ANGSD_DIR),
        pos = '{0}/gls/ngsrelate/{{chrom}}_ngsRelateSNPs_{{site}}_maf{{maf}}.glf.pos.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_gl_forNGSrelate/{chrom}_ngsRelateSNPs_{site}_maf{maf}.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    params:
        out = '{0}/gls/ngsrelate/{{chrom}}_ngsRelateSNPs_{{site}}_maf{{maf}}'.format(ANGSD_DIR),
        max_dp = ANGSD_MAX_DP,
        min_dp_ind = ANGSD_MIN_DP_IND_GL
    threads: 12
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 12000,
        time = '12:00:00'
    wildcard_constraints:
        chrom = 'CM019101.1',
        site='4fold'
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND*80/100 ))
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doGlf 3 \
            -doMajorMinor 1 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -doCounts 1 \
            -setMinDepthInd {params.min_dp_ind} \
            -setMaxDepth {params.max_dp} \
            -baq 2 \
            -ref {input.ref} \
            -minInd $MIN_IND \
            -sites {input.sites} \
            -minQ 20 \
            -minMapQ 30 \
            -minMaf {wildcards.maf} \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """

rule convert_freq_forNGSrelate:
    input:
        rules.angsd_gl_forNGSrelate.output.mafs
    output:
        '{0}/gls/ngsrelate/{{chrom}}_ngsRelate_{{site}}_maf{{maf}}.freqs'.format(ANGSD_DIR)
    log: LOG_DIR + '/convert_freq_forNGSrelate/{chrom}_{site}_maf{maf}_convert_freqs.log'
    wildcard_constraints:
        chrom = 'CM019101.1',
        site='4fold'
    shell:
        """
        zcat {input} | cut -f6 | sed 1d > {output} 2> {log}
        """

rule ngsrelate:
    input:
        bams = rules.create_bam_list_highQualSamples.output,
        gls = rules.angsd_gl_forNGSrelate.output.gls,
        freq = rules.convert_freq_forNGSrelate.output
    output:
        '{0}/{{chrom}}_ngsRelate_{{site}}_maf{{maf}}.out'.format(NGSRELATE_DIR)
    log: LOG_DIR + '/ngsrelate/{chrom}_ngsRelate_{site}_maf{maf}.log'
    container: 'library://james-s-santangelo/ngsrelate/ngsrelate:2.0' 
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '02:00:00'
    wildcard_constraints:
        chrom = 'CM019101.1',
        site='4fold'
    shell:
        """
        N=$( wc -l < {input.bams} );
        ngsRelate -f {input.freq} \
            -O {output} \
            -g {input.gls} \
            -p {threads} \
            -n $N 2> {log}
        """

rule ngsrelate_done:
    input:
        expand(rules.ngsrelate.output, chrom='CM019101.1', maf=['0.05'], site=['4fold'])
    output:
        '{0}/ngsrelate.done'.format(NGSRELATE_DIR)
    shell:
        """
        touch {output}
        """
