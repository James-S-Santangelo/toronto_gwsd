# Rules to estimate relatedness across all 120 samples

###############
#### SETUP ####
###############

rule create_bam_list_highQualSamples:
    input:
        bams = lambda wildcards: expand(rules.subset_bams_degeneracy.output.bam, sample=SAMPLES, site=wildcards.site),
        ref_flag = rules.ref_done.output 
    output:
        '{0}/bam_lists/highQualSamples_{{site}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/create_bam_list/highQualSamples_{site}_bam_list.log'
    run:
        import os
        import logging
        logging.basicConfig(filename=log[0], level=logging.DEBUG)
        try:
            with open(output[0], 'w') as f:
                for bam in input.bams:
                    search = re.search('^(.+)(?=_\w)', os.path.basename(bam))
                    sample = search.group(1) 
                    if sample not in LOWQUAL_SAMPLES_TO_EXCLUDE:
                        f.write('{0}\n'.format(bam))
        except:
            logging.exception("An error occurred!")
            raise

#####################################
#### BINARY GENOTYPE LIKELIHOODS ####
#####################################

rule angsd_gl_forNGSrelate:
    input:
        bams = rules.create_bam_list_highQualSamples.output,
        ref = rules.unzip_reference.output,
        sites = rules.split_random_angsd_sites_byChrom.output,
        idx = rules.index_random_chromosomal_angsd_sites.output
    output:
        gls = temp('{0}/gls/ngsrelate/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_forNGSrelate.glf.gz'.format(ANGSD_DIR)),
        mafs = temp('{0}/gls/ngsrelate/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_forNGSrelate.mafs.gz'.format(ANGSD_DIR)),
        pos = temp('{0}/gls/ngsrelate/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_forNGSrelate.glf.pos.gz'.format(ANGSD_DIR))
    log: LOG_DIR + '/angsd_gl_forNGSrelate/{chrom}_{site}_maf{maf}.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    params:
        out = '{0}/gls/ngsrelate/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_forNGSrelate'.format(ANGSD_DIR),
        max_dp = ANGSD_MAX_DP,
        min_dp_ind = ANGSD_MIN_DP_IND_GL
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '6:00:00'
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

###################
#### NGSRELATE ####
###################

rule convert_freq_forNGSrelate:
    input:
        rules.angsd_gl_forNGSrelate.output.mafs
    output:
        '{0}/gls/ngsrelate/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_forNGSrelate.freqs'.format(ANGSD_DIR)
    log: LOG_DIR + '/convert_freq_forNGSrelate/{chrom}_{site}_maf{maf}_convert_freqs.log'
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
        '{0}/{{chrom}}_{{site}}_maf{{maf}}_NGSrelate.out'.format(NGSRELATE_DIR)
    log: LOG_DIR + '/ngsrelate/{chrom}_{site}_maf{maf}.log'
    container: 'library://james-s-santangelo/ngsrelate/ngsrelate:2.0' 
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '02:00:00'
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
        expand(rules.ngsrelate.output, maf=['0.05'], site=['4fold'], chrom=CHROMOSOMES)
    output:
        '{0}/ngsrelate.done'.format(NGSRELATE_DIR)
    shell:
        """
        touch {output}
        """
