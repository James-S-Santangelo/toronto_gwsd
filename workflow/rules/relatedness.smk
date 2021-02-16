rule bams_list_byPop_multiInd:
    input:
        rules.create_bam_list_highQualSamples.output
    output:
        '{0}/bam_lists/population_bam_lists/{{popu}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/bams_list_byPop_multiInd/{popu}_bam_list.log'
    run:
        import os
        with open(output[0], 'w') as fout:
            with open(input[0], 'r') as fin:
                lines = fin.readlines()
                for line in lines:
                    sline = line.strip()
                    pop = os.path.basename(sline).split('_')[1]
                    if pop == wildcards.popu:
                        fout.write(line)

rule angsd_gl_forNGSrelate:
    input:
        bams = rules.bams_list_byPop_multiInd.output,
        ref = REFERENCE_GENOME
    output:
        gls = '{0}/gls/ngsrelate/pop{{popu}}/{{chrom}}_pop{{popu}}_ngsRelateSNPs_maf{{maf}}.glf.gz'.format(ANGSD_DIR),
        mafs = '{0}/gls/ngsrelate/pop{{popu}}/{{chrom}}_pop{{popu}}_ngsRelateSNPs_maf{{maf}}.mafs.gz'.format(ANGSD_DIR),
        pos = '{0}/gls/ngsrelate/pop{{popu}}/{{chrom}}_pop{{popu}}_ngsRelateSNPs_maf{{maf}}.glf.pos.gz'.format(ANGSD_DIR)
    log: 'logs/angsd_gl_forNGSrelate/{chrom}_{popu}_ngsRelateSNPs_maf{maf}.log'
    conda: '../envs/angsd.yaml'
    params:
        out = '{0}/gls/ngsrelate/pop{{popu}}/{{chrom}}_pop{{popu}}_ngsRelateSNPs_maf{{maf}}'.format(ANGSD_DIR),
        min_dp_ind = ANGSD_MIN_DP_IND_GL
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE,
        time = '12:00:00'
    wildcard_constraints:
        chrom = 'CM019101.1'
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND / 2 ));
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {resources.ntasks} \
            -doGlf 3 \
            -doMajorMinor 1 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -doCounts 1 \
            -setMinDepthInd {params.min_dp_ind} \
            -baq 2 \
            -ref {input.ref} \
            -minInd $MIN_IND \
            -minQ 20 \
            -minMapQ 30 \
            -minMaf {wildcards.maf} \
            -r CM019101.1 \
            -bam {input.bams} 2> {log}
        """

rule convert_freq_forNGSrelate:
    input:
        rules.angsd_gl_forNGSrelate.output.mafs
    output:
        '{0}/gls/ngsrelate/pop{{popu}}/{{chrom}}_pop{{popu}}_ngsRelate_maf{{maf}}.freqs'.format(ANGSD_DIR)
    log: 'logs/convert_freq_forNGSrelate/{chrom}_pop{popu}_maf{maf}_convert_freqs.log'
    wildcard_constraints:
        chrom = 'CM019101.1'
    shell:
        """
        zcat {input} | cut -f6 | sed 1d > {output} 2> {log}
        """

rule ngsrelate:
    input:
        bam = rules.bams_list_byPop_multiInd.output,
        gls = rules.angsd_gl_forNGSrelate.output.gls,
        freq = rules.convert_freq_forNGSrelate.output
    output:
        '{0}/pop{{popu}}/{{chrom}}_pop{{popu}}_ngsRelate_maf{{maf}}.out'.format(NGSRELATE_DIR)
    log: 'logs/ngsrelate/{chrom}_pop{popu}_ngsRelate_maf{maf}.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:ngsrelate_vlatest'
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    wildcard_constraints:
        chrom = 'CM019101.1'
    shell:
        """
        N=$( wc -l < {input.bam} );
        ngsRelate -f {input.freq} \
            -O {output} \
            -g {input.gls} \
            -p {threads} \
            -n $N 2> {log}
        """

rule ngsrelate_done:
    input:
        expand(rules.ngsrelate.output, chrom='CM019101.1', popu=POPS_MULTIPLE_INDS, maf=['0.05'])
    output:
        '{0}/ngsrelate.done'.format(NGSRELATE_DIR)
    shell:
        """
        touch {output}
        """
