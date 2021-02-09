rule bams_list_byPop_multiInd:
    input:
        rules.create_bam_list_varCall.output
    output:
        '{0}/population_bam_lists/{{popu}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
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
        gls = '{0}/ngsRelate_gl/pop{{popu}}/CM019101.1_pop{{popu}}_ngsRelateSNPs_maf{{maf}}.glf.gz'.format(ANGSD_DIR),
        mafs = '{0}/ngsRelate_gl/pop{{popu}}/CM019101.1_pop{{popu}}_ngsRelateSNPs_maf{{maf}}.mafs.gz'.format(ANGSD_DIR),
        pos = '{0}/ngsRelate_gl/pop{{popu}}/CM019101.1_pop{{popu}}_ngsRelateSNPs_maf{{maf}}.glf.pos.gz'.format(ANGSD_DIR)
    log: 'logs/angsd_gl_forNGSrelate/CM019101.1_{popu}_ngsRelateSNPs_maf{maf}.log'
    conda: '../envs/angsd.yaml'
    params:
        out = '{0}/ngsRelate_gl/pop{{popu}}/CM019101.1_pop{{popu}}_ngsRelateSNPs_maf{{maf}}'.format(ANGSD_DIR)
    resources:
        nodes = 1,
        ntasks = CORES_PER_NODE,
        time = '12:00:00'
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
            -setMinDepthInd 3 \
            -baq 2 \
            -ref {input.ref} \
            -minInd $MIN_IND \
            -minQ 20 \
            -minMapQ 30 \
            -minMaf {wildcards.maf} \
            -r CM019101.1 \
            -bam {input.bams} 2> {log}
        """

