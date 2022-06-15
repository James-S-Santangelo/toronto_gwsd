# Rules to estimate within-population diversity and pairwise Fst/PBS using ANGSD

###############
#### SETUP ####
###############

rule create_bam_list_byPop_multiInd:
    input:
        rules.create_bam_list_allFinalSamples.output
    output:
        '{0}/bam_lists/by_population/{{popu}}_{{site}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/create_bam_list/{popu}_{site}_bam_list.log'
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

################################
#### SAF AND SFS ESTIMATION ####
################################

rule angsd_saf_likelihood_byPopulation:
    """
    Generate Site Allele Frequency (SAF) likelihood file for each multi-individual population using ANGSD. Uses only 4fold sites.
    """
    input:
        unpack(get_files_for_saf_estimation_byPopulation)
    output:
        saf = '{0}/sfs/1d/by_population/{{popu}}/{{popu}}_{{site}}.saf.gz'.format(ANGSD_DIR),
        saf_idx = '{0}/sfs/1d/by_population/{{popu}}/{{popu}}_{{site}}.saf.idx'.format(ANGSD_DIR),
        saf_pos = '{0}/sfs/1d/by_population/{{popu}}/{{popu}}_{{site}}.saf.pos.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_saf_likelihood_byPopulation/{popu}_{site}_saf.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    params:
        out = '{0}/sfs/1d/by_population/{{popu}}/{{popu}}_{{site}}'.format(ANGSD_DIR),
        min_dp_ind = ANGSD_MIN_DP_IND_SFS
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '6:00:00'
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND*60/100 ));
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 4 \
            -baq 2 \
            -ref {input.ref} \
            -doCounts 1 \
            -setMinDepthInd {params.min_dp_ind} \
            -minInd $MIN_IND \
            -minQ 20 \
            -minMapQ 30 \
            -doSaf 1 \
            -sites {input.sites} \
            -anc {input.ref} \
            -rf {input.chroms} \
            -bam {input.bams} 2> {log}
        """

rule angsd_byPopulation_done:
    input:
        expand(rules.angsd_saf_likelihood_byPopulation.output, popu=POPS_MULTI_IND, site='4fold'),
    output:
        '{0}/angsd_byPopulation.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """
