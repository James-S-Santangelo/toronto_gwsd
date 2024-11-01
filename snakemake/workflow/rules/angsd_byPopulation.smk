# Rules to estimate within-population diversity and pairwise Fst/PBS using ANGSD

###############
#### SETUP ####
###############

rule create_bam_list_byPop_multiInd:
    """
    Create BAM list with all individuals in a population
    """
    input:
        rules.create_bam_list_allFinalSamples.output
    output:
        '{0}/bam_lists/byPopulation/{{popu}}_{{site}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
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
        saf = '{0}/saf/byPopulation/{{popu}}/{{popu}}_{{site}}.saf.gz'.format(ANGSD_DIR),
        saf_idx = '{0}/saf/byPopulation/{{popu}}/{{popu}}_{{site}}.saf.idx'.format(ANGSD_DIR),
        saf_pos = '{0}/saf/byPopulation/{{popu}}/{{popu}}_{{site}}.saf.pos.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_saf_likelihood_byPopulation/{popu}_{site}_saf.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = '{0}/saf/byPopulation/{{popu}}/{{popu}}_{{site}}'.format(ANGSD_DIR),
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


rule angsd_estimate_sfs_byPopulation:
    """
    Estimate folded SFS separately for each population(i.e., 1D SFS) using realSFS. 
    """
    input:
        saf = rules.angsd_saf_likelihood_byPopulation.output.saf_idx,
        sites = rules.convert_sites_for_angsd.output,
        idx = rules.angsd_index_degenerate_sites.output,
    output:
        '{0}/sfs/1dsfs/byPopulation/{{popu}}/{{popu}}_{{site}}.sfs'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_sfs_byHabitat/{popu}_{site}.sfs.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 6
    wildcard_constraints:
        site='4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        realSFS {input.saf} \
            -sites {input.sites} \
            -P {threads} \
            -fold 1 \
            -maxIter 2000 \
            -seed 42 > {output} 2> {log}
        """

rule angsd_estimate_joint_population_sfs:
    """
    Estimated folded, pairwise population SFS using realSFS. Uses 4fold sites.
    """
    input:
        safs = get_population_saf_files,
        sites = rules.convert_sites_for_angsd.output,
        idx = rules.angsd_index_degenerate_sites.output,
    output:
        '{0}/sfs/2dsfs/byPopulation/{{pop_comb}}_{{site}}.2dsfs'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_population_2dsfs/{pop_comb}_{site}.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '01:00:00'
    shell:
        """
        realSFS {input.safs} \
            -sites {input.sites} \
            -maxIter 2000 \
            -seed 42 \
            -fold 1 \
            -P {threads} > {output} 2> {log}
        """

#########################
#### THETA & FST/PBS ####
#########################


rule angsd_estimate_thetas_byPopulation:
    """
    Generate per-site thetas in each population from 1DSFS
    """
    input:
        saf_idx = rules.angsd_saf_likelihood_byPopulation.output.saf_idx,
        sfs = rules.angsd_estimate_sfs_byPopulation.output
    output:
        idx = '{0}/summary_stats/thetas/byPopulation/{{popu}}/{{popu}}_{{site}}.thetas.idx'.format(ANGSD_DIR),
        thet = '{0}/summary_stats/thetas/byPopulation/{{popu}}/{{popu}}_{{site}}.thetas.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_thetas_byPopulation/{site}_{popu}_thetas.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 4
    wildcard_constraints:
        site='4fold'
    params:
        out = '{0}/summary_stats/thetas/byPopulation/{{popu}}/{{popu}}_{{site}}'.format(ANGSD_DIR)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        realSFS saf2theta {input.saf_idx} \
            -P {threads} \
            -fold 1 \
            -sfs {input.sfs} \
            -outname {params.out} 2> {log}
        """

rule angsd_diversity_neutrality_stats_byPopulation:
    """
    Estimate pi, Waterson's theta, Tajima's D, etc. in each population 
    """
    input:
        rules.angsd_estimate_thetas_byPopulation.output.idx
    output:
       '{0}/summary_stats/thetas/byPopulation/{{popu}}/{{popu}}_{{site}}.thetas.idx.pestPG'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_diversity_neutrality_stats_byPopulation/{site}_{popu}_diversity_neutrality.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    wildcard_constraints:
        site='4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        thetaStat do_stat {input} 2> {log}
        """

rule angsd_population_fst_index:
    """
    Estimate per-site alphas (numerator) and betas (denominator) for population Fst/pbs estimation
    """
    input: 
        saf_files = get_population_saf_files, 
        sfs = rules.angsd_estimate_joint_population_sfs.output
    output:
        fst = temp('{0}/summary_stats/hudson_fst/byPopulation/{{pop_comb}}_{{site}}.fst.gz'.format(ANGSD_DIR)),
        idx = temp('{0}/summary_stats/hudson_fst/byPopulation/{{pop_comb}}_{{site}}.fst.idx'.format(ANGSD_DIR))
    log: LOG_DIR + '/angsd_population_fst_index/{pop_comb}_{site}_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 4
    resources:
        mem_mb = 4000,
        time = '02:00:00'
    params:
        fstout = '{0}/summary_stats/hudson_fst/byPopulation/{{pop_comb}}_{{site}}'.format(ANGSD_DIR)
    shell:
        """
        realSFS fst index {input.saf_files} \
            -sfs {input.sfs} \
            -fold 1 \
            -P {threads} \
            -whichFst 1 \
            -fstout {params.fstout} 2> {log}
        """

rule angsd_population_fst_readable:
    """
    Create readable Fst files. Required due to format of realSFS fst index output files. 
    """
    input:
        idx = rules.angsd_population_fst_index.output.idx,
        fst = rules.angsd_population_fst_index.output.fst
    output:
        '{0}/summary_stats/hudson_fst/byPopulation/{{pop_comb}}_{{site}}_readable.fst'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_population_fst_readable/{pop_comb}_{site}_readable.log'
    resources:
        mem_mb = 1000,
        time = '01:00:00'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    shell:
        """
        realSFS fst print {input.idx} > {output} 2> {log}
        """

##############
#### POST ####
##############

rule angsd_byPopulation_done:
    """
    Create empty flag file signaling completion of per-population ANGSD stats
    """
    input:
        expand(rules.angsd_diversity_neutrality_stats_byPopulation.output, popu=POPS_MULTI_IND, site='4fold'),
        expand(rules.angsd_population_fst_readable.output, pop_comb=POP_COMB_MULTI_IND, site='4fold')
    output:
        '{0}/angsd_byPopulation.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """
