# Rules for estimating SFS (1D and 2D), summary stats, and GLs for urban, rural, and suburban habitats within cities

###############
#### SETUP ####
###############

rule create_bam_list_byHabitat:
    """
    Create text file with paths to BAMS for urban, rural, and suburban samples. BAMs are subsetted around 4fold sites
    """
    input:
        samples = config['samples'],
        bams = rules.create_bam_list_allFinalSamples.output
    output:
        '{0}/bam_lists/{{habitat}}_{{site}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/create_bam_list/{habitat}_{{site}}_bams.log'
    run:
        import os
        import logging
        import pandas as pd
        logging.basicConfig(filename=log[0], level=logging.DEBUG)
        try:
            df = pd.read_table(input.samples, sep = '\t')
            df_sub = df[(df['Habitat'] == wildcards.habitat)]
            samples_habitat = df_sub['Sample'].tolist()
            bams = open(input.bams[0], 'r').readlines()
            with open(output[0], 'w') as f:
                for bam in bams:
                    search = re.search('^(.+)(?=_\w)', os.path.basename(bam))
                    sample = search.group(1) 
                    if sample in samples_habitat:
                        f.write('{0}'.format(bam))
        except:
            logging.exception("An error occured!") 
            raise

################################
#### SAF AND SFS ESTIMATION ####
################################

rule angsd_saf_likelihood_byHabitat:
    """
    Generate Site Allele Frequency (SAF) likelihood file for each habitat using ANGSD. Uses only 4fold sites.
    """
    input:
        unpack(get_files_for_saf_estimation_byHabitat)
    output:
        saf = '{0}/saf/{{habitat}}/{{habitat}}_{{site}}.saf.gz'.format(ANGSD_DIR),
        saf_idx = '{0}/saf/{{habitat}}/{{habitat}}_{{site}}.saf.idx'.format(ANGSD_DIR),
        saf_pos = '{0}/saf/{{habitat}}/{{habitat}}_{{site}}.saf.pos.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_saf_likelihood_byHabitat/{habitat}_{site}_saf.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = '{0}/saf/{{habitat}}/{{habitat}}_{{site}}'.format(ANGSD_DIR),
        min_dp_ind = ANGSD_MIN_DP_IND_SFS
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '12:00:00'
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

rule angsd_estimate_joint_habitat_sfs:
    """
    Estimated folded, two-dimensional urban-rural SFS for each city using realSFS. Uses 4fold sites.
    """
    input:
        safs = get_habitat_saf_files,
        sites = rules.convert_sites_for_angsd.output,
        idx = rules.angsd_index_degenerate_sites.output,
    output:
        '{0}/sfs/2dsfs/byHabitat/{{site}}_{{hab_comb}}.2dsfs'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_habitat_2dsfs/{site}_{hab_comb}.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '03:00:00'
    shell:
        """
        realSFS {input.safs} \
            -sites {input.sites} \
            -tole 1e-6 \
            -maxIter 2000 \
            -seed 42 \
            -fold 1 \
            -P {threads} > {output} 2> {log}
        """

rule angsd_estimate_sfs_byHabitat:
    """
    Estimate folded SFS separately for each habitat (i.e., 1D SFS) using realSFS. 
    """
    input:
        saf = rules.angsd_saf_likelihood_byHabitat.output.saf_idx,
        sites = rules.convert_sites_for_angsd.output,
        idx = rules.angsd_index_degenerate_sites.output,
    output:
        '{0}/sfs/1dsfs/byHabitat/{{site}}_{{habitat}}.sfs'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_sfs_byHabitat/{site}_{habitat}_sfs.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 6
    wildcard_constraints:
        site='4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '03:00:00'
    shell:
        """
        realSFS {input.saf} \
            -sites {input.sites} \
            -P {threads} \
            -tole 1e-6 \
            -fold 1 \
            -maxIter 2000 \
            -seed 42 > {output} 2> {log}
        """

########################
#### FST AND THETAS ####
########################

rule angsd_habitat_fst_index:
    """
    Estimate per-site alphas (numerator) and betas (denominator) for Fst estimation
    """
    input: 
        saf_idx = get_habitat_saf_files,
        joint_sfs = rules.angsd_estimate_joint_habitat_sfs.output
    output:
        fst = temp('{0}/summary_stats/hudson_fst/byHabitat/{{site}}_{{hab_comb}}.fst.gz'.format(ANGSD_DIR)),
        idx = temp('{0}/summary_stats/hudson_fst/byHabitat/{{site}}_{{hab_comb}}.fst.idx'.format(ANGSD_DIR))
    log: LOG_DIR + '/angsd_habitat_fst_index/{site}_{hab_comb}_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 4
    resources:
        mem_mb = 4000,
        time = '02:00:00'
    params:
        fstout = '{0}/summary_stats/hudson_fst/byHabitat/{{site}}_{{hab_comb}}'.format(ANGSD_DIR)
    shell:
        """
        realSFS fst index {input.saf_idx} \
            -sfs {input.joint_sfs} \
            -fold 1 \
            -P {threads} \
            -whichFst 1 \
            -fstout {params.fstout} 2> {log}
        """

rule angsd_habitat_fst_readable:
    """
    Create readable Fst files. Required due to format of realSFS fst index output files. 
    """
    input:
        idx = rules.angsd_habitat_fst_index.output.idx,
        fst = rules.angsd_habitat_fst_index.output.fst
    output:
        '{0}/summary_stats/hudson_fst/byHabitat/{{site}}_{{hab_comb}}_readable.fst'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_habitat_fst_readable/{site}_{hab_comb}_readable.log'
    resources:
        mem_mb = 4000,
        time = '01:00:00'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    shell:
        """
        realSFS fst print {input.idx} > {output} 2> {log}
        """

rule angsd_estimate_thetas_byHabitat:
    """
    Generate per-site thetas in each habitat from 1DSFS
    """
    input:
        saf_idx = rules.angsd_saf_likelihood_byHabitat.output.saf_idx,
        sfs = rules.angsd_estimate_sfs_byHabitat.output
    output:
        idx = temp('{0}/summary_stats/thetas/byHabitat/{{site}}_{{habitat}}.thetas.idx'.format(ANGSD_DIR)),
        thet = temp('{0}/summary_stats/thetas/byHabitat/{{site}}_{{habitat}}.thetas.gz'.format(ANGSD_DIR))
    log: LOG_DIR + '/angsd_estimate_thetas_byHabitat/{site}_{habitat}_thetas.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 4
    wildcard_constraints:
        site='4fold'
    params:
        out = '{0}/summary_stats/thetas/byHabitat/{{site}}_{{habitat}}'.format(ANGSD_DIR)
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

rule angsd_diversity_neutrality_stats_byHabitat:
    """
    Estimate pi, Waterson's theta, Tajima's D, etc. in each habitat
    """
    input:
        idx = rules.angsd_estimate_thetas_byHabitat.output.idx,
        thet = rules.angsd_estimate_thetas_byHabitat.output.thet,
    output:
       '{0}/summary_stats/thetas/byHabitat/{{site}}_{{habitat}}.thetas.idx.pestPG'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_diversity_neutrality_stats_byHabitat/{site}_{habitat}_diversity_neutrality.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    wildcard_constraints:
        site='4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        thetaStat do_stat {input.idx} 2> {log}
        """

##############
#### POST ####
##############

rule angsd_byHabitat_done:
    """
    Generate empty flag file signalling successful completion of SFS and summary stat for habitats
    """
    input:
        expand(rules.angsd_habitat_fst_readable.output, site=['4fold'], hab_comb=HABITAT_COMBOS),
        expand(rules.angsd_diversity_neutrality_stats_byHabitat.output, site=['4fold'], habitat=HABITATS)
    output:
        '{0}/angsd_byHabitat.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """
