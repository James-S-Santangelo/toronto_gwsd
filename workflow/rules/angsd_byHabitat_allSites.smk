# Rules to estimate SAF/SFS by habitat using all sites. Used for detecting selective sweeps

###############
#### SETUP ####
###############

rule create_bam_lists_allFinalSamples_allSites:
    input:
        bams = expand(rules.samtools_markdup.output.bam, sample=SAMPLES)
    output:
        '{0}/bam_lists/allFinalSamples_allSites_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/create_bam_list/allFinalSamples_allSites_bam_list.log'
    run:
        import os
        with open(output[0], 'w') as f:
            for bam in input.bams:
                search = re.search('^(s_\d+_\d+)(?=_\w)', os.path.basename(bam))
                sample = search.group(1) 
                if sample in FINAL_SAMPLES:
                    f.write('{0}\n'.format(bam))


rule create_bam_list_byHabitat_allSites:
    """
    Create text file with paths to BAMS for urban, rural, and suburban samples. BAMs contain all reads genome-wide
    """
    input:
        samples = config['samples'],
        bams = rules.create_bam_lists_allFinalSamples_allSites.output
    output:
        '{0}/bam_lists/{{habitat}}_allSites_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/create_bam_list/{habitat}_allSites_bams.log'
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
                    search = re.search('^(s_\d+_\d+)(?=_\w)', os.path.basename(bam))
                    sample = search.group(1) 
                    if sample in samples_habitat:
                        f.write('{0}'.format(bam))
        except:
            logging.exception("An error occured!") 
            raise

################################
#### SAF AND SFS ESTIMATION ####
################################

rule angsd_saf_likelihood_byHabitat_allSites:
    """
    Generate Site Allele Frequency (SAF) likelihood file for each habitat using ANGSD. Uses only 4fold sites.
    """
    input:
        bams = rules.create_bam_list_byHabitat_allSites.output,
        ref = rules.unzip_reference.output
    output:
        saf = '{0}/sfs/{{habitat}}/allSites/{{chrom}}/{{chrom}}_{{habitat}}_allSites.saf.gz'.format(ANGSD_DIR),
        saf_idx = '{0}/sfs/{{habitat}}/allSites/{{chrom}}/{{chrom}}_{{habitat}}_allSites.saf.idx'.format(ANGSD_DIR),
        saf_pos = '{0}/sfs/{{habitat}}/allSites/{{chrom}}/{{chrom}}_{{habitat}}_allSites.saf.pos.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_saf_likelihood_byHabitat_allSites/{chrom}_{habitat}_allSites_saf.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    params:
        out = '{0}/sfs/{{habitat}}/allSites/{{chrom}}/{{chrom}}_{{habitat}}_allSites'.format(ANGSD_DIR),
        min_dp_ind = ANGSD_MIN_DP_IND_SFS
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
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
            -anc {input.ref} \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """

rule angsd_estimate_joint_habitat_sfs_allSites:
    """
    Estimated folded, two-dimensional urban-rural SFS for each city using realSFS. Uses all sites4
    """
    input:
        safs = get_habitat_saf_files_allSites
    output:
        '{0}/sfs/2dsfs/byHabitat/allSites/{{chrom}}/{{chrom}}_allSites_{{hab_comb}}.2dsfs'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_habitat_2dsfs_allSites/{chrom}_allSites_{hab_comb}.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 30000,
        time = '06:00:00'
    shell:
        """
        realSFS {input.safs} \
            -maxIter 2000 \
            -seed 42 \
            -fold 1 \
            -P {threads} > {output} 2> {log}
        """

rule angsd_estimate_sfs_byHabitat_allSites:
    """
    Estimate folded SFS separately for each habitat (i.e., 1D SFS) using realSFS. 
    """
    input:
        saf = rules.angsd_saf_likelihood_byHabitat_allSites.output.saf_idx
    output:
        '{0}/sfs/1dsfs/byHabitat/allSites/{{chrom}}/{{chrom}}_allSites_{{habitat}}.sfs'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_sfs_byHabitat_allSites/{chrom}_allSites_{habitat}_sfs.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '01:00:00'
    shell:
        """
        realSFS {input.saf} \
            -P {threads} \
            -fold 1 \
            -maxIter 2000 \
            -seed 42 > {output} 2> {log}
        """

########################
#### FST AND THETAS ####
########################

rule angsd_habitat_fst_index_allSites:
    """
    Estimate per-site alphas (numerator) and betas (denominator) for Fst estimation
    """
    input: 
        saf_idx = get_habitat_saf_files_allSites,
        joint_sfs = rules.angsd_estimate_joint_habitat_sfs_allSites.output
    output:
        fst = '{0}/summary_stats/hudson_fst/byHabitat/allSites/{{chrom}}/{{chrom}}_allSites_{{hab_comb}}.fst.gz'.format(ANGSD_DIR),
        idx = '{0}/summary_stats/hudson_fst/byHabitat/allSites/{{chrom}}/{{chrom}}_allSites_{{hab_comb}}.fst.idx'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_habitat_fst_index_allSites/{chrom}_allSites_{hab_comb}_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 4
    resources:
        mem_mb = 4000,
        time = '02:00:00'
    params:
        fstout = '{0}/summary_stats/hudson_fst/byHabitat/allSites/{{chrom}}/{{chrom}}_allSites_{{hab_comb}}'.format(ANGSD_DIR)
    shell:
        """
        realSFS fst index {input.saf_idx} \
            -sfs {input.joint_sfs} \
            -fold 1 \
            -P {threads} \
            -whichFst 1 \
            -fstout {params.fstout} 2> {log}
        """

rule angsd_fst_allSites_readable:
    """
    Create readable Fst files. Required due to format of realSFS fst index output files. 
    """
    input:
        rules.angsd_habitat_fst_index_allSites.output.idx
    output:
        '{0}/summary_stats/hudson_fst/byHabitat/allSites/{{chrom}}/{{chrom}}_allSites_{{hab_comb}}_readable.fst'.format(ANGSD_DIR)
    log: 'logs/angsd_fst_allSites_readable/{chrom}_{hab_comb}_readable_fst.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    shell:
        """
        realSFS fst print {input} > {output} 2> {log}
        """

rule angsd_estimate_thetas_byHabitat_allSites:
    """
    Generate per-site thetas in each habitat from 1DSFS
    """
    input:
        saf_idx = rules.angsd_saf_likelihood_byHabitat_allSites.output.saf_idx,
        sfs = rules.angsd_estimate_sfs_byHabitat_allSites.output
    output:
        idx = '{0}/summary_stats/thetas/byHabitat/allSites/{{chrom}}/{{chrom}}_allSites_{{habitat}}.thetas.idx'.format(ANGSD_DIR),
        thet = '{0}/summary_stats/thetas/byHabitat/allSites/{{chrom}}/{{chrom}}_allSites_{{habitat}}.thetas.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_thetas_byHabitat_allSites/{chrom}_allSites_{habitat}_thetas.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 4
    params:
        out = '{0}/summary_stats/thetas/byHabitat/allSites/{{chrom}}/{{chrom}}_allSites_{{habitat}}'.format(ANGSD_DIR)
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

rule angsd_thetas_allSites_readable:
    """
    Create readable Fst files. Required due to format of realSFS fst index output files. 
    """
    input:
        rules.angsd_estimate_thetas_byHabitat_allSites.output.idx
    output:
        '{0}/summary_stats/thetas/byHabitat/allSites/{{chrom}}/{{chrom}}_allSites_{{habitat}}_readable.thetas'.format(ANGSD_DIR)
    log: 'logs/angsd_thetas_allSites_readable/{chrom}_{habitat}_readable_thetas.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    shell:
        """
        thetaStat print {input} > {output} 2> {log}
        """

###########################
#### WINDOWED ANALYSES ####
###########################

rule windowed_theta:
    input:
        rules.angsd_estimate_thetas_byHabitat_allSites.output.idx
    output:
        "{0}/summary_stats/thetas/byHabitat/allSites/{{chrom}}/{{chrom}}_allSites_{{habitat}}_windowedThetas50.gz.pestPG".format(ANGSD_DIR)
    log: LOG_DIR + '/windowed_theta/{chrom}_{habitat}_windowTheta.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    params:
        out = "{0}/summary_stats/thetas/byHabitat/allSites/{{chrom}}/{{chrom}}_allSites_{{habitat}}_windowedThetas50.gz".format(ANGSD_DIR),
        win = 50000,
        step = 50000
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        thetaStat do_stat {input} -win {params.win} -step {params.step} -outnames {params.out} 2> {log}
        """

rule windowed_fst:
    input:
        rules.angsd_habitat_fst_index_allSites.output.idx
    output:
        "{0}/summary_stats/hudson_fst/byHabitat/allSites/{{chrom}}/{{chrom}}_allSites_{{hab_comb}}_windowed50.fst".format(ANGSD_DIR)
    log: LOG_DIR + '/windowed_fst/{chrom}_{hab_comb}_windowedFst.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    params:
        win = 50000,
        step = 50000
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        realSFS fst stats2 {input} -win {params.win} -step {params.step} > {output} 2> {log}
        """

##############
#### POST ####
##############

rule angsd_byHabitat_allSites_done:
    """
    Generate empty flag file signalling successful completion of SFS and summary stat for habitats
    """
    input:
        expand(rules.windowed_fst.output, chrom=CHROMOSOMES, hab_comb=HABITAT_COMBOS),
        expand(rules.angsd_fst_allSites_readable.output, chrom=CHROMOSOMES, hab_comb=HABITAT_COMBOS),
        expand(rules.windowed_theta.output, chrom=CHROMOSOMES, habitat=HABITATS),
        expand(rules.angsd_thetas_allSites_readable.output, chrom=CHROMOSOMES, habitat=HABITATS)
    output:
        '{0}/angsd_byHabitat_allSites50.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """
