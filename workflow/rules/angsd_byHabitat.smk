# Rules for estimating SFS (1D and 2D), summary stats, and GLs for urban, rural, and suburban habitats within cities

###############################
#### SFS AND SUMMARY STATS ####
###############################

rule create_bam_list_byHabitat:
    """
    Create text file with paths to BAMS for urban, rural, and suburban samples
    """
    input:
        rules.create_bam_list_highQualSamples.output
    output:
        '{0}/bam_lists/{{habitat}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/create_bam_list/{habitat}_bams.log'
    wildcard_constraints:
        habitat='Urban|Rural|Suburban'
    run:
        import os
        import pandas as pd
        df = pd.read_table(config['samples'], sep = '\t')
        df_sub = df[(df['Habitat'] == wildcards.habitat)]
        samples_habitat = df_sub['Sample'].tolist()
        bams = open(input[0], 'r').readlines()
        with open(output[0], 'w') as f:
            for bam in bams:
                sample = os.path.basename(bam).split('_merged')[0]
                if sample in samples_habitat:
                    f.write('{0}'.format(bam))

rule angsd_saf_likelihood_byHabitat:
    """
    Generate Site Allele Frequency (SAF) likelihood file for each habitat using ANGSD. Uses only 4fold sites.
    """
    input:
        unpack(get_files_for_saf_estimation_byHabitat)
    output:
        saf = temp('{0}/sfs/{{habitat}}/{{chrom}}/{{chrom}}_{{habitat}}_allSites.saf.gz'.format(ANGSD_DIR)),
        saf_idx = temp('{0}/sfs/{{habitat}}/{{chrom}}/{{chrom}}_{{habitat}}_allSites.saf.idx'.format(ANGSD_DIR)),
        saf_pos = temp('{0}/sfs/{{habitat}}/{{chrom}}/{{chrom}}_{{habitat}}_allSites.saf.pos.gz'.format(ANGSD_DIR))
    log: 'logs/angsd_saf_likelihood_byHabitat/{chrom}_{habitat}_allSites_saf.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    params:
        out = '{0}/sfs/{{habitat}}/{{chrom}}/{{chrom}}_{{habitat}}_allSites'.format(ANGSD_DIR),
        min_dp_ind = ANGSD_MIN_DP_IND_SFS
    threads: 12
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 12000,
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
            -anc {input.ref} \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """

rule angsd_estimate_joint_habitat_sfs:
    """
    Estimated folded, two-dimensional urban-rural SFS for each city using realSFS. Uses 4fold sites.
    """
    input:
        safs = get_habitat_saf_files,
        sites = rules.split_angsd_sites_byChrom.output,
        sites_idx = rules.angsd_index_sites.output
    output:
        '{0}/sfs/2dsfs/habitat/{{chrom}}/{{chrom}}_Toronto_2dsfs_{{site}}_{{hab_comb}}.2dsfs'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_habitat_2dsfs/{chrom}_Toronto_{site}_{hab_comb}.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 6
    wildcard_constraints:
        site='4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '01:00:00'
    shell:
        """
        realSFS {input.safs} -sites {input.sites} -maxIter 2000 -seed 42 -fold 1 -P {threads} > {output} 2> {log}
        """

rule angsd_habitat_fst_index:
    """
    Estimate per-site alphas (numerator) and betas (denominator) for Fst estimation
    """
    input: 
        saf_idx = get_habitat_saf_files,
        joint_sfs = rules.angsd_estimate_joint_habitat_sfs.output
    output:
        fst = '{0}/summary_stats/hudson_fst/habitat/{{chrom}}/{{chrom}}_Toronto_{{site}}_{{hab_comb}}.fst.gz'.format(ANGSD_DIR),
        idx = '{0}/summary_stats/hudson_fst/habitat/{{chrom}}/{{chrom}}_Toronto_{{site}}_{{hab_comb}}.fst.idx'.format(ANGSD_DIR)
    log: 'logs/angsd_habitat_fst_index/{chrom}_Toronto_{site}_{hab_comb}_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    wildcard_constraints:
        site='4fold'
    threads: 4
    resources:
        mem_mb = 4000,
        time = '01:00:00'
    params:
        fstout = '{0}/summary_stats/hudson_fst/habitat/{{chrom}}/{{chrom}}_Toronto_{{site}}_{{hab_comb}}'.format(ANGSD_DIR)
    shell:
        """
        realSFS fst index {input.saf_idx} -sfs {input.joint_sfs} -fold 1 -P {threads} -whichFst 1 -fstout {params.fstout} 2> {log}
        """

rule angsd_habitat_fst_readable:
    """
    Create readable Fst files. Required due to format of realSFS fst index output files. 
    """
    input:
        rules.angsd_habitat_fst_index.output.idx
    output:
        '{0}/summary_stats/hudson_fst/habitat/{{chrom}}/{{chrom}}_Toronto_{{site}}_{{hab_comb}}_readable.fst'.format(ANGSD_DIR)
    log: 'logs/angsd_habitat_fst_readable/{chrom}_Toronto_{site}_{hab_comb}_readable.log'
    wildcard_constraints:
        site='4fold'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    shell:
        """
        realSFS fst print {input} > {output} 2> {log}
        """

rule angsd_estimate_sfs_byHabitat:
    """
    Estimate folded SFS separately for each habitat (i.e., 1D SFS) using realSFS. 
    """
    input:
        saf = rules.angsd_saf_likelihood_byHabitat.output.saf_idx,
        sites = rules.split_angsd_sites_byChrom.output,
        sides_idx = rules.angsd_index_sites.output
    output:
        '{0}/sfs/1dsfs/habitat/{{chrom}}/{{chrom}}_Toronto_{{site}}_{{habitat}}.sfs'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_sfs_byHabitat/{chrom}_Toronto_{site}_{habitat}_sfs.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 6
    wildcard_constraints:
        site='4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '01:00:00'
    shell:
        """
        realSFS {input.saf} -sites {input.sites} -P {threads} -fold 1 -maxIter 2000 -seed 42 > {output} 2> {log}
        """

rule angsd_estimate_thetas_byHabitat:
    """
    Generate per-site thetas in each habitat from 1DSFS
    """
    input:
        saf_idx = rules.angsd_saf_likelihood_byHabitat.output.saf_idx,
        sfs = rules.angsd_estimate_sfs_byHabitat.output
    output:
        idx = '{0}/summary_stats/thetas/habitat/{{chrom}}/{{chrom}}_Toronto_{{site}}_{{habitat}}.thetas.idx'.format(ANGSD_DIR),
        thet = '{0}/summary_stats/thetas/habitat/{{chrom}}/{{chrom}}_Toronto_{{site}}_{{habitat}}.thetas.gz'.format(ANGSD_DIR)
    log: 'logs/angsd_estimate_thetas_byHabitat/{chrom}_{site}_{habitat}_thetas.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 4
    wildcard_constraints:
        site='4fold'
    params:
        out = '{0}/summary_stats/thetas/habitat/{{chrom}}/{{chrom}}_Toronto_{{site}}_{{habitat}}'.format(ANGSD_DIR)
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
        rules.angsd_estimate_thetas_byHabitat.output.idx
    output:
       '{0}/summary_stats/thetas/habitat/{{chrom}}/{{chrom}}_Toronto_{{site}}_{{habitat}}.thetas.idx.pestPG'.format(ANGSD_DIR)
    log: 'logs/angsd_diversity_neutrality_stats_byHabitat/{chrom}_Toronto_{site}_{habitat}_diversity_neutrality.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    wildcard_constraints:
        site='4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        thetaStat do_stat {input} 2> {log}
        """



##############
###POST ####
##############

rule angsd_byHabitat_done:
    """
    Generate empty flag file signalling successful completion of SFS and summary stat for habitats
    """
    input:
        expand(rules.angsd_habitat_fst_readable.output, site=['4fold'], hab_comb=HABITAT_COMBOS, chrom=CHROMOSOMES),
        expand(rules.angsd_diversity_neutrality_stats.output, site=['4fold'], hab_comb=HABITATS, chrom=CHROMOSOMES)
    output:
        '{0}/angsd_byHabitat.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """
