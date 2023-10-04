# Rules to scan for signatures of selective sweeps

###############
#### SETUP ####
###############

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


rule samples_byHabitat:
    input:
        samples = config['samples']
    output:
        '{0}/selscan/{{chrom}}/{{chrom}}_{{habitat}}.samples'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/samples_byHabitat/{chrom}_{habitat}_samples.log'
    run:
        import os
        import logging
        import pandas as pd
        logging.basicConfig(filename=log[0], level=logging.DEBUG)
        try:
            df = pd.read_table(input.samples, sep = '\t')
            df_sub = df[(df['Habitat'] == wildcards.habitat)]
            samples_habitat = df_sub['Sample'].tolist()
            samples_toWrite = []
            outname = wildcards.chrom + '_' + wildcards.habitat + '_whatshapPhased_shapeitPhased'
            for sample in samples_habitat:
                if sample in FINAL_SAMPLES:
                    samples_toWrite.append(sample)
            with open(output[0], 'w') as f:
                for sample in samples_toWrite:
                    if sample == samples_toWrite[-1]:
                        f.write('{0}\t-\t{1}'.format(sample, outname))
                    else:
                        f.write('{0},'.format(sample))
        except:
            logging.exception("An error occured!") 
            raise
       
rule bcftools_splitVCF_byHabitat:
    input:
        vcf = rules.shapeit_phase.output.vcf,
        samples = lambda w: expand(rules.samples_byHabitat.output, habitat = w.habitat, chrom = w.chrom)
    output:
        vcf = '{0}/vcf/{{chrom}}/{{chrom}}_{{habitat}}_whatshapPhased_shapeitPhased.vcf.gz'.format(FREEBAYES_DIR),
        idx = '{0}/vcf/{{chrom}}/{{chrom}}_{{habitat}}_whatshapPhased_shapeitPhased.vcf.gz.tbi'.format(FREEBAYES_DIR)
    log: LOG_DIR + '/bcftools_splitVCF_byHabitat/{chrom}_{habitat}_split.log'
    conda: '../envs/sweeps.yaml',
    params:
        out =  '{0}/vcf/{{chrom}}/'.format(FREEBAYES_DIR)
    shell:
        """
        ( bcftools +split {input.vcf} \
            --samples-file {input.samples} \
            --output-type z \
            --output {params.out} && tabix {output.vcf} ) 2> {log}
        """

rule genMap_toPlinkFormat:
    input:
        rules.split_genMap.output
    output:
        '{0}/{{chrom}}.map'.format(GENMAP_RESULTS_DIR)
    shell:
        """
        awk -F"\t" '{{print $2 " " $2"_"$1 " " $3 " " $1}}' {input} | tail -n +2 > {output}
        """

##############################
#### ANGSD FST AND THETAS ####
##############################

rule angsd_saf_likelihood_byHabitat_allSites:
    """
    Generate Site Allele Frequency (SAF) likelihood file for each habitat using ANGSD. Uses only 4fold sites.
    """
    input:
        bams = rules.create_bam_list_byHabitat_allSites.output,
        ref = REFERENCE_GENOME
    output:
        saf = '{0}/sfs/{{habitat}}/allSites/{{chrom}}/{{chrom}}_{{habitat}}_allSites.saf.gz'.format(ANGSD_DIR),
        saf_idx = '{0}/sfs/{{habitat}}/allSites/{{chrom}}/{{chrom}}_{{habitat}}_allSites.saf.idx'.format(ANGSD_DIR),
        saf_pos = '{0}/sfs/{{habitat}}/allSites/{{chrom}}/{{chrom}}_{{habitat}}_allSites.saf.pos.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_saf_likelihood_byHabitat_allSites/{chrom}_{habitat}_allSites_saf.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
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
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 8 
    resources:
        mem_mb = 10000,
        time = lambda wildcards, attempt: str(attempt * 24) + ":00:00"
    shell:
        """
        realSFS {input.safs} \
            -maxIter 2000 \
            -seed 42 \
            -tole 1e-6 \
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
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    threads: 8
    resources:
        mem_mb = 10000,
        time = lambda wildcards, attempt: str(attempt * 4) + ":00:00"
    shell:
        """
        realSFS {input.saf} \
            -P {threads} \
            -fold 1 \
            -tole 1e-6 \
            -maxIter 2000 \
            -seed 42 > {output} 2> {log}
        """

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
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
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
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
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
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
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
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    shell:
        """
        thetaStat print {input} > {output} 2> {log}
        """

rule windowed_theta:
    input:
        rules.angsd_estimate_thetas_byHabitat_allSites.output.idx
    output:
        "{0}/summary_stats/thetas/byHabitat/allSites/{{chrom}}/{{chrom}}_allSites_{{habitat}}_windowedThetas.gz.pestPG".format(ANGSD_DIR)
    log: LOG_DIR + '/windowed_theta/{chrom}_{habitat}_windowTheta.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = "{0}/summary_stats/thetas/byHabitat/allSites/{{chrom}}/{{chrom}}_allSites_{{habitat}}_windowedThetas.gz".format(ANGSD_DIR),
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
        "{0}/summary_stats/hudson_fst/byHabitat/allSites/{{chrom}}/{{chrom}}_allSites_{{hab_comb}}_windowed.fst".format(ANGSD_DIR)
    log: LOG_DIR + '/windowed_fst/{chrom}_{hab_comb}_windowedFst.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
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

################
#### XP-NSL ####
################

rule selscan_xpnsl:
    input:
        unpack(selscan_xpnsl_input)
    output:
        temp('{0}/xpnsl/{{chrom}}/{{chrom}}_{{hab_comb}}.xpnsl.out'.format(SWEEPS_DIR))
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    log: LOG_DIR + '/selscan_xpnsl/{chrom}_{hab_comb}_xpnsl.log'
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    threads: 2
    params:
        out = '{0}/xpnsl/{{chrom}}/{{chrom}}_{{hab_comb}}'.format(SWEEPS_DIR) 
    shell:
        """
        selscan --xpnsl \
            --vcf {input.vcf} \
            --vcf-ref {input.vcf_ref} \
            --threads {threads} \
            --out {params.out} 2> {log}
        """

rule norm_xpnsl:
    input:
        lambda w: expand(rules.selscan_xpnsl.output, chrom=CHROMOSOMES, hab_comb=w.hab_comb)
    output:
        expand(f'{SWEEPS_DIR}/xpnsl/{{chrom}}/{{chrom}}_{{hab_comb}}.xpnsl.out.norm', chrom=CHROMOSOMES, allow_missing=True)
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    shell:
        """
        norm --xpnsl --qbins 10 --files {input} 
        """

#############################
#### XP-nSL PERMUTATIONS ####
#############################

rule permuted_samples_byHabitat:
    input:
        samples = config['samples']
    output:
        urb = expand('{0}/selscan/permuted/{{chrom}}/{{chrom}}_Urban_{{n}}.samples'.format(PROGRAM_RESOURCE_DIR), chrom=CHROMOSOMES, allow_missing=True),
        rur = expand('{0}/selscan/permuted/{{chrom}}/{{chrom}}_Rural_{{n}}.samples'.format(PROGRAM_RESOURCE_DIR), chrom=CHROMOSOMES, allow_missing=True)
    log: LOG_DIR + '/permuted_samples_byHabitat/samples_{{n}}.log'
    params:
        samples = FINAL_SAMPLES,
        chroms = CHROMOSOMES
    run:
        import os
        import logging
        import random
        logging.basicConfig(filename=log[0], level=logging.DEBUG)
        try:
            df = pd.read_table(input.samples, sep = '\t')
            samples = df[(df['Habitat'].isin(["Rural", "Urban"]))]["Sample"].tolist()
            samples_toKeep = [x for x in samples if x in params.samples]
            rur_samples = random.sample(samples_toKeep, k=41)
            urb_samples = [x for x in samples_toKeep if x not in rur_samples]
            for chr in params.chroms:
                urban_out = [x for x in output.urb if chr in x][0]
                with open(urban_out, 'w') as f:
                    outname = f"{chr}_Urban_{wildcards.n}_permuted" 
                    for sample in urb_samples:
                        if sample == urb_samples[-1]:
                            f.write('{0}\t-\t{1}'.format(sample, outname))
                        else:
                            f.write('{0},'.format(sample))
                rural_out = [x for x in output.rur if chr in x][0]
                with open(rural_out, 'w') as f:
                    outname = f"{chr}_Rural_{wildcards.n}_permuted" 
                    for sample in rur_samples:
                        if sample == rur_samples[-1]:
                            f.write('{0}\t-\t{1}'.format(sample, outname))
                        else:
                            f.write('{0},'.format(sample))
        except:
            logging.exception("An error occured!") 
            raise
       
rule bcftools_splitVCF_byHabitat_permuted:
    input:
        vcf = rules.shapeit_phase.output.vcf,
        samples = bcftools_splitVCF_permuted_input 
    output:
        vcf = temp('{0}/vcf/permuted/{{chrom}}/{{chrom}}_{{habitat}}_{{n}}_permuted.vcf.gz'.format(FREEBAYES_DIR)),
        idx = temp('{0}/vcf/permuted/{{chrom}}/{{chrom}}_{{habitat}}_{{n}}_permuted.vcf.gz.tbi'.format(FREEBAYES_DIR))
    log: LOG_DIR + '/bcftools_splitVCF_byHabitat_permuted/{chrom}_{habitat}_{n}_split.log'
    conda: '../envs/sweeps.yaml',
    params:
        out =  '{0}/vcf/permuted/{{chrom}}/'.format(FREEBAYES_DIR)
    shell:
        """
        ( bcftools +split {input.vcf} \
            --samples-file {input.samples} \
            --output-type z \
            --output {params.out} && tabix {output.vcf} ) 2> {log}
        """

rule selscan_xpnsl_permuted:
    input:
        vcf_ref = lambda w: expand(rules.bcftools_splitVCF_byHabitat_permuted.output.vcf, chrom=w.chrom, habitat="Rural", n=w.n),
        vcf = lambda w: expand(rules.bcftools_splitVCF_byHabitat_permuted.output.vcf, chrom=w.chrom, habitat="Urban", n=w.n)
    output:
        temp('{0}/xpnsl/permuted/{{chrom}}/{{chrom}}_{{hab_comb}}_{{n}}_permuted.xpnsl.out'.format(SWEEPS_DIR))
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    log: LOG_DIR + '/selscan_xpnsl_permuted/{chrom}_{hab_comb}_{{n}}_xpnsl.log'
    params:
        out = '{0}/xpnsl/permuted/{{chrom}}/{{chrom}}_{{hab_comb}}_{{n}}_permuted'.format(SWEEPS_DIR) 
    shell:
        """
        selscan --xpnsl \
            --vcf {input.vcf} \
            --vcf-ref {input.vcf_ref} \
            --out {params.out} 2> {log}
        """

rule norm_xpnsl_permuted:
    input:
        lambda w: expand(rules.selscan_xpnsl_permuted.output, chrom=CHROMOSOMES, hab_comb=w.hab_comb, n=w.n)
    output:
        expand(f'{SWEEPS_DIR}/xpnsl/permuted/{{chrom}}/{{chrom}}_{{hab_comb}}_{{n}}_permuted.xpnsl.out.norm', chrom=CHROMOSOMES, allow_missing=True)
    log: f"{LOG_DIR}/norm_xpnsl_permuted/{{hab_comb}}_{{n}}.log"
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    shell:
        """
        norm --xpnsl --qbins 10 --files {input} 2> {log} 
        """

rule write_windowed_statistics_permuted:
    input:
        xpnsl = rules.norm_xpnsl_permuted.output
    output:
        xpnsl_df = f'{SWEEPS_DIR}/permuted/windowed_stats/{{hab_comb}}_{{n}}_windowed_xpnsl_permuted.txt'
    params:
        winsize = 50000,
        nSites_xpnsl = 40 
    conda: '../envs/sweeps.yaml'
    script:
        "../scripts/r/write_windowed_statistics_permuted.R"

###############
#### iHH12 ####
###############

rule ihh_OneTwo:
    input:
        vcf = rules.bcftools_splitVCF_byHabitat.output.vcf,
        genMap = rules.genMap_toPlinkFormat.output
    output:
        '{0}/ihh12/{{chrom}}/{{chrom}}_{{habitat}}.ihh12.out'.format(SWEEPS_DIR)
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    log: LOG_DIR + '/selscan_ihh12/{chrom}_{habitat}_ihh12.log'
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    threads: 2
    params:
        out = '{0}/ihh12/{{chrom}}/{{chrom}}_{{habitat}}'.format(SWEEPS_DIR) 
    shell:
        """
        selscan --ihh12 \
            --vcf {input.vcf} \
            --map {input.genMap} \
            --threads {threads} \
            --out {params.out} 2> {log}
        """

rule norm_ihh_OneTwo:
    input:
        lambda w: expand(rules.ihh_OneTwo.output, chrom=CHROMOSOMES, habitat=w.habitat)
    output:
        expand(f'{SWEEPS_DIR}/ihh12/{{chrom}}/{{chrom}}_{{habitat}}.ihh12.out.norm', chrom=CHROMOSOMES, allow_missing=True)
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    shell:
        """
        norm --ihh12 --qbins 10 --files {input} 
        """
        
##################
#### ANALYSES ####
##################

rule write_windowed_statistics:
    input:
        fst = expand(rules.windowed_fst.output, chrom=CHROMOSOMES, hab_comb='Urban_Rural'),
        thetaU = expand(rules.windowed_theta.output, chrom=CHROMOSOMES, habitat='Urban'),
        thetaR = expand(rules.windowed_theta.output, chrom=CHROMOSOMES, habitat='Rural'),
        xpnsl = expand(rules.norm_xpnsl.output, chrom=CHROMOSOMES, hab_comb='Urban_Rural')
    output:
        sfs_df = f'{SWEEPS_DIR}/analyses/windowed_fst_thetas.txt',
        xpnsl_df = f'{SWEEPS_DIR}/analyses/windowed_xpnsl.txt'
    params:
        winsize = 50000,
        nSites_fst = 1500,
        nSites_xpnsl = 40 
    conda: '../envs/sweeps.yaml'
    script:
        "../scripts/r/write_windowed_statistics.R"

rule create_geneToGO_mapfile:
    input:
        GFF_FILE
    output:
        f'{SWEEPS_DIR}/analyses/go/gene2go.map'
    run:
        import re
        with open(output[0], 'w') as fout:
            with open(input[0], 'r') as fin:
                lines = fin.readlines()
                for l in lines:
                    if not l.startswith('#'):
                        sline = l.strip().split('\t')
                        feat = sline[2]
                        if feat == 'mRNA':
                            atts = sline[8]
                            id = re.search('(?<=ID\\=)ACLI19_g\\d+\\.t\\d+(?=;)', atts)[0]
                            transcript = id.split('.')[1]
                            if transcript == 't1':
                                gene = id.split('.')[0]
                                go = re.findall('(GO:\\d+)', atts)
                                go_string = ', '.join(go)
                                fout.write(f'{gene}\t{go_string}\n')


##############
#### POST ####
##############

rule sweeps_done:
    input:
        expand(rules.norm_xpnsl.output, hab_comb=['Urban_Rural','Rural_Suburban']),
        expand(rules.norm_ihh_OneTwo.output, habitat=HABITATS),
        expand(rules.windowed_fst.output, chrom=CHROMOSOMES, hab_comb=HABITAT_COMBOS),
        expand(rules.angsd_fst_allSites_readable.output, chrom=CHROMOSOMES, hab_comb=HABITAT_COMBOS),
        expand(rules.windowed_theta.output, chrom=CHROMOSOMES, habitat=HABITATS),
        expand(rules.angsd_thetas_allSites_readable.output, chrom=CHROMOSOMES, habitat=HABITATS)
    output:
        '{0}/sweeps.done'.format(SWEEPS_DIR)
    shell:
        """
        touch {output}
        """
