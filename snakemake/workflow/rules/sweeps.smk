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
    """
    Create text file withBAM files for samples in each habitat
    """
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
    """
    Convert genetic map to Plink format
    """
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
    """
    Estimate thetas in windows across the genome
    """
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
    """
    Estimate Hudson's Fst in windows across the genome. 
    """
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

#############################
#### PIXY FST AND THETAS ####
#############################

rule create_pixy_popfile:
    """
    Create poulation file for Pixy
    """
    input:
       config['samples']
    output:
        f"{PROGRAM_RESOURCE_DIR}/pixy/popfile.txt"
    run:
        with open(input[0], "r") as fin:
            with open(output[0], "w") as fout:
                lines = fin.readlines()
                for l in lines:
                    sl = l.split('\t')
                    sample = sl[0]
                    pop = sl[1]
                    if sample in FINAL_SAMPLES:
                        fout.write(f"{sample}\t{pop}\n")

rule pixy:
    """
    Run Pixy to estimate windowed stats from VCFs
    """
    input:
        vcf = rules.concat_variant_invariant_sites.output.vcf,
        tbi = rules.concat_variant_invariant_sites.output.tbi,
        popu = rules.create_pixy_popfile.output
    output:
        fst = f"{PIXY_DIR}/{{chrom}}/{{chrom}}_miss{{miss}}_win{{win_size}}_pixy_fst.txt",
        pi = f"{PIXY_DIR}/{{chrom}}/{{chrom}}_miss{{miss}}_win{{win_size}}_pixy_pi.txt",
        dxy = f"{PIXY_DIR}/{{chrom}}/{{chrom}}_miss{{miss}}_win{{win_size}}_pixy_dxy.txt"
    log: f"{LOG_DIR}/pixy/{{chrom}}_miss{{miss}}_win{{win_size}}_pixy.log"
    conda: "../envs/sweeps.yaml"
    params:
        out = f"{PIXY_DIR}/{{chrom}}",
        pref = f"{{chrom}}_miss{{miss}}_win{{win_size}}_pixy"
    shell:
        """
        pixy --stats fst pi dxy \
            --population {input.popu} \
            --vcf {input.vcf} \
            --window_size {wildcards.win_size} \
            --output_folder {params.out} \
            --output_prefix {params.pref} \
            --fst_type hudson &> {log}
        """

rule pixy_perSite:
    """
    Run Pixy to estimate per-site Fst from VCFs
    """
    input:
        vcf = rules.shapeit_phase.output.vcf,
        tbi = rules.shapeit_phase.output.idx,
        popu = rules.create_pixy_popfile.output
    output:
        fst = f"{PIXY_DIR}/{{chrom}}/{{chrom}}_perSite_pixy_fst.txt",
    log: f"{LOG_DIR}/pixy/{{chrom}}_perSite_pixy.log"
    conda: "../envs/sweeps.yaml"
    params:
        out = f"{PIXY_DIR}/{{chrom}}",
        pref = f"{{chrom}}_perSite_pixy"
    threads: 4
    shell:
        """
        pixy --stats fst \
            --population {input.popu} \
            --vcf {input.vcf} \
            --window_size 1 \
            --bypass_invariant_check yes \
            --output_folder {params.out} \
            --output_prefix {params.pref} \
            --n_cores {threads} \
            --fst_type hudson &> {log}
        """

################
#### XP-NSL ####
################

rule selscan_xpnsl:
    """
    Estimate XP-nSL
    """
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
    """
    Normalize XP-nSL
    """
    input:
        lambda w: expand(rules.selscan_xpnsl.output, chrom=CHROMOSOMES, hab_comb=w.hab_comb)
    output:
        expand(f'{SWEEPS_DIR}/xpnsl/{{chrom}}/{{chrom}}_{{hab_comb}}.xpnsl.out.norm', chrom=CHROMOSOMES, allow_missing=True)
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    shell:
        """
        norm --xpnsl --qbins 10 --files {input} 
        """

rule remove_outlier_pops_fromVCF:
    """
    Remove population structure outlier population from VCFs for re-running XP-nSL.
    Urban: Samples from population 40
    Rural: Samples from population 7
    """
    input:
        vcf = lambda w: expand(rules.bcftools_splitVCF_byHabitat.output.vcf, chrom=w.chrom, habitat=w.habitat)
    output:
        vcf = f"{FREEBAYES_DIR}/vcf/{{chrom}}/{{chrom}}_{{habitat}}_outlierPopsRemoved.vcf.gz",
        idx = f"{FREEBAYES_DIR}/vcf/{{chrom}}/{{chrom}}_{{habitat}}_outlierPopsRemoved.vcf.gz.tbi"
    log: f"{LOG_DIR}/remove_outlier_pops_fromVCF/{{chrom}}_{{habitat}}"
    conda: "../envs/sweeps.yaml"
    params:
        to_rem = lambda w: "^s_40_1,s_40_3,s_40_6,s_40_7,s_40_8,s_40_10,s_40_12,s_40_17,s_40_19" if w.habitat == "Urban"
                    else "^s_7_4,s_7_6,s_7_7,s_7_11,s_7_13,s_7_16,s_7_19,s_7_20"
    shell:
        """
        ( bcftools view -O z -o {output.vcf} \
            --samples {params.to_rem} {input.vcf} &&
        sleep 5
        tabix {output.vcf} ) 2> {log}
        """
    
rule selscan_xpnsl_outlierRem:
    """
    Estimate XP-nSL with outlier populations removed
    """
    input:
        vcf = lambda w: expand(rules.remove_outlier_pops_fromVCF.output.vcf, chrom=w.chrom, habitat="Urban"),
        vcf_ref = lambda w: expand(rules.remove_outlier_pops_fromVCF.output.vcf, chrom=w.chrom, habitat="Rural")
    output:
        temp('{0}/xpnsl/{{chrom}}/{{chrom}}_Urban_Rural_outlierRem.xpnsl.out'.format(SWEEPS_DIR))
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    log: LOG_DIR + '/selscan_xpnsl/{chrom}_Urban_Rural_xpnsl_outlierRem.log'
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    threads: 2
    params:
        out = '{0}/xpnsl/{{chrom}}/{{chrom}}_Urban_Rural_outlierRem'.format(SWEEPS_DIR) 
    shell:
        """
        selscan --xpnsl \
            --vcf {input.vcf} \
            --vcf-ref {input.vcf_ref} \
            --threads {threads} \
            --out {params.out} 2> {log}
        """

rule norm_xpnsl_outlierRem:
    """
    Normalize XP-nSL with outlier populations removed
    """
    input:
        lambda w: expand(rules.selscan_xpnsl_outlierRem.output, chrom=CHROMOSOMES)
    output:
        expand(f'{SWEEPS_DIR}/xpnsl/{{chrom}}/{{chrom}}_Urban_Rural_outlierRem.xpnsl.out.norm', chrom=CHROMOSOMES)
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    shell:
        """
        norm --xpnsl --qbins 10 --files {input} 
        """

#############################
#### XP-nSL PERMUTATIONS ####
#############################

rule permuted_samples_byHabitat:
    """
    Generate text files with permuted BAM files
    """
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
    """
    Split Permuted VCFs by habitat
    """
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
    """
    Run selscan to estimate XP-nSL on permuted samples
    """
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
    """
    Normalize permuted XP-nSL scores
    """
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

#############
#### nSL ####
#############

rule nsl:
    """
    Estimate nSL in each habitat
    """
    input:
        vcf = rules.bcftools_splitVCF_byHabitat.output.vcf,
        genMap = rules.genMap_toPlinkFormat.output
    output:
        '{0}/nsl/{{chrom}}/{{chrom}}_{{habitat}}.nsl.out'.format(SWEEPS_DIR)
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    log: LOG_DIR + '/selscan_nsl/{chrom}_{habitat}_nsl.log'
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    threads: 2
    params:
        out = '{0}/nsl/{{chrom}}/{{chrom}}_{{habitat}}'.format(SWEEPS_DIR) 
    shell:
        """
        selscan --nsl \
            --vcf {input.vcf} \
            --map {input.genMap} \
            --threads {threads} \
            --out {params.out} 2> {log}
        """

rule norm_nsl:
    """
    Normalize nSL
    """
    input:
        lambda w: expand(rules.nsl.output, chrom=CHROMOSOMES, habitat=w.habitat)
    output:
        expand(f'{SWEEPS_DIR}/nsl/{{chrom}}/{{chrom}}_{{habitat}}.nsl.out.100bins.norm', chrom=CHROMOSOMES, allow_missing=True)
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    shell:
        """
        norm --nsl --qbins 10 --files {input} 
        """

#############
#### iHs ####
#############

rule ihs:
    """
    Estimate iHS
    """
    input:
        vcf = rules.bcftools_splitVCF_byHabitat.output.vcf,
        genMap = rules.genMap_toPlinkFormat.output
    output:
        '{0}/ihs/{{chrom}}/{{chrom}}_{{habitat}}.ihs.out'.format(SWEEPS_DIR)
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    log: LOG_DIR + '/selscan_ihs/{chrom}_{habitat}_ihs.log'
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    threads: 2
    params:
        out = '{0}/ihs/{{chrom}}/{{chrom}}_{{habitat}}'.format(SWEEPS_DIR) 
    shell:
        """
        selscan --ihs \
            --vcf {input.vcf} \
            --map {input.genMap} \
            --threads {threads} \
            --out {params.out} 2> {log}
        """

rule norm_ihs:
    """
    Normalize iHS
    """
    input:
        lambda w: expand(rules.ihs.output, chrom=CHROMOSOMES, habitat=w.habitat)
    output:
        expand(f'{SWEEPS_DIR}/ihs/{{chrom}}/{{chrom}}_{{habitat}}.ihs.out.100bins.norm', chrom=CHROMOSOMES, allow_missing=True)
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    shell:
        """
        norm --ihs --qbins 10 --files {input} 
        """

###############
#### iHH12 ####
###############

rule ihh_OneTwo:
    """
    Estimate IHH12
    """
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
    """
    Normalize iHH12
    """
    input:
        lambda w: expand(rules.ihh_OneTwo.output, chrom=CHROMOSOMES, habitat=w.habitat)
    output:
        expand(f'{SWEEPS_DIR}/ihh12/{{chrom}}/{{chrom}}_{{habitat}}.ihh12.out.norm', chrom=CHROMOSOMES, allow_missing=True)
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    shell:
        """
        norm --ihh12 --qbins 10 --files {input} 
        """
        
#####################
#### LASSI-PLUS ####
####################

rule generate_salti_spectra:
    """
    Generate haplotype frequency spectra for saltiLassi
    """
    input:
        vcf = rules.bcftools_splitVCF_byHabitat.output.vcf,
        popf = rules.create_pixy_popfile.output
    output:
        spec = '{0}/lassip/{{chrom}}/{{chrom}}_salti.{{habitat}}.lassip.hap.spectra.gz'.format(SWEEPS_DIR)
    container: 'library://james-s-santangelo/lassip/lassip:1.1.1'
    log: LOG_DIR + '/lassip/{chrom}_{habitat}_lassip.log'
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    params:
        out = '{0}/lassip/{{chrom}}/{{chrom}}_salti'.format(SWEEPS_DIR) 
    shell:
        """
        grep '{wildcards.habitat}' {input.popf} > {wildcards.chrom}{wildcards.habitat}_popfile.tmp
        lassip --vcf {input.vcf} \
            --hapstats \
            --calc-spec \
            --winsize 201 \
            --winstep 50 \
            --k 20 \
            --out {params.out} \
            --salti \
            --pop {wildcards.chrom}{wildcards.habitat}_popfile.tmp 2> {log}
        rm {wildcards.chrom}{wildcards.habitat}_popfile.tmp 
        """

rule analyze_salti_spectra:
    """
    Estimate saltiLassi lamda, the number of sweeping haplotypes, and the width of selective sweeps
    """
    input:
        spectra = lambda w: expand(rules.generate_salti_spectra.output.spec, chrom=CHROMOSOMES, habitat=w.habitat),
        maps = expand(rules.genMap_toPlinkFormat.output, chrom=CHROMOSOMES)
    output:
        f'{SWEEPS_DIR}/lassip/{{habitat}}_salti.lassip.hap.out.gz'
    log: LOG_DIR + '/lassip/{habitat}_analyze.log'
    container: 'library://james-s-santangelo/lassip/lassip:1.1.1'
    threads: 16
    params:
        out = '{0}/lassip/{{habitat}}_salti'.format(SWEEPS_DIR) 
    shell:
        """
        cat {input.maps} > {wildcards.habitat}_allChroms.map
        lassip --spectra {input.spectra} \
            --threads {threads} \
            --out {params.out} \
            --salti \
            --k 20 \
            --map {wildcards.habitat}_allChroms.map \
            --dist-type cm &> {log}
        rm {wildcards.habitat}_allChroms.map
        """

#####################
#### PAIRWISE LD ####
#####################

rule plink_pairwise_ld:
    """
    Estimate pairwise LD (as R2) between variants within 50 Kb distance
    """
    input:
        vcf = rules.bcftools_splitVCF_byHabitat.output.vcf
    output:
        f"{PLINK_DIR}/pairwise_ld/{{chrom}}_{{habitat}}.ld.gz"
    log: f"{LOG_DIR}/plink/{{chrom}}_{{habitat}}_plink_ld.log"
    conda: "../envs/sweeps.yaml"
    params:
        out = f"{PLINK_DIR}/pairwise_ld/{{chrom}}_{{habitat}}"
    shell:
        """
        plink --vcf {input.vcf} \
            --allow-extra-chr \
            --set-missing-var-ids @:# \
            --double-id \
            -r2 gz \
            --ld-window 100000\
            --ld-window-kb 50 \
            --ld-window-r2 0 \
            -out {params.out} 2> {log}
        """

##################
#### ANALYSES ####
##################

rule write_windowed_sfs_stats:
    """
    Create file with all windowed sfs stats
    """
    input:
        thetaU = expand(rules.windowed_theta.output, chrom=CHROMOSOMES, habitat='Urban'),
        thetaR = expand(rules.windowed_theta.output, chrom=CHROMOSOMES, habitat='Rural'),
        fst = expand(rules.windowed_fst.output, chrom=CHROMOSOMES, hab_comb='Urban_Rural'),
    output:
        sfs_df = f'{SWEEPS_DIR}/analyses/windowed_sfs_stats.txt',
    conda: '../envs/r.yaml'
    script:
        "../scripts/r/write_windowed_sfs_stats.R"

rule write_windowed_singPop_hapstats:
    """
    Create file with windowed single population haplotype stats
    """
    input:
        norm = get_windowed_singPop_hapstats_input_files
    output:
        hapstats_df = f'{SWEEPS_DIR}/analyses/windowed_{{stat}}.txt'
    params:
        winsize = 50000
    conda: '../envs/r.yaml'
    script:
        "../scripts/r/write_windowed_singPop_hapstats.R"

rule write_windowed_xpnsl:
    """
    Create file with windowed XP-nSL stats
    """
    input:
        norm = get_windowed_xpnsl_input_files
    output:
        hapstats_df = f'{SWEEPS_DIR}/analyses/windowed_{{hab_comb}}_xpnsl.txt'
    params:
        winsize = 50000
    conda: '../envs/r.yaml'
    script:
        "../scripts/r/write_windowed_xpnsl.R"

rule write_windowed_ld:
    """
    Create file with windowed pairwise LD
    """
    input:
        # ld = "/scratch/projects/trifolium/gwsd/results/plink/pairwise_ld/Chr04_Pall_Rural_test.ld.gz"  
        ld = rules.plink_pairwise_ld.output 
    output:
        win_ld = f'{PLINK_DIR}/windowed/{{chrom}}_{{habitat}}_windowed_ld.txt'
        # win_ld = f'{PLINK_DIR}/windowed/Chr04_Pall_Rural_test.txt'
    params:
        winsize = 50000
    conda: '../envs/r.yaml'
    script:
        "../scripts/r/write_windowed_ld.R"

rule write_windowed_xpnsl_outlierRem:
    """
    Create file with windowed XP-nSL stats with outlier populations remove
    """
    input:
        norm = expand(rules.norm_xpnsl_outlierRem.output, chrom=CHROMOSOMES)
    output:
        hapstats_df = f'{SWEEPS_DIR}/analyses/windowed_Urban_Rural_xpnsl_outlierRem.txt'
    params:
        winsize = 50000
    conda: '../envs/r.yaml'
    script:
        "../scripts/r/write_windowed_xpnsl_outlierRem.R"

rule write_windowed_xpnsl_permuted:
    """
    Create file with windowed XP-nSL stats for permuted samples
    """
    input:
        xpnsl = rules.norm_xpnsl_permuted.output
    output:
        xpnsl_df = f'{SWEEPS_DIR}/analyses/permuted_xpnsl/{{hab_comb}}_{{n}}_windowed_xpnsl_permuted.txt'
    params:
        winsize = 50000,
    conda: '../envs/r.yaml'
    script:
        "../scripts/r/write_windowed_xpnsl_permuted.R"

rule create_geneToGO_mapfile:
    """
    Creating mapping file for GO analysis
    """
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

rule compare_observed_permuted_xpnsl:
    """
    Compare observed to permuted XP-nSL distributions
    """
    input:
        obs = expand(rules.write_windowed_singPop_hapstats.output, stat="xpnsl"),
        perm = expand(rules.write_windowed_xpnsl_permuted.output, hab_comb="Urban_Rural", n=[x for x in range(1,1001)])
    output:
        urb_perc = f"{FIGURES_DIR}/selection/xpnsl_perm/urban_percentiles.txt",
        rur_perc = f"{FIGURES_DIR}/selection/xpnsl_perm/rural_percentiles.txt",
    conda: '../envs/r.yaml'
    notebook:
        "../notebooks/compare_observed_permuted_xpnsl.r.ipynb"

rule install_genotype_plot:
    """
    Install genotype_plot into conda env
    """
    output:
        f"{PROGRAM_RESOURCE_DIR}/genotype_plot_install.done"
    conda: "../envs/r.yaml"
    shell:
        """
        R -e 'remotes::install_github("JimWhiting91/genotype_plot")' &&
        R -e 'library(GenotypePlot)' &&
        touch {output}
        """

rule outlier_analysis:
    """
    Perform outlier analysis and generate figures (e.g., Manhattan plot)
    """
    input:
        gt_plot = rules.install_genotype_plot.output,
        chr_lengths = rules.genome_lengths_file.output, 
        win_sfs_fst = rules.write_windowed_sfs_stats.output.sfs_df,
        win_xpnsl_ur = expand(rules.write_windowed_xpnsl.output, hab_comb=['Urban_Rural']),
        win_xpnsl_ur_outRem = expand(rules.write_windowed_xpnsl_outlierRem.output),
        win_xpnsl_sr = expand(rules.write_windowed_xpnsl.output, hab_comb=['Suburban_Rural']),
        win_xpnsl_us = expand(rules.write_windowed_xpnsl.output, hab_comb=['Urban_Suburban']),
        win_xpnsl_perm = expand(rules.write_windowed_xpnsl_permuted.output, hab_comb="Urban_Rural", n=[x for x in range(1,1001)]),
        win_nsl = expand(rules.write_windowed_singPop_hapstats.output, stat="nsl"),
        win_ihh12 = expand(rules.write_windowed_singPop_hapstats.output, stat="ihh12"),
        win_ihs = expand(rules.write_windowed_singPop_hapstats.output, stat="ihs"),
        norm_xpnsl = expand(rules.norm_xpnsl.output, hab_comb=['Urban_Rural']),
        norm_ihh12 = expand(rules.norm_ihh_OneTwo.output, habitat=['Urban', 'Rural']),
        norm_ihs = expand(rules.norm_ihs.output, habitat=['Urban', 'Rural']),
        norm_nsl = expand(rules.norm_nsl.output, habitat=['Urban', 'Rural']),
        gt_win_fst = expand(rules.pixy.output.fst, win_size="50000", miss="0", chrom=CHROMOSOMES),
        gt_win_pi = expand(rules.pixy.output.pi, win_size="50000", miss="0", chrom=CHROMOSOMES),
        win_ld = expand(rules.write_windowed_ld.output, habitat=["Urban", "Rural"], chrom=CHROMOSOMES),
        gen_map = rules.interpolate_genetic_map.output.genMap_interp,
        lassip = expand(rules.analyze_salti_spectra.output, habitat=["Urban", "Rural"]),
        spec = expand(rules.generate_salti_spectra.output, chrom=CHROMOSOMES, habitat=["Urban", "Rural"]),
        vcfs = expand(rules.shapeit_phase.output.vcf, chrom=CHROMOSOMES),
        popmap = rules.create_pixy_popfile.output,
        gff = GFF_FILE 
    output:
        xpnsl_nSites_hist = f'{FIGURES_DIR}/selection/xpnsl_nSites_histogram.pdf',
        xpnsl_manhat_ur = f"{FIGURES_DIR}/selection/manhattan/urban_rural_xpnsl_windowed_manhat.pdf",
        xpnsl_manhat_ur_outRem = f"{FIGURES_DIR}/selection/manhattan/urban_rural_xpnsl_windowed_manhat_outRem.pdf",
        xpnsl_manhat_sr = f"{FIGURES_DIR}/selection/manhattan/suburban_rural_xpnsl_windowed_manhat.pdf",
        xpnsl_manhat_us = f"{FIGURES_DIR}/selection/manhattan/urban_suburban_xpnsl_windowed_manhat.pdf",
        cor_plot = f"{FIGURES_DIR}/selection/xpnsl_perm/observed_permuted_xpnsl_correlation.pdf",
        urb_mean_plot = f"{FIGURES_DIR}/selection/xpnsl_perm/urbanSel_mean.pdf",
        urb_prop_plot = f"{FIGURES_DIR}/selection/xpnsl_perm/urbanSel_prop.pdf",
        rur_mean_plot = f"{FIGURES_DIR}/selection/xpnsl_perm/ruralSel_mean.pdf",
        rur_prop_plot = f"{FIGURES_DIR}/selection/xpnsl_perm/ruralSel_prop.pdf",
        xpnsl_df = f"{FIGURES_DIR}/tables/xpnsl_outliers.txt",
        xpnsl_out_genes = f"{FIGURES_DIR}/tables/xpnsl_outlier_genes.txt",
        rur_nsl_manhat = f"{FIGURES_DIR}/selection/manhattan/rural_nSL_windowed_manhat.pdf",
        urb_nsl_manhat = f"{FIGURES_DIR}/selection/manhattan/urban_nSL_windowed_manhat.pdf",
        nsl_df = f"{FIGURES_DIR}/tables/nsl_outliers.txt",
        rur_ihs_manhat = f"{FIGURES_DIR}/selection/manhattan/rural_iHS_windowed_manhat.pdf",
        urb_ihs_manhat = f"{FIGURES_DIR}/selection/manhattan/urban_iHS_windowed_manhat.pdf",
        ihs_df = f"{FIGURES_DIR}/tables/ihs_outliers.txt",
        rur_ihh12_manhat = f"{FIGURES_DIR}/selection/manhattan/rural_iHH12_windowed_manhat.pdf",
        urb_ihh12_manhat = f"{FIGURES_DIR}/selection/manhattan/urban_iHH12_windowed_manhat.pdf",
        ihh12_df = f"{FIGURES_DIR}/tables/ihh12_outliers.txt",
        gt_fst_manhat = f"{FIGURES_DIR}/selection/manhattan/urban_rural_gt_fst_windowed_manhat.pdf",
        tajima_manhat = f"{FIGURES_DIR}/selection/manhattan/urban_rural_delta_tajima_manhat.pdf",
        gt_fst_df = f"{FIGURES_DIR}/tables/gt_fst_outliers.txt",
        ld_manhat = f"{FIGURES_DIR}/selection/manhattan/urban_rural_ld_manhat.pdf",
        recomb_manhat = f"{FIGURES_DIR}/selection/manhattan/recombination_rate_manhat.pdf",
        xpnsl_vs_recomb_plot = f"{FIGURES_DIR}/selection/xpnsl_vs_recomb_plot.pdf",
        fst_vs_recomb_plot = f"{FIGURES_DIR}/selection/fst_vs_recomb_plot.pdf",
        num_xpnsl_by_num_genes = f"{FIGURES_DIR}/selection/num_xpnsl_by_num_genes.pdf",
        top_hits_genes = f'{FIGURES_DIR}/tables/topHits_selected_regions_genes.txt', 
        top_hits_tbl = f'{FIGURES_DIR}/tables/topHits_selected_regions_urban_rural_table.txt',
        salti_df = f"{FIGURES_DIR}/tables/salti_outliers.txt",
        urban_salti_manhat = f"{FIGURES_DIR}/selection/manhattan/urban_salti_manhat.pdf",
        rural_salti_manhat = f"{FIGURES_DIR}/selection/manhattan/rural_salti_manhat.pdf",
        salti_m_hist = f"{FIGURES_DIR}/selection/salti_m_histogram.pdf",
        salti_sr_hist = f"{FIGURES_DIR}/selection/salti_sr_histogram.pdf",
        logA_perm_plot = f"{FIGURES_DIR}/selection/logA_perm_plot.pdf",
        logA_ur_density = f"{FIGURES_DIR}/selection/logA_ur_density.pdf",
        Chr04_Occ_urb_xpnsl = f"{FIGURES_DIR}/selection/region_plots/Chr04_Occ_urb_xpnsl.pdf",
        Chr04_Occ_urb_ur_haps = f"{FIGURES_DIR}/selection/region_plots/Chr04_Occ_urb_ur_haps.pdf",
        Chr04_Occ_urb_ur_af = f"{FIGURES_DIR}/selection/region_plots/Chr04_Occ_urb_ur_af.pdf",
        Chr04_Occ_urb_ur_pca = f"{FIGURES_DIR}/selection/region_plots/Chr04_Occ_urb_ur_pca.pdf",
        Chr05_Occ_urb_xpnsl = f"{FIGURES_DIR}/selection/region_plots/Chr05_Occ_urb_xpnsl.pdf",
        Chr05_Occ_urb_ur_haps = f"{FIGURES_DIR}/selection/region_plots/Chr05_Occ_urb_ur_haps.pdf",
        Chr05_Occ_urb_ur_af = f"{FIGURES_DIR}/selection/region_plots/Chr05_Occ_urb_ur_af.pdf",
        Chr05_Occ_urb_ur_pca = f"{FIGURES_DIR}/selection/region_plots/Chr05_Occ_urb_ur_pca.pdf",
        Chr04_Occ_rur_xpnsl = f"{FIGURES_DIR}/selection/region_plots/Chr04_Occ_rur_xpnsl.pdf",
        Chr04_Occ_rur_ur_haps = f"{FIGURES_DIR}/selection/region_plots/Chr04_Occ_rur_ur_haps.pdf",
        Chr04_Occ_rur_ur_af = f"{FIGURES_DIR}/selection/region_plots/Chr04_Occ_rur_ur_af.pdf",
        Chr04_Occ_rur_ur_pca = f"{FIGURES_DIR}/selection/region_plots/Chr04_Occ_rur_ur_pca.pdf",
        Chr08_Pall_rur_xpnsl = f"{FIGURES_DIR}/selection/region_plots/Chr08_Pall_rur_xpnsl.pdf",
        Chr08_Pall_rur_ur_haps = f"{FIGURES_DIR}/selection/region_plots/Chr08_Pall_rur_ur_haps.pdf",
        Chr08_Pall_rur_ur_af = f"{FIGURES_DIR}/selection/region_plots/Chr08_Pall_rur_ur_af.pdf",
        Chr08_Pall_rur_ur_pca = f"{FIGURES_DIR}/selection/region_plots/Chr08_Pall_rur_ur_pca.pdf",
        random_unsel_regions_af = f"{FIGURES_DIR}/selection/region_plots/random_unsel_regions_af.pdf",
        random_unsel_regions_pca = f"{FIGURES_DIR}/selection/region_plots/random_unsel_regions_pca.pdf"
    conda: '../envs/r.yaml'
    notebook:
        "../notebooks/outlier_analysis.r.ipynb"

rule generate_dosage_matrix:
    """
    Convert VCFs to dosage matrices for easier calculation of allele frequencies
    """
    input:
        vcf = rules.shapeit_phase.output.vcf,
        tbi = rules.shapeit_phase.output.idx
    output:
        f"{PROGRAM_RESOURCE_DIR}/dosage_matrices/{{chrom}}_dosages.txt"
    container: "docker://ghcr.io/thewanglab/algatr"
    script:
        "../scripts/r/generate_dosage_matrix.R"

rule cline_analysis:
    """
    Examine whether SNPs inselected regions show clina variation
    """
    input:
        top_hits = rules.outlier_analysis.output.top_hits_tbl,
        fst = expand(rules.pixy_perSite.output, chrom=CHROMOSOMES),
        dos = expand(rules.generate_dosage_matrix.output, chrom=CHROMOSOMES)
    output:
        max_fst_df = f"{FIGURES_DIR}/tables/selected_regions_max_fst_sites.txt",
        cline_plot = f"{FIGURES_DIR}/selection/max_fst_sites_clines.pdf",
        selSites_glm_df = f"{FIGURES_DIR}/tables/selected_regions_max_fst_glm.txt",
        freq_byHab_plot = f"{FIGURES_DIR}/selection/selected_regions_max_fst_byHabitat.pdf",
        freq_diff_plot = f"{FIGURES_DIR}/selection/selected_regions_max_fst_frequency_difference.pdf",
    conda: "../envs/r.yaml"
    notebook:
        "../notebooks/cline_analysis.r.ipynb"
        

rule go_enrichment_analysis:
    """
    Perform GO enrichment analysis
    """
    input:
        all_genes = rules.create_geneToGO_mapfile.output,
        top_ten_genes = rules.outlier_analysis.output.top_hits_genes
    output:
        all_go_res = f'{FIGURES_DIR}/tables/all_go_results.txt'
    conda: '../envs/r.yaml'
    notebook:
        "../notebooks/go_enrichment_analysis.r.ipynb"


##############
#### POST ####
##############

rule sweeps_done:
    """
    Create empty file signaling completion of selective sweeps analysis
    """
    input:
        rules.go_enrichment_analysis.output,
        rules.outlier_analysis.output,
        rules.cline_analysis.output,
    output:
        '{0}/sweeps.done'.format(SWEEPS_DIR)
    shell:
        """
        touch {output}
        """
