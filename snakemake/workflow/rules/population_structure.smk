# Rules to infer population structure and admixture proportions

###############
#### SETUP ####
###############

# Estimate LD among 4fold SNPs with MAF > 0.05. 
# Prune these SNPs within 20 Kb using r-squared cutoff of 0.2

rule create_pos_file_for_ngsLD:
    """
    Create text file with site positions for NGSLD
    """
    input:
        rules.angsd_gl_degenerate_allSamples.output.mafs
    output:
        '{0}/ngsld_pos/{{chrom}}_{{site}}_maf{{maf}}.pos'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/create_pos_file_for_ngsLD/{chrom}_{site}_maf{maf}_pos.log'
    shell:
        """
        zcat {input} | cut -f 1,2 | tail -n +2 > {output} 2> {log}
        """

rule ngsLD_degenerateSites:
    """
    Estimate pairwise LD among 4fold sites
    """
    input:
        pos = rules.create_pos_file_for_ngsLD.output,
        gls = rules.angsd_gl_degenerate_allSamples.output.gls
    output:
        '{0}/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}.ld.gz'.format(NGSLD_DIR)
    log: LOG_DIR + '/ngsld/{chrom}_{site}_maf{maf}_calc_ld.log'
    container: 'library://james-s-santangelo/ngsld/ngsld:1.1.1'
    threads: 8
    params:
        n_ind = len(FINAL_SAMPLES)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '1:00:00'
    shell:
        """
        ( NUM_SITES=$(cat {input.pos} | wc -l) &&
          ngsLD --geno {input.gls} \
            --pos {input.pos} \
            --n_ind {params.n_ind} \
            --n_sites $NUM_SITES \
            --probs \
            --n_threads {threads} \
            --max_kb_dist 20 | gzip --best > {output} ) 2> {log}
        """

rule prune_degenerateSNPs_forPopStructure:
    """
    Create text file with IDs for sites to keep post LD pruning
    """
    input:
        rules.ngsLD_degenerateSites.output
    output:
        '{0}/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_pruned.id'.format(NGSLD_DIR)
    log: LOG_DIR + '/prune_degenerateSNP_forPopStructure/{chrom}_{site}_maf{maf}_prune_ld.log'
    container: 'library://james-s-santangelo/ngsld/ngsld:1.1.1'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '3:00:00'
    shell:
        """
        ( zcat {input} | perl /opt/bin/prune_graph.pl \
            --max_kb_dist 20 \
            --min_weight 0.2 | sort -V > {output} ) 2> {log}
        """

rule pruneGLs_degenerateSNPs:
    """
    Prune sites for LD
    """
    input:
        gls = rules.angsd_gl_degenerate_allSamples.output.gls,
        pos = rules.prune_degenerateSNPs_forPopStructure.output
    output:
        '{0}/gls/allSamples/{{site}}/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_pruned.beagle.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/pruneGLs_degenerateSNPs/{chrom}_{site}_maf{maf}_pruneGLs.log'
    params:
        out = '{0}/gls/allSamples/{{site}}/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_pruned.beagle'.format(ANGSD_DIR)
    shell:
        """
        ( zgrep 'marker' {input.gls} > {params.out} &&
                sed 's/:/_/g' {input.pos} | zgrep -w -f - {input.gls} >> {params.out} &&
                gzip {params.out} ) 2> {log}
        """

rule concat_angsd_gl:
    """
    Concatenated GLs from all 16 chromosomes into single file. Done separately for each site type.
    """
    input:
    	lambda wildcards: expand(rules.pruneGLs_degenerateSNPs.output, chrom=CHROMOSOMES, site=wildcards.site, maf=wildcards.maf)
    output:
        '{0}/gls/allSamples/{{site}}/allChroms_{{site}}_maf{{maf}}_pruned.beagle.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/concat_angsd_gl/allSamples_{site}_{maf}_concat.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933' 
    shell:
        """
        first=1
        for f in {input}; do
            if [ "$first"  ]; then
                zcat "$f"
                first=
            else
                zcat "$f"| tail -n +2
            fi
        done | bgzip -c > {output} 2> {log}
        """

#########################
#### PCA & ADMIXTURE ####
#########################

# PCA & admixture analysis using LD-pruned 4fold SNPs from above (MAF > 0.05)

rule pcangsd:
    """
    Estimate variance-covariance matrix of allele frequencies among samples from 4fold sites
    """
    input:
        rules.concat_angsd_gl.output
    output:
        '{0}/pcangsd/allSamples_allChroms_{{site}}_maf{{maf}}_pcangsd.cov'.format(POP_STRUC_DIR)
    log: LOG_DIR + '/pcangsd/allSamples_allChroms_{site}_maf{maf}_pcangsd.log'
    container: 'library://james-s-santangelo/pcangsd/pcangsd:0.99'
    threads: 10
    params:
        out = '{0}/pcangsd/allSamples_allChroms_{{site}}_maf{{maf}}_pcangsd'.format(POP_STRUC_DIR)
    wildcard_constraints:
        site='4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        python3 /opt/pcangsd-v.0.99/pcangsd.py \
            -beagle {input} \
            -o {params.out} \
            -threads {threads} \
            &> {log}
        """

rule ngsadmix:
    """
    Estimate admixture components for all samples
    """
    input:
        rules.concat_angsd_gl.output
    output:
        fopt = '{0}/ngsadmix/K{{k}}/ngsadmix_{{site}}_maf{{maf}}_K{{k}}_seed{{seed}}.fopt.gz'.format(POP_STRUC_DIR),
        qopt = '{0}/ngsadmix/K{{k}}/ngsadmix_{{site}}_maf{{maf}}_K{{k}}_seed{{seed}}.qopt'.format(POP_STRUC_DIR),
        lf = '{0}/ngsadmix/K{{k}}/ngsadmix_{{site}}_maf{{maf}}_K{{k}}_seed{{seed}}.log'.format(POP_STRUC_DIR)
    log: LOG_DIR + '/ngsadmix/{site}_maf{maf}_K{k}_seed{seed}_ngsadmix.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933' 
    threads: 10
    params:
        out = '{0}/ngsadmix/K{{k}}/ngsadmix_{{site}}_maf{{maf}}_K{{k}}_seed{{seed}}'.format(POP_STRUC_DIR)
    wildcard_constraints:
        site = '4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '02:00:00'
    shell:
        """
        NGSadmix -likes {input} \
            -K {wildcards.k} \
            -seed {wildcards.seed} \
            -P {threads} \
            -outfiles {params.out} 2> {log}
        """

rule logfile_for_clumpak:
    """
    Create Inputfile for CLUMPAK containing Log likelihood values of NGSadmix runs for each K
    """
    input:
        expand(rules.ngsadmix.output.lf, site='4fold', maf='0.05', k=NGSADMIX_K, seed=NGSADMIX_SEEDS)
    output:
        '{0}/clumpak/ngsadmix_logfile_for_clumpak.txt'.format(PROGRAM_RESOURCE_DIR)
    run:
        import re
        with open(output[0], 'w') as fout:
            for lf in input:
                # Get K
                m1 = re.search('(?<=_K)(\d+)', lf)
                k = m1.group(1)
                # Get likelihood
                line = open(lf, 'r').readlines()[-1]  # Likelihood always on last line
                m2 = re.search('(?<=like=)(-?\d+.\d+)', line)
                like = m2.group(1)
                fout.write('{0}\t{1}\n'.format(k, like))

rule clumpak_best_k_by_evanno:
    """
    Find optimal K value by city using Evanno method, as implemented in CLUMPAK
    """
    input:
        rules.logfile_for_clumpak.output
    output:
        directory('{0}/bestKbyEvanno'.format(POP_STRUC_DIR))
    log: LOG_DIR + '/clumpak_best_k_by_evanno/evanno.log'
    container: 'library://james-s-santangelo/clumpak/clumpak:1.1'
    params:
        outdir = '{0}/bestKbyEvanno'.format(POP_STRUC_DIR)
    resources:
        mem_mb = 1000,
        time = '01:00:00'
    shell:
        """
        perl /opt/bin/BestKByEvanno.pl --id clumpak_best_k_out \
            --d {params.outdir} \
            --f {input} \
            --inputtype lnprobbyk 2>&1 > {log}
        """

#####################
#### RELATEDNESS ####
#####################

# Estimate pariwise relatedness among all individuals 

rule create_bam_list_highQualSamples:
    """
    Create text files with paths to all BAM files (excluding low quality samples)
    """
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

rule pruned_degenerate_angsd_format:
    """
    Create ANGSD sites-formatted file with position of LD-pruned 4fold sites
    """
    input:
        rules.prune_degenerateSNPs_forPopStructure.output
    output:
        '{0}/angsd_sites/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_pruned.sites'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/pruned_degenerate_angsd_format/{chrom}_{site}_maf{maf}.log'
    shell:
        """
        sed 's/:/\t/g' {input} > {output} 2> {log}
        """

rule angsd_index_prunedSNPs:
    """
    Index LD pruned sites for ANGSD
    """
    input:
        rules.pruned_degenerate_angsd_format.output
    output:
        binary = '{0}/angsd_sites/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_pruned.sites.bin'.format(PROGRAM_RESOURCE_DIR),
        idx = '{0}/angsd_sites/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_pruned.sites.idx'.format(PROGRAM_RESOURCE_DIR)
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    log: LOG_DIR + '/angsd_index_prunedSNPs/{chrom}_{site}_maf{maf}_prunedIndex.log'
    shell:
        """
        angsd sites index {input} 2> {log}
        """
     
rule angsd_gl_forNGSrelate:
    """
    Estimate genotype likelihoods in binary format
    """
    input:
        bams = rules.create_bam_list_highQualSamples.output,
        ref = REFERENCE_GENOME,
        sites = rules.pruned_degenerate_angsd_format.output,
        idx = rules.angsd_index_prunedSNPs.output
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
        time = '3:00:00'
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

rule convert_freq_forNGSrelate:
    """
    Get allele frequencies for NGSrelate
    """
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
    """
    Estimate pairwise relatedness among samples from 4fold sites
    """
    input:
        bams = rules.create_bam_list_highQualSamples.output,
        gls = rules.angsd_gl_forNGSrelate.output.gls,
        freq = rules.convert_freq_forNGSrelate.output
    output:
        '{0}/ngsrelate/{{chrom}}_{{site}}_maf{{maf}}_NGSrelate.out'.format(POP_STRUC_DIR)
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

rule install_ggheatmap:
    """
    Install ggheatmap into Conda environment
    """
    output:
        f"{PROGRAM_RESOURCE_DIR}/ggheatmap_install.done"
    conda: "../envs/r.yaml"
    shell:
        """
        R -e 'install.packages("ggheatmap", repos = "http://cran.us.r-project.org")'
        R -e 'library(ggheatmap)' &&
        touch {output}
        """

rule population_structure_figures:
    """
    Generate population structure figures
    """
    input:
        ggheat_install = rules.install_ggheatmap.output,
        order = expand(rules.extract_sample_angsd.output, site="4fold"),
        cov = expand(rules.pcangsd.output, site="4fold", maf="0.05"),
        evanno = rules.clumpak_best_k_by_evanno.output,
        admix_log = expand(rules.ngsadmix.output.lf, k=NGSADMIX_K, site=['4fold'], maf=['0.05'], seed=NGSADMIX_SEEDS),
        admix_qopt = expand(rules.ngsadmix.output.qopt, k=NGSADMIX_K, site=['4fold'], maf=['0.05'], seed=NGSADMIX_SEEDS),
        pi_byHab = expand(rules.angsd_diversity_neutrality_stats_byHabitat.output, site="4fold", habitat=HABITATS),
        fst_byHab = expand(rules.angsd_habitat_fst_readable.output, site=['4fold'], hab_comb=HABITAT_COMBOS),
        bl = expand(rules.create_bam_list_highQualSamples.output, site="4fold"),
        pi_byPop = expand(rules.angsd_diversity_neutrality_stats_byPopulation.output, popu=POPS_MULTI_IND, site='4fold'),
        fst_byPop = expand(rules.angsd_population_fst_readable.output, pop_comb=POP_COMB_MULTI_IND, site='4fold')
    output:
        "test.txt",
        pca = f"{FIGURES_DIR}/pop_struct/pca_byHabitat.pdf",
        admix_optimal = f"{FIGURES_DIR}/pop_struct/admix_optimal.pdf",
        admix_optimal_minus = f"{FIGURES_DIR}/pop_struct/admix_optimal_minus.pdf",
        admix_optimal_plus = f"{FIGURES_DIR}/pop_struct/admix_optimal_plus.pdf",
        fst_byPop = f"{FIGURES_DIR}/pop_struct/fst_pairwise_population.pdf",
        ibd_plot = f"{FIGURES_DIR}/pop_struct/isolation_by_distance.pdf",
        pi_byHab_df = f"{FIGURES_DIR}/pop_struct/pi_byHab.txt",
        fst_byHab_df = f"{FIGURES_DIR}/pop_struct/fst_byHab.txt"
    conda:'../envs/r.yaml'
    notebook:
        "../notebooks/population_structure.r.ipynb"

##############
#### POST ####
##############

rule pop_structure_done:
    """
    Create empty file signaling completion of population structure
    """
    input:
        expand(rules.pcangsd.output, site=['4fold'], maf=['0.05']),
        expand(rules.ngsadmix.output, k=NGSADMIX_K, site=['4fold'], maf=['0.05'], seed=NGSADMIX_SEEDS),
        expand(rules.clumpak_best_k_by_evanno.output),
        expand(rules.ngsrelate.output, chrom=CHROMOSOMES, site='4fold', maf='0.05'),
        rules.population_structure_figures.output
    output:
        '{0}/population_structure.done'.format(POP_STRUC_DIR)
    shell:
        """
        touch {output}
        """
