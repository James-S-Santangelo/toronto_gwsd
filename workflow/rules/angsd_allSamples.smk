# Rules to estimate genotype likelihoods across all samples for population structure

###############
#### SETUP ####
###############

rule subset_bams_degeneracy:
    """
    Subset BAMs for all samples around 4fold sites to speed up ANGSD computations.
    """
    input:
        unpack(get_subset_bams_degeneracy_input)
    output:
        bam = '{0}/{{site}}/{{sample}}_{{site}}.bam'.format(BAM_DIR),
        idx = '{0}/{{site}}/{{sample}}_{{site}}.bam.bai'.format(BAM_DIR)
    log: LOG_DIR + '/subset_bams_degenerate/{sample}_{site}_subset.log'
    conda: '../envs/ref.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time = '01:00:00'
    shell:
        """
        ( samtools view -bh -L {input.regions} {input.bam} > {output.bam} &&\
            samtools index {output.bam} ) 2> {log}
        """

rule create_bam_list_allFinalSamples:
    input:
        bams = lambda wildcards: expand(rules.subset_bams_degeneracy.output.bam, sample=SAMPLES, site=wildcards.site),
        ref_flag = rules.ref_done.output 
    output:
        '{0}/bam_lists/allFinalSamples_{{site}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/create_bam_list/allFinalSamples_{site}_bam_list.log'
    run:
        import os
        with open(output[0], 'w') as f:
            for bam in input.bams:
                search = re.search('^(.+)(?=_\w)', os.path.basename(bam))
                sample = search.group(1) 
                if sample in FINAL_SAMPLES:
                    f.write('{0}\n'.format(bam))

rule convert_sites_for_angsd:
    input:
        get_bed_to_subset
    output:
        '{0}/angsd_sites/Trepens_{{site}}.sites'.format(PROGRAM_RESOURCE_DIR) 
    log: LOG_DIR + '/convert_sites_for_angsd/convert_{site}_for_angsd.log'
    wildcard_constraints:
        site='0fold|4fold'
    shell:
        """
        awk '{{print $1"\t"$2+1}}' {input} > {output} 2> {log}
        """

rule split_angsd_sites_byChrom:
    input:
        rules.convert_sites_for_angsd.output
    output:
        '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}.sites'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/split_angsd_sites_byChrom/{chrom}/{chrom}_{site}_split_angsd_sites.log'
    shell:
        """
        grep {wildcards.chrom} {input} > {output} 2> {log}
        """

rule angsd_index_sites_byChrom:
    input:
        rules.split_angsd_sites_byChrom.output
    output:
        binary = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}.sites.bin'.format(PROGRAM_RESOURCE_DIR),
        idx = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}.sites.idx'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/angsd_index_sites/{chrom}_{site}_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    shell:
        """
        angsd sites index {input} 2> {log}
        """

rule select_random_degenerate_sites:
    """
    Randomly select `params.nSites` degenerate (i.e., 4fold or 0fold) sites from across the genome
    """
    input:
        rules.convert_sites_for_angsd.output
    output:
        '{0}/angsd_sites/Trepens_{{site}}_random.sites'.format(PROGRAM_RESOURCE_DIR)
    params:
        nSites = 250000
    run:
        import random
        sites = open(input[0], 'r').readlines()
        random.seed(42)
        rand_sites = sorted(random.sample(sites, k = int(params.nSites)))
        with open(output[0], 'w') as fout:
            for site in rand_sites:
                fout.write(site)

rule angsd_index_random_degen_sites:
    """
    Index randomly selected genome-wide degenerate sites for use with ANGSD
    """
    input:
        rules.select_random_degenerate_sites.output
    output:
        binary = '{0}/angsd_sites/Trepens_{{site}}_random.sites.bin'.format(PROGRAM_RESOURCE_DIR),
        idx = '{0}/angsd_sites/Trepens_{{site}}_random.sites.idx'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/angsd_index_random_degen_sites/random_{site}_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    shell:
        """
        angsd sites index {input} 2> {log}
        """

rule split_random_angsd_sites_byChrom:
    """
    Split randomly selected degenerate ANGSD sites file into separate sites files by chromosome. Helps parallelize some computations.
    """
    input:
        rules.select_random_degenerate_sites.output
    output:
        sites = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}_random.sites'.format(PROGRAM_RESOURCE_DIR),
    log: LOG_DIR + '/split_random_angsd_sites_byChrom/{chrom}_{site}_split_angsd_sites_random.log'
    shell:
        """
        grep {wildcards.chrom} {input} > {output.sites} 2> {log}
        """

rule index_random_chromosomal_angsd_sites:
    """
    Index chromosomal ANGSD sites files for use with ANGSD
    """
    input:
        rules.split_random_angsd_sites_byChrom.output
    output:
        binary = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}_random.sites.bin'.format(PROGRAM_RESOURCE_DIR),
        idx = '{0}/angsd_sites/{{chrom}}/{{chrom}}_Trepens_{{site}}_random.sites.idx'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/index_random_chromosomal_angsd_sites/{chrom}_{site}_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933' 
    shell:
        """
        angsd sites index {input} 2> {log}
        """

##############################
#### GENOTYPE LIKELIHOODS ####
##############################

rule angsd_gl_degenerate_allSamples:
    input:
        bams = rules.create_bam_list_allFinalSamples.output,
        ref = rules.unzip_reference.output,
        sites = rules.split_random_angsd_sites_byChrom.output,
        idx = rules.index_random_chromosomal_angsd_sites.output
    output:
        gls = temp('{0}/gls/allSamples/{{site}}/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}.beagle.gz'.format(ANGSD_DIR)),
        mafs = temp('{0}/gls/allSamples/{{site}}/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}.mafs.gz'.format(ANGSD_DIR)),
    log: LOG_DIR + '/angsd_gl_allSamples_degenerate/{chrom}_{site}_maf{maf}_angsd_gl.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 8
    params:
        out = '{0}/gls/allSamples/{{site}}/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}'.format(ANGSD_DIR),
        max_dp = ANGSD_MAX_DP,
        min_dp_ind = ANGSD_MIN_DP_IND_GL
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '6:00:00'
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND*80/100 ))
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doGlf 2 \
            -doMajorMinor 1 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -doCounts 1 \
            -setMinDepthInd {params.min_dp_ind} \
            -setMaxDepth {params.max_dp} \
            -baq 2 \
            -ref {input.ref} \
            -minInd $MIN_IND \
            -minQ 20 \
            -minMapQ 30 \
            -sites {input.sites} \
            -minMaf {wildcards.maf} \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """

rule concat_angsd_gl:
    """
    Concatenated GLs from all 16 chromosomes into single file. Done separately for each site type.
    """
    input:
    	lambda wildcards: expand(rules.angsd_gl_degenerate_allSamples.output.gls, chrom=CHROMOSOMES, site=wildcards.site, maf=wildcards.maf)
    output:
        '{0}/gls/allSamples/{{site}}/allChroms_{{site}}_maf{{maf}}.beagle.gz'.format(ANGSD_DIR)
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

rule concat_angsd_mafs:
    """
    Concatenate MAF files for each of 16 chromosomes into single file. Done separately for each site type.
    """
    input:
    	lambda wildcards: expand(rules.angsd_gl_degenerate_allSamples.output.mafs, chrom=CHROMOSOMES, site=wildcards.site, maf=wildcards.maf)
    output:
        '{0}/gls/allSamples/{{site}}/allChroms_{{site}}_maf{{maf}}.mafs.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/concat_angsd_mafs/allSamples_{site}_{maf}_concat.log'
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


rule extract_sample_angsd:
    input:
        rules.create_bam_list_allFinalSamples.output
    output:
        '{0}/angsd_allFinalSamples_{{site}}_order.txt'.format(PROGRAM_RESOURCE_DIR)
    run:
        with open(output[0], 'w') as fout:
            with open(input[0], 'r') as fin:
                for line in fin:
                    sline = line.strip().split('/')
                    bam = sline[-1]
                    search = re.search('^(.+)(?=_\w)', bam)
                    sample = search.group(1) 
                    fout.write('{0}\n'.format(sample))

rule angsd_allSamples_done:
    input:
        expand(rules.concat_angsd_gl.output, site=['4fold'], maf=['0.05']),
        expand(rules.concat_angsd_mafs.output, site=['4fold'], maf=['0.05']),
        expand(rules.extract_sample_angsd.output, site=['4fold'])
    output:
        '{0}/angsd_allSamples.done'.format(ANGSD_DIR)
    shell:
        "touch {output}"
