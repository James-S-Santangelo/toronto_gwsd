# Rules to infer population structure and admixture proportions

###############
#### SETUP ####
###############

# Estimate LD among 4fold SNPs with MAF > 0.05. 
# Prune these SNPs within 20 Kb using r-squared cutoff of 0.2

rule create_pos_file_for_ngsLD:
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

rule pcangsd:
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

rule pop_structure_done:
    input:
        expand(rules.pcangsd.output, site=['4fold'], maf=['0.05']),
        expand(rules.ngsadmix.output, k=NGSADMIX_K, site=['4fold'], maf=['0.05'], seed=NGSADMIX_SEEDS),
        expand(rules.clumpak_best_k_by_evanno.output)   
    output:
        '{0}/population_structure.done'.format(POP_STRUC_DIR)
    shell:
        """
        touch {output}
        """
