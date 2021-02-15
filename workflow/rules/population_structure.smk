rule pcangsd:
    input:
        rules.concat_angsd_gl.output
    output:
        '{0}/pcangsd/{{sample_set}}/{{sample_set}}_{{site}}_maf{{maf}}_pcangsd.cov'.format(POP_STRUC_DIR),
        '{0}/pcangsd/{{sample_set}}/{{sample_set}}_{{site}}_maf{{maf}}_pcangsd.admix.Q.npy'.format(POP_STRUC_DIR)
    log: 'logs/pcangsd/{sample_set}_{site}_maf{maf}_pcangsd.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:pcangsd_v0.99'
    threads: 10
    params:
        out = '{0}/pcangsd/{{sample_set}}/{{sample_set}}_{{site}}_maf{{maf}}_pcangsd'.format(POP_STRUC_DIR)
    wildcard_constraints:
        site='4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        python3 /opt/pcangsd/pcangsd.py \
            -beagle {input} \
            -o {params.out} \
            -threads {threads} \
            &> {log}
        """

rule ngsadmix:
    input:
        rules.concat_angsd_gl.output
    output:
        '{0}/ngsadmix/{{sample_set}}/K{{k}}/{{sample_set}}_ngsadmix_{{site}}_maf{{maf}}_K{{k}}_seed{{seed}}.fopt.gz'.format(POP_STRUC_DIR),
        '{0}/ngsadmix/{{sample_set}}/K{{k}}/{{sample_set}}_ngsadmix_{{site}}_maf{{maf}}_K{{k}}_seed{{seed}}.qopt'.format(POP_STRUC_DIR)
    log: 'logs/ngsadmix/{sample_set}_{site}_maf{maf}_K{k}_seed{seed}_ngsadmix.log'
    container: 'shub://James-S-Santangelo/singularity-recipes:angsd_v0.933'
    threads: 10
    params:
        out = '{0}/ngsadmix/{{sample_set}}/K{{k}}/{{sample_set}}_ngsadmix_{{site}}_maf{{maf}}_K{{k}}_seed{{seed}}'.format(POP_STRUC_DIR)
    wildcard_constraints:
        site = '4fold',
        sample_set = 'finalSamples_relatedRemoved'
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

rule pop_structure_done:
    input:
        expand(rules.pcangsd.output, site=['4fold'], maf=['0.05'], sample_set=['highQualSamples','finalSamples_relatedRemoved']),
        expand(rules.ngsadmix.output, k=NGSADMIX_K, site=['4fold'], maf=['0.05'], seed=NGSADMIX_SEEDS, sample_set=['finalSamples_relatedRemoved'])
    output:
        '{0}/population_structure.done'.format(POP_STRUC_DIR)
    shell:
        """
        touch {output}
        """


# rule create_pos_file_for_ngsLD:
#     input:
#         rules.angsd_gl_allSites.output.mafs
#     output:
#         '{0}/ngsld_pos/{{chrom}}_angsdGL_withMaf{{maf}}.pos'.format(PROGRAM_RESOURCE_DIR)
#     log: 'logs/create_pos_file_for_ngsLD/{chrom}_withMaf{maf}_pos.log'
#     shell:
#         """
#         zcat {input} | cut -f 1,2 | tail -n +2 > {output} 2> {log}
#         """
# 
# rule calc_ld_angsd_gl:
#     input:
#         pos = rules.create_pos_file_for_ngsLD.output,
#         gls = rules.angsd_gl_allSites.output.gls
#     output:
#         '{0}/{{chrom}}/{{chrom}}_allSamples_withMaf{{maf}}.ld.gz'.format(NGSLD_DIR)
#     log: 'logs/calc_ld_angsd_gl/{chrom}_withMaf{maf}_calc_ld.log'
#     container: 'shub://James-S-Santangelo/singularity-recipes:ngsld_v1.1.1'
#     threads: 16
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 8000,
#         time = '16:00:00'
#     shell:
#         """
#         ( NUM_SITES=$(cat {{input.pos}} | wc -l) &&
#           ngsLD --geno {{input.gls}} \
#             --pos {{input.pos}} \
#             --n_ind {0} \
#             --n_sites $NUM_SITES \
#             --probs \
#             --min_maf {{wildcards.maf}} \
#             --n_threads {{threads}} \
#             --max_kb_dist 100 | gzip --best > {{output}} ) 2> {{log}}
#         """.format(len(SAMPLES))
# 
# rule prune_ld:
#     input:
#         rules.calc_ld_angsd_gl.output
#     output:
#         '{0}/pruned/{{chrom}}/{{chrom}}_withMaf{{maf}}_pruned.id'.format(NGSLD_DIR)
#     log: 'logs/prune_ld/{chrom}_withMaf{maf}_prune_ld.log'
#     container: 'shub://James-S-Santangelo/singularity-recipes:ngsld_v1.1.1'
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 50000,
#         time = '72:00:00'
#     shell:
#         """
#         ( zcat {input} | perl /opt/bin/prune_graph.pl \
#             --max_kb_dist 20 \
#             --min_weight 0.2 | sort -V > {output} ) 2> {log}
#         """
