# Rules to perform demographic modelling using dadi

rule dadi_sfs:
    input:
        unpack(get_dadi_sfs_input_files)
    output:
        '{0}/{{hab_comb}}_dadi_fromAngsd.sfs'.format(DADI_DIR),
    log: LOG_DIR + '/dadi_sfs/{hab_comb}_dadi_sfs.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    threads: 6
    shell:
        """
        realSFS dadi \
            {input.saf_urban} {input.saf_rural} \
            -sfs {input.sfs_urban} \
            -sfs {input.sfs_rural} \
            -ref {input.ref} \
            -anc {input.ref} \
            -seed 42 -P {threads} -maxIter 2000 > {output} 2> {log}
        """

rule format_dadi_sfs:
    input:
        rules.dadi_sfs.output
    output:
        '{0}/{{hab_comb}}_dadi_formatted.sfs'.format(DADI_DIR)
    params:
        pop1_n = 41,
        pop2_n = 41
    shell:
        """
        perl scripts/perl/realsfs2dadi.pl {input} {params.pop1_n} {params.pop2_n} > {output}
        """

rule run_dadi:
    input:
        sfs = rules.format_dadi_sfs.output
    output:
        logf = '{0}/{{hab_comb}}_pop0_pop1_{{rep}}.{{model}}.log.txt'.format(DADI_DIR),
        optimf = '{0}/{{hab_comb}}_pop0_pop1_{{rep}}.{{model}}.optimized.txt'.format(DADI_DIR)
    log: LOG_DIR + '/run_dadi/{hab_comb}_{model}_{rep}.log'
    conda: '../envs/dadi.yaml'
    params:
        prefix = '{0}/{{model}}/'.format(DADI_DIR)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    script:
        "../scripts/python/dadi_pipeline/dadi_Run_2D_Set.py"

rule summarize_dadi_output:
    input:
        expand(rules.run_dadi.output, hab_comb=['Urban_Rural'], model=['no_div', 'no_div_bot', 'no_div_growth', 'no_div_bot_growth'], rep = ['1', '2', '3', '4', '5'])
    output:
        short = '{0}/Results_Summary_Short.txt'.format(DADI_DIR),
        exten = '{0}/Results_Summary_Extended.txt'.format(DADI_DIR)
    params:
        path = '{0}/'.format(DADI_DIR)
    shell:
        """
        python3 scripts/python/dadi_pipeline/Summarize_Outputs.py {params.path}
        """

rule dadi_done:
    input:
        rules.summarize_dadi_output.output
    output:
        '{0}/dadi.done'.format(DADI_DIR)
    shell:
        """
        touch {output}
        """
