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


rule dadi_done:
    input:
        expand(rules.dadi_sfs.output, hab_comb=['Urban_Rural'])
    output:
        '{0}/dadi.done'.format(DADI_DIR)
    shell:
        """
        touch {output}
        """
