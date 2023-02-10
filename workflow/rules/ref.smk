rule samtools_index_reference:
    input:
        REFERENCE_GENOME
    output:
        f'{REFERENCE_GENOME}.fai'
    conda: '../envs/ref.yaml'
    log: LOG_DIR + '/samtools_index_reference/samtools_index_reference.log'
    shell:
        """
        samtools faidx {input} 2> {log}
        """

rule bwa_index_ref:
    input:
        REFERENCE_GENOME
    output:
        multiext(f'{REFERENCE_GENOME}', '.amb', '.ann', '.pac', '.0123', '.bwt.2bit.64')
    conda: '../envs/ref.yaml'
    log: LOG_DIR + '/bwa_index_ref/bwa_index_ref.log'
    resources:
        mem_mb = lambda wildcards, input, attempt: attempt * int(30 * input.size_mb),
        time = '01:00:00'
    shell:
        """
        bwa-mem2 index {input} 2> {log}
        """

rule makeblastdb_fromRef:
    input:
        REFERENCE_GENOME
    output:
        multiext(f'{REFERENCE_GENOME}', '.ndb', '.nhr', '.nin', '.nog', '.nos', '.not', '.nsq', '.ntf', '.nto') 
    conda: '../envs/ref.yaml'
    log: LOG_DIR + '/makeblastdb_fromReb/makeblastdb_fromRef.log'
    shell:
        """
        makeblastdb -in {input} \
            -dbtype nucl \
            -parse_seqids \
            -logfile test.log
        """

# rule clone_degeneracy:
#     """
#     Clone Degeneracy GitHub repo for getting 4fold and 0fold sites.
#     """
#     output:
#         temp(directory('Degeneracy'))
#     log: LOG_DIR + '/clone_degeneracy/clone_degeneracy.log'
#     shell:
#         """
#         git clone https://github.com/James-S-Santangelo/Degeneracy.git
#         """

# rule get_fourfold_zerofold:
#     """
#     Uses get_4fold_sites.sh from Degeneracy to get 4fold and 0fold sites across white clover
#     genome from reference sequence (FASTA) and annotation file (GFF).
#     """
#     input:
#         degen = rules.clone_degeneracy.output,
#         ref = REFERENCE_GENOME,
#         gff = GFF_FILE
#     output:
#         expand('{0}/4fold_0fold/Trepens_{{site}}.bed'.format(REF_DIR), site=['0fold','4fold'])
#     log: LOG_DIR + '/4fold_0fold/get_fourfold_zerofold.log'
#     conda: '../envs/ref.yaml'
#     params:
#         outpath = '{0}/4fold_0fold/'.format(REF_DIR)
#     resources:
#         mem_mb = 4000,
#         time = '06:00:00'
#     shell:
#         """
#         OUT_ABS=$( readlink -f {params.outpath} );
#         REF_ABS=$( readlink -f {input.ref} );
#         GFF_ABS=$( readlink -f {input.gff} )
#         ( cd {input.degen} &&
#         bash get_4fold_sites.sh $GFF_ABS $REF_ABS $OUT_ABS ) 2> {log}
#         """

rule ref_done:
    input:
        rules.samtools_index_reference.output,
        rules.bwa_index_ref.output
    output:
        '{0}/ref.done'.format(REF_DIR)
    shell:
        """
        touch {output}
        """
