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

rule clone_degeneracy:
    """
    Clone Degeneracy GitHub repo for getting 4fold and 0fold sites.
    """
    output:
        temp(directory('Degeneracy'))
    log: LOG_DIR + '/clone_degeneracy/clone_degeneracy.log'
    shell:
        """
        git clone https://github.com/James-S-Santangelo/Degeneracy.git
        """

rule get_fourfold_zerofold:
    """
    Uses get_4fold_sites.sh from Degeneracy to get 4fold and 0fold sites across white clover
    genome from reference sequence (FASTA) and annotation file (GFF).
    """
    input:
        degen = rules.clone_degeneracy.output,
        ref = REFERENCE_GENOME,
        gff = GFF_FILE
    output:
        '{0}/4fold_0fold/Trepens_{{chrom}}_0fold.bed'.format(REF_DIR),
        '{0}/4fold_0fold/Trepens_{{chrom}}_4fold.bed'.format(REF_DIR)
    log: LOG_DIR + '/4fold_0fold/{chrom}_get_fourfold_zerofold.log'
    conda: '../envs/ref.yaml'
    params:
        outpath = '{0}/4fold_0fold/'.format(REF_DIR)
    resources:
        mem_mb = 2000,
        time = '03:00:00'
    shell:
        """
        samtools faidx {input.ref} {wildcards.chrom} > {wildcards.chrom}_tmp.fasta
        grep '{wildcards.chrom}' {input.gff} > {wildcards.chrom}_tmp.gff

        OUT_ABS=$( readlink -f {params.outpath} );
        REF_ABS=$( readlink -f {wildcards.chrom}_tmp.fasta );
        GFF_ABS=$( readlink -f {wildcards.chrom}_tmp.gff )
        ( cd {input.degen} &&
        bash get_4fold_sites.sh $GFF_ABS $REF_ABS $OUT_ABS {wildcards.chrom} ) 2> {log}

        rm {wildcards.chrom}_tmp*
        """

rule concat_degenerate_sites:
    input:
        lambda w: expand(rules.get_fourfold_zerofold.output, chrom=CHROMOSOMES, site=w.site)
    output:
        f"{REF_DIR}/4fold_0fold/Trepens_{{site}}.bed"
    shell:
        """
        cat {input} > {output}
        """

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
