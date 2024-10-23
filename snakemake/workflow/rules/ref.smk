rule samtools_index_reference:
    """
    Indexes reference genome with Samtools
    """
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
    """
    Indexes reference genome with BWA for read alignment
    """
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

rule degenotate:
    """
    Generate text file with degeneracy of every nucleotide in CDSs
    """
    input:
        ref = REFERENCE_GENOME,
        gff = GFF_FILE
    output:
        f'{REF_DIR}/degenotate/degeneracy-all-sites.bed',
    log: f'{LOG_DIR}/degenotate/degenotate.log'
    conda: '../envs/ref.yaml'
    params:
        outpath = f'{REF_DIR}/degenotate/'.format(REF_DIR)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        time = '03:00:00'
    shell:
        """
        sed 's/;$//' {input.gff} > tmp.gff3
        degenotate.py -a tmp.gff3 \
                -g {input.ref} \
                -o {params.outpath} \
                --overwrite &> {log}
        """

rule create_bed_from_degenotate:
    """
    Create BED file for 0fold and 4fold sites
    """
    input:
        rules.degenotate.output
    output:
        f'{REF_DIR}/Trepens_{{site}}.bed'
    params: fold = lambda w: '0' if w.site == '0fold' else '4'
    shell:
        """
        awk '$5 == {params.fold}' {input} | cut -f1,2,3 > {output}
        """

rule ref_done:
    """
    Create empty flag file signaling completion of reference genome setup
    """
    input:
        rules.samtools_index_reference.output,
        rules.bwa_index_ref.output
    output:
        '{0}/ref.done'.format(REF_DIR)
    shell:
        """
        touch {output}
        """
