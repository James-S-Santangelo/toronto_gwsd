rule create_bam_list_allSamples:
    input:
        expand(rules.samtools_markdup.output.bam, sample=SAMPLES)
    output:
        '{0}/bam_lists/allSamples_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/create_bam_list/allSamples_bam_list.log'
    run:
        import os
        with open(output[0], 'w') as f:
            for bam in input:
                sample = os.path.basename(bam).split('_merged')[0]
                f.write('{0}\n'.format(bam))

rule bcftools_chloroplast_gene_variants:
    """
    Call variants in matK and rbcL genes using bcftools. Generates single VCF for each gene. 
    """
    input:
        bams = rules.create_bam_list_allSamples.output,
        ref = rules.unzip_reference.output
    output:
        '{0}/{{gene}}/allSamples_{{gene}}.vcf.gz'.format(SPECIES_ID_DIR)
    log: 'logs/chloroplast_gene_variants/chloroplast_{gene}_variants.log'
    conda: '../envs/species_id.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time = '00:30:00'
    params:
        region = lambda w: 'VCDJ01010680.1:6656-8083' if w.gene == 'rbcl' else 'VCDJ01010680.1:3872-5392'
    shell:
        """
        ( bcftools mpileup \
            --output-type u \
            --fasta-ref {input.ref} \
            --count-orphans \
            --regions {params.region} \
            --bam-list {input.bams} | \
            bcftools call \
            --output-type z \
            --variants-only \
            --multiallelic-caller \
            --ploidy 1 \
            --output {output} ) 2> {log}
        """

rule chloroplast_gene_fasta:
    """
    Extract FASTA sequence for matK and rbcL genes from reference genome using samtools faidx.
    One FASTA for each gene.
    """
    input:
       rules.unzip_reference.output
    output:
       '{0}/{{gene}}/{{gene}}.fna'.format(SPECIES_ID_DIR)
    log: 'logs/chloroplast_gene_fasta/{gene}_fasta.log'
    conda: '../envs/species_id.yaml'
    params:
        region = lambda w: 'VCDJ01010680.1:6656-8083' if w.gene == 'rbcl' else 'VCDJ01010680.1:3872-5392'
    shell:
        """
        samtools faidx {input} {params.region} > {output} 2> {log}
        """

rule index_chloroplast_gene_vcf:
    """
    Index VCF files with variants in each gene. Required for generating consensus
    """
    input:
        rules.bcftools_chloroplast_gene_variants.output
    output:
        '{0}/{{gene}}/allSamples_{{gene}}.vcf.gz.tbi'.format(SPECIES_ID_DIR)
    log: 'logs/index_chloroplast_gene_vcf/index_{gene}_vcf.log'
    conda: '../envs/species_id.yaml'
    shell:
        """
        tabix {input} 2> {log}
        """

rule chloroplast_gene_consensus:
    """
    Use genotype calls in VCFs to generate per-sample matK and rbcL consensus sequences. 
    Outputs two FASTAs per sample (one for each gene)
    """
    input:
        ref = rules.chloroplast_gene_fasta.output,
        vcf = rules.bcftools_chloroplast_gene_variants.output,
        idx = rules.index_chloroplast_gene_vcf.output
    output:
        temp('{0}/{{gene}}/consensus_fasta/{{sample}}_{{gene}}.fna'.format(SPECIES_ID_DIR))
    log: 'logs/{gene}_consensus/{sample}_{gene}_consensus.log'
    conda: '../envs/species_id.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1000,
        time = '00:30:00'
    shell:
        """
        bcftools consensus --fasta-ref {input.ref} \
            --sample {wildcards.sample} \
            {input.vcf} > {output} 2> {log}
        """

rule concat_fasta:
    """
    Concatenate per-sample consensus sequences into single FASTA file. One FASTA for each gene.
    """
    input:
        get_fastas_to_concat
    output:
        '{0}/{{gene}}/consensus_fasta/allSamples_{{gene}}.fasta'.format(SPECIES_ID_DIR)
    log: 'logs/concat_fasta/{gene}_concat_fastas.log'
    run:
        import os
        with open(output[0], 'w') as fout:
            for fasta in input:
                sample = os.path.basename(fasta).split('_{0}'.format(wildcards.gene))[0]
                with open(fasta, 'r') as fin:
                    lines = fin.readlines()
                    seq = ''.join(line.strip() for line in lines[1:])
                    fout.write('>{0};{1}\n{2}\n'.format(sample, wildcards.gene, seq))

checkpoint download_nt_database:
    """
    Download NCBI "nt" database. Uses checkpoint to accomodate changin database size. 
    """
    output:
        temp(directory('{0}/ncbi_nt_database/'.format(PROGRAM_RESOURCE_DIR)))
    log: 'logs/download_nt_database/download_nt_database.log'
    conda: '../envs/species_id.yaml'
    shell:
        """
        ( update_blastdb.pl nt &&
        mv nt* {output} ) 2> {log}
        """

rule extract_nt_database:
    input:
        '{0}/ncbi_nt_database/nt.{{db}}.tar.gz'.format(PROGRAM_RESOURCE_DIR)
    output:
        multiext('{0}/ncbi_nt_database/nt.{{db}}.'.format(PROGRAM_RESOURCE_DIR), 'nhd', 'nhi', 'nhr', 'nin', 'nnd', 'nni', 'nog', 'nsq')
    params:
        outdir = '{0}/ncbi_nt_database'.format(PROGRAM_RESOURCE_DIR)
    shell:
        """
        tar -xzf {input} --directory {params.outdir}
        """

rule ncbi_db_download_done:
    """
    Collects resutls of ncbi database download checkpoint and writes flag file indicating successful completion
    Flag file lists database files
    """
    input:
        aggregate_ncbi_input
    output:
        '{0}/ncbi_nt_database.done'.format(PROGRAM_RESOURCE_DIR)
    shell:
        """
        echo {input} > {output}
        """

rule blast_chloroplast_genes:
    """
    BLAST matK and rbcL FASTA sequences against local NCBI "nt" database using blastn. 
    Outputs single tab-separated text files for each gene.
    """
    input:
        fasta = rules.concat_fasta.output,
        done = rules.ncbi_db_download_done.output
    output:
        '{0}/{{gene}}/{{gene}}_blast_results.txt'.format(SPECIES_ID_DIR)
    log: 'logs/blast_chloroplast_genes/{gene}_blast.log'
    conda: '../envs/species_id.yaml'
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '03:00:00'
    params:
        db = '{0}/ncbi_nt_database/nt'.format(PROGRAM_RESOURCE_DIR),
        outfmt = "'6 qseqid qlen sseqid slen evalue bitscore lenth pident qcovhsp ssciname scomname'"
    shell:
        """
        blastn -query {input.fasta} \
            -db {params.db} \
            -out {output} \
            -num_threads {threads} \
            -max_hsps 1 \
            -max_target_seqs 1 \
            -outfmt {params.outfmt} 2> {log}
        """

rule species_id_done:
    """
    Collect final output files and write flag file signalling successful completion of species ID
    """
    input:
        expand(rules.concat_fasta.output, gene = ['rbcl','matk']),
        expand(rules.blast_chloroplast_genes.output, gene = ['rbcl','matk'])
    output:
        '{0}/species_id.done'.format(SPECIES_ID_DIR)
    shell:
        """
        touch {output}
        """
