rule create_region_files_forFreebayes:
    input:
        ref_idx = rules.samtools_index_reference.output,
        ref = rules.unzip_reference.output,
        bams = rules.create_bam_list_allFinalSamples.output
    output:
        expand('{0}/freebayes_regions/genome.{{chrom}}.region.{{i}}.bed'.format(PROGRAM_RESOURCE_DIR), chrom=CHROMOSOMES, i=FREEBAYES_CHUNKS)
    log: LOG_DIR + '/create_regions_files_forFreebayes/create_region_files_forFreebayes.log'
    params:
        chroms = CHROMOSOMES,
        nchunks = FREEBAYES_NCHUNKS,
        out = '{0}/freebayes_regions/genome'.format(PROGRAM_RESOURCE_DIR)
    shell:
        """
        python3 scripts/python/fasta_generate_regions.py \
            {input.ref_idx} \
            {params.nchunks} \
            --chunks \
            --bed {params.out}\
            --chromosome {params.chroms} 2> {log}
        """

rule freebayes_call_variants:
    input:
        bams = rules.create_bam_list_allFinalSamples.output,
        region = '{0}/freebayes_regions/genome.{{chrom}}.region.{{i}}.bed'.format(PROGRAM_RESOURCE_DIR), 
        ref = rules.unzip_reference.output
    output:
        temp('{0}/vcf/{{chrom}}/{{chrom}}_chunk{{i}}_allSamples.vcf'.format(FREEBAYES_DIR))
    log: LOG_DIR + '/freebayes/{chrom}/{chrom}/{chrom}_chunk{i}_freebayes.log'
    conda: '../envs/freebayes.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = lambda wildcards, attempt: str(attempt * 3) + ":00:00"
    shell:
        """
        freebayes \
            --fasta-reference {input.ref} \
            --targets {input.region} \
            --bam-list {input.bams} \
            --use-best-n-alleles 2 \
            --skip-coverage 15000 \
            --report-monomorphic > {output} 2> {log}
        """
rule concat_vcfs:
    input:
        vcfs = get_vcfs_by_chrom,
        tmp = ancient(rules.create_tmp_dir.output)
    output:
        '{0}/vcf/{{chrom}}/{{chrom}}_allFinalSamples.vcf.gz'.format(FREEBAYES_DIR)
    log: LOG_DIR + '/concat_vcfs/{chrom}_concat_vcfs.log'
    conda: '../envs/freebayes.yaml'
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '03:00:00'
    shell:
        """
        ( bcftools concat --threads {threads} {input.vcfs} |\
            vcfuniq |\
            bgzip -c > {output} ) 2> {log}
        """

rule bcftools_split_variants:
    input:
        vcf = rules.concat_vcfs.output,
        tmp = ancient(rules.create_tmp_dir.output)
    output:
        '{0}/vcf/{{chrom}}/{{chrom}}_allFinalSamples_{{site_type}}.vcf.gz'.format(FREEBAYES_DIR)
    log: LOG_DIR + '/bcftools_split_variants/{chrom}_bcftools_split_variants_{site_type}.log'
    conda: '../envs/freebayes.yaml'
    wildcard_constraints:
        site_type='snps|invariant'
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time = '06:00:00'
    shell:
        """
        if [ {wildcards.site_type} = 'invariant' ]; then
            bcftools view --threads {threads} -O z --include 'N_ALT = 0' {input.vcf} > {output} 2> {log}
        elif [ {wildcards.site_type} = 'snps' ]; then
            ( bcftools view --threads {threads} -O v --types {wildcards.site_type} {input.vcf} |\
            vcfallelicprimitives --keep-info --keep-geno |\
            bcftools view --threads {threads} --types {wildcards.site_type} --min-alleles 2 --max-alleles 2 |\
            bcftools sort -O z -T {input.tmp}/{wildcards.chrom} -o {output} ) 2> {log}
        fi
        """

rule tabix_vcf:
    input:
        rules.bcftools_split_variants.output
    output:
        '{0}/vcf/{{chrom}}/{{chrom}}_allFinalSamples_{{site_type}}.vcf.gz.tbi'.format(FREEBAYES_DIR)
    log: LOG_DIR + '/tabix/{chrom}_tabix_{site_type}.log'
    conda: '../envs/freebayes.yaml'
    shell:
        """
        tabix {input}
        """

rule vcf_to_zarr:
    input:
        rules.bcftools_split_variants.output
    output:
        directory('{0}/zarr_db/{{chrom}}/{{chrom}}_allFinalSamples_{{site_type}}.zarr'.format(FREEBAYES_DIR))
    log: LOG_DIR + '/vcf_to_zarr/{chrom}_vcf_to_zarr_{site_type}.log'
    conda: '../envs/freebayes.yaml'
    wildcard_constraints:
        site_type='snps|invariant'
    threads: 1
    resources:
        mem_mb = 4000,
        time = '02:00:00'
    shell:
        """
        python3 scripts/python/vcf_to_zarr.py {input} {output} 2> {log}
        """

rule freebayes_done:
    input:
        expand(rules.bcftools_split_variants.output, chrom=CHROMOSOMES, site_type=['snps','invariant']),
        expand(rules.tabix_vcf.output, chrom=CHROMOSOMES, site_type=['snps','invariant']),
        expand(rules.vcf_to_zarr.output, chrom=CHROMOSOMES, site_type=['snps','invariant'])
    output:
        '{0}/freebayes.done'.format(FREEBAYES_DIR)
    shell:
        """
        touch {output}
        """
