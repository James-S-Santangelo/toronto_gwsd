rule create_region_files_forFreebayes:
    input:
        ref_idx = rules.samtools_index_reference.output,
        ref = rules.unzip_reference.output,
        bams = rules.create_bam_list_allFinalSamples.output
    output:
        expand('{0}/freebayes_regions/genome.{{chrom}}.region.{{i}}.bed'.format(PROGRAM_RESOURCE_DIR), chrom=CHROMOSOMES, i=FREEBAYES_CHUNKS)
    log: 'logs/create_regions_files_forFreebayes/create_region_files_forFreebayes.log'
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

# rule freebayes_call_variants:
#     input:
#         bams = rules.create_bam_list_allFinalSamples.output,
#         regions = rules.create_region_files_forFreebayes.output,
#         ref = rules.unzip_reference.output
#     output:
#         temp('{0}/vcf/{{chrom}}/{{chrom}}_{{node}}_allSamples.vcf'.format(FREEBAYES_DIR))
#     log: 'logs/freebayes/{chrom}__{node}_freebayes.log'
#     conda: '../envs/freebayes.yaml'
#     resources:
#         nodes = 1,
#         ntasks = CORES_PER_NODE,
#         time = '12:00:00'
#     shell:
#         """
#         ( freebayes-parallel {input.regions} {resources.ntasks} \
#             --fasta-reference {input.ref} \
#             --bam-list {input.bams} \
#             --use-best-n-alleles 2 \
#             --report-monomorphic > {output} ) 2> {log}
#         """
#  
# rule bgzip_vcf:
#     input:
#         rules.freebayes_call_variants.output
#     output:
#         temp('{0}/vcf/{{chrom}}/{{chrom}}_{{node}}_allSamples.vcf.gz'.format(FREEBAYES_DIR))
#     log: 'logs/bgzip/{chrom}_{node}_bgzip.log'
#     conda: '../envs/freebayes.yaml',
#     threads: CORES_PER_NODE
#     resources:
#         nodes = 1,
#         mem_mb = lambda wildcards, attempt: attempt * 500,
#         time = '01:00:00'
#     shell:
#         """
#         bgzip --threads {threads} --force {input} 
#         """
# 
# rule tabix_node_vcf:
#     input: 
#         rules.bgzip_vcf.output
#     output:
#         temp('{0}/vcf/{{chrom}}/{{chrom}}_{{node}}_allSamples.vcf.gz.tbi'.format(FREEBAYES_DIR))
#     log: 'logs/tabix_node_vcf/{chrom}_{node}_tabix.log'
#     conda: '../envs/freebayes.yaml'
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 4000,
#         time = '01:00:00'
#     shell:
#         """
#         tabix --force {input}
#         """
# 
# rule concat_vcfs:
#     input:
#         node_vcfs = get_node_vcfs,
#         node_indices = get_node_tabix_files
#     output:
#         '{0}/vcf/{{chrom}}/{{chrom}}_allSamples.vcf.gz'.format(FREEBAYES_DIR)
#     log: 'logs/concat_vcfs/{chrom}_concat_vcfs.log'
#     conda: '../envs/freebayes.yaml'
#     threads: 8
#     params:
#         nodes = NODES_PER_CHROM
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 4000,
#         time = '02:00:00'
#     shell:
#         """
#         if [[ {params.nodes} -eq 1 ]]
#         then
#             mv {input.node_vcfs} {output} 2> {log}
#         elif [[ {params.nodes} -gt 1 ]]
#         then
#             ( bcftools concat --allow-overlaps \
#                 --threads {threads} \
#                 --output-type z \
#                 --output {output} \
#                 {input.node_vcfs} ) 2> {log}
#         fi
#         """
# 
# rule bcftools_split_variants:
#     input:
#         vcf = rules.concat_vcfs.output,
#         tmp = rules.create_tmp_dir.output
#     output:
#         '{0}/vcf/{{chrom}}/{{chrom}}_allSamples_{{site_type}}_sorted.vcf.gz'.format(FREEBAYES_DIR)
#     log: 'logs/bcftools_split_variants/{chrom}_bcftools_split_variants_{site_type}.log'
#     conda: '../envs/freebayes.yaml'
#     wildcard_constraints:
#         site_type='snps|indels|invariant|mnps|other'
#     threads: 8
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 2000,
#         time = '06:00:00'
#     shell:
#         """
#         if [ {wildcards.site_type} = 'invariant' ]; then
#             bcftools view --threads {threads} -O z --include 'N_ALT = 0' {input.vcf} > {output} 2> {log}
#         elif [ {wildcards.site_type} = 'snps' ]; then
#             ( bcftools view --threads {threads} -O v --types {wildcards.site_type} {input.vcf} |\
#             vcfallelicprimitives --keep-info --keep-geno |\
#             bcftools view --threads {threads} --types {wildcards.site_type} --min-alleles 2 --max-alleles 2 |\
#             bcftools sort -O z -T {inpur.tmp}/{wildcards.chrom} -o {output} ) 2> {log}
#         else
#             bcftools view --threads {threads} -O z --types {wildcards.site_type} {input} > {output} 2> {log}
#         fi
#         """
# 
# rule tabix_vcf:
#     input:
#         rules.bcftools_split_variants.output
#     output:
#         '{0}/vcf/{{chrom}}/{{chrom}}_allSamples_{{site_type}}_sorted.vcf.gz.tbi'.format(FREEBAYES_DIR)
#     log: 'logs/tabix/{chrom}_tabix_{site_type}.log'
#     conda: '../envs/freebayes.yaml'
#     shell:
#         """
#         tabix {input}
#         """
# 
# rule vcf_to_zarr:
#     input:
#         rules.bcftools_split_variants.output
#     output:
#         directory('{0}/zarr_db/{{chrom}}/{{chrom}}_allSamples_{{site_type}}_sorted.zarr'.format(FREEBAYES_DIR))
#     log: 'logs/vcf_to_zarr/{chrom}_vcf_to_zarr_{site_type}.log'
#     conda: '../envs/freebayes.yaml'
#     wildcard_constraints:
#         site_type='snps|invariant'
#     threads: 1
#     resources:
#         time = '02:00:00'
#     script:
#         "../scripts/python/vcf_to_zarr.py"
# 
# rule freebayes_done:
#     input:
#         expand(rules.bcftools_split_variants.output, chrom=CHROMOSOMES, site_type=['snps','indels','invariant','mnps','other']),
#         expand(rules.tabix_vcf.output, chrom=CHROMOSOMES, site_type=['snps','indels','invariant','mnps','other']),
#         expand(rules.vcf_to_zarr.output, chrom=CHROMOSOMES, site_type=['snps','invariant'])
#     output:
#         '{0}/freebayes.done'.format(FREEBAYES_DIR)
#     shell:
#         """
#         touch {output}
#         """
