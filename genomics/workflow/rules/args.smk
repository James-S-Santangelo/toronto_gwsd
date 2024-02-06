##############
### SETUP ####
##############

checkpoint create_regions_file_forARGs:
    input:
        ref_idx = rules.samtools_index_reference.output,
        ref = REFERENCE_GENOME,
        bams = rules.create_bam_lists_allFinalSamples_allSites.output
    output:
        directory(f'{PROGRAM_RESOURCE_DIR}/arg_regions')
    log: f'{LOG_DIR}/create_regions_files_forARGs/create_region_files_forARGs.log'
    params:
        region_size = 1_000_000,
        chroms = CHROMOSOMES,
        out = f'{PROGRAM_RESOURCE_DIR}/arg_regions/genome'
    shell:
        """
        mkdir -p {output}
        python3 scripts/python/fasta_generate_regions.py \
            {input.ref_idx} \
            {params.region_size} \
            --bed {params.out} \
            --chromosome {params.chroms} 2> {log}
        """

rule split_vcf_forARGs:
    input:
        vcf = rules.bcftools_concat_phased_vcfs.output.vcf,
        region = f'{PROGRAM_RESOURCE_DIR}/arg_regions/genome.{{chrom}}.region.{{region_id}}.bed'
    output:
        temp(f"{ARG_DIR}/vcfs/{{chrom}}/{{chrom}}_region{{region_id}}.vcf")
    log: f"{LOG_DIR}/split_vcf_forARGs/{{chrom}}/{{chrom}}_region{{region_id}}_vcf.log"
    conda: "../envs/args.yaml"
    shell:
        """
        tabix -h -R {input.region} {input.vcf} > {output}
        """

def get_vcfs(wildcards):
    ck_output = checkpoints.create_regions_file_forARGs.get(**wildcards).output[0]
    chrom, region_id = glob_wildcards(os.path.join(ck_output, "genome.{chrom}.region.{region_id}.bed"))
    vcfs = expand(f"{ARG_DIR}/vcfs/{{chrom}}/{{chrom}}_region{{region_id}}.vcf", zip, chrom=chrom, region_id=region_id)
    return(vcfs)

checkpoint write_nonempty_vcfs:
    input:
        get_vcfs
    output:
        directory(f"{ARG_DIR}/vcfs/nonempty")
    run:
        import shutil
        os.makedirs(output[0])
        for vcf in input:
            num_sites = sum(1 for l in open(vcf, "r").readlines() if not l.startswith("#"))
            if num_sites > 0:
                shutil.copy(vcf, output[0])

#############################
### SINGER ARG INFERENCE ####
#############################

rule singer_infer_arg:
    input:
        vcf = f"{ARG_DIR}/vcfs/nonempty/{{chrom}}_region{{region_id}}.vcf",
        region = f'{PROGRAM_RESOURCE_DIR}/arg_regions/genome.{{chrom}}.region.{{region_id}}.bed'
    output:
        log = f"{ARG_DIR}/trees/{{chrom}}/region{{region_id}}/region{{region_id}}.log",
    log: f"{LOG_DIR}/singer_infer_arg/{{chrom}}/{{chrom}}_region{{region_id}}_singer.log"
    conda: "../envs/args.yaml"
    params:
        out_prefix = f"{ARG_DIR}/trees/{{chrom}}/region{{region_id}}/region{{region_id}}",
        vcf_prefix = f"{ARG_DIR}/vcfs/{{chrom}}/{{chrom}}_region{{region_id}}",
        n_samples = 100
    shell:
        """
        START=$( cut -f2 {input.region} | cat );
        END=$( cut -f3 {input.region} | cat )
        ~/github-repos/SINGER/releases/singer-0.1.6-beta-linux-x86_64/singer_master \
            -vcf {params.vcf_prefix} \
            -m 1.44e-9 \
            -Ne 1.5e5 \
            -output {params.out_prefix} \
            -start $START \
            -end $END \
            -thin 10 \
            -n {params.n_samples} &> {log}
        """

rule convert_to_tskit:
    input:
        log = rules.singer_infer_arg.output
    output:
        done = f"{ARG_DIR}/trees/{{chrom}}/region{{region_id}}/region{{region_id}}_trees.done",
    log: f"{LOG_DIR}/convert_to_tskit/{{chrom}}/{{chrom}}_region{{region_id}}_tskit_convert.log"
    conda: "../envs/args.yaml"
    params:
        prefix = f"{ARG_DIR}/trees/{{chrom}}/region{{region_id}}/region{{region_id}}",
        n_samples = 100
    shell:
        """
        ( ~/github-repos/SINGER/releases/singer-0.1.6-beta-linux-x86_64/convert_to_tskit \
            -input {params.prefix} \
            -start 0 \
            -end 100 \
            -output {params.prefix} &&
          touch {output.done} ) &> {log}
        """

#####################################
#### FST FROM GENOTYPES AND ARGS ####
#####################################

rule generate_perSite_fst_estimates:
    input:
        trees = rules.convert_to_tskit.output,
        bams = rules.create_bam_lists_allFinalSamples_allSites.output,
        vcf = f"{ARG_DIR}/vcfs/nonempty/{{chrom}}_region{{region_id}}.vcf",
        sfs_fst = expand(rules.angsd_fst_allSites_readable.output, chrom=CHROMOSOMES, hab_comb="Urban_Rural"),
    output:
        f"{ARG_DIR}/fst/{{chrom}}/{{chrom}}_region{{region_id}}_windowed.fst"
    conda: "../envs/args.yaml"
    params:
        arg_path = ARG_DIR,
        n_samples = 100
    notebook:
        "../notebooks/generate_windowed_arg_gt_estimates.py.ipynb"

def get_all_ARGs(wildcards):
    ck_output = checkpoints.write_nonempty_vcfs.get(**wildcards).output[0]
    chrom, region_id = glob_wildcards(os.path.join(ck_output, "{chrom}_region{region_id}.vcf"))
    chroms = []
    region_ids = []
    for i, c in enumerate(chrom):
        if c == "Chr01_Occ":
            chroms.append(c)
            region_ids.append(region_id[i])
    fsts = expand(f"{ARG_DIR}/fst/{{chrom}}/{{chrom}}_region{{region_id}}_windowed.fst", chrom=chroms, region_id=region_ids)
    return fsts 

# rule fst_from_genotypes:
#     input:
#         vcf = rules.split_vcf_forARGs.output, 
#         samples = config["samples"]
#     output:
#         fst = f"{ARG_DIR}/vcftools/region{{region_id}}.weir.fst",
#         log = f"{ARG_DIR}/vcftools/region{{region_id}}.log"
#     conda: "../envs/args.yaml"
#     params:
#         win_size = 10000,
#         out = f"{ARG_DIR}/vcftools/region{{region_id}}" 
#     shell:
#         """
#         grep 'Urban' {input.samples} | cut -f1 > urban{wildcards.n}.txt
#         grep 'Rural' {input.samples} | cut -f1 > rural{wildcards.n}.txt
#         vcftools --vcf {input.vcf} \
#             --weir-fst-pop urban{wildcards.n}.txt \
#             --weir-fst-pop rural{wildcards.n}.txt \
#             --out {params.out} && 
#         rm urban{wildcards.n}.txt
#         rm rural{wildcards.n}.txt
#         """

# rule extract_gt_fst:
#     input:
#         expand(rules.fst_from_genotypes.output.log, n = [x for x in range(1, 363)])
#     output:
#         f"{ARG_DIR}/gt_fst.txt"
#     run:
#         regions = []
#         means = []
#         weighted = []
#         for f in input:
#             region = re.search("(\\d+)", os.path.basename(f)).group()
#             regions.append(region)
#             with open(f, "r") as fin:
#                 lines = fin.readlines()
#                 mean = re.search("(?<=:\\s)(.*$)", lines[-4].strip()).group()
#                 weight = re.search("(?<=:\\s)(.*$)", lines[-3].strip()).group()
#                 means.append(mean)
#                 weighted.append(weight)
#         df = pd.DataFrame(zip(regions, means, weighted), columns = ["regionID", "gt_fst_mean", "gt_fst_weighted"])
#         df.to_csv(output[0], sep="\t", index=False)
 
# rule get_nsites_ntrees_fromARGs:
#     input:
#         trees = lambda w: expand(rules.convert_to_tskit.output, n=w.n),
#     output:
#         f"{ARG_DIR}/nsites_ntrees/region{{region_id}}_nsites_ntrees.txt"
#     conda: "../envs/args.yaml"
#     params:
#         arg_path = ARG_DIR,
#         n_samples = 100,
#     script:
#         "../scripts/python/get_nsites_ntrees_fromARGs.py"

# rule analyse_args:
#     input:
#         arg_fst= expand(rules.calculate_fsts_fromARGs.output, n=[x for x in range(1, 363)]),
#         logs = expand(rules.singer_infer_arg.output.log, n=[x for x in range(1, 363)]),
#         gt_fst = rules.extract_gt_fst.output, 
#         sfs_fst = expand(rules.angsd_fst_allSites_readable.output, chrom=CHROMOSOMES, hab_comb="Urban_Rural"),
#         bams = rules.create_bam_lists_allFinalSamples_allSites.output,
#         regions = rules.create_regions_file_forARGs.output,
#         win_fst = expand(rules.generate_windowed_arg_gt_estimates.output, n=[x for x in range(1, 363)]),
#         nsites = expand(rules.get_nsites_ntrees_fromARGs.output, n=[x for x in range(1, 363)])
#     output:
#         "text.txt"
#     conda: "../envs/args.yaml"
#     notebook:
#         "../notebooks/analyse_args.r.ipynb"

rule args_done:
    input:
        get_all_ARGs
    output:
        f"{ARG_DIR}/args.done"
    shell:
        """
        touch {output}
        """
