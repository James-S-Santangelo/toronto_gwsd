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
            if num_sites >= 100:
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

        rm {params.prefix}_*.txt
        """

#################################
#### SUMMARY STATS FROM ARGS ####
#################################

rule generate_windowed_arg_summary_stats:
    input:
        trees = rules.convert_to_tskit.output,
        bams = rules.create_bam_lists_allFinalSamples_allSites.output,
        region = f'{PROGRAM_RESOURCE_DIR}/arg_regions/genome.{{chrom}}.region.{{region_id}}.bed',
        sfs_fst = expand(rules.angsd_fst_allSites_readable.output, chrom=CHROMOSOMES, hab_comb="Urban_Rural"),
    output:
        f"{ARG_DIR}/summary_stats/{{chrom}}/{{chrom}}_region{{region_id}}_win{{win_size}}_stats.txt"
    conda: "../envs/args.yaml"
    params:
        arg_path = ARG_DIR,
        n_samples = 100
    notebook:
        "../notebooks/generate_windowed_arg_summary_stats.py.ipynb"

######################################
#### SUMMARY STATS FROM GENOTYPES ####
######################################

rule create_pixy_popfile:
    input:
       config['samples']
    output:
        f"{PROGRAM_RESOURCE_DIR}/pixy/popfile.txt"
    run:
        with open(input[0], "r") as fin:
            with open(output[0], "w") as fout:
                lines = fin.readlines()
                for l in lines:
                    sl = l.split('\t')
                    sample = sl[0]
                    pop = sl[1]
                    if sample in FINAL_SAMPLES:
                        fout.write(f"{sample}\t{pop}\n")

def get_pixy_vcf(wildcards):
    if int(wildcards.win_size) > 1:
        vcf = rules.concat_variant_invariant_sites.output.vcf
        tbi = rules.concat_variant_invariant_sites.output.tbi
    else:
        vcf = f"{ARG_DIR}/vcfs/nonempty/{{chrom}}_region{{region_id}}.vcf"
    return vcf


rule pixy:
    input:
        vcf = get_pixy_vcf,
        popu = rules.create_pixy_popfile.output,
        region = f'{PROGRAM_RESOURCE_DIR}/arg_regions/genome.{{chrom}}.region.{{region_id}}.bed'
    output:
        fst = f"{PIXY_DIR}/{{chrom}}/{{chrom}}_miss{{miss}}_region{{region_id}}_win{{win_size}}_pixy_fst.txt",
        pi = f"{PIXY_DIR}/{{chrom}}/{{chrom}}_miss{{miss}}_region{{region_id}}_win{{win_size}}_pixy_pi.txt",
        dxy = f"{PIXY_DIR}/{{chrom}}/{{chrom}}_miss{{miss}}_region{{region_id}}_win{{win_size}}_pixy_dxy.txt"
    log: f"{LOG_DIR}/pixy/{{chrom}}_miss{{miss}}_region{{region_id}}_win{{win_size}}_pixy.log"
    conda: "../envs/pixy.yaml"
    params:
        start = lambda wildcards, input: int(open(input['region'], "r").readlines()[0].split("\t")[1]) + 1,
        end = lambda wildcards, input: int(open(input['region'], "r").readlines()[0].split("\t")[2].strip()) + 1,
        tmp_vcf = f"{{chrom}}{{miss}}{{region_id}}{{win_size}}_tmp.vcf.gz",
        tmp_sites = f"{{chrom}}{{miss}}{{region_id}}{{win_size}}.sites",
        out = f"{PIXY_DIR}/{{chrom}}",
        pref = f"{{chrom}}_miss{{miss}}_region{{region_id}}_win{{win_size}}_pixy"
    shell:
        """
        if [ {wildcards.win_size} = '1' ]; then
            bgzip -c {input.vcf} > {params.tmp_vcf};
            tabix {params.tmp_vcf};
            zgrep -v '#' {params.tmp_vcf} | cut -f1,2 > {params.tmp_sites}; 
            pixy --stats fst pi dxy \
                --population {input.popu} \
                --vcf {params.tmp_vcf} \
                --bypass_invariant_check yes \
                --window_size {wildcards.win_size} \
                --output_folder {params.out} \
                --output_prefix {params.pref} \
                --sites_file {params.tmp_sites} \
                --fst_type hudson &> {log}
            rm {params.tmp_vcf}*
            rm {params.tmp_sites}
        else
            pixy --stats fst pi dxy \
                --population {input.popu} \
                --vcf {input.vcf} \
                --interval_start {params.start} \
                --interval_end {params.end} \
                --window_size {wildcards.win_size} \
                --output_folder {params.out} \
                --output_prefix {params.pref} \
                --fst_type hudson &> {log}
        fi
        """

def get_all_ARGs(wildcards):
    ck_output = checkpoints.write_nonempty_vcfs.get(**wildcards).output[0]
    chrom, region_id = glob_wildcards(os.path.join(ck_output, "{chrom}_region{region_id}.vcf"))
    chroms = []
    regions = []
    for i, c in enumerate(chrom):
        if c == wildcards.chrom:
            chroms.append(c)
            regions.append(region_id[i])
    win_sizes = [wildcards.win_size for x in range(len(chroms))]
    fsts = expand(f"{ARG_DIR}/summary_stats/{{chrom}}/{{chrom}}_region{{region_id}}_win{{win_size}}_stats.txt",
                  zip, chrom=chroms, region_id=regions, win_size=win_sizes)
    return fsts 

def get_all_ARG_region_files(wildcards):
    ck_output = checkpoints.create_regions_file_forARGs.get(**wildcards).output[0]
    chrom, region_id = glob_wildcards(os.path.join(ck_output, "genome.{chrom}.region.{region_id}.bed"))
    chroms = []
    regions = []
    for i, c in enumerate(chrom):
        if c == wildcards.chrom:
            chroms.append(c)
            regions.append(region_id[i])
    regions = expand(f'{PROGRAM_RESOURCE_DIR}/arg_regions/genome.{{chrom}}.region.{{region_id}}.bed',
                     zip, chrom=chroms, region_id=regions)
    return regions 

def get_all_pixy_files(wildcards):
    ck_output = checkpoints.write_nonempty_vcfs.get(**wildcards).output[0]
    chrom, region_id = glob_wildcards(os.path.join(ck_output, "{chrom}_region{region_id}.vcf"))
    chroms = []
    regions = []
    for i, c in enumerate(chrom):
        if c == wildcards.chrom:
            chroms.append(c)
            regions.append(region_id[i])
    win_sizes = [wildcards.win_size for x in range(len(chroms))]
    misses = [wildcards.miss for x in range(len(chroms))]
    fsts = expand(f"{PIXY_DIR}/{{chrom}}/{{chrom}}_miss{{miss}}_region{{region_id}}_win{{win_size}}_pixy_fst.txt",
                  zip, chrom=chroms, region_id=regions, miss=misses, win_size=win_sizes)
    return fsts 

##################
#### ANALYSES ####
##################

rule write_all_fsts:
    input:
        arg_fst = get_all_ARGs,
        regions = get_all_ARG_region_files,
        gt_fsts = get_all_pixy_files,
        sfs_fsts = lambda w: expand(rules.angsd_fst_allSites_readable.output, chrom=w.chrom, hab_comb="Urban_Rural")
    output:
        f"{ARG_DIR}/summary_stats/{{chrom}}/{{chrom}}_all_fsts_miss{{miss}}_win{{win_size}}.txt"    
    conda: "../envs/args.yaml"
    script:
        "../scripts/r/write_all_fsts.R"

def get_all_ntree_nsite_files(wildcards):
    ck_output = checkpoints.write_nonempty_vcfs.get(**wildcards).output[0]
    chrom, region_id = glob_wildcards(os.path.join(ck_output, "{chrom}_region{region_id}.vcf"))
    files = expand(f"{ARG_DIR}/nsites_ntrees/{{chrom}}/{{chrom}}_region{{region_id}}_nsites_ntrees.txt",
                  zip, chrom=chrom, region_id=region_id)
    return files 

rule plot_arg_gt_fst_correlations:
    input:
        nSites_nTrees = get_all_ntree_nsite_files,
        all_fsts = expand(rules.write_all_fsts.output, chrom=CHROMOSOMES, miss="0", win_size=["1", "10000"])
    output:
        site_fst_winSize_hist = f"{ARG_DIR}/figures/site_fst_hist_by_winSize.pdf",
        site_fst_winSize_hist_filt = f"{ARG_DIR}/figures/site_fst_hist_by_winSize_filtered.pdf",
        branch_fst_winSize_hist_filt = f"{ARG_DIR}/figures/branch_fst_hist_by_winSize_filtered.pdf",
        site_gt_fst_cor_by_winSize= f"{ARG_DIR}/figures/site_gt_fst_cor_by_winSize.pdf",
        branch_gt_fst_cor_by_winSize= f"{ARG_DIR}/figures/branch_gt_fst_cor_by_winSize.pdf",
        site_sfs_fst_cor_by_winSize= f"{ARG_DIR}/figures/site_sfs_fst_cor_by_winSize.pdf",
        branch_sfs_fst_cor_by_winSize= f"{ARG_DIR}/figures/branch_sfs_fst_cor_by_winSize.pdf",
        branch_gt_cor_hist_by_winSize = f"{ARG_DIR}/figures/branch_gt_cor_hist_by_winSize.pdf",
        site_gt_cor_hist_by_winSize = f"{ARG_DIR}/figures/site_gt_cor_hist_by_winSize.pdf",
        branch_sfs_cor_hist_by_winSize = f"{ARG_DIR}/figures/branch_sfs_cor_hist_by_winSize.pdf",
        site_sfs_cor_hist_by_winSize = f"{ARG_DIR}/figures/site_sfs_cor_hist_by_winSize.pdf",
        nTrees_nSites_hist = f"{ARG_DIR}/figures/nSites_nTrees_hist.pdf"
    conda: "../envs/args.yaml"
    notebook:
        "../notebooks/plot_arg_gt_fst_correlations.r.ipynb"

rule get_nsites_ntrees_fromARGs:
    input:
        trees = lambda w: expand(rules.convert_to_tskit.output, chrom=w.chrom, region_id=w.region_id)
    output:
        f"{ARG_DIR}/nsites_ntrees/{{chrom}}/{{chrom}}_region{{region_id}}_nsites_ntrees.txt"
    conda: "../envs/args.yaml"
    params:
        arg_path = ARG_DIR,
        n_samples = 100,
    script:
        "../scripts/python/get_nsites_ntrees_fromARGs.py"

rule args_done:
    input:
        rules.plot_arg_gt_fst_correlations.output
    output:
        f"{ARG_DIR}/args.done"
    shell:
        """
        touch {output}
        """
