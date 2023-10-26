rule create_regions_file_forARGs:
    input:
        xpnsl = rules.write_windowed_statistics.output.xpnsl_df
    output:
        arg_regions = f'{PROGRAM_RESOURCE_DIR}/args/args_regions.txt'
    conda: '../envs/args.yaml'
    script:
        "../scripts/r/write_arg_regions.R"

rule split_arg_regions:
    input:
        rules.create_regions_file_forARGs.output
    output:
        f"{PROGRAM_RESOURCE_DIR}/args/region{{n}}.bed"
    run:
        df = pd.read_csv(input[0], sep = "\t")
        df_filt = df[df["regionID"] == int(wildcards.n)]
        chr = df_filt["Chr"].iloc[0]
        start = df_filt["start"].iloc[0]
        end = df_filt["end"].iloc[0]
        with open(output[0], "w") as fout:
            fout.write(f"{chr}\t{start}\t{end}")

rule split_vcf_forARGs:
    input:
        vcf = rules.bcftools_concat_phased_vcfs.output.vcf,
        region = rules.split_arg_regions.output
    output:
        f"{ARG_DIR}/vcfs/region{{n}}.vcf"
    log: f"{LOG_DIR}/split_vcf_forARGs/region{{n}}_vcf.log"
    conda: "../envs/args.yaml"
    shell:
        """
        tabix -h -R {input.region} {input.vcf} > {output} 2> {log}
        """

ITERS = [i for i in range(0, 500)]
rule singer_infer_arg:
    input:
        vcf = rules.split_vcf_forARGs.output,
        region = rules.split_arg_regions.output
    output:
        log = f"{ARG_DIR}/trees/region{{n}}/region{{n}}.log",
        recombs = expand(f"{ARG_DIR}/trees/region{{n}}/region{{n}}_recombs_{{i}}.txt", i=ITERS, allow_missing=True),
        muts = expand(f"{ARG_DIR}/trees/region{{n}}/region{{n}}_muts_{{i}}.txt", i=ITERS, allow_missing=True),
        nodes = expand(f"{ARG_DIR}/trees/region{{n}}/region{{n}}_nodes_{{i}}.txt", i=ITERS, allow_missing=True),
    log: f"{LOG_DIR}/singer_infer_arg/region{{n}}_singer.log"
    params:
        out_prefix = f"{ARG_DIR}/trees/region{{n}}/region{{n}}",
        vcf_prefix = f"{ARG_DIR}/vcfs/region{{n}}"
    shell:
        """
        START=$( cut -f2 {input.region} | cat );
        END=$( cut -f3 {input.region} | cat )
        ~/github-repos/SINGER/previous_releases/singer_master \
            -vcf {params.vcf_prefix} \
            -m 1.8e-8 \
            -Ne 2.5e5 \
            -output {params.out_prefix} \
            -start $START \
            -end $END \
            -thin 1 \
            -n 500 &> {log}
        """

rule convert_singer_arg_to_tskit:
    input:
        rules.singer_infer_arg.output.log
    output:
        trees = expand(f"{ARG_DIR}/trees/region{{n}}/region{{n}}_{{i}}.trees", i=ITERS, allow_missing=True),
    log: f"{LOG_DIR}/convert_singer_arg_to_tskit/region{{n}}_convert.log"
    conda: "../envs/args.yaml"
    params:
        prefix = f"{ARG_DIR}/trees/region{{n}}/region{{n}}"
    shell:
        """
        ~/github-repos/SINGER/previous_releases/convert_to_tskit \
            -log {params.prefix} \
            -output {params.prefix} 2> {log}
        """
                
FAILED_SAMPLES = [37, 145, 317, 38, 146, 362, 69, 35, 343, 270]
rule calculate_fsts_fromARGs:
    input:
        trees = lambda w: expand(rules.convert_singer_arg_to_tskit.output, n=w.n),
        bams = rules.create_bam_lists_allFinalSamples_allSites.output,
    output:
        f"{ARG_DIR}/fst/region{{n}}.fst"
    conda: "../envs/args.yaml"
    script:
        "../scripts/python/calculate_fsts_fromARGs.py"

rule fst_from_genotypes:
    input:
        vcf = rules.split_vcf_forARGs.output,
        samples = config["samples"]
    output:
        fst = f"{ARG_DIR}/vcftools/region{{n}}.weir.fst",
        log = f"{ARG_DIR}/vcftools/region{{n}}.log"
    conda: "../envs/args.yaml"
    params:
        out = f"{ARG_DIR}/vcftools/region{{n}}" 
    shell:
        """
        grep 'Urban' {input.samples} | cut -f1 > urban{wildcards.n}.txt
        grep 'Rural' {input.samples} | cut -f1 > rural{wildcards.n}.txt
        vcftools --vcf {input.vcf} \
            --weir-fst-pop urban{wildcards.n}.txt \
            --weir-fst-pop rural{wildcards.n}.txt \
            --out {params.out} && 
        rm urban{wildcards.n}.txt
        rm rural{wildcards.n}.txt
        """

rule extract_gt_fst:
    input:
        expand(rules.fst_from_genotypes.output.log, n = [x for x in range(1, 363) if x not in FAILED_SAMPLES])
    output:
        f"{ARG_DIR}/gt_fst.txt"
    run:
        regions = []
        means = []
        weighted = []
        for f in input:
            print(f)
            region = re.search("(\\d+)", os.path.basename(f)).group()
            regions.append(region)
            with open(f, "r") as fin:
                lines = fin.readlines()
                mean = re.search("(?<=:\\s)(.*$)", lines[-4].strip()).group()
                weight = re.search("(?<=:\\s)(.*$)", lines[-3].strip()).group()
                means.append(mean)
                weighted.append(weight)
        df = pd.DataFrame(zip(regions, means, weighted), columns = ["regionID", "gt_fst_mean", "gt_fst_weighted"])
        df.to_csv(output[0], sep="\t", index=False)

rule analyse_args:
    input:
        arg_fst= expand(rules.calculate_fsts_fromARGs.output, n=[x for x in range(1, 363) if x not in FAILED_SAMPLES]),
        logs = expand(rules.singer_infer_arg.output.log, n=[x for x in range(1, 363) if x not in FAILED_SAMPLES]),
        gt_fst = rules.extract_gt_fst.output, 
        sfs_fst = expand(rules.angsd_fst_allSites_readable.output, chrom=CHROMOSOMES, hab_comb="Urban_Rural"),
        bams = rules.create_bam_lists_allFinalSamples_allSites.output,
        regions = rules.create_regions_file_forARGs.output
    output:
        "test.txt"
    conda: "../envs/args.yaml"
    notebook:
        "../notebooks/analyse_args.r.ipynb"

rule args_done:
    input:
        expand(rules.calculate_fsts_fromARGs.output, n=[x for x in range(1, 363) if x not in FAILED_SAMPLES]),
    output:
        f"{ARG_DIR}/args.done"
    shell:
        """
        touch {output}
        """
