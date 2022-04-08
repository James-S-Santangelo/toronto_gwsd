# Rules to scan for signatures of selective sweeps

rule samples_byHabitat:
    input:
        samples = config['samples']
    output:
        '{0}/selscan/{{chrom}}/{{chrom}}_{{habitat}}.samples'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/samples_byHabitat/{chrom}_{habitat}_samples.log'
    run:
        import os
        import logging
        import pandas as pd
        logging.basicConfig(filename=log[0], level=logging.DEBUG)
        try:
            df = pd.read_table(input.samples, sep = '\t')
            df_sub = df[(df['Habitat'] == wildcards.habitat)]
            samples_habitat = df_sub['Sample'].tolist()
            samples_toWrite = []
            outname = wildcards.chrom + '_' + wildcards.habitat + '_whatshapPhased_shapeitPhased'
            for sample in samples_habitat:
                if sample in FINAL_SAMPLES:
                    samples_toWrite.append(sample)
            with open(output[0], 'w') as f:
                for sample in samples_toWrite:
                    if sample == samples_toWrite[-1]:
                        f.write('{0}\t-\t{1}'.format(sample, outname))
                    else:
                        f.write('{0},'.format(sample))
        except:
            logging.exception("An error occured!") 
            raise
       
rule bcftools_splitVCF_byHabitat:
    input:
        vcf = rules.shapeit_phase.output.vcf,
        samples = lambda w: expand(rules.samples_byHabitat.output, habitat = w.habitat, chrom = w.chrom)
    output:
        vcf = '{0}/vcf/{{chrom}}/{{chrom}}_{{habitat}}_whatshapPhased_shapeitPhased.vcf.gz'.format(FREEBAYES_DIR),
        idx = '{0}/vcf/{{chrom}}/{{chrom}}_{{habitat}}_whatshapPhased_shapeitPhased.vcf.gz.tbi'.format(FREEBAYES_DIR)
    log: LOG_DIR + '/bcftools_splitVCF_byHabitat/{chrom}_{habitat}_split.log'
    conda: '../envs/sweeps.yaml',
    params:
        out =  '{0}/vcf/{{chrom}}/'.format(FREEBAYES_DIR)
    shell:
        """
        ( bcftools +split {input.vcf} \
            --samples-file {input.samples} \
            --output-type z \
            --output {params.out} && tabix {output.vcf} ) 2> {log}
        """

rule format_genMap_forSelscan:
    input:
        rules.split_genMap.output
    output:
        '{0}/{{chrom}}_forSelscan.map'.format(GENMAP_RESULTS_DIR)
    shell:
        """
        awk -F"\t" '{{print $2 " " "rs"$1 " " $3 " " $1}}' {input} | tail -n +2 > {output}
        """

rule selscan_xpehh:
    input:
        vcf = lambda w: expand(rules.bcftools_splitVCF_byHabitat.output.vcf, chrom=w.chrom, habitat='Urban'),
        vcf_ref = lambda w: expand(rules.bcftools_splitVCF_byHabitat.output.vcf, chrom=w.chrom, habitat='Rural'),
        genMap = rules.format_genMap_forSelscan.output
    output:
        temp('{0}/{{chrom}}/{{chrom}}_UR.xpehh.out'.format(SWEEPS_DIR))
    conda: '../envs/sweeps.yaml'
    log: LOG_DIR + '/selscan_xpehh/{chrom}_UR_xpehh.log'
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    threads: 6
    params:
        out = '{0}/{{chrom}}/{{chrom}}_UR'.format(SWEEPS_DIR) 
    shell:
        """
        selscan --xpehh \
            --vcf {input.vcf} \
            --vcf-ref {input.vcf_ref} \
            --map {input.genMap} \
            --threads {threads} \
            --out {params.out} 2> {log}
        """

rule norm_xpehh:
    input:
        rules.selscan_xpehh.output
    output:
        norm = '{0}/{{chrom}}/{{chrom}}_UR.xpehh.out.norm'.format(SWEEPS_DIR),
        logf = '{0}/{{chrom}}/{{chrom}}_UR.xpehh.out.log'.format(SWEEPS_DIR)
    log: LOG_DIR + '/norm_xpehh/{chrom}_xpehh_norm.log'
    conda: '../envs/sweeps.yaml'
    shell:
        """
        norm --xpehh --files {input} --log {output.logf} 2> {log}
        """
    

rule sweeps_done:
    input:
        expand(rules.norm_xpehh.output, chrom = CHROMOSOMES, habitat = HABITATS)
    output:
        '{0}/sweeps.done'.format(SWEEPS_DIR)
    shell:
        """
        touch {output}
        """
