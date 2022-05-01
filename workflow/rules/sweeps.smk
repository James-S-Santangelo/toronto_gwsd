# Rules to scan for signatures of selective sweeps

###############
#### SETUP ####
###############

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

rule genMap_toPlinkFormat:
    input:
        rules.split_genMap.output
    output:
        '{0}/{{chrom}}.map'.format(GENMAP_RESULTS_DIR)
    shell:
        """
        awk -F"\t" '{{print $2 " " $2"_"$1 " " $3 " " $1}}' {input} | tail -n +2 > {output}
        """

################
#### XP-NSL ####
################

rule selscan_xpnsl:
    input:
        unpack(selscan_xpnsl_input)
    output:
        temp('{0}/xpnsl/{{chrom}}/{{chrom}}_{{hab_comb}}.xpnsl.out'.format(SWEEPS_DIR))
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    log: LOG_DIR + '/selscan_xpnsl/{chrom}_{hab_comb}_xpnsl.log'
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    threads: 2
    params:
        out = '{0}/xpnsl/{{chrom}}/{{chrom}}_{{hab_comb}}'.format(SWEEPS_DIR) 
    shell:
        """
        selscan --xpnsl \
            --vcf {input.vcf} \
            --vcf-ref {input.vcf_ref} \
            --threads {threads} \
            --out {params.out} 2> {log}
        """

rule norm_xpnsl:
    input:
        lambda w: expand(rules.selscan_xpnsl.output, chrom=CHROMOSOMES, hab_comb=w.hab_comb)
    output:
        done = '{0}/xpnsl/{{hab_comb}}_xpnsl_normalization.done'.format(SWEEPS_DIR)
    log: LOG_DIR + '/norm_xpnsl/{hab_comb}_xpnsl_norm.log'
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    params:
        winsize = 50000
    shell:
        """
        norm --xpnsl --bp-win --winsize {params.winsize} --qbins 10 --files {input} 2> {log} &&
        touch {output}
        """

################
#### XP-CLR ####
################

rule xpclr:
    input:
       unpack(xpclr_input)
    output:
        '{0}/xpclr/{{chrom}}/{{chrom}}_{{hab_comb}}_xpclr.out'.format(SWEEPS_DIR)
    log: LOG_DIR + '/xpclr/{chrom}_{hab_comb}_xpcl.log'
    container: 'library://james-s-santangelo/xpclr/xpclr:1.2.1'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    params:
        size = 50000,
        step = 50000,
        maxsnps = 600,
        ld = 0.95
    shell:
        """
        POP1_SAMPLES=$( cut -f1 {input.pop1s} )
        POP2_SAMPLES=$( cut -f1 {input.pop2s} )
        xpclr --format vcf \
            --input {input.vcf} \
            --samplesA $POP1_SAMPLES \
            --samplesB $POP2_SAMPLES \
            --out {output} \
            --size {params.size} \
            --step {params.step} \
            --maxsnps {params.maxsnps} \
            --ld {params.ld} \
            --gdistkey CM \
            --phased \
            --chr {wildcards.chrom} 2> {log}
        """

###############
#### iHH12 ####
###############

rule ihh_OneTwo:
    input:
        vcf = rules.bcftools_splitVCF_byHabitat.output.vcf,
        genMap = rules.genMap_toPlinkFormat.output
    output:
        '{0}/ihh12/{{chrom}}/{{chrom}}_{{habitat}}.ihh12.out'.format(SWEEPS_DIR)
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    log: LOG_DIR + '/selscan_ihh12/{chrom}_{habitat}_ihh12.log'
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    threads: 2
    params:
        out = '{0}/ihh12/{{chrom}}/{{chrom}}_{{habitat}}'.format(SWEEPS_DIR) 
    shell:
        """
        selscan --ihh12 \
            --vcf {input.vcf} \
            --map {input.genMap} \
            --threads {threads} \
            --out {params.out} 2> {log}
        """

rule norm_ihh_OneTwo:
    input:
        lambda w: expand(rules.ihh_OneTwo.output, chrom=CHROMOSOMES, habitat=w.habitat)
    output:
        done = '{0}/ihh12/{{habitat}}_ihh12_normalization.done'.format(SWEEPS_DIR)
    log: LOG_DIR + '/norm_ihh12/{habitat}_ihh12_norm.log'
    container: 'library://james-s-santangelo/selscan/selscan:1.3.0'
    params:
        winsize = 50000
    shell:
        """
        norm --ihh12 --bp-win --winsize {params.winsize} --qbins 10 --files {input} 2> {log} &&
        touch {output}
        """

###############
#### RAiSD ####
###############

rule raisd:
    input:
        vcf = rules.bcftools_splitVCF_byHabitat.output.vcf
    output:
        info = '{0}/raisd/{{chrom}}/RAiSD_Info.{{chrom}}_{{habitat}}'.format(SWEEPS_DIR),
        report = '{0}/raisd/{{chrom}}/RAiSD_Report.{{chrom}}_{{habitat}}'.format(SWEEPS_DIR)
    container: 'library://james-s-santangelo/raisd/raisd:2.9'
    log: LOG_DIR + '/raisd/{chrom}_{habitat}_raisd.log'
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    params:
        name = '{0}/raisd/{{chrom}}/RAiSD_Info.{{chrom}}_{{habitat}}'.format(SWEEPS_DIR)
    shell:
        """
        RAiSD -n {wildcards.chrom}_{wildcards.habitat} \
            -I {input} \
            -R \
            -a 42 2> {log}
        """
        

rule sweeps_done:
    input:
        expand(rules.norm_xpnsl.output, hab_comb=['Urban_Rural', 'Rural_Suburban']),
        expand(rules.norm_ihh_OneTwo.output, habitat=HABITATS),
        #expand(rules.raisd.output, chrom='CM019101.1', habitat=HABITATS),
        expand(rules.xpclr.output, chrom=CHROMOSOMES, hab_comb=['Urban_Rural', 'Rural_Suburban'])
    output:
        '{0}/sweeps.done'.format(SWEEPS_DIR)
    shell:
        """
        touch {output}
        """
