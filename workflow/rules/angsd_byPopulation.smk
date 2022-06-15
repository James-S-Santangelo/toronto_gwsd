# Rules to estimate within-population diversity and pairwise Fst/PBS using ANGSD

rule create_bam_list_byPop_multiInd:
    input:
        rules.create_bam_list_allFinalSamples.output
    output:
        '{0}/bam_lists/by_population/{{popu}}_{{site}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: 'logs/create_bam_list/{popu}_{site}_bam_list.log'
    run:
        import os
        with open(output[0], 'w') as fout:
            with open(input[0], 'r') as fin:
                lines = fin.readlines()
                for line in lines:
                    sline = line.strip()
                    pop = os.path.basename(sline).split('_')[1]
                    if pop == wildcards.popu:
                        fout.write(line)

rule angsd_byPopulation_done:
    input:
        expand(rules.create_bam_list_byPop_multiInd.output, popu=POPS_MULTI_IND, site='4fold'),
    output:
        '{0}/angsd_byPopulation.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """
