# Rules for fitting and interpolating recombination maps
# Relies on genetic positions from Olsen et al. 2021 and GBS marker sequences sent by Olsen

rule blast_markers:
    input:
        markers = '{0}/{{map_pop}}_tags.fa'.format(GENMAP_RESOURCE_DIR),
        db_flag = rules.makeblastdb_fromRef.output,
        ref = rules.unzip_reference.output
    output:
        '{0}/{{map_pop}}_marker_blast.txt'.format(GENMAP_RESULTS_DIR)
    conda: '../envs/recombination_map.yaml'
    log: LOG_DIR + '/blast_markers/{map_pop}_marker_blast.log'
    params: 
        outfmt = "'6 qseqid sseqid pident length mismatch gamap_popen qstart qend qlen sstart send slen evalue bitscore qcovs qcovhsp'"
    shell:
        """
        blastn -db {input.ref} \
            -query {input.markers} \
            -out {output} \
            -outfmt {params.outfmt} \
            -num_threads {threads} \
            -evalue 1e-10 \
            -max_hsps 5 \
            -max_target_seqs 5 2> {log}  
        """

rule recombination_map_done:
    input:
        expand(rules.blast_markers.output, map_pop = ['SG', 'DG'])
    output:
        '{0}/recombination_map_done'.format(GENMAP_RESULTS_DIR)
    shell:
        """
        touch {output}
        """
