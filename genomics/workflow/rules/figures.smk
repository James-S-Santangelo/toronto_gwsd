# rule multiqc_stats_notebook:
#     input:
#         rules.multiqc.output
#     output:
#         '{0}/supplemental/qc/qualimap_general_error_rate_histogram.pdf'.format(FIGURES_DIR),
#         '{0}/supplemental/qc/qualimap_aligned_vs_coverage_by_category.pdf'.format(FIGURES_DIR)
#     conda: '../envs/notebooks.yaml'
#     log: notebook = 'logs/notebooks/multiqc_stats_notebook.r.ipynb'
#     notebook:
#         "../notebooks/multiqc_stats_notebook.r.ipynb"

rule manhattan_plots:
    input:
        fst_win = rules.write_windowed_statistics.output.sfs_df,
        xpnsl_win = rules.write_windowed_statistics.output.xpnsl_df,
        xpnsl_raw = expand(rules.norm_xpnsl.output, hab_comb='Urban_Rural'),
        top_ten = rules.write_selected_regions.output.top_ten_tbl
    output:
        fst_manhat= f"{FIGURES_DIR}/manhattan/fst_allChroms.pdf",
        xpnsl_manhat= f"{FIGURES_DIR}/manhattan/xpnsl_allChroms.pdf",
        chr04_Occ_urb =f"{FIGURES_DIR}/manhattan/chr04_Occ_urb.pdf",
        chr05_Occ_urb =f"{FIGURES_DIR}/manhattan/chr05_Occ_urb.pdf",
        chr04_Occ_rur =f"{FIGURES_DIR}/manhattan/chr04_Occ_rur.pdf",
        chr05_Pall_rur =f"{FIGURES_DIR}/manhattan/chr05_Pall_rur.pdf",
        chr07_Occ_rur =f"{FIGURES_DIR}/manhattan/chr07_Occ_rur.pdf",
        chr08_Pall_rur1 =f"{FIGURES_DIR}/manhattan/chr08_Pall_rur1.pdf",
        chr08_Pall_rur2 =f"{FIGURES_DIR}/manhattan/chr08_Pall_rur2.pdf",
    conda: '../envs/figures.yaml'
    notebook:
        "../notebooks/manhattan_plots.r.ipynb"

rule population_structure_figures:
    input:
        order = expand(rules.extract_sample_angsd.output, site="4fold"),
        cov = expand(rules.pcangsd.output, site="4fold", maf="0.05"),
        evanno = rules.clumpak_best_k_by_evanno.output,
        admix_log = expand(rules.ngsadmix.output.lf, k=NGSADMIX_K, site=['4fold'], maf=['0.05'], seed=NGSADMIX_SEEDS),
        admix_qopt = expand(rules.ngsadmix.output.qopt, k=NGSADMIX_K, site=['4fold'], maf=['0.05'], seed=NGSADMIX_SEEDS),
        pi_byHab = expand(rules.angsd_diversity_neutrality_stats_byHabitat.output, site="4fold", habitat=HABITATS),
        fst_byHab = expand(rules.angsd_habitat_fst_readable.output, site=['4fold'], hab_comb=HABITAT_COMBOS),
        relate = expand(rules.ngsrelate.output, maf=['0.05'], site=['4fold'], chrom=CHROMOSOMES),
        bl = expand(rules.create_bam_list_highQualSamples.output, site="4fold"),
        pi_byPop = expand(rules.angsd_diversity_neutrality_stats_byPopulation.output, popu=POPS_MULTI_IND, site='4fold'),
        fst_byPop = expand(rules.angsd_population_fst_readable.output, pop_comb=POP_COMB_MULTI_IND, site='4fold')
    output:
        "test.txt",
        pca = f"{FIGURES_DIR}/pca_byHabitat.pdf",
        umap = f"{FIGURES_DIR}/umap_byHabitat.pdf",
        admix_optimal = f"{FIGURES_DIR}/admix_optimal.pdf",
        admix_optimal_minus = f"{FIGURES_DIR}/admix_optimal_minus.pdf",
        admix_optimal_plus = f"{FIGURES_DIR}/admix_optimal_plus.pdf",
        pi_byHab_df = f"{FIGURES_DIR}/pi_byHab.txt",
        fst_byHab_df = f"{FIGURES_DIR}/fst_byHab.txt",
        relate_byHabComb_df = f"{FIGURES_DIR}/pairwise_habitat_relatedness.txt",
        relate_bySampleComb_df = f"{FIGURES_DIR}/pairwise_sample_relatedness.txt",
        relate_byHabComb = f"{FIGURES_DIR}/relatedness_byHabitatCombination.pdf"
    conda:'../envs/figures.yaml'
    notebook:
        "../notebooks/population_structure.r.ipynb"
