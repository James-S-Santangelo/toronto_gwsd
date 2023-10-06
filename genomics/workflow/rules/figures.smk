# Snakefile with rules to generate figures and some summary text files

rule compare_observed_permuted_xpnsl:
    input:
        obs = rules.write_windowed_statistics.output.xpnsl_df,
        perm = expand(rules.write_windowed_statistics_permuted.output, hab_comb="Urban_Rural", n=[x for x in range(1,1001)])
    output:
        cor_plot = f"{FIGURES_DIR}/selection/xpnsl_perm/observed_permuted_xpnsl_correlation.pdf",
        urb_mean_plot = f"{FIGURES_DIR}/selection/xpnsl_perm/urbanSel_mean.pdf",
        urb_prop_plot = f"{FIGURES_DIR}/selection/xpnsl_perm/urbanSel_prop.pdf",
        rur_mean_plot = f"{FIGURES_DIR}/selection/xpnsl_perm/ruralSel_mean.pdf",
        rur_prop_plot = f"{FIGURES_DIR}/selection/xpnsl_perm/ruralSel_prop.pdf",
        urb_perc = f"{FIGURES_DIR}/selection/xpnsl_perm/urban_percentiles.txt",
        rur_perc = f"{FIGURES_DIR}/selection/xpnsl_perm/rural_percentiles.txt",
    conda: '../envs/sweeps.yaml'
    notebook:
        "../notebooks/compare_observed_permuted_xpnsl.r.ipynb"

rule write_selected_regions:
    input:
        fst = rules.write_windowed_statistics.output.sfs_df,
        xpnsl = rules.write_windowed_statistics.output.xpnsl_df,
        urb_perc = rules.compare_observed_permuted_xpnsl.output.urb_perc,
        rur_perc = rules.compare_observed_permuted_xpnsl.output.rur_perc,
        gff = GFF_FILE 
    output:
        top_ten_genes = f'{FIGURES_DIR}/selection/top10_selected_regions_genes.txt', 
        top_ten_tbl = f'{FIGURES_DIR}/selection/top10_selected_regions_urban_rural_table.txt',
        all_xpnsl_sel = f'{FIGURES_DIR}/selection/all_selected_regions_genes.txt'
    conda: '../envs/sweeps.yaml'
    notebook:
        "../notebooks/write_top_selected_regions.r.ipynb"

rule go_enrichment_analysis:
    input:
        all_genes = rules.create_geneToGO_mapfile.output,
        all_sel = rules.write_selected_regions.output.all_xpnsl_sel,
        top_ten_genes = rules.write_selected_regions.output.top_ten_genes
    output:
        all_go_res = f'{FIGURES_DIR}/selection/all_go_results.txt'
    conda: '../envs/sweeps.yaml'
    notebook:
        "../notebooks/go_enrichment_analysis.r.ipynb"

rule manhattan_plots:
    input:
        fst_win = rules.write_windowed_statistics.output.sfs_df,
        xpnsl_win = rules.write_windowed_statistics.output.xpnsl_df,
        xpnsl_raw = expand(rules.norm_xpnsl.output, hab_comb='Urban_Rural'),
        top_ten = rules.write_selected_regions.output.top_ten_tbl
    output:
        fst_manhat= f"{FIGURES_DIR}/selection/manhattan/fst_allChroms.pdf",
        xpnsl_manhat= f"{FIGURES_DIR}/selection/manhattan/xpnsl_allChroms.pdf",
        chr04_Occ_urb =f"{FIGURES_DIR}/selection/manhattan/chr04_Occ_urb.pdf",
        chr05_Occ_urb =f"{FIGURES_DIR}/selection/manhattan/chr05_Occ_urb.pdf",
        chr04_Occ_rur =f"{FIGURES_DIR}/selection/manhattan/chr04_Occ_rur.pdf",
        chr05_Pall_rur =f"{FIGURES_DIR}/selection/manhattan/chr05_Pall_rur.pdf",
        chr07_Occ_rur =f"{FIGURES_DIR}/selection/manhattan/chr07_Occ_rur.pdf",
        chr08_Pall_rur1 =f"{FIGURES_DIR}/selection/manhattan/chr08_Pall_rur1.pdf",
        chr08_Pall_rur2 =f"{FIGURES_DIR}/selection/manhattan/chr08_Pall_rur2.pdf",
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
        pca = f"{FIGURES_DIR}/pop_struct/pca_byHabitat.pdf",
        umap = f"{FIGURES_DIR}/pop_struct/umap_byHabitat.pdf",
        admix_optimal = f"{FIGURES_DIR}/pop_struct/admix_optimal.pdf",
        admix_optimal_minus = f"{FIGURES_DIR}/pop_struct/admix_optimal_minus.pdf",
        admix_optimal_plus = f"{FIGURES_DIR}/pop_struct/admix_optimal_plus.pdf",
        pi_byHab_df = f"{FIGURES_DIR}/pop_struct/pi_byHab.txt",
        fst_byHab_df = f"{FIGURES_DIR}/pop_struct/fst_byHab.txt",
        relate_byHabComb_df = f"{FIGURES_DIR}/pop_struct/pairwise_habitat_relatedness.txt",
        relate_bySampleComb_df = f"{FIGURES_DIR}/pop_struct/pairwise_sample_relatedness.txt",
        relate_byHabComb = f"{FIGURES_DIR}/pop_struct/relatedness_byHabitatCombination.pdf"
    conda:'../envs/figures.yaml'
    notebook:
        "../notebooks/population_structure.r.ipynb"

rule figures_done:
    input:
        rules.population_structure_figures.output,
        rules.manhattan_plots.output,
        rules.go_enrichment_analysis.output,
        rules.compare_observed_permuted_xpnsl.output
    output:
        f"{FIGURES_DIR}/figures.done"
    shell:
        """
        touch {output}
        """
