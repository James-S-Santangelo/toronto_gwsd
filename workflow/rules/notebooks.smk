rule multiqc_stats_notebook:
    input:
        rules.multiqc.output
    output:
        '{0}/supplemental/qc/qualimap_general_error_rate_histogram.pdf'.format(FIGURES_DIR),
        '{0}/supplemental/qc/qualimap_aligned_vs_coverage_by_category.pdf'.format(FIGURES_DIR)
    conda: '../envs/notebooks.yaml'
    log: notebook = 'logs/notebooks/multiqc_stats_notebook.r.ipynb'
    notebook:
        "../notebooks/multiqc_stats_notebook.r.ipynb"

rule population_structure_notebook:
    input:
        rules.pop_structure_done.output
    output:
        pca_plot = '{0}/population_structure/pca_allFinalSamples.pdf'.format(FIGURES_DIR),
        admix_plot = '{0}/population_structure/admixture_allFinalSamples_bestK.pdf'.format(FIGURES_DIR)
    conda:'../envs/notebooks.yaml'
    notebook:
        "../notebooks/population_structure.r.ipynb"
