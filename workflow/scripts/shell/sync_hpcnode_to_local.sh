# ANGSD
rsync -vuar -P \
    --prune-empty-dirs \
    --include='*/' \
    --include='*allChroms.sfs' \
    --include='*01.1_sum*' \
    --include='*GL*_baq*.sfs' \
    --include='allSamples*.pestPG' \
    --exclude='*' \
    santang3@hpcnode1.utm.utoronto.ca:/scratch/research/projects/trifolium/toronto_gwsd/results/angsd/ \
    ../../../data

# NGSRELATE
rsync -vuar -P \
    --prune-empty-dirs \
    --include='*/' \
    --include='*ngsRelate*' \
    --exclude='*' \
    santang3@hpcnode1.utm.utoronto.ca:/scratch/research/projects/trifolium/toronto_gwsd/results/ngsrelate \
    ../../../data

# POPULATION STRUCTURE
rsync -vuar -P \
    --prune-empty-dirs \
    --include='*/' \
    --include='*.cov' \
    --include='*seed1.qopt' \
    --exclude='*' \
    santang3@hpcnode1.utm.utoronto.ca:/scratch/research/projects/trifolium/toronto_gwsd/results/population_structure/ \
    ../../../data

# PROGRAM RESOURCES
rsync -vuar -P \
    --prune-empty-dirs \
    --include='*/' \
    --include='angsd_sample_order.txt' \
    --exclude='*' \
    santang3@hpcnode1.utm.utoronto.ca:/scratch/research/projects/trifolium/toronto_gwsd/results/program_resources/ \
    ../../../data
