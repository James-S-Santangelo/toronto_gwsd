rsync -vuar -P \
    santang3@nia-datamover1.scinet.utoronto.ca:/scratch/n/nessrobe/santang3/toronto_gwsd/workflow/logs/ \
    ../../logs && 

rsync -vuar -P \
    santang3@nia-datamover1.scinet.utoronto.ca:/scratch/n/nessrobe/santang3/toronto_gwsd/workflow/slurm_logs/ \
    ../../slurm_logs && 

rsync -vuar -P \
    santang3@nia-datamover1.scinet.utoronto.ca:/scratch/n/nessrobe/santang3/toronto_gwsd/results/angsd/ \
    ../../../results/angsd
