# rsync -vuar -P --delete \
    # santang3@gra-dtn1.computecanada.ca:/home/santang3/scratch/gwsd/logs/ \
    # /scratch/projects/trifolium/gwsd/logs

# rsync -vuar -P --delete \
    # santang3@gra-dtn1.computecanada.ca:/home/santang3/scratch/gwsd/slurm_logs/ \
    # /scratch/projects/trifolium/gwsd/slurm_logs

rsync -vuar -P \
    santang3@gra-dtn1.computecanada.ca:/home/santang3/scratch/gwsd/results/ \
    /scratch/projects/trifolium/gwsd/results
