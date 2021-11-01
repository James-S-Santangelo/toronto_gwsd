rsync -vuar -P \
    santang3@gra-dtn1.computecanada.ca:/home/santang3/scratch/gwsd/toronto_gwsd/workflow/logs/ \
    ../../logs

rsync -vuar -P \
    santang3@gra-dtn1.computecanada.ca:/home/santang3/scratch/gwsd/toronto_gwsd/workflow/slurm_logs/ \
    ../../../slurm_logs

rsync -vuar -P \
    santang3@gra-dtn1.computecanada.ca:/home/santang3/scratch/gwsd/results/ \
    ../../../results
