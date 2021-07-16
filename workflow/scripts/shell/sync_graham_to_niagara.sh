rsync -vuar -P \
    --prune-empty-dirs \
    --include '*/' \
    --include 'program_resources/**' \
    --exclude '*' \
    ../../../results/ \
    santang3@nia-datamover1.scinet.utoronto.ca:/scratch/n/nessrobe/santang3/toronto_gwsd/results
