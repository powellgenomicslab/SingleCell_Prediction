# Run jobs

for i in 1 2 3 4; do qsub -v CELL_TYPE=$i -N pred_$i run_prediction.pbs ; done