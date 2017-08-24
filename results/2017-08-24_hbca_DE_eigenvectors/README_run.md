# Run jobs

for i in DC1 DC2 DC3 DC4 DC5 DC6 Mono1 Mono2 Mono3 Mono4; do qsub -v CELL_TYPE=$i -N de_$i run_prediction.pbs ; done
