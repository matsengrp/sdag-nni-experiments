#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cjennin2@fredhutch.org


cd ../../
tree_file="data/flu100/_output/short_mcmc/dsflu100.t"
for i in 1000 10000 100000 1000000
do
    input_file="check_mb_sdag_timing/flu100/${i}_lines.t"
    output_file="check_mb_sdag_timing/flu100/${i}_lines.results.csv"
    head -n $i $tree_file > $input_file
    python python_processing_of_mb_file.py ${input_file} 1 $output_file
    rm $input_file
    rm $output_file
one
