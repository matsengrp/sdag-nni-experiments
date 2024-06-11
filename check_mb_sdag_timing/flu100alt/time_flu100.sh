#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cjennin2@fredhutch.org



root=1

cd ../../
tree_file="data/flu100/_output/short_mcmc/dsflu100.t"
input_file="check_mb_sdag_timing/flu100alt/1000000_lines.t"
output_file="check_mb_sdag_timing/flu100alt/1000000_lines.results.csv"
head -n 1000000 $tree_file > $input_file

echo $(date)

awk '$1~/tree/ {print $NF}' ${input_file} | nw_topology - | nw_reroot - ${root} \
  | nw_order - | uniq -c | sed -e "s/^[ ]*//" -e "s/[ ]/\t/" \
  > ${output_file}

echo $(date)

rm $input_file
rm $output_file