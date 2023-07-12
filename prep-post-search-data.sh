#!/bin/bash

declare -a ds_datasets=( 1 3 4 5 6 7 8 )
declare -a prior_types=( "_output" "_exp_output" )
declare -a search_types=( "tp" "gp" )

function prep_post_search_data() {
    ds=$1
    prior_type=$2
    search=$3
    max_iter=500    

    data_path="data/ds${ds}"
    results_path="${data_path}/${prior_type}"
    fasta=${data_path}/ds${ds}.fasta
    posterior_newick=${results_path}/ds${ds}.mb-trees.with-fake_branches.nwk
    pp_csv=${results_path}/ds${ds}.mb-pp.csv
    search_results=${results_path}/ds${ds}.results.${search}.${max_iter}.csv
    out_path1=${results_path}/ds${ds}.found-trees-dag-stats.${search}.credible.csv
    out_path2=${results_path}/ds${ds}.found-trees-dag-stats.${search}.posterior.csv

    python found_tree_dag_stats.py $fasta $posterior_newick $pp_csv $search_results $out_path1

    python found_tree_dag_stats.py $fasta $posterior_newick $pp_csv $search_results $out_path2 \
      --restrict_to_credible_set False
}

for ds in ${ds_datasets[@]}; do
  for prior_type in ${prior_types[@]}; do
    for search_type in ${search_types[@]}; do
      prep_post_search_data $ds $prior_type $search_type
    done
  done
done
