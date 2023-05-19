#!/bin/bash
# set -x

dss=( 1 )
result_types=( "_output" )
search_types=( "tp" "gp" "pcsp" )

function nni_search() {
    ds=$1
    result_type=$2
    search=$3
    max_iter=75
    echo "### nni_search ds=${ds} result_type=${result_type} search=${search}..."

    data_path="/home/drich/matsen-lab/sdag-nni-experiments/data/ds${ds}/"
    results_path="${data_path}/${result_type}"

    python nni_search.py nni-search ${data_path}/ds${ds}.fasta ${results_path}/ds${ds}.top1.nwk ${results_path}/ds${ds}.credible.with-branches.rerooted.nwk ${results_path}/ds${ds}.mb-pp.csv ${results_path}/ds${ds}.pcsp-pp.csv --${search} --iter-max ${max_iter} -o ${results_path}/ds${ds}.results.${search}.${max_iter}.csv
}

for ds in ${dss[@]}; do
  for result_type in ${result_types[@]}; do
    for search_type in ${search_types[@]}; do
      nni_search $ds $result_type $search_type
    done
  done
done
