#!/bin/bash
#SBATCH --job-name=nni_search_single
#SBATCH --ntasks=1
#SBATCH --time=08:00:00
#SBATCH --output=nni_search_single.%j.%t.log


dss=( 1 3 4 5 6 7 8 )
result_types=( "_output" "_exp_output" )
search_types=( "tp" "gp" "pcsp" )

function nni_search() {
    ds=$1
    result_type=$2
    search=$3
    max_iter=$4
    echo "### nni_search ds=${ds} result_type=${result_type} search=${search}, max_iter=${max_iter}..."

    data_path="/home/drich/matsen-lab/sdag-nni-experiments/data/ds${ds}/"
    results_path="${data_path}/${result_type}"

    python nni_search.py nni-search ${data_path}/ds${ds}.fasta ${results_path}/ds${ds}.top1.nwk ${results_path}/ds${ds}.credible.with-branches.rerooted.nwk ${results_path}/ds${ds}.mb-pp.csv ${results_path}/ds${ds}.pcsp-pp.csv --${search} --iter-max ${max_iter} -o ${results_path}/ds${ds}.results.${search}.${max_iter}.csv
}

if [ "$#" != "4" ]; then
    echo "ERROR: Incorrect number of args."
    echo "Usage: <ds> <result_type> <search_type> <max_iter>"
    exit
fi

ds=$1
result_type=$2
search_type=$3
max_iter=$4
nni_search $ds $result_type $search_type $max_iter
