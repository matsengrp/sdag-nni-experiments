#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err

ds=5
result_type="_output"
search_type="tp"

function nni_search() {
    ds=$1
    result_type=$2
    search=$3
    max_iter=2000
    echo "### nni_search ds=${ds} result_type=${result_type} search=${search}..."

    data_path="data/ds${ds}"
    results_path="${data_path}/${result_type}"

    fasta=${data_path}/ds${ds}.fasta
    seed_newick=${results_path}/ds${ds}.top1.nwk

    credible_newick=${results_path}/ds${ds}.credible.with-fake-branches.nwk
    posterior_newick=${results_path}/ds${ds}.mb-trees.with-fake_branches.nwk
    pp_csv=${results_path}/ds${ds}.mb-pp.csv
    pcsp_pp_csv=${results_path}/ds${ds}.pcsp-pp.csv
    subsplit_csv=${results_path}/ds${ds}.subsplit.csv
    output_path=${results_path}/ds${ds}.results.${search}.${max_iter}.csv

    python nni_search.py nni-search $fasta $seed_newick $credible_newick \
      $posterior_newick $pp_csv $pcsp_pp_csv $subsplit_csv --${search} --iter-max \
      ${max_iter} -o $output_path
}

nni_search $ds $result_type $search_type



