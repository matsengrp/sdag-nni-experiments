#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cjennin2@fredhutch.org


declare -A ds_roots=( [3]=7 )
result_types=( "_output" )
search_types=( "tp" )

function nni_search() {
    ds=$1
    root=$2
    result_type=$3
    search=$4
    max_iter=300
    echo "### nni_search ds=${ds} result_type=${result_type} search=${search}..."

    data_path="data/ds${ds}"
    results_path="${data_path}/${result_type}"
    fasta=${data_path}/ds${ds}.fasta
    cred_nwk="${results_path}/ds${ds}.credible.nwk"
    seed_newick=${results_path}/ds${ds}.raxml_start.nwk
    top_topology_nwk="multiple_trees_data/raxml_tests/ds${ds}/ds${ds}.raxml.unique.nwk"
    top_tree_nwk="${results_path}/ds${ds}.top_raxml.tree.nwk" 


    echo "# Obtaining optimal branch lengths for top trees via wtch-branch-optimization.py..."
    python wtch-branch-optimization.py $top_topology_nwk $fasta $top_tree_nwk --sort=True
    nw_reroot $top_tree_nwk $root > $seed_newick 

    credible_newick=${results_path}/ds${ds}.credible.with-fake-branches.nwk
    posterior_newick=${results_path}/ds${ds}.mb-trees.with-fake_branches.nwk
    pp_csv=${results_path}/ds${ds}.mb-pp.csv
    pcsp_pp_csv=${results_path}/ds${ds}.pcsp-pp.csv
    subsplit_csv=${results_path}/ds${ds}.subsplit.csv
    output_path=${results_path}/ds${ds}.results.${search}.${max_iter}.top_raxml.csv


    python nni_search.py nni-search $fasta $seed_newick $credible_newick $posterior_newick $pp_csv \
      $pcsp_pp_csv $subsplit_csv --${search} --iter-max ${max_iter} -o $output_path
}


for ds in "${!ds_roots[@]}"; do
  root=${ds_roots[$ds]}
  for result_type in ${result_types[@]}; do
    for search_type in ${search_types[@]}; do
	    nni_search $ds $root $result_type $search_type
    done
  done
done