#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cjennin2@fredhutch.org


### This script runs top pruning on the ds-datadsets using the top N trees of the 
### posterior for the initial trees, with N=1,2,..,10.



function nni_search() {
    root=1
    result_type="_output"
    search="tp"
    start_tree_count=10 
    max_iter=300

    data_path="data/flu100"
    results_path="${data_path}/${result_type}"
    fasta=${data_path}/flu100.fasta
    cred_nwk="${results_path}/flu100.credible.nwk"
    seed_newick=${results_path}/flu100.multiple_start${start_tree_count}.nwk
    top_topology_nwk="${results_path}/flu100.top${start_tree_count}.topology.nwk"
    top_tree_nwk="${results_path}/flu100.top${start_tree_count}.tree.nwk" 

    echo "# Obtaining optimal branch lengths for top trees via wtch-branch-optimization.py..."
    head -n $start_tree_count $cred_nwk > $top_topology_nwk
    python wtch-branch-optimization.py $top_topology_nwk $fasta $top_tree_nwk --sort=True
    nw_reroot $top_tree_nwk $root > $seed_newick 

    credible_newick=${results_path}/flu100.credible.with-fake-branches.nwk
    posterior_newick=${results_path}/flu100.mb-trees.with-fake_branches.nwk
    pp_csv=${results_path}/flu100.mb-pp.csv
    pcsp_pp_csv=${results_path}/flu100.pcsp-pp.csv
    subsplit_csv=${results_path}/flu100.subsplit.csv
    output_path=${results_path}/flu100.results.${search}.${max_iter}.top${start_tree_count}.csv

    python nni_search.py nni-search $fasta $seed_newick $credible_newick $posterior_newick $pp_csv \
      $pcsp_pp_csv $subsplit_csv --${search} --iter-max ${max_iter} -o $output_path
}


nni_search

