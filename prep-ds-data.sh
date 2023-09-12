#!/bin/bash


declare -A ds_roots=( [1]=15 [3]=7 [4]=8 [5]=1 [6]=6 [7]=7 [8]=4 )
declare -a prior_types=( "_output" "_exp_output" )

function prep_ds_data() {
    ds=$1
    root=$2
    prior_type=$3
    echo "### Preparing data files for ds${ds} in ${prior_type}..."

    dest_path="data/ds${ds}/${prior_type}"
    golden_path="/fh/fast/matsen_e/shared/vip/ds-golden/"

    pre_trprobs="${golden_path}/ds"
    post_trprobs="/${prior_type}/ds"$ds".trprobs"
    fasta="${golden_path}/ds"$ds"/${prior_type}/ds"$ds".fasta"
    post_pkl="${dest_path}/posterior.pkl"
    
    cred_nwk="${dest_path}/ds${ds}.credible.nwk"
    
    top_topology_nwk="${dest_path}/ds${ds}.top1.topology.nwk"
    top_tree_nwk="${dest_path}/ds${ds}.top1.tree.nwk" 
    first_nwk="${dest_path}/ds${ds}.top1.nwk" 
    
    cred_with_fake_branches_nwk="${dest_path}/ds${ds}.credible.with-fake-branches.nwk" 
    mb_trees_with_fake_branches_nwk="${dest_path}/ds${ds}.mb-trees.with-fake_branches.nwk"    
    
    mb_trees="${dest_path}/ds${ds}.mb-trees.nwk"
    mb_pp="${dest_path}/ds${ds}.mb-pp.csv"
    pcsp_pp="${dest_path}/ds${ds}.pcsp-pp.csv"
    subsplit_path="${dest_path}/ds${ds}.subsplit.csv"

    mkdir -p $dest_path

    echo "# Processing trprobs file via wtch-process-trprobs.py..."
    wtch-process-trprobs.py $pre_trprobs $ds $ds $post_trprobs $root $post_pkl

    echo "# Unpickling cumulative density via wtch-unpickle-cdf.py..."
    wtch-unpickle-cdf.py $post_pkl $cred_nwk $mb_trees $mb_pp

    echo "# Obtaining optimal branch lengths for top tree via wtch-branch-optimization.py..."
    head -n 1 $cred_nwk > $top_topology_nwk
    wtch-branch-optimization.py $top_topology_nwk $fasta $top_tree_nwk --sort=False

    echo "# Rerooting optimized top tree via nw_reroot..."
    nw_reroot $top_tree_nwk $root > $first_nwk

    echo "# Adding dummy branch lengths to Mr. Bayes credible set and posterior trees..."
    python add_dummy_branch_lengths.py $mb_trees $mb_trees_with_fake_branches_nwk
    python add_dummy_branch_lengths.py $cred_nwk $cred_with_fake_branches_nwk

    echo "# Calculating pcsp posterior weights via nni_search.py build-pcsp-map..."
    python nni_search.py build-pcsp-map $fasta $mb_trees $mb_pp $first_nwk -o $pcsp_pp

    echo "# Calculating credible and posterior subsplits via nni_search.py build-subsplit-map..."
    python nni_search.py build-subsplit-map $fasta $mb_trees $mb_pp $first_nwk -o $subsplit_path

    echo "# Copying fasta file to local directory..."
    fasta="${golden_path}/ds"$ds"/${prior_type}/ds"$ds".fasta"
    cp $fasta data/ds${ds}/ds${ds}.fasta
    
    echo "# Finished processing files for ds${ds} in $prior_type."
}

#scripts_path="/home/drich/matsen-lab/sdag-nni-experiments/scripts"
#cd $scripts_path

for ds in "${!ds_roots[@]}"; do
  root=${ds_roots[$ds]}
  for prior_type in ${prior_types[@]}; do
    prep_ds_data $ds $root $prior_type
  done
done

