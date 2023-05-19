#!/bin/bash
# set -x

# declare -A ds_roots=( [1]=15 [3]=7 [4]=8 [5]=1 [6]=6 [7]=7 [8]=4 )
# result_types=( "_output" "_exp_output" )

declare -A ds_roots=( [7]=7 [8]=4 )
result_types=( "_exp_output" )

function get_branch_lengths() {
    ds=$1
    root=$2
    result_type=$3
    echo "### get-branch-lengths ds=${ds} root=${root}..."

    src_path="/home/drich/matsen-lab/sdag-nni-experiments/ds-data/ds${ds}"
    dest_path="/home/drich/matsen-lab/sdag-nni-experiments/data/ds${ds}/${result_type}"
    golden_path="/fh/fast/matsen_e/shared/vip/ds-golden/"

    pre_trprobs="${golden_path}/ds"
    post_trprobs="/_exp_output/ds"$ds".trprobs"
    fasta="${golden_path}/ds"$ds"/${result_type}/ds"$ds".fasta"
    post_pkl="${dest_path}/posterior.pkl"
    cred_nwk="${dest_path}/ds${ds}.credible.nwk"
    cred_w_branches_nwk="${dest_path}/ds${ds}.credible.with-branches.nwk"
    mb_trees="${dest_path}/ds${ds}.mb-trees.nwk"
    mb_pp="${dest_path}/ds${ds}.mb-pp.csv"

    mkdir $dest_path

    echo "# wtch-process-trprobs.py..."
    python wtch-process-trprobs.py $pre_trprobs $ds $ds $post_trprobs $root $post_pkl
    echo "# wtch-unpickle-cdf.py..."
    python wtch-unpickle-cdf.py $post_pkl $cred_nwk $mb_trees $mb_pp
    echo "# wtch-branch-optimization.py..."
    python wtch-branch-optimization.py $cred_nwk $fasta $cred_w_branches_nwk --sort=False
}

scripts_path="/home/drich/matsen-lab/sdag-nni-experiments/scripts"
cd $scripts_path

for ds in "${!ds_roots[@]}"; do
  root=${ds_roots[$ds]}
  for result_type in $result_types; do
    get_branch_lengths $ds $root $result_type
  done
done
