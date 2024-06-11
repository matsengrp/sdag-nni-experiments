#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err



function prep_flu_data() {
    root="1"
    echo "### Preparing data files for flu data set."

    dest_path="data/flu100/_output"
    fasta="data/flu100/flu100.fasta"
    post_pkl="${dest_path}/posterior.pkl"
    cred_nwk="${dest_path}/flu100.credible.nwk"
    top_topology_nwk="${dest_path}/flu100.top1.topology.nwk"
    top_tree_nwk="${dest_path}/flu100.top1.tree.nwk" 
    first_nwk="${dest_path}/flu100.top1.nwk" 
    cred_with_fake_branches_nwk="${dest_path}/flu100.credible.with-fake-branches.nwk" 
    mb_trees_with_fake_branches_nwk="${dest_path}/flu100.mb-trees.with-fake_branches.nwk"    
    mb_trees="${dest_path}/flu100.mb-trees.nwk"
    mb_pp="${dest_path}/flu100.mb-pp.csv"
    pcsp_pp="${dest_path}/flu100.pcsp-pp.csv"
    subsplit_path="${dest_path}/flu100.subsplit.csv"

    mkdir -p $dest_path

    echo "# Processing trprobs file via wtch-process-trprobs.py..."
    wtch-process-trprobs.py "data/flu10" "0" "0" "flu100.trprobs" $root $post_pkl

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
    python nni_search.py build-pcsp-map $fasta $mb_trees $mb_pp -o $pcsp_pp

    echo "# Calculating credible and posterior subsplits via nni_search.py build-subsplit-map..."
    python nni_search.py build-subsplit-map $fasta $mb_trees $mb_pp -o $subsplit_path
}


function prep_mb_data() {
    root=1
    prior_str="unconstrained:uniform(0,1)"

    echo "### Preparing Mr. Bayes data files."

    fasta_path="data/flu100/flu100.fasta"
    data_path="data/flu100/_output"
    posterior_pickle_path="${data_path}/posterior.pkl" 
    first_tree_path="${data_path}/flu100.top1.nwk" 
    posterior_newick_path="${data_path}/flu100.mb-trees.with-fake_branches.nwk"
    pp_csv="${data_path}/flu100.mb-pp.csv"
    pcsp_pp_csv="${data_path}/flu100.pcsp-pp.csv"
    subsplit_path="${data_path}/flu100.subsplit.csv"
    mb_run_dir="${data_path}/short_mcmc"
    mb_run_path="${mb_run_dir}/run.mb"
    rerooted_topologies_path="${mb_run_dir}/rerooted-topology-sequence.tab"
    out_path="${mb_run_dir}/mcmc_search_sdag_stats.csv"    

    mkdir -p $mb_run_dir

    # Create the Mr. Bayes file.
    python make_mb_file.py "flu100" $first_tree_path $prior_str $mb_run_path \
	    --generations 1000000


    echo "execute ../../flu100.n.nex;" > "${mb_run_path}.temp"
    tail +2 $mb_run_path >> "${mb_run_path}.temp"
    rm $mb_run_path
    mv "${mb_run_path}.temp" $mb_run_path

    echo "### Running Mr. Bayes..."
    cd $mb_run_dir
    mb run.mb | tee mb.log

    echo "### Processing Mr. Bayes trees file..."
    awk '$1~/tree/ {print $NF}' dsflu100.t | nw_topology - | nw_reroot - ${root} \
      | nw_order - | uniq -c | sed -e "s/^[ ]*//" -e "s/[ ]/\t/" \
      > rerooted-topology-sequence.tab
    cd ../../../..

    echo "### Calculating statistics..."
    python mb_comparison_stats.py $posterior_pickle_path $rerooted_topologies_path \
      $fasta_path $first_tree_path $posterior_newick_path $pp_csv $pcsp_pp_csv \
      $subsplit_path $out_path

    echo "### Results written to mb_run_dir="${data_path}/short_mcmc""

}



function nni_search() {
    search="tp"
    max_iter=2000
    echo "### nni_search..."

    data_path="data/flu100"
    results_path="${data_path}/_output"

    fasta=${data_path}/flu100.fasta
    seed_newick=${results_path}/flu100.top1.nwk

    credible_newick=${results_path}/flu100.credible.with-fake-branches.nwk
    posterior_newick=${results_path}/flu100.mb-trees.with-fake_branches.nwk
    pp_csv=${results_path}/flu100.mb-pp.csv
    pcsp_pp_csv=${results_path}/flu100.pcsp-pp.csv
    subsplit_csv=${results_path}/flu100.subsplit.csv
    output_path=${results_path}/flu100.results.${search}.${max_iter}.csv

    python nni_search.py nni-search $fasta $seed_newick $credible_newick \
      $posterior_newick $pp_csv $pcsp_pp_csv $subsplit_csv --${search} --iter-max \
      ${max_iter} -o $output_path
}

prep_flu_data
prep_mb_data
nni_search





