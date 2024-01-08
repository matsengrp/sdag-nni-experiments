#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err


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
    output_path=${results_path}/flu100.results.${search}.${max_iter}.time.csv

    python nni_search.py nni-search $fasta $seed_newick $credible_newick \
      $posterior_newick $pp_csv $pcsp_pp_csv $subsplit_csv --${search} --iter-max \
      ${max_iter} -o $output_path --log-time-only
}

nni_search





