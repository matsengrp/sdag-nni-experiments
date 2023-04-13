import sys
import numpy as np
import pandas as pd
from pathlib import Path
import pprint
import click

import bito
import bito.beagle_flags as beagle_flags
import bito.phylo_model_mapkeys as model_keys


pp = pprint.PrettyPrinter(indent=4)
SIMPLE_SPECIFICATION = bito.PhyloModelSpecification(
    substitution="JC69", site="constant", clock="none"
)


def parse_args():
    arg_dict = {}
    if len(sys.argv) != 6:
        print("Usage: <i:fasta> <i:seed_newick> <i:full_newick> <i:full_posterior> <i:max_iter> <i:temp_folder> <o:results>")
        exit()
    arg_dict["fasta_path"] = Path(sys.argv[1])
    arg_dict["seed_newick_path"] = Path(sys.argv[2])
    arg_dict["full_newick_path"] = Path(sys.argv[3])
    arg_dict["full_posterior_path"] = Path(sys.argv[4])
    arg_dict["max_iter"] = int(sys.argv[5])
    arg_dict["temp_folder"] = Path(sys.argv[6])
    arg_dict["results"] = Path(sys.argv[7])
    return arg_dict


def default_args():
    arg_dict = {
        "fasta_path": Path("_data/ds1.fasta"),
        "seed_newick_path": Path("_data/ds1.credible.rerooted.nwk"),
        "full_newick_path": Path("_data/ds1.mb-trees.rerooted.nwk"),
        "full_posterior": Path("_data/ds1.mb-trees.rerooted.nwk"),
        "max_iter": 1,
        "temp_folder": Path("_ignore"),
        "results": Path("results.nni-search.csv")
    }
    return arg_dict


def create_inst(fasta_path, newick_path, temp_folder):
    inst = bito.gp_instance(str(temp_folder) + "/mmap_pvs.data")
    inst.read_fasta_file(str(fasta_path))
    inst.read_newick_file(str(newick_path))
    inst.make_dag()
    return inst


def create_inst_with_all_engines(fasta_path, newick_path, temp_folder):
    inst = create_inst(fasta_path, newick_path, temp_folder)
    inst.make_gp_engine()
    inst.make_tp_engine()
    inst.make_nni_engine()
    return inst


def nni_engine_accept_top_nni(nni_engine):
    # nni_engine.run_main_loop()
    nni_engine.graft_adjacent_nnis_to_dag()
    nni_engine.filter_pre_update()
    nni_engine.filter_eval_adjacent_nnis()
    nni_engine.filter_post_update()
    nni_engine.filter_process_adjacent_nnis()
    nni_engine.remove_all_graft_nnis_from_dag()
    nni_engine.add_accepted_nnis_to_dag()

    # nni_engine.run_post_loop()
    nni_engine.update_adjacent_nnis()
    nni_engine.update_accepted_nnis()
    nni_engine.update_rejected_nnis()
    nni_engine.update_scored_nnis()


############
### MAIN ###
############

if __name__ == "__main__":
    args = parse_args()
    seed_inst = create_inst(args["fasta_path"], args["seed_newick_path"], args)
    # full_inst = create_inst(args["fasta_path"], args["full_newick_path"], args)

    seed_inst.make_gp_engine()
    seed_inst.make_tp_engine()
    seed_inst.make_nni_engine()

    nni_engine = seed_inst.get_nni_engine()
    nni_engine.run_init()

    itr = 0
    while itr < args["max_iter"]:
        print("iter:", itr, ", adj_nni_count:",
              nni_engine.adjacent_nni_count())

        itr += 1
