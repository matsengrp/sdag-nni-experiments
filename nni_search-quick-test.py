import sys
from pathlib import Path
import bito
import numpy as np
import pprint
import bito.beagle_flags as beagle_flags
import bito.phylo_model_mapkeys as model_keys
import pandas as pd

SIMPLE_SPECIFICATION = bito.PhyloModelSpecification(
    substitution="JC69", site="constant", clock="none"
)

pp = pprint.PrettyPrinter(indent=4)


def print_dag_stats(inst):
    dag = inst.get_dag()
    print("dag: {} {} {}".format(dag.node_count(),
          dag.edge_count(), dag.taxon_count()))


def print_graft_dag_stats(inst):
    graft_dag = inst.get_nni_engine().get_graft_dag()
    print(
        "graft_dag: {} {}".format(
            graft_dag.graft_node_count(), graft_dag.graft_edge_count()
        )
    )


def print_nni_engine_stats(inst):
    nni_engine = inst.get_nni_engine()
    dag = inst.get_dag()
    print(
        "iter: {}, adj_count: {}, accepted_count: {} {}, rejected_count: {} {}, topology_count: {}".format(
            nni_engine.iter_count(),
            nni_engine.adjacent_nni_count(),
            nni_engine.accepted_nni_count(),
            nni_engine.past_accepted_nni_count(),
            nni_engine.rejected_nni_count(),
            nni_engine.past_rejected_nni_count(),
            dag.topology_count(),
        )
    )


if len(sys.argv) < 6:
    print(
        "Usage: <i:fasta_path> <i:newick_path> <i:temp_dir> <i:test_newick_path> <i:test_posteriors>"
    )
    exit()
fasta_path = Path(sys.argv[1])
seed_newick_path = Path(sys.argv[2])
temp_dir = Path(sys.argv[3])
mmap_path = Path(temp_dir / "mmapped_pv.test.data")
truth_newick_path = Path(sys.argv[4])
truth_posterior_path = Path(sys.argv[5])
truth_mmap_path = Path(temp_dir / "mmapped_pv.truth.data")

gp_cutoff_threshold = -20000
tp_like_cutoff_threshold = -12960
tp_pars_cutoff_threshold = 0

if not fasta_path.exists():
    print("<i:fasta_path> does not exist.")
    exit()
if not seed_newick_path.exists():
    print("<i:newick_path> does not exist.")
    exit()
if not temp_dir.exists():
    print("<i:temp_dir> does not exist.")
    exit()
if not temp_dir.is_dir():
    print("<i:temp_dir> is not a directory.")
    exit()
if not truth_newick_path.exists():
    print("<i:test_newick_path> does not exist.")
    exit()


def create_gp_inst(fasta_path_, newick_path_, mmap_path_):
    # build gp_instance from paths
    inst = bito.gp_instance(str(mmap_path_))
    inst.read_fasta_file(str(fasta_path_))
    inst.read_newick_file(str(newick_path_))
    return inst


def init_gp_inst(inst):
    inst.make_dag()
    inst.make_gp_engine()
    inst.make_tp_engine()
    inst.make_nni_engine()
    return

### MAIN ###


inst = create_gp_inst(fasta_path, seed_newick_path, mmap_path)
init_gp_inst(inst)
dag = inst.get_dag()

gp_engine = inst.get_gp_engine()
inst.take_first_branch_length()

tp_engine = inst.get_tp_engine()
inst.tp_engine_set_branch_lengths_by_taking_first()

loaded_trees = inst.currently_loaded_trees_with_gp_branch_lengths()
loaded_tree = loaded_trees.trees[0]
top_tree = tp_engine.get_top_tree_with_edge(bito.edge_id(1))

print("loaded_tree:", dag.tree_to_newick_tree(loaded_tree))
print("top_tree:", dag.tree_to_newick_tree(top_tree))
