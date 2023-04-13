import sys
from pathlib import Path
import bito
import numpy as np
import pprint
import bito.beagle_flags as beagle_flags
import bito.phylo_model_mapkeys as model_keys

SIMPLE_SPECIFICATION = bito.PhyloModelSpecification(
    substitution="JC69", site="constant", clock="none"
)

pp = pprint.PrettyPrinter(indent=4)


def print_dag_stats(inst):
    dag = inst.get_dag()
    print("dag: {} {} {}".format(dag.node_count(), dag.edge_count(), dag.taxon_count()))


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


if len(sys.argv) < 4:
    print("Usage: <i:fasta_path> <i:newick_path> <i:temp_dir>")
    exit()
fasta_path = Path(sys.argv[1])
newick_path = Path(sys.argv[2])
temp_dir = Path(sys.argv[3])
mmap_path = Path(temp_dir / "mmapped_pv.data")

gp_cutoff_threshold = -12960
tp_like_cutoff_threshold = -12960
tp_pars_cutoff_threshold = 15

if not fasta_path.exists():
    print("<i:fasta_path> does not exist.")
    exit()
if not newick_path.exists():
    print("<i:newick_path> does not exist.")
    exit()
if not temp_dir.exists():
    print("<i:temp_dir> does not exist.")
    exit()
if not temp_dir.is_dir():
    print("<i:temp_dir> is not a directory.")
    exit()

inst = bito.gp_instance(str(mmap_path))
inst.read_fasta_file(str(fasta_path))
inst.read_newick_file(str(newick_path))
inst.make_dag()
print_dag_stats(inst)
