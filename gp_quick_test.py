import sys
import itertools
import pandas as pd
import bito


def build_dag_inst(name, fasta_path, newick_path, mmap_path):
    inst = bito.gp_instance(mmap_path)
    inst.read_fasta_file(fasta_path)
    inst.read_newick_file(newick_path)
    inst.make_dag()
    return inst


def build_dag_with_all_adj_nnis(name, fasta_path, newick_path, mmap_path):
    inst = build_dag_inst(name, fasta_path, newick_path, mmap_path)
    dag = inst.get_dag()
    inst.make_nni_engine()
    nni_engine = inst.get_nni_engine()
    nni_engine.sync_adjacent_nnis_with_dag()
    adj_nnis = nni_engine.adjacent_nnis()
    for nni in adj_nnis:
        dag.add_node_pair(nni.get_parent(), nni.get_child())
    return dag


def build_tree_map(name, fasta_path, newick_path):
    inst_mb = bito.rooted_instance(name)
    inst_mb.read_fasta_file(fasta_path)
    inst_mb.read_newick_file(newick_path)
    trees = inst_mb.tree_collection.trees
    tree_map = {}
    tree_id = 0
    for tree in trees:
        tree_map[tree_id] = tree
    return tree_map


def build_tree_pp_map(tree_map, pp_path):
    fp = open(pp_path)
    lines = fp.readlines()
    print("# num_lines: ", len(lines), "num_trees: ", len(tree_map))
    tree_pp_map = {}
    for (tree_id, line) in zip(tree_map, lines):
        words = line.split()
        pp = float(words[0])
        tree_pp_map[tree_id] = pp
    fp.close()
    return tree_pp_map


def build_pcsp_pp_map_over_dag_pcsps(dag, tree_map, tree_pp_map):
    ''' Build a map of PCSP posteriors over all PCSP contained in DAG. '''
    pcsp_pp_map = {}
    dag_pcsps = dag.build_sorted_vector_of_edge_bitsets()
    for pcsp in dag_pcsps:
        pcsp_pp_map[pcsp] = 0
    for tree_id in tree_map:
        tree = tree_map[tree_id]
        tree_pcsps = tree.topology().build_vector_of_pcsps()
        for pcsp in tree_pcsps:
            if pcsp in pcsp_pp_map:
                pcsp_pp_map[pcsp] += tree_pp_map[pcsp]
    return pcsp_pp_map


def accumulate_pp_over_trees(dag, tree_map, tree_pp_map):
    dag_pp = 0
    dag_tree_count = 0
    for tree_id in tree_pp_map:
        tree = tree_map[tree_id]
        if dag.contains_tree(tree):
            dag_pp += tree_pp_map[tree_id]
            dag_tree_count += 1
    return dag_pp, dag_tree_count


### MAIN ###


req_args = 5
if (len(sys.argv) != (req_args + 1)):
    print("Usage: <i:fasta> <i:seed_trees> <i:mb_trees/mb_credible_trees> <i:mb_pp> <o:results>")
    exit(0)

mmap_path1 = "_ignore/mmapped_pv.data"
mmap_path2 = "_ignore/mmapped_pv.data"
mmap_path3 = "_ignore/mmapped_pv.data"
fasta_path = sys.argv[1]
seed_trees_path = sys.argv[2]
mb_trees_path = sys.argv[3]
mb_pp_path = sys.argv[4]
out_path = sys.argv[5]

print("# building seed dag...")
seed_inst = build_dag_inst("seed_dag", fasta_path, seed_trees_path, mmap_path1)
seed_dag = seed_inst.get_dag()
seed_inst.make_nni_engine()
nni_engine = seed_inst.get_nni_engine()

print("# building mb dag...")
mb_tree_map = build_tree_map("mb_trees", fasta_path, mb_trees_path)

# build map from tree in posterior to posterior density
print("# building tree posterior map...")
mb_tree_pp_map = build_tree_pp_map(mb_tree_map, mb_pp_path)

# get posterior probability of dag with each nni.
nni_dag_pp_map = {}
nni_engine.sync_adjacent_nnis_with_dag()
adj_nnis = nni_engine.adjacent_nnis()
for nni_id, nni in enumerate(adj_nnis):
    print("# finding pp for dag with nni: ", nni)
    nni_inst = build_dag_inst("nni_dag", fasta_path,
                              seed_trees_path, mmap_path3)
    nni_dag = nni_inst.get_dag()
    nni_dag.add_node_pair(nni.get_parent(), nni.get_child())
    nni_inst.make_gp_engine()
    nni_inst.estimate_branch_lengths(1e-10, 2)
    nni_inst.populate_plvs()
    nni_inst.compute_likelihoods()
    llhs = nni_inst.get_per_pcsp_log_likelihoods()
    print("# llhs: ", llhs)
    nni_dag_pp, nni_dag_tree_count = accumulate_pp_over_trees(
        nni_dag, mb_tree_map, mb_tree_pp_map)
    nni_dag_pp_map[nni_id] = (nni, nni_dag_pp, nni_dag_tree_count)
    print("# nni_dag_pp: ", nni_dag_pp,
          "nni_dag_tree_count: ", nni_dag_tree_count)

print("### RESULTS ###")
print("# nni_id nni_parent nni_child nni_pp nni_tree_count")
fp = open(out_path, "w")
pcsps = []
parent_subsplits = []
child_subsplits = []
for nni_id in nni_dag_pp_map:
    nni, pp, tree_count = nni_dag_pp_map[nni_id]
    out = f"{nni_id} {nni.get_parent().to_subsplit()} {nni.get_child().to_subsplit()} {pp} {tree_count}"
    print(out)
    fp.write(out + "\n")
fp.close()
print("# ...done")
