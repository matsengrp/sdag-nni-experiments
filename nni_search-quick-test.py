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


def load_test_inst():
    inst = create_gp_inst(fasta_path, seed_newick_path, mmap_path)
    init_gp_inst(inst)
    return inst


def load_truth_inst():
    truth_inst = create_gp_inst(fasta_path, truth_newick_path, truth_mmap_path)
    init_gp_inst(truth_inst)
    return truth_inst


def load_truth_sbn_inst():
    truth_sbn_inst = bito.rooted_instance("truth")
    truth_sbn_inst.read_fasta_file(str(fasta_path))
    truth_sbn_inst.read_newick_file(str(truth_newick_path))
    return truth_sbn_inst


def load_truth_posterior():
    posteriors = []
    with truth_posterior_path.open("r") as fp:
        for line in fp:
            posteriors.append(float(line))
    return posteriors


def run_nni_search(inst, truth_sbn_inst, truth_inst):
    print("=== RUN SEARCH ===")
    dag = inst.get_dag()
    nni_engine = inst.get_nni_engine()
    gp_engine = inst.get_gp_engine()
    tp_engine = inst.get_tp_engine()
    graft_dag = nni_engine.get_graft_dag()
    nni_engine.run_init()

    print("INITIAL:")
    print_dag_stats(inst)
    print_graft_dag_stats(inst)
    print_nni_engine_stats(inst)

    max_iter = 1
    itr = 0
    while itr < max_iter:
        print("iter:", itr, ", adj_nni_count:",
              nni_engine.adjacent_nni_count())
        # print(nni_engine.adjacent_nnis())
        is_adj_valid = []
        for adj_nni in nni_engine.adjacent_nnis():
            is_adj_valid.append(
                graft_dag.get_host_dag().is_valid_add_node_pair(
                    adj_nni.parent(), adj_nni.child()
                )
            )
        print("is_adj_valid:", is_adj_valid)
        if nni_engine.adjacent_nni_count() <= 0:
            break

        print("== MAIN ==")
        nni_engine.run_main_loop()
        print_nni_engine_stats(inst)
        print_dag_stats(inst)
        print_graft_dag_stats(inst)

        print("NNIS:")
        pp.pprint(nni_engine.scored_nnis())

        print("== POST ==")
        nni_engine.run_post_loop()

        print("== FIND TREES ==")
        find_trees_in_dag(
            inst.get_dag(), truth_inst.currently_loaded_trees_with_gp_branch_lengths()
        )

        print("== RESYNC ==")
        nni_engine.sync_adjacent_nnis_with_dag()
        print_nni_engine_stats(inst)
        itr += 1

    print("== FINAL ==")
    print_dag_stats(inst)
    print_graft_dag_stats(inst)
    print_nni_engine_stats(inst)


def run_test_results(inst, truth_sbn_inst, truth_inst):
    # Check for trees in DAG
    dag = inst.get_dag()
    print("== FIND TREES ==")
    posteriors = load_truth_posterior()
    contains_trees_in_dag = []
    tree_count = 0
    tree_found = 0
    posterior_found = 0.0
    for tree in truth_sbn_inst.tree_collection.trees:
        # print("Tree:", tree.to_newick())
        contains_tree_in_dag = dag.contains_tree(tree, True)
        contains_trees_in_dag.append((tree, contains_tree_in_dag))
        # print("Contains Tree:", contains_tree_in_dag)
        if contains_tree_in_dag:
            tree_found += 1
            # posterior_found += posteriors[tree_count]
        tree_count += 1

    print(
        "Percent of Support Trees Found:",
        tree_found,
        tree_count,
        float(tree_found) / float(tree_count),
    )
    print("contains_tree_in_dag:", contains_trees_in_dag)

    # Check for unused nodes and edges in DAG
    truth_inst.make_dag()
    truth_dag = truth_inst.get_dag()

    # Check that taxons mapping matches
    taxon = dag.get_taxon_map()
    truth_taxon = truth_dag.get_taxon_map()
    print("Truth Taxons:", truth_taxon)
    print("Test Taxons:", taxon)
    for key in taxon.keys():
        print("Compare Taxons:", key, str(taxon[key]) == str(truth_taxon[key]))

    def node_vector_to_str(my_dag):
        set_of_nodes = []
        for node in my_dag.build_sorted_vector_of_node_bitsets():
            set_of_nodes.append(node.to_string())
        return set(set_of_nodes)

    nodes = node_vector_to_str(dag)
    truth_nodes = node_vector_to_str(truth_dag)
    extra_nodes = nodes.difference(truth_nodes)
    missed_nodes = truth_nodes.difference(nodes)
    print("node counts:", dag.node_count(), truth_dag.node_count())
    print("extra_nodes:", len(extra_nodes))
    print("missed_nodes:", len(missed_nodes))

    def edge_vector_to_str(my_dag):
        set_of_edges = []
        for edge in my_dag.build_sorted_vector_of_edge_bitsets():
            set_of_edges.append(edge.to_string())
        return set(set_of_edges)

    edges = edge_vector_to_str(dag)
    truth_edges = edge_vector_to_str(truth_dag)
    extra_edges = edges.difference(truth_edges)
    missed_edges = truth_edges.difference(edges)
    print("edge counts:", dag.edge_count(), truth_dag.edge_count())
    print("extra_edges:", len(extra_edges))
    print("missed_edges:", len(missed_edges))


def calc_top_trees(fasta_path_, newick_path_, mmap_path_):
    print("=== build engines ===")
    inst = create_gp_inst(fasta_path_, newick_path_, mmap_path_)
    # inst.estimate_branch_lengths(1e-4, 100)
    inst.take_first_branch_length()
    dag = inst.get_dag()
    tp_engine = inst.get_tp_engine()
    # like_tree_engine = inst.get_likelihood_tree_engine()
    # pars_tree_engine = inst.get_parsimony_tree_engine()
    print("=== tree scores [BEGIN] ===")
    print("edge_count:", tp_engine.edge_count())
    for i in range(tp_engine.edge_count()):
        edge_id = bito.edge_id(i)
        tree = tp_engine.get_top_tree_with_edge(edge_id)
        print("branch_lengths:", np.array(tree.branch_lengths))
        likelihood_score = 0
        # likelihood_score = inst.compute_tree_likelihood(tree)
        parsimony_score = inst.compute_tree_parsimony(tree)
        print(
            "edge_id:",
            edge_id,
            "likelihood_score:",
            likelihood_score,
            "parsimony_score:",
            parsimony_score,
        )
    print("=== tree scores [END] ===")


def calc_ranked_list(fasta_path_, newick_path_, mmap_path_):
    print("=== build engines ===")
    inst = create_gp_inst(fasta_path_, newick_path_, mmap_path_)
    # inst.estimate_branch_lengths(1e-4, 100)
    inst.take_first_branch_length()
    tree_collection = inst.currently_loaded_trees_with_gp_branch_lengths()
    print("branch_lengths:", np.array(inst.get_branch_lengths))
    tree_counter = 0
    for tree in tree_collection.trees:
        # print("tree:", tree_counter, tree.to_newick())
        # print("branch_lengths:", np.array(tree.branch_lengths))
        likelihood_score = inst.compute_tree_likelihood(tree)
        parsimony_score = inst.compute_tree_parsimony(tree)
        print(
            "tree_counter:",
            tree_counter,
            "likelihood_score:",
            likelihood_score,
            "parsimony_score:",
            parsimony_score,
        )
        tree_counter += 1


def find_trees_in_dag(dag, tree_collection):
    tree_counter = 0
    trees_found = 0
    for tree in tree_collection.trees:
        contains_tree = dag.contains_tree(tree, True)
        if contains_tree:
            trees_found += 1
        # print("tree_counter:", tree_counter, "contains_tree:", contains_tree)
        tree_counter += 1
    print("trees_found:", trees_found, "of", tree_counter)
    return trees_found / tree_counter


def run_nni_expansion_until_all_trees_found(fasta_path_, newick_path_, mmap_path_, final_newick_path_, final_mmap_path_):
    print("=== RUN_NNI_EXPANSION ===")
    print("=== load ===")
    inst = create_gp_inst(fasta_path_, newick_path_, mmap_path_)
    init_gp_inst(inst)
    nni_engine = inst.get_nni_engine()
    nni_engine.set_no_filter(True)
    tp_engine = inst.get_tp_engine()
    dag = inst.get_dag()
    node_bitsets = dag.build_sorted_vector_of_node_bitsets()
    edge_bitsets = dag.build_sorted_vector_of_edge_bitsets()

    truth_inst = create_gp_inst(
        fasta_path_, final_newick_path_, final_mmap_path_)
    init_gp_inst(truth_inst)
    truth_dag = truth_inst.get_dag()
    truth_node_bitsets = truth_dag.build_sorted_vector_of_node_bitsets()
    truth_edge_bitsets = truth_dag.build_sorted_vector_of_edge_bitsets()

    truth_sbn_inst = bito.rooted_instance("truth")
    truth_sbn_inst.read_fasta_file(str(fasta_path_))
    truth_sbn_inst.read_newick_file(str(final_newick_path_))

    print("node_bitsets:", dag.node_count(), truth_dag.node_count())
    cover_node_count = 0
    for bitset in truth_node_bitsets:
        contains = dag.contains_node(bitset)
        if contains:
            cover_node_count += 1
    print("nodes_covered:", cover_node_count, "of", truth_dag.node_count())

    print("edge_bitsets:", dag.edge_count(), truth_dag.edge_count())
    cover_edge_count = 0
    for bitset in truth_edge_bitsets:
        contains = dag.contains_edge(bitset)
        if contains:
            cover_edge_count += 1
    print("edges_covered:", cover_edge_count, "of", truth_dag.edge_count())

    iterations = 0
    found_all_trees = False
    trees_added = []

    print("=== begin search ===")
    while not found_all_trees:
        iterations += 1
        print("dag:", dag.node_count(), dag.edge_count())
        print("iterations:", iterations)
        nni_engine.sync_adjacent_nnis_with_dag()
        print("adj_nnis:", nni_engine.adjacent_nni_count())
        nni_engine.filter_eval_adjacent_nnis()
        accepted_nni_cnt = 0
        accepted_nnis = {}
        adj_nnis = nni_engine.adjacent_nnis()

        for nni in adj_nnis:
            contains_nni = truth_dag.contains_nni(nni)
            contains_parent = truth_dag.contains_node(nni.parent())
            contains_child = truth_dag.contains_node(nni.child())
            contains_nni = contains_parent and contains_child
            accepted_nni_cnt += contains_nni
            accepted_nnis[nni] = contains_nni

        print("found_nnis:", accepted_nni_cnt)
        if accepted_nni_cnt == 0:
            for nni in adj_nnis:
                accepted_nnis[nni] = True

        nni_engine.test_eval_adjacent_nnis(accepted_nnis)

        print("accepted_nnis:", nni_engine.accepted_nni_count())
        nni_engine.add_accepted_nnis_to_dag()
        # nni_engine.test_add_adjacent_nnis_to_dag()
        print("check_trees:")
        perc_trees_found = find_trees_in_dag(
            dag, truth_sbn_inst.tree_collection)
        print("perc_trees_found:", perc_trees_found)
        if (perc_trees_found == 1.0):
            found_all_trees = True
        top_trees = tp_engine.build_vector_of_unique_top_trees()
        print("num_top_trees: ", len(top_trees))
        tree_count = 0
        fp = open("test_output.iter-{}.txt".format(iterations), 'w')
        for top_tree in top_trees:
            print(tree_count, dag.tree_to_newick_topology(top_tree))
            fp.write(dag.tree_to_newick_topology(top_tree))
            tree_count += 1
        fp.close()
    print("=== end search ===")

    fp = open("test_output.final.txt", 'w')
    top_tree_cnt = 0
    for top_tree in tp_engine.build_vector_of_unique_top_trees():
        likelihood = inst.compute_tree_likelihood(top_tree)
        parsimony = inst.compute_tree_parsimony(top_tree)

        line = "{}, {}, {}, {}".format(
            top_tree_cnt, likelihood, parsimony, dag.tree_to_newick_topology(
                top_tree)
        )
        fp.write(line)
        print(line)
        top_tree_cnt += 1
    fp.close()


def quick_test(fasta_path_, newick_path_, mmap_path_, final_newick_path_, final_mmap_path_):
    print("=== RUN_NNI_EXPANSION ===")
    print("=== load ===")
    inst = create_gp_inst(fasta_path_, newick_path_, mmap_path_)
    init_gp_inst(inst)
    inst.take_first_branch_length()
    nni_engine = inst.get_nni_engine()
    tp_engine = inst.get_tp_engine()
    dag = inst.get_dag()
    nni_engine.set_tp_likelihood_cutoff_filtering_scheme(
        tp_like_cutoff_threshold)
    nni_engine.set_no_filter(True)
    seed_tree_set = inst.currently_loaded_trees_with_gp_branch_lengths()
    seed_trees
    for tree in seed_tree_set.trees:


    inst_3 = create_gp_inst(fasta_path_, final_newick_path_, final_mmap_path_)
    init_gp_inst(inst_3)
    inst_3.take_first_branch_length()
    dag_3 = inst_3.make_dag()
    dag_3 = inst_3.get_dag()
    final_trees_set = inst_3.currently_loaded_trees_with_gp_branch_lengths()
    final_trees = []
    for tree in final_trees_set.trees:
        final_trees.append(tree)

    inst_2 = create_gp_inst(fasta_path_, newick_path_, final_mmap_path_)
    init_gp_inst(inst_2)
    inst_2.take_first_branch_length()
    nni_engine_2 = inst_2.get_nni_engine()
    tp_engine_2 = inst_2.get_tp_engine()
    dag_2 = inst_2.get_dag()
    nni_engine_2.set_tp_likelihood_cutoff_filtering_scheme(
        tp_like_cutoff_threshold)
    nni_engine_2.set_no_filter(True)

    nni_engine_2.run_init()
    nni_engine_2.run_main_loop()
    nni_engine_2.run_post_loop()

    nni_engine.run_init()
    nni_engine.sync_adjacent_nnis_with_dag()
    adj_nnis = nni_engine.adjacent_nnis()
    nni_engine.filter_pre_update()
    nni_engine.filter_eval_adjacent_nnis()
    nni_engine.filter_post_update()

    data = []

    for nni in adj_nnis:
        sc1 = nni_engine.get_score_by_nni(nni)
        sc2 = nni_engine_2.get_score_by_nni(nni)
        edge_id = dag_2.get_edge_id(nni)
        top_tree = tp_engine_2.get_top_tree_with_edge(edge_id)
        llh = inst_2.compute_tree_likelihood(top_tree)
        prs = inst_2.compute_tree_parsimony(top_tree)
        contains_nni = dag_3.contains_nni(nni)
        contains_parent = dag_3.contains_node(nni.parent())
        contains_child = dag_3.contains_node(nni.child())
        contains_nodes = contains_parent and contains_child

        contains_tree = -np.inf
        topo_cnt = 0
        for final_tree in final_trees:
            if top_tree.compare_by_topology(final_tree):
                contains_tree = topo_cnt
                break
            topo_cnt += 1

        data.append([edge_id, sc1, sc2, llh, prs, contains_nni,
                    contains_parent, contains_child, contains_nodes, contains_tree, dag_3.tree_to_newick_topology(top_tree)])
        print("data:", data)

    df = pd.DataFrame(data)
    print("DataFrame")
    fp = open("test_output.csv", 'w')
    fp.write(df.to_csv())
    fp.close()
    return data


test_data = quick_test(
    fasta_path, seed_newick_path, mmap_path, truth_newick_path, truth_mmap_path)
exit()


print("=== LOAD DATA ===")
inst = load_test_inst()
truth_sbn_inst = load_truth_sbn_inst()
truth_inst = load_truth_inst()
posteriors = load_truth_posterior()
dag = inst.get_dag()
nni_engine = inst.get_nni_engine()
tp_engine = inst.get_tp_engine()

print("=== SET BRANCH LENGTHS ===")
inst.take_first_branch_length()
truth_inst.take_first_branch_length()

nni_engine = inst.get_nni_engine()
nni_engine.set_gp_likelihood_cutoff_filtering_scheme(gp_cutoff_threshold)
nni_engine.set_tp_likelihood_cutoff_filtering_scheme(tp_like_cutoff_threshold)
# nni_engine.set_tp_parsimony_cutoff_filtering_scheme(tp_pars_cutoff_threshold)

print("=== FIND TREES ===")
perc_trees_found = find_trees_in_dag(
    inst.get_dag(), truth_inst.currently_loaded_trees_with_gp_branch_lengths()
)
print("perc_trees_found:", perc_trees_found)

dag = inst.get_dag()
nni_engine = inst.get_nni_engine()
tp_engine = inst.get_tp_engine()
# nni_engine.run_init()
scored_nnis = nni_engine.past_scored_nnis()
print("dag:", dag.node_count(), dag.edge_count())
for i in range(dag.edge_count()):
    edge_id = bito.edge_id(i)
    top_tree = tp_engine.get_top_tree_with_edge(edge_id)
    likelihood = inst.compute_tree_likelihood(top_tree)
    parsimony = inst.compute_tree_parsimony(top_tree)
    score = tp_engine.get_top_tree_likelihood_with_edge(edge_id)
    print(
        "nni: edge_id: {}, score: {}, likelihood: {}, parsimony: {}".format(
            edge_id, score, likelihood, parsimony
        )
    )

# print("=== NNI SEARCH ===")
# run_nni_search(inst, truth_sbn_inst, truth_inst)
# run_test_results(inst, truth_sbn_inst, truth_inst)

print("=== FIND TREES ===")
find_trees_in_dag(
    inst.get_dag(), truth_inst.currently_loaded_trees_with_gp_branch_lengths()
)

unique_trees = tp_engine.build_vector_of_unique_top_trees()
print("unique_trees:", len(unique_trees))
counter = 0
for tree in unique_trees:
    print(tree.to_newick_topology())

loaded_trees = truth_inst.currently_loaded_trees_with_gp_branch_lengths().trees
print("loaded_trees:", len(loaded_trees))
for tree in loaded_trees:
    likelihood = inst.compute_tree_likelihood(tree)
    parsimony = inst.compute_tree_parsimony(tree)
    print("likelihood:", likelihood, "parsimony:", parsimony)
    # print(tree.to_newick_topology())

print("loaded_trees_found:", trees_found)

print("=== COMPARE SCORES ===")
# nni_engine.run_init()
scored_nnis = nni_engine.past_scored_nnis()

print("dag:", dag.node_count(), dag.edge_count())
for i in range(dag.edge_count()):
    edge_id = bito.edge_id(i)
    nni = dag.get_nni(edge_id)
    for other_nni in scored_nnis.keys():
        if nni == other_nni:
            top_tree = tp_engine.get_top_tree_with_edge(edge_id)
            likelihood = inst.compute_tree_likelihood(top_tree)
            parsimony = inst.compute_tree_parsimony(top_tree)
            score = tp_engine.get_top_tree_likelihood_with_edge(edge_id)
            # adj_score = scored_nnis[nni]
            print(
                "nni: edge_id: {}, score: {}, likelihood: {}, parsimony: {}".format(
                    edge_id, score, likelihood, parsimony
                )
            )
