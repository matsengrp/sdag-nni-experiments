import sys
import itertools
import numpy as np
import scipy.stats as ss
import pandas as pd
import bito


def parse_args__build_pcsp_pp_map():
    req_args = 3
    if (len(sys.argv) != (req_args + 1)):
        print("Usage: <i:fasta> <i:mb_credible_newick> <i:pp_csv>")
        exit()
    args = {}
    args["fasta"] = sys.argv[1]
    args["credible_newick"] = sys.argv[2]
    args["pp_csv"] = sys.argv[3]
    return args


def parse_args__nni_search():
    req_args = 5
    if (len(sys.argv) != (req_args + 1)):
        print("Usage: <i:fasta> <i:seed_newick> <i:credible_newick> <i:pp_csv> <i:pcsp_pp_csv>")
        exit()

    args = {}
    args["fasta"] = sys.argv[1]
    args["seed_newick"] = sys.argv[2]
    args["credible_newick"] = sys.argv[3]
    args["pp_csv"] = sys.argv[4]
    args["pcsp_pp_csv"] = sys.argv[5]
    return args


def build_pcsp_weight_table(pcsp_weight_path):
    pcsp_weight_path = sys.argv[1]
    pcsp_weight_header = ["pcsp", "pp"]
    pcsp_wt_df = pd.read_csv(pcsp_weight_path, names=pcsp_weight_header)
    pcsp_wt_df = pcsp_wt_df[1:]

    parent_subsplits = []
    child_subsplits = []
    for pcsp in pcsp_wt_df["pcsp"]:
        pcsp_bitset = bito.bitset(pcsp)
        parent_subsplits.append(
            pcsp_bitset.pcsp_get_parent_subsplit().subsplit_to_string())
        child_subsplits.append(
            pcsp_bitset.pcsp_get_child_subsplit().subsplit_to_string())
    pcsp_wt_df["parent"] = parent_subsplits
    pcsp_wt_df["child"] = child_subsplits
    pcsp_wt_df = pcsp_wt_df.drop(["pcsp"], axis=1)
    return pcsp_wt_df


def load_dag(fasta_path, newick_path):
    dag_inst = bito.gp_instance("_ignore/mmap.data")
    dag_inst.read_fasta_file(fasta_path)
    dag_inst.read_newick_file(newick_path)
    dag_inst.make_dag()
    dag = dag_inst.get_dag()
    return dag_inst, dag


def load_trees(fasta_path, newick_path):
    tree_inst = bito.rooted_instance("trees")
    tree_inst.read_fasta_file(fasta_path)
    tree_inst.read_newick_file(newick_path)
    trees = tree_inst.tree_collection.trees
    return tree_inst, trees


def save_pcsp_pp_map(pcsp_pp_map, fileout_path):
    fp = open(fileout_path, 'w')
    for pcsp in pcsp_pp_map:
        parent = pcsp.pcsp_get_parent_subsplit()
        child = pcsp.pcsp_get_child_subsplit()
        pp = pcsp_pp_map[pcsp]
        line = f"{parent.subsplit_to_string()} {child.subsplit_to_string()} {pp}\n"
        print(line)
        fp.write(line)
    fp.close()


def load_pps(pp_path):
    pps = []
    with open(pp_path, 'r') as fp:
        for line in fp.readlines():
            pps.append(float(line))
    return pps


def load_pcsp_pp_map(pcsp_pp_path):
    pcsp_pp_map = {}
    fp = open(pcsp_pp_path, 'r')
    for line in fp.readlines():
        fields = line.strip().split(" ")
        parent = fields[0].split("|")
        parent = bito.subsplit(parent[0], parent[1])
        child = fields[1].split("|")
        child = bito.subsplit(child[0], child[1])
        pcsp = bito.pcsp(parent, child)
        pp = float(fields[2])
        pcsp_pp_map[pcsp] = pp
    fp.close()
    return pcsp_pp_map


def build_tree_id_map(trees):
    tree_id_map = {}
    for id, tree in enumerate(trees):
        tree_id_map[id] = tree
    return tree_id_map


def build_tree_pp_map(tree_id_map, pps):
    tree_pp_map = {}
    for (tree_id, pp) in zip(tree_id_map, pps):
        tree_pp_map[tree_id] = pp
    return tree_pp_map


def build_pcsp_pp_map(dag, tree_id_map, tree_pp_map):
    dag_pcsps = dag.build_sorted_vector_of_edge_bitsets()
    pcsp_pp_map = {}
    for pcsp in dag_pcsps:
        pcsp_pp_map[pcsp] = 0.0
    for tree_id in tree_pp_map:
        tree = tree_id_map[tree_id]
        pp = tree_pp_map[tree_id]
        tree_pcsps = tree.build_vector_of_pcsps()
        for pcsp in tree_pcsps:
            pcsp_pp_map[pcsp] += pp
    return pcsp_pp_map


def get_dag_pp(dag, tree_id_map, tree_pp_map):
    dag_pp = 0.0
    for tree_id in tree_id_map:
        tree = tree_id_map[tree_id]
        if (dag.contains_tree(tree)):
            pp = tree_pp_map[tree_id]
            dag_pp += pp
    return dag_pp


def get_pcsp_pp(nni, pcsp_pp_map):
    pcsp = nni.get_central_edge_pcsp()
    if pcsp in pcsp_pp_map.keys():
        return pcsp_pp_map[pcsp]
    else:
        return 0.0


def build_ranked_list(list):
    total_list = [len(list)] * len(list)
    rank_list = total_list - ss.rankdata(list, method='ordinal')
    return rank_list


def get_dag_pp_for_all_nnis():
    print("# parse args...")
    args = parse_args__nni_search()

    print("# load trees...")
    tree_inst, trees = load_trees(args["fasta"], args["credible_newick"])
    print("# load pps...")
    pps = load_pps(args["pp_csv"])
    print("# build maps...")
    tree_id_map = build_tree_id_map(trees)
    tree_pp_map = build_tree_pp_map(tree_id_map, pps)
    pcsp_pp_map = load_pcsp_pp_map(args["pcsp_pp_csv"])

    added_nnis = []
    dfs = {}

    iter_count = 0
    dfs = {}
    for iter_count in range(10):
        print("# iter_count:", iter_count)
        dag_inst, _ = load_dag(args["fasta"], args["seed_newick"])
        dag = dag_inst.get_dag()
        for nni in added_nnis:
            dag.add_node_pair(nni.get_parent(), nni.get_child())
        print("# dag:", dag.node_count(), dag.edge_count())
        print("# dag_pp:", get_dag_pp(dag, tree_id_map, tree_pp_map))

        dag_inst.make_tp_engine()
        dag_inst.make_nni_engine()
        dag_inst.tp_engine_set_branch_lengths_by_taking_first()
        dag_inst.tp_engine_set_choice_map_by_taking_first()
        nni_engine = dag_inst.get_nni_engine()

        nni_engine.set_tp_likelihood_cutoff_filtering_scheme(0.0)
        nni_engine.set_top_n_score_filtering_scheme(1)
        nni_engine.run_init(True)
        nni_engine.graft_adjacent_nnis_to_dag()
        nni_engine.filter_pre_update()
        nni_engine.filter_eval_adjacent_nnis()
        nni_engine.filter_post_update()
        nni_engine.filter_process_adjacent_nnis()
        nni_engine.remove_all_graft_nnis_from_dag()
        scored_nni_map = nni_engine.scored_nnis()
        scored_nnis = sorted(scored_nni_map.items(),
                             key=lambda x: x[1], reverse=True)

        print("# accepted_nnis:", nni_engine.accepted_nni_count())
        for nni in nni_engine.accepted_nnis():
            print(nni.get_parent().subsplit_to_string(),
                  nni.get_child().subsplit_to_string(), scored_nni_map[nni])
        # print("# scored_nnis:", len(scored_nnis))
        # for (nni, score) in scored_nnis:
        #     print(nni.get_parent().subsplit_to_string(),
        #           nni.get_child().subsplit_to_string(), score)

        my_dict = {
            'iter': [iter_count] * len(scored_nnis),
            'nni_id': [],
            'parent': [],
            'child': [],
            'tp_llh': [],
            'tree_pp': [],
            'pcsp_pp': [],
            'tp_llh_rank': [],
            'tree_pp_rank': [],
            'pcsp_pp_rank': []
        }
        for nni_id, (nni, llh_score) in enumerate(scored_nnis):
            _, nni_dag = load_dag(args["fasta"], args["seed_newick"])
            for old_nni in added_nnis:
                nni_dag.add_node_pair(
                    old_nni.get_parent(), old_nni.get_child())
            nni_dag.add_node_pair(nni.get_parent(), nni.get_child())
            nni_dag_pp = get_dag_pp(nni_dag, tree_id_map, tree_pp_map)
            nni_pcsp_pp = get_pcsp_pp(nni, pcsp_pp_map)
            my_dict['nni_id'].append(nni_id)
            my_dict['parent'].append(nni.get_parent().subsplit_to_string())
            my_dict['child'].append(nni.get_child().subsplit_to_string())
            my_dict['tree_pp'].append(nni_dag_pp)
            my_dict['pcsp_pp'].append(nni_pcsp_pp)
            my_dict['tp_llh'].append(llh_score)

        tmp = [len(scored_nnis)] * len(scored_nnis)
        my_dict['tree_pp_rank'] = build_ranked_list(my_dict['tree_pp'])
        my_dict['pcsp_pp_rank'] = build_ranked_list(my_dict['pcsp_pp'])
        my_dict['tp_llh_rank'] = build_ranked_list(my_dict['tp_llh'])
        df = pd.DataFrame(my_dict)
        dfs[iter_count] = df
        print("# dataframe:")
        print(df)
        df.to_csv(f"_ignore/ds1.top10.results.{iter_count}.csv")

        for nni in nni_engine.accepted_nnis():
            added_nnis.append(nni)

    print("added_nnis:")
    for (id, nni) enumerate(added_nnis):
        print(id, nni)
    return


def build_and_save_pcsp_pp_map():
    print("# parse args...")
    args = parse_args__build_pcsp_pp_map()
    print("# load dag...")
    dag_inst, dag = load_dag(args["fasta"], args["credible_newick"])
    print("# load trees...")
    tree_inst, trees = load_trees(args["fasta"], args["credible_newick"])
    print("# load pps...")
    pps = load_pps(args["pp_csv"])
    print("# build maps...")
    tree_id_map = build_tree_id_map(trees)
    tree_pp_map = build_tree_pp_map(tree_id_map, pps)
    pcsp_pp_map = build_pcsp_pp_map(dag, tree_id_map, tree_pp_map)


def nni_search__add_top_nni():
    print("# parse args...")
    args = parse_args__nni_search()
    print("# load dag...")
    cred_inst, cred_dag = load_dag(args["fasta"], args["credible_newick"])


def dag_nni_search():
    print("# parse args...")
    args = parse_args__nni_search()
    print("# load dag...")
    dag_inst, dag = load_dag(args["fasta"], args["seed_newick"])
    print("# load trees...")
    tree_inst, trees = load_trees(args["fasta"], args["credible_newick"])
    print("# load pps...")
    pps = load_pps(args["pp_csv"])
    print("# build maps...")
    tree_id_map = build_tree_id_map(trees)
    tree_pp_map = build_tree_pp_map(tree_id_map, pps)
    # pcsp_pp_map = load_pcsp_pp_map(args["pcsp_pp_csv"])

    dag_pp = get_dag_pp(dag, tree_id_map, tree_pp_map)
    print("dag_pp:", dag_pp)
    added_nnis = []

    print("# nni search...")
    dag_inst.make_tp_engine()
    dag_inst.make_nni_engine()
    nni_engine = dag_inst.get_nni_engine()

    # nni_engine.set_tp_likelihood_cutoff_filtering_scheme(0)
    nni_engine.run_init(False)
    # nni_engine.reset_all_nnis()
    # nni_engine.sync_adjacent_nnis_with_dag()
    # nni_engine.prep_eval_engine()
    dag_inst.tp_engine_set_branch_lengths_by_taking_first()
    dag_inst.tp_engine_set_choice_map_by_taking_first(False)
    nni_engine.set_top_n_score_filtering_scheme(1, True)
    # nni_engine.filter_init()

    # print("Adjacent NNIs:", nni_engine.adjacent_nnis())

    for i in range(10):
        nni_engine.run_main_loop()
        print("Scored NNIs:", nni_engine.scored_nnis())
        nni_engine.run_post_loop()

    return dag_inst, dag

    ############
    ### MAIN ###
    ############


if __name__ == "__main__":
    print("# build_pp_ranking_by_pcsp...")
    # build_and_save_pcsp_pp_map()
    get_dag_pp_for_all_nnis()
    print("# ...done")
