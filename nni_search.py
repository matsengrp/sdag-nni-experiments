import sys
import argparse
import itertools
import numpy as np
import scipy.stats as ss
import pandas as pd
import bito
import cProfile

verbose = True


def print_v(*args):
    if verbose:
        print(*args)


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


def load_pps(pp_path):
    pps = []
    with open(pp_path, 'r') as fp:
        for line in fp.readlines():
            pps.append(float(line))
    return pps


def load_pcsp_pp_map(pcsp_pp_path):
    pcsp_pp_map = {}
    df = pd.read_csv(pcsp_pp_path)
    for index, row in df.iterrows():
        parent = row["parent"].split("|")
        parent = bito.subsplit(parent[0], parent[1])
        child = row["child"].split("|")
        child = bito.subsplit(child[0], child[1])
        pcsp = bito.pcsp(parent, child)
        pp = float(row["pcsp_pp"])
        pcsp_pp_map[pcsp] = pp
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


def get_tree_pp(dag, tree_id_map, tree_pp_map):
    dag_pp = 0.0
    for tree_id in tree_id_map:
        tree = tree_id_map[tree_id]
        if (dag.contains_tree(tree)):
            pp = tree_pp_map[tree_id]
            dag_pp += pp
    return dag_pp


def get_pcsp_pp(nni, pcsp_pp_map):
    if type(nni) == bito.nni_op:
        pcsp = nni.get_central_edge_pcsp()
    else:
        pcsp = nni
    if pcsp in pcsp_pp_map.keys():
        return pcsp_pp_map[pcsp]
    else:
        return 0.0


def get_pcsp_pp_rank(best_nni, scored_nnis, pcsp_pp_map):
    best_pcsp_pp = get_pcsp_pp(best_nni, pcsp_pp_map)
    pcsp_pp_rank = 0
    for nni in scored_nnis:
        nni_pcsp = get_pcsp_pp(nni, pcsp_pp_map)
        if nni_pcsp > best_pcsp_pp:
            pcsp_pp_rank += 1
    return pcsp_pp_rank


def get_credible_edge_count(dag, pcsp_pp_map):
    cred_edge_count = 0
    noncred_edge_count = 0
    pcsps = dag.build_sorted_vector_of_edge_bitsets()
    for pcsp in pcsps:
        if pcsp in pcsp_pp_map:
            cred_edge_count += 1
        else:
            noncred_edge_count += 1
    return (cred_edge_count, noncred_edge_count)


def build_ranked_list(list):
    total_list = [len(list)] * len(list)
    rank_list = total_list - ss.rankdata(list, method='ordinal')
    return rank_list


def init_engine_for_gp_search(dag_inst, args):
    dag_inst.make_gp_engine()
    dag_inst.make_nni_engine()
    dag_inst.take_first_branch_length()
    nni_engine = dag_inst.get_nni_engine()
    nni_engine.set_include_rootsplits(args.include_rootsplits)
    nni_engine.set_gp_likelihood_cutoff_filtering_scheme(0.0)
    nni_engine.set_top_n_score_filtering_scheme(1)
    dag_inst.estimate_branch_lengths(1e-5, 3)


def init_engine_for_tp_search(dag_inst, args):
    dag_inst.make_tp_engine()
    dag_inst.make_nni_engine()
    dag_inst.tp_engine_set_branch_lengths_by_taking_first()
    dag_inst.tp_engine_set_choice_map_by_taking_first()
    nni_engine = dag_inst.get_nni_engine()
    nni_engine.set_include_rootsplits(args.include_rootsplits)
    nni_engine.set_tp_likelihood_cutoff_filtering_scheme(0.0)
    nni_engine.set_top_n_score_filtering_scheme(1)


def build_and_save_pcsp_pp_map(args):
    print_v("# load dag...")
    dag_inst, dag = load_dag(args.fasta, args.credible_newick)
    print_v("# load trees...")
    tree_inst, trees = load_trees(args.fasta, args.credible_newick)
    print_v("# load pps...")
    pps = load_pps(args.pp_csv)
    print_v("# build maps...")
    tree_id_map = build_tree_id_map(trees)
    tree_pp_map = build_tree_pp_map(tree_id_map, pps)
    pcsp_pp_map = build_pcsp_pp_map(dag, tree_id_map, tree_pp_map)
    print_v("pcsp_pp_map:", len(pcsp_pp_map), pcsp_pp_map)
    my_dict = {
        'parent': [],
        'child': [],
        'pcsp_pp': []
    }
    for pcsp in pcsp_pp_map:
        parent = pcsp.pcsp_get_parent_subsplit().subsplit_to_string()
        child = pcsp.pcsp_get_child_subsplit().subsplit_to_string()
        pcsp_pp = get_pcsp_pp(pcsp, pcsp_pp_map)
        my_dict["parent"].append(parent)
        my_dict["child"].append(child)
        my_dict["pcsp_pp"].append(pcsp_pp)
    df = pd.DataFrame(my_dict)
    df.to_csv(args.output)
    return pcsp_pp_map


def nni_search(args):
    print_v("# load trees...")
    tree_inst, trees = load_trees(args.fasta, args.credible_newick)
    print_v("# load pps...")
    pps = load_pps(args.pp_csv)
    print_v("# build maps...")
    tree_id_map = build_tree_id_map(trees)
    tree_pp_map = build_tree_pp_map(tree_id_map, pps)
    pcsp_pp_map = load_pcsp_pp_map(args.pcsp_pp_csv)
    final_dict = final_data_init()

    print_v("# load dag...")
    dag_inst, _ = load_dag(args.fasta, args.seed_newick)
    dag = dag_inst.get_dag()
    print_v("# init engine...")
    if args.tp:
        init_engine_for_tp_search(dag_inst, args)
    if args.gp:
        init_engine_for_gp_search(dag_inst, args)
    nni_engine = dag_inst.get_nni_engine()
    nni_engine.run_init(True)

    iter_count = 0
    prev_nni_count = 0
    llhs_computed = 0
    while iter_count < args.iter_max:
        # run iteration of search
        print_v(f"# iter_count: {iter_count} of {args.iter_max}...")
        print_v("# dag:", dag.node_count(), dag.edge_count())
        tree_pp = get_tree_pp(dag, tree_id_map, tree_pp_map)
        print_v("# dag_tree_pp:", tree_pp)
        if args.gp:
            nni_engine.run_init(True)
        nni_engine.graft_adjacent_nnis_to_dag()
        nni_engine.filter_pre_update()
        nni_engine.filter_eval_adjacent_nnis()
        nni_engine.filter_post_update()
        nni_engine.filter_process_adjacent_nnis()
        nni_engine.remove_all_graft_nnis_from_dag()
        nni_engine.add_accepted_nnis_to_dag(False)
        scored_nnis = nni_engine.scored_nnis()
        accepted_nni_count = len(nni_engine.accepted_nnis())
        print_v("# scored_nnis:", len(scored_nnis), scored_nnis)
        print_v("# accepted_nnis:", nni_engine.accepted_nni_count())

        # add entry to final data
        new_nni_count = len(scored_nnis) - (prev_nni_count - 1)
        if args.tp:
            llhs_computed += new_nni_count
        if args.gp:
            llhs_computed += len(scored_nnis)
        tree_count = dag.topology_count()
        cred_edge_count, noncred_edge_count = get_credible_edge_count(
            dag, pcsp_pp_map)
        for (nni_id, nni) in enumerate(nni_engine.accepted_nnis()):
            pcsp_pp = get_pcsp_pp(nni, pcsp_pp_map)
            pcsp_pp_rank = get_pcsp_pp_rank(
                nni, scored_nnis, pcsp_pp_map)
            final_dict['iter'].append(iter_count)
            final_dict['acc_nni_id'].append(nni_id)
            final_dict['acc_nni_count'].append(accepted_nni_count)
            final_dict['score'].append(scored_nnis[nni])
            final_dict['tree_pp'].append(tree_pp)
            final_dict['pcsp_pp'].append(pcsp_pp)
            final_dict['pcsp_pp_rank'].append(pcsp_pp_rank)
            final_dict['node_count'].append(dag.node_count())
            final_dict['edge_count'].append(dag.edge_count())
            final_dict['cred_edge_count'].append(cred_edge_count)
            final_dict['tree_count'].append(tree_count)
            final_dict['adj_nni_count'].append(nni_engine.adjacent_nni_count())
            final_dict['new_nni_count'].append(new_nni_count)
            final_dict['llhs_computed'].append(llhs_computed)
            final_dict['parent'].append(nni.get_parent().subsplit_to_string())
            final_dict['child'].append(nni.get_child().subsplit_to_string())
        df = pd.DataFrame(final_dict)
        df.to_csv(args.output)

        # add entry for additional NNI info
        if args.nni_info:
            nni_info_df = nni_df_add_entry(**locals())

        # bookkeeping post iteration
        if len(nni_engine.accepted_nnis()) == 0:
            print_v("# NO ACCEPTED NNIS")
            break
        nni_engine.run_post_loop()
        nni_engine.sync_adjacent_nnis_with_dag()
        prev_nni_count = len(scored_nnis)
        iter_count += 1

    print_v("# final dataframe:")
    df = pd.DataFrame(final_dict)
    pd.set_option('display.max_colwidth', None)
    print_v(df)
    df.to_csv(args.output)
    return final_dict


def nni_search__choose_by_best_pcsp_pp(args):
    print_v("# parse args...")
    args = parse_args__nni_search()
    print_v("# load trees...")
    tree_inst, trees = load_trees(args.fasta, args.credible_newick)
    print_v("# load pps...")
    pps = load_pps(args.pp_csv)
    print_v("# build maps...")
    tree_id_map = build_tree_id_map(trees)
    tree_pp_map = build_tree_pp_map(tree_id_map, pps)
    pcsp_pp_map = load_pcsp_pp_map(args.pcsp_pp_csv)

    final_dict = {
        'iter': [],
        'tree_pp': [],
        'added_nni': [],
        'parent': [],
        'child': [],
    }

    print_v("# load dag...")
    dag_inst, _ = load_dag(args.fasta, args.seed_newick)
    dag = dag_inst.get_dag()
    dag_inst.make_nni_engine()
    nni_engine = dag_inst.get_nni_engine()
    nni_engine.set_include_rootsplits(include_rootsplits)

    for iter_count in range(iter_max):
        print_v(f"# iter_count: {iter_count} of {iter_max}...")
        print_v("# dag:", dag.node_count(), dag.edge_count())
        tree_pp = get_tree_pp(dag, tree_id_map, tree_pp_map)
        print_v("# dag_tree_pp:", tree_pp)
        nni_pcsp_map = {}

        best_pcsp = -np.inf
        nni_engine.sync_adjacent_nnis_with_dag()
        adj_nnis = nni_engine.adjacent_nnis()
        print_v("adj_nnis:", adj_nnis)
        for nni_id, nni in enumerate(adj_nnis):
            nni_pcsp_pp = get_pcsp_pp(nni, pcsp_pp_map)
            nni_pcsp_map[nni] = nni_pcsp_pp
            if nni_pcsp_pp > best_pcsp:
                best_pcsp = nni_pcsp_pp
                best_nni = nni
        print_v("# nni_pcsp_map:", nni_pcsp_map)
        dag.add_node_pair(best_nni.get_parent(), best_nni.get_child())
        final_dict['iter'].append(iter_count)
        final_dict['tree_pp'].append(tree_pp)
        final_dict['added_nni'].append(best_nni)
        final_dict['parent'].append(best_nni.get_parent().subsplit_to_string())
        final_dict['child'].append(best_nni.get_child().subsplit_to_string())

    del final_dict['added_nni']
    print_v("# final dataframe:")
    df = pd.DataFrame(final_dict)
    print_v(df)
    df.to_csv(f"_ignore/test.best_by_pcsp_pp.final.csv")
    return final_dict


def nni_df_add_entry(**kwargs):
    info_dict = {
        'iter': [iter_count] * len(scored_nnis),
        'nni_id': [],
        'llh': [],
        'tree_pp': [],
        'pcsp_pp': [],
        'llh_rank': [],
        'tree_pp_rank': [],
        'pcsp_pp_rank': [],
        'parent': [],
        'child': [],
    }
    for nni_id, (nni, llh_score) in enumerate(scored_nnis):
        print_v("# nni_id", nni_id, "of", len(scored_nnis), "...")
        _, nni_dag = load_dag(args.fasta, args.seed_newick)
        for old_nni in final_dict['added_nni']:
            nni_dag.add_node_pair(
                old_nni.get_parent(), old_nni.get_child())
        nni_dag.add_node_pair(nni.get_parent(), nni.get_child())
        if include_tree_pp:
            nni_tree_pp = get_tree_pp(
                nni_dag, tree_id_map, tree_pp_map)
        else:
            nni_tree_pp = -np.inf
        nni_pcsp_pp = get_pcsp_pp(nni, pcsp_pp_map)
        info_dict['nni_id'].append(nni_id)
        info_dict['parent'].append(
            nni.get_parent().subsplit_to_string())
        info_dict['child'].append(nni.get_child().subsplit_to_string())
        info_dict['tree_pp'].append(nni_tree_pp)
        info_dict['pcsp_pp'].append(nni_pcsp_pp)
        info_dict['llh'].append(llh_score)

    print_v("# build data...")
    if args.include_tree_pp:
        info_dict['tree_pp_rank'] = build_ranked_list(my_dict['tree_pp'])
    else:
        del info_dict['tree_pp']
        del info_dict['tree_pp_rank']
    info_dict['pcsp_pp_rank'] = build_ranked_list(my_dict['pcsp_pp'])
    info_dict['llh_rank'] = build_ranked_list(my_dict['llh'])
    df = pd.DataFrame(my_dict)
    dfs[iter_count] = df
    df.to_csv(f"{args.output}.{iter_count}")
    return df


def final_data_init():
    final_data = {
        'iter': [],
        'acc_nni_id': [],
        'acc_nni_count': [],
        'score': [],
        'tree_pp': [],
        'pcsp_pp': [],
        'pcsp_pp_rank': [],
        'node_count': [],
        'edge_count': [],
        'cred_edge_count': [],
        'tree_count': [],
        'adj_nni_count': [],
        'new_nni_count': [],
        'llhs_computed': [],
        'parent': [],
        'child': []
    }
    return final_data

############
### MAIN ###
############


def main_arg_parse(args):
    parser = argparse.ArgumentParser(
        description='Tools for performing NNI systematic search.')
    parser.add_argument('-v', '--verbose',
                        help='verbose', type=int, default=1)
    parser.add_argument('-p', '--profiler',
                        action='store_true', help='profile program')
    subparsers = parser.add_subparsers(title='programs', dest='program')

    # nni search
    subparser1 = subparsers.add_parser(
        'nni_search', help='Perform iterative NNI search.')
    subparser1.add_argument('fasta', help='fasta file', type=str)
    subparser1.add_argument(
        'seed_newick', help='newick file for initial trees in DAG', type=str)
    subparser1.add_argument(
        'credible_newick', help='newick file for trees in credible posterior', type=str)
    subparser1.add_argument(
        'pp_csv', help='csv file containing the posterior weights of the trees from credible_trees', type=str)
    subparser1.add_argument(
        'pcsp_pp_csv', help='csv file containing the per-PCSP posterior weights', type=str)
    group = subparser1.add_mutually_exclusive_group(required=True)
    group.add_argument('--gp', action='store_true',
                       help='Use generalized pruning.')
    group.add_argument('--tp', action='store_true',
                       help='Use top pruning via log likelihood.')
    group.add_argument('--pcsp', action='store_true',
                       help='Use by selecting the NNI with the best per-PCSP podsterior.')
    # options
    subparser1.add_argument('-o', '--output', help='output file', type=str,
                            default='results.nni_search.csv')
    subparser1.add_argument(
        '--iter-max', help='number of NNI search iterations', type=int, default=10)
    subparser1.add_argument(
        '--nni-info', help='Give additional NNI info per iteration', type=str)
    subparser1.add_argument(
        '--tree-pp', help='Compute DAG pre-tree posterior every iteration', type=str)
    subparser1.add_argument("--include-rootsplits", action='store_true',
                            help='Whether to include rootsplits in NNI search')

    # pcsp map builder
    subparser2 = subparsers.add_parser(
        'build_pcsp_map', help='Build per-PCSP map.')
    subparser2.add_argument('fasta', help='fasta file', type=str)
    subparser2.add_argument(
        'credible_newick', help='newick file for trees in credible posterior', type=str)
    subparser2.add_argument(
        'pp_csv', help='csv file containing the posterior weights of the trees from credible_trees', type=str)
    # options
    subparser2.add_argument('-o', '--output', help='output file', type=str,
                            default='results.pcsp_pp_map.csv')
    # run parser
    parsed_args = parser.parse_args(args)
    args_dict = vars(parsed_args)
    return parsed_args


if __name__ == "__main__":
    print_v("# begin...")
    args = main_arg_parse(sys.argv[1:])
    verbose = (args.verbose > 0)
    if (args.program == 'nni_search'):
        if (args.tp):
            nni_search(args)
        if (args.gp):
            nni_search(args)
        if (args.pcsp):
            nni_search__choose_by_best_pcsp_pp(args)
    if (args.program == 'build_pcsp_map'):
        build_and_save_pcsp_pp_map(args)
    print_v("# ...done")
