import sys
import argparse
import numpy as np
import scipy.stats as ss
import pandas as pd
import bito
import tempfile
from collections import defaultdict
import os
import time


verbose = False


def print_v(*args):
    """
    Optional verbose printing.
    """
    if verbose:
        print(*args)


def load_dag(fasta_path, newick_path, temp_dir=None):
    """
    Create and return a new bito gp_instance and dag from the given fasta and newick
    files.

    This version is susceptible to the bito taxon label issue when using multiple rooted
    tree instances. For now, call load_dag_with_fake_first instead, which does call this
    method.

    Parameters:
        fasta_path (str): The path for the fasta file.
        newick_path (str): The path for the newick file.
        temp_dir (str): The path of the directory to use when writing the mmap data file
            of the gp_instance. When None, the data is written to the directory
            "_ignore". If multiple instances of nni_search run at the same time, they
            require different paths.
    """
    if temp_dir is not None:
        temp_data_path = os.path.join(temp_dir, "mmap.dat")
    else:
        temp_data_path = "_ignore/mmap.data"
    dag_inst = bito.gp_instance(temp_data_path)
    dag_inst.read_fasta_file(fasta_path)
    dag_inst.read_newick_file(newick_path)
    dag_inst.make_dag()
    dag = dag_inst.get_dag()
    return dag_inst, dag


def load_dag_with_fake_first(fasta_path, newick_path, fake_first_path, temp_dir=None):
    """
    Create and return a new bito gp_instance and dag from the given fasta and newick
    files.

    This method creates a new newick_file by appending the newick string from
    fake_first_path to newick_path, so that trees from different instances are
    comparable. However, if this topology is not in the sdag spanned by the topologies
    of newick path, then this dag is not the one you want.

    Parameters:
        fasta_path (str): The path for the fasta file.
        newick_path (str): The path for the newick file.
        fake_first_path (str):  The path for the first newick entry.
        temp_dir (str): The path of the directory to use when writing the mmap data file
            of the gp_instance. When None, the data is written to the directory
            "_ignore". If multiple instances of nni_search run at the same time, they
            require different paths.
    """
    with tempfile.TemporaryDirectory() as another_temp_dir:
        temp_newick_path = f"{another_temp_dir}/temp.nwk"
        with open(temp_newick_path, "w") as temp_newick_file:
            with open(fake_first_path) as the_file:
                temp_newick_file.write(the_file.readline() + "\n")
            with open(newick_path) as the_file:
                for line in the_file:
                    temp_newick_file.write(line + "\n")
        dag_inst, dag = load_dag(fasta_path, temp_newick_path, temp_dir)
    return dag_inst, dag


def load_trees(fasta_path, newick_path, temp_dir=None):
    """
    Create and return a new bito rooted_instance and list of trees from the given fasta
    and newick files.

    This version is susceptible to the bito taxon label issue when using multiple rooted
    tree instances. For now, call load_trees_with_fake_first instead, which does call this
    method.

    Parameters:
        fasta_path (str): The path for the fasta file.
        newick_path (str): The path for the newick file.
        temp_dir (str): The path of the directory to use as part of the name of the
            rooted_istance. This might be unnecessary, if there's no actual file.
    """
    if temp_dir is not None:
        temp_data_path = os.path.join(temp_dir, "trees")
    else:
        temp_data_path = "trees"
    tree_inst = bito.rooted_instance(temp_data_path)
    tree_inst.read_fasta_file(fasta_path)
    tree_inst.read_newick_file(newick_path)
    trees = tree_inst.tree_collection.trees
    return tree_inst, trees


def load_trees_with_fake_first(fasta_path, newick_path, fake_first_path, temp_dir=None):
    """
    Create and return a new bito rooted_instance and list of trees from the given fasta
    and newick files.

    This method creates a new newick_file by appending the newick string from
    fake_first_path to newick_path, so that trees from different instances are
    comparable. This topology is ommitted from the list of trees that is returned.

    Parameters:
        fasta_path (str): The path for the fasta file.
        newick_path (str): The path for the newick file.
        fake_first_path (str):  The path for the first newick entry.
        temp_dir (str): The path of the directory to use as part of the name of the
            rooted_istance. This might be unnecessary, if there's no actual file.
    """
    with tempfile.TemporaryDirectory() as another_temp_dir:
        temp_newick_path = f"{another_temp_dir}/temp.nwk"
        with open(temp_newick_path, "w") as temp_newick_file:
            with open(fake_first_path) as the_file:
                temp_newick_file.write(the_file.readline() + "\n")
            with open(newick_path) as the_file:
                for line in the_file:
                    temp_newick_file.write(line + "\n")
        tree_inst, trees = load_trees(fasta_path, temp_newick_path, temp_dir)
    return tree_inst, trees[1:]


def build_and_save_pcsp_pp_map(
    fasta_path, posterior_newick_path, posterior_probs_path, output_path
):
    """
    Calculate and write to output_path the posterior probabilities, based on those in
    posterior_probs_path, for the PCSPs in the sDAG spanned by the topologies in
    posterior_newick_path.

    This version is susceptible to the bito taxon label issue when using multiple rooted
    tree instances and the general issue of taxon order in bito bitstrings. For now,
    call build_and_save_pcsp_pp_map_with_fake_first instead, which does not call this
    method.

    Parameters:
        fasta_path (str): The path for the fasta file.
        posterior_newick_path (str): The path for the newick file of posterior trees.
        posterior_probs_path (str): The of the csv file of probabilies for the trees of
            posterior_newick_path.
        output_path (str): The path to write out the csv of pcsp posterior
            probabilities.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        print_v("# load dag...")
        dag_inst, dag = load_dag(fasta_path, posterior_newick_path, temp_dir)
        print_v("# load trees...")
        tree_inst, trees = load_trees(fasta_path, posterior_newick_path, temp_dir)
        print_v("# load pps...")
        pps = load_pps(posterior_probs_path)

        print_v("# build maps...")
        tree_id_map, tree_pp_map, tree_cred_map, _ = build_tree_dicts(trees, pps)
        pcsp_pp_map, pcsp_cred_map = build_pcsp_dicts(
            dag, tree_id_map, tree_pp_map, tree_cred_map
        )

        print_v("pcsp_pp_map:", len(pcsp_pp_map), pcsp_pp_map)
        pcsp_dict = {"parent": [], "child": [], "pcsp_pp": [], "in_cred_set": []}
        for pcsp in pcsp_pp_map:
            parent = pcsp.pcsp_get_parent_subsplit().subsplit_to_string()
            child = pcsp.pcsp_get_child_subsplit().subsplit_to_string()
            pcsp_pp = get_pcsp_pp(pcsp, pcsp_pp_map)
            is_pcsp_in_cred_set = pcsp_cred_map[pcsp]
            pcsp_dict["parent"].append(parent)
            pcsp_dict["child"].append(child)
            pcsp_dict["pcsp_pp"].append(pcsp_pp)
            pcsp_dict["in_cred_set"].append(is_pcsp_in_cred_set)
        df = pd.DataFrame(pcsp_dict)
        df.to_csv(output_path)
    return pcsp_pp_map, pcsp_cred_map


def build_and_save_pcsp_pp_map_with_fake_first(
    fasta_path,
    posterior_newick_path,
    posterior_probs_path,
    fake_first_path,
    output_path,
):
    """
    Calculate and write to output_path the posterior probabilities, based on those in
    posterior_probs_path, for the PCSPs in the sDAG spanned by the topologies in
    posterior_newick_path.

    This version creates the sDAG with the additional topology in fake_first_path, so
    there may be more PCSPs than expected. However, the reported posterior probabilities
    of the expected PCSPs are correct.

    Parameters:
        fasta_path (str): The path for the fasta file.
        posterior_newick_path (str): The path for the newick file of posterior trees.
        posterior_probs_path (str): The of the csv file of probabilies for the trees of
            posterior_newick_path.
        fake_first_path (str): The path for the first newick entry.
        output_path (str): The path to write out the csv of pcsp posterior
            probabilities.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        print_v("# load dag...")
        dag_inst, dag = load_dag_with_fake_first(
            fasta_path, posterior_newick_path, fake_first_path, temp_dir
        )
        print_v("# load trees...")
        tree_inst, trees = load_trees_with_fake_first(
            fasta_path, posterior_newick_path, fake_first_path, temp_dir
        )
        print_v("# load pps...")
        pps = load_pps(posterior_probs_path)

        print_v("# build maps...")
        tree_id_map, tree_pp_map, tree_cred_map, _ = build_tree_dicts(trees, pps)
        pcsp_pp_map, pcsp_cred_map = build_pcsp_dicts(
            dag, tree_id_map, tree_pp_map, tree_cred_map
        )

        print_v("pcsp_pp_map:", len(pcsp_pp_map), pcsp_pp_map)
        pcsp_dict = {"parent": [], "child": [], "pcsp_pp": [], "in_cred_set": []}
        for pcsp in pcsp_pp_map:
            parent = pcsp.pcsp_get_parent_subsplit().subsplit_to_string()
            child = pcsp.pcsp_get_child_subsplit().subsplit_to_string()
            pcsp_pp = get_pcsp_pp(pcsp, pcsp_pp_map)
            is_pcsp_in_cred_set = pcsp_cred_map[pcsp]
            pcsp_dict["parent"].append(parent)
            pcsp_dict["child"].append(child)
            pcsp_dict["pcsp_pp"].append(pcsp_pp)
            pcsp_dict["in_cred_set"].append(is_pcsp_in_cred_set)
        df = pd.DataFrame(pcsp_dict)
        df.to_csv(output_path)
    return pcsp_pp_map, pcsp_cred_map


def build_and_save_subsplit_map_with_fake_first(
    fasta_path,
    posterior_newick_path,
    posterior_probs_path,
    fake_first_path,
    output_path,
):
    """
    Calculate and write to output_path the credibility of subsplits, based on those in
    posterior_probs_path, for the subsplits in the sDAG spanned by the topologies in
    posterior_newick_path.

    This version creates the sDAG with the additional topology in fake_first_path, so
    there may be more subsplits than expected. However, the reported membership in the
    credible set of the expected subsplits are correct. Additionally, this issue does
    not exist when the topology of fake_first_path is contained in the sDAG (this is
    always the case under the current design).

    Parameters:
        fasta_path (str): The path for the fasta file.
        posterior_newick_path (str): The path for the newick file of posterior trees.
        posterior_probs_path (str): The of the csv file of probabilies for the trees of
            posterior_newick_path.
        fake_first_path (str): The path for the first newick entry.
        output_path (str): The path to write out the csv of pcsp posterior
            probabilities.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        print_v("# load dag...")
        dag_inst, dag = load_dag_with_fake_first(
            fasta_path, posterior_newick_path, fake_first_path, temp_dir
        )
        print_v("# load trees...")
        tree_inst, trees = load_trees_with_fake_first(
            fasta_path, posterior_newick_path, fake_first_path, temp_dir
        )
        print_v("# load pps...")
        pps = load_pps(posterior_probs_path)

        print_v("# build maps...")
        tree_id_map, _, tree_cred_map, _ = build_tree_dicts(trees, pps)
        subsplit_cred_map = build_subsplit_dict(dag, tree_id_map, tree_cred_map)

        with open(output_path, "w") as the_file:
            header = "subsplit,in_cred_set\n"
            the_file.write(header)
            subsplit_lines = (
                f"{subsplit.subsplit_to_string()},{is_credible}\n"
                for subsplit, is_credible in subsplit_cred_map.items()
            )
            the_file.writelines(subsplit_lines)

    return subsplit_cred_map



def load_pps(pp_path):
    """
    Returns the list of posterior probabilities in pp_path.
    """
    with open(pp_path, "r") as fp:
        pps = [float(line) for line in fp.readlines()]
    return pps


def load_pcsp_pp_map(pcsp_pp_path):
    """
    Reads the PCSP data in pcsp_pp_path and returns three default dictionaries. The
    first is a dictionary mapping a PCSP (as a bito bitset object) to its posterior
    probability, with default 0.0. The second is a dictionary mapping a PCSP to the
    truth value for the PCSP appearing in a credible set topology. The third is a
    default dictionary, with default False, mapping a PCSP to the truth value of the
    PCSP being found by a search method.

    The purpose of the final dictionary is to take advantage of the fact that the search
    methods enlarge the sDAG of the previous iteration. Once an edge is found in the
    sDAG of an iteration of the search, the edge is present in the sDAGs of all future
    iterations.

    The csv file pcsp_pp_path is expected to have columns "parent", "child", "pcsp_pp",
    and "in_cred_set", as prepared by the methods build_and_save_pcsp_pp_map and
    build_and_save_pcsp_pp_map_with_fake_first. Note that since this is loading the data
    from file, the credible set (95% by default) cannot be changed. I.e., if you want
    edges with a different cutoff for the cumulative density of the credible set
    topologies, then you need to rerun build build-pcsp-map.
    """
    to_subsplit = lambda clades: bito.subsplit(*clades)
    to_pcsp = lambda subsplits: bito.pcsp(*subsplits)

    df = pd.read_csv(pcsp_pp_path)
    parent_clades = df["parent"].str.split(pat="|", expand=True)
    child_clades = df["child"].str.split(pat="|", expand=True)
    df["parent"] = parent_clades.apply(to_subsplit, axis=1)
    df["child"] = child_clades.apply(to_subsplit, axis=1)
    df["pcsp"] = df[["parent", "child"]].apply(to_pcsp, axis=1)
    df.set_index("pcsp", inplace=True)

    pcsp_pp_map = defaultdict(float, df["pcsp_pp"].to_dict())
    pcsp_cred_map = defaultdict(bool, df["in_cred_set"].to_dict())
    pcsp_found_map = defaultdict(bool, {pcsp: False for pcsp in pcsp_pp_map})
    return pcsp_pp_map, pcsp_cred_map, pcsp_found_map


def load_subsplit_map(subsplit_path):
    """
    Reads the subsplit data in subsplit_path and returns two default dictionaries. The
    first is a dictionary mapping a subsplit (as a bito bitset object) to the truth 
    value for the subsplit appearing in a credible set topology. The second is a
    default dictionary, with default False, mapping a subsplit to the truth value of the
    subsplit being found by a search method.

    The purpose of the final dictionary is to take advantage of the fact that the search
    methods enlarge the sDAG of the previous iteration. Once a subsplit is found in the
    sDAG of an iteration of the search, the subsplit is present in the sDAGs of all future
    iterations.

    The csv file should be prepared by the build_and_save_subsplit_map_with_fake_first
    method. Note that since this is loading the data from file, the credible set (95% by 
    default) cannot be changed. I.e., if you want subsplits with a different cutoff for 
    the cumulative density of the credible set topologies, then you need to rerun 
    build-subsplit-map.
    """
    to_subsplit = lambda x: bito.subsplit(*x.split("|"))
    with open(subsplit_path) as the_file:
        # Skip the header.
        the_file.readline()
        subsplit_cred_map = {
            to_subsplit((split := line.split(","))[0]): split[1][0] == "T"
            for line in the_file.readlines()
        }

    subsplit_cred_map = defaultdict(bool, subsplit_cred_map)
    subsplit_found_map = defaultdict(
        bool, {subsplit: False for subsplit in subsplit_cred_map}
    )
    return subsplit_cred_map, subsplit_found_map


def build_tree_dicts(trees, pps, cred_pp=0.95):
    """
    Given a list of bito RootedTrees and list of posterior probabilities of those trees,
    returns four dictionaries. The first is simply a dictionary mapping the list index
    to the corresponding RootedTree, with the list index considered as the tree_id.
    The second is a default dictionary mapping a tree_id to the posterior probability,
    with default 0.0. The third is a default dictionary mapping a tree_id to the truth
    value of the tree being in the credible set, with default False. This credible set
    is the 95% credible set by default, but the user may specify another value. The
    fourth is a default dictionary, with default False, mapping a tree_id to the truth
    value of the tree being found by a search method.

    The purpose of the final dictionary is to take advantage of the fact that the search
    methods enlarge the sDAG of the previous iteration. Once a tree is found in the sDAG
    of an iteration of the search, the tree is present in the sDAGs of all future
    iterations.

    Parameters:
        trees (list): The list of bito RootedTrees.
        pps (list): The list of posterior probabilities.
        cred_pp (float): The credible set cumulative density.
    """
    tree_id_map = dict(enumerate(trees))
    tree_pp_map = defaultdict(float, zip(tree_id_map, pps))
    in_credible_set = np.cumsum(pps) <= cred_pp
    tree_cred_map = defaultdict(bool, enumerate(in_credible_set))
    tree_found_map = defaultdict(bool, {j: False for j in tree_id_map})

    return tree_id_map, tree_pp_map, tree_cred_map, tree_found_map


def build_pcsp_dicts(dag, tree_id_map, tree_pp_map, tree_cred_map):
    """
    Given a bito dag and tree dictionaries, this method constructs and returns the two
    dictionaries mapping the PCSPs (as bito bitsets) of the dag to posterior
    probabilities and truth value of membership of the credible set.

    Parameters:
        dag (bito.dag): The sDAG.
        tree_id_map (dict): The dictionary mapping tree_id to bito RootedTree.
        tree_pp_map (dict): The dictionary mapping tree_id to posterior probability.
        tree_cred_map (dict): The dictionary mapping tree_id to truth value of
            membership of the credible set.
    """
    dag_pcsps = dag.build_set_of_edge_bitsets()
    pcsp_pp_map = defaultdict(float, {pcsp: 0.0 for pcsp in dag_pcsps})
    pcsp_cred_map = defaultdict(bool, {pcsp: False for pcsp in dag_pcsps})

    n_taxa = dag.taxon_count()
    top_subsplit = bito.subsplit("0" * n_taxa, "1" * n_taxa)
    make_root_pcsp = lambda tree: bito.pcsp(top_subsplit, tree.build_subsplit())

    for tree_id, pp in tree_pp_map.items():
        tree = tree_id_map[tree_id]
        is_credible = tree_cred_map[tree_id]
        tree_pcsps = tree.build_vector_of_pcsps()
        tree_pcsps.append(make_root_pcsp(tree))
        for pcsp in tree_pcsps:
            pcsp_pp_map[pcsp] += pp
            if is_credible:
                pcsp_cred_map[pcsp] = True
    return pcsp_pp_map, pcsp_cred_map


def build_subsplit_dict(dag, tree_id_map, tree_cred_map):
    """
    Given a bito dag and tree dictionaries, this method constructs and returns the
    dictionary mapping the subsplits (as bito bitsets) of the dag to truth values for
    membership of the credible set.

    Parameters:
        dag (bito.dag): The sDAG.
        tree_id_map (dict): The dictionary mapping tree_id to bito RootedTree.
        tree_cred_map (dict): The dictionary mapping tree_id to truth value of
            membership of the credible set.
    """

    sdag_nodes = dag.build_set_of_node_bitsets()
    subsplit_cred_map = defaultdict(bool, {node: False for node in sdag_nodes})

    n_taxa = dag.taxon_count()
    root = bito.subsplit("0" * n_taxa, "1" * n_taxa)
    subsplit_cred_map[root] = True

    for tree_id, is_credible in tree_cred_map.items():
        if is_credible:
            tree = tree_id_map[tree_id]
            tree_nodes = tree.build_set_of_subsplits()
            for node in tree_nodes:
                subsplit_cred_map[node] = True
    return subsplit_cred_map


def get_pcsp_pp(nni, pcsp_pp_map):
    """Returns the posterior probability of the PCSP associated to the NNI."""
    pcsp = nni.get_central_edge_pcsp() if type(nni) == bito.nni_op else nni
    return pcsp_pp_map[pcsp]


def get_pcsp_pp_rank(best_nni, scored_nnis, pcsp_pp_map):
    best_pcsp_pp = get_pcsp_pp(best_nni, pcsp_pp_map)
    pcsp_pp_rank = 0
    for nni in scored_nnis:
        nni_pcsp = get_pcsp_pp(nni, pcsp_pp_map)
        if nni_pcsp > best_pcsp_pp:
            pcsp_pp_rank += 1
    return pcsp_pp_rank


def get_credible_edge_count(pcsp_cred_map, pcsp_found_map):
    """
    Returns the number of sdag edges that are in topologies of the credible set and marked as found.
    """
    return sum((pcsp_cred_map[pcsp] for pcsp, found in pcsp_found_map.items() if found))


def get_posterior_edge_count(pcsp_found_map):
    """
    Returns the number of sdag edges that are in topologies of the posterior and marked as found.
    """
    return sum(pcsp_found_map.values())


def get_credible_subsplit_count(subsplit_cred_map, subsplit_found_map):
    """
    Returns the number of sdag nodes or subsplits that are in topologies of the credible set and marked as found.
    """
    return sum(
        (subsplit_cred_map[node] for node, found in subsplit_found_map.items() if found)
    )


def get_posterior_subsplit_count(subsplit_found_map):
    """
    Returns the number of sdag nodes or subsplits that are in topologies of the posterior and marked as found.
    """
    return sum(subsplit_found_map.values())


def get_tree_pp(tree_pp_map, tree_found_map):
    """
    Returns the cumulative posterior probabilities of trees marked as found.
    """
    dag_pp = sum(
        (tree_pp_map[tree_id] for tree_id, found in tree_found_map.items() if found)
    )
    return dag_pp


def get_credible_tree_count(tree_cred_map, tree_found_map):
    """
    Returns the count of topologies in the credible set that are marked as found.
    """
    cred_tree_count = sum(
        (
            1
            for tree_id, cred in tree_cred_map.items()
            if cred and tree_found_map[tree_id]
        )
    )
    return cred_tree_count


def update_found_trees(dag, tree_id_map, tree_found_map):
    """
    Updates the tree_found_map dictionary to mark as found the trees of tree_id_map
    contained in the dag. Additionally, this method returns a list of tree_ids for the
    trees found in the dag that were not found previously.

    Parameters:
        dag (bito.dag): The sDAG.
        tree_id_map (dict): The dictionary mapping tree_id to bito RootedTree.
        tree_found_map (dict): The dictionary mapping tree_id to truth value of the
            tree being already found.
    """
    in_dag = lambda tree_id: dag.contains_tree(tree_id_map[tree_id])
    newly_found_tree_ids = [
        tree_id
        for tree_id, found in tree_found_map.items()
        if not found and in_dag(tree_id)
    ]
    for tree_id in newly_found_tree_ids:
        tree_found_map[tree_id] = True
    return newly_found_tree_ids


def update_found_edges(dag, pcsp_found_map):
    """
    Updates the pcsp_found_map dictionary to mark as found the edges contained in the
    dag.

    Parameters:
        dag (bito.dag): The sDAG.
        pcsp_found_map (dict): The dictionary mapping a PCSP to truth value of the PCSP
            being already found.
    """
    pcsps = dag.build_set_of_edge_bitsets()
    newly_found_pcsps = (
        pcsp for pcsp, found in pcsp_found_map.items() if not found and pcsp in pcsps
    )
    for pcsp in newly_found_pcsps:
        pcsp_found_map[pcsp] = True
    return None


def update_found_nodes(dag, subsplit_found_map):
    """
    Updates the subsplit_found_map dictionary to mark as found the subsplits contained
    in the dag.

    Parameters:
        dag (bito.dag): The sDAG.
        subsplit_found_map (dict): The dictionary mapping a subsplit to truth value of
            the subsplit being already found.
    """
    nodes = dag.build_set_of_node_bitsets()
    newly_found_nodes = (
        node
        for node, found in subsplit_found_map.items()
        if not found and node in nodes
    )
    for node in newly_found_nodes:
        subsplit_found_map[node] = True
    return None


def build_ranked_list(list):
    total_list = [len(list)] * len(list)
    rank_list = total_list - ss.rankdata(list, method="ordinal")
    return rank_list


def init_engine_for_gp_search(dag_inst, args):
    """
    Initialize the generalized-pruning search engine based on command line arguments.
    """
    dag_inst.make_gp_engine()
    dag_inst.make_nni_engine()
    dag_inst.take_first_branch_length()
    nni_engine = dag_inst.get_nni_engine()
    nni_engine.set_include_rootsplits(args.include_rootsplits)
    nni_engine.set_gp_likelihood_cutoff_filtering_scheme(0.0)
    nni_engine.set_top_n_score_filtering_scheme(1)
    if args.use_cutoff:
        nni_engine.set_gp_likelihood_cutoff_filtering_scheme(args.threshold)
    if args.use_dropoff:
        nni_engine.set_gp_likelihood_drop_filtering_scheme(args.threshold)
    if args.use_top_n:
        nni_engine.set_top_n_score_filtering_scheme(args.top_n)
    dag_inst.estimate_branch_lengths(1e-5, 3, True)


def init_engine_for_tp_search(dag_inst, args):
    """Initialize the top-pruning search engine based on command line arguments."""
    dag_inst.make_tp_engine()
    dag_inst.make_nni_engine()
    dag_inst.tp_engine_set_branch_lengths_by_taking_first()
    dag_inst.tp_engine_set_choice_map_by_taking_first()
    nni_engine = dag_inst.get_nni_engine()
    nni_engine.set_include_rootsplits(args.include_rootsplits)
    nni_engine.set_tp_likelihood_cutoff_filtering_scheme(0.0)
    nni_engine.set_top_n_score_filtering_scheme(1)
    if args.use_cutoff:
        nni_engine.set_tp_likelihood_cutoff_filtering_scheme(args.threshold)
    if args.use_dropoff:
        nni_engine.set_tp_likelihood_drop_filtering_scheme(args.threshold)
    if args.use_top_n:
        nni_engine.set_top_n_score_filtering_scheme(args.top_n)


def init_results_file(
    file_path,
    dag,
    tree_id_map,
    tree_pp_map,
    tree_cred_map,
    tree_found_map,
    pcsp_cred_map,
    pcsp_found_map,
    subsplit_cred_map,
    subsplit_found_map,
):
    """
    Create a csv file at file_path and writes to this file the header line and sdag
    stats for iteration 0 of an nni-search.
    """
    tree_count = int(dag.topology_count())
    node_count = dag.node_count()
    edge_count = dag.edge_count()
    posterior_tree_ids = update_found_trees(dag, tree_id_map, tree_found_map)
    posterior_tree_ids = f'"{posterior_tree_ids}"'
    tree_pp = get_tree_pp(tree_pp_map, tree_found_map)
    cred_tree_count = get_credible_tree_count(tree_cred_map, tree_found_map)

    update_found_edges(dag, pcsp_found_map)
    cred_edge_count = get_credible_edge_count(pcsp_cred_map, pcsp_found_map)
    posterior_edge_count = get_posterior_edge_count(pcsp_found_map)

    update_found_nodes(dag, subsplit_found_map)
    cred_node_count = get_credible_subsplit_count(subsplit_cred_map, subsplit_found_map)
    posterior_node_count = get_posterior_subsplit_count(subsplit_found_map)

    header = "iter,acc_nni_id,acc_nni_count,score,tree_pp,pcsp_pp,pcsp_pp_rank,"
    header += "node_count,cred_node_count,posterior_node_count,edge_count,"
    header += "cred_edge_count,posterior_edge_count,tree_count,cred_tree_count,"
    header += "posterior_tree_ids,adj_nni_count,new_nni_count,llhs_computed,parent,"
    header += "child\n"
    first_row = f"0,,0,,{tree_pp},,,{node_count},{cred_node_count},"
    first_row += f"{posterior_node_count},{edge_count},{cred_edge_count},"
    first_row += f"{posterior_edge_count},{tree_count},{cred_tree_count},"
    first_row += f"{posterior_tree_ids},0,0,0,,\n"
    with open(file_path, "w") as the_file:
        the_file.write(header)
        the_file.write(first_row)
    return None


def write_results_line(file_path, *args):
    """
    Appends a line of sdag stats to the csv at file_path.
    """
    with open(file_path, "a") as the_file:
        the_file.write(",".join(map(str, args)) + "\n")
    return None


def nni_search(args):
    """
    Perform an nni-search based on command line arguments.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        start_time = time.time()

        print("nni_search")
        print_v("# load trees...")
        # tree_inst, trees = load_trees(args.fasta, args.posterior_newick)
        tree_inst, trees = load_trees_with_fake_first(
            args.fasta, args.posterior_newick, args.seed_newick, temp_dir
        )

        print_v("# load pps...")
        pps = load_pps(args.pp_csv)

        print_v("# build maps...")
        tree_id_map, tree_pp_map, tree_cred_map, tree_found_map = build_tree_dicts(
            trees, pps
        )
        pcsp_pp_map, pcsp_cred_map, pcsp_found_map = load_pcsp_pp_map(args.pcsp_pp_csv)
        subsplit_cred_map, subsplit_found_map = load_subsplit_map(args.subsplit_csv)

        print_v("# load dag...")
        dag_inst, dag = load_dag(args.fasta, args.seed_newick, temp_dir)

        print_v("# init engine...")
        if args.tp:
            init_engine_for_tp_search(dag_inst, args)
        if args.gp:
            init_engine_for_gp_search(dag_inst, args)
        nni_engine = dag_inst.get_nni_engine()
        nni_engine.run_init(True)

        init_results_file(
            args.output,
            dag,
            tree_id_map,
            tree_pp_map,
            tree_cred_map,
            tree_found_map,
            pcsp_cred_map,
            pcsp_found_map,
            subsplit_cred_map,
            subsplit_found_map,
        )

        log_iteration_frequency = args.log_freq

        prev_nni_count = 0
        llhs_computed = 0
        current_time = time.time() - start_time
        iteration_times = [(0, current_time)]
        for iter_count in range(1, args.iter_max + 1):
            # run iteration of search
            log_this_iter = (iter_count % log_iteration_frequency) == 0

            print_v(f"\n# iter_count: {iter_count} of {args.iter_max}...")
            print_v(f"# dag: {dag.node_count()} nodes, {dag.edge_count()} edges")

            if args.gp:
                # Generalized-pruning specific code.
                init_engine_for_gp_search(dag_inst, args)
                nni_engine = dag_inst.get_nni_engine()
                nni_engine.run_init(True)
            if not args.pcsp:
                # Top- and generalized-pruning code.
                nni_engine.graft_adjacent_nnis_to_dag()
                nni_engine.filter_pre_update()
                nni_engine.filter_eval_adjacent_nnis()
                nni_engine.filter_post_update()
                nni_engine.filter_process_adjacent_nnis()
                nni_engine.remove_all_graft_nnis_from_dag()
                nni_engine.add_accepted_nnis_to_dag(False)
                scored_nnis = nni_engine.scored_nnis()
                accepted_nni_count = len(nni_engine.accepted_nnis())
                print_v("# scored_nnis:", len(scored_nnis))
                print_v("# accepted_nnis:", nni_engine.accepted_nni_count())
            else:
                # Old code block for choosing the highest posterior density nni. Need to
                # clean this up at some point.
                for nni_id, nni in enumerate(nni_engine.adjacent_nnis()):
                    nni_pcsp_pp = get_pcsp_pp(nni, pcsp_pp_map)
                    nni_pcsp_map[nni] = nni_pcsp_pp
                    if nni_pcsp_pp > best_pcsp:
                        best_pcsp = nni_pcsp_pp
                        best_nni = nni
                dag.add_node_pair(best_nni.get_parent(), best_nni.get_child())

            # add entry to final data
            new_nni_count = len(scored_nnis) - (prev_nni_count - 1)
            if args.tp:
                llhs_computed += new_nni_count
            if args.gp:
                llhs_computed += len(scored_nnis)

            if log_this_iter:
                if args.log_time_only:
                    current_time = time.time() - start_time
                    iteration_times.append((iter_count, current_time))
                else:
                    tree_count = int(dag.topology_count())
                    node_count = dag.node_count()
                    edge_count = dag.edge_count()
                    posterior_tree_ids = update_found_trees(
                        dag, tree_id_map, tree_found_map
                    )
                    posterior_tree_ids = f'"{posterior_tree_ids}"'
                    tree_pp = get_tree_pp(tree_pp_map, tree_found_map)
                    cred_tree_count = get_credible_tree_count(
                        tree_cred_map, tree_found_map
                    )
                    update_found_edges(dag, pcsp_found_map)
                    cred_edge_count = get_credible_edge_count(
                        pcsp_cred_map, pcsp_found_map
                    )
                    posterior_edge_count = get_posterior_edge_count(pcsp_found_map)
                    update_found_nodes(dag, subsplit_found_map)
                    cred_node_count = get_credible_subsplit_count(
                        subsplit_cred_map, subsplit_found_map
                    )
                    posterior_node_count = get_posterior_subsplit_count(
                        subsplit_found_map
                    )

                    adjacent_nni_count = nni_engine.adjacent_nni_count()
                    print_v("# dag_tree_pp:", tree_pp)

                    for nni_id, nni in enumerate(nni_engine.accepted_nnis()):
                        pcsp_pp = get_pcsp_pp(nni, pcsp_pp_map)
                        pcsp_pp_rank = get_pcsp_pp_rank(nni, scored_nnis, pcsp_pp_map)
                        parent = nni.get_parent().subsplit_to_string()
                        child = nni.get_child().subsplit_to_string()

                        write_results_line(
                            args.output,
                            iter_count,
                            nni_id,
                            accepted_nni_count,
                            scored_nnis[nni],
                            tree_pp,
                            pcsp_pp,
                            pcsp_pp_rank,
                            node_count,
                            cred_node_count,
                            posterior_node_count,
                            edge_count,
                            cred_edge_count,
                            posterior_edge_count,
                            tree_count,
                            cred_tree_count,
                            posterior_tree_ids,
                            adjacent_nni_count,
                            new_nni_count,
                            llhs_computed,
                            parent,
                            child,
                        )

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
        print_v(f"# final dataframe written to {args.output}")

        if args.log_time_only:
            with open(f"{args.output}.time.csv", "w") as the_file:
                the_file.write("iter,time\n")
                for iter, seconds in iteration_times:
                    the_file.write(f"{iter},{seconds}\n")
    return None


def nni_df_add_entry(**kwargs):
    """...old code for debugging..."""
    info_dict = {
        "iter": [iter_count] * len(scored_nnis),
        "nni_id": [],
        "llh": [],
        "tree_pp": [],
        "pcsp_pp": [],
        "llh_rank": [],
        "tree_pp_rank": [],
        "pcsp_pp_rank": [],
        "parent": [],
        "child": [],
    }
    for nni_id, (nni, llh_score) in enumerate(scored_nnis):
        print_v("# nni_id", nni_id, "of", len(scored_nnis), "...")
        _, nni_dag = load_dag(args.fasta, args.seed_newick)
        for old_nni in final_dict["added_nni"]:
            nni_dag.add_node_pair(old_nni.get_parent(), old_nni.get_child())
        nni_dag.add_node_pair(nni.get_parent(), nni.get_child())
        if include_tree_pp:
            nni_tree_pp = get_tree_pp(nni_dag, tree_id_map, tree_pp_map)
        else:
            nni_tree_pp = -np.inf
        nni_pcsp_pp = get_pcsp_pp(nni, pcsp_pp_map)
        info_dict["nni_id"].append(nni_id)
        info_dict["parent"].append(nni.get_parent().subsplit_to_string())
        info_dict["child"].append(nni.get_child().subsplit_to_string())
        info_dict["tree_pp"].append(nni_tree_pp)
        info_dict["pcsp_pp"].append(nni_pcsp_pp)
        info_dict["llh"].append(llh_score)

    print_v("# build data...")
    if args.include_tree_pp:
        info_dict["tree_pp_rank"] = build_ranked_list(my_dict["tree_pp"])
    else:
        del info_dict["tree_pp"]
        del info_dict["tree_pp_rank"]
    info_dict["pcsp_pp_rank"] = build_ranked_list(my_dict["pcsp_pp"])
    info_dict["llh_rank"] = build_ranked_list(my_dict["llh"])
    df = pd.DataFrame(my_dict)
    dfs[iter_count] = df
    df.to_csv(f"{args.output}.{iter_count}")
    return df


def score_info(nni_engine):
    """...old debug code?..."""
    scored_nnis = nni_engine.scored_nnis()
    n = len(scored_nnis)
    the_max = None if n == 0 else max(scored_nnis.values())
    print(f"There are {n} scored NNIs, max score {the_max}.")

    adjacent_nnis = nni_engine.adjacent_nnis()
    n = len(adjacent_nnis)
    # the_max = None if n == 0 else max(scored_nnis.values())
    print(f"There are {n} adjacent NNIs.")

    return None


############
### MAIN ###
############


def main_arg_parse(args):
    parser = argparse.ArgumentParser(
        description="Tools for performing NNI systematic search."
    )
    parser.add_argument("-v", "--verbose", help="verbose", type=int, default=1)
    parser.add_argument("-p", "--profiler", action="store_true", help="profile program")
    subparsers = parser.add_subparsers(title="programs", dest="program")

    # nni search
    subparser1 = subparsers.add_parser(
        "nni-search", help="Perform systematic NNI search."
    )
    subparser1.add_argument("fasta", help="fasta file", type=str)
    subparser1.add_argument(
        "seed_newick", help="newick file for initial trees in DAG", type=str
    )
    subparser1.add_argument(
        "credible_newick", help="newick file for trees in credible posterior", type=str
    )
    subparser1.add_argument(
        "posterior_newick",
        help="newick file for trees in empirical posterior",
        type=str,
    )
    subparser1.add_argument(
        "pp_csv",
        help="csv file containing the posterior weights of the trees from credible_trees",
        type=str,
    )
    subparser1.add_argument(
        "pcsp_pp_csv",
        help="csv file containing the per-PCSP posterior weights",
        type=str,
    )
    subparser1.add_argument(
        "subsplit_csv",
        help="csv file containing the posterior and credible subsplits",
        type=str,
    )
    # search method group
    group = subparser1.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--gp",
        action="store_true",
        help="Selects best NNI according to Generalized Pruning.",
    )
    group.add_argument(
        "--tp",
        action="store_true",
        help="Selects best NNI according to Top Pruning via Likelihood.",
    )
    group.add_argument(
        "--pcsp",
        action="store_true",
        help="Selects best NNI according to per-PCSP cumulative posterior.",
    )
    # search scheme group
    group = subparser1.add_mutually_exclusive_group(required=False)
    group.add_argument(
        "--use-top-n", action="store_true", help="Selects the top N scoring NNIs."
    )
    group.add_argument(
        "--use-cutoff",
        action="store_true",
        help="Selects all NNIs scoring above a given threshold.",
    )
    group.add_argument(
        "--use-dropoff",
        action="store_true",
        help="Selects all NNIs scoring above a given dropoff below best NNI.",
    )
    # search scheme arguments
    subparser1.add_argument(
        "--top-n", help="number of NNIs to accept per iteration", type=int, default=1
    )
    subparser1.add_argument(
        "--threshold", help="cutoff/dropoff threshold", type=float, default=0.0
    )
    # options
    subparser1.add_argument(
        "-o", "--output", help="output file", type=str, default="results.nni_search.csv"
    )
    subparser1.add_argument(
        "--iter-max", help="number of NNI search iterations", type=int, default=10
    )
    subparser1.add_argument(
        "--log-freq", help="frequency to log stats for the search", type=int, default=1
    )
    subparser1.add_argument(
        "--log-time-only",
        help="log run-time and no other statistics",
        action="store_true",
    )
    subparser1.add_argument(
        "--nni-info", help="Give additional NNI info per iteration", type=str
    )
    # subparser1.add_argument(
    #     '--tree-pp', action='store_true' help='Compute DAG pre-tree posterior every iteration')
    subparser1.add_argument(
        "--include-rootsplits",
        action="store_true",
        help="Whether to include rootsplits in NNI search",
    )

    # pcsp map builder
    subparser2 = subparsers.add_parser("build-pcsp-map", help="Build per-PCSP map.")
    subparser2.add_argument("fasta", help="fasta file", type=str)
    subparser2.add_argument(
        "posterior_newick",
        help="newick file for trees in empirical posterior",
        type=str,
    )
    subparser2.add_argument(
        "pp_csv",
        help="csv file containing the posterior weights of the trees from posterior_newick",
        type=str,
    )
    subparser2.add_argument(
        "seed_newick",
        help="newick file for initial trees in DAG, needed now to avoid a taxon label problem",
        type=str,
    )
    # options
    subparser2.add_argument(
        "-o",
        "--output",
        help="output file",
        type=str,
        default="results.pcsp_pp_map.csv",
    )

    # subsplit map builder
    subparser3 = subparsers.add_parser("build-subsplit-map", help="Build subsplit map.")
    subparser3.add_argument("fasta", help="fasta file", type=str)
    subparser3.add_argument(
        "posterior_newick",
        help="newick file for trees in empirical posterior",
        type=str,
    )
    subparser3.add_argument(
        "pp_csv",
        help="csv file containing the posterior weights of the trees from posterior_newick",
        type=str,
    )
    subparser3.add_argument(
        "seed_newick",
        help="newick file for initial trees in DAG, needed now to avoid a taxon label problem",
        type=str,
    )
    # options
    subparser3.add_argument(
        "-o",
        "--output",
        help="output file",
        type=str,
        default="results.subsplit_node_map.csv",
    )
    # run parser
    parsed_args = parser.parse_args(args)
    args_dict = vars(parsed_args)
    return parsed_args


if __name__ == "__main__":
    print_v("# begin...")
    args = main_arg_parse(sys.argv[1:])
    verbose = args.verbose > 0
    if args.program == "nni-search":
        if args.tp:
            nni_search(args)
        if args.gp:
            nni_search(args)
        if args.pcsp:
            nni_search(args)
    if args.program == "build-pcsp-map":
        build_and_save_pcsp_pp_map_with_fake_first(
            args.fasta,
            args.posterior_newick,
            args.pp_csv,
            args.seed_newick,
            args.output,
        )
    if args.program == "build-subsplit-map":
        build_and_save_subsplit_map_with_fake_first(
            args.fasta,
            args.posterior_newick,
            args.pp_csv,
            args.seed_newick,
            args.output,
        )
    print_v("# ...done")
