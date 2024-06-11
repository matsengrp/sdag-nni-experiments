#!/usr/bin/env python

import click
import pandas as pd
import pickle
import subprocess
import tempfile
import os
from ete3 import Tree
from collections import OrderedDict, defaultdict
from Bio import Phylo
from io import StringIO



def mcmc_treeprob(filename):
    """
    Given a Phylo parsable file, in nexus format, extracts the trees and their weights.
    Returns a triple for a dictionary mapping tree names to ete3 Tree objects based on
    the newick representations, a list of tree names, and a list of tree weights. The
    two lists are in matching order.

    The expected use is with filename being a trprobs file from a MrBayes run.
    """
    data_type = "nexus"
    tree_format = "newick"
    mcmc_samp_tree_stats = Phylo.parse(filename, data_type)
    mcmc_samp_tree_dict = OrderedDict()
    mcmc_samp_tree_name = []
    mcmc_samp_tree_wts = []

    for tree in mcmc_samp_tree_stats:
        handle = StringIO()
        Phylo.write(tree, handle, tree_format)
        mcmc_samp_tree_dict[tree.name] = Tree(handle.getvalue().strip())
        handle.close()
        mcmc_samp_tree_name.append(tree.name)
        mcmc_samp_tree_wts.append(tree.weight)

    return mcmc_samp_tree_dict, mcmc_samp_tree_name, mcmc_samp_tree_wts


def combine_trprobs_files(file_paths):
    """
    Given a collection of trprobs files listed in file_paths, find the average weight
    for each tree among the trprobs files. This method returns a triple of a dictionary
    mapping a tree name to a ete3 Tree object, a list of tree names, and a list of tree
    weights. The two lists are in matching order. The tree names will not match the
    tree names in the trprobs files. This method assumes it makes sense to simply
    average the weight across the trprobs files.
    """
    number_of_trpobs = len(file_paths)
    tree_dict_total = OrderedDict()
    tree_dict_map_total = defaultdict(float)
    tree_names_total = []
    tree_wts_total = []
    n_samp_tree = 0
    for file_path in file_paths:
        tree_dict_rep, tree_name_rep, tree_wts_rep = mcmc_treeprob(file_path)

        for j, name in enumerate(tree_name_rep):
            tree_id = tree_dict_rep[name].get_topology_id()
            if tree_id not in tree_dict_map_total:
                n_samp_tree += 1
                tree_names_total.append("tree_{}".format(n_samp_tree))
                tree_dict_total[tree_names_total[-1]] = tree_dict_rep[name]

            tree_dict_map_total[tree_id] += tree_wts_rep[j]

    for key in tree_dict_map_total:
        tree_dict_map_total[key] /= number_of_trpobs

    tree_wts_total = [
        tree_dict_map_total[tree_dict_total[name].get_topology_id()]
        for name in tree_names_total
    ]

    return tree_dict_total, tree_names_total, tree_wts_total


def reroot(trees, reroot_number):
    rooted_trees = []
    with tempfile.TemporaryDirectory() as tmpdir:
        unrooted_path = os.path.join(tmpdir, "unrooted")
        rooted_path = os.path.join(tmpdir, "rooted")
        with open(unrooted_path, "w") as temp_file:
            temp_file.writelines(trees)
        subprocess.check_call(
            f"nw_reroot {unrooted_path} {reroot_number} | nw_order - > {rooted_path}",
            shell=True,
        )
        with open(rooted_path, "r") as temp_file:
            for line in temp_file:
                rooted_trees.append(line.strip())
    return rooted_trees


@click.command()
@click.argument("base_file_path")
@click.argument("var_path_start", type=int)
@click.argument("var_path_stop", type=int)
@click.argument("file_name")
@click.argument("reroot_number")
@click.argument("out_pickle_path")
def run(
    base_file_path,
    var_path_start,
    var_path_stop,
    file_name,
    reroot_number,
    out_pickle_path,
):
    """
    Given a collection of trprobs located at:
        base_file_path{j}/file_name     for var_path_start <= j <= var_path_stop
    take the average weight of the trees in these files and write this back to file.
    This writes to out_file_path each tree with its posterior probablity; one entry per
    line with the tree in newick string format ending in a semicolon, following by the
    posterior probability. This dumps a pickle file to out_pickle_path containing a
    pair of the list of trees in the 95% credible set and a dictionary for tree to
    posterior probability.
    """
    file_paths = [
        f"{base_file_path}{j}/{file_name}"
        for j in range(var_path_start, 1 + var_path_stop)
    ]

    tree_names_to_objects, tree_names, tree_weights = combine_trprobs_files(file_paths)
    newick_strings = [
        tree_names_to_objects[tree].write(format=9) for tree in tree_names
    ]
    newick_strings = reroot(newick_strings, reroot_number)
    newick_to_pp = dict(zip(newick_strings, tree_weights))

    df = pd.DataFrame(data={"tree": newick_strings, "weight": tree_weights})
    df.sort_values("weight", ascending=False, inplace=True)
    df["total_pp"] = df["weight"].cumsum()
    cred_set = df[df["total_pp"] <= 0.95]["tree"].to_list()

    with open(out_pickle_path, "wb") as the_file:
        pickle.dump((newick_to_pp, cred_set), the_file)


if __name__ == "__main__":
    run()
