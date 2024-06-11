#!/usr/bin/env python
import pickle
import os
from collections import namedtuple
import tempfile
import subprocess
import numpy as np
import click
from asyncio.subprocess import DEVNULL

TreeData = namedtuple("TreeData", "pp_dict tree_set")


def tree_data_of_path(tree_pickle_path):
    pp_dict, tree_ci_list = pickle.load(open(tree_pickle_path, "rb"))
    return TreeData(pp_dict, set(tree_ci_list))


def optimize_branch_lengths(topology_set, sequence_file_path, sort=True):
    """
    Returns the list of trees in topology_set with optimal branch lengths, optionally
    ordered by likelihood (highest likelihood first).
            Parameters:
                    topology_set (set): A set of topologies represented as strings of
                        their Newick tree format (without branch lengths).
                    sequence_file_path (string): The file containing the sequencing
                        data for the tree tips.
            Returns:
                    optimized_trees (list): The list of trees from topology_set. Each
                    tree is represented as a string of their Newick tree format (with
                    optimal branch lengths). This list is optionally ordered according
                    to the log-likelihood from iqtree, with maximum likelihood first.
    """
    topology_file_name = "topology.nwk"
    optimized_file_name = "optimized-trees.nwk"
    likelihood_file_name = "tree-likelihoods.temp"
    with tempfile.TemporaryDirectory() as tmpdir:
        topology_path = os.path.join(tmpdir, topology_file_name)
        optimized_path = os.path.join(tmpdir, optimized_file_name)
        likelihood_path = os.path.join(tmpdir, likelihood_file_name)

        for topology in topology_set:
            with open(topology_path, "w") as fp:
                fp.write(topology + "\n")

            # Run iqtree on the tree. The command below may throw a runtime warning:
            #       OMP: Info #271: omp_set_nested routine deprecated, please use
            #       omp_set_max_active_levels instead.
            # This is an issue with the installed version of openMP. Since it is only a
            # warning, people say just ignore it. We redirect to suppress the message.
            subprocess.check_call(
                f"iqtree --redo --quiet -s {sequence_file_path} -te {topology_path} "
                + "-m jc69",
                shell=True,
                stderr=DEVNULL,
            )

            # Add the tree with optimized branch lengths somewhere temporary, same for
            # the likelihoods. Grabbing the likelihood requires some magic to find it
            # in the full report.
            subprocess.check_call(
                f"cat {sequence_file_path}.treefile >> {optimized_path}", shell=True
            )
            subprocess.check_call(
                f'grep "Log-likelihood of the tree: " {sequence_file_path}.iqtree | '
                + "sed -e 's/^Log\\-likelihood of the tree: \\(.*\\) (.*$/\\1/;'"
                + f">> {likelihood_path}",
                shell=True,
            )

        # Load all the tree info into variables.
        with open(optimized_path, "r") as the_file:
            optimized_trees = the_file.read().splitlines()

        if sort:
            tree_likelihoods = np.loadtxt(likelihood_path, delimiter=";", dtype=float)
            indices_for_sort = np.flip(np.argsort(tree_likelihoods))
            optimized_trees = [optimized_trees[j] for j in indices_for_sort]

        # Clean up local files.
        for extension in ["ckp.gz", "iqtree", "log", "treefile"]:
            subprocess.check_call(f"rm {sequence_file_path}.{extension}", shell=True)

        return optimized_trees


@click.command()
@click.argument("topology_path")
@click.argument("fasta_path")
@click.argument("output_path")
@click.option("--sort", default=True)
def wrapper_for_tree_optimizing(topology_path, fasta_path, output_path, sort=True):
    with open(topology_path, "r") as the_file:
        topology_data = the_file.read().splitlines()

    optimized_trees = optimize_branch_lengths(topology_data, fasta_path, sort)

    with open(output_path, "w") as the_output_file:
        for tree in optimized_trees:
            the_output_file.write(f"{tree}\n")


if __name__ == "__main__":
    wrapper_for_tree_optimizing()
