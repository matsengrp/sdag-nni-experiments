#!/usr/bin/env python
import click
from ete3 import Tree


@click.command()
@click.argument("newick_input_path")
@click.argument("out_path")
def run_it(newick_input_path, out_path):
    """
    Some code expects trees, but doesn't actually use the branchg lengths. This script
    adds branch lengths for that purpose.
    """

    dummy_length = 0.1
    with open(newick_input_path, "r") as the_topology_file:
        with open(out_path, "w") as the_tree_file:
            for line in the_topology_file:
                tree = Tree(line)
                for node in tree.traverse():
                    node.dist = dummy_length
                the_tree_file.write(tree.write(format=5) + "\n")
    return None


if __name__ == "__main__":
    run_it()
