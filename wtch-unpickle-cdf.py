#!/usr/bin/env python

import pickle
import os
import click


@click.command()
@click.argument("read_path")
@click.argument("credible_tree_write_path")
@click.argument("pp_tree_write_path")
@click.argument("pp_value_write_path")
def cli(read_path, credible_tree_write_path, pp_tree_write_path, pp_value_write_path):
    """
    Given a pickle file of a tuple from a MrBayes run, write out three files for the
    credible set topologies, the topologies with pp-values, and the pp-values of the
    those topologies.
    """

    pp_dict, cred_tree_list = pickle.load(open(read_path, "rb"))
    with open(credible_tree_write_path, "wt") as out_file:
        for tree in cred_tree_list:
            out_file.write(tree + "\n")
    # The use of 7 digits is specific to the current setup of the default 6 digits from
    # MrBayes and the extra digit due to the 10 runs.  
    formatter = lambda tree : "{:.7f}".format(pp_dict[tree])
    with open(pp_tree_write_path, "wt") as tree_out_file:
        with open(pp_value_write_path, "wt") as pp_out_file:
            for tree in pp_dict.keys():
                tree_out_file.write(tree + "\n")
                pp_out_file.write(formatter(tree) + "\n")


if __name__ == "__main__":
    cli()
