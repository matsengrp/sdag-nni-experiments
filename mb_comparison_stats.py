import scipy.stats as ss
import pickle
import bito
import pandas as pd
import tempfile
import os
import pathlib
import time
from collections import namedtuple
import click
from nni_search import (
    build_tree_dicts,
    get_credible_edge_count,
    get_credible_subsplit_count,
    get_credible_tree_count,
    get_posterior_edge_count,
    get_posterior_subsplit_count,
    get_tree_pp,
    load_pcsp_pp_map,
    load_subsplit_map,
    load_pps,
    load_trees,
    update_found_edges,
    update_found_nodes,
    update_found_trees,
)


GoldenData = namedtuple("GoldenData", "pp_dict credible_set")


def golden_data_of_path(golden_pickle_path):
    pp_dict, tree_ci_list = pickle.load(open(golden_pickle_path, "rb"))
    return GoldenData(pp_dict, set(tree_ci_list))


def topology_set_to_path(topology_set, topologies_path):
    "Write a set of topologies."
    with open(topologies_path, "w") as topologies_file:
        for topology in topology_set:
            topologies_file.write(topology + "\n")


def topology_set_of_path(topologies_path):
    "Read a set of topologies."
    with open(topologies_path) as topologies_file:
        return {t.strip() for t in topologies_file}


def mcmc_df_of_topology_sequence(topology_seq_path, seen_topology_dir, golden):
    pathlib.Path(seen_topology_dir).mkdir(exist_ok=True)
    df = pd.read_csv(topology_seq_path, sep="\t", names=["dwell_count", "topology"])
    # The set of topologies seen so far.
    seen = set()
    seen_list = []
    # For writing the topologies to file, we use a list, so that the order is
    # consistent. This avoids issues with bito loading multiple newick files with
    # different first lines.
    first_time = []

    for topology in df["topology"]:
        if topology in seen:
            first_time.append(False)
        else:
            first_time.append(True)
            seen.add(topology)
            seen_list.append(topology)

    write_path = seen_topology_dir + f"/topologies-seen.{len(seen_list)}.nwk"
    topology_set_to_path(seen_list, write_path)
    df["first_time"] = first_time
    df["support_size"] = df["first_time"].cumsum()
    df["mcmc_iters"] = df.dwell_count.cumsum().shift(1, fill_value=0)
    df["pp"] = df["topology"].apply(lambda t: golden.pp_dict.get(t, 0.0))
    df["total_pp"] = (df["pp"] * df["first_time"]).cumsum()
    df["in_credible_set"] = df["topology"].apply(golden.credible_set.__contains__)
    df["credible_set_found"] = (df["in_credible_set"] & df["first_time"]).cumsum()
    df["credible_set_frac"] = df["credible_set_found"] / len(golden.credible_set)
    df.drop_duplicates(subset="support_size", keep="first", inplace=True)
    df.drop(
        columns=["dwell_count", "topology", "first_time", "pp", "in_credible_set"],
        inplace=True,
    )
    return df


def sdag_results_df_of(
    seen_trees,
    max_topology_count,
    all_seen_path,
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
    Calculate sdag stats for topologies from Mr. Bayes.

    Parameters:
        seen_trees (list): The list of bito RootedTrees found by the short run of Mr.
            Bayes, in the order of seen topologies.
        max_topology_count (int): The maximum number of topologies for which to compute
            sDAG statistics.
        all_seen_path (str): The path of the newick files for all seen topologies,
            sorted by sorder order.
        tree_id_map (dict): The dictionary mapping tree_id to bito RootedTree.
        tree_pp_map (dict): The dictionary mapping tree_id to posterior probability.
        tree_cred_map (dict): The dictionary mapping tree_id to truth value of
            membership of the credible set.
        tree_found_map (dict): The dictionary mapping tree_id to truth value of the
            tree being already found.
        pcsp_cred_map (dict): The dictionary mapping a PCSP (as a bito bitset object) to
            the truth value of the PCSP appearing in a topology of the credible set.
        pcsp_found_map (dict): The dictionary mapping a PCSP (as a bito bitset object)
            to the truth value of the PCSP being already found.
        subsplit_cred_map (dict): The dictionary mapping a subsplit (as a bito bitset object) to
            the truth value of the subsplit appearing in a topology of the credible set.
        subsplit_found_map (dict): The dictionary mapping a subsplit (as a bito bitset object)
            to the truth value of the subsplit being already found.
    """
    rows_for_df = []
    with tempfile.TemporaryDirectory() as tmpdir:
        temp_data_path = os.path.join(tmpdir, "mmap.dat")
        tree_increases_dag = True
        seen_path = os.path.join(tmpdir, "current_seen.nwk")
        with open(all_seen_path) as all_seen, open(seen_path, "a") as seen_file:
            for topology in range(1, max_topology_count + 1):
                current_line = next(all_seen)
                print(f"Tree {topology} is necessary: {tree_increases_dag}")
                if tree_increases_dag:
                    seen_file.write(current_line + "\n")
                    seen_file.flush()
                    start_time = time.time()
                    inst = bito.gp_instance(temp_data_path)
                    inst.read_newick_file(seen_path)
                    inst.make_dag()
                    dag = inst.get_dag()
                    dag_construction_time = time.time() - start_time

                    update_found_trees(dag, tree_id_map, tree_found_map)
                    tree_count = int(dag.topology_count())
                    node_count = dag.node_count()
                    edge_count = dag.edge_count()
                    tree_pp = get_tree_pp(tree_pp_map, tree_found_map)
                    cred_tree_count = get_credible_tree_count(
                        tree_cred_map, tree_found_map
                    )
                    update_found_edges(dag, pcsp_found_map)
                    posterior_edge_count = get_posterior_edge_count(pcsp_found_map)
                    credible_edge_count = get_credible_edge_count(
                        pcsp_cred_map, pcsp_found_map
                    )
                    update_found_nodes(dag, subsplit_found_map)
                    posterior_subsplit_count = get_posterior_subsplit_count(
                        subsplit_found_map
                    )
                    credible_subsplit_count = get_credible_subsplit_count(
                        subsplit_cred_map, subsplit_found_map
                    )

                row = [topology, node_count, edge_count, tree_count, cred_tree_count]
                row.extend([tree_pp, posterior_edge_count, credible_edge_count])
                row.extend([posterior_subsplit_count, credible_subsplit_count])
                row.append(dag_construction_time)
                rows_for_df.append(row)
                if topology < max_topology_count:
                    next_tree = seen_trees[topology]
                    tree_increases_dag = not dag.contains_tree(next_tree)

    the_df = pd.DataFrame(
        rows_for_df,
        columns=[
            "support_size",
            "sdag_node_count",
            "sdag_edge_count",
            "sdag_topos_total",
            "sdag_topos_in_credible",
            "sdag_total_pp",
            "sdag_edges_in_posterior",
            "sdag_edges_in_credible",
            "sdag_nodes_in_posterior",
            "sdag_nodes_in_credible",
            "sdag_build_time",
        ],
    )
    credible_count = sum(tree_cred_map.values())
    the_df["sdag_credible_set_frac"] = the_df.sdag_topos_in_credible / credible_count
    return the_df


def restrict_to_dag(newick_path, tree_dicts):
    """
    Given a list of tree dictionaries, the first of which is a tree_id_map (integer to
    bito RootedTree), remove the entries of the dictionaries for trees that are not in
    the sDAG spanned by the topologies in newick_path. The modification to the
    dictionaries is done in-place
    """
    tree_id_map = tree_dicts[0]
    with tempfile.TemporaryDirectory() as tmpdir:
        print("Building maximal sdag.")
        inst = bito.gp_instance(os.path.join(tmpdir, "mmap.dat"))
        inst.read_newick_file(newick_path)
        inst.make_dag()
        max_dag = inst.get_dag()
        print(f"Checking trees for membership.")
        trees_to_remove = [
            tree_id
            for tree_id, tree in tree_id_map.items()
            if not max_dag.contains_tree(tree)
        ]
    tree_count = len(tree_id_map)
    remove_count = len(trees_to_remove)
    print(f"Removing {remove_count} of the {tree_count} posterior trees.")

    for tree_id in trees_to_remove:
        for tree_dict in tree_dicts:
            tree_dict.pop(tree_id)

    return None


@click.command()
@click.argument("golden_pickle_path", type=str)
@click.argument("topology_sequence_path", type=str)
@click.argument("fasta_path", type=str)
@click.argument("seed_newick_path", type=str)
@click.argument("posterior_newick_path", type=str)
@click.argument("pp_csv", type=str)
@click.argument("pcsp_pp_csv", type=str)
@click.argument("subsplit_csv", type=str)
@click.argument("out_path", type=str)
@click.option("--skip_sdag_stats", is_flag=True)
def run(
    golden_pickle_path,
    topology_sequence_path,
    fasta_path,
    seed_newick_path,
    posterior_newick_path,
    pp_csv,
    pcsp_pp_csv,
    subsplit_csv,
    out_path,
    skip_sdag_stats=False,
):
    """
    Calculate posterior and sdag stats for a short Mr. Bayes run. Optionally skip the
    sdag stats, which take longer to compute.

    Parameters:
        golden_pickle_path (str): Path to the pickled Mr. Bayes empirical posterior.
        topology_sequence_path (str): Path to the rerooted topologies of the short Mr.
            Bayes run.
        fasta_path (str): Path to the fasta file.
        seed_newick_path (str): Path to the file with the newick of the common topology
            for multiple bito object instances.
        posterior_newick_path (str): Path to the newick file of the posterior trees from
            the long run of Mr. Bayes.
        pp_csv (str): Path to the csv file with the posterior probababilities of the
            trees in posterior_newick_path.
        pcsp_pp_csv (str): Path to the csv file with the posterior probababilities of
            the edges of the sDAG spanned by the topologies in posterior_newick_path.
        subsplit_csv (str): Path to the csv file with subsplits and their membership
            in the credible set, for the subsplits of the sDAG spanned by the
            topoloiogies in posterior_newick_path.
        out_path (str): Write out path for the data.
    """
    # Gather the topologies visited by Mr Bayes.
    golden = golden_data_of_path(golden_pickle_path)

    seen_topology_dir = os.path.dirname(topology_sequence_path) + "/topologies-seen"
    accumulation_df = mcmc_df_of_topology_sequence(
        topology_sequence_path, seen_topology_dir, golden
    )
    total_seen_count = accumulation_df.support_size.max()
    print(f"Accumulation dataframe built. Seen topology count: {total_seen_count}")

    if skip_sdag_stats:
        accumulation_df.to_csv(out_path)
        return None

    # Get all posterior trees, then restrict to those in the largest sdag of MCMC trees.
    # The rest are never visited and don't contribute to any of the stats.
    with tempfile.TemporaryDirectory() as temp_dir:
        tree_inst, trees = load_trees(fasta_path, posterior_newick_path, temp_dir)
        pps = load_pps(pp_csv)
        tree_id_map, tree_pp_map, tree_cred_map, tree_found_map = build_tree_dicts(
            trees, pps
        )

        all_seen_path = f"{seen_topology_dir}/topologies-seen.{total_seen_count}.nwk"
        restrict_to_dag(
            all_seen_path, (tree_id_map, tree_pp_map, tree_cred_map, tree_found_map)
        )
        _, pcsp_cred_map, pcsp_found_map = load_pcsp_pp_map(pcsp_pp_csv)
        subsplit_cred_map, subsplit_found_map = load_subsplit_map(subsplit_csv)

        # Construct sdags from MCMC trees and get stats.
        seen_tree_inst, seen_trees = load_trees(fasta_path, all_seen_path, temp_dir)
        sdag_results_df = sdag_results_df_of(
            seen_trees,
            total_seen_count,
            all_seen_path,
            tree_id_map,
            tree_pp_map,
            tree_cred_map,
            tree_found_map,
            pcsp_cred_map,
            pcsp_found_map,
            subsplit_cred_map,
            subsplit_found_map,
        )

    final_df = accumulation_df.merge(sdag_results_df, on="support_size")

    final_df.to_csv(out_path)
    return None


if __name__ == "__main__":
    run()
