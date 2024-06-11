import bito
import tempfile
import os
from ete3 import Tree
from collections import defaultdict


def check_raxml(ds):
    path = f"{ds}/{ds}.fasta.raxml.mlTrees"
    with open(path) as the_file:
        trees = [Tree(line) for line in the_file]
    line_count = len(trees)
    input_topology_count = len({tree.get_topology_id() for tree in trees})
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_data_path = os.path.join(temp_dir, "mmap.dat")
        dag_inst = bito.gp_instance(temp_data_path)
        dag_inst.read_newick_file(path)
        dag_inst.make_dag()
        dag = dag_inst.get_dag()
        tree_count = int(dag.topology_count())
        node_count = dag.node_count()
        edge_count = dag.edge_count()
    print(
        f"tree searchs: {line_count}\ndistinct input topologies: {input_topology_count},\n"
        + f"sDAG: {tree_count} topologies, {node_count} nodes, {edge_count} edges."
    )


def write_unique(ds):
    path = f"{ds}/{ds}.fasta.raxml.mlTrees"
    with open(path) as the_file:
        trees = [Tree(line) for line in the_file]

    id_to_trees = defaultdict(list)
    for tree in trees:
        id_to_trees[tree.get_topology_id()].append(tree)

    print(f"There are {len(id_to_trees)} distinct topologies.")
    for the_id in id_to_trees:
        print(f"{the_id}: {len(id_to_trees[the_id])} trees")

    with open(f"{ds}/{ds}.raxml.unique.nwk", "w") as the_file:
        for the_id in id_to_trees:
            the_file.write(id_to_trees[the_id][0].write(format=5) + "\n")


if __name__ == "__main__":
    for ds in ("ds1", "ds3", "ds4", "ds5", "ds6", "ds7", "ds8", "flu100"):
        print(ds)
        check_raxml(ds)
        write_unique(ds)
        print("#" * 30)
