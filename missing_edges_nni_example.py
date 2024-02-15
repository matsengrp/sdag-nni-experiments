import bito
import tempfile
import os
import subprocess
from ete3 import Tree
import unittest

# This script verifies an example of an sDAG that does not have edges between all
# compatible subsplits and enlarging the sDAG by an NNI yields...
# 1) A topology that is not in the sDAG and is not an NNI of any topology in the sDAG.
# 2) A topology that is not in the sDAG and does not use the central edge of the NNI.


class test_utils_methods(unittest.TestCase):
    def setUp(self):
        self.make_paths()
        self.make_newicks()
        self.write_data_files()
        self.write_nni_neighbors()
        self.make_sdags()
        self.make_subsplits_and_edges()
        self.make_rooted_trees()
        return None

    def make_paths(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.dir = self.temp_dir.name
        self.fasta_path = os.path.join(self.dir, "0")
        self.pre_nni_sdag_gen_path = os.path.join(self.dir, "1")
        self.post_nni_sdag_gen_path = os.path.join(self.dir, "2")
        self.one_nni_nwk_path = os.path.join(self.dir, "3")
        self.missed_by_nni_nwk_path = os.path.join(self.dir, "4")
        self.no_central_edge_nwk_path = os.path.join(self.dir, "5")
        self.all_nni_nwk_path = os.path.join(self.dir, "6")
        self.dummy0 = os.path.join(self.dir, "-1")
        self.dummy1 = os.path.join(self.dir, "-2")
        self.dummy2 = os.path.join(self.dir, "-3")
        self.dummy3 = os.path.join(self.dir, "-4")
        self.dummy4 = os.path.join(self.dir, "-5")
        return None

    def make_newicks(self):
        self.tau0 = "(0,(((((1,(2,3)),(4,5)),(6,(7,8))),9),10));"
        self.tau1 = "(0,(((1,((2,6),((3,7),8))),(4,5)),(9,10)));"
        self.tau2 = "(0,(((((((1,2),3),((6,7),8)),4),5),10),9));"
        self.tau3 = "(0,(((((1,(2,3)),(6,(7,8))),(4,5)),9),10));"
        self.tau4 = "(0,(((((1,2),3),((6,7),8)),(4,5)),(9,10)));"
        self.tau5 = "(0,((((1,((2,6),((3,7),8))),(4,5)),9),10));"
        self.pre_nni_generators = [self.tau0, self.tau1, self.tau2]
        self.post_nni_generators = [self.tau0, self.tau1, self.tau2, self.tau3]
        return None

    def write_data_files(self):
        with open(self.fasta_path, "w") as the_file:
            the_file.write("".join((f">{j}\nAGCT\n" for j in range(11))))
        with open(self.pre_nni_sdag_gen_path, "w") as the_file:
            for tau in self.pre_nni_generators:
                the_file.write(tau + "\n")
        with open(self.post_nni_sdag_gen_path, "w") as the_file:
            for tau in self.post_nni_generators:
                the_file.write(tau + "\n")
        with open(self.one_nni_nwk_path, "w") as the_file:
            the_file.write(self.tau3 + "\n")
        with open(self.missed_by_nni_nwk_path, "w") as the_file:
            the_file.write(self.tau4 + "\n")
        with open(self.no_central_edge_nwk_path, "w") as the_file:
            the_file.write(self.tau5 + "\n")
        return None

    def write_nni_neighbors(self):
        temp_file = tempfile.NamedTemporaryFile(dir=self.dir)
        temp = temp_file.name
        with open(self.pre_nni_sdag_gen_path) as the_file:
            for line in the_file:
                subprocess.check_call(
                    f"echo '{line[:-1]}' | spr_neighbors --nni >> {temp}",
                    shell=True,
                )
        subprocess.call(
            f"cat {temp} | nw_reroot - 0 | nw_order - | sort -R | uniq > "
            + f"{self.all_nni_nwk_path}",
            shell=True,
        )
        temp_file.close()
        return None

    def make_sdags(self):
        self.pre_nni_sdag, self.pre_nni_sdag_inst = self.dag_from_path(
            self.dummy0, self.pre_nni_sdag_gen_path
        )
        self.post_nni_sdag, self.post_nni_sdag_inst = self.dag_from_path(
            self.dummy1, self.post_nni_sdag_gen_path
        )
        return None

    def make_subsplits_and_edges(self):
        # The order of taxa in bitstrings (left to right) is 0,1,10,2,3,4,5,6,7,8,9.
        # NOT 0,1,2,3,4,5,6,7,8,9,10.
        self.t_subsplit = bito.subsplit("01011110000", "00000001110")
        self.s_subsplit = bito.subsplit("01011000000", "00000110000")
        self.t_to_s = bito.pcsp(self.t_subsplit, self.s_subsplit)
        self.tp_subsplit = bito.subsplit("01011001110", "00000110000")
        self.sp_subsplit = bito.subsplit("01011000000", "00000001110")
        self.tp_to_sp = bito.pcsp(self.tp_subsplit, self.sp_subsplit)
        return None

    def make_rooted_trees(self):
        _, dag_inst = self.dag_from_path(self.dummy2, self.one_nni_nwk_path)
        self.rooted_tree3 = dag_inst.generate_complete_rooted_tree_collection().trees[0]
        _, dag_inst = self.dag_from_path(self.dummy3, self.missed_by_nni_nwk_path)
        self.rooted_tree4 = dag_inst.generate_complete_rooted_tree_collection().trees[0]
        _, dag_inst = self.dag_from_path(self.dummy4, self.no_central_edge_nwk_path)
        self.rooted_tree5 = dag_inst.generate_complete_rooted_tree_collection().trees[0]
        return None

    def dag_from_path(self, data_path, tree_path):
        dag_inst = bito.gp_instance(data_path)
        dag_inst.read_newick_file(tree_path)
        dag_inst.read_fasta_file(self.fasta_path)
        dag_inst.make_gp_engine()
        dag_inst.make_dag()
        dag = dag_inst.get_dag()
        return dag, dag_inst

    def tearDown(self):
        self.temp_dir.cleanup()
        return None

    def test_sdag_construction(self):
        # The pre-nni-sdag should consist of exactly three topologies. In particular,
        # the three choices are at the root and not next to the subplit t nor s, so that
        # applying the sDAG-NNI is the same as adding a single tree-nni. Also the
        # post-nni-sdag is larger only by edges, not subsplits.
        self.assertEqual(self.pre_nni_sdag.topology_count(), 3)
        post_nni_subsplits = self.post_nni_sdag.build_set_of_node_bitsets()
        post_nni_edges = self.post_nni_sdag.build_set_of_edge_bitsets()
        self.assertTrue(all(map(self.pre_nni_sdag.contains_node, post_nni_subsplits)))
        self.assertFalse(all(map(self.pre_nni_sdag.contains_edge, post_nni_edges)))
        self.assertEqual(
            self.post_nni_sdag.edge_count(), self.pre_nni_sdag.edge_count() + 4
        )
        nodes = (self.t_subsplit, self.s_subsplit, self.tp_subsplit, self.sp_subsplit)
        for sdag in (self.pre_nni_sdag, self.post_nni_sdag):
            for node in nodes:
                self.assertTrue(sdag.contains_node(node))
        self.assertTrue(self.pre_nni_sdag.contains_edge(self.t_to_s))
        self.assertFalse(self.pre_nni_sdag.contains_edge(self.tp_to_sp))
        self.assertTrue(self.post_nni_sdag.contains_edge(self.t_to_s))
        self.assertTrue(self.post_nni_sdag.contains_edge(self.tp_to_sp))
        return None

    def test_nnis(self):
        # The topology tau3 should an NNI of tau1, tau2, or tau3. The topology tau4
        # should not.
        with open(self.all_nni_nwk_path) as the_file:
            nni_newicks = list(map(str.strip, the_file))
        nni_ete_ids = [Tree(newick).get_topology_id() for newick in nni_newicks]
        self.assertIn(self.tau3, nni_newicks)
        self.assertIn(Tree(self.tau3).get_topology_id(), nni_ete_ids)
        self.assertNotIn(self.tau4, nni_newicks)
        self.assertNotIn(Tree(self.tau4).get_topology_id(), nni_ete_ids)
        return None

    def test_central_edge(self):
        # The edge t'->s' should be in tau3 and tau4, but not tau5.
        self.assertIn(self.tp_to_sp, self.rooted_tree3.build_set_of_pcsps())
        self.assertIn(self.tp_to_sp, self.rooted_tree4.build_set_of_pcsps())
        self.assertNotIn(self.tp_to_sp, self.rooted_tree5.build_set_of_pcsps())
        return None

    def test_sdag_membership(self):
        # All of tau3, tau4, tau5 should not be in the pre-nni-sdag, should be in the
        # post-nni-sdag, and should have all of their nodes in the pre-nni-sdag.
        for tree in (self.rooted_tree3, self.rooted_tree4, self.rooted_tree5):
            edges = tree.build_set_of_pcsps()
            subsplits = tree.build_set_of_subsplits()
            self.assertFalse(self.pre_nni_sdag.contains_tree(tree))
            self.assertFalse(all(map(self.pre_nni_sdag.contains_edge, edges)))
            self.assertTrue(all(map(self.pre_nni_sdag.contains_node, subsplits)))
            self.assertTrue(self.post_nni_sdag.contains_tree(tree))
        return None


if __name__ == "__main__":
    unittest.main()
