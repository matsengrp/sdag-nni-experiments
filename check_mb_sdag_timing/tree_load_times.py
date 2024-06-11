import scipy.stats as ss
import sys

sys.path.append("../")
import bito
from nni_search import load_trees
import tempfile
import time


def load_times():
    for ds in [1, 3, 4, 5, 6, 7, 8]:
        fasta_path = f"../data/ds{ds}/ds{ds}.fasta"
        for line_count in [1000, 10000, 100000]:
            newick_path = f"ds{ds}/{line_count}_lines.results.csv"
            with tempfile.TemporaryDirectory() as temp_dir:
                start = time.time()
                tree_inst, trees = load_trees(fasta_path, newick_path, temp_dir)
                runtime = time.time() - start
                tree_count = len(trees)
                print(f"Loading {tree_count} trees from ds{ds} took {runtime} seconds.")


if __name__ == "__main__":
    load_times()
