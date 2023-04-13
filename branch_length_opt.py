import sys
from pathlib import Path
import bito
import numpy as np
import pprint
import bito.beagle_flags as beagle_flags
import bito.phylo_model_mapkeys as model_keys

SIMPLE_SPECIFICATION = bito.PhyloModelSpecification(
    substitution="JC69", site="constant", clock="none"
)

pp = pprint.PrettyPrinter(indent=4)

if len(sys.argv) < 6:
    print(
        "Usage: <i:fasta_path> <i:newick_path> <i:temp_dir> <o:newick_with_branches> <o:likelihoods>"
    )
    exit()
fasta_path = Path(sys.argv[1])
newick_path = Path(sys.argv[2])
temp_dir = Path(sys.argv[3])
mmap_path = Path(temp_dir / "mmapped_pv.data")
newick_out_path = Path(sys.argv[4])
likelihood_out_path = Path(sys.argv[5])

if not fasta_path.exists():
    print("<i:fasta_path> does not exist.")
    exit()
if not newick_path.exists():
    print("<i:newick_path> does not exist.")
    exit()
if not temp_dir.exists():
    print("<i:temp_dir> does not exist.")
    exit()
if not newick_out_path.parents[0].exists():
    print("<o:newick_with_branches> containing folder does not exist.")
    exit()
if not likelihood_out_path.parents[0].exists():
    print("<o:likelihoods> containing folder does not exist.")
    exit()

temp_newicks = []

# write each line into a separate file.
with newick_path.open("r") as fp:
    for line in fp:
        temp_newick = (temp_dir / fasta_path.name).with_suffix(
            ".{}.nwk".format(len(temp_newicks))
        )
        if False:
            temp_fp = temp_newick.open("w")
            temp_fp.write(line.strip() + "\n")
            temp_fp.close()
        temp_newicks.append(temp_newick)
        print(len(temp_newicks), line)

print(temp_newicks)

fp_newick_out = newick_out_path.open("w")
fp_likelihood_out = likelihood_out_path.open("w")

# for each tree in newick file, optimize branch lengths.
min_likelihood = np.inf
step_size = 10
for i in range(0, len(temp_newicks), step_size):
    temp_newick = temp_newicks[i]
    print("writing out: ", i, temp_newick)
    inst = bito.gp_instance(str(temp_dir))
    inst.read_fasta_file(str(fasta_path))
    inst.read_newick_file(str(temp_newick))
    inst.make_dag()
    inst.make_gp_engine()
    # inst.estimate_branch_lengths(tol=1e-4, max_iter=100, quiet=True)
    # dag = inst.get_dag()
    # branch_lengths = np.array(inst.get_branch_lengths())
    # tree_collection = inst.generate_complete_rooted_tree_collection()
    # newick = tree_collection.newick()
    inst.populate_plvs()
    inst.compute_likelihoods()
    log_likelihoods = np.array(inst.get_per_pcsp_log_likelihoods())
    likelihood = log_likelihoods[0]
    if min_likelihood > likelihood:
        min_likelihood = likelihood
    print("log_likelihood: ", i, likelihood)

    # fp_newick_out.write(newick.strip() + "\n")
    fp_likelihood_out.write(str(likelihood) + "\n")

fp_newick_out.close()
fp_likelihood_out.close()

print("min_likelihood:", min_likelihood)
