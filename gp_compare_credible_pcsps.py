import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import bito

fasta_path = "data/ds1.fasta"
full_path = "data/ds1.mb-trees.rerooted.nwk"
cred_path = "data/ds1.credible.rerooted.nwk"
outdata_path = "_output/gp_cred_vs_noncred_results.csv"

# build inst from mb-trees in the credible set
print("# build cred inst...")
cred_inst = bito.gp_instance("_ignore/cred-mmap.data")
cred_inst.read_fasta_file(fasta_path)
cred_inst.read_newick_file(cred_path)
cred_inst.make_dag()
cred_dag = cred_inst.get_dag()
# cred_inst.make_gp_engine()
# cred_inst.take_first_branch_length()
# cred_inst.populate_plvs()
# cred_inst.compute_likelihoods()
# cred_llhs = cred_inst.get_per_pcsp_log_likelihoods()

# build inst from all mb-trees
print("# build full inst...")
full_inst = bito.gp_instance("_ignore/full-mmap.data")
full_inst.read_fasta_file(fasta_path)
full_inst.read_newick_file(full_path)
full_inst.make_dag()
full_dag = full_inst.get_dag()
full_inst.make_gp_engine()
full_inst.take_first_branch_length()
full_inst.populate_plvs()
full_inst.compute_likelihoods()
full_llhs = full_inst.get_per_pcsp_log_likelihoods()


# build map from pcsp to edge_id.
full_pcsps = full_dag.build_sorted_vector_of_edge_bitsets()
full_edge_map = {}
for pcsp in full_pcsps:
    edge_id = full_dag.get_edge_id(pcsp)
    full_edge_map[edge_id.value()] = pcsp

# build map of llhs
full_pcsp_llh_map = {}
for edge_id in full_edge_map:
    full_pcsp_llh_map[edge_id] = full_llhs[edge_id]

# get total likelihood
total_llh = 0
for llh in full_llhs:
    total_llh += llh

# build map of probabilities
full_pcsp_prob_map = {}
for edge_id in full_edge_map:
    full_pcsp_prob_map[edge_id] = full_llhs[edge_id] / total_llh

# find which edges exist in the credible
cred_pcsps = set(cred_dag.build_sorted_vector_of_edge_bitsets())
is_pcsp_cred_map = {}
noncred_pcsps = []
for edge_id in full_edge_map:
    pcsp = full_edge_map[edge_id]
    if pcsp in cred_pcsps:
        is_pcsp_cred_map[edge_id] = True
    else:
        is_pcsp_cred_map[edge_id] = False

# build map of smaller clade size counts.
smaller_clade_size_map = {}
for edge_id in full_edge_map:
    pcsp = full_edge_map[edge_id]
    subsplit = pcsp.pcsp_get_parent_subsplit()
    smaller_clade_size = min(subsplit.subsplit_get_clade(0).clade_get_count(),
                             subsplit.subsplit_get_clade(1).clade_get_count())
    smaller_clade_size_map[edge_id] = smaller_clade_size

print("# build data...")
pcsp_data = []
for edge_id in full_edge_map:
    pcsp = full_edge_map[edge_id].pcsp_to_string()
    llh = full_pcsp_llh_map[edge_id]
    cred = is_pcsp_cred_map[edge_id]
    smaller_clade_size = smaller_clade_size_map[edge_id]
    pcsp_data.append([pcsp, llh, cred, smaller_clade_size])
# build data frame
df = pd.DataFrame(pcsp_data, columns=[
                  "pcsp", "llh", "cred", "smaller_clade_size"])

df.to_csv(outdata_path)

print("# ...done.")
