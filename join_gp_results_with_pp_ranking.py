import sys
import itertools
import pandas as pd
import bito

req_args = 3
if (len(sys.argv) != (req_args + 1)):
    print("Usage: <i:llh_scores> <i:pp_scores> <o:joined_results>")
    exit(0)

llh_filename = sys.argv[1]
pp_filename = sys.argv[2]
results_filename = sys.argv[3]

names1 = ["parent", "child", "llh"]
df1 = pd.read_csv(llh_filename, sep=" ",
                  names=names1, index_col=False)

print("# read in pp_ranking_by_tree...")
names2 = ["iter", "parent", "child", "tree_pp", "pcsp_pp"]
df2 = pd.read_csv(pp_filename, sep=" ", names=names2)

results_df = pd.merge(df1, df2, on=["parent", "child"])

results_df.to_csv(results_filename)
print("# ...done.")
