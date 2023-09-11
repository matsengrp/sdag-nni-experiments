import pandas as pd


def write_edge_counts_to_file(out_path):
    """
    Get the number of edges from credible and posterior trees based on the csv files for
    PCSPs (prepared by running nni-search with build-pcsp-map) and write to out_path.
    """
    the_data = [
        (ds, prior, *get_posterior_edge_counts(ds, prior))
        for ds in [1, 3, 4, 5, 6, 7, 8]
        for prior in ["uniform", "exponential"]
    ]
    header = "ds,prior,file_edge_count,posterior_edge_count,credible_edge_count\n"
    with open(out_path, "w") as the_file:
        the_file.write(header)
        for line in the_data:
            the_file.write(",".join(map(str, line)) + "\n")
    return None


def get_posterior_edge_counts(ds, prior="uniform"):
    prior_path = {"uniform": "_output", "exponential": "_exp_output"}
    pcsp_pp_path = f"data/ds{ds}/{prior_path[prior]}/ds{ds}.pcsp-pp.csv"
    edge_df = pd.read_csv(pcsp_pp_path)
    edge_count = len(edge_df)
    posterior_edge_count = sum(edge_df.pcsp_pp > 0)
    credible_edge_count = sum(edge_df.in_cred_set)
    return edge_count, posterior_edge_count, credible_edge_count


if __name__ == "__main__":
    write_edge_counts_to_file("data/posterior_edge_counts.csv")
