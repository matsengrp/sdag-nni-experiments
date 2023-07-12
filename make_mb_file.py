import re
import click


@click.command()
@click.argument("ds", type=str)
@click.argument("initial_tree_path", type=str)
@click.argument("branch_length_prior", type=str)
@click.argument("out_path", type=str)
def prepare_mb_run_file(ds, initial_tree_path, branch_length_prior, out_path):
    """
    Use the mb template file and the specified settings to create a Mr Baye file that is
    ready to run.

    Parameters:
        ds (string): The dataset number (without the prefix "ds").
        initial_tree_path (string): The file path for the newick string, with branch
            lengths, of the starting tree.
        branch_length_prior (string): The prior on the branch lengths, in the format of
            Mr. Bayes.
    """
    template_path = "data/mb_template.mb"
    with open(template_path) as mb_template_file:
        file_lines = mb_template_file.readlines()
    with open(initial_tree_path) as initial_tree_file:
        initial_tree_newick = initial_tree_file.readline().strip()

    replacement_dict = {
        "{{ds}}": ds,
        "{{starttree}}": initial_tree_newick,
        "{{brlenspr}}": branch_length_prior,
    }
    replacement_match = "(" + ")|(".join(replacement_dict.keys()) + ")"
    replacement_sub = lambda match: replacement_dict[match.group()]

    formatted_file_lines = [
        re.sub(replacement_match, replacement_sub, line) for line in file_lines
    ]

    print(f"Trying to write to {out_path}")
    with open(out_path, "w") as mb_file:
        for line in formatted_file_lines:
            mb_file.write(line)


if __name__ == "__main__":
    prepare_mb_run_file()
