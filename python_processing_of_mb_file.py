import click
import subprocess
import time


@click.command()
@click.argument("mb_t_file_path", type=str)
@click.argument("root", type=str)
@click.argument("out_path", type=str)
def process_mb_t_file(mb_t_file_path, root, out_path):
    """
    Performs the same functionality as the bash script to prepare the
    rerooted-topology-sequence.tab from the Mr. Bayes .t file, but through python to
    handle extra large files.
    """
    start = time.time()
    previous_topology = ""
    previous_topology_count = 0
    with open(mb_t_file_path) as the_file:
        for line in the_file:
            line = line.strip()
            if line.startswith("tree"):
                topology = process_topology(line[line.index("(") :], root)
                if topology == previous_topology:
                    previous_topology_count += 1
                else:
                    if previous_topology_count != 0:
                        with open(out_path, "a") as the_out_file:
                            the_out_file.write(
                                f"{previous_topology_count}\t{previous_topology}"
                            )
                    previous_topology = topology
                    previous_topology_count = 1

        with open(out_path, "a") as the_out_file:
            the_out_file.write(f"{previous_topology_count}\t{previous_topology}")
    stop = time.time()
    print(f"Run-time: {stop-start} seconds")


def process_topology(newick, root):
    command = f"echo '{newick}' | nw_topology - | nw_reroot - {root} | nw_order -"
    return subprocess.check_output(command, text=True, shell=True)


if __name__ == "__main__":
    process_mb_t_file()
