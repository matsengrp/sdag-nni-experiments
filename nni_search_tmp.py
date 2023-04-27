import sys
import argparse

############
### MAIN ###
############


def main_arg_parse(args):
    parser = argparse.ArgumentParser(
        description='Tools for performing NNI systematic search.')
    parser.add_argument('-v', '--verbose',
                        action='store_true', help='verbose output')
    parser.add_argument('-o', '--output', help='output file', type=str,
                        default='nni_search.results.txt')
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand')

    # nni search
    subparser1 = subparsers.add_parser(
        'nni_search', help='Perform iterative NNI search.')
    subparser1.add_argument('fasta', help='fasta file', type=str)
    subparser1.add_argument(
        'seed_newick', help='newick file for initial trees in DAG', type=str)
    subparser1.add_argument(
        'credible_newick', help='newick file for trees in credible posterior', type=str)
    subparser1.add_argument(
        'pp_csv', help='csv file containing the posterior weights of the trees from credible_trees', type=str)
    subparser1.add_argument(
        'pcsp_pp_csv', help='csv file containing the per-PCSP posterior weights', type=str)
    subparser1.add_argument(
        '--iter', help='number of NNI search iterations', type=int, default=10)
    group = subparser1.add_mutually_exclusive_group(required=True)
    group.add_argument('--gp', action='store_true',
                       help='Use generalized pruning.')
    group.add_argument('--tp', action='store_true',
                       help='Use top pruning via log likelihood.')

    # pcsp map builder
    subparser2 = subparsers.add_parser(
        'build_pcsp_map', help='Build per-PCSP map.')
    subparser2.add_argument('fasta', help='fasta file', type=str)
    subparser2.add_argument(
        'credible_newick', help='newick file for trees in credible posterior', type=str)
    subparser2.add_argument(
        'pp_csv', help='csv file containing the posterior weights of the trees from credible_trees', type=str)

    return parser.parse_args(args)


if __name__ == "__main__":
    print("# begin...")
    print("# begin args:", sys.argv)
    print("# begin args:", sys.argv[1:])
    args = main_arg_parse(sys.argv[1:])
    print("# end args:", args)
    print("# ...done")
