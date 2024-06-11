# This script demonstrates a fundamental design flaw in the variant of top-pruning where
# we add (at most) 5 edges per NNI and allow only the two NNIs at an edge using the
# choice map tree through the edge: there can be unobtainable edges and subsplits.

# I messed up the algorithm described in bito issue #440 at step 7. We cannot add the
# NNIs at the sister and parent edges of u->t', because those edges alreeady exist and
# their choice maps are set. There may be other mistakes...


class sdag:
    """A simple class for the topological structure and choice maps of an sDAG."""

    def __init__(self, root_split):
        self.root_split = root_split
        self.edges = []
        self.subsplits = [root_split]
        self.taxon_set = root_split.left_clade.union(root_split.right_clade)
        self.parent_choice = {}
        self.sibling_choice = {}
        self.left_child_choice = {}
        self.right_child_choice = {}

    def add_subsplit(self, the_subsplit):
        """
        Adds the subsplit to this sdag, if the subsplit is not already present. Returns
        the truth value of the subsplit being new to the sdag.
        """
        if self.subsplit_exists(the_subsplit):
            return False
        else:
            self.subsplits.append(the_subsplit)
            return True

    def add_edge(self, the_edge):
        """
        Adds the edge to this sdag, if the edge is not already present. Returns the
        truth value of the edge being new to the sdag.
        """
        if self.edge_exists(the_edge):
            return False
        else:
            self.edges.append(the_edge)
            return True

    def subsplit_exists(self, the_subsplit):
        """Is the subsplit in this sdag?"""
        return the_subsplit in self.subsplits

    def edge_exists(self, the_edge):
        """Is the edge in this sdag?"""
        return the_edge in self.edges

    def valid_nni_edge(self, the_edge):
        """
        Is the the edge a valid edge to apply an NNI (i.e., neither leaf nor root)?
        """
        return not (the_edge.is_leaf_edge() or the_edge.is_root_edge(self))

    def set_choices_to_first_available(self):
        """
        Initialize the choice maps for the edges of this sdag, using the first viable
        edge in self.edges.
        """
        for e in self.edges:
            if e.is_root_edge(self):
                self.parent_choice[e] = None
            else:
                parents = [parent for parent in self.edges if e.valid_parent(parent)]
                self.parent_choice[e] = parents[0]
            siblings = [sibling for sibling in self.edges if e.valid_sibling(sibling)]
            self.sibling_choice[e] = siblings[0]
            if e.is_leaf_edge():
                self.left_child_choice[e] = None
                self.right_child_choice[e] = None
            else:
                lefts = [child for child in self.edges if e.valid_left_child(child)]
                rights = [child for child in self.edges if e.valid_right_child(child)]
                self.left_child_choice[e] = lefts[0]
                self.right_child_choice[e] = rights[0]
        return None

    def NNI_child_edge(self, the_edge, swap_right):
        """
        Performs the NNI at the_edge of the choice map topology through the edge,
        swapping the sibling edge with either the right or left child edge. The
        resulting subsplits and edges are added to this sdag, if they do not already
        exist. For edges that do not already exist, we set their choice maps on this
        sdag.

        Returned is a tuple of tuples, containing the 5 edges and truth values for
        them being new.
        """
        if not self.edge_exists(the_edge):
            raise ValueError("The given edge is not in the sdag.")
        if the_edge.is_leaf_edge():
            raise ValueError("The given edge ends with a leaf node.")
        if the_edge.is_root_edge(self):
            raise ValueError("The given edge begins with a rootsplit")

        swap_left = not swap_right
        old_central = the_edge
        old_parent = self.parent_choice[old_central]
        old_sibling = self.sibling_choice[old_central]
        old_left_child = self.left_child_choice[old_central]
        old_right_child = self.right_child_choice[old_central]
        old_parent_subsplit = old_central.parent_subsplit
        old_child_subsplit = old_central.child_subsplit

        X, Y = old_child_subsplit.left_clade, old_child_subsplit.right_clade
        # Determine if the old parent subsplit is XY|Z or Z|XY.
        old_child_subsplit_is_on_left = old_central.is_child_subsplit_on_left()
        old_child_subsplit_is_on_right = not old_child_subsplit_is_on_left
        if old_child_subsplit_is_on_left:
            Z = old_parent_subsplit.right_clade
        else:
            Z = old_parent_subsplit.left_clade

        if swap_right:
            parent_subsplit = subsplit(X.union(Z), Y)
            child_subsplit = subsplit(X, Z)
        else:
            # Determine if the new child subsplit is Y|Z or Z|Y
            old_sibling_clade_goes_on_right = min(Y) < min(Z)
            old_sibling_clade_goes_on_left = not old_sibling_clade_goes_on_right
            parent_subsplit = subsplit(X, Y.union(Z))
            child_subsplit = subsplit(Y, Z)
        self.add_subsplit(parent_subsplit)
        self.add_subsplit(child_subsplit)

        parent = edge(old_parent.parent_subsplit, parent_subsplit)
        central = edge(parent_subsplit, child_subsplit)
        if swap_right:
            sibling = edge(parent_subsplit, old_right_child.child_subsplit)
            if old_child_subsplit_is_on_left:
                left_child = edge(child_subsplit, old_left_child.child_subsplit)
                right_child = edge(child_subsplit, old_sibling.child_subsplit)
            else:
                left_child = edge(child_subsplit, old_sibling.child_subsplit)
                right_child = edge(child_subsplit, old_left_child.child_subsplit)
        else:
            sibling = edge(parent_subsplit, old_left_child.child_subsplit)
            if old_sibling_clade_goes_on_left:
                left_child = edge(child_subsplit, old_sibling.child_subsplit)
                right_child = edge(child_subsplit, old_right_child.child_subsplit)
            else:
                left_child = edge(child_subsplit, old_right_child.child_subsplit)
                right_child = edge(child_subsplit, old_sibling.child_subsplit)
        parent_is_new = self.add_edge(parent)
        central_is_new = self.add_edge(central)
        sibling_is_new = self.add_edge(sibling)
        left_is_new = self.add_edge(left_child)
        right_is_new = self.add_edge(right_child)

        if parent_is_new:
            self.parent_choice[parent] = self.parent_choice[old_parent]
            self.sibling_choice[parent] = self.sibling_choice[old_parent]
            if swap_right or old_child_subsplit_is_on_right:
                self.left_child_choice[parent] = central
                self.right_child_choice[parent] = sibling
            else:
                self.left_child_choice[parent] = sibling
                self.right_child_choice[parent] = central
        if central_is_new:
            self.parent_choice[central] = parent
            self.sibling_choice[central] = sibling
            self.left_child_choice[central] = left_child
            self.right_child_choice[central] = right_child
        if sibling_is_new:
            self.parent_choice[sibling] = parent
            self.sibling_choice[sibling] = central
            pre_nni_edge = old_right_child if swap_right else old_left_child
            self.left_child_choice[sibling] = self.left_child_choice[pre_nni_edge]
            self.right_child_choice[sibling] = self.right_child_choice[pre_nni_edge]
        if left_is_new:
            self.parent_choice[left_child] = central
            self.sibling_choice[left_child] = right_child
            if swap_right and old_child_subsplit_is_on_left:
                pre_nni_edge = old_left_child
            elif swap_left and old_sibling_clade_goes_on_right:
                pre_nni_edge = old_right_child
            else:
                pre_nni_edge = old_sibling
            self.left_child_choice[left_child] = self.left_child_choice[pre_nni_edge]
            self.right_child_choice[left_child] = self.right_child_choice[pre_nni_edge]
        if right_is_new:
            self.parent_choice[right_child] = central
            self.sibling_choice[right_child] = left_child
            if swap_right and old_child_subsplit_is_on_right:
                pre_nni_edge = old_left_child
            elif swap_left and old_sibling_clade_goes_on_left:
                pre_nni_edge = old_right_child
            else:
                pre_nni_edge = old_sibling
            self.left_child_choice[right_child] = self.left_child_choice[pre_nni_edge]
            self.right_child_choice[right_child] = self.right_child_choice[pre_nni_edge]

        return (
            (parent, central, sibling, left_child, right_child),
            (parent_is_new, central_is_new, sibling_is_new, left_is_new, right_is_new),
        )


class edge:
    """A simple class for sdag edges as pairs of subsplits."""

    def __init__(self, parent_subsplit, child_subsplit):
        self.parent_subsplit = parent_subsplit
        self.child_subsplit = child_subsplit
        self.hash_value = hash(str(self))

    def __eq__(self, other):
        """
        Equality of edges is in terms of pairs of subsplits.
        """
        return (
            self.parent_subsplit == other.parent_subsplit
            and self.child_subsplit == other.child_subsplit
        )

    def __str__(self):
        return str(self.parent_subsplit) + " ---> " + str(self.child_subsplit)

    def __hash__(self):
        return self.hash_value

    def is_leaf_edge(self):
        """Is the child subsplit of this edge a leaf node?"""
        return self.child_subsplit.is_leaf_subsplit()

    def is_root_edge(self, the_sdag):
        """Is the parent subsplit of this edge a rootsplit for the taxa of the_sdag?"""
        return self.parent_subsplit.is_root_split(the_sdag)

    def is_child_subsplit_on_left(self):
        """
        Is the child subsplit of this edge bipartitioning the left subsplit clade of
        the parent subsplit?
        """
        return self.parent_subsplit.left_clade == self.child_subsplit.taxa

    def is_child_subsplit_on_right(self):
        """
        Is the child subsplit of this edge bipartitioning the right subsplit clade of
        the parent subsplit?
        """
        return self.parent_subsplit.right_clade == self.child_subsplit.taxa

    def valid_parent(self, other):
        """Is the other edge a valid choice for parent edge of self?"""
        return other.valid_child(self)

    def valid_sibling(self, other):
        """Is the other edge a valid choice for sibling edge of self?"""
        return (
            self.parent_subsplit == other.parent_subsplit
            and self.is_child_subsplit_on_left() == other.is_child_subsplit_on_right()
        )

    def valid_left_child(self, other):
        """Is the other edge a valid choice for left child edge of self?"""
        return (
            self.valid_child(other)
            and self.child_subsplit.left_clade == other.child_subsplit.taxa
        )

    def valid_right_child(self, other):
        """Is the other edge a valid choice for right child edge of self?"""
        return (
            self.valid_child(other)
            and self.child_subsplit.right_clade == other.child_subsplit.taxa
        )

    def valid_child(self, other):
        """Is the other edge a valid choice for child edge of self?"""
        return self.child_subsplit == other.parent_subsplit


class subsplit:
    """A simple class for subsplits as pairs of sets of taxa, with an enforced order."""

    def __init__(self, clade1, clade2):
        """Create a subsplit from clade1 and clade2, which should be sets of integers."""
        if (len(clade1) + len(clade2) == 0) or len(clade1.intersection(clade2)) != 0:
            raise ValueError("The clades do not form a subsplit.")
        c1_first = len(clade2) == 0 or (len(clade1) != 0 and min(clade1) < min(clade2))
        self.left_clade = clade1 if c1_first else clade2
        self.right_clade = clade2 if c1_first else clade1
        self.taxa = clade1.union(clade2)

    def __eq__(self, other):
        """Equality of subsplits is equality of both clades."""
        return (
            self.left_clade == other.left_clade
            and self.right_clade == other.right_clade
        )

    def __str__(self):
        clade_string = lambda clade: ",".join(map(str, sorted(clade)))
        return f"[{clade_string(self.left_clade)}|{clade_string(self.right_clade)}]"

    def is_root_split(self, the_sdag):
        """Is this subsplit a rootsplit for the taxa of the_sdag?"""
        return len(self.taxa) == len(the_sdag.taxon_set)

    def is_leaf_subsplit(self):
        """Is this subsplit a leaf node?"""
        return len(self.left_clade) == 1 and len(self.right_clade) == 0


def build_caterpillar(N):
    """
    Returns the sdag constructed from the topology where 0 splits off first, then 1
    splits off, then 2 splits off, ..., and lastly N-2 and N-1 split.
    """
    parent_subsplit = subsplit({0}, set(range(1, N)))
    the_sdag = sdag(parent_subsplit)
    for i in range(N - 1):
        left_subsplit = subsplit({i}, set())
        right_subsplit = subsplit({i + 1}, set(range(i + 2, N)))
        left_edge = edge(parent_subsplit, left_subsplit)
        right_edge = edge(parent_subsplit, right_subsplit)
        the_sdag.add_subsplit(left_subsplit)
        the_sdag.add_subsplit(right_subsplit)
        the_sdag.add_edge(left_edge)
        the_sdag.add_edge(right_edge)
        parent_subsplit = right_subsplit
    return the_sdag


def check_sdag_from_caterpillar(N):
    """
    This method constructs an sDAG from a caterpillar tree with N taxa, taking the taxon
    0 as the outgroup for the rootsplit. NNIs are then repeatedly applied to this sDAG,
    following the "just 5 edges" version of top-pruning. Unlike the actual top-pruning
    algorithm, we do not rank and add NNIs by likelihood. Instead we add them in order
    of a breadth first search. However, like top-pruning we allow only two NNIs per
    sDAG edge and the NNIs are based on the choice maps at the edges. We check for
    missing edges in the sDAG after applying all allowed NNIs.
    """
    the_sdag = build_caterpillar(N)
    the_sdag.set_choices_to_first_available()
    edges_to_nni = [e for e in the_sdag.edges if the_sdag.valid_nni_edge(e)]
    while len(edges_to_nni) != 0:
        current_edge = edges_to_nni.pop(0)
        for swap_right in (False, True):
            edges, are_edges_new = the_sdag.NNI_child_edge(current_edge, swap_right)
            for e, is_new in zip(edges, are_edges_new):
                if is_new and the_sdag.valid_nni_edge(e):
                    edges_to_nni.append(e)

    node_count = len(the_sdag.subsplits)
    missing_node_count = max_node_count(N - 1) + 1 - node_count
    edge_count = len(the_sdag.edges)
    missing_edge_count = max_edge_count(N - 1) + 1 - edge_count

    print(
        f"The sdag has {node_count} nodes, including the rootsplit and leaf 0, missing "
        + f"{missing_node_count}.\nThe sdag has {edge_count} edges, including those "
        + f"from the rootsplit, missing {missing_edge_count}."
    )
    if missing_edge_count != 0:
        print("A few of the missing edges are:")
        find_missing_edges(the_sdag)
    print()


def find_missing_edges(the_sdag, how_many=10):
    """
    Prints the edges that are not in the_sdag, but are valid edges between subsplits in
    the_sdag.
    """
    count = 0
    for parent in the_sdag.subsplits:
        if not parent.is_leaf_subsplit():
            for child in the_sdag.subsplits:
                if child.taxa == parent.left_clade or child.taxa == parent.right_clade:
                    if not the_sdag.edge_exists(edge(parent, child)):
                        print(f"The sdag is missing the edge from {parent} to {child}.")
                        count += 1
                        if count == how_many:
                            return None
    return None


def max_edge_count(N):
    """
    Returns the maximum number of edges in an sDAG with N taxa, including edges from the
    universal ancestor and edges to leaf nodes (bito style).
    """
    return int(2 ** (N - 1) * (2**N + N + 4) - 3 * (3**N + 1) / 2 - N)


def max_node_count(N):
    """
    Returns the maximum number of subsplits in an sDAG with N taxa, including the
    universal ancestor and leaf nodes (bito style).
    """
    return int(N + (3**N + 1) / 2 - 2**N + 1)


if __name__ == "__main__":
    # Missing edges start at N=7, missing subsplits start at N=9.
    for N in range(3, 10):
        check_sdag_from_caterpillar(N)
