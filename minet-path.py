#!/usr/bin/env python3

# Import modules
import networkx as nx
import MineClient3 as mc
import sys
import argparse
import pickle
import multiprocessing as mp
from multiprocessing.managers import BaseManager
from collections import defaultdict as dd
from itertools import repeat
import time

# Define functions
def Chunks(lst,n):
    return [ lst[i::n] for i in range(n) ]


def CountRxns(network):
    """Count the number of unique reactions in a network."""
    rxns = set()
    for node in network.nodes(data=True):
        if node[1]['type'] != 'c':
            rxns.add(node[1]['mid'])
    return len(rxns)

def test_CountRxns():
    p1 = nx.DiGraph()
    p1.add_path(range(1,5))
    for node_data in enumerate([('c','C1'), ('rf','R1'), ('pf','R1'), ('c','C2')]):
        n = node_data[0] + 1
        p1.node[n]['type'] = node_data[1][0]
        p1.node[n]['mid'] = node_data[1][1]

    p2 = nx.DiGraph()
    p2.add_path(range(1,5))
    for node_data in enumerate([('c','C3'), ('rr','R2'), ('pr','R2'), ('c','C4')]):
        n = node_data[0] + 1
        p2.node[n]['type'] = node_data[1][0]
        p2.node[n]['mid'] = node_data[1][1]

    p3 = nx.DiGraph()
    p3.add_node(1, type='c', mid='C5')

    p4 = [('c','C6'), ('rf','R3'), ('pf','R3'), ('c','C7'), ('rr','R3'), ('rp','R3'), ('c','C8')]
    p4 = nx.DiGraph()
    p4.add_path(range(1,8))
    for node_data in enumerate([('c','C6'), ('rf','R3'), ('pf','R3'), ('c','C7'), ('rr','R3'), ('rp','R3'), ('c','C8')]):
        n = node_data[0] + 1
        p4.node[n]['type'] = node_data[1][0]
        p4.node[n]['mid'] = node_data[1][1]

    assert CountRxns(p1) == 1
    assert CountRxns(p2) == 1
    assert CountRxns(p3) == 0
    assert CountRxns(p4) == 1


def FindStartCompNodes(network):
    """Returns a list starting compound nodes in a MINE network."""
    start_comp_nodes = []
    for node in network.nodes():
        if network.node[node]['type'] == 'c' and network.node[node]['start']:
            start_comp_nodes.append(node)
    return set(start_comp_nodes)

def test_FindStartCompNodes():
    G = nx.DiGraph()
    G.add_node(1,type='c',start=True)
    G.add_node(2,type='rf')
    G.add_node(3,type='pf')
    G.add_node(4,type='c',start=False)
    G.add_path([1,2,3,4])
    assert FindStartCompNodes(G) == set([1])


def FindValidReactantNodes(network, proc_num=1, comp_node_set=set(), force_parallel=False):
    """Find and return all reactant nodes that are valid given the supplied compound node set."""
    # If the set is empty, use starting compounds as the compound node set
    if comp_node_set == set():
        comp_node_set = FindStartCompNodes(network)

    # Define the Worker
    def Worker(work):
        results = set()
        for node in work:
            if network.node[node]['type'] in {'rf','rr'}:
                c = network.node[node]['c']
                if c.issubset(comp_node_set):
                    results.add(node)
        output.extend(list(results))

    # Only go parallel if there are 5k or more items
    if len(network.nodes()) < 5000 and not force_parallel:
        output = []
        Worker(network.nodes())
        valid_reactant_nodes = set(output)

    else:
        with mp.Manager() as manager:
            # Initialize output list in manager
            output = manager.list()

            # Initialize processes
            procs = []
            for work in Chunks(network.nodes(), proc_num):
                p = mp.Process(target=Worker, args=(work,))
                procs.append(p)
                p.start()

            # Stop workers
            for p in procs:
                p.join()

            # Get the results
            valid_reactant_nodes = set(output)

    return valid_reactant_nodes

def test_FindValidReactantNodes():
    G = nx.DiGraph()
    G.add_node(1,type='c',start=True)
    G.add_node(2,type='rf',c={1})
    G.add_node(3,type='pf',c={4})
    G.add_node(4,type='c',start=False)
    G.add_path([1,2,3,4])
    G.add_node(9,type='rr',c={4})
    G.add_node(10,type='pr',c={1})
    G.add_path([4,9,10,1])
    G.add_node(5,type='c',start=False)
    G.add_node(6,type='rf',c={1,5})
    G.add_node(7,type='pf',c={8})
    G.add_node(8,type='c',start=False)
    G.add_path([1,6,7,8])
    G.add_edge(5,6)
    G.add_node(11,type='rr',c={8})
    G.add_node(12,type='pr',c={1,5})
    G.add_path([8,11,12,1])
    G.add_edge(12,5)
    start_comps = FindStartCompNodes(G)

    assert FindValidReactantNodes(G, 4, start_comps) == set([2])
    assert FindValidReactantNodes(G, 4, start_comps, force_parallel=True) == set([2])

    assert FindValidReactantNodes(G, 4) == set([2])
    assert FindValidReactantNodes(G, 4, force_parallel=True) == set([2])

    assert FindValidReactantNodes(G, 4, set([1,5])) == set([2,6])
    assert FindValidReactantNodes(G, 4, set([1,5]), force_parallel=True) == set([2,6])

    assert FindValidReactantNodes(G, 4, set([8])) == set([11])
    assert FindValidReactantNodes(G, 4, set([8]), force_parallel=True) == set([11])


def ExpandValidCompoundSet(network, proc_num=1, valid_reactant_nodes=set(), comp_node_set=set(), force_parallel=False):
    """Expands the valid compound set with products of valid reactant nodes."""
    if comp_node_set == set():
        comp_node_set = FindStartCompNodes(network)
    else:
        # Define the Worker
        def Worker(work):
            compounds = set()
            for r_node in work:
                p_node = network.successors(r_node)[0]
                products = network.node[p_node]['c']
                compounds = compounds.union(products)
            output.extend(compounds)

        # Only go parallel if there are 5k or more items
        if len(valid_reactant_nodes) < 5000  and not force_parallel:
            output = []
            Worker(valid_reactant_nodes)
            comp_node_set = comp_node_set.union(set(output))

        else:
            with mp.Manager() as manager:
                # Initialize output list in manager
                output = manager.list()

                # Initialize processes
                procs = []
                for work in Chunks(list(valid_reactant_nodes), proc_num):
                    p = mp.Process(target=Worker, args=(work,))
                    procs.append(p)
                    p.start()

                # Stop workers
                for p in procs:
                    p.join()

                # Get results
                comp_node_set = comp_node_set.union(set(output))

    return comp_node_set

def test_ExpandValidCompoundSet():
    G = nx.DiGraph()
    G.add_node(1,type='c',start=True)
    G.add_node(2,type='rf',c={1})
    G.add_node(3,type='pf',c={4})
    G.add_node(4,type='c',start=False)
    G.add_path([1,2,3,4])
    G.add_node(9,type='rr',c={4})
    G.add_node(10,type='pr',c={1})
    G.add_path([4,9,10,1])
    G.add_node(5,type='c',start=False)
    G.add_node(6,type='rf',c={1,5})
    G.add_node(7,type='pf',c={8})
    G.add_node(8,type='c',start=False)
    G.add_path([1,6,7,8])
    G.add_edge(5,6)
    G.add_node(11,type='rr',c={8})
    G.add_node(12,type='pr',c={1,5})
    G.add_path([8,11,12,1])
    G.add_edge(12,5)

    assert ExpandValidCompoundSet(G) == set([1])
    assert ExpandValidCompoundSet(G, force_parallel=True) == set([1])

    assert ExpandValidCompoundSet(G, 4, set([2,9]), set([1,4])) == set([1,4])
    assert ExpandValidCompoundSet(G, 4, set([2,9]), set([1,4]), force_parallel=True) == set([1,4])

    assert ExpandValidCompoundSet(G, 4, FindValidReactantNodes(G, 4, set([8])), set([8])) == set([1,5,8])
    assert ExpandValidCompoundSet(G, 4, FindValidReactantNodes(G, 4, set([8]), force_parallel=True), set([8]), force_parallel=True) == set([1,5,8])

    assert ExpandValidCompoundSet(G, 4, set([2,6,11]), ExpandValidCompoundSet(G, 4, set([11]), set([8]))) == set([1,4,5,8])
    assert ExpandValidCompoundSet(G, 4, set([2,6,11]), ExpandValidCompoundSet(G, 4, set([11]), set([8]), force_parallel=True), force_parallel=True) == set([1,4,5,8])


def DistanceToOrigin(network, proc_num=1, N=-1):
    """
    Calculates the shortest distance (number of reactions) from the starting compounds (origin) to
    every node up to distance N. Set N to -1 to exhaustively calculate the minimum distance to
    every node that is reachable.

    Returns two sets in a tuple: Valid compound nodes and valid reactant nodes.
    """

    sys.stdout.write("\nCalculating minimum distance of nodes to origin...\n\n")
    sys.stdout.flush()

    time_start = time.time()

    # Number of nodes for formatting
    L = len(network.nodes())
    l = len(str(L))

    # Set up counters
    n = 0
    c = 0
    rf = 0
    pf = 0
    rr = 0
    pr = 0

    # Start with no valid reactant or compound nodes
    valid_reactant_nodes = set([])
    valid_compound_nodes = set([])

    # The "previous" lists are also empty
    prev_vrn = list(valid_reactant_nodes)
    prev_vcn = list(valid_compound_nodes)

    # Start with no new valid reactant nodes
    new_vrn = set([])

    while True:

        # Valid product nodes will be connected at the start of an expansion cycle
        # They are however in principle identified in the previous cycle via valid reactant nodes
        for r_node in new_vrn:
            p_node = network.successors(r_node)[0]
            network.node[p_node]['dist'] = n
            node_type = network.node[p_node]['type']
            if node_type == 'pf':
                pf += 1
            if node_type == 'pr':
                pr += 1

        # Expand the valid compound set
        # When n = 0, this means the starting compounds
        # When n > 0, the valid compound set will be expanded
        new_vcn = ExpandValidCompoundSet(network, proc_num, new_vrn, valid_compound_nodes) - valid_compound_nodes
        valid_compound_nodes = new_vcn.union(valid_compound_nodes)

        new_vrn = FindValidReactantNodes(network, proc_num, valid_compound_nodes) - valid_reactant_nodes
        valid_reactant_nodes = new_vrn.union(valid_reactant_nodes)

        for node in new_vcn:
            network.node[node]['dist'] = n
            c += 1
        for node in new_vrn:
            network.node[node]['dist'] = n
            node_type = network.node[node]['type']
            if node_type == 'rf':
                rf += 1
            if node_type == 'rr':
                rr += 1

        # Nicely (hopefully) formatted progress output
        output = '{0:<%s} {1:>%s} {2:>%s} {3:>%s} {4:>%s} {5:>%s}' % (str(l+6), str(l+5), str(l+5), str(l+5), str(l+5), str(l+5))
        print(output.format('Step ' + str(n) + ':', str(c) + ' c', str(rf) + ' rf', str(pf) + ' pf', str(rr) + ' rr', str(pr) + ' pr'))

        n += 1

        if set(prev_vrn) == valid_reactant_nodes and set(prev_vcn) == valid_compound_nodes:
            # When no new valid compound or reactant nodes have been identified, it is time to stop
            break
        else:
            if n > N and N !=-1:
                # n starts at 0 and increments by one before each round dealing
                # with that particular step n - stop when n exceeds the limit
                break
        prev_vrn = list(valid_reactant_nodes)
        prev_vcn = list(valid_compound_nodes)

    total_time = time.time() - time_start

    sys.stdout.write("\nDone in %ss.\n" %str(total_time))
    sys.stdout.flush()

    return (valid_compound_nodes, valid_reactant_nodes)

def test_DistanceToOrigin():
    G = nx.DiGraph()
    G.add_node(1,type='c',start=False)
    G.add_node(2,type='rf',c={1})
    G.add_node(3,type='pf',c={4})
    G.add_node(4,type='c',start=False)
    G.add_path([1,2,3,4])
    G.add_node(9,type='rr',c={4})
    G.add_node(10,type='pr',c={1})
    G.add_path([4,9,10,1])
    G.add_node(5,type='c',start=False)
    G.add_node(6,type='rf',c={1,5})
    G.add_node(7,type='pf',c={8})
    G.add_node(8,type='c',start=True) # Compound 8 is now the start
    G.add_path([1,6,7,8])
    G.add_edge(5,6)
    G.add_node(11,type='rr',c={8})
    G.add_node(12,type='pr',c={1,5})
    G.add_path([8,11,12,1])
    G.add_edge(12,5)

    output_0 = DistanceToOrigin(G.copy(), 4, 0) # Creating a copy, since the function modifies the network it is given
    output_1 = DistanceToOrigin(G.copy(), 4, 1)
    output_2 = DistanceToOrigin(G, 4, 2) # Not creating a copy in order to check modification capabilities

    assert output_0 == (set([8]), set([11])) # Compound and reactant nodes reachable within 0 reaction steps, respectively
    assert output_1 == (set([8,1,5]), set([11,6,2]))
    assert output_2 == (set([8,1,5,4]), set([11,6,2,9]))

    assert set([G.node[n]['dist'] for n in [8,11]]) == set([0])
    assert set([G.node[n]['dist'] for n in [1,2,5,6,12]]) == set([1])
    assert set([G.node[n]['dist'] for n in [3,4,7,9]]) == set([2])

    # Test for detection of non-reachable nodes
    # Only nodes 1, 2, 3, 4, 5 and 6 should be reachable and will receive a 'dist' value - other nodes do not
    Y = nx.DiGraph()
    Y.add_node(1,type='c',start=True)
    Y.add_node(2,type='rf',c={1})
    Y.add_node(3,type='pf',c={6})
    Y.add_node(6,type='c',start=False)
    Y.add_path([1,2,3,6])
    Y.add_node(4,type='rr',c={6})
    Y.add_node(5,type='pr',c={1})
    Y.add_path([6,4,5,1])

    Y.add_node(7,type='c',start=False)
    Y.add_node(8,type='rf',c={6,7})
    Y.add_node(9,type='pf',c={12})
    Y.add_node(12,type='c',start=False)
    Y.add_path([7,8,9,12])
    Y.add_edge(6,8)
    Y.add_node(10,type='rr',c={12})
    Y.add_node(11,type='pr',c={6,7})
    Y.add_path([6,8,9,12])
    Y.add_edge(12,10)

    output_Y = DistanceToOrigin(Y, 4, -1)

    z = 0
    for node in Y.nodes():
        try:
            x = Y.node[node]['dist']
        except KeyError:
            z += 1

    assert z == 6
    assert [Y.node[n]['dist'] for n in range(1,7)] == [0,0,1,1,2,1]


def PruneNetwork(network):
    """Remove all nodes that are 'unreachable' defined as lacking a 'dist' data key."""
    for node in network.nodes():
        try:
            x = network.node[node]['dist']
        except KeyError:
            network.remove_node(node)

def test_PruneNetwork():
    # Nodes that were not reached in DistanceToOrigin are going to lack the 'dist' data key
    # Such nodes are expected to be removed
    G = nx.DiGraph()
    G.add_node(1,type='c',start=False)
    G.add_node(2,type='rf',c={1})
    G.add_node(3,type='pf',c={4})
    G.add_node(4,type='c',start=False)
    G.add_path([1,2,3,4])
    G.add_node(9,type='rr',c={4})
    G.add_node(10,type='pr',c={1})
    G.add_path([4,9,10,1])
    G.add_node(5,type='c',start=False)
    G.add_node(6,type='rf',c={1,5})
    G.add_node(7,type='pf',c={8})
    G.add_node(8,type='c',start=True) # Compound 8 is now the start
    G.add_path([1,6,7,8])
    G.add_edge(5,6)
    G.add_node(11,type='rr',c={8})
    G.add_node(12,type='pr',c={1,5})
    G.add_path([8,11,12,1])
    G.add_edge(12,5)

    H = G.copy()

    output = DistanceToOrigin(G, 4, 1)
    output = DistanceToOrigin(H, 4, 2)

    Y = G.subgraph([1,2,5,6,8,11,12])
    Z = H.subgraph([1,2,3,4,5,6,7,8,9,11,12])

    PruneNetwork(G)
    PruneNetwork(H)

    assert nx.is_isomorphic(G,Y)
    assert G.nodes(data=True) == Y.nodes(data=True)
    assert set(G.edges()) == set(Y.edges())

    assert nx.is_isomorphic(H,Z)
    assert H.nodes(data=True) == Z.nodes(data=True)
    assert set(H.edges()) == set(Z.edges())


def FindPaths(network, start_node, target_node, reaction_limit):
    """
    Find all simple paths from an origin reactant node to a target compound node,
    limiting the total number of reactions.
    """

    paths = []

    if network.node[start_node]['type'] not in {'rf','rr'} or network.node[target_node]['type'] != 'c':
        sys.stderr.write("Warning: Incorrect type detected in start and target nodes '%s' and '%s'.\n" % (str(start_node), str(target_node)))
        sys.stderr.flush()
        return paths

    if nx.has_path(network, start_node, target_node):
        for path in nx.shortest_simple_paths(network, start_node, target_node):
            if len(path)/3 > reaction_limit:
                # Dividing the path length by three gives the number of reaction steps
                break
            if nx.is_directed_acyclic_graph(network.subgraph(path)):
                # Only paths that represent an acyclic sub-network are allowed
                # Paths may not fold back onto themselves
                paths.append(path)

    return paths

def test_FindPaths(capsys):

    # Set up testing network
    G = nx.DiGraph()
    G.add_nodes_from(range(1,5), type='c')
    G.add_nodes_from(range(5,18,4), type='rf')
    G.add_nodes_from(range(6,19,4), type='pf')
    G.add_nodes_from(range(7,20,4), type='rr')
    G.add_nodes_from(range(8,21,4), type='pr')
    G.add_path([1,5,6,2,9,10,13,14,4,15,16,3,11,12,2,8,7,1])
    G.add_path([2,17,18,4,19,20,2])

    # Perform testing
    assert FindPaths(G, 5, 1, 5) == [] # Cyclicity
    assert FindPaths(G, 5, 4, 1) == [] # Path length limit 1
    assert FindPaths(G, 5, 4, 2) == [[5,6,2,17,18,4]] # Path length limit 2
    assert len(FindPaths(G, 5, 4, 3)) == 2 # Path length limit 3
    assert FindPaths(G, 13, 2, 2) == [[13,14,4,19,20,2]] # Different starting point + Cyclicity
    assert FindPaths(G, 1, 4, 5) == []
    out, err = capsys.readouterr()
    assert err == "Warning: Incorrect type detected in start and target nodes '1' and '4'.\n"


def SortPaths(paths):
    """Sort paths into bins based on their root - the product node at position -2."""
    path_dict = {}
    for path in paths:
        root_node = path[-2]
        try:
            path_dict[root_node].append(path)
        except KeyError:
            path_dict[root_node] = [path]
    return path_dict

def test_SortPaths():
    paths = [
    [1,2,3,4,5,6,7,8,9445],
    [10,20,30,40,50,9445],
    [60,70,80,7,8,9445],
    [125,92,119,3810,393,291,40,50,9445],
    [111,222,333,444,555,666,777,888,9445]
    ]

    assert [len(b) for b in [SortPaths(paths)[node] for node in [8,50,888]]] == [2,2,1]
    assert [paths[0],paths[2],paths[1],paths[3],paths[4]] == [x for y in [SortPaths(paths)[node] for node in [8,50,888]] for x in y]
    assert SortPaths(paths)[888][0] == paths[4]
    assert SortPaths(paths)[8][0] == paths[0]
    assert len(SortPaths(paths)) == 3

# Main code block
def main(infile_name, compound, reaction_limit, n_procs, simple, outfile_name):
    # Load and trim the network (data for every compound and reaction is lost)
    sys.stdout.write("\nLoading network pickle...")
    sys.stdout.flush()
    minet = pickle.load(open(infile_name, 'rb'))
    sys.stdout.write(" Done.\n")
    sys.stdout.flush()
    sys.stdout.write("Identifying starting compounds...")
    sys.stdout.flush()
    sys.stdout.write(" Done.\n")
    # Check for simple flag
    if simple:
        if n_procs != 1:
            sys.stdout.write("Performing network pruning using %s processes...\n" % n_procs)
        else:
            sys.stdout.write("Performing network pruning using one process...\n")
        dummy = DistanceToOrigin(minet, n_procs, reaction_limit)
        PruneNetwork(minet)
        sys.stdout.write("\nDone.\n")
        sys.stdout.flush()
    else:
        sys.exit('Not implemented. Try --simple.')
    sys.stdout.write("\nWriting results to pickle...")
    sys.stdout.flush()
    pickle.dump(minet, open(outfile_name, 'wb'))
    sys.stdout.write(" Done.\n")
    sys.stdout.flush()


if __name__ == "__main__":
    # Read arguments from the commandline
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='Read minet network pickle.')
    parser.add_argument('outfile', help='Write identified pathways as a list of graphs to output pickle.')
    parser.add_argument('-c', '--compound', type=str, help='Target compound.', required=True)
    parser.add_argument('-r', '--reactions', type=int, default=10, help='Maximum number of reactions.')
    parser.add_argument('-p', '--processes', type=int, default=1, help='Number of parallell processes to run.')
    parser.add_argument('--simple', action="store_true", help='Only calculate simple, direct paths (no branching).')
    args = parser.parse_args()
    main(args.infile, args.compound, args.reactions, args.processes, args.simple, args.outfile)
