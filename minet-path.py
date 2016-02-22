#!/usr/bin/env python3

# Import modules
import networkx as nx
import MineClient3 as mc
import sys
import argparse
import pickle
import multiprocessing as mp
from collections import defaultdict as dd
from itertools import repeat

# Define functions
def CountRxns(path):
    """Count the number of unique reactions in a path graph."""
    rxns = set()
    for node in path.nodes(data=True):
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


def TrimNetwork(network):
    """Trims the fat off a network in order to improve multiprocessing memory usage."""
    del network.graph['mine_data']

def test_TrimNetwork():
    # Set up testing networks
    G = nx.DiGraph()
    Y = nx.DiGraph()
    p1 = [1,2,3]
    p2 = [1,5,6]
    p3 = [5,7,4,2]
    for p in [p1,p2,p3]:
        G.add_path(p)
        Y.add_path(p)
    for n in G.nodes():
        G.node[n]['mid'] = 'dummy'
        Y.node[n]['mid'] = 'dummy'

    G.graph['mine_data'] = {'R1':{'_id':'R1'}}
    Y.graph['mine_data'] = {'R1':{'_id':'R1'}}

    # Trim G but not Y
    TrimNetwork(G)

    # Topography should be the same
    assert nx.is_isomorphic(G, Y)

    # Data should be gone from G
    try:
        x = G.graph['mine_data']
        d = True
    except KeyError:
        d = False

    assert d == False


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


def FindValidReactantNodes(network, comp_node_set=set([])):
    """Find and return all reactant nodes that are valid given the supplied compound node set."""
    # If the set is empty, use starting compounds as the compound node set
    if comp_node_set == set([]):
        comp_node_set = FindStartCompNodes(network)

    valid_reactant_nodes = set([])

    for node in network.nodes():
        if network.node[node]['type'] in {'rf','rr'}:
            c = network.node[node]['c']
            if c.issubset(comp_node_set):
                valid_reactant_nodes.add(node)

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
    assert FindValidReactantNodes(G, start_comps) == set([2])
    assert FindValidReactantNodes(G) == set([2])
    assert FindValidReactantNodes(G, set([1,5])) == set([2,6])
    assert FindValidReactantNodes(G, set([8])) == set([11])


def ExpandValidCompoundSet(network, valid_reactant_nodes=set([]), comp_node_set=set([])):
    """Expands the valid compound set with products of valid reactant nodes."""
    if comp_node_set == set([]):
        comp_node_set = FindStartCompNodes(network)
    else:
        for r_node in valid_reactant_nodes:
            p_node = network.successors(r_node)[0]
            comp_node_set = comp_node_set.union(network.node[p_node]['c'])
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
    assert ExpandValidCompoundSet(G, set([2,9]), set([1,4])) == set([1,4])
    assert ExpandValidCompoundSet(G, FindValidReactantNodes(G, set([8])), set([8])) == set([1,5,8])
    assert ExpandValidCompoundSet(G, set([2,6,11]), ExpandValidCompoundSet(G, set([11]), set([8]))) == set([1,4,5,8])


def DistanceToOrigin(network, N=-1):
    """
    Calculates the shortest distance (number of reactions) from the starting compounds (origin) to
    every node up to distance N. Set N to -1 to exhaustively calculate the minimum distance to
    every node that is reachable.

    Returns two sets in a tuple: Valid compound nodes and valid reactant nodes.
    """

    sys.stdout.write("\nCalculating minimum distance of nodes to origin...\n\n")
    sys.stdout.flush()

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

    while True:

        # Valid product nodes will be connected at the start of an expansion cycle
        # They are however in principle identified in the previous cycle via valid reactant nodes
        try:
            for r_node in new_vrn:
                p_node = network.successors(r_node)[0]
                network.node[p_node]['dist'] = n
                node_type = network.node[p_node]['type']
                if node_type == 'pf':
                    pf += 1
                if node_type == 'pr':
                    pr += 1
        except NameError:
            # new_vrn is not present in the first round
            pass

        # Expand the valid compound set
        # When n = 0, this means the starting compounds
        # When n > 0, the valid compound set will be expanded
        new_vcn = ExpandValidCompoundSet(network, valid_reactant_nodes, valid_compound_nodes) - valid_compound_nodes
        valid_compound_nodes = new_vcn.union(valid_compound_nodes)

        new_vrn = FindValidReactantNodes(network, valid_compound_nodes) - valid_reactant_nodes
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

        if N == -1:
            if set(prev_vrn) == valid_reactant_nodes and set(prev_vcn) == valid_compound_nodes:
                # When no new valid compound or reactant nodes have been identified, it is time to stop
                break
        else:
            if n > N:
                # n starts at 0 and increments by one before each round dealing
                # with that particular step n - stop when n exceeds the limit
                break
        prev_vrn = list(valid_reactant_nodes)
        prev_vcn = list(valid_compound_nodes)

    sys.stdout.write("\nDone.\n")
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

    output_0 = DistanceToOrigin(G.copy(), 0) # Creating a copy, since the function modifies the network it is given
    output_1 = DistanceToOrigin(G.copy(), 1)
    output_2 = DistanceToOrigin(G, 2) # Not creating a copy in order to check modification capabilities

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

    output_Y = DistanceToOrigin(Y, -1)

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

    output = DistanceToOrigin(G, 1)
    output = DistanceToOrigin(H, 2)

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


def DetectLoops(path):
    """Detect reaction loops, which induce 'bootstrapping', in a path."""
    rxns = dd(int)
    for node in path:
        if node[0] != 'c':
            rxns[node[1]] +=1
    for n in rxns.values():
        if n > 2:
            return True
    return False

def test_DetectLoops():
    p1 = [('c','C1'), ('rf','R1'), ('pf','R1'), ('c','C2')]
    p2 = [('c','C1'), ('rf','R1'), ('pf','R1'), ('c','C2'), ('rr','R1'), ('pr','R1'), ('c','C3')]
    assert DetectLoops(p2)
    assert not DetectLoops(p1)


def FindPaths(network, start_id, target_id, reaction_limit):
    """Find all simple paths from a starting compound to a target compound,
    limiting the total number of reactions."""
    start_node = ('c', start_id)
    target_node = ('c', target_id)
    paths = []
    if nx.has_path(network, start_node, target_node):
        for path in nx.shortest_simple_paths(network, start_node, target_node, weight='weight'):
            if CountRxns(path) > reaction_limit:
                break
            if not DetectLoops(path):
                # Paths that loop through reactions are not allowed
                paths.append(path)
    return paths

def test_FindPaths():

    # Set up testing network
    p1 = [('c','C1'), ('rf','R1'), ('pf','R1'), ('c','C2')]
    p2 = [('c','C2'), ('rf','R2'), ('pf','R2'), ('c','C3')]
    p3 = [('c','C3'), ('rr','R2'), ('pr','R2'), ('c','C2')]
    p4 = [('c','C2'), ('rr','R1'), ('pr','R1'), ('c','C1')]
    p5 = [('c','C1'), ('rf','R3'), ('pf','R3'), ('c','C4')]
    p6 = [('c','C4'), ('rf','R4'), ('pf','R4'), ('c','C2')]
    p7 = [('c','C4'), ('rr','R3'), ('pr','R3'), ('c','C1')]
    p8 = [('c','C2'), ('rr','R4'), ('pr','R4'), ('c','C4')]
    G = nx.DiGraph()
    for p in [p1,p2,p3,p4,p5,p6,p7,p8]:
        G.add_path(p)
    for e in G.edges():
        if e[0][0] == 'c' or e[1][0] == 'c':
            G.edge[e[0]][e[1]]['weight'] = 0
        else:
            G.edge[e[0]][e[1]]['weight'] = 1

    # Perform testing
    assert FindPaths(G, 'C1', 'C3', 3)[1][3] == ('c','C4')
    assert len(FindPaths(G, 'C1', 'C3', 10)) == 2
    assert len(FindPaths(G, 'C1', 'C3', 2)) == 1


def GetStartCompIds(network):
    """Lists the IDs of all starting compounds in the network."""
    start_comp_ids = []
    for node in network.nodes(data=True):
        if node[0][0] == 'c':
            if node[1]['start']:
                start_comp_ids.append(node[0][1])
    return start_comp_ids

def test_GetStartCompIds():
    G = nx.DiGraph()
    G.add_nodes_from([('c','C1'),('c','C2'),('c','C3')], start=False)
    G.add_nodes_from([('rf','R1'),('pf','R1')])
    G.node[('c','C2')]['start'] = True
    assert len(GetStartCompIds(G)) == 1
    assert GetStartCompIds(G)[0] == 'C2'


def CheckDependence(path, network):
    """Returns IDs of non-starting compounds among any of the reactants along the path."""
    dependencies = []
    for n in path:
        if n[0] in {'rf','rr'}:
            for c in network.node[n]['c']:
                try:
                    # The compound is a dependency if it is not in the direct path and not a starting compound
                    if not network.node[('c',c)]['start'] and ('c',c) not in path:
                        dependencies.append(c)
                except KeyError:
                    continue
    return dependencies


def test_CheckDependence():

    # Set up testing network
    p1 = [('c','C1'), ('rf','R1'), ('pf','R1'), ('c','C2')]
    p1b = [('c','C5'), ('rf','R1'), ('pf','R1'), ('c','C2')]
    p2 = [('c','C2'), ('rr','R1'), ('pr','R1'), ('c','C1')]
    p2b = [('c','C2'), ('rr','R1'), ('pr','R1'), ('c','C5')]
    p3 = [('c','C1'), ('rf','R3'), ('pf','R3'), ('c','C4')]
    p4 = [('c','C4'), ('rr','R3'), ('pr','R3'), ('c','C1')]
    p5 = [('c','C4'), ('rf','R4'), ('pf','R4'), ('c','C6')]
    p5b = [('c','C6'), ('rr','R4'), ('pr','R4'), ('c','C4')]
    G = nx.DiGraph()
    for p in [p1,p1b,p2,p2b,p3,p4,p5,p5b]:
        G.add_path(p)
    G.node[('rf','R1')]['c'] = set(['C1','C5'])
    G.node[('pf','R1')]['c'] = set(['C2'])
    G.node[('rr','R1')]['c'] = set(['C2'])
    G.node[('pr','R1')]['c'] = set(['C1','C5'])
    G.node[('rf','R3')]['c'] = set(['C1'])
    G.node[('pf','R3')]['c'] = set(['C3'])
    G.node[('rr','R3')]['c'] = set(['C3'])
    G.node[('pr','R3')]['c'] = set(['C1'])
    G.node[('rf','R4')]['c'] = set(['C4'])
    G.node[('pf','R4')]['c'] = set(['C6'])
    G.node[('rr','R4')]['c'] = set(['C6'])
    G.node[('pr','R4')]['c'] = set(['C4'])

    for n in G.nodes():
        if n[0] == 'c':
            G.node[n]['start'] = False

    G.node[('c','C1')]['start'] = True

    assert set(CheckDependence(p1, G)) == set(['C5'])
    assert set(CheckDependence(p2, G)) == set([])
    assert set(CheckDependence(p3, G)) == set([])
    assert set(CheckDependence(p3[0:3] + p5, G)) == set([])


def EnumeratePaths(network, start_id, target_id, reaction_limit):
    """Enumerate branched pathways from the starting compound to the target compound,
    limiting the total number of reactions. Returns a graph describing the pathway."""
    start_comp_ids = GetStartCompIds(network)
    # Recursively? find paths
    path_graphs = []
    for path in FindPaths(network, start_id, target_id, reaction_limit):
        reaction_number = 0
        path_len = CountRxns(path)
        reaction_number += path_len

    return None

def test_EnumeratePaths():

    # Construct testing network
    G = nx.DiGraph()

    p1 = [('c','C1'),('rf','R1'),('pf','R1'),('c','C3'),('rf','R3'),('pf','R3'),('c','C5')]
    p2 = [('c','C5'),('rr','R3'),('pr','R3'),('c','C3'),('rr','R1'),('pr','R1'),('c','C1')]
    p3 = [('c','C2'),('rf','R2'),('pf','R2'),('c','C4'),('rf','R3'),('pf','R3'),('c','C5')]
    p4 = [('c','C5'),('rr','R3'),('pr','R3'),('c','C4'),('rr','R2'),('pr','R2'),('c','C2')]
    p5 = [('c','C6'),('rf','R4'),('pf','R4'),('c','C5')]
    p6 = [('c','C5'),('rr','R4'),('pr','R4'),('c','C6')]
    for p in [p1,p2,p3,p4,p5,p6]:
        G.add_path(p)

    G.node[('rf','R1')]['c'] = set(['C1'])
    G.node[('pf','R1')]['c'] = set(['C3'])
    G.node[('rr','R1')]['c'] = set(['C3'])
    G.node[('pr','R1')]['c'] = set(['C1'])

    G.node[('rf','R2')]['c'] = set(['C2'])
    G.node[('pf','R2')]['c'] = set(['C4'])
    G.node[('rr','R2')]['c'] = set(['C4'])
    G.node[('pr','R2')]['c'] = set(['C2'])

    G.node[('rf','R3')]['c'] = set(['C3','C4'])
    G.node[('pf','R3')]['c'] = set(['C5'])
    G.node[('rr','R3')]['c'] = set(['C5'])
    G.node[('pr','R3')]['c'] = set(['C3','C4'])

    G.node[('rf','R4')]['c'] = set(['C6'])
    G.node[('pf','R4')]['c'] = set(['C5'])
    G.node[('rr','R4')]['c'] = set(['C5'])
    G.node[('pr','R4')]['c'] = set(['C6'])

    G.node[('c','C1')]['start'] = True
    G.node[('c','C2')]['start'] = True
    G.node[('c','C3')]['start'] = False
    G.node[('c','C4')]['start'] = False
    G.node[('c','C5')]['start'] = False
    G.node[('c','C6')]['start'] = False

    for e in G.edges():
        if e[0][0] == 'c' or e[1][0] == 'c':
            G.edge[e[0]][e[1]]['weight'] = 0
        else:
            G.edge[e[0]][e[1]]['weight'] = 1

    for n in G.nodes():
        G.node[n]['data'] = {'_id':n[1]}

    path_graph_C1_C5 = nx.compose(nx.subgraph(G,p1), nx.subgraph(G,p3))

    path_graph_C2_C5 = path_graph_C1_C5 # Both starting compounds should yield the same graph

    path_graph_C1_C6 = nx.compose(path_graph_C1_C5, nx.subgraph(G,p6))

    assert nx.is_isomorphic(EnumeratePaths(G, 'C1', 'C5', 5)[0], path_graph_C1_C5)
    assert nx.is_isomorphic(EnumeratePaths(G, 'C2', 'C5', 5)[0], path_graph_C2_C5)
    assert nx.is_isomorphic(EnumeratePaths(G, 'C1', 'C6', 5)[0], path_graph_C1_C6)
    assert len(EnumeratePaths(G, 'C1', 'C6', 4)) == 1
    assert len(EnumeratePaths(G, 'C1', 'C6', 3)) == 0


def GetDirectIndependentPaths(network, start_id, target_id, reaction_limit):
    """Lists all paths that directly connect the starting compound to a target."""
    path_graphs = []
    rejected = '!'
    for path in FindPaths(network, start_id, target_id, reaction_limit):
        if len(CheckDependence(path, network)) == 0:
            path_graphs.append(nx.subgraph(network, path))
            n = CountRxns(path)
            if n < 10:
                num = str(n)
            else:
                num = '*'
            sys.stdout.write(num)
            sys.stdout.flush()
        else:
            sys.stdout.write(rejected)
            sys.stdout.flush()
        rejected = '.'
    return path_graphs

def test_GetDirectIndependentPaths():
    # Construct testing network
    G = nx.DiGraph()

    p1 = [('c','C1'),('rf','R1'),('pf','R1'),('c','C3'),('rf','R3'),('pf','R3'),('c','C5')]
    p2 = [('c','C5'),('rr','R3'),('pr','R3'),('c','C3'),('rr','R1'),('pr','R1'),('c','C1')]
    p3 = [('c','C2'),('rf','R2'),('pf','R2'),('c','C4'),('rf','R3'),('pf','R3'),('c','C5')]
    p4 = [('c','C5'),('rr','R3'),('pr','R3'),('c','C4'),('rr','R2'),('pr','R2'),('c','C2')]
    p5 = [('c','C6'),('rf','R4'),('pf','R4'),('c','C5')]
    p6 = [('c','C5'),('rr','R4'),('pr','R4'),('c','C6')]
    for p in [p1,p2,p3,p4,p5,p6]:
        G.add_path(p)

    G.node[('rf','R1')]['c'] = set(['C1'])
    G.node[('pf','R1')]['c'] = set(['C3'])
    G.node[('rr','R1')]['c'] = set(['C3'])
    G.node[('pr','R1')]['c'] = set(['C1'])

    G.node[('rf','R2')]['c'] = set(['C2'])
    G.node[('pf','R2')]['c'] = set(['C4'])
    G.node[('rr','R2')]['c'] = set(['C4'])
    G.node[('pr','R2')]['c'] = set(['C2'])

    G.node[('rf','R3')]['c'] = set(['C3','C4'])
    G.node[('pf','R3')]['c'] = set(['C5'])
    G.node[('rr','R3')]['c'] = set(['C5'])
    G.node[('pr','R3')]['c'] = set(['C3','C4'])

    G.node[('rf','R4')]['c'] = set(['C6'])
    G.node[('pf','R4')]['c'] = set(['C5'])
    G.node[('rr','R4')]['c'] = set(['C5'])
    G.node[('pr','R4')]['c'] = set(['C6'])

    G.node[('c','C1')]['start'] = True
    G.node[('c','C2')]['start'] = True
    G.node[('c','C3')]['start'] = False
    G.node[('c','C4')]['start'] = False
    G.node[('c','C5')]['start'] = False
    G.node[('c','C6')]['start'] = False

    for e in G.edges():
        if e[0][0] == 'c' or e[1][0] == 'c':
            G.edge[e[0]][e[1]]['weight'] = 0
        else:
            G.edge[e[0]][e[1]]['weight'] = 1

    for n in G.nodes():
        G.node[n]['data'] = {'_id':n[1]}

    path_C1_C3 = nx.subgraph(G, p1[0:4])
    path_C1_C4 = nx.subgraph(G, p1 + p4[1:4]) # "Bootstrapped" path; C5 cannot exist without C4
    path_C5_C6 = nx.subgraph(G, p6) # C5 is non-starting, but included in path

    assert nx.is_isomorphic(path_C1_C3, GetDirectIndependentPaths(G, 'C1', 'C3', 3)[0])
    assert GetDirectIndependentPaths(G, 'C1', 'C4', 3) == []
    assert nx.is_isomorphic(path_C5_C6, GetDirectIndependentPaths(G, 'C5', 'C6', 5)[0])
    assert set(path_C5_C6.nodes()) == set(GetDirectIndependentPaths(G, 'C5', 'C6', 5)[0].nodes())



# Main code block
def main(infile_name, compound, reaction_limit, n_procs, simple, outfile_name):
    # Load and trim the network (data for every compound and reaction is lost)
    sys.stdout.write("\nLoading network pickle...")
    sys.stdout.flush()
    minetwork = pickle.load(open(infile_name, 'rb'))
    TrimNetwork(minetwork)
    sys.stdout.write(" Done.\n")
    sys.stdout.flush()
    sys.stdout.write("Identifying starting compounds...")
    sys.stdout.flush()
    start_comp_ids = GetStartCompIds(minetwork)
    sys.stdout.write(" Done.\n")
    # Check for simple flag
    if simple:
        sys.stdout.write("Performing pathfinding using %s processes...\n" % n_procs)
        pool = mp.Pool(processes=n_procs)
        M = len(start_comp_ids)
        arguments = zip(repeat(minetwork, M), start_comp_ids, repeat(compound, M), repeat(reaction_limit, M))
        results_0 = pool.starmap_async(GetDirectIndependentPaths, arguments)
        results = results_0.get()
        #n = 0
        #sys.stdout.write(" Done.\n\nRetrieving results...\n")
        #sys.stdout.flush()
        #for result in results_0:
        #    result = p.get()
        #    results.append(result)
        #    n += 1
        #    sys.stdout.write("%0.2f%\r" % n / M * 100)
        #    sys.stdout.flush()
        sys.stdout.write("\nDone.\n")
        sys.stdout.flush()
    else:
        sys.exit('Not implemented. Try --simple.')
    sys.stdout.write("\nWriting results to pickle...")
    sys.stdout.flush()
    pickle.dump(results, open(outfile_name, 'wb'))
    sys.stdout.write(" Done.\n")
    sys.stdout.flush()
    return None


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
