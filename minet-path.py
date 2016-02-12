#!/usr/bin/env python3

# Import modules
import networkx as nx
import MineClient3 as mc
import sys
import argparse
import pickle
import multiprocessing as mp

# Define functions
def CountRxns(path):
    """Count the number of unique reactions in a path or list of nodes."""
    rxns = set()
    for node in path:
        if node[0] != 'c':
            rxns.add(node[1])
    return len(rxns)

def test_CountRxns():
    p1 = [('c','C1'), ('rf','R1'), ('pf','R1'), ('c','C2')]
    p2 = [('c','C3'), ('rr','R2'), ('pr','R2'), ('c','C4')]
    p3 = [('c','C5')]
    p4 = [('c','C6'), ('rf','R3'), ('pf','R3'), ('c','C7'), ('rr','R3'), ('rp','R3'), ('c','C8')]
    assert CountRxns(p1) == 1
    assert CountRxns(p2) == 1
    assert CountRxns(p3) == 0
    assert CountRxns(p4) == 1


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
                if not network.node[('c',c)]['start']:
                    dependencies.append(c)
    return dependencies


def test_CheckDependence():

    # Set up testing network
    p1 = [('c','C1'), ('rf','R1'), ('pf','R1'), ('c','C2')]
    p1b = [('c','C5'), ('rf','R1'), ('pf','R1'), ('c','C2')]
    p2 = [('c','C2'), ('rr','R1'), ('pr','R1'), ('c','C1')]
    p2b = [('c','C2'), ('rr','R1'), ('pr','R1'), ('c','C5')]
    p3 = [('c','C1'), ('rf','R3'), ('pf','R3'), ('c','C4')]
    p4 = [('c','C4'), ('rr','R3'), ('pr','R3'), ('c','C1')]
    G = nx.DiGraph()
    for p in [p1,p1b,p2,p2b,p3,p4]:
        G.add_path(p)
    G.node[('rf','R1')]['c'] = set(['C1','C5'])
    G.node[('pf','R1')]['c'] = set(['C2'])
    G.node[('rr','R1')]['c'] = set(['C2'])
    G.node[('pr','R1')]['c'] = set(['C1','C5'])
    G.node[('rf','R3')]['c'] = set(['C1'])
    G.node[('pf','R3')]['c'] = set(['C3'])
    G.node[('rr','R3')]['c'] = set(['C3'])
    G.node[('pr','R3')]['c'] = set(['C1'])

    for n in G.nodes():
        if n[0] == 'c':
            G.node[n]['start'] = False

    G.node[('c','C1')]['start'] = True

    assert set(CheckDependence(p1, G)) == set(['C5'])
    assert set(CheckDependence(p2, G)) == set(['C2'])
    assert set(CheckDependence(p3, G)) == set([])


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
        path_graph = 
        while


        if reaction_number >= reaction_limit:
            break




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


# Main code block
def main(infile_name, compound, reaction_limit, n_threads, outfile_name):
    # Code here
    return None


if __name__ == "__main__":
    # Read arguments from the commandline
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='Read minet network pickle.')
    parser.add_argument('outfile', help='')
    parser.add_argument('-c', type=str, help='')
    parser.add_argument('-r', type=int, default=10, help='')
    parser.add_argument('-N', type=int, default=1, help='')
    args = parser.parse_args()
    main()
