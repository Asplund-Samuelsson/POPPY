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


def TrimNetwork(network):
    """Trims the fat off a network in order to improve multiprocessing memory usage."""
    for node in network.nodes():
        del network.node[node]['data']
    return network


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
    assert set(CheckDependence(p2, G)) == set(['C2'])
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
    """Lists all paths that directly connects the starting compound to a target."""
    path_graphs = []
    for path in FindPaths(network, start_id, target_id, reaction_limit):
        if len(CheckDependence(path, network)) == 0:
            path_graphs.append(nx.subgraph(network, path))
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

    assert nx.is_isomorphic(path_C1_C3, GetDirectIndependentPaths(G, 'C1', 'C3', 3)[0])
    assert GetDirectIndependentPaths(G, 'C1', 'C4', 3) == []
    assert GetDirectIndependentPaths(G, 'C5', 'C6', 5) == [] # C5 is first in the path, but is non-starting and thus a dependency



# Main code block
def main(infile_name, compound, reaction_limit, n_procs, simple, outfile_name):
    # Load the network
    sys.stdout.write("\nLoading network pickle...")
    sys.stdout.flush()
    minetwork = pickle.load(open(infile_name, 'rb'))
    minetwork_lite = TrimNetwork(minetwork)
    sys.stdout.write(" Done.\n")
    sys.stdout.flush()
    sys.stdout.write("Identifying starting compounds...")
    sys.stdout.flush()
    start_comp_ids = GetStartCompIds(minetwork)
    sys.stdout.write(" Done.\n")
    # Check for simple flag
    if simple:
        sys.stdout.write("Performing pathfinding using %s processes..." % n_procs)
        pool = mp.Pool(processes=n_procs)
        M = len(start_comp_ids)
        arguments = zip(repeat(minetwork_lite, M), start_comp_ids, repeat(compound, M), repeat(reaction_limit, M))
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
