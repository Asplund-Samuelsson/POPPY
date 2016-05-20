#!/usr/bin/env python3

# Import modules
import networkx as nx
import sys
import argparse
import pickle
import multiprocessing as mp
import time
from datetime import timedelta as delta
import threading
import re
from itertools import product

# Import scripts
from progress import Progress
import mineclient3 as mc
from mepmap_origin_helpers import *
from mepmap_helpers import *

# Define functions
def count_reactions(network):
    """Count the number of unique reactions in a network."""
    rxns = set()
    for n in network.nodes():
        if network.node[n]['type'] != 'c':
            rxns.add(network.node[n]['mid'])
    return len(rxns)

def test_count_reactions():
    p1 = nx.DiGraph()
    p1.add_path(range(1,5))
    for node_data in enumerate([
            ('c','C1'), ('rf','R1'), ('pf','R1'), ('c','C2')
        ]):
        n = node_data[0] + 1
        p1.node[n]['type'] = node_data[1][0]
        p1.node[n]['mid'] = node_data[1][1]

    p2 = nx.DiGraph()
    p2.add_path(range(1,5))
    for node_data in enumerate([
            ('c','C3'), ('rr','R2'), ('pr','R2'), ('c','C4')
        ]):
        n = node_data[0] + 1
        p2.node[n]['type'] = node_data[1][0]
        p2.node[n]['mid'] = node_data[1][1]

    p3 = nx.DiGraph()
    p3.add_node(1, type='c', mid='C5')

    p4 = [
        ('c','C6'), ('rf','R3'), ('pf','R3'),
        ('c','C7'), ('rr','R3'), ('rp','R3'),
        ('c','C8')
    ]
    p4 = nx.DiGraph()
    p4.add_path(range(1,8))
    for node_data in enumerate([
            ('c','C6'), ('rf','R3'), ('pf','R3'),
            ('c','C7'), ('rr','R3'), ('rp','R3'),
            ('c','C8')
        ]):
        n = node_data[0] + 1
        p4.node[n]['type'] = node_data[1][0]
        p4.node[n]['mid'] = node_data[1][1]

    assert count_reactions(p1) == 1
    assert count_reactions(p2) == 1
    assert count_reactions(p3) == 0
    assert count_reactions(p4) == 1


def find_paths(network, reactant_node, compound_node, reaction_limit):
    """Find all simple paths from an origin reactant node to a target
    compound node, limiting the total number of reactions."""

    paths = []

    if nx.has_path(network, reactant_node, compound_node):
        for path in nx.shortest_simple_paths(network, reactant_node, compound_node):
            if len(path)/3 > reaction_limit:
                # Dividing path length by three yields reaction step number
                break
            if nx.is_directed_acyclic_graph(network.subgraph(path)):
                # Only paths that represent an acyclic sub-network are allowed
                # Paths may not fold back onto themselves
                paths.append(path)

    return paths

def test_find_paths():

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

    # Cyclicity
    assert find_paths(G, 5, 1, 5) == []
    # Path length limit 1
    assert find_paths(G, 5, 4, 1) == []
    # Path length limit 2
    assert find_paths(G, 5, 4, 2) == [[5,6,2,17,18,4]]
    # Path length limit 3
    assert len(find_paths(G, 5, 4, 3)) == 2
    # Different starting point + Cyclicity
    assert find_paths(G, 13, 2, 2) == [[13,14,4,19,20,2]]


def generate_paths(network, target_node, reaction_limit, n_procs=1, quiet=False):
    """Generate a list of paths to a target node from origin reactant nodes."""

    if not quiet:
        s_out("Generating paths...\n")

    # Define the worker
    def worker():
        while True:
            origin_node = Work.get()
            if origin_node is None:
                break
            paths = find_paths(network, origin_node, target_node, reaction_limit)
            output.extend(paths)
            with lock:
                n_work_done.value += 1

    # Set up reporter thread
    def reporter():
        t_start = time.time()
        n_left = Work.qsize()
        while n_left > 0:
            # Check progress
            n_made = len(output)
            n_left = Work.qsize()
            n_done = n_work - n_left
            # Calculate speed and time left
            speed = n_done / (time.time()-t_start)
            if n_done != 0:
                t_left = round(n_left / speed)
                t_left_str = str(delta(seconds=t_left))
            else:
                t_left_str = '.:..:..'
            status_format = "{0:<10} {1:<25} {2:<25}"
            progress = "%0.1f%%" % float(n_done / n_work * 100)
            found = str(n_made) + " paths found."
            t_out = "Time left: " + t_left_str
            status = status_format.format(progress, found, t_out)
            s_out('\r' + status)
            time.sleep(1)

    with mp.Manager() as manager:
        # Initialize Work queue in manager
        Work = manager.Queue()

        for origin_node in find_valid_reactant_nodes(network):
            Work.put(origin_node)

        # Place stop signals on queue
        for i in range(n_procs):
            Work.put(None)

        n_work = Work.qsize()

        # Initialize output list in manager
        output = manager.list()

        # Initialize number of tasks done counter
        n_work_done = mp.Value('i', 0)
        lock = mp.Lock()

        # Start reporter
        if not quiet:
            reporter = threading.Thread(target=reporter)
            reporter.start()

        # Start processes
        procs = []
        for i in range(n_procs):
            p = mp.Process(target=worker)
            procs.append(p)
            p.start()

        # Wait until all work is done
        while n_work_done.value != n_work - n_procs:
            time.sleep(1)

        # Terminate the processes
        for p in procs:
            if p.is_alive():
                p.terminate()

        # Stop reporter
        if not quiet:
            reporter.join()
            s_out("\nDone.\n")

        return list(output)


def test_generate_paths():
    G = nx.DiGraph()

    G.add_nodes_from([1,4,9,14,22,27,32,37,42,47,52], type='c', start=False)
    for node in [1,32,37,42,52]: G.node[node]['start'] = True

    rf = [2,7,20,12,101,25,30,35,40,45,55,50]
    rr = [5,10,23,18,15,28,33,38,43,48,57,53]
    pf = [x + 1 for x in rf]
    pr = [x + 1 for x in rr]

    assert set(rf).intersection(set(rr)) == set()
    assert len(set(rf)) == len(set(rr)) == 12
    assert set(G.nodes()).intersection(set(rf)) == set()
    assert set(G.nodes()).intersection(set(rr)) == set()
    assert set(G.nodes()).intersection(set(pf)) == set()
    assert set(G.nodes()).intersection(set(pr)) == set()

    G.add_nodes_from(rf, type='rf')
    G.add_nodes_from(pf, type='pf')
    G.add_nodes_from(rr, type='rr')
    G.add_nodes_from(pr, type='pr')

    G.add_path([1,2,3,4,7,8,9,20,21,22])
    G.add_path([22,23,24,9,10,11,4,5,6,1])
    G.add_path([4,12,13,14,101,102,9])
    G.add_path([9,18,19,14,15,16,4])
    G.add_path([32,33,34,27,28,29,22])
    G.add_path([22,25,26,27,30,31,32])
    G.add_path([22,35,36,37,40,41,42,45,46,47,55,56,22])
    G.add_path([22,57,58,47,48,49,42,43,44,37,38,39,22])
    G.add_path([52,53,54,47,50,51,52])

    for node in rf + rr: G.node[node]['c'] = set(G.predecessors(node))

    for node in pf + pr: G.node[node]['c'] = set(G.successors(node))

    only_one_reactant = True

    for node in G.nodes():
        if G.node[node]['type'] in {'rf','rr'}:
            if len(G.node[node]['c']) != 1:
                only_one_reactant = False

    assert only_one_reactant


    paths = generate_paths(G, 22, 5, n_procs=4)

    assert len(paths) == 9

    assert [2,3,4,7,8,9,20,21,22] in paths
    assert [2,3,4,12,13,14,101,102,9,20,21,22] in paths
    assert [33,34,27,28,29,22] in paths
    assert [38,39,22] in paths
    assert [43,44,37,38,39,22] in paths
    assert [53,54,47,48,49,42,43,44,37,38,39,22] in paths
    assert [53,54,47,55,56,22] in paths
    assert [45,46,47,55,56,22] in paths
    assert [40,41,42,45,46,47,55,56,22] in paths

    paths_2 = generate_paths(G, 9, 2, n_procs=4)

    assert len(paths_2) == 2

    assert [2,3,4,7,8,9] in paths_2
    assert [38,39,22,23,24,9] in paths_2

    assert generate_paths(G, 9, 1) == []


def nodes_being_produced(network):
    """Create set of compound nodes produced by the reactions in the network."""
    produced_nodes = set()
    for node in network.nodes():
        if network.node[node]['type'] in {'pf','pr'}:
            produced_nodes = produced_nodes.union(network.node[node]['c'])
    return produced_nodes

def test_nodes_being_produced():
    G = nx.DiGraph()
    G.add_nodes_from([1,2,3], type='c')
    G.add_nodes_from([10,20], type='pf')
    G.add_nodes_from([30,40], type='pr')
    G.add_node(50, type='rf', c=set([1]))
    G.node[10]['c'] = set([1,2])
    G.node[20]['c'] = set([2,3])
    G.node[30]['c'] = set([4])
    G.node[40]['c'] = set([5])

    assert nodes_being_produced(G) == set([1,2,3,4,5])


def nodes_being_consumed(network):
    """Create set of compound nodes consumed by the reactions in the network."""
    consumed_nodes = set()
    for node in network.nodes():
        if network.node[node]['type'] in {'rf','rr'}:
            consumed_nodes = consumed_nodes.union(network.node[node]['c'])
    return consumed_nodes

def test_nodes_being_consumed():
    G = nx.DiGraph()
    G.add_nodes_from([1,2,3], type='c')
    G.add_nodes_from([10,20], type='rf')
    G.add_nodes_from([30,40], type='rr')
    G.add_node(50, type='pf', c=set([1]))
    G.node[10]['c'] = set([1,2])
    G.node[20]['c'] = set([2,3])
    G.node[30]['c'] = set([4])
    G.node[40]['c'] = set([5])

    assert nodes_being_consumed(G) == set([1,2,3,4,5])


def remove_incomplete_reactions(network):
    """
    Removes reactions that require reactants not provided as start compounds or
    as products in other reactions of the network.

    Removes compound nodes if they are connected only to the affected reaction
    and are not a starting compound.

    The process is iterative, as removal of one reaction may make
    additional reactions 'incomplete' in terms of reactant presence.
    """

    start_comp_nodes = find_start_comp_nodes(network)

    while True:

        # Stop iteration if no nodes were removed in the previous round
        try:
            if prev_node_count == len(network.nodes()):
                break
        except NameError:
            # First round yields a NameError
            pass

        # Determine what compounds are available
        available = start_comp_nodes.union(nodes_being_produced(network))

        # Determine what reactions should be removed based on
        # unfulfilled reactant requirements
        nodes_to_remove = set()
        for node in network.nodes():
            if network.node[node]['type'] in {'rf','rr'}:
                if not network.node[node]['c'].issubset(available):
                    nodes_to_remove.add(node)
                    nodes_to_remove.add(network.successors(node)[0])
                    # Reactant nodes have one successor, i.e. a product node
                    # Go through and remove compounds directly downstream of the
                    # reaction, if they are connected only to this reaction
                    downstream = network.successors(network.successors(node)[0])
                    for comp_node in downstream:
                        if network.node[comp_node]['type'] != 'c':
                            s_err("Warning: '" + str(comp_node) + \
                            "' is not a compound node as expected (successor" +\
                            " of " + str(network.successors(node)[0]) + ")")
                            continue
                        else:
                            dependent = network.predecessors(comp_node) == \
                            network.successors(node)
                            start = network.node[comp_node]['start']
                            if dependent and not start:
                                nodes_to_remove.add(comp_node)

        # Count the number of nodes before purging
        prev_node_count = len(network.nodes())

        # Remove the nodes
        for node in nodes_to_remove:
            network.remove_node(node)

def test_remove_incomplete_reactions():
    G = nx.DiGraph()

    # Reactions tier 1
    G.add_node(1,type='rf',c=set([101,102])) # Origin reaction
    G.add_node(2,type='pf',c=set([103]))
    G.add_node(3,type='rr',c=set([104])) # Origin reaction
    G.add_node(4,type='pr',c=set([105]))

    # Reactions tier 2

    # Non-origin reaction that is complete:
    G.add_node(5,type='rf',c=set([103]))
    G.add_node(6,type='pf',c=set([106,107]))
    # Non-origin reaction that is incomplete:
    G.add_node(7,type='rf',c=set([105,108]))
    G.add_node(8,type='pf',c=set([109]))

    # Reactions tier 3
    G.add_node(9,type='rr',c=set([107,109]))
    G.add_node(10,type='pr',c=set([110]))

    # Compounds (only start compounds and those in
    # the direct path from origin to target)
    G.add_nodes_from([101,102,104], type='c', start=True)
    G.add_nodes_from([103,105,107,109,110], type='c', start=False)

    # Add edges
    G.add_path([101,1,2,103,5,6,107,9,10,110])
    G.add_path([104,3,4,105,7,8,109,9,10,110])
    G.add_edge(102,1)

    remove_incomplete_reactions(G)

    assert set(G.nodes()) == set([101,102,103,104,105,107,1,2,3,4,5,6])


    H = nx.DiGraph()

    H.add_node(101, type='rf', c=set([1]))
    H.add_node(102, type='pf', c=set([2]))
    H.add_node(201, type='rr', c=set([3]))
    H.add_node(202, type='pf', c=set([2,4]))

    H.add_node(1, type='c', start=True)
    H.add_nodes_from([2,3,4], type='c', start=False)

    H.add_path([1,101,102,2])
    H.add_path([3,201,202,2])
    H.add_edge(202,4)

    remove_incomplete_reactions(H)

    assert set(H.nodes()) == set([1,101,102,2,3])


def find_branch_nodes(network, severed=False):
    """Identifies reactant nodes that represent branch points in the pathway
    network.

    IMPORTANT: This function assumes that non-start predecessors leading to a
    severed branch are not present in the network, as would be the case for
    compounds not found in a direct path from origin to target."""

    branch_nodes = set()

    for node in network.nodes():

        # Branch nodes are reactant nodes
        if network.node[node]['type'] in {'rf','rr'}:
            S = 0

            # Iterate over predecessors, which are compound nodes
            for predecessor in network.predecessors(node):
                if network.node[predecessor]['start']:
                    S += 1

            # Start compound nodes should always be included in the network
            # Severed non-start compound nodes might not be included
            N = len(network.node[node]['c']) - S

            # Branch reactant nodes have more than one predecessor
            if len(network.node[node]['c']) > 1:
                if not severed:
                    branch_nodes.add(node)
                else:
                    # A severed branch node lacks 1 or more predecessors
                    # compared to what is listed in the 'c' set
                    pred_num = len(network.predecessors(node))
                    c_num = len(network.node[node]['c'])
                    if pred_num < c_num:
                        branch_nodes.add(node)

    return branch_nodes

def test_find_branch_nodes():
    # Set up testing network
    G = nx.DiGraph()

    # Add Compounds
    G.add_nodes_from([1,2,4,5,6,10], type='c', start=True)
    G.add_nodes_from([3,7,8,9,11,12,13], type='c', start=False)

    # Add reactions
    rf = [101,103,105,107,109,111,113,115]
    G.add_nodes_from(rf, type='rf')
    G.add_nodes_from([x+1 for x in rf], type='pf')
    G.node[107]['type'] = 'rr' # Don't forget reverse reactions
    G.node[108]['type'] = 'pr'

    # Add reactant sets
    G.node[101]['c'] = set([1])
    G.node[103]['c'] = set([2])
    G.node[105]['c'] = set([3])
    G.node[107]['c'] = set([4])
    G.node[109]['c'] = set([5])
    G.node[111]['c'] = set([6,7])
    G.node[113]['c'] = set([8,9])
    G.node[115]['c'] = set([10,11,12])

    # Add product sets
    G.node[102]['c'] = set([3])
    G.node[104]['c'] = set([3])
    G.node[106]['c'] = set([7])
    G.node[108]['c'] = set([8])
    G.node[110]['c'] = set([9])
    G.node[112]['c'] = set([11])
    G.node[114]['c'] = set([12])
    G.node[116]['c'] = set([13])

    # Add paths
    G.add_path([1,101,102,3,105,106,7,111,112,11,115,116,13])
    G.add_path([4,107,108,8,113,114,12,115])
    G.add_path([5,109,110,9,113])
    G.add_path([2,103,104,3])
    G.add_edge(6,111)
    G.add_edge(10,115)

    # Ensure branch node identification is working
    assert find_branch_nodes(G) == set([111,113,115])

    # Sever branches and then identify them
    G.remove_nodes_from([107,108,8])
    assert find_branch_nodes(G, severed=True) == set([113])
    G.remove_nodes_from([109,110,9])
    assert find_branch_nodes(G, severed=True) == set([113])
    G.remove_nodes_from([105,106,7])
    assert find_branch_nodes(G, severed=True) == set([111,113])

    H = G.copy()
    G.remove_nodes_from([111,112,11])
    H.remove_nodes_from([113,114,12])
    assert find_branch_nodes(G, severed=True) == set([113,115])
    assert find_branch_nodes(H, severed=True) == set([111,115])

    G.remove_nodes_from([113,114,12])
    assert find_branch_nodes(G, severed=True) == set([115])


def find_switch_nodes(network):
    """
    Identifies and returns a list of 'switch nodes' in the network.

    Switch nodes are compound nodes with multiple predecessors and thereby
    multiple paths of synthesis.
    """
    switch_nodes = set()
    for node in network.nodes():
        if network.node[node]['type'] == 'c':
            if len(network.predecessors(node)) > 1:
                switch_nodes.add(node)
    return switch_nodes

def test_find_switch_nodes():
    # Switch nodes are compound nodes that have multiple predecessors
    G = nx.DiGraph()
    G.add_nodes_from([1,2], type='c')
    G.add_nodes_from([101,201,301], type='pf')
    G.add_edges_from([(101,1),(201,2),(301,2)])

    assert find_switch_nodes(G) == set([2])


def generate_termini(network):
    """
    Generates termini and eliminates start-compound induced branching by
    cutting all start compound to reactant node edges.
    """
    for node in network.nodes():
        if network.node[node]['type'] == 'c':
            if network.node[node]['start']:
                network.remove_edges_from(network.out_edges(node))

def test_generate_termini():
    G = nx.DiGraph()
    G.add_nodes_from([1,3], type='c', start=True)
    G.add_nodes_from([2,4], type='c', start=False)
    G.add_node(101, type='rf')
    G.add_node(102, type='pf')
    G.add_node(201, type='rr')
    G.add_node(202, type='pr')
    G.add_path([1,101,102,2])
    G.add_path([2,201,202,4])
    G.add_edge(3,201)

    H = G.copy()

    generate_termini(G)

    assert len(G.nodes()) == len(H.nodes())
    assert set(G.edges()) == {(101,102),(102,2),(2,201),(201,202),(202,4)}
    assert nx.has_path(G, 101, 4)


def digraph_connected_component(network, target_node):
    """
    Performs a reverse depth-first search to find all nodes that have a path
    leading to the target node. Then returns that component of the network.
    """
    return set(nx.dfs_preorder_nodes(network.reverse(), target_node))

def test_digraph_connected_component():
    G = nx.DiGraph()
    G.add_path([1,2,3,4])
    G.add_path([1,5,6,7])
    G.add_path([8,9,10,4])
    G.add_path([11,12,13,14])

    H = nx.DiGraph()
    H.add_path([1,2,3,4])
    H.add_path([8,9,10,4])

    Y = G.subgraph(digraph_connected_component(G,4))

    assert nx.is_isomorphic(H, Y)
    assert set(H.nodes()) == set(Y.nodes())


def subnetwork_from_paths(network, paths, target_node):

    s_out("\nConstructing sub-network from identified paths...")

    # Acquire basic data
    path_nodes = set([n for p in paths for n in p])
    start_comp_nodes = find_start_comp_nodes(network)

    # Construct a sub-network of all paths,
    # start compounds and produced compounds
    subnet = network.subgraph(path_nodes.union(start_comp_nodes))
    subnet = network.subgraph(
        set(subnet.nodes()).union(nodes_being_produced(subnet))
    )

    n_removed = 1
    while n_removed:
        n_before = len(subnet)
        # Identify the incomplete reactions and remove them
        remove_incomplete_reactions(subnet)
        # Reduce to the connected component
        subnet = subnet.subgraph(digraph_connected_component(subnet, target_node))
        n_removed = n_before - len(subnet)
        # Exit with an error message if the target node was removed
        if target_node not in subnet.nodes():
            sys.exit("\nError: Target node cannot be produced by the sub-network.\n")

    s_out(" Done.\n")

    return subnet

def test_subnetwork_from_paths():
    # Test case has the following in common with test_generate_paths:
    G = nx.DiGraph()

    G.add_nodes_from([1,4,9,14,22,27,32,37,42,47,52], type='c', start=False)
    for node in [1,32,37,42,52]: G.node[node]['start'] = True

    rf = [2,7,20,12,101,25,30,35,40,45,55,50]
    rr = [5,10,23,18,15,28,33,38,43,48,57,53]
    pf = [x + 1 for x in rf]
    pr = [x + 1 for x in rr]

    G.add_nodes_from(rf, type='rf')
    G.add_nodes_from(pf, type='pf')
    G.add_nodes_from(rr, type='rr')
    G.add_nodes_from(pr, type='pr')

    G.add_path([1,2,3,4,7,8,9,20,21,22])
    G.add_path([22,23,24,9,10,11,4,5,6,1])
    G.add_path([4,12,13,14,101,102,9])
    G.add_path([9,18,19,14,15,16,4])
    G.add_path([32,33,34,27,28,29,22])
    G.add_path([22,25,26,27,30,31,32])
    G.add_path([22,35,36,37,40,41,42,45,46,47,55,56,22])
    G.add_path([22,57,58,47,48,49,42,43,44,37,38,39,22])
    G.add_path([52,53,54,47,50,51,52])

    # Now spice it with a compound that cannot be produced by the network
    G.add_node(100, type='c', start=False)
    G.add_edge(100, 28)
    G.add_edge(26, 100)

    # Finally, list the reactant and product nodes
    for node in rf + rr: G.node[node]['c'] = set(G.predecessors(node))
    for node in pf + pr: G.node[node]['c'] = set(G.successors(node))

    paths = generate_paths(G, 22, 4, n_procs=4)
    subnet = subnetwork_from_paths(G, paths, 22)

    assert len(subnet.nodes()) == 33
    assert len(subnet.edges()) == 36
    assert set(subnet.nodes()) == set([
        2,3,4,7,8,12,13,14,101,
        102,9,20,21,22,56,55,
        47,54,53,46,45,39,38,
        1,37,42,52,
        44,43,48,49,40,41
    ])


def format_graphml(network, subnet):
    # Need to label the origin reactant nodes
    origins = find_valid_reactant_nodes(network)

    for node in subnet.nodes():
        if subnet.node[node]['type'] in {'rf','rr'}:
            if node in origins:
                subnet.node[node]['origin'] = True
            else:
                subnet.node[node]['origin'] = False

    # Also need to create compound name labels
    for node in subnet.nodes():
        if subnet.node[node]['type'] == 'c':
            mid = network.node[node]['mid']
            try:
                common_name = network.graph['mine_data'][mid]['Names'][0]
            except KeyError:
                common_name = network.graph['mine_data'][mid]['Formula']
            subnet.node[node]['common_name'] = common_name


    # Furthermore, label reaction nodes with the reactants or products
    for node in subnet.nodes():
        if subnet.node[node]['type'] != 'c':
            c_names = []
            for c_node in subnet.node[node]['c']:
                mid = network.node[c_node]['mid']
                try:
                    common_name = network.graph['mine_data'][mid]['Names'][0]
                except KeyError:
                    common_name = network.graph['mine_data'][mid]['Formula']
                c_names.append(common_name)
            subnet.node[node]['common_name'] = " + ".join(c_names)

    # Copy and remove the incompatible stuff from nodes and graph
    outnet = subnet.copy()

    for node in outnet.nodes():
        if outnet.node[node]['type'] != 'c':
            del outnet.node[node]['c']

    data_keys = list(outnet.graph.keys())
    for key in data_keys:
        del outnet.graph[key]

    return outnet


def paths_to_pathways(network, paths, target_node):
    """Enumerate complete branched pathways capable of producing the target"""

    # Function for determining whether a path is part of a network
    def path_in_network(path, network):
        for i in range(1, len(path)):
            n = path[i]
            if n in network.nodes():
                predecessors = network.predecessors(n)
            else:
                predecessors = []
            if path[i-1] not in predecessors:
                return False
        return True

    # Construct a subnetwork
    subnet = subnetwork_from_paths(network, paths, target_node)

    # Determine compounds nodes available in the subnetwork
    start_comp_nodes = find_start_comp_nodes(network) # Start compounds
    prod_comp_nodes = nodes_being_produced(subnet)
    tot_comp = start_comp_nodes.union(prod_comp_nodes)

    # Filter the supplied paths to those that are present in the subnetwork
    paths_filtered = []
    p = Progress(design = 'pt', max_val = len(paths))
    n = 0
    for path in paths:
        n += 1
        s_out("\rFiltering paths... %s" % p.to_string(n))
        if path_in_network(path, subnet):
            paths_filtered.append(path)

    print("")

    # Generate a dictionary with path segments producing the key node
    segments = {}
    p = Progress(max_val = len(paths_filtered))
    m = 0

    # Go through all filtered paths
    for path in paths_filtered:

        # Report progress
        m += 1
        s_out("\rGenerating path segments... %s" % p.to_string(m))

        # Iterate over the elements of the path
        for element in enumerate(path):
            i = element[0]
            n = element[1]

            # Check if product node
            if subnet.node[n]['type'] in {'pf','pr'}:

                # Construct segments
                path_segs = []
                offset = 0
                while True:
                    path_segs.append(path[i - (1 + offset):i + 1])
                    try:
                        # If the upstream compound node is not start
                        if not subnet.node[path[i - (2 + offset)]]['start']:
                            # Stop construction
                            break
                    # If there are no more nodes upstream
                    except IndexError:
                        # Stop construction
                        break
                    offset += 3

                # Iterate over products
                for segment in path_segs:
                    for c_node in subnet.node[n]['c']:

                        # Save the segment(s) under every produced node
                        try:
                            segments[c_node].add(tuple(segment + [c_node]))
                        except KeyError:
                            segments[c_node] = set([tuple(segment + [c_node])])

    print("")

    # Pathway enumeration
    print("\nEnumerating pathways...")

    # Storage container for finished pathways
    finished_pathways = set()
    detected_cycles = set()
    unfinished_pathways = set([frozenset(p) for p in segments[target_node]])

    # Progress setup
    max_length = 0
    min_length = count_reactions(subnet)
    p = Progress(design='s')
    p_form = '{0} Finished: {1:<12} Unfinished: {2:<10} Reactions (min/max): {3:^3}/{4:^3}'

    # Iterate through unfinished pathways
    while unfinished_pathways:

        # Pop a path off the set of unfinished pathways
        path = set(unfinished_pathways.pop())

        # Check if the pathway has any known cycles
        if not sum([c.issubset(pathnet.nodes()) for c in detected_cycles]):

            # Extract the paths sub-network
            pathnet = subnet.subgraph(path)

            # Determine what compounds are produced by the path
            cpd_prod = nodes_being_produced(pathnet)

            # Determine what compounds are consumed by the path
            cpd_cons = nodes_being_consumed(pathnet)

            # Check if the path is complete on its own
            if cpd_cons.issubset(cpd_prod.union(start_comp_nodes)):

                # Make a copy for cycle identification
                path_c = pathnet.copy()

                # Cut off start compounds to allow benign cycles
                generate_termini(path_c)

                # Identify cycles
                cycles = set([frozenset(c) for c in nx.simple_cycles(path_c)])
                # Cycles are indicative of compounds that cannot be produced
                # de novo by the pathway, but enter due to a bootstrap effect
                if cycles:
                    detected_cycles = detected_cycles.union(cycles)
                else:
                    max_length = max(max_length, count_reactions(pathnet))
                    min_length = min(min_length, count_reactions(pathnet))
                    finished_pathways.add(frozenset(path))

            # If not, add all combinations of segments that might complement it
            else:

                # Identify the missing compound nodes
                missing = cpd_cons - cpd_prod - start_comp_nodes

                # Construct sets of complementary segments that will satisfy the
                # current missing compound nodes
                try:
                    for complement in product(*[segments[i] for i in missing]):
                        unfinished_pathways.add(
                            frozenset(set(path).union(*complement))
                        )
                except KeyError:
                    # If a complement cannot be found, discard the pathway
                    pass

        # Report progress
        n_done = len(finished_pathways)
        n_left = len(unfinished_pathways)
        p_msg = "\r" + p_form.format(p.to_string(), n_done, n_left, min_length, max_length)
        s_out(p_msg)

    print("")

    return finished_pathways

def test_paths_to_pathways():
    # Set up testing network - Same as for Identifyfind_branch_nodes
    G = nx.DiGraph()

    # Add Compounds
    G.add_nodes_from([1,2,4,5,6,10], type='c', start=True)
    G.add_nodes_from([3,7,8,9,11,12,13], type='c', start=False)

    # Add reactions
    rf = [101,103,105,107,109,111,113,115]
    G.add_nodes_from(rf, type='rf')
    G.add_nodes_from([x+1 for x in rf], type='pf')
    G.node[107]['type'] = 'rr' # Don't forget reverse reactions
    G.node[108]['type'] = 'pr'

    # Add reactant sets
    G.node[101]['c'] = set([1])
    G.node[103]['c'] = set([2])
    G.node[105]['c'] = set([3])
    G.node[107]['c'] = set([4])
    G.node[109]['c'] = set([5])
    G.node[111]['c'] = set([6,7])
    G.node[113]['c'] = set([8,9])
    G.node[115]['c'] = set([10,11,12])

    # Add product sets
    G.node[102]['c'] = set([3])
    G.node[104]['c'] = set([3])
    G.node[106]['c'] = set([7])
    G.node[108]['c'] = set([8])
    G.node[110]['c'] = set([9])
    G.node[112]['c'] = set([11])
    G.node[114]['c'] = set([12])
    G.node[116]['c'] = set([13])

    # Add paths
    G.add_path([1,101,102,3,105,106,7,111,112,11,115,116,13])
    G.add_path([4,107,108,8,113,114,12,115])
    G.add_path([5,109,110,9,113])
    G.add_path([2,103,104,3])
    G.add_edge(6,111)
    G.add_edge(10,115)

    # Add three extra paths to increase complexity

    # One leading through 6 (a starting compound)
    G.add_path([20,201,202,6])
    G.node[20]['type'] = 'c'
    G.node[20]['start'] = True
    G.node[201]['type'] = 'rr'
    G.node[201]['c'] = set([20])
    G.node[202]['type'] = 'pr'
    G.node[202]['c'] = set([6])

    # One adding an additional option to produce 7
    G.add_path([30,301,302,7])
    G.node[30]['type'] = 'c'
    G.node[30]['start'] = True
    G.node[301]['type'] = 'rf'
    G.node[301]['c'] = set([20])
    G.node[302]['type'] = 'pf'
    G.node[302]['c'] = set([7])

    # One adding an additional root node and a route beginning in 5
    G.add_path([5,401,402,40,403,404,13])
    G.node[401]['type'] = 'rr'
    G.node[401]['c'] = set([5])
    G.node[402]['type'] = 'pr'
    G.node[402]['c'] = set([40])
    G.node[40]['type'] = 'c'
    G.node[40]['start'] = False
    G.node[403]['type'] = 'rr'
    G.node[403]['c'] = set([40])
    G.node[404]['type'] = 'pr'
    G.node[404]['c'] = set([13])

    # Also adding a shortcut from 40 to 12 that will result in an invalid rxn
    G.add_path([40,501,502,12])
    G.add_edge(50,501)
    G.node[50]['type'] = 'c'
    G.node[50]['start'] = False # This compound should yield an incomplete rxn
    G.node[501]['type'] = 'rr'
    G.node[501]['c'] = set([40,50])
    G.node[502]['type'] = 'pr'
    G.node[502]['c'] = set([12])

    # Add 'mids'
    mids = [
        (107,'R1'), (108,'R1'), (109,'R2'), (110,'R2'), (401,'R3'), (402,'R3'),
        (201,'R4'), (202,'R4'), (301,'R5'), (302,'R5'), (101,'R6'), (102,'R6'),
        (103,'R7'), (104,'R7'), (113,'R8'), (114,'R8'), (403,'R9'), (404,'R9'),
        (501,'Ra'), (502,'Ra'), (105,'Rb'), (106,'Rb'), (111,'Rc'), (112,'Rc'),
        (115,'Rd'), (116,'Rd')
    ]
    for m in mids:
        G.node[m[0]]['mid'] = m[1]

    # These are the complete sets of nodes representing branched pathways
    req_1 = set([7,111,112,11,115,116,13,12,114,113,9,8,109,110,107,108])
    expected_branched_paths = {
    frozenset(set([101,102,3,105,106]).union(req_1)),
    frozenset(set([103,104,3,105,106]).union(req_1)),
    frozenset(set([301,302]).union(req_1)),
    frozenset(set([201,202,6,101,102,3,105,106]).union(req_1)),
    frozenset(set([201,202,6,103,104,3,105,106]).union(req_1)),
    frozenset(set([201,202,6,301,302]).union(req_1)),
    frozenset([401,402,40,403,404,13])
    }

    # Letting the automated functions produce a result
    paths = generate_paths(G, 13, 5, quiet=True)
    output_branched_paths = paths_to_pathways(G, paths, 13)

    paths_equal = True
    missing = []
    unexpected = []
    for path in expected_branched_paths:
        if path not in output_branched_paths:
            paths_equal = False
            missing.append(sorted(list(path)))
    for path in output_branched_paths:
        if path not in expected_branched_paths:
            paths_equal = False
            unexpected.append(sorted(list(path)))

    if not paths_equal:
        s_err("Missing: " + str(missing) + "\n")
        s_err("Unexpected: " + str(unexpected) + "\n")

    assert paths_equal
    assert len(expected_branched_paths) == len(output_branched_paths)

    # Test for parallel paths detection and acceptance
    N = nx.DiGraph()
    N.add_nodes_from([1,7], type='c', start=True)
    N.add_nodes_from([2,3,4,5,6], type='c', start=False)
    N.add_nodes_from([101,301,501], type='rf')
    N.add_nodes_from([102,302,502], type='pf')
    N.add_nodes_from([103,303,503], type='rr')
    N.add_nodes_from([104,304,504], type='pr')
    N.add_nodes_from([201,401], type='rr')
    N.add_nodes_from([202,402], type='pr')
    N.add_nodes_from([203,403], type='rf')
    N.add_nodes_from([204,404], type='pf')
    N.add_path([1,101,102,2,201,202,5,401,402,6])
    N.add_path([6,403,404,5,203,204,2,103,104,1])
    N.add_path([1,101,102,3,301,302,4,401,402,6])
    N.add_path([6,403,404,4,303,304,3,103,104,1])
    N.add_path([7,501,502,4,503,504,7])
    for node in N.nodes():
        if N.node[node]['type'] in {'rf','rr'}:
            N.node[node]['c'] = set(N.predecessors(node))
        if N.node[node]['type'] in {'pf','pr'}:
            N.node[node]['c'] = set(N.successors(node))

    # Add 'mids'
    mids = [
        (101,'R1'), (102,'R1'), (103,'R2'), (104,'R2'), (201,'R3'), (202,'R3'),
        (203,'R4'), (204,'R4'), (301,'R5'), (302,'R5'), (303,'R6'), (304,'R6'),
        (401,'R7'), (402,'R7'), (403,'R8'), (404,'R8'), (501,'R9'), (502,'R9'),
        (503,'Ra'), (504,'Ra')
    ]
    for m in mids:
        N.node[m[0]]['mid'] = m[1]

    N_expected_paths = {
    frozenset([101,102,2,201,202,5,401,402,6,3,301,302,4]),
    frozenset([101,102,2,201,202,5,401,402,6,501,502,4])
    }

    paths = generate_paths(N, 6, 5, quiet=True)
    output_branched_paths = paths_to_pathways(N, paths, 6)

    paths_equal = True
    missing = []
    unexpected = []
    for path in N_expected_paths:
        if path not in output_branched_paths:
            paths_equal = False
            missing.append(sorted(list(path)))
    for path in output_branched_paths:
        if path not in N_expected_paths:
            paths_equal = False
            unexpected.append(sorted(list(path)))

    if not paths_equal:
        s_err("Missing: " + str(missing) + "\n")
        s_err("Unexpected: " + str(unexpected) + "\n")

    assert paths_equal
    assert len(N_expected_paths) == len(output_branched_paths)

    # Test for pathways going through starting compounds
    H = nx.DiGraph()
    rf_nodes = [101,201,301,111]
    pf_nodes = [102,202,302,112]
    H.add_nodes_from([1,2,3,4,5], type='c', start=True) # All start compounds
    H.add_nodes_from(rf_nodes, type='rf')
    H.add_nodes_from(pf_nodes, type='pf')
    H.add_path([1,101,102,2,201,202,3,301,302,4])
    H.add_path([5,111,112,2])
    for n in rf_nodes:
        H.node[n]['c'] = set(H.predecessors(n))
    for n in pf_nodes:
        H.node[n]['c'] = set(H.successors(n))
    for i in range(len(rf_nodes)):
        mid = 'R' + str(i)
        H.node[rf_nodes[i]]['mid'] = mid
        H.node[pf_nodes[i]]['mid'] = mid

    H_expected_pathways = {
    frozenset([301,302,4]),
    frozenset([201,202,3,301,302,4]),
    frozenset([101,102,2,201,202,3,301,302,4]),
    frozenset([111,112,2,201,202,3,301,302,4])
    }

    paths = generate_paths(H, 4, 5, quiet=True)
    output_branched_paths = paths_to_pathways(H, paths, 4)

    paths_equal = True
    missing = []
    unexpected = []
    for path in H_expected_pathways:
        if path not in output_branched_paths:
            paths_equal = False
            missing.append(sorted(list(path)))
    for path in output_branched_paths:
        if path not in H_expected_pathways:
            paths_equal = False
            unexpected.append(sorted(list(path)))

    if not paths_equal:
        s_err("Missing: " + str(missing) + "\n")
        s_err("Unexpected: " + str(unexpected) + "\n")

    assert paths_equal
    assert len(H_expected_pathways) == len(output_branched_paths)


def parse_compound(compound, network):
    """Determines the type of compound identifier and returns the node."""

    node = None

    mid_match = re.match('^[CX]{1}[0-9,a-f]{40}$', compound)
    kegg_match = re.match('^C{1}[0-9]{5}$', compound)

    if mid_match:
        try:
            node = network.graph['cmid2node'][compound]
            if not node in network.nodes():
                s_err("Error: MINE ID '" + compound + \
                "' appears to not be available in the network.\n")
                node = None
        except KeyError:
            s_err("Error: MINE ID '" + compound + \
            "' appears to not be available in the network.\n")
    elif kegg_match:
        try:
            nodes = network.graph['kegg2nodes'][compound]
            if len(nodes) > 1:
                nodes = "'\n'".join(
                    sorted([network.node[n]['mid'] for n in nodes])
                )
                s_err("Error: '" + compound + "' refers to multiple nodes." + \
                " Use --exact_comp_id and --compound followed by one of:\n'" + \
                nodes + "'\n")
            else:
                node = list(nodes)[0]
                if not node in network.nodes():
                    s_err("Error: KEGG ID '" + compound + \
                    "' appears to not be available in the network.\n")
                    node = None
        except KeyError:
            s_err("Error: KEGG ID '" + compound + \
            "' appears to not be available in the network.\n")
    else:
        try:
            nodes = network.graph['name2nodes'][compound]
            if len(nodes) > 1:
                nodes = "'\n'".join(
                    sorted([network.node[n]['mid'] for n in nodes])
                )
                s_err("Error: '" + compound + "' refers to multiple nodes." + \
                " Use --exact_comp_id and --compound followed by one of:\n'" + \
                nodes + "'\n")
            else:
                node = list(nodes)[0]
                if not node in network.nodes():
                    s_err("Error: Name '" + compound + \
                    "' appears to not be available in the network.\n")
                    node = None
        except KeyError:
            s_err("Error: Name '" + compound + \
            "' appears to not be available in the network.\n")

    return node

def test_parse_compound(capsys):
    G = nx.DiGraph()
    G.graph['cmid2node'] = {}
    G.graph['kegg2nodes'] = {}
    G.graph['name2nodes'] = {}

    G.graph['cmid2node'] = {
    'Cf647c96ae2e66c3b6ab160faa1d8498be5112fe4':1,
    'C6a6f4d5234ea2b14b42c391eb760d6311afa8388':2,
    'C31b986bd97ef9152d7534436e847379077d4553c':3,
    'Cbef44c34f09ee9b80ad4b9daaa32afe088b41507':4,
    'C5a2e2841cff1008380531689c8c45b6dbecd04b6':5,
    'C0736b9ad03e466caa7698fbd3fccf06f6654fb53':6,
    'C0d3cb8256055b59b846cd5930d69c266bb13deb1':7,
    'C12c16f3e8910911f982fe6fcd541c35bca59119e':8,
    'Ce4a58113b67f1e7edb22e28123f300f36b763903':9,
    'C31890':10,
    'C00291':11,
    'C67559':12,
    'C67560':13,
    'C99999':14
    }

    G.graph['kegg2nodes'] = {
    'C31890':set([1]), 'C00291':set([2]), 'C67559':set([3]),
    'C67560':set([3]), 'C99999':set([14])
    }

    G.graph['name2nodes'] = {
    'Alpha':set([7]),
    'Beta':set([8]),
    'Gamma':set([9]),
    'n-Alpha':set([7]),
    'n-Beta':set([8]),
    'n-Gamma':set([9]),
    'Twin':set([5,8])
    }

    G.add_nodes_from(G.graph['cmid2node'].values())
    node2cmid = {v: k for k, v in G.graph['cmid2node'].items()}
    for n in G.nodes():
        G.node[n]['mid'] = node2cmid[n]

    # Test KEGG IDs
    assert parse_compound('C31890',G) == 1
    assert parse_compound('C00291',G) == 2
    assert parse_compound('C67559',G) == 3
    assert parse_compound('C67560',G) == 3

    # Test MINE IDs
    res = [parse_compound(c,G) for c in [node2cmid[n] for n in range(1,10)]]
    assert res == list(range(1,10))

    # Test Names
    assert parse_compound('Alpha',G) == parse_compound('n-Alpha',G) == 7
    assert parse_compound('Beta',G) == parse_compound('n-Beta',G) == 8
    assert parse_compound('Gamma',G) == parse_compound('n-Gamma',G) == 9
    assert parse_compound('Twin',G) == None

    # Test error output
    assert parse_compound('Cffffffffffffffffffffffffffffffffffffffff',G) == None
    assert parse_compound('C00001', G) == None
    assert parse_compound('Delta', G) == None

    exp_err = "Error: 'Twin' refers to multiple nodes. Use --exact_comp_id" + \
    " and --compound followed by one of:\n" + \
    "'C12c16f3e8910911f982fe6fcd541c35bca59119e'\n" + \
    "'C5a2e2841cff1008380531689c8c45b6dbecd04b6'\n"
    exp_err = exp_err + "Error: MINE ID " + \
    "'Cffffffffffffffffffffffffffffffffffffffff' appears to not be " + \
    "available in the network.\n"
    exp_err = exp_err + "Error: KEGG ID 'C00001' appears to not be " + \
    "available in the network.\n"
    exp_err = exp_err + "Error: Name 'Delta' appears to not be " + \
    "available in the network.\n"

    out, err = capsys.readouterr()

    assert err == exp_err


def update_start_compounds(network, start_comp_ids):
    """Update the start compound status for each compound node."""
    S = set(start_comp_ids)
    for n in network.nodes():
        if network.node[n]['type'] == 'c':
            if network.node[n]['mid'] in S:
                start = True
            else:
                start = False
            network.node[n]['start'] = start

def test_update_start_compounds():
    G = nx.DiGraph()

    G.add_node(1,type='c',mid='C1',start=True) # Will be changed
    G.add_node(2,type='c',mid='C2',start=False) # Will be changed
    G.add_node(3,type='c',mid='C3',start=True) # Will not be changed
    G.add_node(4,type='c',mid='C4',start=False) # Will not be changed
    G.add_node(5,type='rf',mid='R1')
    G.add_node(6,type='pf',mid='R1')
    G.add_node(7,type='rr',mid='R1')
    G.add_node(8,type='pr',mid='R1')

    G.add_path([1,5,6,2])

    H = G.copy()
    H.node[1]['start'] = False
    H.node[2]['start'] = True

    start_comp_ids = ['C2','C3']

    update_start_compounds(G, start_comp_ids)

    assert G.nodes(data=True) == H.nodes(data=True)
    assert nx.is_isomorphic(G, H)


def format_pathway_text(network, pathways, target_node):
    """Create pathways record in text format"""
    pathway_lines = []

    # Go through pathways in length order
    for pathway in sorted(list(pathways), key=len):

        # Create a subnetwork
        subnet = network.subgraph(pathway)

        # Identify the reactant nodes in the pathway
        reactant_nodes = [
            n[0] for n in subnet.nodes(data=True) if n[1]['type'] in {'rf','rr'}
        ]

        # Find the order of the reactions
        rxn_nodes = sorted(
            reactant_nodes,
            key = lambda rn : len(nx.shortest_path(subnet, rn, target_node)),
            reverse = True)

        # Add reactions in the detected order
        for n in rxn_nodes:
            rxn_type = subnet.node[n]['type']
            rxn_id = subnet.node[n]['mid']
            rxn = network.graph['mine_data'][rxn_id]
            if rxn_type == 'rf':
                r = rxn['Reactants']
                p = rxn['Products']
            else:
                r = rxn['Products']
                p = rxn['Reactants']

            # Flatten lists and add plus characters
            r = [a for b in list(joinit(r, ['+'])) for a in b]
            p = [a for b in list(joinit(p, ['+'])) for a in b]

            # Combine and add reaction arrow
            rxn_elements = r + ['<=>'] + p

            # Remove "1" coefficients
            rxn_elements = list(filter(lambda x : x is not 1, rxn_elements))

            # Create text
            pathway_lines.append(
                rxn_id + "\t" + " ".join([str(x) for x in rxn_elements])
            )

        # Add pathway divider
        pathway_lines.append("//")

    # Return formatted pathways text
    return "\n".join(pathway_lines) + "\n"

def test_format_pathway_text():

    # Set up the testing network
    N = nx.DiGraph()
    N.graph['mine_data'] = {
        'R1':{'Reactants':[[1,'C4']],
              'Products':[[2,'C1']]},
        'R2':{'Reactants':[[1,'C2'],[1,'C3']],
              'Products':[[1,'C5']]},
        'R3':{'Reactants':[[1,'C4'],[1,'C5']],
              'Products':[[1,'C6'],[1,'X0']]},
        'R4':{'Reactants':[[1,'X9'],[2,'C7']],
              'Products':[[1,'C6'],[1,'X8']]}
    }
    N.add_nodes_from([0,1,2,3,4,5,6,7,8,9], type='c')
    N.add_nodes_from([11,21,31,41], type='rf')
    N.add_nodes_from([12,22,32,42], type='pf')
    N.add_nodes_from([13,23,33,43], type='rr')
    N.add_nodes_from([14,24,34,44], type='pr')
    for n in [11,12,13,14,21,22,23,24,31,32,33,34,41,42,43,44]:
        N.node[n]['mid'] = 'R' + str(n)[0]
    N.add_path([1,13,14,4,31,32,6,43,44,7])
    N.add_path([7,41,42,6,33,34,4,11,12,1])
    N.add_edges_from([
        (42,8),(8,43),(44,9),(9,41),
        (3,21),(24,3),(34,5),(5,31),
        (32,0),(0,33)
    ])
    N.add_path([2,21,22,5])
    N.add_path([5,23,24,2])

    # Define the pathways
    pathways = {
    frozenset([13,14,4,31,32,6,43,44,7,21,22,5]),
    frozenset([31,32,6,43,44,7,21,22,5])
    }

    # Describe the expected text output
    exp_pathway_text = "\n".join([
    "R2\tC2 + C3 <=> C5",
    "R3\tC4 + C5 <=> C6 + X0",
    "R4\tC6 + X8 <=> X9 + 2 C7",
    "//",
    "R1\t2 C1 <=> C4",
    "R2\tC2 + C3 <=> C5",
    "R3\tC4 + C5 <=> C6 + X0",
    "R4\tC6 + X8 <=> X9 + 2 C7",
    "//",
    ]) + "\n"

    assert format_pathway_text(N, pathways, 7) == exp_pathway_text


# Main code block
def main(infile_name, compound, start_comp_id_file, exact_comp_id,
    reaction_limit, n_procs, sub_network_out, outfile_name, pathway_text):

    # Default results are empty
    results = {}

    # Load the network
    s_out("\nLoading network pickle...")
    network = pickle.load(open(infile_name, 'rb'))
    s_out(" Done.\n")

    # Update starting compounds
    if start_comp_id_file:
        S = [L.rstrip() for L in open(start_comp_id_file, 'r').readlines()]
        update_start_compounds(network, S)

    # Pathway enumeration
    if not exact_comp_id:
        target_node = parse_compound(compound, network)
    else:
        try:
            target_node = network.graph['cmid2node'][compound]
        except KeyError:
            target_node = None

    if target_node == None:
        sys.exit("Error: Target node was not found. Check compound '" + \
        compound + "'.\n")

    paths = generate_paths(network, target_node, reaction_limit, n_procs)

    # Save sub-network graphml
    if sub_network_out:
        subnet = subnetwork_from_paths(network, paths, target_node)
        s_out("\nWriting sub-network to graphml...")
        subnet = format_graphml(network, subnet)
        nx.write_graphml(subnet, sub_network_out)
        s_out(" Done.\n")

    # Enumerate pathways and save results
    if outfile_name or pathway_text:
        pathways = paths_to_pathways(network, paths, target_node)
        if outfile_name:
            s_out("\nWriting pathways to pickle...")
            pickle.dump(pathways, open(outfile_name, 'wb'))
            s_out(" Done.\n")
        if pathway_text:
            s_out("\nWriting pathways to text...")
            with open(pathway_text, 'w') as txt:
                txt.write(format_pathway_text(network, pathways, target_node))
            s_out(" Done.\n")


if __name__ == "__main__":
    # Read arguments from the commandline
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'infile',
        help='Read mepmap network pickle.'
    )
    parser.add_argument(
        'compound', type=str, default=False,
        help='Target compound.'
    )
    parser.add_argument(
        '-S', '--start_comp_ids', type=str,
        help='Provide new starting compounds in a file.'
    )
    parser.add_argument(
        '-o', '--outfile', type=str, default=False,
        help='Save identified pathways in pickle.'
    )
    parser.add_argument(
        '-t', '--pathway_text', type=str, default=False,
        help='Save identified pathways as a text file.'
    )
    parser.add_argument(
        '-s', '--sub_network', type=str, default=False,
        help='Save sub-network as graphml (requires -c).'
    )
    parser.add_argument(
        '-e', '--exact_comp_id', action='store_true',
        help='Look for exact compound ID.'
    )
    parser.add_argument(
        '-r', '--reactions', type=int, default=5,
        help='Maximum number of reactions.'
    )
    parser.add_argument(
        '-p', '--processes', type=int, default=1,
        help='Number of parallel processes to run.'
    )
    args = parser.parse_args()
    main(args.infile, args.compound, args.start_comp_ids, args.exact_comp_id,
    args.reactions, args.processes, args.sub_network, args.outfile,
    args.pathway_text)
