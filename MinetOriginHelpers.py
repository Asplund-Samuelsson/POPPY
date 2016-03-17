# Minet Origin node identification functions

# Import modules
import networkx as nx
import multiprocessing as mp

# Import scripts
from MinetHelpers import *

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
