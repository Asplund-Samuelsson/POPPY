# Minet Origin node identification functions

# Import modules
import networkx as nx
import multiprocessing as mp

# Import scripts
from poppy_helpers import *

def find_start_comp_nodes(network):
    """Returns a list starting compound nodes in a MINE network."""
    start_comp_nodes = []
    for node in network.nodes():
        if network.node[node]['type'] == 'c' and network.node[node]['start']:
            start_comp_nodes.append(node)
    return set(start_comp_nodes)


def find_valid_reactant_nodes(network, proc_num=1, comp_node_set=set(),
    force_parallel=False):
    """Find valid reactant nodes given the supplied compound node set."""
    # If the set is empty, use starting compounds as the compound node set
    if comp_node_set == set():
        comp_node_set = find_start_comp_nodes(network)

    # Define the worker
    def worker(work):
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
        worker(network.nodes())
        valid_reactant_nodes = set(output)

    else:
        with mp.Manager() as manager:
            # Initialize output list in manager
            output = manager.list()

            # Initialize processes
            procs = []
            for work in chunks(network.nodes(), proc_num):
                p = mp.Process(target=worker, args=(work,))
                procs.append(p)
                p.start()

            # Stop workers
            for p in procs:
                p.join()

            # Get the results
            valid_reactant_nodes = set(output)

    return valid_reactant_nodes
