# Minet Origin node identification functions

# Add repository root to the path
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

# Import the script to be tested
from poppy_origin_helpers import *

# Define tests
def test_find_start_comp_nodes():
    G = nx.DiGraph()
    G.add_node(1,type='c',start=True)
    G.add_node(2,type='rf')
    G.add_node(3,type='pf')
    G.add_node(4,type='c',start=False)
    G.add_path([1,2,3,4])
    assert find_start_comp_nodes(G) == set([1])


def test_find_valid_reactant_nodes():
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
    start_comps = find_start_comp_nodes(G)

    assert find_valid_reactant_nodes(G, 4, start_comps) == set([2])
    assert find_valid_reactant_nodes(G, 4, start_comps, \
    force_parallel=True) == set([2])

    assert find_valid_reactant_nodes(G, 4) == set([2])
    assert find_valid_reactant_nodes(G, 4, force_parallel=True) == set([2])

    assert find_valid_reactant_nodes(G, 4, set([1,5])) == set([2,6])
    assert find_valid_reactant_nodes(G, 4, set([1,5]), \
    force_parallel=True) == set([2,6])

    assert find_valid_reactant_nodes(G, 4, set([8])) == set([11])
    assert find_valid_reactant_nodes(G, 4, set([8]), \
    force_parallel=True) == set([11])
