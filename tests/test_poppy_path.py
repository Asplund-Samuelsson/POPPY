#!/usr/bin/env python3

# Add repository root to the path
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

# Import the script to be tested
from poppy_path import *

# Define tests
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


def test_find_switch_nodes():
    # Switch nodes are compound nodes that have multiple predecessors
    G = nx.DiGraph()
    G.add_nodes_from([1,2], type='c')
    G.add_nodes_from([101,201,301], type='pf')
    G.add_edges_from([(101,1),(201,2),(301,2)])

    assert find_switch_nodes(G) == set([2])


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


def test_has_cycles():
    # Set up testing network
    Z = nx.DiGraph()

    # Add Compounds
    Z.add_nodes_from([1,2,7], type='c', start=True)
    Z.add_nodes_from([3,4,5,6], type='c', start=False)

    # Add reactions
    rf = [101,201,301,401]
    Z.add_nodes_from(rf, type='rf')
    Z.add_nodes_from([x+1 for x in rf], type='pf')
    Z.node[201]['type'] = 'rr' # Don't forget reverse reactions
    Z.node[202]['type'] = 'pr'

    # Add reactant sets
    Z.node[101]['c'] = set([1])
    Z.node[201]['c'] = set([2,6])
    Z.node[301]['c'] = set([3,4])
    Z.node[401]['c'] = set([7])

    # Add product sets
    Z.node[102]['c'] = set([3])
    Z.node[202]['c'] = set([4])
    Z.node[302]['c'] = set([5,6])
    Z.node[402]['c'] = set([2])

    # Add paths
    Z.add_path([1,101,102,3,301,302,1])
    Z.add_path([7,401,402,2,201,202,4,301])
    Z.add_edge(302,5)
    Z.add_path([302,6,201])

    # Current cycles are...
    # 1>101>102>3>301>302>1 (okay, 1 is start)
    # 6>201>202>4>301>302>6 (not okay, 6 is a bootstrap compound)

    Z_path_1 = Z.subgraph({1,2,101,102,201,202,3,4,301,302,5,6})
    assert has_cycles(Z_path_1, Z)

    Z_path_2 = Z.subgraph({1,2,101,102,201,202,3,4,301,302,5,6})
    Z_path_2.remove_edges_from([(302,6),(6,201)])
    assert has_cycles(Z_path_2, Z)

    Z_path_3 = Z.subgraph({7,401,402,2,201,202,4,301,302,5})
    assert has_cycles(Z_path_3, Z)

    Z_path_4 = Z.subgraph({1,101,102,3,301,302,5})
    assert not has_cycles(Z_path_4, Z)


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
    expected_branched_paths_shallow = {
    frozenset(set([101,102,3,105,106]).union(req_1)),
    frozenset(set([103,104,3,105,106]).union(req_1)),
    frozenset(set([301,302]).union(req_1)),
    frozenset([401,402,40,403,404,13])
    }

    # Letting the automated functions produce a result
    paths = generate_paths(G, 13, 5, quiet=True)

    # Standard procedure
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

    # Shallow enumeration
    shallow_pathways = paths_to_pathways(G, paths, 13, shallow=True)

    paths_equal = True
    missing = []
    unexpected = []
    for path in expected_branched_paths_shallow:
        if path not in shallow_pathways:
            paths_equal = False
            missing.append(sorted(list(path)))
    for path in shallow_pathways:
        if path not in expected_branched_paths_shallow:
            paths_equal = False
            unexpected.append(sorted(list(path)))

    if not paths_equal:
        s_err("Missing: " + str(missing) + "\n")
        s_err("Unexpected: " + str(unexpected) + "\n")

    assert paths_equal
    assert len(expected_branched_paths_shallow) == len(shallow_pathways)

    # Check the reaction limit
    exp_limited_paths = {
        frozenset(set([301,302]).union(req_1)),
        frozenset([401,402,40,403,404,13])
    }
    assert paths_to_pathways(G, paths, 13, 6) == exp_limited_paths

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

    exp_err = "Error: 'Twin' refers to multiple IDs. Use --exact_comp_id" + \
    " and one of:\n" + \
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

    # Test the return of sets
    assert parse_compound('Twin', G, return_set=True) == {5,8}
    assert parse_compound('C99999', G, return_set=True) == {14}


def test_update_start_compounds():
    G = nx.DiGraph()

    G.graph['cmid2node'] = {}
    G.graph['kegg2nodes'] = {}
    G.graph['name2nodes'] = {}

    G.graph['cmid2node'] = {
    'C1ae5f6786e5a8f65a865e865f68a5e68f5a86e8a':1,
    'C2ae5f6786e5a8f65a865e865f68a5e68f5a86e8a':2,
    'C3ae5f6786e5a8f65a865e865f68a5e68f5a86e8a':3,
    'C4ae5f6786e5a8f65a865e865f68a5e68f5a86e8a':4,
    'C00004':11,
    'Cf5dc8599a48d0111a3a5f618296752e1b53c8d30':12,
    'Xf5dc8599a48d0111a3a5f618296752e1b53c8d30':13
    }

    G.graph['kegg2nodes'] = {
    'C31890':set([1]), 'C00291':set([2]), 'C67559':set([3]),'C67560':set([3]),
    'C00004':set([11]), 'C30184':set([12])
    }

    G.graph['name2nodes'] = {
    'Alpha':set([1]),
    'Beta':set([2]),
    'Gamma':set([3]),
    'n-Alpha':set([1]),
    'n-Beta':set([2]),
    'n-Gamma':set([3]),
    'Twin':set([1,4]),
    'NAD':set([11]),
    'NADX':set([12])
    }

    G.graph['mine_data'] = {
    'C1ae5f6786e5a8f65a865e865f68a5e68f5a86e8a':{
    'Names':['Alpha','n-Alpha','Twin'],
    'DB_links':{'KEGG':['C31890']},
    },
    'C2ae5f6786e5a8f65a865e865f68a5e68f5a86e8a':{
    'Names':['Beta','n-Beta'],
    'DB_links':{'KEGG':['C00291']},
    },
    'C3ae5f6786e5a8f65a865e865f68a5e68f5a86e8a':{
    'Names':['Gamma','n-Gamma'],
    'DB_links':{'KEGG':['C67559','C67560']},
    },
    'C4ae5f6786e5a8f65a865e865f68a5e68f5a86e8a':{
    'Names':['Twin'],
    },
    'C00004':{
    'Names':['NAD','NADX','DPN'],
    'DB_links':{'KEGG':['C00004','C30184']},
    },
    'Cf5dc8599a48d0111a3a5f618296752e1b53c8d30':{
    'Names':['NAD','NADX','DPN'],
    'DB_links':{'KEGG':['C00004','C30184']},
    },
    'Xf5dc8599a48d0111a3a5f618296752e1b53c8d30':{
    'Names':['NAD','NADX','DPN'],
    'DB_links':{'KEGG':['C00004','C30184']},
    }
    }

    G.add_node(1,type='c',\
    mid='C1ae5f6786e5a8f65a865e865f68a5e68f5a86e8a',start=True)
    G.add_node(2,type='c',\
    mid='C2ae5f6786e5a8f65a865e865f68a5e68f5a86e8a',start=False)
    G.add_node(3,type='c',\
    mid='C3ae5f6786e5a8f65a865e865f68a5e68f5a86e8a',start=True)
    G.add_node(4,type='c',\
    mid='C4ae5f6786e5a8f65a865e865f68a5e68f5a86e8a',start=False)
    G.add_node(5,type='rf',mid='R1')
    G.add_node(6,type='pf',mid='R1')
    G.add_node(7,type='rr',mid='R1')
    G.add_node(8,type='pr',mid='R1')
    G.add_node(11,type='c',mid='C00004',start=False)
    G.add_node(12,type='c',\
    mid='Cf5dc8599a48d0111a3a5f618296752e1b53c8d30',start=False)
    G.add_node(13,type='c',\
    mid='Xf5dc8599a48d0111a3a5f618296752e1b53c8d30',start=False)

    G.add_path([1,5,6,2])

    H = G.copy()
    H.node[1]['start'] = False
    H.node[2]['start'] = True

    start_comp_ids = [
    'C2ae5f6786e5a8f65a865e865f68a5e68f5a86e8a',
    'C3ae5f6786e5a8f65a865e865f68a5e68f5a86e8a'
    ]

    update_start_compounds(G, start_comp_ids)

    assert G.nodes(data=True) == H.nodes(data=True)
    assert nx.is_isomorphic(G, H)

    H.node[1]['start'] = True
    H.node[2]['start'] = True
    H.node[3]['start'] = False

    update_start_compounds(G, ['Alpha', 'Beta', 'n-Beta'])

    assert G.nodes(data=True) == H.nodes(data=True)
    assert nx.is_isomorphic(G, H)

    H.node[1]['start'] = True
    H.node[2]['start'] = True
    H.node[3]['start'] = True
    H.node[4]['start'] = True

    update_start_compounds(G, \
    ['C67559', 'C67560', 'C2ae5f6786e5a8f65a865e865f68a5e68f5a86e8a', 'Twin'])

    assert G.nodes(data=True) == H.nodes(data=True)
    assert nx.is_isomorphic(G, H)

    for n in [1,2,3,4,11,12,13]:
        G.node[n]['start'] = False

    for n in [1,2,3,4]:
        H.node[n]['start'] = False
    for n in [11,12,13]:
        H.node[n]['start'] = True

    update_start_compounds(G, ['NAD'])

    assert G.nodes(data=True) == H.nodes(data=True)
    assert nx.is_isomorphic(G, H)

    for n in [1,2,3,4,11,12,13]:
        G.node[n]['start'] = False

    update_start_compounds(G, ['C00004','FAKE'])

    assert G.nodes(data=True) == H.nodes(data=True)
    assert nx.is_isomorphic(G, H)

    for n in [1,2,3,4,11,12,13]:
        G.node[n]['start'] = False

    update_start_compounds(G, ['NADX', 'C30184'])

    assert G.nodes(data=True) == H.nodes(data=True)
    assert nx.is_isomorphic(G, H)


def test_format_reaction_text():
    reactions = [
        {'Reactants':[[1,'C4']],
         'Products':[[2,'C1']]},
        {'Reactants':[[1,'C2'],[1,'C3']],
         'Products':[[1,'C5']]},
        {'Reactants':[[1,'C4'],[1,'C5']],
         'Products':[[1,'C6'],[1,'X0']]},
        {'Reactants':[[1,'X9'],[2,'C7']],
         'Products':[[1,'C6'],[1,'X8']]}
    ]
    exp_txt = [
        "C4 <=> 2 C1", "2 C1 <=> C4",
        "C2 + C3 <=> C5", "C5 <=> C2 + C3",
        "C4 + C5 <=> C6 + X0", "C6 + X0 <=> C4 + C5",
        "X9 + 2 C7 <=> C6 + X8", "C6 + X8 <=> X9 + 2 C7"
    ]
    for R in enumerate(reactions):
        assert format_reaction_text(R[1]) == exp_txt[R[0]*2]
        assert format_reaction_text(R[1], reverse=True) == exp_txt[R[0]*2+1]


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


def test_format_mdf_summary():
    # Set up testing "network"
    network = nx.DiGraph()
    network.graph['mine_data'] = {
        'R1':{'Reactants':[[1,'A']],'Products':[[2,'B']]},
        'R2':{'Reactants':[[1,'B']],'Products':[[1,'C']]},
        'R3':{'Reactants':[[1,'C']],'Products':[[2,'D']]},
        'R4':{'Reactants':[[2,'B']],'Products':[[1,'D'],[1,'E']]},
        'R5':{'Reactants':[[1,'D']],'Products':[[1,'B']]},
        'R6':{'Reactants':[[1,'D'],[1,'X']],'Products':[[1,'B']]}
    }

    # Set up testing pathways in text format
    pw_1 = "\n".join([
        'R1' + "\t" + format_reaction_text(network.graph['mine_data']['R1']),
        'R2' + "\t" + format_reaction_text(network.graph['mine_data']['R2']),
        'R3' + "\t" + format_reaction_text(network.graph['mine_data']['R3'])
    ])
    pw_2 = "\n".join([
        'R1' + "\t" + format_reaction_text(network.graph['mine_data']['R1']),
        'R4' + "\t" + format_reaction_text(network.graph['mine_data']['R4'])
    ])
    pw_3 = "\n".join([
        'R1' + "\t" + format_reaction_text(network.graph['mine_data']['R1']),
        'R5' + "\t" + format_reaction_text(network.graph['mine_data']['R5'], 1)
    ])
    pw_4 = "\n".join([
        'R1' + "\t" + format_reaction_text(network.graph['mine_data']['R1']),
        'R6' + "\t" + format_reaction_text(network.graph['mine_data']['R6'], 1)
    ])

    # Set up the testing MDF dictionary
    mdf_dict = {pw_1 : 5.0006, pw_2 : 5.0006, pw_3 : 4.3, pw_4 : None}

    # The expected output
    exp_summary_txt = "\n".join([
        "\t".join(x) for x in [
            ['pathway', 'MDF', 'length', 'reactions'],
            ['47afbbd9e7', '5.001', '2', 'R1,R4'],
            ['8d415da2b3', '5.001', '3', 'R1,R2,R3'],
            ['a3943b79d8', '4.300', '2', 'R1,-R5'],
            ['560b193f2a', 'NA', '2', 'R1,-R6']]
    ]) + "\n"

    assert exp_summary_txt == format_mdf_summary(mdf_dict, network)[1]


def test_format_pathway_html():
    # Set up testing "network"
    network = nx.DiGraph()
    network.graph['mine_data'] = {
        'R1':{'Reactants':[[1,'A']],
              'Products':[[2,'B']],
              'Operators':['M:1.2.3.a','1.2.3.54']},
        'R2':{'Reactants':[[1,'B']],
              'Products':[[1,'C']],
              'Operators':['3.2.1.-']},
        'R3':{'Reactants':[[1,'C']],
              'Products':[[2,'D']],
              'Operators':['2.2.2.22','M:2.3.-4.c']},
        'R4':{'Reactants':[[2,'B']],
              'Products':[[1,'D'],[1,'E']],
              'Operators':['M:4.1.-1.4c','4.1.1.1']},
        'R5':{'Reactants':[[1,'D']],
              'Products':[[1,'B']],
              'Operators':['1.3.1.78','1.3.1.79']},
        'RM6':{'Reactants':[[1,'D'],[1,'X']],
               'Products':[[1,'B']],
               'Operators':['M:3.1.1.b']},
        'A':{'Names':['cA','cA1']},
        'B':{'Names':['cB','cB1']},
        'C':{'Names':['cC']},
        'D':{'Names':['cD']},
        'E':{'Names':['cE','cE1','cE2']},
        'X':{'Names':['cX','cXy']}
    }
    network.graph['cmid2node'] = {'A':1, 'B':2, 'C':3, 'D':4, 'E':5, 'X':6}
    network.add_node(1, mid='A', type='c', start=True)
    network.add_node(2, mid='B', type='c', start=False)
    network.add_node(3, mid='C', type='c', start=False)
    network.add_node(4, mid='D', type='c', start=False)
    network.add_node(5, mid='E', type='c', start=False)
    network.add_node(6, mid='X', type='c', start=True)

    # Set up testing pathways in text format
    pw_1 = "\n".join([
        'R1' + "\t" + format_reaction_text(network.graph['mine_data']['R1']),
        'R2' + "\t" + format_reaction_text(network.graph['mine_data']['R2']),
        'R3' + "\t" + format_reaction_text(network.graph['mine_data']['R3'])
    ])
    pw_2 = "\n".join([
        'R1' + "\t" + format_reaction_text(network.graph['mine_data']['R1']),
        'R4' + "\t" + format_reaction_text(network.graph['mine_data']['R4'])
    ])
    pw_3 = "\n".join([
        'R1' + "\t" + format_reaction_text(network.graph['mine_data']['R1']),
        'R5' + "\t" + format_reaction_text(network.graph['mine_data']['R5'], 1)
    ])
    pw_4 = "\n".join([
        'R1' + "\t" + format_reaction_text(network.graph['mine_data']['R1']),
        'RM6' + "\t" + format_reaction_text(network.graph['mine_data']['RM6'], 1)
    ])

    # Set up the testing MDF dictionary
    mdf_dict = {pw_1 : 5.0006, pw_2 : 5.0006, pw_3 : 4.3, pw_4 : None}

    # Set up pathways dataframe
    pw_df, pw_txt = format_mdf_summary(mdf_dict, network)

    timestamp = time.strftime('%Y-%m-%d %H:%M:%S')

    exp_html = [
        '<!DOCTYPE html>',
        '<html>',
        '<head>',
        '  <link rel="stylesheet" href="style.css">',
        '<base target="_blank">',
        '</head>',
        '<body>',
        '',
        '<table border=0>',
        '<tr>',
        '  <td>',
        '  <h3><b>Report: 4 cD pathways</b></h3>',
        '  </td>',
        '  <td>',
        '    <h5><b>MDF</b></h5>',
        '  </td>',
        '  <td>',
        '    <h5><b>Pathway length</b></h5>',
        '  </td>',
        '</tr>',
        '<tr>',
        '  <td>',
        '  <li><h5><b>Report generated on:</b> %s</h5>' % timestamp,
        '  <li><h5><b>Target compound:</b> cD (D)</h5>',
        '  <li><h5><b>Reaction depth:</b> 3</h5>',
        '  <li><h5><b>Reaction limit:</b> 4</h5>',
        '  <li><h5><b>Number of pathways:</b> 4</h5>',
        '  <li><h5><b>Feasible/infeasible:</b> 3/1</h5>',
        '  </td>',
        '  <td>',
        '  <img src="mdf_summary.png" style="width: 200px;"/>',
        '  </td>',
        '  <td>',
        '  <img src="length_summary.png" style="width: 200px;"/>',
        '  </td>',
        '</tr>',
        '</table>',
        '',
        '<hr>',
        '',
        '<h5>Pathway "47afbbd9e7": 5.001 kJ/mol, 2 reactions</h5>',
        '',
        '<table border=0>',
        '<tr>',
        '  <td><a href="http://www.genome.jp/dbget-bin/www_bget?R1">R1</a></td>',
        '  <td style="background-color:#c7eae5"><a href="http://www.genome.jp/dbget-bin/www_bget?A">cA</a></td>',
        '  <td></td>',
        '  <td style="background-color:#f6e8c3"><a href="http://www.genome.jp/dbget-bin/www_bget?B">cB</a></td>',
        '</tr>',
        '<tr>',
        '  <td><p>',
        '    <a href="http://enzyme.expasy.org/EC/1.2.3.-">1.2.3.a*</a><br>',
        '    <a href="http://enzyme.expasy.org/EC/1.2.3.54">1.2.3.54</a><br>',
        '  </p></td>',
        '  <td><img src="cpd_png/A.png" style="width: 75px;"/></td>',
        '  <td><=></td>',
        '  <td>2 <img src="cpd_png/B.png" style="width: 75px;"/></td>',
        '</tr>',
        '</table>',
        '',
        '<table border=0>',
        '<tr>',
        '  <td><a href="http://www.genome.jp/dbget-bin/www_bget?R4">R4</a></td>',
        '  <td style="background-color:#f6e8c3"><a href="http://www.genome.jp/dbget-bin/www_bget?B">cB</a></td>',
        '  <td></td>',
        '  <td style="background-color:#c2a5cf"><a href="http://www.genome.jp/dbget-bin/www_bget?D">cD</a></td>',
        '  <td></td>',
        '  <td style="background-color:#f6e8c3"><a href="http://www.genome.jp/dbget-bin/www_bget?E">cE</a></td>',
        '</tr>',
        '<tr>',
        '  <td><p>',
        '    <a href="http://enzyme.expasy.org/EC/4.1.1.-">4.1.-1.4c*</a><br>',
        '    <a href="http://enzyme.expasy.org/EC/4.1.1.1">4.1.1.1</a><br>',
        '  </p></td>',
        '  <td>2 <img src="cpd_png/B.png" style="width: 75px;"/></td>',
        '  <td><=></td>',
        '  <td><img src="cpd_png/D.png" style="width: 75px;"/></td>',
        '  <td>+</td>',
        '  <td><img src="cpd_png/E.png" style="width: 75px;"/></td>',
        '</tr>',
        '</table>',
        '',
        '<hr>',
        '',
        '<h5>Pathway "8d415da2b3": 5.001 kJ/mol, 3 reactions</h5>',
        '',
        '<table border=0>',
        '<tr>',
        '  <td><a href="http://www.genome.jp/dbget-bin/www_bget?R1">R1</a></td>',
        '  <td style="background-color:#c7eae5"><a href="http://www.genome.jp/dbget-bin/www_bget?A">cA</a></td>',
        '  <td></td>',
        '  <td style="background-color:#f6e8c3"><a href="http://www.genome.jp/dbget-bin/www_bget?B">cB</a></td>',
        '</tr>',
        '<tr>',
        '  <td><p>',
        '    <a href="http://enzyme.expasy.org/EC/1.2.3.-">1.2.3.a*</a><br>',
        '    <a href="http://enzyme.expasy.org/EC/1.2.3.54">1.2.3.54</a><br>',
        '  </p></td>',
        '  <td><img src="cpd_png/A.png" style="width: 75px;"/></td>',
        '  <td><=></td>',
        '  <td>2 <img src="cpd_png/B.png" style="width: 75px;"/></td>',
        '</tr>',
        '</table>',
        '',
        '<table border=0>',
        '<tr>',
        '  <td><a href="http://www.genome.jp/dbget-bin/www_bget?R2">R2</a></td>',
        '  <td style="background-color:#f6e8c3"><a href="http://www.genome.jp/dbget-bin/www_bget?B">cB</a></td>',
        '  <td></td>',
        '  <td style="background-color:#f6e8c3"><a href="http://www.genome.jp/dbget-bin/www_bget?C">cC</a></td>',
        '</tr>',
        '<tr>',
        '  <td><p>',
        '    <a href="http://enzyme.expasy.org/EC/3.2.1.-">3.2.1.-</a><br>',
        '  </p></td>',
        '  <td><img src="cpd_png/B.png" style="width: 75px;"/></td>',
        '  <td><=></td>',
        '  <td><img src="cpd_png/C.png" style="width: 75px;"/></td>',
        '</tr>',
        '</table>',
        '',
        '<table border=0>',
        '<tr>',
        '  <td><a href="http://www.genome.jp/dbget-bin/www_bget?R3">R3</a></td>',
        '  <td style="background-color:#f6e8c3"><a href="http://www.genome.jp/dbget-bin/www_bget?C">cC</a></td>',
        '  <td></td>',
        '  <td style="background-color:#c2a5cf"><a href="http://www.genome.jp/dbget-bin/www_bget?D">cD</a></td>',
        '</tr>',
        '<tr>',
        '  <td><p>',
        '    <a href="http://enzyme.expasy.org/EC/2.2.2.22">2.2.2.22</a><br>',
        '    <a href="http://enzyme.expasy.org/EC/2.3.4.-">2.3.-4.c*</a><br>',
        '  </p></td>',
        '  <td><img src="cpd_png/C.png" style="width: 75px;"/></td>',
        '  <td><=></td>',
        '  <td>2 <img src="cpd_png/D.png" style="width: 75px;"/></td>',
        '</tr>',
        '</table>',
        '',
        '<hr>',
        '',
        '<h5>Pathway "a3943b79d8": 4.300 kJ/mol, 2 reactions</h5>',
        '',
        '<table border=0>',
        '<tr>',
        '  <td><a href="http://www.genome.jp/dbget-bin/www_bget?R1">R1</a></td>',
        '  <td style="background-color:#c7eae5"><a href="http://www.genome.jp/dbget-bin/www_bget?A">cA</a></td>',
        '  <td></td>',
        '  <td style="background-color:#f6e8c3"><a href="http://www.genome.jp/dbget-bin/www_bget?B">cB</a></td>',
        '</tr>',
        '<tr>',
        '  <td><p>',
        '    <a href="http://enzyme.expasy.org/EC/1.2.3.-">1.2.3.a*</a><br>',
        '    <a href="http://enzyme.expasy.org/EC/1.2.3.54">1.2.3.54</a><br>',
        '  </p></td>',
        '  <td><img src="cpd_png/A.png" style="width: 75px;"/></td>',
        '  <td><=></td>',
        '  <td>2 <img src="cpd_png/B.png" style="width: 75px;"/></td>',
        '</tr>',
        '</table>',
        '',
        '<table border=0>',
        '<tr>',
        '  <td><a href="http://www.genome.jp/dbget-bin/www_bget?R5">R5</a></td>',
        '  <td style="background-color:#f6e8c3"><a href="http://www.genome.jp/dbget-bin/www_bget?B">cB</a></td>',
        '  <td></td>',
        '  <td style="background-color:#c2a5cf"><a href="http://www.genome.jp/dbget-bin/www_bget?D">cD</a></td>',
        '</tr>',
        '<tr>',
        '  <td><p>',
        '    <a href="http://enzyme.expasy.org/EC/1.3.1.78">1.3.1.78</a><br>',
        '    <a href="http://enzyme.expasy.org/EC/1.3.1.79">1.3.1.79</a><br>',
        '  </p></td>',
        '  <td><img src="cpd_png/B.png" style="width: 75px;"/></td>',
        '  <td><=></td>',
        '  <td><img src="cpd_png/D.png" style="width: 75px;"/></td>',
        '</tr>',
        '</table>',
        '',
        '<hr>',
        '',
        '<h5>Pathway "560b193f2a": MDF FAILED, 2 reactions</h5>',
        '',
        '<table border=0>',
        '<tr>',
        '  <td><a href="http://www.genome.jp/dbget-bin/www_bget?R1">R1</a></td>',
        '  <td style="background-color:#c7eae5"><a href="http://www.genome.jp/dbget-bin/www_bget?A">cA</a></td>',
        '  <td></td>',
        '  <td style="background-color:#f6e8c3"><a href="http://www.genome.jp/dbget-bin/www_bget?B">cB</a></td>',
        '</tr>',
        '<tr>',
        '  <td><p>',
        '    <a href="http://enzyme.expasy.org/EC/1.2.3.-">1.2.3.a*</a><br>',
        '    <a href="http://enzyme.expasy.org/EC/1.2.3.54">1.2.3.54</a><br>',
        '  </p></td>',
        '  <td><img src="cpd_png/A.png" style="width: 75px;"/></td>',
        '  <td><=></td>',
        '  <td>2 <img src="cpd_png/B.png" style="width: 75px;"/></td>',
        '</tr>',
        '</table>',
        '',
        '<table border=0>',
        '<tr>',
        '  <td>RM6</td>',
        '  <td style="background-color:#f6e8c3"><a href="http://www.genome.jp/dbget-bin/www_bget?B">cB</a></td>',
        '  <td></td>',
        '  <td style="background-color:#c2a5cf"><a href="http://www.genome.jp/dbget-bin/www_bget?D">cD</a></td>',
        '  <td></td>',
        '  <td style="background-color:#c7eae5"><a href="http://www.genome.jp/dbget-bin/www_bget?X">cX</a></td>',
        '</tr>',
        '<tr>',
        '  <td><p>',
        '    <a href="http://enzyme.expasy.org/EC/3.1.1.-">3.1.1.b*</a><br>',
        '  </p></td>',
        '  <td><img src="cpd_png/B.png" style="width: 75px;"/></td>',
        '  <td><=></td>',
        '  <td><img src="cpd_png/D.png" style="width: 75px;"/></td>',
        '  <td>+</td>',
        '  <td><img src="cpd_png/X.png" style="width: 75px;"/></td>',
        '</tr>',
        '</table>',
        ''
    ]

    auto_html = format_pathway_html(pw_df, network, 4, 3, 4)

    for line in enumerate(auto_html.split("\n")):
        print(line[1])
        print(exp_html[line[0]])
        assert line[1] == exp_html[line[0]]


def test_disconnect_reactants_products():

    # Set up testing network
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

    # Add dictionaries
    N.graph['cmid2node'] = {
    'Cf647c96ae2e66c3b6ab160faa1d8498be5112fe4':1,
    'C31890':2,
    'C00291':3,
    'C6a6f4d5234ea2b14b42c391eb760d6311afa8388':4,
    'C67559':5,
    'C67560':6,
    'C99999':7
    }

    N.graph['kegg2nodes'] = {
    'C31890':set([2]), 'C00291':set([3]), 'C67559':set([5]),
    'C67560':set([6]), 'C99999':set([7])
    }

    N.graph['name2nodes'] = {
    'Alpha':set([1]),
    'Beta':set([2]),
    'Gamma':set([3]),
    'Twin':set([4,5]),
    'n-Alpha':set([6]),
    'n-Beta':set([7])
    }

    # Select compounds and reaction nodes to ban
    ban_reacs = {'Alpha','C31890'}
    ban_prods = {'Twin','Cf647c96ae2e66c3b6ab160faa1d8498be5112fe4','C99999'}

    ban_reacs_n = set()
    for node in N.nodes():
        if N.node[node]['type'] in {'rr', 'rf'}:
            if N.node[node]['c'].intersection({1, 2}):
                ban_reacs_n.add(node)

    ban_prods_n = set()
    for node in N.nodes():
        if N.node[node]['type'] in {'pr', 'pf'}:
            if N.node[node]['c'].intersection({4, 5, 1, 7}):
                ban_prods_n.add(node)

    assert ban_reacs_n == {201, 101, 103}
    assert ban_prods_n == {104, 202, 302, 404, 502, 504}

    # Make copy of network and ban in three different ways
    G = N.copy()
    disconnect_reactants_products(G, reactants=ban_reacs)
    assert not nx.is_isomorphic(G, N)
    assert set(N.nodes()) - ban_reacs_n == set(G.nodes())

    G = N.copy()
    disconnect_reactants_products(G, products=ban_prods)
    assert not nx.is_isomorphic(G, N)
    assert set(N.nodes()) - ban_prods_n == set(G.nodes())

    G = N.copy()
    disconnect_reactants_products(G, ban_reacs, ban_prods)
    assert not nx.is_isomorphic(G, N)
    assert set(N.nodes()) - ban_reacs_n - ban_prods_n == set(G.nodes())

    # Make sure that supplying empty sets does not affect the network
    G = N.copy()
    disconnect_reactants_products(G, set(), set())
    assert nx.is_isomorphic(G, N)
