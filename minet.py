#!/usr/bin/env python3

# Import modules
import networkx as nx
import MineClient3 as mc
import re
import sys
import argparse
import pickle

# Define functions
def ReadCompounds(filename):
    """Read a file with KEGG compound IDs."""
    sys.stdout.write("Reading compound ID file...")
    sys.stdout.flush()
    compounds = [line.rstrip() for line in open(filename, 'r')]
    for c in compounds:
        if re.fullmatch("^C[0-9]{5}$", c) == None:
            msg = "Warning: The supplied string '", c, "' is not a valid KEGG compound ID."
            sys.exit(msg)
    print(" Done.")
    return compounds

def test_ReadCompounds():

    import pytest
    import tempfile

    t1 = str.encode("C10000\nC40055\nC13482\n")
    t2 = str.encode("C13854\nR10309\nC33190\n")

    f1 = tempfile.NamedTemporaryFile()
    f1.write(t1)
    f1.flush()
    f2 = tempfile.NamedTemporaryFile()
    f2.write(t2)
    f2.flush()

    assert set(ReadCompounds(f1.name)) == set(['C10000','C40055','C13482'])
    with pytest.raises(SystemExit): ReadCompounds(f2.name)

    f1.close()
    f2.close()


def KeggToMineId(kegg_ids):
    """Translate KEGG IDs to MINE IDs."""
    sys.stdout.write("Translating from KEGG IDs to MINE IDs...")
    sys.stdout.flush()
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    kegg_id_dict = {}
    for kegg_id in kegg_ids:
        try:
            kegg_id_dict[kegg_id] = con.quick_search("KEGGexp2", kegg_id)[0]['_id']
        except:
            continue
    print(" Done.")
    return kegg_id_dict

def test_KeggToMineId():
    assert KeggToMineId(['C15667', 'C16519', 'C00130']) == {'C15667':'C023e725c27a385fd9057c8171ca1992a32a3e9a4',
    'C16519':'Cfa3c20a209e12ac4a70f10530cf43bea2b4422fe',
    'C00130':'Cae5be0a0a2bb6b6baef29e4d4c6b3f4c1a67ad19'}
    # C00003 is not present in the database
    assert KeggToMineId(['C00003', 'C16519', 'C00130']) == {'C16519':'Cfa3c20a209e12ac4a70f10530cf43bea2b4422fe',
    'C00130':'Cae5be0a0a2bb6b6baef29e4d4c6b3f4c1a67ad19'}


def GetRawNetwork(comp_id_list, step_limit=10, comp_limit=100000):
    """Download connected reactions and compounds up to the limits."""

    sys.stdout.write("Downloading raw network data...")
    sys.stdout.flush()

    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    # Set up output dictionaries
    comp_dict = {}
    rxn_dict = {}

    # Set up counters
    steps = 0
    comps = 0

    # First add the starting compounds
    for comp_id in comp_id_list:
        comps += 1
        comp = con.get_comps(db, [comp_id])[0]
        comp_dict[comp_id] = comp # Add compound to dict

    extended_comp_ids = set()

    # Perform stepwise expansion of downloaded data
    while steps < step_limit and comps < comp_limit:
        steps += 1
        unextended_comp_ids = set(comp_dict.keys()) - extended_comp_ids
        for comp_id in unextended_comp_ids:
            comp = comp_dict[comp_id] # New compounds are always in the dictionary
            rxn_id_list = []
            try:
                rxn_id_list.extend(comp['Reactant_in'])
            except KeyError:
                pass
            try:
                rxn_id_list.extend(comp['Product_of'])
            except KeyError:
                pass
            for rxn_id in rxn_id_list:
                # Only download new reactions
                try:
                    rxn = rxn_dict[rxn_id]
                except KeyError:
                    rxn = con.get_rxns(db, [rxn_id])[0]
                    rxn_dict[rxn_id] = rxn # Add new reaction
                rxn_comp_ids = [x[1] for x in rxn['Products']] + [x[1] for x in rxn['Reactants']]
                for rxn_comp_id in rxn_comp_ids:
                    # Only download new compounds
                    try:
                        rxn_comp = comp_dict[rxn_comp_id]
                    except KeyError:
                        rxn_comp = con.get_comps(db, [rxn_comp_id])[0]
                        comp_dict[rxn_comp_id] = rxn_comp # Add new compound
                        comps += 1
            extended_comp_ids.add(comp_id)
    print(" Done.")
    return (comp_dict, rxn_dict)

def test_GetRawNetwork():

    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    db = "KEGGexp2"
    con = mc.mineDatabaseServices(server_url)

    rxn_id_list = ['R04759e864c86cfd0eaeb079404d5f18dae6c7227', 'Re598257045ae3ce45dabf450b57708d84e642558']

    c_id_1 = [a for b in [[y[1] for y in rxn['Reactants']] + [y[1] for y in rxn['Products']] for rxn in con.get_rxns(db,rxn_id_list)] for a in b]

    for comp in con.get_comps(db, c_id_1):
        try:
            rxn_id_list.extend(comp['Reactant_in'])
        except KeyError:
            pass
        try:
            rxn_id_list.extend(comp['Product_of'])
        except KeyError:
            pass

    c_id_2 = [a for b in [[y[1] for y in rxn['Reactants']] + [y[1] for y in rxn['Products']] for rxn in con.get_rxns(db,rxn_id_list)] for a in b]

    rxn_id_list = list(set(rxn_id_list))
    c_id_2 = list(set(c_id_2))

    a = con.get_comps(db, c_id_2[0:50])
    b = con.get_comps(db, c_id_2[50:100])
    c = con.get_comps(db, c_id_2[100:])
    comps = a + b + c

    comp_dict = dict([(comp['_id'], comp) for comp in comps])
    rxn_dict = dict([(rxn['_id'], rxn) for rxn in con.get_rxns(db, rxn_id_list)])

    compound_ids = ['Cefbaa83ea06e7c31820f93c1a5535e1378aba42b','C38a97a9f962a32b984b1702e07a25413299569ab']

    assert GetRawNetwork(compound_ids, step_limit=2) == (comp_dict, rxn_dict)


def AddQuadReactionNode(graph, rxn):
    """
    Adds a "Quad Reaction Node" (QRN) group of nodes to a graph.

    The QRN consists of two nodes constituting the intended forward direction
    of the reaction and two nodes constituting the reverse direction. Each pair
    of nodes is connected by an edge in the direction of the reaction. Each node
    represents a group of compounds on one side of the reaction equation.
    """

    reactants_f = set([x[1] for x in rxn['Reactants']])
    products_f = set([x[1] for x in rxn['Products']])
    reactants_r = products_f
    products_r = reactants_f

    graph.add_node(('rf', rxn['_id']), data=rxn, c=reactants_f)
    graph.add_node(('pf', rxn['_id']), data=rxn, c=products_f)
    graph.add_edge(('rf', rxn['_id']), ('pf', rxn['_id']), weight=1)

    graph.add_node(('rr', rxn['_id']), data=rxn, c=reactants_r)
    graph.add_node(('pr', rxn['_id']), data=rxn, c=products_r)
    graph.add_edge(('rr', rxn['_id']), ('pr', rxn['_id']), weight=1)

    return graph

def test_AddQuadReactionNode():
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    db = "KEGGexp2"
    con = mc.mineDatabaseServices(server_url)

    c = 'C38a97a9f962a32b984b1702e07a25413299569ab'
    rxn = con.get_rxns(db, ['R04759e864c86cfd0eaeb079404d5f18dae6c7227'])[0]

    reactants = set(['Caf6fc55862387e5fd7cd9635ef9981da7f08a531', 'X25a9fafebc1b08a0ae0fec015803771c73485a61'])
    products = set(['Cefbaa83ea06e7c31820f93c1a5535e1378aba42b', 'Xf729c487f9b991ec6f645c756cf34b9a20b9e8a4'])

    G1 = nx.DiGraph()
    G1.add_node(('c', c), data=con.get_comps(db, [c])[0])

    rxn_id = 'R04759e864c86cfd0eaeb079404d5f18dae6c7227'

    G1.add_node(('rf', rxn_id), data=rxn, c=reactants) # Forward (intended) direction reactants
    G1.add_node(('pf', rxn_id), data=rxn, c=products) # Forward (intended) direction products
    G1.add_node(('rr', rxn_id), data=rxn, c=products) # Reverse direction reactants
    G1.add_node(('pr', rxn_id), data=rxn, c=reactants) # Reverse direction products
    G1.add_edge(('rf', rxn_id), ('pf', rxn_id), weight=1) # Directed edge for the forward reaction
    G1.add_edge(('rr', rxn_id), ('pr', rxn_id), weight=1) # Directed edge for the reverse reaction

    G2 = nx.DiGraph()
    G2.add_node(('c', c), data=con.get_comps(db, [c])[0])

    assert nx.is_isomorphic(AddQuadReactionNode(G2, rxn), G1)


def AddCompoundNode(graph, compound):
    """Adds a compound node to the graph."""
    node_data = compound
    node = ('c', compound['_id'])
    graph.add_node(node, data=node_data)
    return graph

def test_AddCompoundNode():
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    db = "KEGGexp2"
    con = mc.mineDatabaseServices(server_url)

    G1 = nx.DiGraph()
    G2 = nx.DiGraph()
    comp1 = con.get_comps(db, ['Xf5dc8599a48d0111a3a5f618296752e1b53c8d30'])[0]
    comp2 = con.get_comps(db, ['C38a97a9f962a32b984b1702e07a25413299569ab'])[0]

    G2.add_node(('c', comp1['_id']), data=comp1)
    G2.add_node(('c', comp2['_id']), data=comp2)

    assert nx.is_isomorphic(AddCompoundNode(AddCompoundNode(G1, comp1), comp2), G2)


def ConstructNetwork(comp_dict, rxn_dict, start_comp_ids=[]):
    """Constructs a directed graph (network) from the compound and reaction
    dictionaries produced by GetRawNetwork."""

    sys.stdout.write("Constructing network...")
    sys.stdout.flush()

    start_comp_ids = set(start_comp_ids)

    def CheckConnection(minetwork, c_node, r_node):
        """Checks that the compound-to-reaction node connection is valid."""
        if c_node[1] not in minetwork.node[r_node]['c']:
            msg_dict = {
            'rf':'forward reactants',
            'pf':'forward products',
            'rr':'reverse reactants',
            'pr':'reverse products'
            }
            set_name = msg_dict[r_node[0]]
            node_set = ", ".join(minetwork.node[r_node]['c'])
            message = "Warning: Compound %s is not found in the %s of reaction %s (%s). Connection not created.\n" % (c_node[1], set_name, r_node[1], node_set)
            sys.stderr.write(message)
            return False
        else:
            return True

    # Initialise directed graph
    minetwork = nx.DiGraph()

    # Add all compounds
    for comp in comp_dict.values():
        minetwork = AddCompoundNode(minetwork, comp)
        # Add an attribute to the compound specifying whether it is a starting compound
        if comp['_id'] in start_comp_ids:
            minetwork.node[('c',comp['_id'])]['start'] = True
        else:
            minetwork.node[('c',comp['_id'])]['start'] = False

    # Add all reactions
    for rxn in rxn_dict.values():
        minetwork = AddQuadReactionNode(minetwork, rxn)

    # Iterate over compounds, connecting them to the correct reaction nodes
    for comp in comp_dict.values():
        c_node = ('c', comp['_id']) # The compound node to connect
        try:
            # Connect to reactions in which the compound is a reactant
            rxn_ids = comp['Reactant_in']
            # Only connect reactions that are present in the network/dictionary
            rxn_ids = set.intersection(set(rxn_dict.keys()), set(rxn_ids))
            for rxn_id in rxn_ids:
                rf_node = ('rf', rxn_id)
                if CheckConnection(minetwork, c_node, rf_node):
                    minetwork.add_edge(c_node, rf_node, weight=0) # Connect c -> forward reactants
                pr_node = ('pr', rxn_id)
                if CheckConnection(minetwork, c_node, pr_node):
                    minetwork.add_edge(pr_node, c_node, weight=0) # Connect reverse products -> c
        except KeyError:
            pass
        try:
            # Connect to reactions in which the compound is a product
            rxn_ids = comp['Product_of']
            # Only connect reactions that are present in the network/dictionary
            rxn_ids = set.intersection(set(rxn_dict.keys()), set(rxn_ids))
            for rxn_id in rxn_ids:
                pf_node = ('pf', rxn_id)
                if CheckConnection(minetwork, c_node, pf_node):
                    minetwork.add_edge(pf_node, c_node, weight=0) # Connect c -> forward products
                rr_node = ('rr', rxn_id)
                if CheckConnection(minetwork, c_node, rr_node):
                    minetwork.add_edge(c_node, rr_node, weight=0) # Connect reverse reactants -> c
        except KeyError:
            pass
    print(" Done.")
    return minetwork

def test_ConstructNetwork(capsys):
    comp_dict = {}
    comp_dict['C1'] = {'_id':'C1', 'Reactant_in':['R99']}
    comp_dict['C2'] = {'_id':'C2', 'Reactant_in':['R1e'], 'Product_of':['R99']}
    comp_dict['C3'] = {'_id':'C3', 'Reactant_in':['Rcd'], 'Product_of':['R99','Rc3']}
    comp_dict['C4'] = {'_id':'C4', 'Product_of':['R1e']}
    comp_dict['C5'] = {'_id':'C5', 'Reactant_in':['Rc3','R2f'], 'Product_of':['Rcd','R1e']}
    comp_dict['C6'] = {'_id':'C6', 'Product_of':['Rcd']}
    comp_dict['C7'] = {'_id':'C7', 'Product_of':['R2f','R7f']} # Seeding with non-expanded reaction R7f
    comp_dict['C8'] = {'_id':'C8', 'Reactant_in':['Rb7'], 'Product_of':['R2f']} # Seeding with non-expanded reaction Rb7

    rxn_dict = {}
    rxn_dict['R99'] = {'_id':'R99', 'Products':[[1,'C2'],[1,'C3']], 'Reactants':[[1,'C1']]}
    rxn_dict['R1e'] = {'_id':'R1e', 'Products':[[1,'C4'],[1,'C5']], 'Reactants':[[1,'C2']]}
    rxn_dict['Rc3'] = {'_id':'Rc3', 'Products':[[1,'C3']], 'Reactants':[[1,'C5']]}
    rxn_dict['Rcd'] = {'_id':'Rcd', 'Products':[[1,'C5'],[1,'C6']], 'Reactants':[[1,'C3']]}
    rxn_dict['R2f'] = {'_id':'R2f', 'Products':[[1,'C7'],[1,'C8']], 'Reactants':[[1,'C5']]}

    G = nx.DiGraph()

    for comp in comp_dict.values():
        G = AddCompoundNode(G, comp)

    for rxn in rxn_dict.values():
        G = AddQuadReactionNode(G, rxn)

    for node in G.nodes():
        if node[0] == 'c':
            G.node[node]['start'] = False
    G.node[('c','C1')]['start'] = True # C1 is the starting compound

    # C1 edges
    c = ('c','C1')
    G.add_edge(c, ('rf','R99'), weight=0)
    G.add_edge(('pr','R99'), c, weight=0)

    # C2 edges
    c = ('c','C2')
    G.add_edge(('pf','R99'), c, weight=0)
    G.add_edge(c, ('rr','R99'), weight=0)
    G.add_edge(c, ('rf','R1e'), weight=0)
    G.add_edge(('pr','R1e'), c, weight=0)

    # C3 edges
    c = ('c','C3')
    G.add_edge(('pf','R99'), c, weight=0)
    G.add_edge(c, ('rr','R99'), weight=0)
    G.add_edge(c, ('rr','Rc3'), weight=0)
    G.add_edge(('pf','Rc3'), c, weight=0)
    G.add_edge(c, ('rf','Rcd'), weight=0)
    G.add_edge(('pr','Rcd'), c, weight=0)

    # C4 edges
    c = ('c','C4')
    G.add_edge(('pf','R1e'), c, weight=0)
    G.add_edge(c, ('rr','R1e'), weight=0)

    # C5 edges
    c = ('c','C5')
    G.add_edge(('pf','R1e'), c, weight=0)
    G.add_edge(c, ('rr','R1e'), weight=0)
    G.add_edge(('pr','Rc3'), c, weight=0)
    G.add_edge(c, ('rf','Rc3'), weight=0)
    G.add_edge(('pf','Rcd'), c, weight=0)
    G.add_edge(c, ('rr','Rcd'), weight=0)
    G.add_edge(c, ('rf','R2f'), weight=0)
    G.add_edge(('pr','R2f'), c, weight=0)
    # C7 edges

    # C6 edges
    c = ('c','C6')
    G.add_edge(('pf','Rcd'), c, weight=0)
    G.add_edge(c, ('rr','Rcd'), weight=0)

    c = ('c','C7')
    G.add_edge(('pf','R2f'), c, weight=0)
    G.add_edge(c, ('rr','R2f'), weight=0)

    # C8 edges
    c = ('c','C8')
    G.add_edge(('pf','R2f'), c, weight=0)
    G.add_edge(c, ('rr','R2f'), weight=0)

    assert nx.is_isomorphic(ConstructNetwork(comp_dict,rxn_dict,['C1']), G)

    # Test contents node by node
    t = True
    for node in ConstructNetwork(comp_dict,rxn_dict,['C1']).nodes(data=True):
        if G.node[node[0]] != node[1]:
            t = False
            break
    assert t

    # Test output of CheckConnection
    comp_dict = {'c':{'_id':'c', 'Reactant_in':['r']}}
    rxn_dict = {'r':{'_id':'r', 'Products':[[1,'a'],[1,'b']], 'Reactants':[[1,'z']]}}
    ConstructNetwork(comp_dict, rxn_dict, ['C1'])
    out, err = capsys.readouterr()
    assert err == """Warning: Compound c is not found in the forward reactants of reaction r (z). Connection not created.\nWarning: Compound c is not found in the reverse products of reaction r (z). Connection not created.\n"""


# Main code block
def main(infile_name, step_limit, comp_limit, outfile_name):
    # Get starting compound MINE IDs
    start_ids = list(set(KeggToMineId(ReadCompounds(infile_name)).values()))
    # Create the network
    minetwork = ConstructNetwork(*GetRawNetwork(start_ids, step_limit, comp_limit), start_ids)
    # Save to Pickle
    pickle.dump(minetwork, open(outfile_name, 'wb'))


if __name__ == "__main__":
    # Read arguments from the commandline
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='Read KEGG compound identifiers from text file.')
    parser.add_argument('outfile', help='Write MINE network to Python Pickle file.')
    parser.add_argument('-r', type=int, default=10, help='Maximum number of reaction steps to download.')
    parser.add_argument('-c', type=int, default=100000, help='Maximum number of compounds to download.')
    args = parser.parse_args()
    main(args.infile, args.r, args.c, args.outfile)
