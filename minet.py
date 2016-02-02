#!/usr/bin/env python3

# Import modules
import networkx as nx
import MineClient3 as mc
import re
import sys

# Define functions
def ReadCompounds(filename):
    """Read a file with KEGG compound IDs."""
    compounds = [line.rstrip() for line in open(filename, 'r')]
    for c in compounds:
        if re.fullmatch("^C[0-9]{5}$", c) == None:
            msg = "Warning: The supplied string '", c, "' is not a valid KEGG compound ID."
            sys.exit(msg)
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
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    kegg_id_dict = {}
    for kegg_id in kegg_ids:
        try:
            kegg_id_dict[kegg_id] = con.quick_search("KEGGexp2", kegg_id)[0]['_id']
        except:
            continue
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
            print(comp_id)
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
                print(rxn_id)
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
                        print(rxn_comp_id)
            extended_comp_ids.add(comp_id)
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


def ConstructNetwork(comp_dict, rxn_dict):
    # Code here
    return minetwork

def test_ConstructNetwork():
    comp_dict = {}
    comp_dict['C1'] = {'_id':'C1', 'Reactant_in':['R99']}
    comp_dict['C2'] = {'_id':'C2', 'Reactant_in':['R1e'], 'Product_of':['R99']}
    comp_dict['C3'] = {'_id':'C3', 'Reactant_in':['Rcd'], 'Product_of':['R99','Rc3']}
    comp_dict['C4'] = {'_id':'C4', 'Product_of':['R1e']}
    comp_dict['C5'] = {'_id':'C5', 'Reactant_in':['Rc3','R2f'], 'Product_of':['Rcd','R1e']}
    comp_dict['C6'] = {'_id':'C6', 'Product_of':['Rcd']}
    comp_dict['C7'] = {'_id':'C7', 'Product_of':['R2f']}
    comp_dict['C8'] = {'_id':'C8', 'Product_of':['R2f']}

    rxn_dict = {}
    rxn_dict['R99'] = {'_id':'R99', 'Products':[[1,'C2'],[1,'C3']], 'Reactants':[[1,'C1']]}
    rxn_dict['R1e'] = {'_id':'R1e', 'Products':[[1,'C4'],[1,'C5']], 'Reactants':[[1,'C2']]}
    rxn_dict['Rc3'] = {'_id':'Rc3', 'Products':[[1,'C3']], 'Reactants':[[1,'C5']]}
    rxn_dict['Rcd'] = {'_id':'Rcd', 'Products':[[1,'C5'],[1,'C6']], 'Reactants':[[1,'C3']]}
    rxn_dict['R2f'] = {'_id':'R2f', 'Products':[[1,'C7'],[1,'C8']], 'Reactants':[[1,'C5']]}

    G = nx.DiGraph()


    assert ConstructNetwork({},{}) == G


# Main code block
def main():
    # Code here
    pass

def test_main():
    assert main() == ""

if __name__ == "__main__":
    main(options)
