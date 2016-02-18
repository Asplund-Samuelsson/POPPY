#!/usr/bin/env python3

# Import modules
import networkx as nx
import MineClient3 as mc
import re
import sys
import argparse
import pickle
import time

# Define functions
def QuickSearch(con, db, query):
    """Wrapper for MineClient3 quick_search() with reconnect functionality."""
    n = 0
    results = []
    while True:
        try:
            results = con.quick_search(db, query)
            return results
        except mc.ServerError:
            return results
        except:
            # Server not responding, try again
            n += 1
            if n % 5 == 0:
                message = "Warning: Server not responding after %s attempts ('%s').\n" % (str(n), query)
                sys.stderr.write(message)
                sys.stderr.flush()
            if n >= 36:
                message = "Warning: Connection attempt limit reached. Returning empty list.\n"
                sys.stderr.write(message)
                sys.stderr.flush()
                return results
            if n <= 12:
                time.sleep(10)
            if n > 12:
                time.sleep(30)

def test_QuickSearch():

    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    assert QuickSearch(con, db, 'C00022')[0]['Names'][0] == 'Pyruvate'
    assert QuickSearch(con, db, 'random_query') == []


def GetComps(con, db, comp_id_list):
    """Wrapper for MineClient3 get_comps() with reconnect functionality."""
    n = 0
    results = []
    while True:
        try:
            results = con.get_comps(db, comp_id_list)
            return results
        except mc.ServerError:
            return results
        except:
            # Server not responding, try again
            n += 1
            if n % 5 == 0:
                message = "Warning: Server not responding after %s attempts ('%s').\n" % (str(n), str(comp_id_list))
                sys.stderr.write(message)
                sys.stderr.flush()
            if n >= 36:
                message = "Warning: Connection attempt limit reached. Returning empty list.\n"
                sys.stderr.write(message)
                sys.stderr.flush()
                return results
            if n <= 12:
                time.sleep(10)
            if n > 12:
                time.sleep(30)

def test_GetComps():
    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    assert GetComps(con, db, ['Cc93137cc81324a5b2872b0bf1c77866c234d66e1'])[0]['Formula'] == 'C7H15O10P'
    assert GetComps(con, db, ['Cc93137cc81324a5b2872b0bf1c77866c234d66e1'])[0]['dG_error'] == 1.02079
    assert GetComps(con, db, ['not_a_comp_id']) == [None]


def GetRxns(con, db, rxn_id_list):
    """Wrapper for MineClient3 get_rxns() with reconnect functionality."""
    n = 0
    results = []
    while True:
        try:
            results = con.get_rxns(db, rxn_id_list)
            return results
        except mc.ServerError:
            return results
        except:
            # Server not responding, try again
            n += 1
            if n % 5 == 0:
                message = "Warning: Server not responding after %s attempts ('%s').\n" % (str(n), rxn_id_list)
                sys.stderr.write(message)
                sys.stderr.flush()
            if n >= 36:
                message = "Warning: Connection attempt limit reached. Returning empty list.\n"
                sys.stderr.write(message)
                sys.stderr.flush()
                return results
            if n <= 12:
                time.sleep(10)
            if n > 12:
                time.sleep(30)

def test_GetRxns():

    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    rxn_id = 'Re598257045ae3ce45dabf450b57708d84e642558'
    rxn_op = '1.14.13.e'
    rxn_rlen = 4

    assert len(GetRxns(con, db, [rxn_id])) == 1
    assert GetRxns(con, db, [rxn_id])[0]['Operators'] == [rxn_op]
    assert len(GetRxns(con, db, [rxn_id])[0]['Reactants']) == 4
    assert GetRxns(con, db, ['random_reaction']) == [None]


def ReadCompounds(filename):
    """Read a file with KEGG compound IDs."""
    sys.stdout.write("\nReading compound ID file...\n")
    sys.stdout.flush()
    compounds = [line.rstrip() for line in open(filename, 'r')]
    for c in compounds:
        if re.fullmatch("^C[0-9]{5}$", c) == None:
            msg = "Warning: The supplied string '", c, "' is not a valid KEGG compound ID."
            sys.exit(msg)
    print("Done.")
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
    sys.stdout.write("\nTranslating from KEGG IDs to MINE IDs...\n")
    sys.stdout.flush()
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"
    kegg_id_dict = {}
    for kegg_id in kegg_ids:
        try:
            kegg_id_dict[kegg_id] = QuickSearch(con, db, kegg_id)[0]['_id']
        except IndexError:
            sys.stderr.write("Warning: '%s' is not present in the database.\n" % kegg_id)
            sys.stderr.flush()
            continue
    print("Done.")
    return kegg_id_dict

def test_KeggToMineId():
    assert KeggToMineId(['C15667', 'C16519', 'C00130']) == {'C15667':'C023e725c27a385fd9057c8171ca1992a32a3e9a4',
    'C16519':'Cfa3c20a209e12ac4a70f10530cf43bea2b4422fe',
    'C00130':'Cae5be0a0a2bb6b6baef29e4d4c6b3f4c1a67ad19'}
    # C00003 is not present in the database
    assert KeggToMineId(['C00003', 'C16519', 'C00130']) == {'C16519':'Cfa3c20a209e12ac4a70f10530cf43bea2b4422fe',
    'C00130':'Cae5be0a0a2bb6b6baef29e4d4c6b3f4c1a67ad19'}


def ExtractReactionCompIds(rxn):
    """Extracts all compound IDs (reactants and products) from a MINE reaction object."""

    rxn_comp_ids = []

    # Get reaction ID and test if the reaction is valid
    try:
        rxn_id = rxn['_id']
    except KeyError:
        sys.stderr.write("Warning: '%s' does not have a reaction ID.\n" % str(rxn))
        sys.stderr.flush()
        rxn_id = 'UnknownReaction'
    except TypeError:
        sys.stderr.write("Warning: '%s' is not a valid reaction.\n" % str(rxn))
        sys.stderr.flush()
        return rxn_comp_ids

    # Try to get the reactants
    try:
        rxn_p = rxn['Reactants']
        try:
            rxn_comp_ids.extend([x[1] for x in rxn_p])
        except IndexError:
            sys.stderr.write("Warning: The reactant list of '%s' is not valid.\n" % rxn_id)
            sys.stderr.flush()
    except KeyError:
        sys.stderr.write("Warning: '%s' does not list its reactants.\n" % rxn_id)
        sys.stderr.flush()

    # Try to get the products
    try:
        rxn_p = rxn['Products']
        try:
            rxn_comp_ids.extend([x[1] for x in rxn_p])
        except IndexError:
            sys.stderr.write("Warning: The product list of '%s' is not valid.\n" % rxn_id)
            sys.stderr.flush()
    except KeyError:
        sys.stderr.write("Warning: '%s' does not list its products.\n" % rxn_id)
        sys.stderr.flush()

    return rxn_comp_ids

def test_ExtractReactionCompIds(capsys):
    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    rxn1 = {'_id':'R1','Products':[[1,'C1'],[1,'C2']],'Reactants':[[1,'X1'],[1,'X2']]}
    rxn2_id = 'Re598257045ae3ce45dabf450b57708d84e642558'
    rxn2 = con.get_rxns(db, [rxn2_id])[0]
    rxn3 = {'_id':'R3','Reactants':[['XZ']]}

    assert set(ExtractReactionCompIds(rxn1)) == set(['C1', 'C2', 'X1', 'X2'])
    assert set(ExtractReactionCompIds(rxn2)) == set([x[1] for x in rxn2['Products']] + [x[1] for x in rxn2['Reactants']])

    ExtractReactionCompIds(rxn3)
    out, err = capsys.readouterr()
    assert err == "Warning: The reactant list of 'R3' is not valid.\nWarning: 'R3' does not list its products.\n"


def LimitCarbon(con, db, rxn, C_limit=25):
    """Returns True if a compound in the reaction exceeds the carbon atom limit."""
    regex = re.compile('C{1}[0-9]*')
    comps = GetComps(con, db, ExtractReactionCompIds(rxn))
    try:
        rxn_id = rxn['_id']
    except KeyError:
        sys.stderr.write("Warning: Reaction '%s' will be not be considered as it lacks an ID." % str(rxn))
        sys.stderr.flush()
        return False
    for comp in comps:
        if type(comp) != dict:
            # Skip compound if it is not a dictionary (i.e. None for lack of db entry)
            continue
        try:
            formula = comp['Formula']
        except KeyError:
            try:
                comp_id = comp['_id']
            except KeyError:
                comp_id = 'UnknownCompound'
            sys.stderr.write("Warning: Compound '%s' in reaction '%s' lacks a formula and will not be evaluated for C content." % (comp_id, rxn_id))
            sys.stderr.flush()
            continue
        match = re.search(regex, formula)
        if match:
            try:
                C_count = int(match.group(0).split('C')[1])
            except ValueError:
                C_count = 1
        else:
            C_count = 0
        if C_count > C_limit:
            return True

def test_LimitCarbon():
    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    rxn = con.get_rxns(db, ['R180569d0b4cec9c8392f78015bf8d5341ca05c66'])[0]
    assert LimitCarbon(con, db, rxn, 25)
    assert not LimitCarbon(con, db, rxn, 50)

    rxn = con.get_rxns(db, ['R25b1c5f3ec86899ccbd244413c5e53140c626646'])[0]
    assert not LimitCarbon(con, db, rxn)
    assert LimitCarbon(con, db, rxn, 20)


def GetRawNetwork(comp_id_list, step_limit=10, comp_limit=100000, C_limit=25):
    """Download connected reactions and compounds up to the limits."""

    sys.stdout.write("\nDownloading raw network data...\n")
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
        sys.stdout.write("\rStep %s: Compound %s ('%s')..." % (str(steps), str(comps), comp_id))
        sys.stdout.flush()
        try:
            comp = GetComps(con, db, [comp_id])[0]
            comp_dict[comp_id] = comp # Add compound to dict
        except IndexError:
            sys.stderr.write("Warning: '%s' could not be retrieved from the database.\n" % comp_id)
            sys.stderr.flush()

    extended_comp_ids = set()

    # Perform stepwise expansion of downloaded data
    while steps < step_limit:
        steps += 1
        print("")
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
                    try:
                        rxn = GetRxns(con, db, [rxn_id])[0]
                        if rxn != None:
                            if not LimitCarbon(con, db, rxn, C_limit):
                                rxn_dict[rxn_id] = rxn # Add new reaction if it exists and passes the C limit
                            else:
                                continue
                        else:
                            sys.stderr.write("Warning: '%s' appears to have no record in the database.\n" % rxn_id)
                            sys.stderr.flush()
                            continue
                    except IndexError:
                        sys.stderr.write("Warning: '%s' could not be retrieved from the database.\n" % rxn_id)
                        sys.stderr.flush()
                        continue
                rxn_comp_ids = ExtractReactionCompIds(rxn)
                for rxn_comp_id in rxn_comp_ids:
                    # Only download new compounds
                    try:
                        rxn_comp = comp_dict[rxn_comp_id]
                    except KeyError:
                        try:
                            rxn_comp = GetComps(con, db, [rxn_comp_id])[0]
                            if rxn_comp != None:
                                comp_dict[rxn_comp_id] = rxn_comp # Add new compound
                                comps += 1
                                sys.stdout.write("\rStep %s: Compound %s ('%s')..." % (str(steps), str(comps), rxn_comp_id))
                                sys.stdout.flush()
                            else:
                                sys.stderr.write("Warning: '%s' in '%s' appears to have no record in the database.\n" % (rxn_comp_id, rxn_id))
                                sys.stderr.flush()
                            # Break here if comp_limit is reached
                            if comps >= comp_limit:
                                print("Done.")
                                return (comp_dict, rxn_dict)
                        except IndexError:
                            sys.stderr.write("Warning: '%s' in '%s' could not be retrieved from the database.\n" % (rxn_id, rxn_id))
                            sys.stderr.flush()
            extended_comp_ids.add(comp_id)
    print("Done.")
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

    assert GetRawNetwork(compound_ids, step_limit=2, C_limit=500) == (comp_dict, rxn_dict)

    # NAD+ should not be connected via reactions
    nad_plus = 'Xf5dc8599a48d0111a3a5f618296752e1b53c8d30'
    nad_comp = con.get_comps(db, [nad_plus])[0]

    assert GetRawNetwork([nad_plus], C_limit=500) == ({nad_plus : nad_comp}, {})

    # Huge compounds are not allowed to grow
    huge = 'Caf6fc55862387e5fd7cd9635ef9981da7f08a531'
    huge_comp = con.get_comps(db, [huge])[0]

    assert GetRawNetwork([huge]) == ({huge : huge_comp}, {})


def AddCompoundNode(graph, compound, start_comp_ids):
    """Adds a compound node to the graph."""
    N = len(graph.nodes()) + 1
    try:
        mid = compound['_id']
    except:
        sys.stderr.write("Warning: Compound '%s' is malformed and will not be added to the network.\n" % str(compound))
        sys.stderr.flush()
        return graph
    if mid in start_comp_ids:
        start = True
    else:
        start = False
    graph.add_node(N, type='c', mid=mid, start=start)
    return graph

def test_AddCompoundNode():
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    db = "KEGGexp2"
    con = mc.mineDatabaseServices(server_url)

    G1 = nx.DiGraph()
    G2 = nx.DiGraph()
    comp1 = con.get_comps(db, ['Xf5dc8599a48d0111a3a5f618296752e1b53c8d30'])[0]
    comp2 = con.get_comps(db, ['C38a97a9f962a32b984b1702e07a25413299569ab'])[0]

    G2.add_node(1, type='c', mid=comp1['_id'], start=True)
    G2.add_node(2, type='c', mid=comp2['_id'], start=False)

    sids = set(['Xf5dc8599a48d0111a3a5f618296752e1b53c8d30'])

    assert nx.is_isomorphic(AddCompoundNode(AddCompoundNode(G1, comp1, sids), comp2, sids), G2)
    assert G1.node[1]['mid'] == G2.node[1]['mid']
    assert G1.node[2]['mid'] == G2.node[2]['mid']
    assert G1.node[1]['start'] == G2.node[1]['start'] == True
    assert G1.node[2]['start'] == G2.node[2]['start'] == False
    assert G1.nodes(data=True) == G2.nodes(data=True)


def CheckConnection(minetwork, c_node, r_node):
    """Checks that the compound-to-reaction node connection is valid."""

    con_check = False

    c_mid = minetwork.node[c_node]['mid']
    r_mid = minetwork.node[r_node]['mid']
    r_type = minetwork.node[r_node]['type']

    if r_type in {'rf','pr'}:
        try:
            if r_mid in minetwork.graph['mine_data'][c_mid]['Reactant_in']:
                con_check = True
        except KeyError:
            pass

    if r_type in {'pf','rr'}:
        try:
            if r_mid in minetwork.graph['mine_data'][c_mid]['Product_of']:
                con_check = True
        except KeyError:
            pass

    return con_check

def test_CheckConnection(capsys):
    G = nx.DiGraph(
    mine_data = {'C1':{'_id':'C1','Reactant_in':['R1']}, 'C2':{'_id':'C2','Reactant_in':['R1']}, 'R1':{'_id':'R1', 'Reactants':[[1,'C1'],[1,'C2']], 'Products':[[1,'C3']]}, 'C3':{'_id':'C3','Product_of':['R1']}},
    )
    G.add_node(1, type='c', mid='C1')
    G.add_node(2, type='c', mid='C2')
    G.add_node(3, type='rf', mid='R1', c=set([1,2]))
    G.add_node(4, type='pf', mid='R1', c=set([5]))
    G.add_node(5, type='c', mid='C3')
    G.add_node(6, type='rr', mid='R1', c=set([5]))
    G.add_node(7, type='pr', mid='R1', c=set([1,2]))

    assert CheckConnection(G, 1, 3)
    assert CheckConnection(G, 5, 6)

    assert not CheckConnection(G, 2, 4)


def AddQuadReactionNode(graph, rxn):
    """
    Adds a "Quad Reaction Node" (QRN) group of nodes to a graph, and connects
    them to the correct compound nodes.

    The QRN consists of two nodes constituting the intended forward direction
    of the reaction and two nodes constituting the reverse direction. Each pair
    of nodes is connected by an edge in the direction of the reaction. Each node
    represents a group of compounds on one side of the reaction equation.
    """

    # Make sure the reaction is in good shape

    rxn_malformed = False

    try:
        rxn_id = rxn['_id']
    except:
        rxn_malformed = True

    try:
        reactants_f = set([x[1] for x in rxn['Reactants']])
        products_f = set([x[1] for x in rxn['Products']])
        reactants_r = products_f
        products_r = reactants_f
    except:
        rxn_malformed = True

    if rxn_malformed:
        sys.stderr.write("Warning: Reaction '%s' is malformed and will not be added to the network.\n" % str(rxn))
        sys.stderr.flush()
        return graph

    # Find the compound nodes of the reactants and the products
    rf = set([])
    pf = set([])
    rr = set([])
    pr = set([])

    for node in graph.nodes():
        if graph.node[node]['type'] != 'c':
            continue
        if graph.node[node]['mid'] in reactants_f:
            rf.add(node)
            pr.add(node)
        if graph.node[node]['mid'] in products_f:
            pf.add(node)
            rr.add(node)

    # Create the reaction nodes
    N = len(graph.nodes()) + 1

    graph.add_node(N, type='rf', mid=rxn_id, c=rf)
    for c_node in rf:
        if CheckConnection(graph, c_node, N):
            graph.add_edge(c_node, N)

    N += 1

    graph.add_node(N, type='pf', mid=rxn_id, c=pf)
    for c_node in pf:
        if CheckConnection(graph, c_node, N):
            graph.add_edge(N, c_node)

    graph.add_edge(N-1, N) # Forward reaction edge

    N += 1

    graph.add_node(N, type='rr', mid=rxn_id, c=rr)
    for c_node in rr:
        if CheckConnection(graph, c_node, N):
            graph.add_edge(c_node, N)

    N += 1

    graph.add_node(N, type='pr', mid=rxn_id, c=pr)
    for c_node in pr:
        if CheckConnection(graph, c_node, N):
            graph.add_edge(N, c_node)

    graph.add_edge(N-1, N) # Reverse reaction edge

    return graph

def test_AddQuadReactionNode():
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    db = "KEGGexp2"
    con = mc.mineDatabaseServices(server_url)

    rxn = con.get_rxns(db, ['R04759e864c86cfd0eaeb079404d5f18dae6c7227'])[0]

    r_mid = ['Caf6fc55862387e5fd7cd9635ef9981da7f08a531', 'X25a9fafebc1b08a0ae0fec015803771c73485a61']
    p_mid = ['Cefbaa83ea06e7c31820f93c1a5535e1378aba42b', 'Xf729c487f9b991ec6f645c756cf34b9a20b9e8a4']
    r_node_c = set([1,2])
    p_node_c = set([3]) # Xf729c487f9... is here not a compound node in the network

    G1 = nx.DiGraph()
    G1.add_node(1, type='c', mid=con.get_comps(db, [r_mid[0]])[0]['_id'], start=False) # r_mid[0] is connected
    G1.add_node(2, type='c', mid=con.get_comps(db, [r_mid[1]])[0]['_id'], start=True) # r_mid[1] is ATP; not connected
    G1.add_node(3, type='c', mid=con.get_comps(db, [p_mid[0]])[0]['_id'], start=False) # p_mid[0] is connected

    rxn_id = 'R04759e864c86cfd0eaeb079404d5f18dae6c7227'

    G1.add_node(4, type='rf', mid=rxn_id, c=r_node_c) # Forward (intended) direction reactants
    G1.add_node(5, type='pf', mid=rxn_id, c=p_node_c) # Forward (intended) direction products
    G1.add_node(6, type='rr', mid=rxn_id, c=p_node_c) # Reverse direction reactants
    G1.add_node(7, type='pr', mid=rxn_id, c=r_node_c) # Reverse direction products
    G1.add_edge(4, 5) # Directed edge for the forward reaction
    G1.add_edge(6, 7) # Directed edge for the reverse reaction

    # Edges connecting compound and reaction nodes
    G1.add_edge(1, 4)
    #G1.add_edge(2, 4) # ATP should not be connected
    G1.add_edge(5, 3)
    G1.add_edge(3, 6)
    G1.add_edge(7, 1)
    #G1.add_edge(7, 2) # ATP should not be connected

    G2 = nx.DiGraph(mine_data={
    r_mid[0] : con.get_comps(db, [r_mid[0]])[0],
    r_mid[1] : con.get_comps(db, [r_mid[1]])[0],
    p_mid[0] : con.get_comps(db, [p_mid[0]])[0]
    })
    G2.add_node(1, type='c', mid=con.get_comps(db, [r_mid[0]])[0]['_id'], start=False)
    G2.add_node(2, type='c', mid=con.get_comps(db, [r_mid[1]])[0]['_id'], start=True)
    G2.add_node(3, type='c', mid=con.get_comps(db, [p_mid[0]])[0]['_id'], start=False)
    G2 = AddQuadReactionNode(G2, rxn)

    assert nx.is_isomorphic(G1, G2)
    assert G1.nodes(data=True) == G2.nodes(data=True)

    r1 = {'_id':'R1','Reactants':[[1,'C1'],[1,'X1']],'Products':[[1,'C2'],[1,'X2']]}
    r2 = {'_id':'R2','Reactants':[[1,'C2']],'Products':[[1,'C3']]}
    c1 = {'_id':'C1','Reactant_in':['R1']}
    c2 = {'_id':'C2','Product_of':['R1'],'Reactant_in':['R2']}
    c3 = {'_id':'C3','Product_of':['R2']}
    x1 = {'_id':'X1'}
    G3 = nx.DiGraph(mine_data={'R1':r1,'R2':r2,'C1':c1,'C2':c2,'C3':c3})
    G3.add_node(1,mid='C1',type='c')
    G3.add_node(2,mid='C2',type='c')
    G3.add_node(3,mid='C3',type='c')
    G3.add_node(4,mid='X1',type='c')
    G3 = AddQuadReactionNode(G3, r1)
    G3 = AddQuadReactionNode(G3, r2)

    assert len(G3.edges()) == 12
    assert len(G3.nodes()) == 12


def ConstructNetwork(comp_dict, rxn_dict, start_comp_ids=[]):
    """Constructs a directed graph (network) from the compound and reaction
    dictionaries produced by GetRawNetwork."""

    sys.stdout.write("\nConstructing network...\n")
    sys.stdout.flush()

    start_comp_ids = set(start_comp_ids)

    # Initialise directed graph
    minetwork = nx.DiGraph(mine_data={**comp_dict, **rxn_dict})

    # Add all compounds
    for comp_id in sorted(comp_dict.keys()):
        comp = comp_dict[comp_id]
        minetwork = AddCompoundNode(minetwork, comp, start_comp_ids)

    # Add all reactions
    for rxn_id in sorted(rxn_dict.keys()):
        rxn = rxn_dict[rxn_id]
        minetwork = AddQuadReactionNode(minetwork, rxn)

    print("Done.")
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
    rxn_dict['R1e'] = {'_id':'R1e', 'Products':[[1,'C4'],[1,'C5']], 'Reactants':[[1,'C2']]} #9
    rxn_dict['R2f'] = {'_id':'R2f', 'Products':[[1,'C7'],[1,'C8']], 'Reactants':[[1,'C5']]} #13
    rxn_dict['R99'] = {'_id':'R99', 'Products':[[1,'C2'],[1,'C3']], 'Reactants':[[1,'C1']]} #17
    rxn_dict['Rc3'] = {'_id':'Rc3', 'Products':[[1,'C3']], 'Reactants':[[1,'C5']]} #21
    rxn_dict['Rcd'] = {'_id':'Rcd', 'Products':[[1,'C5'],[1,'C6']], 'Reactants':[[1,'C3']]} #25

    G = nx.DiGraph(mine_data = {**comp_dict, **rxn_dict})

    start_comp_ids = set(['C1'])

    for comp_id in sorted(comp_dict.keys()):
        comp = comp_dict[comp_id]
        G = AddCompoundNode(G, comp, start_comp_ids)

    N = 8

    for rxn_id in sorted(rxn_dict.keys()):
        rxn = rxn_dict[rxn_id]
        reactants = set([int(x[1][1]) for x in rxn['Reactants']])
        products = set([int(x[1][1]) for x in rxn['Products']])
        N += 1
        G.add_node(N, type='rf', mid=rxn_id, c=reactants)
        N += 1
        G.add_node(N, type='pf', mid=rxn_id, c=products)
        G.add_edge(N-1, N)
        N += 1
        G.add_node(N, type='rr', mid=rxn_id, c=products)
        N += 1
        G.add_node(N, type='pr', mid=rxn_id, c=reactants)
        G.add_edge(N-1, N)

    # C1 edges
    G.add_edge(1, 17)
    G.add_edge(20, 1)

    # C2 edges
    G.add_edge(18, 2)
    G.add_edge(2, 19)
    G.add_edge(2, 9)
    G.add_edge(12, 2)

    # C3 edges
    G.add_edge(18, 3)
    G.add_edge(3, 19)
    G.add_edge(3, 23)
    G.add_edge(22, 3)
    G.add_edge(3, 25)
    G.add_edge(28, 3)

    # C4 edges
    G.add_edge(10, 4)
    G.add_edge(4, 11)

    # C5 edges
    G.add_edge(10, 5)
    G.add_edge(5, 11)
    G.add_edge(24, 5)
    G.add_edge(5, 21)
    G.add_edge(26, 5)
    G.add_edge(5, 27)
    G.add_edge(5, 13)
    G.add_edge(16, 5)

    # C6 edges
    G.add_edge(26, 6)
    G.add_edge(6, 27)

    # C7 edges
    G.add_edge(14, 7)
    G.add_edge(7, 15)

    # C8 edges
    G.add_edge(14, 8)
    G.add_edge(8, 15)

    assert nx.is_isomorphic(ConstructNetwork(comp_dict,rxn_dict,['C1']), G)

    # Test contents node by node
    t = True
    for node in ConstructNetwork(comp_dict,rxn_dict,['C1']).nodes(data=True):
        if G.node[node[0]] != node[1]:
            t = False
            break
    assert t


# Main code block
def main(infile_name, step_limit, comp_limit, C_limit, outfile_name):
    # Get starting compound MINE IDs
    start_ids = list(set(KeggToMineId(ReadCompounds(infile_name)).values()))
    # Create the network
    minetwork = ConstructNetwork(*GetRawNetwork(start_ids, step_limit, comp_limit, C_limit), start_ids)
    # Save to Pickle
    pickle.dump(minetwork, open(outfile_name, 'wb'))


if __name__ == "__main__":
    # Read arguments from the commandline
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='Read KEGG compound identifiers from text file.')
    parser.add_argument('outfile', help='Write MINE network to Python Pickle file.')
    parser.add_argument('-r', type=int, default=10, help='Maximum number of reaction steps to download.')
    parser.add_argument('-c', type=int, default=100000, help='Maximum number of compounds to download.')
    parser.add_argument('-C', type=int, default=25, help='Maximum number of C atoms per molecule for following a reaction.')
    args = parser.parse_args()
    main(args.infile, args.r, args.c, args.C, args.outfile)
