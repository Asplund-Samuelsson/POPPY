#!/usr/bin/env python3

# Import modules
import networkx as nx
import MineClient3 as mc
import re
import sys
import argparse
import pickle
import time
import queue
import threading

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


def ThreadedQuickSearch(con, db, query_list):
    """Threaded implementation of QuickSearch, taking a list of queries as input."""
    def Worker():
        while True:
            query = work.get()
            if query is None:
                break
            sys.stdout.write("\rHandling quick query '%s'." % str(query))
            sys.stdout.flush()
            output.put(QuickSearch(con, db, query))
            work.task_done()

    work = queue.Queue()
    output = queue.Queue()

    threads = []
    num_workers = 128

    for i in range(num_workers):
        t = threading.Thread(target=Worker)
        t.start()
        threads.append(t)

    for query in query_list:
        work.put(query)

    # Block until all work is done
    work.join()

    # Stop workers
    for i in range(num_workers):
        work.put(None)
    for t in threads:
        t.join()

    # Get the results
    results = []

    while not output.empty():
        results.extend(output.get())

    return results

def test_ThreadedQuickSearch():
    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    assert ThreadedQuickSearch(con, db, ['C00022'])[0]['Names'][0] == 'Pyruvate'
    assert ThreadedQuickSearch(con, db, ['random_query']) == []
    assert len(ThreadedQuickSearch(con, db, ['C00022','C01719','C13842','C00231'])) == 4


def GetComp(con, db, comp_id):
    """Wrapper for MineClient3 get_comps() with reconnect functionality."""
    n = 0
    while True:
        try:
            results = con.get_comps(db, [comp_id])
            break
        except mc.ServerError:
            results = None
        except:
            # Server not responding, try again
            n += 1
            if n % 5 == 0:
                message = "Warning: Server not responding after %s attempts ('%s').\n" % (str(n), comp_id)
                sys.stderr.write(message)
                sys.stderr.flush()
            if n >= 36:
                message = "Warning: Connection attempt limit reached. Results negative.\n"
                sys.stderr.write(message)
                sys.stderr.flush()
                results = None
            if n <= 12:
                time.sleep(10)
            if n > 12:
                time.sleep(30)
    try:
        results = results[0]
    except IndexError or TypeError:
        results = None
    if results == None:
        sys.stderr.write("Warning: '%s' could not be retrieved from the database.\n" % comp_id)
        sys.stderr.flush()
    return results

def test_GetComp():
    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    assert GetComp(con, db, 'Cc93137cc81324a5b2872b0bf1c77866c234d66e1')['Formula'] == 'C7H15O10P'
    assert GetComp(con, db, 'Cc93137cc81324a5b2872b0bf1c77866c234d66e1')['dG_error'] == 1.02079
    assert GetComp(con, db, 'not_a_comp_id') == None


def ThreadedGetComps(con, db, comp_id_list):
    """Threaded implementation of GetComp, taking a list of compound ids as input."""
    def Worker():
        while True:
            comp_id = work.get()
            if comp_id is None:
                break
            sys.stdout.write("\rHandling compound query '%s'." % str(comp_id))
            sys.stdout.flush()
            output.put(GetComp(con, db, comp_id))
            work.task_done()

    work = queue.Queue()
    output = queue.Queue()

    threads = []
    num_workers = 128

    for i in range(num_workers):
        t = threading.Thread(target=Worker)
        t.start()
        threads.append(t)

    for comp_id in comp_id_list:
        work.put(comp_id)

    # Block until all work is done
    work.join()

    # Stop workers
    for i in range(num_workers):
        work.put(None)
    for t in threads:
        t.join()

    # Get the results
    comps = []

    while not output.empty():
        comps.append(output.get())

    return comps

def test_ThreadedGetComps():
    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    comp_ids = ['C1bb250660ea917ddaa2b2777b4773facd6bebb33',
    'C9effc25891ed5be2d4e0804f72e7c78f24e08825',
    'Ce0b888f73c8eabf45289f3fd8e564ff0a92f0014',
    'Cee14c71f197998d923eefb144761a1626a87b738',
    'C6efa5f2bc583af46e2f0c53f112c875abc916d37']

    comps = [con.get_comps(db, [comp_id])[0] for comp_id in comp_ids]

    comps_t = ThreadedGetComps(con, db, comp_ids)

    elements_identical = True

    for e in comps:
        if not e in comps_t:
            elements_identical = False
    for e in comps_t:
        if not e in comps:
            elements_identical = False

    assert elements_identical


def GetRxn(con, db, rxn_id):
    """Wrapper for MineClient3 get_rxns() with reconnect functionality."""
    n = 0
    while True:
        try:
            results = con.get_rxns(db, [rxn_id])
            break
        except mc.ServerError:
            results = None
        except:
            # Server not responding, try again
            n += 1
            if n % 5 == 0:
                message = "Warning: Server not responding after %s attempts ('%s').\n" % (str(n), rxn_id)
                sys.stderr.write(message)
                sys.stderr.flush()
            if n >= 36:
                message = "Warning: Connection attempt limit reached. Results negative.\n"
                sys.stderr.write(message)
                sys.stderr.flush()
                results = None
            if n <= 12:
                time.sleep(10)
            if n > 12:
                time.sleep(30)
    try:
        results = results[0]
    except IndexError or TypeError:
        results = None
    if results == None:
        sys.stderr.write("Warning: '%s' could not be retrieved from the database.\n" % rxn_id)
        sys.stderr.flush()
    return results


def test_GetRxn():

    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    rxn_id = 'Re598257045ae3ce45dabf450b57708d84e642558'
    rxn_op = '1.14.13.e'
    rxn_rlen = 4

    assert type(GetRxn(con, db, rxn_id)) == dict
    assert GetRxn(con, db, rxn_id)['Operators'] == [rxn_op]
    assert len(GetRxn(con, db, rxn_id)['Reactants']) == 4
    assert GetRxn(con, db, 'random_reaction') == None


def ThreadedGetRxns(con, db, rxn_id_list):
    """Threaded implementation of GetRxn, taking a list of reaction ids as input."""
    def Worker():
        while True:
            rxn_id = work.get()
            if rxn_id is None:
                break
            sys.stdout.write("\rHandling reaction query '%s'." % str(rxn_id))
            sys.stdout.flush()
            output.put(GetRxn(con, db, rxn_id))
            work.task_done()

    work = queue.Queue()
    output = queue.Queue()

    threads = []
    num_workers = 128

    for i in range(num_workers):
        t = threading.Thread(target=Worker)
        t.start()
        threads.append(t)

    for rxn_id in rxn_id_list:
        work.put(rxn_id)

    # Block until all work is done
    work.join()

    # Stop workers
    for i in range(num_workers):
        work.put(None)
    for t in threads:
        t.join()

    # Get the results
    rxns = []

    while not output.empty():
        rxns.append(output.get())

    return rxns

def test_ThreadedGetRxns():
    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    rxn_ids = ['R39b3f701d4b0949c38469e31ef675cb7ca1b0fde',
    'R62503c9b0dab64629bea90753f557c451ad5a9b1',
    'Rec4bc0816e3e97672c93e81bee581f6710eac00f',
    'Rc72c92a8ea137cdcc3ada34dc2589553f94faf20',
    'Rc1015bf465307226440d0692919c708e8d38cfb1',
    'R47a4684b398ad812c44c5eae69b34972f8a4b624',
    'R1d52cfafb75c8fc3f5dbdbc681c623a02b4014f7']

    rxns = [con.get_rxns(db, [rxn_id])[0] for rxn_id in rxn_ids]

    rxns_t = ThreadedGetRxns(con, db, rxn_ids)

    elements_identical = True

    for e in rxns:
        if not e in rxns_t:
            elements_identical = False
    for e in rxns_t:
        if not e in rxns:
            elements_identical = False

    assert elements_identical


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
    for kegg_comp in ThreadedGetComps(con, db, [x['_id'] for x in ThreadedQuickSearch(con, db, kegg_ids)]):
        for kegg_id in kegg_comp['DB_links']['KEGG']:
            kegg_id_dict[kegg_id] = kegg_comp['_id']
    for kegg_id in kegg_ids:
        try:
            kegg_comp = kegg_id_dict[kegg_id]

        except KeyError:
            sys.stderr.write("Warning: '%s' is not present in the database.\n" % kegg_id)
            sys.stderr.flush()
            continue
    print("\nDone.\n")
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


def LimitCarbon(comp, C_limit=25):
    """Returns True if the compound exceeds the carbon atom limit, otherwise False."""
    regex = re.compile('C{1}[0-9]*')
    try:
        formula = comp['Formula']
    except KeyError:
        try:
            comp_id = comp['_id']
        except KeyError:
            comp_id = 'UnknownCompound'
        sys.stderr.write("Warning: Compound '%s' lacks a formula and will pass the C limit." % (comp_id, rxn_id))
        sys.stderr.flush()
        return False
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
    else:
        return False

def test_LimitCarbon():
    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    rxn = con.get_rxns(db, ['R180569d0b4cec9c8392f78015bf8d5341ca05c66'])[0]

    test_1_25 = []
    test_1_50 = []
    for comp in [con.get_comps(db, [comp_id])[0] for comp_id in ExtractReactionCompIds(rxn)]:
        test_1_25.append(LimitCarbon(comp, 25))
        test_1_50.append(LimitCarbon(comp, 50))

    assert True in test_1_25
    assert True not in test_1_50

    rxn = con.get_rxns(db, ['R25b1c5f3ec86899ccbd244413c5e53140c626646'])[0]

    test_2_def = []
    test_2_20 = []
    for comp in [con.get_comps(db, [comp_id])[0] for comp_id in ExtractReactionCompIds(rxn)]:
        test_2_def.append(LimitCarbon(comp))
        test_2_20.append(LimitCarbon(comp, 20))

    assert True not in test_2_def
    assert True in test_2_20


def ExtractCompReactionIds(comp):
    """Extracts all reaction IDs from a MINE compound object."""
    rxn_id_list = []
    try:
        rxn_id_list.extend(comp['Reactant_in'])
    except KeyError:
        pass
    try:
        rxn_id_list.extend(comp['Product_of'])
    except KeyError:
        pass
    return rxn_id_list


def test_ExtractCompReactionIds():
    C1 = {'_id':'C1', 'Reactant_in':['R1']}
    C2 = {'_id':'C1', 'Reactant_in':['R2'], 'Product_of':['R3']}
    C3 = {'_id':'C1', 'Product_of':['R3', 'R4']}
    C4 = {'_id':'C4'}
    assert ExtractCompReactionIds(C1) == ['R1']
    assert ExtractCompReactionIds(C2) == ['R2','R3']
    assert ExtractCompReactionIds(C3) == ['R3','R4']
    assert ExtractCompReactionIds(C4) == []


def GetRawNetwork(comp_id_list, step_limit=10, comp_limit=100000, C_limit=25):
    """Download connected reactions and compounds up to the limits."""

    sys.stdout.write("\nDownloading raw network data...\n\n")
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
    for comp in ThreadedGetComps(con, db, comp_id_list):
        if comp == None:
            continue
        try:
            comp_id = comp['_id']
        except KeyError:
            sys.stderr.write("Warning: '%s' is not a valid compound.\n" % str(comp))
            sys.stderr.flush()
            continue
        if not LimitCarbon(comp, C_limit):
            comp_dict[comp_id] = comp # Add compound to dict
            comps += 1
        else:
            sys.stderr.write("Warning: Starting compound '%s' exceeds the C limit and is excluded.\n" % comp_id)
            sys.stderr.flush()

    sys.stdout.write("\nStep %s finished at %s compounds.\n" % (str(steps), str(comps)))
    sys.stdout.flush()

    extended_comp_ids = set()
    rxn_exceeding_C_limit = set()
    comp_exceeding_C_limit = set()
    comp_cache = {}

    # Perform stepwise expansion of downloaded data
    while steps < step_limit:
        # A new step begins
        steps += 1
        print("")

        # Get the unexplored compounds by subtracting explored from all that are stored
        unextended_comp_ids = set(comp_dict.keys()) - extended_comp_ids

        # Go through each unexplored compound and get a list of the reactions that need to be downloaded
        rxn_ids_to_download = set()

        for comp_id in unextended_comp_ids:
            comp = comp_dict[comp_id] # New compounds are always in the dictionary
            # Get a list of the reactions that the compound is involved in
            rxn_id_list = ExtractCompReactionIds(comp)
            # Go through each reaction
            for rxn_id in rxn_id_list:
                # Reactions that are not in the reaction dictionary
                # and do not exceed the C limit will be downloaded and further explored
                if rxn_id not in rxn_dict.keys() and rxn_id not in rxn_exceeding_C_limit:
                    rxn_ids_to_download.add(rxn_id)

        # Download new rxns
        new_rxns = ThreadedGetRxns(con, db, list(rxn_ids_to_download))
        print("")

        # Go through the downloaded reactions and get a list of compounds to download
        comp_ids_to_download = set()

        for rxn in new_rxns:
            if rxn == None: continue
            rxn_comp_ids = ExtractReactionCompIds(rxn)
            for rxn_comp_id in rxn_comp_ids:
                # Compounds that are not in the reaction dictionary
                # and do not exceed the C limit will be downloaded and further explored
                if rxn_comp_id not in comp_dict.keys() and rxn_comp_id not in comp_cache.keys() and rxn_comp_id not in comp_exceeding_C_limit:
                    comp_ids_to_download.add(rxn_comp_id)

        # Download new compounds
        new_comps = ThreadedGetComps(con, db, list(comp_ids_to_download))

        # Expand the comp_cache with the new compounds
        for comp in new_comps:
            if comp == None: continue
            try:
                comp_id = comp['_id']
            except KeyError:
                sys.stderr.write("Warning: Compound '%s' lacks an ID and will be skipped.\n" % str(comp))
                sys.stderr.flush()
                continue
            comp_cache[comp_id] = comp

        # Go through each new reaction and its compounds
        for rxn in new_rxns:
            new_rxn_comp_ids = set()
            if rxn == None: continue
            try:
                rxn_id = rxn['_id']
            except KeyError:
                sys.stderr.write("Warning: Reaction '%s' lacks an ID and will be skipped.\n" % str(rxn))
                sys.stderr.flush()
                continue
            rxn_comp_ids = ExtractReactionCompIds(rxn)
            for rxn_comp_id in rxn_comp_ids:
                if rxn_comp_id not in comp_dict.keys():
                    # The compound has not been added to the compound dict
                    try:
                        rxn_comp = comp_cache[rxn_comp_id]
                        if LimitCarbon(rxn_comp, C_limit):
                            # The compound and the reaction both exceed the C limit
                            comp_exceeding_C_limit.add(rxn_comp_id)
                            rxn_exceeding_C_limit.add(rxn_id)
                            # We don't want to explore this reaction further and thus break
                            break
                        # The compound passed the C limit
                        new_rxn_comp_ids.add(rxn_comp_id)
                    except KeyError:
                        # The compound was never downloaded
                        continue
            # We've made it through the compounds of the reaction
            if rxn_id in rxn_exceeding_C_limit:
                continue
            # The reaction did not exceed the C limit, so let's harvest the new compounds
            rxn_dict[rxn_id] = rxn # The reaction should also be placed in the reaction dictionary
            for new_rxn_comp_id in new_rxn_comp_ids:
                comp_dict[new_rxn_comp_id] = comp_cache[new_rxn_comp_id]
                comps += 1
                # Stop at compound limit here
                if comps >= comp_limit:
                    sys.stdout.write("\nStep %s finished at %s compounds.\n" % (str(steps), str(comps)))
                    sys.stdout.flush()
                    print("\nDone.")
                    return (comp_dict, rxn_dict)
        # All reactions in the current step have been explored
        extended_comp_ids = extended_comp_ids.union(unextended_comp_ids)
        sys.stdout.write("\nStep %s finished at %s compounds.\n" % (str(steps), str(comps)))
        sys.stdout.flush()
    print("\nDone.")
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

    assert GetRawNetwork([huge]) == ({}, {})

    # Using octanol to test the carbon limit
    octanol = 'Cf6baa9f91035ac294770d5e0bfbe039e5ab67261'
    C24_comp = 'C479f661686a597fa18f69c533438aa7bf0e1fd89' # This is connected by 1 step

    net_C25 = GetRawNetwork([octanol], 1)
    net_C20 = GetRawNetwork([octanol], 1, C_limit=20)

    try:
        x = net_C25[0][C24_comp]['_id']
    except KeyError:
        x = 1
    try:
        y = net_C20[0][C24_comp]['_id']
    except KeyError:
        y = 1

    assert x == C24_comp
    assert y == 1


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
    elif 'C' + mid[1:] in start_comp_ids:
        start = True
    elif mid[0] == 'X' and not LimitCarbon(compound, 0):
        # Make sure that inorganic X compounds are treated as starting compounds
        start = True
    else:
        start = False
    graph.add_node(N, type='c', mid=mid, start=start)
    # Also add a mid to node entry in the graph dictionary
    try:
        graph.graph['cmid2node'][mid] = N
    except KeyError:
        graph.graph['cmid2node'] = {mid : N}
    return graph

def test_AddCompoundNode():
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    db = "KEGGexp2"
    con = mc.mineDatabaseServices(server_url)

    G1 = nx.DiGraph()
    G2 = nx.DiGraph()
    comp1 = con.get_comps(db, ['Xf5dc8599a48d0111a3a5f618296752e1b53c8d30'])[0]
    comp2 = con.get_comps(db, ['C38a97a9f962a32b984b1702e07a25413299569ab'])[0]
    comp3 = con.get_comps(db, ['X96ff2c653c25b4f3c6fab12b241ec78bff13a751'])[0] # Phosphate - not listed as start, no carbon
    comp4 = con.get_comps(db, ['C89b394fd02e5e5e60ae1e167780ea7ab3276288e'])[0]

    G2.add_node(1, type='c', mid=comp1['_id'], start=True)
    G2.add_node(2, type='c', mid=comp2['_id'], start=False)
    G2.add_node(3, type='c', mid=comp3['_id'], start=True)
    G2.add_node(4, type='c', mid=comp4['_id'], start=True)

    sids = set(['Cf5dc8599a48d0111a3a5f618296752e1b53c8d30', 'C89b394fd02e5e5e60ae1e167780ea7ab3276288e']) # Note that the X has been replaced with C

    for comp in [comp1, comp2, comp3, comp4]:
        AddCompoundNode(G1, comp, sids)

    assert nx.is_isomorphic(G1, G2)

    assert G1.node[1]['mid'] == G2.node[1]['mid']
    assert G1.node[2]['mid'] == G2.node[2]['mid']
    assert G1.node[3]['mid'] == G2.node[3]['mid']
    assert G1.node[4]['mid'] == G2.node[4]['mid']

    assert G1.node[1]['start'] == G2.node[1]['start'] == True
    assert G1.node[2]['start'] == G2.node[2]['start'] == False
    assert G1.node[3]['start'] == G2.node[3]['start'] == True
    assert G1.node[4]['start'] == G2.node[4]['start'] == True

    assert G1.nodes(data=True) == G2.nodes(data=True)

    assert G1.graph['cmid2node'] == {
        'Xf5dc8599a48d0111a3a5f618296752e1b53c8d30':1,
        'C38a97a9f962a32b984b1702e07a25413299569ab':2,
        'X96ff2c653c25b4f3c6fab12b241ec78bff13a751':3,
        'C89b394fd02e5e5e60ae1e167780ea7ab3276288e':4
        }


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

    for c_mid in reactants_f:
        try:
            node = graph.graph['cmid2node'][c_mid]
            rf.add(node)
            pr.add(node)
        except KeyError:
            # If a reactant is missing, the reaction should not be added
            sys.stderr.write("Warning: Compound '%s' in reaction '%s' is missing. Reaction nodes were not added to the network.\n" % (c_mid, rxn_id))
            sys.stderr.flush()
            return graph
    for c_mid in products_f:
        try:
            node = graph.graph['cmid2node'][c_mid]
            pf.add(node)
            rr.add(node)
        except KeyError:
            # If a product is missing, the reaction should not be added
            sys.stderr.write("Warning: Compound '%s' in reaction '%s' is missing. Reaction nodes were not added to the network.\n" % (c_mid, rxn_id))
            sys.stderr.flush()
            return graph

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

def test_AddQuadReactionNode(capsys):
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    db = "KEGGexp2"
    con = mc.mineDatabaseServices(server_url)

    rxn = con.get_rxns(db, ['R04759e864c86cfd0eaeb079404d5f18dae6c7227'])[0]

    r_mid = ['Caf6fc55862387e5fd7cd9635ef9981da7f08a531', 'X25a9fafebc1b08a0ae0fec015803771c73485a61']
    p_mid = ['Cefbaa83ea06e7c31820f93c1a5535e1378aba42b', 'Xf729c487f9b991ec6f645c756cf34b9a20b9e8a4']
    r_node_c = set([1,2])
    p_node_c = set([3,4])

    G1 = nx.DiGraph()
    G1.add_node(1, type='c', mid=r_mid[0], start=False) # r_mid[0] is connected
    G1.add_node(2, type='c', mid=r_mid[1], start=True) # r_mid[1] is ATP; not connected
    G1.add_node(3, type='c', mid=p_mid[0], start=False) # p_mid[0] is connected
    G1.add_node(4, type='c', mid=p_mid[1], start=False) # r_mid[2] is ADP; not connected

    G1.graph['cmid2node'] = {}
    for node in G1.nodes():
        G1.graph['cmid2node'][G1.node[node]['mid']] = node

    rxn_id = 'R04759e864c86cfd0eaeb079404d5f18dae6c7227'

    G1.add_node(5, type='rf', mid=rxn_id, c=r_node_c) # Forward (intended) direction reactants
    G1.add_node(6, type='pf', mid=rxn_id, c=p_node_c) # Forward (intended) direction products
    G1.add_node(7, type='rr', mid=rxn_id, c=p_node_c) # Reverse direction reactants
    G1.add_node(8, type='pr', mid=rxn_id, c=r_node_c) # Reverse direction products
    G1.add_edge(5, 6) # Directed edge for the forward reaction
    G1.add_edge(7, 8) # Directed edge for the reverse reaction

    # Edges connecting compound and reaction nodes
    G1.add_edge(1, 5)
    #G1.add_edge(2, 5) # ATP should not be connected
    G1.add_edge(6, 3)
    #G1.add_edge(6, 4) # ADP should not be connected
    G1.add_edge(3, 7)
    #G1.add_edge(4, 7) # ADP should not be connected
    G1.add_edge(8, 1)
    #G1.add_edge(8, 2) # ATP should not be connected

    G2 = nx.DiGraph(mine_data={
    r_mid[0] : con.get_comps(db, [r_mid[0]])[0],
    r_mid[1] : con.get_comps(db, [r_mid[1]])[0],
    p_mid[0] : con.get_comps(db, [p_mid[0]])[0],
    p_mid[1] : con.get_comps(db, [p_mid[1]])[0]
    })
    G2.add_node(1, type='c', mid=r_mid[0], start=False)
    G2.add_node(2, type='c', mid=r_mid[1], start=True)
    G2.add_node(3, type='c', mid=p_mid[0], start=False)
    G2.add_node(4, type='c', mid=p_mid[1], start=False)

    G2.graph['cmid2node'] = {}
    for node in G2.nodes():
        G2.graph['cmid2node'][G2.node[node]['mid']] = node

    G2 = AddQuadReactionNode(G2, rxn)

    assert nx.is_isomorphic(G1, G2)
    assert G1.nodes(data=True) == G2.nodes(data=True)

    r1 = {'_id':'R1','Reactants':[[1,'C1'],[1,'X1']],'Products':[[1,'C2'],[1,'X2']]}
    r2 = {'_id':'R2','Reactants':[[1,'C2']],'Products':[[1,'C3']]}
    c1 = {'_id':'C1','Reactant_in':['R1']}
    c2 = {'_id':'C2','Product_of':['R1'],'Reactant_in':['R2']}
    c3 = {'_id':'C3','Product_of':['R2']}
    x1 = {'_id':'X1'}
    x2 = {'_id':'X2'}
    G3 = nx.DiGraph(mine_data={'R1':r1,'R2':r2,'C1':c1,'C2':c2,'C3':c3,'X1':x1,'X2':x2})
    G3.add_node(1,mid='C1',type='c')
    G3.add_node(2,mid='C2',type='c')
    G3.add_node(3,mid='C3',type='c')
    G3.add_node(4,mid='X1',type='c')
    G3.add_node(5,mid='X2',type='c')

    G3.graph['cmid2node'] = {}
    for node in G3.nodes():
        G3.graph['cmid2node'][G3.node[node]['mid']] = node

    G3 = AddQuadReactionNode(G3, r1)
    G3 = AddQuadReactionNode(G3, r2)

    assert len(G3.edges()) == 12
    assert len(G3.nodes()) == 13


    G4 = nx.DiGraph()
    G4.add_node(1, type='c', mid=r_mid[0], start=False)
    G4.add_node(2, type='c', mid=r_mid[1], start=True)
    G4.add_node(3, type='c', mid=p_mid[0], start=False)
    # G4.add_node(4, type='c', mid=p_mid[1], start=False) # Skip this node

    G4.graph['cmid2node'] = {}
    for node in G4.nodes():
        G4.graph['cmid2node'][G4.node[node]['mid']] = node

    G4 = AddQuadReactionNode(G4, rxn)

    out, err = capsys.readouterr()
    assert err == "Warning: Compound 'Xf729c487f9b991ec6f645c756cf34b9a20b9e8a4' in reaction 'R04759e864c86cfd0eaeb079404d5f18dae6c7227' is missing. Reaction nodes were not added to the network.\n"
    assert len(G4.edges()) == 0


def ExpandStartCompIds(comp_dict, start_comp_ids, extra_kegg_ids=[]):
    """
    Expands a set of start compound IDs with those of compounds sharing
    the same KEGG ID.
    """
    start_kegg_ids = set(extra_kegg_ids)
    for start_comp_id in start_comp_ids:
        start_comp = comp_dict[start_comp_id]
        try:
            start_kegg_ids = start_kegg_ids.union(set(start_comp['DB_links']['KEGG']))
        except KeyError:
            pass
    for comp_id in comp_dict.keys():
        comp = comp_dict[comp_id]
        try:
            if len(set(comp['DB_links']['KEGG']).intersection(start_kegg_ids)):
                start_comp_ids.add(comp_id)
        except KeyError:
            pass
    return start_comp_ids

def test_ExpandStartCompIds():
    comp_dict = {
        'S1':{'_id':'S1','DB_links':{'KEGG':['C00001']}},
        'S2':{'_id':'S2','DB_links':{'KEGG':['C00002','C00003']}},
        'S3':{'_id':'S3','DB_links':{}},
        'S4':{'_id':'S4'},
        'C1':{'_id':'C1','DB_links':{'KEGG':['C00001']}},
        'C2':{'_id':'C2','DB_links':{}},
        'C3':{'_id':'C3'},
        'C4':{'_id':'C4','DB_links':{'KEGG':['C00002']}},
        'X1':{'_id':'X1','DB_links':{'KEGG':['C00002','C10284']}},
        'X5':{'_id':'X5','DB_links':{'KEGG':['C00006','C00007']}},
        'X6':{'_id':'X6','DB_links':{'KEGG':['C11111']}}
    }
    assert ExpandStartCompIds(comp_dict, set(['S1','S2','S3','S4']), ['C11111']) == set(['S1','S2','S3','S4','C1','C4','X1','X6'])


def ConstructNetwork(comp_dict, rxn_dict, start_comp_ids=[], extra_kegg_ids=[]):
    """Constructs a directed graph (network) from the compound and reaction
    dictionaries produced by GetRawNetwork."""

    sys.stdout.write("\nConstructing network...\n")
    sys.stdout.flush()

    # ExpandStartCompIds catches "unlisted" compounds with the same KEGG ID
    start_comp_ids = ExpandStartCompIds(comp_dict, set(start_comp_ids), extra_kegg_ids=extra_kegg_ids)

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
    start_kegg_ids = ReadCompounds(infile_name)
    start_ids = list(set(KeggToMineId(start_kegg_ids).values()))
    # Create the network
    minetwork = ConstructNetwork(*GetRawNetwork(start_ids, step_limit, comp_limit, C_limit), start_ids, extra_kegg_ids=start_kegg_ids)
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
