#!/usr/bin/env python3

# Import modules
import networkx as nx
import re
import sys
import argparse
import pickle
import time
import queue
import threading
from requests import get as rget
from rdkit import Chem
from copy import deepcopy
from itertools import repeat, product

# Import scripts
import mineclient3 as mc
from mepmap_helpers import *
from mepmap_origin_helpers import *
from mepmap_KEGG_helpers import *
from progress import Progress

# Define functions
def allow_reaction_listing(kegg_comp, kegg_rxn):
    # Is the compound inorganic or CO2?
    if 'Formula' in kegg_comp.keys():
        if not limit_carbon(kegg_comp, 0) or kegg_comp['_id'] == "C00011":
            return False
    # Is the compound CoA or ACP?
    if kegg_comp['_id'] in {"C00010", "C00229"}:
        return False
    # Is the compound a cofactor or part of RP00003?
    if "RPair" in kegg_rxn.keys():
        for rp in kegg_rxn['RPair'].items():
            if kegg_comp['_id'] in rp[1][0] and rp[1][1] == "cofac":
                return False
            if kegg_comp['_id'] in rp[1][0] and rp[0] == "RP00003":
                return False
    # If the compound passed all the tests, the reaction is free to be listed
    return True

def test_allow_reaction_listing():
    # C1 is not a cofactor
    cpd = {"_id":"C1", "Reactions":['R1','R2'], "Formula":"C10"}
    rxn = {"_id":"R1", "RPair":{"RP1":("C1_C2","main"),"RP2":("C3_C4","cofac")}}
    assert allow_reaction_listing(cpd, rxn)

    # C1 is a cofactor
    rxn = {"_id":"R2", "RPair":{"RP3":("C1_C5","cofac"),"RP4":("C6_C7","main")}}
    assert not allow_reaction_listing(cpd, rxn)

    # C1 is not in an RPair
    rxn = {"_id":"R1", "RPair":{"RP1":("C5_C2","main"),"RP2":("C3_C4","cofac")}}
    assert allow_reaction_listing(cpd, rxn)

    # ATP/ADP reaction pair should not be listed (RP00003)
    cpd = {"_id":"C00002", "Reactions":['R1','R2']}
    rxn = {"_id":"R1", "RPair":{"RP00003":("C00002_C00008","ligase")}}
    assert not allow_reaction_listing(cpd, rxn)

    # ATP might be involved in other reactions
    cpd = format_KEGG_compound(get_KEGG_text("C00002"))
    rxn = format_KEGG_reaction(get_KEGG_text("R00085")) # ATP -> AMP
    assert allow_reaction_listing(cpd, rxn)

    # CoA should not be listed
    cpd = format_KEGG_compound(get_KEGG_text("C00010"))
    rxn = format_KEGG_reaction(get_KEGG_text(cpd['Reactions'][201]))
    assert not allow_reaction_listing(cpd, rxn)

    # ACP should not be listed
    cpd = format_KEGG_compound(get_KEGG_text("C00229"))
    rxn = format_KEGG_reaction(get_KEGG_text(cpd['Reactions'][15]))
    assert not allow_reaction_listing(cpd, rxn)

    # Water is often a cofactor and should not be listed
    cpd = format_KEGG_compound(get_KEGG_text("C00001"))
    rxn = format_KEGG_reaction(get_KEGG_text(cpd['Reactions'][45]))
    assert not allow_reaction_listing(cpd, rxn)

    # Inorganic compounds need to be disallowed
    cpd = format_KEGG_compound(get_KEGG_text("C00009"))
    rxn = format_KEGG_reaction(get_KEGG_text(cpd['Reactions'][167]))
    assert not allow_reaction_listing(cpd, rxn)

    # ...as well as CO2
    cpd = format_KEGG_compound(get_KEGG_text("C00011"))
    rxn = format_KEGG_reaction(get_KEGG_text(cpd['Reactions'][89]))
    assert not allow_reaction_listing(cpd, rxn)

    # Keep single-carbon reduced compounds though
    cpd = format_KEGG_compound(get_KEGG_text("C00132")) # Methanol
    rxn = format_KEGG_reaction(get_KEGG_text(cpd['Reactions'][23]))
    assert allow_reaction_listing(cpd, rxn)


def sort_KEGG_reactions(kegg_comp_dict, kegg_rxn_dict, verbose=False):
    """
    Re-organizes reactions of a KEGG compound into 'Reactant_in' and
    'Product_of' categories based on the information contained in the reactions
    dictionary.
    """
    # Go through all compounds
    for kegg_comp_id in kegg_comp_dict.keys():
        kegg_comp = kegg_comp_dict[kegg_comp_id]
        # Go through its reactions
        if "Reactions" in kegg_comp.keys():
            for rxn_id in kegg_comp["Reactions"]:
                if rxn_id in kegg_rxn_dict.keys():
                    rxn = kegg_rxn_dict[rxn_id]
                else:
                    if verbose:
                        s_err("Warning: '" + kegg_comp_id + \
                        "' lists missing reaction '" + rxn_id + "'.\n")
                    continue
                # Check if a reaction listing is allowed
                if not allow_reaction_listing(kegg_comp, rxn):
                    continue
                # Add to Reactant_in list
                if "Reactants" in rxn.keys():
                    if kegg_comp_id in [x[1] for x in rxn["Reactants"]]:
                        try:
                            kegg_comp_dict[kegg_comp_id]\
                            ['Reactant_in'].append(rxn["_id"])
                        except KeyError:
                            kegg_comp_dict[kegg_comp_id]\
                            ['Reactant_in'] = [rxn["_id"]]
                # Add to Product_of list
                if "Products" in rxn.keys():
                    if kegg_comp_id in [x[1] for x in rxn["Products"]]:
                        try:
                            kegg_comp_dict[kegg_comp_id]\
                            ['Product_of'].append(rxn["_id"])
                        except KeyError:
                            kegg_comp_dict[kegg_comp_id]\
                            ['Product_of'] = [rxn["_id"]]

def test_sort_KEGG_reactions():
    # C1 is cofactor in one reaction, not in another
    # C2 is inorganic
    # C3 is a reactant in one reaction, product in another
    # C4 is a product of C3, and lists a reaction that doesn't exist
    # C5 lists a reaction in which it is not listed
    # C6 does not list reactions
    kegg_comp_dict = {
        "C1":{"_id":"C1","Reactions":["R1","R2"],"Formula":"C10H18O2"},
        "C2":{"_id":"C2","Reactions":["R3","R4"],"Formula":"XeF4"},
        "C3":{"_id":"C3","Reactions":["R5","R6"],"Formula":"C10H12O3"},
        "C4":{"_id":"C4","Reactions":["R5","RX"],"Formula":"C2H5O"},
        "C5":{"_id":"C5","Reactions":["R7"],"Formula":"CH3O"},
        "C6":{"_id":"C6","Formula":"C12"}
    }
    kegg_rxn_dict = {
        "R1":{"_id":"R1","Reactants":[[1,"C1"]],"Products":[[1,"X1"]],
              "RPair":{"RP1":("C1_X1","main")}},
        "R2":{"_id":"R2","Reactants":[[1,"C1"],[1,"C100"]],
              "Products":[[1,"C101"],[2,"C10"]],"RPair":{"RP2":("C100_C101",
              "main"),"RP3":("C1_C10","cofac")}},
        "R3":{"_id":"R3","Reactants":[[1,"C2"]],"Products":[[1,"X2"]],
              "RPair":{"RP4":("C2_X2","main")}},
        "R4":{"_id":"R3","Reactants":[[1,"Z2"]],"Products":[[1,"C2"]],
              "RPair":{"RP5":("Z2_C2","main")}},
        "R5":{"_id":"R5","Reactants":[[1,"C3"],[1,"Z9"]],"Products":[[1,"C4"]],
              "RPair":{"RP6":("C3_C4","main"),"RP7":("Z9_C4","trans")}},
        "R6":{"_id":"R6","Reactants":[[1,"C9"]],"Products":[[1,"C8"],[1,"C3"]],
              "RPair":{"RP8":("C9_C3","main")}},
        "R7":{"_id":"R7","Reactants":[[1,"X4"]],"Products":[[1,"Z4"]],
              "RPair":{"RP9":("X4_Z4","main")}}
    }
    expected_comp_dict = {
        "C1":{"_id":"C1","Reactions":["R1","R2"],"Formula":"C10H18O2",
              "Reactant_in":["R1"]},
        "C2":{"_id":"C2","Reactions":["R3","R4"],"Formula":"XeF4"},
        "C3":{"_id":"C3","Reactions":["R5","R6"],"Formula":"C10H12O3",
              "Reactant_in":["R5"],"Product_of":["R6"]},
        "C4":{"_id":"C4","Reactions":["R5","RX"],"Formula":"C2H5O",
              "Product_of":["R5"]},
        "C5":{"_id":"C5","Reactions":["R7"],"Formula":"CH3O"},
        "C6":{"_id":"C6","Formula":"C12"}
    }
    sort_KEGG_reactions(kegg_comp_dict, kegg_rxn_dict) # Edits dict. directly
    assert kegg_comp_dict == expected_comp_dict

    # How about a real example?
    # Butanol (C06142)
    kegg_comp_ids = [
        "C01412","C00005","C00080",
        "C06142","C00006","C00004",
        "C00003"
    ]
    kegg_rxn_ids = ["R03545","R03544"]
    kegg_comp_dict = dict(zip(kegg_comp_ids,
                              [format_KEGG_compound(get_KEGG_text(x)) \
                              for x in kegg_comp_ids]))
    kegg_rxn_dict = dict(zip(kegg_rxn_ids,
                             [format_KEGG_reaction(get_KEGG_text(x)) \
                             for x in kegg_rxn_ids]))

    expected_comp_dict = deepcopy(kegg_comp_dict)
    expected_comp_dict['C06142']['Product_of'] = ['R03544','R03545']
    expected_comp_dict['C01412']['Reactant_in'] = ['R03544','R03545']

    sort_KEGG_reactions(kegg_comp_dict, kegg_rxn_dict)
    assert kegg_comp_dict == expected_comp_dict



def get_raw_KEGG(kegg_comp_ids=[], kegg_rxn_ids=[], krest="http://rest.kegg.jp",
    n_threads=128, test_limit=0):
    """
    Downloads all KEGG compound (C) and reaction (R) records and formats them
    as MINE database compound or reaction entries. The final output is a tuple
    containing a compound dictionary and a reaction dictionary.

    Alternatively, downloads only a supplied list of compounds and reactions.
    """

    s_out("\nDownloading KEGG data via %s/...\n" % krest)

    # Acquire list of KEGG compound IDs
    if not len(kegg_comp_ids):
        s_out("Downloading KEGG compound list...")
        r = rget("/".join([krest,"list","compound"]))
        if r.status_code == 200:
            for line in r.text.split("\n"):
                if line == "": break # The end
                kegg_comp_id = line.split()[0].split(":")[1]
                kegg_comp_ids.append(kegg_comp_id)
        else:
            msg = "Error: Unable to download KEGG rest compound list.\n"
            sys.exit(msg)
        s_out(" Done.\n")

    # Acquire list of KEGG reaction IDs
    if not len(kegg_rxn_ids):
        s_out("Downloading KEGG reaction list...")
        r = rget("/".join([krest,"list","reaction"]))
        if r.status_code == 200:
            for line in r.text.split("\n"):
                if line == "": break # The end
                kegg_rxn_id = line.split()[0].split(":")[1]
                kegg_rxn_ids.append(kegg_rxn_id)
        else:
            msg = "Error: Unable to download KEGG rest reaction list.\n"
            sys.exit(msg)
        s_out(" Done.\n")

    # Limit download length, for testing only
    if test_limit:
        kegg_comp_ids = kegg_comp_ids[0:test_limit]
        kegg_rxn_ids = kegg_rxn_ids[0:test_limit]

    # Download compounds (threaded)
    kegg_comp_dict = {}
    print("Downloading KEGG compounds...")
    for comp in get_KEGG_comps(kegg_comp_ids):
        if comp == None:
            continue
        try:
            kegg_comp_dict[comp['_id']] = comp
        except KeyError:
            s_err("Warning: '" + str(comp) + \
            "' lacks an ID and will be discarded.\n")
            continue

    print("")

    # Download reactions (threaded)
    kegg_rxn_dict = {}
    print("Downloading KEGG reactions...")
    for rxn in get_KEGG_rxns(kegg_rxn_ids):
        if rxn == None:
            continue
        try:
            kegg_rxn_dict[rxn['_id']] = rxn
        except KeyError:
            s_err("Warning: '" + str(rxn) + \
            "' lacks an ID and will be discarded.\n")
            continue

    print("")

    # Re-organize compound reaction listing, taking cofactor role into account
    s_out("Organizing reaction lists...")
    sort_KEGG_reactions(kegg_comp_dict, kegg_rxn_dict)
    s_out(" Done.\n")

    s_out("KEGG download completed.\n")
    return (kegg_comp_dict, kegg_rxn_dict)

def test_get_raw_KEGG_1():
    # Butanol (C06142)
    kegg_comp_ids = [
        "C01412","C00005","C00080",
        "C06142","C00006","C00004",
        "C00003"
    ]
    kegg_rxn_ids = ["R03545","R03544"]
    kegg_comp_dict = dict(zip(kegg_comp_ids,
                              [format_KEGG_compound(get_KEGG_text(x)) \
                              for x in kegg_comp_ids]))
    kegg_rxn_dict = dict(zip(kegg_rxn_ids,
                             [format_KEGG_reaction(get_KEGG_text(x)) \
                             for x in kegg_rxn_ids]))

    kegg_comp_dict['C06142']['Product_of'] = ['R03544','R03545']
    kegg_comp_dict['C01412']['Reactant_in'] = ['R03544','R03545']

    expected = (kegg_comp_dict, kegg_rxn_dict)

    assert get_raw_KEGG(kegg_comp_ids, kegg_rxn_ids) == expected

def test_get_raw_KEGG_2():
    # Random sample
    random_comp_ids = [
    "C14978","C01268","C09868","C05562","C08104",
    "C15636","C14337","C00988","C08400","C19305",
    "C07495","C09986","C04144","C06578","C00508",
    "C17617","C10048","C16549","C04299","C18093"
    ]

    random_rxn_ids = []

    for comp_id in random_comp_ids:
        try:
            random_rxn_ids.extend(KEGG_rest_dict(get_KEGG_text(comp_id))\
            ['REACTION'])
        except KeyError:
            continue

    random_comp_dict = dict(zip(random_comp_ids,
                                [format_KEGG_compound(get_KEGG_text(x)) \
                                for x in random_comp_ids]))
    random_rxn_dict = dict(zip(random_rxn_ids,
                               [format_KEGG_reaction(get_KEGG_text(x)) \
                               for x in random_rxn_ids]))
    sort_KEGG_reactions(random_comp_dict, random_rxn_dict)

    expected = (random_comp_dict, random_rxn_dict)

    assert get_raw_KEGG(random_comp_ids, random_rxn_ids) == expected

def test_get_raw_KEGG_3():
    # First 20 compounds and reactions
    first_comp_ids = [x.split("\t")[0].split(":")[1] \
    for x in rget("http://rest.kegg.jp/list/compound").text.split("\n")[0:20]]
    first_rxn_ids = [x.split("\t")[0].split(":")[1] \
    for x in rget("http://rest.kegg.jp/list/reaction").text.split("\n")[0:20]]

    first_comp_dict = dict(zip(first_comp_ids,
                               [format_KEGG_compound(get_KEGG_text(x)) \
                               for x in first_comp_ids]))
    first_rxn_dict = dict(zip(first_rxn_ids,
                              [format_KEGG_reaction(get_KEGG_text(x)) \
                              for x in first_rxn_ids]))
    sort_KEGG_reactions(first_comp_dict, first_rxn_dict)

    assert get_raw_KEGG(test_limit=20) == (first_comp_dict, first_rxn_dict)


def quicksearch(con, db, query):
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
                s_err("Warning: Server not responding after " + \
                "%s attempts ('%s').\n" % (str(n), query))
            if n >= 36:
                s_err("Warning: Connection attempt limit reached for '" + \
                query + "'.\n")
                return results
            if n <= 12:
                time.sleep(10)
            if n > 12:
                time.sleep(30)

def test_quicksearch():

    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    assert quicksearch(con, db, 'C00022')[0]['Names'][0] == 'Pyruvate'
    assert quicksearch(con, db, 'random_query') == []


def threaded_quicksearch(con, db, query_list):
    """Threaded implementation of quicksearch."""
    def worker():
        while True:
            query = work.get()
            if query is None:
                work.task_done()
                break
            output.put(quicksearch(con, db, query))
            work.task_done()

    work = queue.Queue()
    output = queue.Queue()

    threads = []
    num_workers = 128

    for i in range(num_workers):
        t = threading.Thread(target=worker)
        t.start()
        threads.append(t)

    for query in query_list:
        work.put(query)

    # Report progress
    M = len(query_list)
    p = Progress(design='pbct', max_val=M)
    while True:
        if not work.qsize():
            n = M - work.qsize()
            p.write(n)
            break
        else:
            n = M - work.qsize()
            p.write(n)
            time.sleep(1)
    print("")

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

def test_threaded_quicksearch():
    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    res = threaded_quicksearch(con, db, ['C00022'])[0]['Names'][0]
    assert res == 'Pyruvate'

    res = threaded_quicksearch(con, db, ['random_query'])
    assert res == []

    res = threaded_quicksearch(con, db, ['C00022','C01719','C13842','C00231'])
    assert len(res) == 4


def getcomp(con, db, comp_id):
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
                s_err("Warning: Server not responding after " + \
                "%s attempts ('%s').\n" % (str(n), comp_id))
            if n >= 36:
                s_err("Warning: Connection attempt limit reached for '" + \
                comp_id + "'.\n")
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
        s_err("Warning: '" + comp_id + \
        "' could not be retrieved from the database.\n")
    return results

def test_getcomp():
    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    res = getcomp(con, db, 'Cc93137cc81324a5b2872b0bf1c77866c234d66e1')
    assert res['Formula'] == 'C7H15O10P'

    res = getcomp(con, db, 'Cc93137cc81324a5b2872b0bf1c77866c234d66e1')
    assert res['dG_error'] == 1.02079

    assert getcomp(con, db, 'not_a_comp_id') == None


def threaded_getcomps(con, db, comp_id_list):
    """Threaded implementation of getcomp."""
    def worker():
        while True:
            comp_id = work.get()
            if comp_id is None:
                work.task_done()
                break
            output.put(getcomp(con, db, comp_id))
            work.task_done()

    work = queue.Queue()
    output = queue.Queue()

    threads = []
    num_workers = 128

    for i in range(num_workers):
        t = threading.Thread(target=worker)
        t.start()
        threads.append(t)

    for comp_id in comp_id_list:
        work.put(comp_id)

    # Report progress
    M = len(comp_id_list)
    p = Progress(design='pbct', max_val=M)
    while True:
        if not work.qsize():
            n = M - work.qsize()
            p.write(n)
            break
        else:
            n = M - work.qsize()
            p.write(n)
            time.sleep(1)
    print("")

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

def test_threaded_getcomps():
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

    comps_t = threaded_getcomps(con, db, comp_ids)

    elements_identical = True

    for e in comps:
        if not e in comps_t:
            elements_identical = False
    for e in comps_t:
        if not e in comps:
            elements_identical = False

    assert elements_identical


def getrxn(con, db, rxn_id):
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
                s_err("Warning: Server not responding after " + \
                "%s attempts ('%s').\n" % (str(n), rxn_id))
            if n >= 36:
                s_err("Warning: Connection attempt limit reached for '" + \
                rxn_id + "'.\n")
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
        s_err("Warning: '" + rxn_id + \
        "' could not be retrieved from the database.\n")
    return results


def test_getrxn():

    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    rxn_id = 'Re598257045ae3ce45dabf450b57708d84e642558'
    rxn_op = '1.14.13.e'
    rxn_rlen = 4

    assert type(getrxn(con, db, rxn_id)) == dict
    assert getrxn(con, db, rxn_id)['Operators'] == [rxn_op]
    assert len(getrxn(con, db, rxn_id)['Reactants']) == 4
    assert getrxn(con, db, 'random_reaction') == None


def threaded_getrxn(con, db, rxn_id_list):
    """Threaded implementation of getrxn."""
    def worker():
        while True:
            rxn_id = work.get()
            if rxn_id is None:
                work.task_done()
                break
            output.put(getrxn(con, db, rxn_id))
            work.task_done()

    work = queue.Queue()
    output = queue.Queue()

    threads = []
    num_workers = 128

    for i in range(num_workers):
        t = threading.Thread(target=worker)
        t.start()
        threads.append(t)

    for rxn_id in rxn_id_list:
        work.put(rxn_id)

    # Report progress
    M = len(rxn_id_list)
    p = Progress(design='pbct', max_val=M)
    while True:
        if not work.qsize():
            n = M - work.qsize()
            p.write(n)
            break
        else:
            n = M - work.qsize()
            p.write(n)
            time.sleep(1)
    print("")

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

def test_threaded_getrxn():
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

    rxns_t = threaded_getrxn(con, db, rxn_ids)

    elements_identical = True

    for e in rxns:
        if not e in rxns_t:
            elements_identical = False
    for e in rxns_t:
        if not e in rxns:
            elements_identical = False

    assert elements_identical


def read_compounds(filename):
    """Read a file with KEGG compound IDs."""
    s_out("\nReading compound ID file...")
    compounds = [line.rstrip() for line in open(filename, 'r')]
    for c in compounds:
        if re.fullmatch("^C[0-9]{5}$", c) == None:
            msg = "".join(["\nWarning: The supplied string '",
                           c,
                           "' is not a valid KEGG compound ID."])
            sys.exit(msg)
    s_out(" Done.")
    return compounds

def test_read_compounds():

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

    assert set(read_compounds(f1.name)) == set(['C10000','C40055','C13482'])
    with pytest.raises(SystemExit): read_compounds(f2.name)

    f1.close()
    f2.close()


def KEGG_to_MINE_id(kegg_ids):
    """Translate KEGG IDs to MINE IDs."""
    s_out("\nTranslating from KEGG IDs to MINE IDs...\n")
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"
    kegg_id_dict = {}
    for kegg_comp in threaded_getcomps(con, db, [x['_id'] \
    for x in threaded_quicksearch(con, db, kegg_ids)]):
        for kegg_id in kegg_comp['DB_links']['KEGG']:
            kegg_id_dict[kegg_id] = kegg_comp['_id']
    for kegg_id in kegg_ids:
        try:
            kegg_comp = kegg_id_dict[kegg_id]
        except KeyError:
            s_err("Warning: '%s' is not present in the database.\n" % kegg_id)
            continue
    print("\nDone.\n")
    return kegg_id_dict

def test_KEGG_to_MINE_id():
    assert KEGG_to_MINE_id(['C15667', 'C16519', 'C00130']) == {
        'C15667':'C023e725c27a385fd9057c8171ca1992a32a3e9a4',
        'C16519':'Cfa3c20a209e12ac4a70f10530cf43bea2b4422fe',
        'C00130':'Cae5be0a0a2bb6b6baef29e4d4c6b3f4c1a67ad19'
    }
    # C00003 is not present in the database
    assert KEGG_to_MINE_id(['C00003', 'C16519', 'C00130']) == {
        'C16519':'Cfa3c20a209e12ac4a70f10530cf43bea2b4422fe',
        'C00130':'Cae5be0a0a2bb6b6baef29e4d4c6b3f4c1a67ad19'
    }


def extract_reaction_comp_ids(rxn):
    """Extract reactant and product IDs from a MINE reaction object."""

    rxn_comp_ids = []

    # Get reaction ID and test if the reaction is valid
    try:
        rxn_id = rxn['_id']
    except KeyError:
        s_err("Warning: '%s' does not have a reaction ID.\n" % str(rxn))
        rxn_id = 'UnknownReaction'
    except TypeError:
        s_err("Warning: '%s' is not a valid reaction.\n" % str(rxn))
        return rxn_comp_ids

    # Try to get the reactants
    try:
        rxn_p = rxn['Reactants']
        try:
            rxn_comp_ids.extend([x[1] for x in rxn_p])
        except IndexError:
            s_err("Warning: The reactant list of '%s' is not valid.\n" % rxn_id)
    except KeyError:
        s_err("Warning: '%s' does not list its reactants.\n" % rxn_id)

    # Try to get the products
    try:
        rxn_p = rxn['Products']
        try:
            rxn_comp_ids.extend([x[1] for x in rxn_p])
        except IndexError:
            s_err("Warning: The product list of '%s' is not valid.\n" % rxn_id)
    except KeyError:
        s_err("Warning: '%s' does not list its products.\n" % rxn_id)

    return rxn_comp_ids

def test_extract_reaction_comp_ids(capsys):
    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    rxn1 = {
        '_id':'R1',
        'Products':[[1,'C1'],[1,'C2']],
        'Reactants':[[1,'X1'],[1,'X2']]
    }
    rxn2_id = 'Re598257045ae3ce45dabf450b57708d84e642558'
    rxn2 = con.get_rxns(db, [rxn2_id])[0]
    rxn3 = {'_id':'R3','Reactants':[['XZ']]}

    exp1 = set(['C1', 'C2', 'X1', 'X2'])
    exp2 = set([x[1] for x in rxn2['Products']] + \
    [x[1] for x in rxn2['Reactants']])

    assert set(extract_reaction_comp_ids(rxn1)) == exp1
    assert set(extract_reaction_comp_ids(rxn2)) == exp2

    extract_reaction_comp_ids(rxn3)
    out, err = capsys.readouterr()
    assert err == "\n".join(["Warning: The reactant list of 'R3' is not valid.",
                             "Warning: 'R3' does not list its products.\n"])


def limit_carbon(comp, C_limit=25):
    """True if the compound exceeds the carbon atom limit, otherwise False."""
    regex = re.compile('[A-Z]{1}[a-z]*[0-9]*')
    try:
        formula = comp['Formula']
    except KeyError:
        try:
            comp_id = comp['_id']
        except KeyError:
            comp_id = 'UnknownCompound'
        s_err("Warning: '%s' has no formula and passes the C-limit." % (comp_id))
        return False
    # Find all elements in the formula
    match = re.findall(regex, formula)
    C_count = 0
    for element in match:
        # Check if the element is carbon
        if re.search('C{1}[^a-z]*$', element):
            # Determine carbon count
            try:
                C_count = int(element.split('C')[1])
            except ValueError:
                C_count = 1
            # Assume that the first C encountered is carbon
            # Later C's might be part of e.g. ACP
            break
    if C_count > C_limit:
        return True
    else:
        return False

def test_limit_carbon():
    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    rxn = con.get_rxns(db, ['R180569d0b4cec9c8392f78015bf8d5341ca05c66'])[0]
    rxn_cids = extract_reaction_comp_ids(rxn)

    test_1_25 = []
    test_1_50 = []
    for comp in [con.get_comps(db, [comp_id])[0] for comp_id in rxn_cids]:
        test_1_25.append(limit_carbon(comp, 25))
        test_1_50.append(limit_carbon(comp, 50))

    assert True in test_1_25
    assert True not in test_1_50

    rxn = con.get_rxns(db, ['R25b1c5f3ec86899ccbd244413c5e53140c626646'])[0]
    rxn_cids = extract_reaction_comp_ids(rxn)

    test_2_def = []
    test_2_20 = []
    for comp in [con.get_comps(db, [comp_id])[0] for comp_id in rxn_cids]:
        test_2_def.append(limit_carbon(comp))
        test_2_20.append(limit_carbon(comp, 20))

    assert True not in test_2_def
    assert True in test_2_20

    assert not limit_carbon({'_id':'X1','Formula':'HCl'}, C_limit=0)
    assert not limit_carbon({'_id':'X2','Formula':'CsI'}, C_limit=0)
    assert limit_carbon({'_id':'X3','Formula':'CH2Cl2'}, C_limit=0)


def extract_comp_reaction_ids(comp):
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


def test_extract_comp_reaction_ids():
    C1 = {'_id':'C1', 'Reactant_in':['R1']}
    C2 = {'_id':'C1', 'Reactant_in':['R2'], 'Product_of':['R3']}
    C3 = {'_id':'C1', 'Product_of':['R3', 'R4']}
    C4 = {'_id':'C4'}
    assert extract_comp_reaction_ids(C1) == ['R1']
    assert extract_comp_reaction_ids(C2) == ['R2','R3']
    assert extract_comp_reaction_ids(C3) == ['R3','R4']
    assert extract_comp_reaction_ids(C4) == []


def get_raw_MINE(comp_id_list, step_limit=10, comp_limit=100000, C_limit=25):
    """Download connected reactions and compounds up to the limits."""

    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    s_out("\nDownloading MINE data via %s/...\n\n" % server_url)

    # Set up output dictionaries
    comp_dict = {}
    rxn_dict = {}

    # Set up counters
    steps = 0
    comps = 0

    # First add the starting compounds
    for comp in threaded_getcomps(con, db, comp_id_list):
        if comp == None:
            continue
        try:
            comp_id = comp['_id']
        except KeyError:
            s_err("Warning: '%s' is not a valid compound.\n" % str(comp))
            continue
        if not limit_carbon(comp, C_limit):
            comp_dict[comp_id] = comp # Add compound to dict
            comps += 1
        else:
            s_err("Warning: Starting compound '" + comp_id + \
            "' exceeds the C limit.\n")

    s_out("\nStep %s finished at %s compounds.\n" % (str(steps), str(comps)))

    extended_comp_ids = set()
    rxn_exceeding_C_limit = set()
    comp_exceeding_C_limit = set()
    comp_cache = {}

    # Perform stepwise expansion of downloaded data
    while steps < step_limit:
        # A new step begins
        steps += 1
        print("")

        # Get the unexplored compounds by subtracting
        # explored from all that are stored
        unextended_comp_ids = set(comp_dict.keys()) - extended_comp_ids

        # Go through each unexplored compound and get a list
        # of the reactions that need to be downloaded
        rxn_ids_to_download = set()

        for comp_id in unextended_comp_ids:
            # New compounds are always in the dictionary
            comp = comp_dict[comp_id]
            # Get a list of the reactions that the compound is involved in
            rxn_id_list = extract_comp_reaction_ids(comp)
            # Go through each reaction
            for rxn_id in rxn_id_list:
                # Reactions that are not in the reaction dictionary
                # and do not exceed the C limit will be downloaded
                # and further explored
                req1 = rxn_id not in rxn_dict.keys()
                req2 = rxn_id not in rxn_exceeding_C_limit
                if req1 and req2:
                    rxn_ids_to_download.add(rxn_id)

        # Download new rxns
        new_rxns = threaded_getrxn(con, db, list(rxn_ids_to_download))
        print("")

        # Go through the downloaded reactions and get a list
        # of compounds to download
        comp_ids_to_download = set()

        for rxn in new_rxns:
            if rxn == None: continue
            rxn_comp_ids = extract_reaction_comp_ids(rxn)
            for rxn_comp_id in rxn_comp_ids:
                # Compounds that are not in the reaction dictionary
                # and do not exceed the C limit will be downloaded
                # and further explored
                req1 = rxn_comp_id not in comp_dict.keys()
                req2 = rxn_comp_id not in comp_cache.keys()
                req3 = rxn_comp_id not in comp_exceeding_C_limit
                if req1 and req2 and req3:
                    comp_ids_to_download.add(rxn_comp_id)

        # Download new compounds
        new_comps = threaded_getcomps(con, db, list(comp_ids_to_download))

        # Expand the comp_cache with the new compounds
        for comp in new_comps:
            if comp == None: continue
            try:
                comp_id = comp['_id']
            except KeyError:
                s_err("Warning: Compound '%s' lacks an ID " % str(comp) + \
                "and will be skipped.\n")
                continue
            comp_cache[comp_id] = comp

        # Go through each new reaction and its compounds
        for rxn in new_rxns:
            new_rxn_comp_ids = set()
            if rxn == None: continue
            try:
                rxn_id = rxn['_id']
            except KeyError:
                s_err("Warning: Reaction '%s' lacks an ID " % str(rxn) + \
                "and will be skipped.\n")
                continue
            rxn_comp_ids = extract_reaction_comp_ids(rxn)
            for rxn_comp_id in rxn_comp_ids:
                if rxn_comp_id not in comp_dict.keys():
                    # The compound has not been added to the compound dict
                    try:
                        rxn_comp = comp_cache[rxn_comp_id]
                        if limit_carbon(rxn_comp, C_limit):
                            # Both compound reaction exceed the C limit
                            comp_exceeding_C_limit.add(rxn_comp_id)
                            rxn_exceeding_C_limit.add(rxn_id)
                            # We don't want to explore this reaction further
                            break
                        # The compound passed the C limit
                        new_rxn_comp_ids.add(rxn_comp_id)
                    except KeyError:
                        # The compound was never downloaded
                        continue
            # We've made it through the compounds of the reaction
            if rxn_id in rxn_exceeding_C_limit:
                continue
            # The reaction did not exceed the C limit,
            # so let's harvest the new compounds
            # The reaction should also be placed in the reaction dictionary
            rxn_dict[rxn_id] = rxn
            for new_rxn_comp_id in new_rxn_comp_ids:
                comp_dict[new_rxn_comp_id] = comp_cache[new_rxn_comp_id]
                comps += 1
                # Stop at compound limit here
                if comps >= comp_limit:
                    s_out("".join([
                        "\nStep ", str(steps), " finished at ",
                        str(comps), " compounds.\n"
                    ]))
                    print("\nDone.")
                    return (comp_dict, rxn_dict)
        # All reactions in the current step have been explored
        extended_comp_ids = extended_comp_ids.union(unextended_comp_ids)
        s_out("".join([
            "\nStep ", str(steps), " finished at ", str(comps), " compounds.\n"
        ]))
    print("\nDone.")
    return (comp_dict, rxn_dict)

def test_get_raw_MINE():

    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    db = "KEGGexp2"
    con = mc.mineDatabaseServices(server_url)

    rxn_id_list = [
        'R04759e864c86cfd0eaeb079404d5f18dae6c7227',
        'Re598257045ae3ce45dabf450b57708d84e642558'
    ]

    c_id_1 = [
        a for b in [[y[1] for y in rxn['Reactants']] + \
        [y[1] for y in rxn['Products']] for rxn in \
        con.get_rxns(db, rxn_id_list)] for a in b
    ]

    for comp in con.get_comps(db, c_id_1):
        try:
            rxn_id_list.extend(comp['Reactant_in'])
        except KeyError:
            pass
        try:
            rxn_id_list.extend(comp['Product_of'])
        except KeyError:
            pass

    c_id_2 = [
        a for b in [[y[1] for y in rxn['Reactants']] + \
        [y[1] for y in rxn['Products']] for rxn in \
        con.get_rxns(db, rxn_id_list)] for a in b
    ]

    rxn_id_list = list(set(rxn_id_list))
    c_id_2 = list(set(c_id_2))

    a = con.get_comps(db, c_id_2[0:50])
    b = con.get_comps(db, c_id_2[50:100])
    c = con.get_comps(db, c_id_2[100:])
    comps = a + b + c

    rxns = con.get_rxns(db, rxn_id_list)
    comp_dict = dict([(comp['_id'], comp) for comp in comps])
    rxn_dict = dict([(rxn['_id'], rxn) for rxn in rxns])

    compound_ids = [
        'Cefbaa83ea06e7c31820f93c1a5535e1378aba42b',
        'C38a97a9f962a32b984b1702e07a25413299569ab'
    ]

    res = get_raw_MINE(compound_ids, step_limit=2, C_limit=500)
    assert res == (comp_dict, rxn_dict)

    # NAD+ should not be connected via reactions
    nad_plus = 'Xf5dc8599a48d0111a3a5f618296752e1b53c8d30'
    nad_comp = con.get_comps(db, [nad_plus])[0]

    assert get_raw_MINE([nad_plus], C_limit=500) == ({nad_plus : nad_comp}, {})

    # Huge compounds are not allowed to grow
    huge = 'Caf6fc55862387e5fd7cd9635ef9981da7f08a531'
    huge_comp = con.get_comps(db, [huge])[0]

    assert get_raw_MINE([huge]) == ({}, {})

    # Using octanol to test the carbon limit
    octanol = 'Cf6baa9f91035ac294770d5e0bfbe039e5ab67261'
    C24_comp = 'C479f661686a597fa18f69c533438aa7bf0e1fd89' # Connected by 1 step

    net_C25 = get_raw_MINE([octanol], 1)
    net_C20 = get_raw_MINE([octanol], 1, C_limit=20)

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


def add_compound_node(graph, compound, start_comp_ids):
    """Adds a compound node to the graph."""
    N = len(graph.nodes()) + 1
    try:
        mid = compound['_id']
    except:
        s_err("Warning: Compound '%s' is malformed.\n" % str(compound))
        return graph
    if mid in start_comp_ids:
        start = True
    elif 'C' + mid[1:] in start_comp_ids:
        start = True
    elif not limit_carbon(compound, 0):
        # Make sure that inorganic compounds are treated as starting compounds
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

def test_add_compound_node():
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    db = "KEGGexp2"
    con = mc.mineDatabaseServices(server_url)

    G1 = nx.DiGraph()
    G2 = nx.DiGraph()
    comp1 = con.get_comps(db, ['Xf5dc8599a48d0111a3a5f618296752e1b53c8d30'])[0]
    comp2 = con.get_comps(db, ['C38a97a9f962a32b984b1702e07a25413299569ab'])[0]
    # Phosphate below - not listed as start, no carbon
    comp3 = con.get_comps(db, ['X96ff2c653c25b4f3c6fab12b241ec78bff13a751'])[0]
    comp4 = con.get_comps(db, ['C89b394fd02e5e5e60ae1e167780ea7ab3276288e'])[0]

    G2.add_node(1, type='c', mid=comp1['_id'], start=True)
    G2.add_node(2, type='c', mid=comp2['_id'], start=False)
    G2.add_node(3, type='c', mid=comp3['_id'], start=True)
    G2.add_node(4, type='c', mid=comp4['_id'], start=True)

    sids = set([
        'Cf5dc8599a48d0111a3a5f618296752e1b53c8d30',
        'C89b394fd02e5e5e60ae1e167780ea7ab3276288e'
    ]) # Note that the X has been replaced with C

    for comp in [comp1, comp2, comp3, comp4]:
        add_compound_node(G1, comp, sids)

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


def check_connection(network, c_node, r_node):
    """Checks that the compound-to-reaction node connection is valid."""

    con_check = False

    c_mid = network.node[c_node]['mid']
    r_mid = network.node[r_node]['mid']
    r_type = network.node[r_node]['type']

    if r_type in {'rf','pr'}:
        try:
            if r_mid in network.graph['mine_data'][c_mid]['Reactant_in']:
                con_check = True
        except KeyError:
            pass

    if r_type in {'pf','rr'}:
        try:
            if r_mid in network.graph['mine_data'][c_mid]['Product_of']:
                con_check = True
        except KeyError:
            pass

    return con_check

def test_check_connection(capsys):
    G = nx.DiGraph(
    mine_data = {
        'C1':{'_id':'C1','Reactant_in':['R1']},
        'C2':{'_id':'C2','Reactant_in':['R1']},
        'R1':{'_id':'R1',
              'Reactants':[[1,'C1'],[1,'C2']],
              'Products':[[1,'C3']]},
        'C3':{'_id':'C3','Product_of':['R1']}
        },
    )
    G.add_node(1, type='c', mid='C1')
    G.add_node(2, type='c', mid='C2')
    G.add_node(3, type='rf', mid='R1', c=set([1,2]))
    G.add_node(4, type='pf', mid='R1', c=set([5]))
    G.add_node(5, type='c', mid='C3')
    G.add_node(6, type='rr', mid='R1', c=set([5]))
    G.add_node(7, type='pr', mid='R1', c=set([1,2]))

    assert check_connection(G, 1, 3)
    assert check_connection(G, 5, 6)

    assert not check_connection(G, 2, 4)


def add_quad_reaction_node(graph, rxn):
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
        s_err("Warning: Reaction '%s' is malformed.\n" % str(rxn))
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
            s_err("Warning: Compound '" + c_mid + "' in reaction '" + rxn_id + \
            "' is missing. Reaction nodes were not added to the network.\n")
            return graph
    for c_mid in products_f:
        try:
            node = graph.graph['cmid2node'][c_mid]
            pf.add(node)
            rr.add(node)
        except KeyError:
            # If a product is missing, the reaction should not be added
            s_err("Warning: Compound '" + c_mid + "' in reaction '" + rxn_id + \
            "' is missing. Reaction nodes were not added to the network.\n")
            return graph

    # Create the reaction nodes
    N = len(graph.nodes()) + 1

    graph.add_node(N, type='rf', mid=rxn_id, c=rf)
    for c_node in rf:
        if check_connection(graph, c_node, N):
            graph.add_edge(c_node, N)

    N += 1

    graph.add_node(N, type='pf', mid=rxn_id, c=pf)
    for c_node in pf:
        if check_connection(graph, c_node, N):
            graph.add_edge(N, c_node)

    graph.add_edge(N-1, N) # Forward reaction edge

    N += 1

    graph.add_node(N, type='rr', mid=rxn_id, c=rr)
    for c_node in rr:
        if check_connection(graph, c_node, N):
            graph.add_edge(c_node, N)

    N += 1

    graph.add_node(N, type='pr', mid=rxn_id, c=pr)
    for c_node in pr:
        if check_connection(graph, c_node, N):
            graph.add_edge(N, c_node)

    graph.add_edge(N-1, N) # Reverse reaction edge

    return graph

def test_add_quad_reaction_node(capsys):
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    db = "KEGGexp2"
    con = mc.mineDatabaseServices(server_url)

    rxn = con.get_rxns(db, ['R04759e864c86cfd0eaeb079404d5f18dae6c7227'])[0]

    r_mid = [
        'Caf6fc55862387e5fd7cd9635ef9981da7f08a531',
        'X25a9fafebc1b08a0ae0fec015803771c73485a61'
    ]
    p_mid = [
        'Cefbaa83ea06e7c31820f93c1a5535e1378aba42b',
        'Xf729c487f9b991ec6f645c756cf34b9a20b9e8a4'
    ]
    r_node_c = set([1,2])
    p_node_c = set([3,4])

    G1 = nx.DiGraph()

    # r_mid[0] is connected
    G1.add_node(1, type='c', mid=r_mid[0], start=False)
    # r_mid[1] is ATP; not connected
    G1.add_node(2, type='c', mid=r_mid[1], start=True)
    # p_mid[0] is connected
    G1.add_node(3, type='c', mid=p_mid[0], start=False)
    # r_mid[2] is ADP; not connected
    G1.add_node(4, type='c', mid=p_mid[1], start=False)

    G1.graph['cmid2node'] = {}
    for node in G1.nodes():
        G1.graph['cmid2node'][G1.node[node]['mid']] = node

    rxn_id = 'R04759e864c86cfd0eaeb079404d5f18dae6c7227'

    # Forward (intended) direction reactants
    G1.add_node(5, type='rf', mid=rxn_id, c=r_node_c)
    # Forward (intended) direction products
    G1.add_node(6, type='pf', mid=rxn_id, c=p_node_c)
    # Reverse direction reactants
    G1.add_node(7, type='rr', mid=rxn_id, c=p_node_c)
    # Reverse direction products
    G1.add_node(8, type='pr', mid=rxn_id, c=r_node_c)
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

    G2 = add_quad_reaction_node(G2, rxn)

    assert nx.is_isomorphic(G1, G2)
    assert G1.nodes(data=True) == G2.nodes(data=True)

    r1 = {
        '_id':'R1',
        'Reactants':[[1,'C1'],[1,'X1']],
        'Products':[[1,'C2'],[1,'X2']]
    }
    r2 = {'_id':'R2','Reactants':[[1,'C2']],'Products':[[1,'C3']]}
    c1 = {'_id':'C1','Reactant_in':['R1']}
    c2 = {'_id':'C2','Product_of':['R1'],'Reactant_in':['R2']}
    c3 = {'_id':'C3','Product_of':['R2']}
    x1 = {'_id':'X1'}
    x2 = {'_id':'X2'}
    G3 = nx.DiGraph(mine_data={
        'R1':r1,'R2':r2,'C1':c1,'C2':c2,'C3':c3,'X1':x1,'X2':x2
        })
    G3.add_node(1,mid='C1',type='c')
    G3.add_node(2,mid='C2',type='c')
    G3.add_node(3,mid='C3',type='c')
    G3.add_node(4,mid='X1',type='c')
    G3.add_node(5,mid='X2',type='c')

    G3.graph['cmid2node'] = {}
    for node in G3.nodes():
        G3.graph['cmid2node'][G3.node[node]['mid']] = node

    G3 = add_quad_reaction_node(G3, r1)
    G3 = add_quad_reaction_node(G3, r2)

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

    G4 = add_quad_reaction_node(G4, rxn)

    out, err = capsys.readouterr()
    assert err == "".join([
        "Warning: Compound ",
        "'Xf729c487f9b991ec6f645c756cf34b9a20b9e8a4' in reaction ",
        "'R04759e864c86cfd0eaeb079404d5f18dae6c7227' is missing. ",
        "Reaction nodes were not added to the network.\n"
        ])
    assert len(G4.edges()) == 0


def expand_start_comp_ids(comp_dict, start_comp_ids, extra_kegg_ids=[]):
    """
    Expands a set of start compound IDs with those of compounds sharing
    the same KEGG ID.
    """
    start_kegg_ids = set(extra_kegg_ids)
    for start_comp_id in start_comp_ids:
        try:
            start_comp = comp_dict[start_comp_id]
        except KeyError:
            # Missing start compound IDs missing is not optimal
            # By-passing them here
            continue
        try:
            kegg_ids = set(start_comp['DB_links']['KEGG'])
            start_kegg_ids = start_kegg_ids.union(kegg_ids)
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

def test_expand_start_comp_ids():
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
    start_comp_ids = set(['S1','S2','S3','S4'])
    res = expand_start_comp_ids(comp_dict, start_comp_ids, ['C11111'])
    assert res == set(['S1','S2','S3','S4','C1','C4','X1','X6'])


def construct_network(comp_dict, rxn_dict, start_comp_ids=[], extra_kegg_ids=[]):
    """
    Constructs a directed graph (network) from the compound and reaction
    dictionaries produced by get_raw_MINE and/or get_raw_KEGG.
    """

    s_out("\nConstructing network...\n")

    # expand_start_comp_ids catches "unlisted" compounds with the same KEGG ID
    start_comp_ids = expand_start_comp_ids(comp_dict, set(start_comp_ids),
        extra_kegg_ids=extra_kegg_ids)

    # Initialise directed graph
    network = nx.DiGraph(mine_data={**comp_dict, **rxn_dict})

    # Add all compounds
    p = Progress(max_val=len(comp_dict), design='pct')
    n_done = 0
    for comp_id in sorted(comp_dict.keys()):
        comp = comp_dict[comp_id]
        network = add_compound_node(network, comp, start_comp_ids)
        progress = p.to_string(n_done)
        s_out("\rAdding compounds... %s" % progress)
        n_done += 1

    # Add all reactions
    print("")
    p = Progress(max_val=len(rxn_dict), design='pct')
    n_done = 0
    for rxn_id in sorted(rxn_dict.keys()):
        rxn = rxn_dict[rxn_id]
        network = add_quad_reaction_node(network, rxn)
        progress = p.to_string(n_done)
        s_out("\rAdding reactions... %s" % progress)
        n_done += 1

    print("\nDone.")
    return network

def test_construct_network(capsys):
    comp_dict = {}
    comp_dict['C1'] = {'_id':'C1', 'Reactant_in':['R99']}
    comp_dict['C2'] = {'_id':'C2', 'Reactant_in':['R1e'], 'Product_of':['R99']}
    comp_dict['C3'] = {'_id':'C3', 'Reactant_in':['Rcd'],
                       'Product_of':['R99','Rc3']}
    comp_dict['C4'] = {'_id':'C4', 'Product_of':['R1e']}
    comp_dict['C5'] = {'_id':'C5', 'Reactant_in':['Rc3','R2f'],
                       'Product_of':['Rcd','R1e']}
    comp_dict['C6'] = {'_id':'C6', 'Product_of':['Rcd']}
    # Seeding with non-expanded reaction R7f:
    comp_dict['C7'] = {'_id':'C7', 'Product_of':['R2f','R7f']}
    # Seeding with non-expanded reaction Rb7:
    comp_dict['C8'] = {'_id':'C8', 'Reactant_in':['Rb7'], 'Product_of':['R2f']}

    rxn_dict = {}
    rxn_dict['R1e'] = {
        '_id':'R1e','Products':[[1,'C4'],[1,'C5']], 'Reactants':[[1,'C2']]
    } #9
    rxn_dict['R2f'] = {
        '_id':'R2f', 'Products':[[1,'C7'],[1,'C8']], 'Reactants':[[1,'C5']]
    } #13
    rxn_dict['R99'] = {
        '_id':'R99', 'Products':[[1,'C2'],[1,'C3']], 'Reactants':[[1,'C1']]
    } #17
    rxn_dict['Rc3'] = {
        '_id':'Rc3', 'Products':[[1,'C3']], 'Reactants':[[1,'C5']]
    } #21
    rxn_dict['Rcd'] = {
        '_id':'Rcd', 'Products':[[1,'C5'],[1,'C6']], 'Reactants':[[1,'C3']]
    } #25

    G = nx.DiGraph(mine_data = {**comp_dict, **rxn_dict})

    start_comp_ids = set(['C1'])

    for comp_id in sorted(comp_dict.keys()):
        comp = comp_dict[comp_id]
        G = add_compound_node(G, comp, start_comp_ids)

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

    assert nx.is_isomorphic(construct_network(comp_dict,rxn_dict,['C1']), G)

    # Test contents node by node
    t = True
    for node in construct_network(comp_dict,rxn_dict,['C1']).nodes(data=True):
        if G.node[node[0]] != node[1]:
            t = False
            break
    assert t


def is_connected_MINE_comp(comp_id, network):
    """Determines if the MINE compound is connected."""
    con_keys = set(['Reactant_in','Product_of'])
    if con_keys.intersection(network.graph['mine_data'][comp_id].keys()):
        return True
    else:
        return False

def test_is_connected_MINE_comp():
    G = nx.DiGraph()
    G.graph['mine_data'] = {
        'C1':{'Reactant_in':[]},
        'X2':{},
        'C3':{'Reactant_in':[],'Product_of':[]},
        'C4':{'Product_of':[]},
        'C5':{}
    }
    res = [is_connected_MINE_comp(mid,G) for mid in ['C1','X2','C3','C4','C5']]
    assert res == [True,False,True,True,False]


def KEGG_MINE_integration(network):
    """
    Integrates KEGG and MINE sub-networks

    Transfers incoming and outgoing edges of KEGG nodes to their matching MINE
    nodes. Integrated KEGG nodes are then removed.

    Does not alter the 'mine_data' dictionary of reactions and compounds.
    """

    s_out("\nPerforming KEGG/MINE integration...\n")

    # Identify KEGG compound nodes
    s_out("Identifying KEGG compound nodes...")
    kegg_comp_nodes = []
    for node in network.nodes():
        if network.node[node]['type'] == 'c':
            if re.match('^C{1}[0-9]{5}$', network.node[node]['mid']):
                kegg_comp_nodes.append(node)
    s_out(" Done.\n")

    # Set up dictionary that lists MINE IDs for KEGG IDs
    s_out("Setting up KEGG to MINE ID dictionary...")
    kegg_to_mine = {}
    c_node_count = 0
    mc_node_count = 0
    for node in network.nodes():
        if network.node[node]['type'] == 'c':
            c_node_count += 1
            mine_id = network.node[node]['mid']
            mid_match = re.match('^[CX]{1}[0-9,a-f]{40}$', mine_id)
            if mid_match:
                mc_node_count += 1
                try:
                    kegg_ids = set(network.graph['mine_data']\
                    [mine_id]['DB_links']['KEGG'])
                except KeyError:
                    continue
                for kegg_id in kegg_ids:
                    try:
                        kegg_to_mine[kegg_id].add(mine_id)
                    except KeyError:
                        kegg_to_mine[kegg_id] = set([mine_id])
    s_out(" Done.\n")

    # Set up dictionary that lists disconnected nodes for KEGG nodes
    s_out("Setting up disconnected KEGG node dictionary...")
    kegg_node_dis_nodes = {}
    for node in network.nodes():
        rxn_id = network.node[node]['mid']
        kegg_match = re.match('^R[0-9]{5}$', rxn_id)
        if kegg_match:
            if network.node[node]['type'] in {'rf','rr'}:
                for c_node in network.node[node]['c']:
                    if not network.has_edge(c_node, node):
                        try:
                            kegg_node_dis_nodes[c_node].add(node)
                        except KeyError:
                            kegg_node_dis_nodes[c_node] = set([node])
            if network.node[node]['type'] in {'pf','pr'}:
                for c_node in network.node[node]['c']:
                    if not network.has_edge(node, c_node):
                        try:
                            kegg_node_dis_nodes[c_node].add(node)
                        except KeyError:
                            kegg_node_dis_nodes[c_node] = set([node])

    s_out(" Done.\n")

    # Go through all KEGG compound nodes and perform transplantation
    kc_node_count = c_node_count - mc_node_count
    n = 0
    for kc_node in kegg_comp_nodes:
        n += 1
        progress = float(100 * n / kc_node_count)
        s_out("\rTransferring edges... %0.1f%%" % progress)
        kegg_id = network.node[kc_node]['mid']
        try:
            mine_ids = kegg_to_mine[kegg_id]
        except KeyError:
            # No mine_ids found associated with this KEGG node
            continue
        # Ensure that there is no more than one connected MINE compound
        # node and one disconnected MINE compound node
        n_con = 0
        n_dis = 0
        for mine_id in mine_ids:
            if is_connected_MINE_comp(mine_id, network):
                n_con += 1
            else:
                n_dis += 1
        # Get connected reactant nodes (downstream)
        con_r_nodes = network.successors(kc_node)
        # Get connected product nodes (upstream)
        con_p_nodes = network.predecessors(kc_node)
        # Get disconnected reactant and product nodes
        try:
            dis_nodes = kegg_node_dis_nodes[kc_node]
        except KeyError:
            dis_nodes = set()
        # Go through the corresponding MINE nodes
        for mine_id in mine_ids:
            mine_node = network.graph['cmid2node'][mine_id]
            # Transfer connections to the MINE twin
            if n_con <= 1:
                if is_connected_MINE_comp(mine_id, network):
                    # Incoming edges, from product nodes
                    network.add_edges_from(zip(
                        con_p_nodes, repeat(mine_node,len(con_p_nodes))
                    ))
                    network.remove_edges_from(zip(
                        con_p_nodes, repeat(kc_node, len(con_p_nodes))
                    ))
                    # Outgoing edges, to reactant nodes
                    network.add_edges_from(zip(
                        repeat(mine_node,len(con_r_nodes)), con_r_nodes
                    ))
                    network.remove_edges_from(zip(
                        repeat(kc_node, len(con_r_nodes)), con_r_nodes
                    ))
                    # Update KEGG reactant nodes with the MINE replacement
                    for con_r_node in con_r_nodes:
                        network.node[con_r_node]['c'].add(mine_node)
                        network.node[con_r_node]['c'].discard(kc_node)
                    # Update KEGG product nodes with the MINE replacement
                    for con_p_node in con_p_nodes:
                        network.node[con_p_node]['c'].add(mine_node)
                        network.node[con_p_node]['c'].discard(kc_node)
            else:
                s_err("Warning: KEGG ID '%s' has more than one " % kegg_id + \
                "connected MINE twin.\n")
            if n_dis <= 1:
                if not is_connected_MINE_comp(mine_id, network):
                    # Product and reactant nodes will be updated,
                    # but no edges will be transferred
                    for dis_node in dis_nodes:
                        network.node[dis_node]['c'].add(mine_node)
                        network.node[dis_node]['c'].discard(kc_node)
            else:
                s_err("Warning: KEGG ID '%s' has more than one " % kegg_id + \
                "disconnected MINE twin.\n")
            # A KEGG node with connections,
            # without a twin with connections, is kept as it is
            pass
            # A KEGG node that is disconnected,
            # without a twin that is disconnected, is kept as it is
            pass

    s_out("\nIntegration completed.\n")

def test_KEGG_MINE_integration():
    G = nx.DiGraph()

    # Add KEGG compound nodes
    G.add_node(1,type='c',mid='C00001',start=True) # 'X490c4e...'     A
    G.add_node(2,type='c',mid='C00002',start=False) # 'C683de2...'    B
    G.add_node(3,type='c',mid='C00003',start=False) # 'C069ca5...'    C
    G.add_node(4,type='c',mid='C00004',start=False) # 'C123097...'    D

    # Add MINE compound nodes
    G.add_node(5,type='c',mid='X490c4e9c5d9c3b903bab41ff596eca62ed06130d',start=True) #  A
    G.add_node(6,type='c',mid='C683de2716dd472f4da0a144683d31a10e48a45fc',start=False) # B
    G.add_node(7,type='c',mid='C069ca544492566919b8c9d20984e55b37a9f79a8',start=False) # C
    G.add_node(8,type='c',mid='C123097ef07e00abcd707e873bbd09783da730a38',start=False) # D

    # Add KEGG reaction nodes (A<->C and C<->B)
    G.add_node(9,type='rf',mid='R00001',c=set([1]))
    G.add_node(10,type='pf',mid='R00001',c=set([3]))
    G.add_node(11,type='rr',mid='R00001',c=set([3]))
    G.add_node(12,type='pr',mid='R00001',c=set([1]))
    G.add_path([1,9,10,3])
    G.add_path([3,11,12,1])

    G.add_node(13,type='rf',mid='R00002',c=set([3]))
    G.add_node(14,type='pf',mid='R00002',c=set([2]))
    G.add_node(15,type='rr',mid='R00002',c=set([2]))
    G.add_node(16,type='pr',mid='R00002',c=set([3]))
    G.add_path([3,13,14,2])
    G.add_path([2,15,16,3])

    # Add MINE reaction nodes (Disconnected A<->C and B<->D)
    G.add_node(17,type='rf',mid='Rf2279c67b1b433641502020c3ddd46b911827b88',c=set([5]))
    G.add_node(18,type='pf',mid='Rf2279c67b1b433641502020c3ddd46b911827b88',c=set([7]))
    G.add_node(19,type='rr',mid='Rf2279c67b1b433641502020c3ddd46b911827b88',c=set([7]))
    G.add_node(20,type='pr',mid='Rf2279c67b1b433641502020c3ddd46b911827b88',c=set([5]))
    G.add_path([17,18,7])
    G.add_path([7,19,20])

    G.add_node(21,type='rf',mid='Re9283748451e3dc8254bcd45342926db929b2176',c=set([6]))
    G.add_node(22,type='pf',mid='Re9283748451e3dc8254bcd45342926db929b2176',c=set([8]))
    G.add_node(23,type='rr',mid='Re9283748451e3dc8254bcd45342926db929b2176',c=set([8]))
    G.add_node(24,type='pr',mid='Re9283748451e3dc8254bcd45342926db929b2176',c=set([6]))
    G.add_path([6,21,22,8])
    G.add_path([8,23,24,6])

    # Add mine_data dictionary to network
    G.graph['mine_data'] = {
        'C00001':{"DB_links":{'KEGG':['C00001']},'Reactant_in':['R00001']},
        'X490c4e9c5d9c3b903bab41ff596eca62ed06130d':{"DB_links":{'KEGG':['C00001']}},
        'C00002':{"DB_links":{'KEGG':['C00002']},'Product_of':['R00002']},
        'C683de2716dd472f4da0a144683d31a10e48a45fc':{
            "DB_links":{'KEGG':['C00002']},
            'Reactant_in':['Re9283748451e3dc8254bcd45342926db929b2176']
        },
        'C00003':{
            "DB_links":{'KEGG':['C00003']},
            'Reactant_in':['R00002'],'Product_of':['R00001']
        },
        'C069ca544492566919b8c9d20984e55b37a9f79a8':{
            "DB_links":{'KEGG':['C00003']},
            'Product_of':['Rf2279c67b1b433641502020c3ddd46b911827b88']
        },
        'C00004':{"DB_links":{'KEGG':['C00004']}},
        'C123097ef07e00abcd707e873bbd09783da730a38':{
            "DB_links":{'KEGG':['C00004']},
            'Product_of':['Re9283748451e3dc8254bcd45342926db929b2176']
        },
        'R00001':{},
        'Rf2279c67b1b433641502020c3ddd46b911827b88':{},
        'R00002':{},
        'Re9283748451e3dc8254bcd45342926db929b2176':{}
    }

    # Add the cmid2node dictionary to the network
    G.graph['cmid2node'] = dict(zip(
        [G.node[n]['mid'] for n in range(1,9)],range(1,9)
    ))

    # Copy and integrate
    H = G.copy()
    KEGG_MINE_integration(H)

    # A new path should have been introduced
    assert not nx.has_path(G, 6, 7)
    assert nx.has_path(H, 6, 7)

    # Node 5 should stay disconnected
    assert len(G.predecessors(5)) == len(G.successors(5)) == \
    len(H.predecessors(5)) == len(H.successors(5)) == 0

    # Check the expected connection status of the KEGG nodes
    assert (H.predecessors(1), H.successors(1)) == ([12],[9])
    assert (H.predecessors(2), H.successors(2)) == ([],[])
    assert (H.predecessors(3), H.successors(3)) == ([],[])
    assert (H.predecessors(4), H.successors(4)) == ([],[])

    # Every reaction node should have the correct c node set
    c_sets = [set([x]) for x in [1,7,7,1,7,6,6,7,5,7,7,5,6,8,8,6]]
    res = [H.node[en[0]+9]['c'] == en[1] for en in enumerate(c_sets)]
    assert sum(res) == 16

    # All connections must be transferred correctly
    expected_edges = set([
        (9,10), (11,12), (13,14), (15,16), (17,18), (19,20), (21,22), (23,24),
        (1,9), (12,1),
        (10,7), (7,11),
        (7,13), (14,6), (6,15), (16,7),
        (18,7), (7,19),
        (6,21), (22,8), (8,23), (24,6)
    ])
    assert set(H.edges()) == expected_edges

    # The function should not do anything with non-KEGG reaction nodes
    # or the raw reaction/compound dictionary
    assert G.graph['mine_data'] == H.graph['mine_data']
    assert G.graph['cmid2node'] == H.graph['cmid2node']
    for node in range(17,25):
        assert G.node[node] == H.node[node]

    # Test multiple MINE node situations

    # Compound dict
    comp_dict = {
        'C00001':{
            '_id':'C00001',
            'Reactant_in':['R00001'],'DB_links':{'KEGG':['C00001']}
        },
        'C00002':{
            '_id':'C00002',
            'Product_of':['R00001'],'DB_links':{'KEGG':['C00002']}
        },
        'C00003':{
            '_id':'C00003',
            'Product_of':['R00002'],'DB_links':{'KEGG':['C00003']}
        },
        'C123097ef07e00abcd707e873bbd09783da730a38':{
            '_id':'C123097ef07e00abcd707e873bbd09783da730a38',
            'Reactant_in':['R12f097ef07e00abcd707e873bbd09783da730a38'],
            'DB_links':{'KEGG':['C00001']}
        },
        'X123097ef07e00abcd707e873bbd09783da730a38':{
            '_id':'X123097ef07e00abcd707e873bbd09783da730a38',
            'DB_links':{'KEGG':['C00001']}
        },
        'C31095054707709e8798fbbd89707d0987c8d897c':{
            '_id':'C31095054707709e8798fbbd89707d0987c8d897c',
            'Product_of':[
                'Rfeb4b35607e00abcd707e873bbd09783da730a38',
                'R12f097ef07e00abcd707e873bbd09783da730a38'
                ]
        }
    }

    # Reaction dict
    rxn_dict = {
        'R00001':{
            '_id':'R00001','Reactants':[[1,'C00001']],'Products':[[1,'C00002']]
        },
        'R00002':{
            '_id':'R00002','Reactants':[[1,'C00001']],'Products':[[1,'C00003']]
        },
        'R12f097ef07e00abcd707e873bbd09783da730a38':{
            '_id':'R12f097ef07e00abcd707e873bbd09783da730a38',
            'Reactants':[[1,'C123097ef07e00abcd707e873bbd09783da730a38']],
            'Products':[[1,'C31095054707709e8798fbbd89707d0987c8d897c']]
        },
        'Rfeb4b35607e00abcd707e873bbd09783da730a38':{
            '_id':'Rfeb4b35607e00abcd707e873bbd09783da730a38',
            'Reactants':[[1,'X123097ef07e00abcd707e873bbd09783da730a38']],
            'Products':[[1,'C31095054707709e8798fbbd89707d0987c8d897c']]
        }
    }

    # Construct network
    Y = construct_network(comp_dict, rxn_dict)
    X = Y.copy()

    # Integrate X but not Y
    KEGG_MINE_integration(X)

    # C00001 is connected in R00001 and disconnected in R00002.
    # R00001 should be connected to C123097ef07e00abcd707e873bbd09783da730a38
    # R00002 should list X123097ef07e00abcd707e873bbd09783da730a38 as a reactant,
    # but remain disconnected
    for node in X.nodes():
        if X.node[node]['mid'] == 'R00001' and X.node[node]['type'] == 'rf':
            assert X.node[node]['c'] == \
            set([X.graph['cmid2node']['C123097ef07e00abcd707e873bbd09783da730a38']])
            assert X.predecessors(node) == \
            [X.graph['cmid2node']['C123097ef07e00abcd707e873bbd09783da730a38']]
        if X.node[node]['mid'] == 'R00002' and X.node[node]['type'] == 'rf':
            assert X.node[node]['c'] == \
            set([X.graph['cmid2node']['X123097ef07e00abcd707e873bbd09783da730a38']])
            assert X.predecessors(node) == []

    # Test functinality when multiple KEGG IDs point to the same compound
    G = nx.DiGraph()

    # Add compound nodes
    G.add_node(1,type='c',mid='C10000')
    G.add_node(2,type='c',mid='C20000')
    G.add_node(3,type='c',mid='C123097ef07e00abcd707e873bbd09783da730a38')

    # Add reaction nodes
    G.add_node(10,type='pf',mid='R10000',c=set([1]))
    G.add_node(20,type='rr',mid='R20000',c=set([2]))
    G.add_node(30,type='rf',mid='R123097ef07e00abcd707e873bbd09783da730a38',c=set([3]))

    # Add edges
    G.add_edge(10,1)
    G.add_edge(2,20)
    G.add_edge(3,30)

    # Add the cmid2node and mine_data dictionaries
    G.graph['cmid2node'] = {
        'C10000':1,'C20000':2,'C123097ef07e00abcd707e873bbd09783da730a38':3
    }
    G.graph['mine_data'] = {
        'C10000':{'Product_of':['R10000'],'DB_links':{'KEGG':['C10000']}},
        'C20000':{'Reactant_in':['R20000'],'DB_links':{'KEGG':['C20000']}},
        'C123097ef07e00abcd707e873bbd09783da730a38':{
            'Reactant_in':['R123097ef07e00abcd707e873bbd09783da730a38'],
            'DB_links':{'KEGG':['C10000','C20000']}
        }
    }

    # Integrate and make sure that connections have been transferred
    KEGG_MINE_integration(G)

    assert G.node[10]['c'] == G.node[20]['c'] == G.node[30]['c'] == set([3])
    assert set(G.edges()) == set([(10,3),(3,20),(3,30)])


def expand_valid_compound_set(network, proc_num=1, valid_reactant_nodes=set(),
    comp_node_set=set(), force_parallel=False):
    """Expands the valid compound set with products of valid reactant nodes."""
    if comp_node_set == set():
        comp_node_set = find_start_comp_nodes(network)
    else:
        # Define the worker
        def worker(work):
            compounds = set()
            for r_node in work:
                p_node = network.successors(r_node)[0]
                products = network.node[p_node]['c']
                compounds = compounds.union(products)
            output.extend(compounds)

        # Only go parallel if there are 5k or more items
        if len(valid_reactant_nodes) < 5000 and not force_parallel:
            output = []
            worker(valid_reactant_nodes)
            comp_node_set = comp_node_set.union(set(output))

        else:
            with mp.Manager() as manager:
                # Initialize output list in manager
                output = manager.list()

                # Initialize processes
                procs = []
                for work in chunks(list(valid_reactant_nodes), proc_num):
                    p = mp.Process(target=worker, args=(work,))
                    procs.append(p)
                    p.start()

                # Stop workers
                for p in procs:
                    p.join()

                # Get results
                comp_node_set = comp_node_set.union(set(output))

    return comp_node_set

def test_expand_valid_compound_set():
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

    assert expand_valid_compound_set(G) == set([1])
    assert expand_valid_compound_set(G, force_parallel=True) == set([1])

    assert expand_valid_compound_set(G, 4, set([2,9]), set([1,4])) == set([1,4])
    assert expand_valid_compound_set(G, 4, set([2,9]), set([1,4]), \
    force_parallel=True) == set([1,4])

    assert expand_valid_compound_set(G, 4, find_valid_reactant_nodes(G, 4, \
    set([8])), set([8])) == set([1,5,8])
    assert expand_valid_compound_set(G, 4, find_valid_reactant_nodes(G, 4, \
    set([8]), force_parallel=True), set([8]), force_parallel=True) == set([1,5,8])

    assert expand_valid_compound_set(G, 4, set([2,6,11]), \
    expand_valid_compound_set(G, 4, set([11]), set([8]))) == set([1,4,5,8])
    assert expand_valid_compound_set(G, 4, set([2,6,11]), \
    expand_valid_compound_set(G, 4, set([11]), set([8]), force_parallel=True), \
    force_parallel=True) == set([1,4,5,8])


def distance_to_origin(network, proc_num=1, N=-1):
    """
    Calculates the shortest distance (number of reactions) from the starting
    compounds (origin) to every node up to distance N. Set N to -1 to
    exhaustively calculate the minimum distance to every node that is reachable.

    Returns two sets in a tuple: Valid compound nodes and valid reactant nodes.
    """

    s_out("\nCalculating minimum distance of nodes to origin...\n\n")

    time_start = time.time()

    # Number of nodes for formatting
    L = len(network.nodes())
    l = len(str(L))

    # Set up counters
    n = 0
    c = 0
    rf = 0
    pf = 0
    rr = 0
    pr = 0

    # Start with no valid reactant or compound nodes
    valid_reactant_nodes = set([])
    valid_compound_nodes = set([])

    # The "previous" lists are also empty
    prev_vrn = list(valid_reactant_nodes)
    prev_vcn = list(valid_compound_nodes)

    # Start with no new valid reactant nodes
    new_vrn = set([])

    while True:

        # Valid product nodes will be connected at the start of an expansion
        # cycle. They are however in principle identified in the previous cycle
        # via valid reactant nodes
        for r_node in new_vrn:
            p_node = network.successors(r_node)[0]
            network.node[p_node]['dist'] = n
            node_type = network.node[p_node]['type']
            if node_type == 'pf':
                pf += 1
            if node_type == 'pr':
                pr += 1

        # Expand the valid compound set
        # When n = 0, this means the starting compounds
        # When n > 0, the valid compound set will be expanded
        new_vcn = expand_valid_compound_set(network, proc_num, new_vrn, \
        valid_compound_nodes) - valid_compound_nodes

        valid_compound_nodes = new_vcn.union(valid_compound_nodes)

        new_vrn = find_valid_reactant_nodes(network, proc_num, \
        valid_compound_nodes) - valid_reactant_nodes

        valid_reactant_nodes = new_vrn.union(valid_reactant_nodes)

        for node in new_vcn:
            network.node[node]['dist'] = n
            c += 1
        for node in new_vrn:
            network.node[node]['dist'] = n
            node_type = network.node[node]['type']
            if node_type == 'rf':
                rf += 1
            if node_type == 'rr':
                rr += 1

        # Nicely formatted progress output
        output = ''.join([
            '{0:<', str(l+6),
            '} {1:>', str(l+5),
            '} {2:>', str(l+5),
            '} {3:>', str(l+5),
            '} {4:>', str(l+5),
            '} {5:>', str(l+5),
            '}'
        ])
        print(output.format('Step ' + str(n) + ':', str(c) + ' c', str(rf) + \
        ' rf', str(pf) + ' pf', str(rr) + ' rr', str(pr) + ' pr'))

        n += 1

        no_new_vrn = set(prev_vrn) == valid_reactant_nodes
        no_new_vcn = set(prev_vcn) == valid_compound_nodes

        if no_new_vrn and no_new_vcn:
            # When no new valid compound or reactant nodes have been identified,
            # it is time to stop
            break
        else:
            if n > N and N !=-1:
                # n starts at 0 and increments by one before each round dealing
                # with that particular step n - stop when n exceeds the limit
                break
        prev_vrn = list(valid_reactant_nodes)
        prev_vcn = list(valid_compound_nodes)

    total_time = time.time() - time_start

    s_out("\nDone in %ss.\n" %str(total_time))

    return (valid_compound_nodes, valid_reactant_nodes)

def test_distance_to_origin():
    G = nx.DiGraph()
    G.add_node(1,type='c',start=False)
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
    G.add_node(8,type='c',start=True) # Compound 8 is now the start
    G.add_path([1,6,7,8])
    G.add_edge(5,6)
    G.add_node(11,type='rr',c={8})
    G.add_node(12,type='pr',c={1,5})
    G.add_path([8,11,12,1])
    G.add_edge(12,5)

    # Creating a copy, since the function modifies the network it is given
    output_0 = distance_to_origin(G.copy(), 4, 0)
    output_1 = distance_to_origin(G.copy(), 4, 1)
    # Not creating a copy in order to check modification capabilities
    output_2 = distance_to_origin(G, 4, 2)

    # Compound and reactant nodes reachable within 0 reaction steps
    assert output_0 == (set([8]), set([11]))

    assert output_1 == (set([8,1,5]), set([11,6,2]))
    assert output_2 == (set([8,1,5,4]), set([11,6,2,9]))

    assert set([G.node[n]['dist'] for n in [8,11]]) == set([0])
    assert set([G.node[n]['dist'] for n in [1,2,5,6,12]]) == set([1])
    assert set([G.node[n]['dist'] for n in [3,4,7,9]]) == set([2])

    # Test for detection of non-reachable nodes
    # Only nodes 1, 2, 3, 4, 5 and 6 should be reachable
    # and will receive a 'dist' value - other nodes do not
    Y = nx.DiGraph()
    Y.add_node(1,type='c',start=True)
    Y.add_node(2,type='rf',c={1})
    Y.add_node(3,type='pf',c={6})
    Y.add_node(6,type='c',start=False)
    Y.add_path([1,2,3,6])
    Y.add_node(4,type='rr',c={6})
    Y.add_node(5,type='pr',c={1})
    Y.add_path([6,4,5,1])

    Y.add_node(7,type='c',start=False)
    Y.add_node(8,type='rf',c={6,7})
    Y.add_node(9,type='pf',c={12})
    Y.add_node(12,type='c',start=False)
    Y.add_path([7,8,9,12])
    Y.add_edge(6,8)
    Y.add_node(10,type='rr',c={12})
    Y.add_node(11,type='pr',c={6,7})
    Y.add_path([6,8,9,12])
    Y.add_edge(12,10)

    output_Y = distance_to_origin(Y, 4, -1)

    z = 0
    for node in Y.nodes():
        try:
            x = Y.node[node]['dist']
        except KeyError:
            z += 1

    assert z == 6
    assert [Y.node[n]['dist'] for n in range(1,7)] == [0,0,1,1,2,1]


def prune_network(network, remove_cfm=True):
    """
    Remove all nodes that are 'unreachable' (lacking a 'dist' data key).

    Also removes CFM spectra from the mine_data dictionary by default.
    """
    for node in network.nodes():
        try:
            x = network.node[node]['dist']
        except KeyError:
            network.remove_node(node)
    if remove_cfm:
        for mid in network.graph['mine_data'].keys():
            if 'Neg_CFM_spectra' in network.graph['mine_data'][mid].keys():
                del network.graph['mine_data'][mid]['Neg_CFM_spectra']
            if 'Pos_CFM_spectra' in network.graph['mine_data'][mid].keys():
                del network.graph['mine_data'][mid]['Pos_CFM_spectra']
    network.graph['pruned'] = True

def test_prune_network():
    # Nodes that were not reached in distance_to_origin
    # are going to lack the 'dist' data key
    # Such nodes are expected to be removed
    G = nx.DiGraph()
    G.add_node(1,type='c',start=False)
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
    G.add_node(8,type='c',start=True) # Compound 8 is now the start
    G.add_path([1,6,7,8])
    G.add_edge(5,6)
    G.add_node(11,type='rr',c={8})
    G.add_node(12,type='pr',c={1,5})
    G.add_path([8,11,12,1])
    G.add_edge(12,5)

    G.graph['mine_data'] = {
    'C1':{'_id':'C1','Neg_CFM_spectra':{'Dummy'},'Pos_CFM_spectra':{'Dummy'}},
    'C4':{'_id':'C4','Neg_CFM_spectra':{'Dummy'},'Pos_CFM_spectra':{'Dummy'}},
    'C5':{'_id':'C5','Neg_CFM_spectra':{'Dummy'},'Pos_CFM_spectra':{'Dummy'}},
    'C8':{'_id':'C1','Neg_CFM_spectra':{'Dummy'},'Pos_CFM_spectra':{'Dummy'}}
    }

    H = G.copy()

    output = distance_to_origin(G, 4, 1)
    output = distance_to_origin(H, 4, 2)

    Y = G.subgraph([1,2,5,6,8,11,12])
    Z = H.subgraph([1,2,3,4,5,6,7,8,9,11,12])

    prune_network(G)
    prune_network(H)

    assert nx.is_isomorphic(G,Y)
    assert G.nodes(data=True) == Y.nodes(data=True)
    assert set(G.edges()) == set(Y.edges())

    assert nx.is_isomorphic(H,Z)
    assert H.nodes(data=True) == Z.nodes(data=True)
    assert set(H.edges()) == set(Z.edges())

    assert G.graph['mine_data'] == H.graph['mine_data'] == {
    'C1':{'_id':'C1'},
    'C4':{'_id':'C4'},
    'C5':{'_id':'C5'},
    'C8':{'_id':'C1'}
    }


def prepare_dictionaries(network):
    """Prepare dictionaries for translating KEGG IDs and Names to MINE IDs."""
    network.graph['kegg2nodes'] = {}
    network.graph['name2nodes'] = {}
    for node in network.nodes():
        # Only consider nodes that are reachable, i.e. have predecessor(s)
        if not network.predecessors(node):
            continue
        # Get associated KEGG IDs
        try:
            mid = network.node[node]['mid']
            kegg_ids = network.graph['mine_data'][mid]['DB_links']['KEGG']
        except KeyError:
            kegg_ids = []
        # Get associated names
        try:
            mid = network.node[node]['mid']
            names = network.graph['mine_data'][mid]['Names']
        except KeyError:
            names = []
        # Add node to set of nodes for each KEGG ID
        for kegg_id in kegg_ids:
            try:
                network.graph['kegg2nodes'][kegg_id].add(node)
            except KeyError:
                network.graph['kegg2nodes'][kegg_id] = set([node])
        # Add node to set of nodes for each name
        for name in names:
            try:
                network.graph['name2nodes'][name].add(node)
            except KeyError:
                network.graph['name2nodes'][name] = set([node])

def test_prepare_dictionaries():
    G = nx.DiGraph()
    G.graph['mine_data'] = {
        'C1' : {
            '_id':'C1', 'Names':['Something','Anything'],
            'DB_links':{'KEGG':['C93102']}
        },
        'X1' : {
            '_id':'X1', 'Names':['Something','Anything'],
            'DB_links':{'KEGG':['C93102']}
        },
        'C2' : {
            '_id':'C2', 'Names':['Whatever'],
            'DB_links':{'KEGG':['C33391','C33392']}
        },
        'C3' : {'_id':'C3','Names':['Bit','Bob']},
        'C4' : {'_id':'C4'},
        'X5' : {'_id':'X5','Names':['Whatever']},
        'C7' : {'_id':'C7','DB_links':{'KEGG':['C00001']}},
        'C8' : {'_id':'C8','Names':['Twin'],'DB_links':{'KEGG':['C11011']}},
        'C9' : {'_id':'C9','Names':['Twin'],'DB_links':{'KEGG':['C11011']}},
        'R1' : {'_id':'R1'}, 'R2' : {'_id':'R2'}, 'R3' : {'_id':'R3'},
        'R4' : {'_id':'R4'}, 'R8' : {'_id':'R8'}, 'R9' : {'_id':'R9'}
    }

    # Add nodes
    nodes = [
        (1,{'mid':'C1'}), (2,{'mid':'X1'}), (3,{'mid':'C2'}),
        (4,{'mid':'C3'}), (5,{'mid':'C4'}), (6,{'mid':'X5'}),
        (7,{'mid':'C7'}), (8,{'mid':'C8'}), (9,{'mid':'C9'}),
        (11,{'mid':'R1'}), (12,{'mid':'R2'}), (13,{'mid':'R3'}),
        (14,{'mid':'R4'}), (17,{'mid':'R7'}), (18,{'mid':'R8'}),
        (19,{'mid':'R9'})
    ]
    G.add_nodes_from(nodes)

    # Connect nodes to simulate connectedness
    # Note that (7,17) means 7 is only connected outwards and cannot be reached
    G.add_edges_from([(11,1),(12,3),(13,4),(14,5),(7,17),(18,8),(19,9)])

    expected_kegg2nodes = {
        'C93102':set([1]), 'C33391':set([3]),
        'C33392':set([3]), 'C11011':set([8,9])
    }

    expected_name2nodes = {
        'Something':set([1]), 'Anything':set([1]), 'Whatever':set([3]),
        'Bit':set([4]), 'Bob':set([4]), 'Twin':set([8,9])
    }

    prepare_dictionaries(G)

    assert G.graph['kegg2nodes'] == expected_kegg2nodes
    assert G.graph['name2nodes'] == expected_name2nodes


def create_SMILES_to_KEGG_dict(KEGG_dict):
    """Create a dictionary for translating SMILES to KEGG IDs."""
    SMILES_to_KEGG = {}
    for KEGG_id in KEGG_dict.keys():
        KEGG_comp = KEGG_dict[KEGG_id]
        try:
            SMILES = KEGG_comp['SMILES']
        except KeyError:
            continue
        try:
            SMILES_to_KEGG[SMILES].add(KEGG_id)
        except KeyError:
            SMILES_to_KEGG[SMILES] = {KEGG_id}
    return SMILES_to_KEGG

def test_create_SMILES_to_KEGG_dict():
    KEGG_dict = {
        'C80000':{'SMILES':'OCC1OC(O)C(O)C(O)C1O'},
        'C80001':{'SMILES':'CCO'},
        'C80002':{'SMILES':'CCO'},
        'C80003':{'_id':'C80002'}
    }
    exp_SMILES_to_KEGG = {
        'OCC1OC(O)C(O)C(O)C1O':{'C80000'},
        'CCO':{'C80001', 'C80002'}
    }
    assert exp_SMILES_to_KEGG == create_SMILES_to_KEGG_dict(KEGG_dict)


def MINE_comps_KEGG_filter(comps, SMILES_to_KEGG):
    """Remove compounds that have no or cannot be assigned a KEGG ID."""
    comps_filtered = []
    for comp in comps:
        try:
            KEGG_ids = comp['DB_links']['KEGG']
        except KeyError:
            KEGG_ids = []
        try:
            SMILES = comp['SMILES']
        except KeyError:
            SMILES = ''
        if KEGG_ids:
            comps_filtered.append(comp)
            continue
        if SMILES in SMILES_to_KEGG.keys():
            KEGG_ids = sorted(list(SMILES_to_KEGG[SMILES]))
            try:
                comp['DB_links']['KEGG'] = KEGG_ids
            except KeyError:
                comp['DB_links'] = {'KEGG' : KEGG_ids}
            comps_filtered.append(comp)
    return comps_filtered

def test_MINE_comps_KEGG_filter():
    comps = [
        {'_id':'C1', # Should be kept
        'DB_links':{'KEGG':['C00001']},
        'Reactant_in':['R1','R2'],
        'Product_of':['R3','R4'],
        'SMILES':'O=CO'},
        {'_id':'C2', # Should be kept
        'DB_links':{'KEGG':['C00002', 'C00003']},
        'Reactant_in':['R4'],
        'SMILES':'CCCOC'},
        {'_id':'C3', # Should be kept (via SMILES_to_KEGG)
        'DB_links':{'CUD':['CUD_5']},
        'Reactant_in':['R5'],
        'Product_of':['R6', 'R7'],
        'SMILES':'CC(N)OC'},
        {'_id':'C4', # Should be kept (via SMILES_to_KEGG)
        'SMILES':'C(C)C=(O)OP=([O])(O)[O-]'},
        {'_id':'C5', # Should be removed
        'Product_of':['R9']},
        {'_id':'C6', # Should be removed
        'DB_links':{'CUD':['CUD_8']},
        'Product_of':['R10']}
    ]
    SMILES_to_KEGG = {
        'CC(N)OC' : {'C00004'},
        'C(C)C=(O)OP=([O])(O)[O-]' : {'C00005', 'C00006'}
    }
    exp_comps = comps[:4]
    exp_comps[2]['DB_links'] = {'CUD':['CUD_5'], 'KEGG':['C00004']}
    exp_comps[3]['DB_links'] = {'KEGG':['C00005','C00006']}
    assert MINE_comps_KEGG_filter(comps, SMILES_to_KEGG) == exp_comps

def operators_identical(op1, op2):
    """Checks if two BNICE EC operators are eachother's reverse/identical."""
    for e1, e2 in zip(op1.split('.'), op2.split('.')):
        if e1.lstrip('-') != e2.lstrip('-'):
            return False
    return True

def test_operators_identical():
    assert operators_identical('1.1.1.-', '1.1.1.-')
    assert operators_identical('2.42.1.a', '2.42.-1.a')
    assert not operators_identical('1.2.3.-', '1.2.4.-')


def extract_ints(str_list):
    int_list = []
    for s in str_list:
        try:
            int_list.append(int(s))
        except ValueError:
            continue
    return int_list

def test_extract_ints():
    L1 = ['1','2','23','c']
    L2 = ['4','42','-1','-']
    L3 = ['1','1','1','22']
    assert extract_ints(L1) == [1,2,23]
    assert extract_ints(L2) == [4,42,-1]
    assert extract_ints(L3) == [1,1,1,22]


def remove_redundant_MINE_rxns(rxns):
    """Removes identical but reversed MINE reactions."""

    discarded_rxns = set()

    # Set up progress indicator
    p = Progress(max_val = len(rxns)**2, design='pt')
    n = 0

    # Iterate over all reaction pairs
    for rp in product(enumerate(rxns), enumerate(rxns)):
        # Report progress
        n += 1
        s_out("\rRemoving redundant MINE reactions... %s" % p.to_string(n))

        # Don't compare a reaction to itself
        if rp[0][0] == rp[1][0]:
            continue

        # Don't perform further comparisons for discarded reactions
        if set([rp[0][0], rp[1][0]]).intersection(discarded_rxns):
            continue

        # Compare the Products and reactants
        try:
            r1p = set([tuple(c) for c in rp[0][1]['Products']])
        except KeyError:
            r1p = set()
        try:
            r2p = set([tuple(c) for c in rp[1][1]['Products']])
        except KeyError:
            r2p = set()
        try:
            r1r = set([tuple(c) for c in rp[0][1]['Reactants']])
        except KeyError:
            r1r = set()
        try:
            r2r = set([tuple(c) for c in rp[1][1]['Reactants']])
        except KeyError:
            r2r = set()

        if r1p == r2r and r2p == r1r:
            are_mutual_reverse = True
        else:
            are_mutual_reverse = False

        # Compare the sets of operators
        ops1 = rp[0][1]['Operators']
        ops2 = rp[1][1]['Operators']
        n_identical = 0
        for op_pair in product(ops1, ops2):
            if operators_identical(*op_pair):
                n_identical += 1

        # If the reactions have the same operators and are eachothers reverse,
        # determine which one to discard
        if n_identical == len(ops1) == len(ops2) and are_mutual_reverse:

            # One reaction needs to be discarded; Find which one to keep
            ops1_int = (extract_ints(op.split('.')) for op in ops1)
            ops2_int = (extract_ints(op.split('.')) for op in ops2)
            ops1_neg = sum(any(x < 0 for x in op) for op in ops1_int)
            ops2_neg = sum(any(x < 0 for x in op) for op in ops2_int)

            # Discard the reaction with most negative operators
            if ops1_neg < ops2_neg:
                discarded_rxns.add(rp[1][0])

            elif ops1_neg > ops2_neg:
                discarded_rxns.add(rp[0][0])

            # Otherwise discard the second reaction
            else:
                if rp[0][0] < rp[1][0]:
                    discarded_rxns.add(rp[1][0])
                else:
                    discarded_rxns.add(rp[0][0])

    # Return reactions that were not discarded
    print("")
    return [rxns[i] for i in range(len(rxns)) if i not in discarded_rxns]

def test_remove_redundant_MINE_rxns():
    rxns = [
        {'_id':'R1', # Should be removed; reverse of R2
        'Operators':['2.7.-1.a'],
        'Products':[[1,'C1'],[1,'X1']],
        'Reactants':[[1,'X2'],[1,'C2']]},
        {'_id':'R2', # Should remain
        'Operators':['2.7.1.a'],
        'Reactants':[[1,'C1'],[1,'X1']],
        'Products':[[1,'C2'],[1,'X2']]},
        {'_id':'R3', # Should remain; different Reactants
        'Operators':['2.7.1.a'],
        'Reactants':[[1,'C3'],[1,'X1']],
        'Products':[[1,'C2'],[1,'X2']]},
        {'_id':'R4', # Should remain; different Operators
        'Operators':['3.12.2.-'],
        'Reactants':[[1,'C1'],[1,'X1']],
        'Products':[[1,'C2'],[1,'X2']]},
        {'_id':'R5', # Should remain; different Operators
        'Operators':['2.7.-1.a', '4.1.1.g2'],
        'Products':[[1,'C1'],[1,'X1']],
        'Reactants':[[1,'C2'],[1,'X2']]},
        {'_id':'R6', # Should remain
        'Operators':['1.1.1.a', '2.2.2.-'],
        'Reactants':[[1,'C4']],
        'Products':[[2,'C5'],[1,'X3']]},
        {'_id':'R7', # Should be removed; reverse of R6
        'Operators':['2.2.-2.-', '1.1.-1.a'],
        'Products':[[1,'C4']],
        'Reactants':[[2,'C5'],[1,'X3']]},
        {'_id':'R8', # Should remain
        'Operators':['1.1.1.a', '2.2.2.-', '3.3.-3.g'],
        'Reactants':[[1,'C4']],
        'Products':[[2,'C5'],[1,'X4']]},
        {'_id':'R9', # Should be removed; reverse of R8 and most 'negative' ops.
        'Operators':['2.2.-2.-', '1.1.-1.a', '3.3.3.g'],
        'Products':[[1,'C4']],
        'Reactants':[[2,'C5'],[1,'X4']]},
        {'_id':'R10', # Should remain
        'Operators':['1.2.3.a', '4.5.-6.-'],
        'Reactants':[[1,'C4']],
        'Products':[[2,'C5'],[1,'X3']]},
        {'_id':'R11', # Should be removed; reverse of R11 and second encountered
        'Operators':['1.2.-3.a', '4.5.6.-'],
        'Products':[[1,'C4']],
        'Reactants':[[2,'C5'],[1,'X3']]},
        {'_id':'R12', # Should be removed; reverse of R11 and third encountered
        'Operators':['1.2.-3.a', '4.5.6.-'],
        'Products':[[1,'C4']],
        'Reactants':[[2,'C5'],[1,'X3']]},
        {'_id':'R13', # Should be removed; reverse of R14
        'Operators':['2.8.-1.a'],
        'Products':[[1,'C8'],[1,'X8']],
        'Reactants':[[1,'C9'],[1,'X9']]},
        {'_id':'R14', # Should remain
        'Operators':['2.8.1.a'],
        'Reactants':[[1,'C8'],[1,'X8']],
        'Products':[[1,'C9'],[1,'X9']]},
    ]
    exp_rxns = [rxns[1]] + rxns[2:6] + [rxns[7]] + [rxns[9]] + [rxns[-1]]
    filtered_rxns = remove_redundant_MINE_rxns(rxns)
    for exp_rxn in exp_rxns:
        assert exp_rxn in filtered_rxns
    for rxn in filtered_rxns:
        assert rxn in exp_rxns
    assert exp_rxns == filtered_rxns


def remove_non_KEGG_MINE_rxns(rxns, comps):
    """Return reactions with only KEGG compounds."""
    filtered_rxns = []
    allowed_ids = set([c['_id'] for c in comps])
    p = Progress(max_val = len(rxns), design = 'p')
    n = 0
    for rxn in rxns:
        n += 1
        s_out("\rRemoving non-KEGG MINE reactions... %s" % p.to_string(n))
        if set(extract_reaction_comp_ids(rxn)).issubset(allowed_ids):
            filtered_rxns.append(rxn)
    print("")
    return filtered_rxns

def test_remove_non_KEGG_MINE_rxns():
    rxns = [
        {'_id':'R1', # All KEGG; keep
        'Reactants':[[1,'C1'],[1,'C2']],
        'Products':[[1,'C3'],[1,'C4']]},
        {'_id':'R2', # One non-KEGG reactant; remove
        'Reactants':[[1,'C5'],[1,'X1']],
        'Products':[[1,'X2'],[2,'C1']]},
        {'_id':'R3', # One non-KEGG product; remove
        'Reactants':[[1,'C3'],[1,'X3']],
        'Products':[[1,'X4'],[2,'C6']]},
        {'_id':'R4', # All non-KEGG; remove
        'Reactants':[[1,'C8'],[1,'X8']],
        'Products':[[1,'X9'],[1,'C9']]},
        {'_id':'R5', # All KEGG; keep
        'Reactants':[[1,'C3'],[1,'C2']],
        'Products':[[1,'C3'],[1,'X2']]},
        {'_id':'R6', # All KEGG; keep
        'Reactants':[[2,'C1']],
        'Products':[[1,'X4'],[1,'C2']]}
    ]
    comps = [
        {'_id':'C1'}, {'_id':'C2'}, {'_id':'C3'}, {'_id':'C4'},
        {'_id':'X1'}, {'_id':'X2'}, {'_id':'X3'}, {'_id':'X4'},
    ]
    assert remove_non_KEGG_MINE_rxns(rxns, comps) == [rxns[0]] + rxns[4:]


def KEGG_rxns_from_MINE_rxns(rxns, comps):
    """Produce reactions with KEGG IDs instead of MINE IDs."""

    # Set up a MINE ID to KEGG ID dictionary
    M2K = dict([(c['_id'], c['DB_links']['KEGG']) for c in comps])

    # Set up a MINE ID to Cofactor status dictionary
    def is_not_connected(comp):
        if set(['Reactant_in','Product_of']).intersection(comp.keys()):
            return False
        else:
            return True

    is_cof = dict([(c['_id'], is_not_connected(c)) for c in comps])

    # Produce KEGG-labeled reactions
    KEGG_rxns = []
    p = Progress(max_val=len(rxns), design='pt')
    n = 0
    for rxn in rxns:
        # Report progress
        n += 1
        s_out("\rProducing MINE reactions with KEGG IDs... %s" % p.to_string(n))

        # Determine positions of cofactors
        r_cofacs = []
        for e in enumerate([is_cof[c[1]] for c in rxn['Reactants']]):
            if e[1]:
                r_cofacs.append(e[0])
        p_cofacs = []
        for e in enumerate([is_cof[c[1]] for c in rxn['Products']]):
            if e[1]:
                p_cofacs.append(e[0])

        # Create combinations of KEGG IDs
        r_combos = product(*[M2K[c[1]] for c in rxn['Reactants']])
        p_combos = product(*[M2K[c[1]] for c in rxn['Products']])
        r_coeffs = [c[0] for c in rxn['Reactants']]
        p_coeffs = [c[0] for c in rxn['Products']]
        rp_pairs = product(r_combos, p_combos)

        # Iterate over the new reactant/product pairs and produce reactions
        n = 0

        for rp_pair in rp_pairs:

            # Create a copy of the old reaction
            new_rxn = deepcopy(rxn)

            # Add an index to the end of the ID
            new_rxn['_id'] = rxn['_id'] + '_' + str(n)

            # Update Reactants
            new_rxn['Reactants'] = [list(i) for i in zip(r_coeffs, rp_pair[0])]

            # Update Products
            new_rxn['Products'] = [list(i) for i in zip(p_coeffs, rp_pair[1])]

            # Create RPair dictionary
            rp_dict = {}

            # Define a function for building up the RPair entries
            def build_RPair_dict(rp_dict, en_comp, cof_indices):
                # The compound could be a cofactor...
                if en_comp[0] in cof_indices:
                    try:
                        entry = rp_dict['cofac'][0] + '_' + en_comp[1]
                        rp_dict['cofac'][0] = entry
                    except KeyError:
                        rp_dict['cofac'] = [en_comp[1], 'cofac']

                # ...or a main participant
                else:
                    try:
                        entry = rp_dict['main'][0] + '_' + en_comp[1]
                        rp_dict['main'][0] = entry
                    except KeyError:
                        rp_dict['main'] = [en_comp[1], 'main']

                # Return the expanded RPair dictionary
                return rp_dict

            # Build up RPair entries with reactants
            for en_r in enumerate(rp_pair[0]):
                build_RPair_dict(rp_dict, en_r, r_cofacs)

            # Build up RPair entries with products
            for en_p in enumerate(rp_pair[1]):
                build_RPair_dict(rp_dict, en_p, p_cofacs)

            # Convert to tuples
            for item in rp_dict.items():
                rp_dict[item[0]] = tuple(item[1])

            # Add the RPair entries to the reaction
            new_rxn['RPair'] = rp_dict

            # Finally, add the reaction to the list of new reactions
            KEGG_rxns.append(new_rxn)

            # ...and count one step up
            n += 1

    return KEGG_rxns

def test_KEGG_rxns_from_MINE_rxns():
    rxns = [
        {'_id':'R1', 'Operators':['1.1.1.a'],
        'Reactants':[[1,'C1'],[1,'C2']],
        'Products':[[1,'C3'],[1,'C4']]},
        {'_id':'R5', 'Operators':['1.2.1.a'],
        'Reactants':[[1,'X3'],[1,'C2']],
        'Products':[[1,'C3'],[1,'X2']]},
        {'_id':'R6', 'Operators':['1.3.1.a'],
        'Reactants':[[2,'C1'],[1,'X1']],
        'Products':[[1,'X4'],[1,'C2']]},
        {'_id':'R9', 'Operators':['1.1.4.a', '3.4.2.-'],
        'Reactants':[[1,'C2']],
        'Products':[[1,'C4']]}
    ]
    comps = [
        {'_id':'C1', 'DB_links':{'KEGG':['C10000','C10001']},
        'Reactant_in':['R1','R6']},
        {'_id':'C2', 'DB_links':{'KEGG':['C20000']},
        'Reactant_in':['R1','R5','R9'], 'Product_of':['R6']},
        {'_id':'C3', 'DB_links':{'KEGG':['C30000','C30001']},
        'Product_of':['R1','R5']},
        {'_id':'C4', 'DB_links':{'KEGG':['C40000']},
        'Product_of':['R1','R9']},
        {'_id':'X1', 'DB_links':{'KEGG':['C10000']}},
        {'_id':'X2', 'DB_links':{'KEGG':['C00002']}},
        {'_id':'X3', 'DB_links':{'KEGG':['C00003']}},
        {'_id':'X4', 'DB_links':{'KEGG':['C40000','C40001']}},
    ]
    exp_rxns = [
        {'_id':'R1_0', 'Operators':['1.1.1.a'],
        'Reactants':[[1,'C10000'],[1,'C20000']],
        'Products':[[1,'C30000'],[1,'C40000']],
        'RPair':{
            'main':('C10000_C20000_C30000_C40000','main')
        }},
        {'_id':'R1_1', 'Operators':['1.1.1.a'],
        'Reactants':[[1,'C10000'],[1,'C20000']],
        'Products':[[1,'C30001'],[1,'C40000']],
        'RPair':{
            'main':('C10000_C20000_C30001_C40000','main')
        }},
        {'_id':'R1_2', 'Operators':['1.1.1.a'],
        'Reactants':[[1,'C10001'],[1,'C20000']],
        'Products':[[1,'C30000'],[1,'C40000']],
        'RPair':{
            'main':('C10001_C20000_C30000_C40000','main')
        }},
        {'_id':'R1_3', 'Operators':['1.1.1.a'],
        'Reactants':[[1,'C10001'],[1,'C20000']],
        'Products':[[1,'C30001'],[1,'C40000']],
        'RPair':{
            'main':('C10001_C20000_C30001_C40000','main')
        }},
        {'_id':'R5_0', 'Operators':['1.2.1.a'],
        'Reactants':[[1,'C00003'],[1,'C20000']],
        'Products':[[1,'C30000'],[1,'C00002']],
        'RPair':{
            'main':('C20000_C30000','main'),
            'cofac':('C00003_C00002', 'cofac')
        }},
        {'_id':'R5_1', 'Operators':['1.2.1.a'],
        'Reactants':[[1,'C00003'],[1,'C20000']],
        'Products':[[1,'C30001'],[1,'C00002']],
        'RPair':{
            'main':('C20000_C30001','main'),
            'cofac':('C00003_C00002', 'cofac')
        }},
        {'_id':'R6_0', 'Operators':['1.3.1.a'],
        'Reactants':[[2,'C10000'],[1,'C10000']],
        'Products':[[1,'C40000'],[1,'C20000']],
        'RPair':{
            'main':('C10000_C20000','main'),
            'cofac':('C10000_C40000', 'cofac')
        }},
        {'_id':'R6_1', 'Operators':['1.3.1.a'],
        'Reactants':[[2,'C10000'],[1,'C10000']],
        'Products':[[1,'C40001'],[1,'C20000']],
        'RPair':{
            'main':('C10000_C20000','main'),
            'cofac':('C10000_C40001', 'cofac')
        }},
        {'_id':'R6_2', 'Operators':['1.3.1.a'],
        'Reactants':[[2,'C10001'],[1,'C10000']],
        'Products':[[1,'C40000'],[1,'C20000']],
        'RPair':{
            'main':('C10001_C20000','main'),
            'cofac':('C10000_C40000', 'cofac')
        }},
        {'_id':'R6_3', 'Operators':['1.3.1.a'],
        'Reactants':[[2,'C10001'],[1,'C10000']],
        'Products':[[1,'C40001'],[1,'C20000']],
        'RPair':{
            'main':('C10001_C20000','main'),
            'cofac':('C10000_C40001', 'cofac')
        }},
        {'_id':'R9_0', 'Operators':['1.1.4.a', '3.4.2.-'],
        'Reactants':[[1,'C20000']],
        'Products':[[1,'C40000']],
        'RPair':{
            'main':('C20000_C40000','main')
        }}
    ]
    new_rxns = KEGG_rxns_from_MINE_rxns(rxns, comps)
    assert len(exp_rxns) == len(new_rxns)
    assert exp_rxns == new_rxns


def add_MINE_rxns_to_KEGG_comps(comps, rxns):
    """Add MINE reactions to KEGG compound 'Reactant_in'/'Product_of' lists."""

    # Initialize list of updated compounds
    new_comps = []

    # Iterate over the KEGG compounds
    p = Progress(max_val=len(comps), design='pt')
    n = 0
    for comp in comps:
        # Report progress
        n += 1
        s_out("\rAdding MINE reactions to KEGG compounds... %s" % p.to_string(n))

        # Initialize a new compound
        new_comp = deepcopy(comp)

        # Check if the compound should list any of the reactions
        for rxn in rxns:

            # Check if the compound may list the reaction
            if not allow_reaction_listing(comp, rxn):
                continue

            # Add the reaction to Reactant_in if applicable
            if comp['_id'] in [c[1] for c in rxn['Reactants']]:
                try:
                    new_comp['Reactant_in'].append(rxn['_id'])
                except KeyError:
                    new_comp['Reactant_in'] = [rxn['_id']]

            # Add the reaction to Product_of if applicable
            if comp['_id'] in [c[1] for c in rxn['Products']]:
                try:
                    new_comp['Product_of'].append(rxn['_id'])
                except KeyError:
                    new_comp['Product_of'] = [rxn['_id']]

        # Add the updated compound to the new compounds list
        new_comps.append(new_comp)

    return new_comps

def test_add_MINE_rxns_to_KEGG_comps():
    comps = [
        {'_id':'C00001', 'Reactant_in':['R00001']},
        {'_id':'C00002', 'Product_of':['R00001'], 'Reactant_in':['R00002']},
        {'_id':'C00003'},
        {'_id':'C00004', 'Product_of':['R00004']}
    ]
    rxns = [
        {'_id':'R1_0',
        'Reactants':[[1,'C00002'],[1,'C00004']],
        'Products':[[1,'C00003'],[1,'C00001']],
        'RPair':{
            'main':('C00002_C00003','main'),
            'cofac':('C00004_C00001','cofac')
        }},
        {'_id':'R2_2',
        'Reactants':[[1,'C00003']],
        'Products':[[1,'C00001']],
        'RPair':{
            'main':('C00003_C00001','main')
        }},
        {'_id':'R3_12',
        'Reactants':[[1,'C00001'],[1,'C00003'],[1,'C00004']],
        'Products':[[1,'C00002']],
        'RPair':{
            'main':('C00001_C00003_C00004_C00002','main')
        }}
    ]
    exp_comps = [
        {'_id':'C00001',
            'Reactant_in':['R00001','R3_12'],
            'Product_of':['R2_2']},
        {'_id':'C00002',
            'Product_of':['R00001','R3_12'],
            'Reactant_in':['R00002','R1_0']},
        {'_id':'C00003',
            'Reactant_in':['R2_2','R3_12'],
            'Product_of':['R1_0']},
        {'_id':'C00004',
            'Product_of':['R00004'],
            'Reactant_in':['R3_12']}
    ]
    new_comps = add_MINE_rxns_to_KEGG_comps(comps, rxns)
    assert len(exp_comps) == len(new_comps)
    assert exp_comps == new_comps


def enhance_KEGG_with_MINE(KEGG_comp_dict, KEGG_rxn_dict):
    """Enhance a raw KEGG reaction network with MINE reactions."""

    print("\nEnhancing KEGG reaction network data with MINE reactions...\n")

    # Set up MINE connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    # Create a list of KEGG compound IDs
    KEGG_comp_ids = list(KEGG_comp_dict.keys())

    # Create a list of IDs for MINE compounds to download
    print("Identifying KEGG compounds in MINE...")
    MINE_comp_ids = set()
    for MINE_query_result in threaded_quicksearch(con, db, KEGG_comp_ids):
        try:
            MINE_comp_ids.add(MINE_query_result['_id'])
        except KeyError:
            continue

    # Download the MINE compounds corresponding to KEGG IDs
    print("\nDownloading MINE compounds...")
    MINE_comps = threaded_getcomps(con, db, list(MINE_comp_ids))
    print("")

    # Create a list of IDs for MINE reactions to download
    MINE_rxn_ids = set()
    p = Progress(max_val=len(MINE_comps), design='p')
    n = 0
    for comp in MINE_comps:
        n += 1
        s_out("\rIdentifying MINE reactions to download... %s" % p.to_string(n))
        for MINE_rxn_id in extract_comp_reaction_ids(comp):
            MINE_rxn_ids.add(MINE_rxn_id)
    s_out("\rIdentifying MINE reactions to download... Done. \n")

    # Download the MINE reactions listed for the MINE compounds
    s_out("Downloading MINE reactions...\n")
    MINE_rxns = list(filter(None, threaded_getrxn(con, db, list(MINE_rxn_ids))))

    # Download the 'X' MINE compounds listed for the reactions
    s_out("\nIdentifying cofactors...")
    MINE_X_comp_ids = set()
    for rxn in MINE_rxns:
        cids = set(extract_reaction_comp_ids(rxn))
        for cid in cids:
            if cid.startswith('X'):
                MINE_X_comp_ids.add(cid)
    s_out(" Done.\n")
    s_out("\nDownloading MINE cofactors...\n")
    MINE_comps = MINE_comps + threaded_getcomps(con, db, list(MINE_X_comp_ids))

    # Create a SMILES to KEGG dictionary from the KEGG compound dictionary
    s_out("\nSetting up SMILES to KEGG translation...")
    SMILES_to_KEGG = create_SMILES_to_KEGG_dict(KEGG_comp_dict)
    s_out(" Done.\n")

    # Filter the MINE compounds for KEGG equivalents
    s_out("Filtering MINE compounds...")
    MINE_comps = filter(None, MINE_comps) # Remove None results
    MINE_comps = MINE_comps_KEGG_filter(MINE_comps, SMILES_to_KEGG)
    s_out(" Done.\n")

    # Remove MINE reactions that have non-KEGG compounds
    MINE_rxns = remove_non_KEGG_MINE_rxns(MINE_rxns, MINE_comps)

    # Remove redundant MINE reactions
    MINE_rxns = remove_redundant_MINE_rxns(MINE_rxns)

    # Create a list of MINE reactions in terms of KEGG compounds
    MINE_rxns = KEGG_rxns_from_MINE_rxns(MINE_rxns, MINE_comps)

    # Add the new MINE reactions to the KEGG compounds
    KEGG_comps = list(KEGG_comp_dict.values())
    KEGG_comps = add_MINE_rxns_to_KEGG_comps(KEGG_comps, MINE_rxns)

    # Add the new MINE reactions to the KEGG reactions dictionary
    s_out("Updating raw KEGG reaction network data...")
    KEGG_rxn_dict.update(dict([(rxn['_id'], rxn) for rxn in MINE_rxns]))

    # Create a new KEGG compound dictionary
    KEGG_comp_dict = dict([(comp['_id'], comp) for comp in KEGG_comps])
    s_out(" Done.\n")

    # Return the new KEGG compound dict. and the updated KEGG reaction dict.
    return (KEGG_comp_dict, KEGG_rxn_dict)

def test_enhance_KEGG_with_MINE():
    # THE TESTING SYSTEM
    # Compounds:
    # C00001    H2O
    # C00003    NAD+
    # C00004    NADH
    # C00011    CO2
    # C00080    H+
    # C01412    Butanal
    # C02804    5-Hydroxypentanoate
    # C06142    Butanol
    #
    # KEGG reaction R03544:
    # DEFINITION  Butanal + NADH + H+ <=> 1-Butanol + NAD+
    # EQUATION    C01412 + C00004 + C00080 <=> C06142 + C00003
    #
    # MINE reaction:
    # DEFINITION  5-Hydroxypentanoate <=> 1-Butanol + CO2
    # EQUATION    C02804 <=> C06142 + C00011

    # Set up MINE connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    # Download KEGG data
    KEGG_comp_ids = ['C00001', 'C00003', 'C00004', 'C00011',
                     'C00080', 'C01412', 'C02804', 'C06142']
    KEGG_rxn_ids = ['R03544']
    KEGG_comp_dict, KEGG_rxn_dict = get_raw_KEGG(KEGG_comp_ids, KEGG_rxn_ids)

    # This is the expected MINE data (butanol from butanal or hydroxypentanoate)
    MINE_rxns = [
    {'Rxn_Hash': '011bbc8f35061fb116cc4cbe46285b9db340031f3d3449f02ddf228666c686ac-5ef24b4a29298e80df3e410e68928b6e-b6560c09a82d56b15c28787b9dc1b017',
    'Products': [[1, 'C01412'], [1, 'C00004'], [1, 'C00080']],
    'Reactants': [[1, 'C06142'], [1, 'C00003']],
    'RPair': {
        'main': ('C06142_C01412', 'main'),
        'cofac': ('C00003_C00004_C00080', 'cofac')},
    'Operators': ['1.1.1.a'],
    '_id': 'R86e63dd5a75dedff25511e9535e77e2316e4c7af_0'
    },
    {'Rxn_Hash': '1b67a3e6edc9cb29d69e8c21af86cf42129d566f728ff721bfd50914cbae3357-5ef24b4a29298e80df3e410e68928b6e-19772c3a68c961ee657f9e206e21bd88',
    'Products': [[1, 'C06142'], [1, 'C00011']],
    'Reactants': [[1, 'C02804']],
    'RPair': {
        'cofac': ('C00011', 'cofac'),
        'main': ('C02804_C06142', 'main')},
    'Operators': ['4.1.1.k'],
    '_id': 'R766e3bcdbdf841f470c146b7ad9b74ca35d5c3e6_0'}
    ]

    # The expected KEGG reaction dictionary
    exp_KEGG_rxn_dict = deepcopy(KEGG_rxn_dict)
    exp_KEGG_rxn_dict[MINE_rxns[0]['_id']] = MINE_rxns[0]
    exp_KEGG_rxn_dict[MINE_rxns[1]['_id']] = MINE_rxns[1]

    # The expected KEGG compound dictionary
    exp_KEGG_comp_dict = deepcopy(KEGG_comp_dict)

    exp_KEGG_comp_dict['C06142']['Product_of'].append(MINE_rxns[1]['_id'])
    exp_KEGG_comp_dict['C06142']['Reactant_in'] = [MINE_rxns[0]['_id']]

    exp_KEGG_comp_dict['C01412']['Product_of'] = [MINE_rxns[0]['_id']]

    exp_KEGG_comp_dict['C02804']['Reactant_in'] = [MINE_rxns[1]['_id']]

    # Produce enhanced KEGG dictionaries
    enh_KEGG_comp_dict, enh_KEGG_rxn_dict = enhance_KEGG_with_MINE(
        KEGG_comp_dict, KEGG_rxn_dict
    )

    # Assert equality
    assert exp_KEGG_comp_dict == enh_KEGG_comp_dict
    assert exp_KEGG_rxn_dict == enh_KEGG_rxn_dict


# Main code block
def main(infile, mine, kegg, step_limit,
    comp_limit, C_limit, enhance, outfile_name):

    # Exit if a database choice has not been specified
    if not mine and not kegg:
        msg = (
        "\nPlease choose one or more databases for network construction:\n" +
        "-M, --mine for MINE\n" +
        "-K, --kegg for KEGG\n"
        )
        sys.exit(msg)

    # Get starting compounds
    if infile:
        start_kegg_ids = read_compounds(infile)
    else:
        start_kegg_ids = []

    # Acquire raw KEGG dictionaries
    kegg_comp_dict = {} # Default
    kegg_rxn_dict = {} # Default
    if kegg:
        kegg_comp_dict, kegg_rxn_dict = get_raw_KEGG()

    # Acquire raw MINE dictionaries
    start_ids = [] # Default
    mine_comp_dict = {} # Default
    mine_rxn_dict = {} # Default
    if mine and not enhance:
        start_ids = list(set(KEGG_to_MINE_id(start_kegg_ids).values()))
        raw_mine = get_raw_MINE(start_ids, step_limit, comp_limit, C_limit)
        mine_comp_dict = raw_mine[0]
        mine_rxn_dict = raw_mine[1]

    # Combine KEGG and MINE dictionaries
    mine_comp_dict.update(kegg_comp_dict)
    mine_rxn_dict.update(kegg_rxn_dict)

    # Perform KEGG enhancement
    if not mine and kegg and enhance:
        mine_comp_dict, mine_rxn_dict = enhance_KEGG_with_MINE(
            mine_comp_dict, mine_rxn_dict
        )

    # Create the network
    network = construct_network(mine_comp_dict, mine_rxn_dict,
    start_ids, extra_kegg_ids=start_kegg_ids)

    # Integrate KEGG with MINE
    if mine and kegg and not enhance:
        KEGG_MINE_integration(network)

    # Prune the network
    if not enhance:
        s_out("\nPruning network...\n")
        dummy = distance_to_origin(network, 2, -1)
        prune_network(network)
        s_out("\nDone.\n")
        network_pruned = True

    # Prepare dictionaries
    s_out("\nPreparing compound dictionaries...")
    prepare_dictionaries(network)
    s_out(" Done.\n")

    # Save to Pickle
    pickle.dump(network, open(outfile_name, 'wb'))


if __name__ == "__main__":
    # Read arguments from the commandline
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'outfile',
        help='Write network to Python Pickle file.'
    )
    parser.add_argument(
        '-i', '--infile',
        help='Read KEGG compound identifiers from text file.'
    )
    parser.add_argument(
        '-M', '--mine', action='store_true',
        help='Use MINE for network construction.'
    )
    parser.add_argument(
        '-K', '--kegg', action='store_true',
        help='Use KEGG for network construction.'
    )
    parser.add_argument(
        '-r', type=int, default=10,
        help='Maximum number of MINE reaction steps to download.'
    )
    parser.add_argument(
        '-c', type=int, default=400000,
        help='Maximum number of MINE compounds to download.'
    )
    parser.add_argument(
        '-C', type=int, default=25,
        help='Maximum number of C atoms per MINE molecule.'
    )
    parser.add_argument(
        '-e', '--enhance', action='store_true',
        help='Produce a KEGG network enhanced with MINE reactions.'
    )

    args = parser.parse_args()

    main(args.infile, args.mine, args.kegg, args.r, \
    args.c, args.C, args.enhance, args.outfile)
