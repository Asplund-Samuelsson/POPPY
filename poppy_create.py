#!/usr/bin/env python3

# Import modules
import networkx as nx
import re
import sys
import os
import argparse
import pickle
import time
import queue
import threading
import multiprocessing as mp
from requests import get as rget
from rdkit import Chem
from copy import deepcopy
from itertools import repeat, product, combinations

# Import scripts
import mineclient3 as mc
from poppy_helpers import *
from poppy_origin_helpers import *
from poppy_KEGG_helpers import *
from progress import Progress
from equilibrator_query import *

# Specify path to repository
global repo_dir
repo_dir = os.path.dirname(__file__)

# Define functions
def allow_reaction_listing(kegg_comp, kegg_rxn):
    """Determine whether to allow the compound to list the reaction"""

    # Is the compound inorganic?
    inorganic_C = ["C00011","C00288","C01353","C00237"]
    if 'Formula' in kegg_comp.keys():
        if not limit_carbon(kegg_comp, 0) or kegg_comp['_id'] in inorganic_C:
            return False
    else:
        # Does the compound have a formula?
        return False

    # Is the compound CoA, ACP or SAM?
    if kegg_comp['_id'] in {"C00010", "C00229", "C00019"}:
        return False

    # Is the compound a cofactor?
    cofactor_pairs = {
        "C00004" : ("C00004", "C00003"), # NADH/NAD+
        "C00003" : ("C00004", "C00003"), # NADH/NAD+
        "C00005" : ("C00005", "C00006"), # NADPH/NADP+
        "C00006" : ("C00005", "C00006"), # NADPH/NADP+
        "C01352" : ("C01352", "C00016"), # FADH2/FAD
        "C00016" : ("C01352", "C00016"), # FADH2/FAD
        "C00113" : ("C00113", "C01359"), # PQQ/PQQH2
        "C01359" : ("C00113", "C01359"), # PQQ/PQQH2
        "C00019" : ("C00019", "C00021"), # SAM/S-Adenosyl-L-homocysteine
        "C00021" : ("C00019", "C00021"), # SAM/S-Adenosyl-L-homocysteine
        "C00342" : ("C00342", "C00343"), # Thioredoxin/Thiredoxin disulfide
        "C00343" : ("C00342", "C00343"), # Thioredoxin/Thiredoxin disulfide
        "C00399" : ("C00399", "C00390"), # Ubiquinone/Ubiquinol
        "C00390" : ("C00399", "C00390"), # Ubiquinone/Ubiquinol
        "C00138" : ("C00138", "C00139"), # Ferredoxin (red/ox)
        "C00139" : ("C00138", "C00139")  # Ferredoxin (red/ox)
    }
    try:
        C = cofactor_pairs[kegg_comp['_id']]
        R = [c[1] for c in kegg_rxn['Reactants']]
        P = [c[1] for c in kegg_rxn['Products']]
        if C[0] in R and C[1] in P:
            return False
        if C[0] in P and C[1] in R:
            return False
    except KeyError:
        pass

    # Is the compound a nucleotide involved in a reaction with nucleotides on
    # both reactant and product sides?
    nucleotides = [
        'C00002','C00008','C00020', # ATP, ADP, AMP
        'C00075','C00015','C00105', # UTP, UDP, UMP
        'C00044','C00035','C00144', # GTP, GDP, GMP
        'C00063','C00112','C00055', # CTP, CDP, CMP
        'C00054','C03850','C18344'  # AXP variants
    ]
    if kegg_comp['_id'] in nucleotides:
        R = [c[1] for c in kegg_rxn['Reactants']]
        P = [c[1] for c in kegg_rxn['Products']]
        r_nucs = [k for k in R if k in nucleotides]
        p_nucs = [k for k in P if k in nucleotides]
        if r_nucs and p_nucs:
            return False

    # If the compound passed all the tests, the reaction is free to be listed
    return True


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
    # Fix missing KEGG IDs
    mid_to_kegg = {
        "X71306b6c4efe11bc7c485fbc71932f3deb14fa2c":"C00080",
        "X08a914cde05039694ef0194d9ee79ff9a79dde33":"C00001"
    }
    if comp_id in mid_to_kegg.keys():
        kegg_id = mid_to_kegg[comp_id]
        try:
            results['DB_links']['KEGG'] = [kegg_id]
        except KeyError:
            results['DB_links'] = {'KEGG':[kegg_id]}
    # Return compound record
    return results


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


def KEGG_to_MINE_id(kegg_ids):
    """Translate KEGG IDs to MINE IDs."""
    s_out("\nTranslating from KEGG IDs to MINE IDs...\n")
    server_url = "http://modelseed.org/services/mine-database"
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
        s_err("Warning: '%s' has no formula and passes the C-limit.\n" % (comp_id))
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


def get_raw_MINE(comp_id_list, step_limit=10, comp_limit=100000, C_limit=25):
    """Download connected reactions and compounds up to the limits."""

    # Set up connection
    server_url = "http://modelseed.org/services/mine-database"
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


def is_connected_MINE_comp(comp_id, network):
    """Determines if the MINE compound is connected."""
    try:
        pred = network.successors(network.graph['cmid2node'][comp_id])
        succ = network.predecessors(network.graph['cmid2node'][comp_id])
    except KeyError:
        return False
    if len(pred) + len(succ) > 0:
        return True
    else:
        return False


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


def operators_identical(op1, op2):
    """Checks if two BNICE EC operators are eachother's reverse/identical."""
    for e1, e2 in zip(op1.split('.'), op2.split('.')):
        if e1.lstrip('-') != e2.lstrip('-'):
            return False
    return True


def extract_ints(str_list):
    int_list = []
    for s in str_list:
        try:
            int_list.append(int(s))
        except ValueError:
            continue
    return int_list


def remove_redundant_MINE_rxns(rxns):
    """Removes identical but reversed MINE reactions."""

    discarded_rxns = set()

    # Identify all operators
    all_operators = set()
    p = Progress(max_val=len(rxns), design='p')
    n = 0
    for rxn in rxns:
        n += 1
        s_out("\rIdentifying operators... %s" % p.to_string(n))
        for operator in rxn['Operators']:
            all_operators.add(operator)
    s_out("\rIdentifying operators... Done. \n")

    # Identify operators that have a reverse
    operators_with_reverse = set()
    for op1 in all_operators:
        for op2 in all_operators - set([op1]):
            if operators_identical(op1, op2):
                operators_with_reverse.add(op1)

    # Reduce the reactions to those in which all operators have a reverse
    rxns_red = []
    p = Progress(max_val=len(rxns), design='p')
    n = 0
    for rxn in enumerate(rxns):
        n += 1
        s_out("\rIdentifying redundancy candidates... %s" % p.to_string(n))
        add_rxn = True
        for operator in rxn[1]['Operators']:
            if operator not in operators_with_reverse:
                add_rxn = False
                break
        if add_rxn:
            rxns_red.append(rxn)
    s_out("\rIdentifying possibly redundant reactions... Done. \n")

    # Set up progress indicator
    p = Progress(max_val = len(rxns_red)**2, design='pt')
    n = 0

    # Iterate over all reaction pairs
    for rp in product(rxns_red, rxns_red):
        # Report progress
        n += 1
        s_out("\rRemoving redundant MINE reactions... %s" % p.to_string(n))

        # Don't compare a reaction to itself
        if rp[0][0] == rp[1][0]:
            continue

        # Don't perform further comparisons for discarded reactions
        if rp[0][0] in discarded_rxns or rp[1][0] in discarded_rxns:
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


def KEGG_rxns_from_MINE_rxns(rxns, comps, KEGG_comp_ids):
    """Produce reactions with KEGG IDs instead of MINE IDs."""

    # Set up a MINE ID to KEGG ID dictionary - filter KEGG IDs
    K = set(KEGG_comp_ids)
    M2K = dict(
    [(c['_id'], [k for k in c['DB_links']['KEGG'] if k in K]) for c in comps]
    )

    # Produce KEGG-labeled reactions
    KEGG_rxns = []
    p = Progress(max_val=len(rxns), design='p')
    n = 0
    for rxn in rxns:
        # Report progress
        n += 1
        s_out("\rProducing MINE reactions with KEGG IDs... %s" % p.to_string(n))

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

            # Finally, add the reaction to the list of new reactions
            KEGG_rxns.append(new_rxn)

            # ...and count one step up
            n += 1

    s_out("\rProducing MINE reactions with KEGG IDs... Done. \n")

    return KEGG_rxns


def add_MINE_rxns_to_KEGG_comps(comps, rxns):
    """Add MINE reactions to KEGG compound 'Reactant_in'/'Product_of' lists."""

    # Create a KEGG ID to reaction dictionary
    K2R = {}
    p = Progress(max_val=len(rxns), design='p')
    n = 0
    # Go through the rxns and store a set of valid indices for each cpd
    for rxn in enumerate(rxns):
        n += 1
        s_out("\rListing MINE reactions per KEGG compound... %s" % p.to_string(n))
        for comp_id in extract_reaction_comp_ids(rxn[1]):
            try:
                K2R[comp_id].add(rxn[0])
            except KeyError:
                K2R[comp_id] = set([rxn[0]])

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

        # Now go through only the reactions in which the compound takes part
        try:
            comp_rxns = [rxns[i] for i in sorted(list(K2R[comp['_id']]))]
        except KeyError:
            comp_rxns = []

        # Check if the compound should list any of the reactions
        for rxn in comp_rxns:

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

    print("")

    return new_comps


def enhance_KEGG_with_MINE(KEGG_comp_dict, KEGG_rxn_dict):
    """Enhance a raw KEGG reaction network with MINE reactions."""

    print("\nEnhancing KEGG reaction network data with MINE reactions...\n")

    # Set up MINE connection
    server_url = "http://modelseed.org/services/mine-database"
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
    #MINE_rxns = pickle.load(open('/ssd/common/db/mine/MINE_rxns.pickle', 'rb'))

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

    # Create a list of MINE reactions in terms of KEGG compounds
    MINE_rxns = KEGG_rxns_from_MINE_rxns(MINE_rxns, MINE_comps, KEGG_comp_ids)

    # Create KEGG reaction list and remove redundant MINE reactions
    KEGG_rxns = list(KEGG_rxn_dict.values())
    MINE_rxns, KEGG_rxns = merge_MINE_KEGG_rxns(MINE_rxns, KEGG_rxns)

    # Construct a new KEGG reaction dictionary from the KEGG reaction list
    KEGG_rxn_dict = dict(zip([rxn['_id'] for rxn in KEGG_rxns], KEGG_rxns))

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


def KEGG_rxns_Equilibrator_filter(rxns):
    """Remove Equilibrator-incompatible reactions"""

    print("\nChecking Equilibrator compatibility...")

    # Extract the compound IDs
    comp_ids = set()
    for rxn in rxns.values():
        for comp_id in extract_reaction_comp_ids(rxn):
            comp_ids.add(comp_id)

    # Query Equilibrator
    #eq_results = threaded_equilibrator_gibbf([(cid,) for cid in comp_ids])

    # Determine valid compounds
    #valid_comp_ids = set()
    #for query in eq_results.keys():
    #    if eq_results[query]:
    #        valid_comp_ids.add(query[0])

    valid_file = os.path.join(repo_dir, 'data/eq_KEGG_IDs.pkl')
    valid_comp_ids = pickle.load(open(valid_file,'rb'))

    # Identify invalid reactions
    invalid_reactions = set()
    p = Progress(max_val=len(rxns))
    n = 0
    for rxn in rxns.values():
        n += 1
        s_out("\rIdentifying Equilibrator-incompatible reactions... " + \
        p.to_string(n))
        for comp_id in extract_reaction_comp_ids(rxn):
            if comp_id not in valid_comp_ids:
                invalid_reactions.add(rxn['_id'])
                break

    # Remove invalid reactions
    for rxn_id in invalid_reactions:
        del rxns[rxn_id]

    s_out('\rIdentifying Equilibrator-incompatible reactions... Done. \n')


def merge_MINE_KEGG_rxns(MINE_rxns, KEGG_rxns):
    """Merge equivalent MINE reactions into KEGG reactions or new reactions"""

    # Define reaction comparison function
    def stoichiometric_set(rxn):
        try:
            p = frozenset([tuple(c) for c in rxn['Products']])
        except KeyError:
            p = frozenset()
        try:
            r = frozenset([tuple(c) for c in rxn['Reactants']])
        except KeyError:
            r = frozenset()
        return frozenset([r, p])

    # First merge MINE reactions into KEGG reactions when possible
    KEGG_indices = range(len(KEGG_rxns))
    MINE_indices = range(len(MINE_rxns))

    # Set up progress reporting
    N = len(KEGG_indices)*2 + len(MINE_indices)*2
    p = Progress(max_val = N)
    n = 0

    # Group reactions by identical stoichiometry
    S2i = {}

    # KEGG reactions
    for iK in KEGG_indices:
        n += 1
        s_out("\rMerging reactions... %s" % p.to_string(n))
        try:
            set_dict = S2i[stoichiometric_set(KEGG_rxns[iK])]
        except KeyError:
            set_dict = {}
        try:
            set_dict['K'].add(iK)
        except KeyError:
            set_dict['K'] = set([iK])
        S2i[stoichiometric_set(KEGG_rxns[iK])] = set_dict

    # MINE reactions
    for iM in MINE_indices:
        n += 1
        s_out("\rMerging reactions... %s" % p.to_string(n))
        try:
            set_dict = S2i[stoichiometric_set(MINE_rxns[iM])]
        except KeyError:
            set_dict = {}
        try:
            set_dict['M'].add(iM)
        except KeyError:
            set_dict['M'] = set([iM])
        S2i[stoichiometric_set(MINE_rxns[iM])] = set_dict

    # Extend KEGG Operator lists
    for iK in KEGG_indices:
        n += 1
        s_out("\rMerging reactions... %s" % p.to_string(n))
        try:
            iMs = list(S2i[stoichiometric_set(KEGG_rxns[iK])]['M'])
        except KeyError:
            continue
        new_ops = []
        for iM in iMs:
            new_ops.extend(['M:' + op for op in MINE_rxns[iM]['Operators']])
        try:
            new_ops.extend(KEGG_rxns[iK]['Operators'])
        except KeyError:
            pass
        KEGG_rxns[iK]['Operators'] = sorted(list(set(new_ops)))

    # Determine what MINE reactions were merged into KEGG reactions
    MINE_discard = set()
    for iM in MINE_indices:
        n += 1
        s_out("\rMerging reactions... %s" % p.to_string(n))
        try:
            iKs = list(S2i[stoichiometric_set(MINE_rxns[iM])]['K'])
            MINE_discard.add(iM)
        except KeyError:
            pass

    # Now merge MINE reactions
    # Extract equivalent reaction index sets and sort them
    eq_MINE_rxn_sets = []
    for iM in [i for i in MINE_indices if i not in MINE_discard]:
        eq_MINE_rxn_sets.append(
            frozenset(S2i[stoichiometric_set(MINE_rxns[iM])]['M'])
        )

    eq_MINE_rxn_sets = sorted(
        list(set(eq_MINE_rxn_sets)), key = lambda x : min(x)
    )

    # Merge reactions
    new_MINE_rxns = []
    n = 0
    for equivalent_reaction_set in eq_MINE_rxn_sets:

        # The new MINE reaction entry is based on the one with the lowest index
        base_rxn_index = min(equivalent_reaction_set)
        n += 1
        new_MINE_rxn = {
            '_id' : 'RM' + str(n),
            'Operators' : [],
            'MINE_id' : [],
            'Products' : MINE_rxns[base_rxn_index]['Products'],
            'Reactants' : MINE_rxns[base_rxn_index]['Reactants']
        }

        # Add the equivalent MINE reaction IDs and operators
        for i in sorted(list(equivalent_reaction_set)):
            rxn = MINE_rxns[i]
            new_MINE_rxn['MINE_id'].append(rxn['_id'])
            for op in rxn['Operators']:
                new_MINE_rxn['Operators'].append('M:' + op)

        # Make the operator list unique and sorted
        new_MINE_rxn['Operators'] = sorted(list(set(new_MINE_rxn['Operators'])))

        # Append to the list of new MINE reactions
        new_MINE_rxns.append(new_MINE_rxn)

    print("")

    # Return MINE reactions and KEGG reactions
    return (new_MINE_rxns, KEGG_rxns)


def formula_to_dict(formula, H=False):
    """Create a dictionary with atom counts from a formula string

    ARGUMENTS

    formula : string
        A chemical formula.
    H : bool
        Set to True to include H(ydrogen) in the output.

    RETURNS

    formula_dict : dictionary
        A dictionary with the counts of each element in the formula.
    """

    # Setup element-matching regular expressions
    regex = re.compile('[A-Z]{1}[a-z]*[0-9]*')

    # Prepare output dictionary
    formula_dict = {}

    # Find all elements in the formula
    for element in re.findall(regex, formula):
        atom = re.match('[A-Z,a-z]*', element).group()
        try:
            count = next(filter(None, re.findall('[0-9]*', element)))
        except StopIteration:
            count = 1
        if atom != 'H' or H:
            try:
                formula_dict[atom] += float(count)
            except KeyError:
                formula_dict[atom] = float(count)

    return formula_dict


def is_balanced(reaction, comp_dict):
    """Determines whether the reaction is balanced"""
    def increment_elements(elements, n, comp_id):
        try:
            formula = formula_to_dict(comp_dict[comp_id]['Formula'], H=True)
            for element in formula:
                try:
                    elements[element] += float(n) * formula[element]
                except KeyError:
                    elements[element] = float(n) * formula[element]
        except KeyError:
            pass
        return elements

    reactant_elements = {}
    for n, reactant in reaction['Reactants']:
        reactant_elements = increment_elements(reactant_elements, n, reactant)

    product_elements = {}
    for n, product in reaction['Products']:
        product_elements = increment_elements(product_elements, n, product)

    return reactant_elements == product_elements


# Main code block
def main(outfile_name, infile, mine, kegg, step_limit,
    comp_limit, C_limit, enhance, eq_filter):

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
        #kegg_file = '/ssd/common/db/kegg/KEGG_cpd_rxn.pickle'
        #kegg_comp_dict, kegg_rxn_dict = pickle.load(open(kegg_file, 'rb'))
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

    # Filter to equilibrator compatible reactions
    if eq_filter and enhance:
        KEGG_rxns_Equilibrator_filter(mine_rxn_dict)

    # Filter to balanced reactions
    s_out("\nRemoving unbalanced reactions... ")
    unbalanced_rxn_ids = set()
    for rxn_id in mine_rxn_dict:
        if not is_balanced(mine_rxn_dict[rxn_id], mine_comp_dict):
            unbalanced_rxn_ids.add(rxn_id)
    for rxn_id in unbalanced_rxn_ids:
        del mine_rxn_dict[rxn_id]
    s_out(" Done.\n")

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
    parser.add_argument(
        '-E', '--equilibrator_filter', action='store_true',
        help='Remove equilibrator incompatible reactions.'
    )

    args = parser.parse_args()

    main(args.outfile, args.infile, args.mine, args.kegg, args.r, \
    args.c, args.C, args.enhance, args.equilibrator_filter)
