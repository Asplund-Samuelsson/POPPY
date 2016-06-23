#!/usr/bin/env python3

# Minet KEGG access functions

# Import modules
import re
import sys
import time
import queue
import threading
from time import sleep
from requests import get as rget
from rdkit import Chem

# Define functions
def sWrite(string):
    sys.stdout.write(string)
    sys.stdout.flush()


def sError(string):
    sys.stderr.write(string)
    sys.stderr.flush()


def kegg_get(kegg_id, server="http://rest.kegg.jp"):
    """
    Downloads the raw REST text entry for the provided KEGG ID,
    via the KEGG rest API @ http://rest.kegg.jp/
    """
    # Construct query
    rest_address = "/".join([server,"get",kegg_id])
    # Set up con_attempts counter
    con_attempts = 0
    # Try connecting
    while con_attempts < 5:
        con_attempts += 1
        r = rget(rest_address)
        if r.status_code == 200:
            return r.text
        else:
            # Server not returning a result, try again
            time.sleep(2)
    # The connection attempt limit was reached
    sError("Warning: Unable to download KEGG data for '%s'.\n" % str(kegg_id))
    return None


def threaded_kegg_get(queries):
    """Threaded implementation of kegg_get."""

    def worker():
        while True:
            query = work.get()
            if query is None:
                break
            result = kegg_get(query)
            output.put((query, result))
            work.task_done()

    # Initialise queues
    work = queue.Queue()
    output = queue.Queue()

    for query in queries:
        work.put(query)

    # Start threads
    threads = []
    for i in range(16):
        t = threading.Thread(target=worker)
        t.start()
        threads.append(t)

    # Report on progress
    while True:
        if len(queries) == 0:
            progress = 100.0
        else:
            progress = float(output.qsize() / len(queries) * 100)
        sys.stdout.write("\rQuerying KEGG... %0.1f%%" % progress)
        sys.stdout.flush()
        if output.qsize() == len(queries):
            print("")
            break
        sleep(0.5)

    # Join work queue
    work.join()

    # Stop workers
    for i in range(len(threads)):
        work.put(None)
    for t in threads:
        t.join()

    # Get results
    results = {}
    while not output.empty():
        result = output.get()
        results[result[0]] = result[1]
    return results


def create_kegg_dict(kegg_text):
    """
    Parses a KEGG REST text record into a dictionary. Accepts a single record,
    as it stops after the first '///'. All lines under a key are presented as a
    list of lists.
    """
    # Initialize dictionary
    kegg_dict = {}
    # Iterate over lines in raw text
    for line in kegg_text.split("\n"):
        if line == "///":
            # The first entry ends here
            break
        line = re.split(" +", line.rstrip())
        if line[0] != "":
            key = line[0]
            line = line[1:]
            try:
                kegg_dict[key].append(line)
            except KeyError:
                kegg_dict[key] = [line]
        else:
            kegg_dict[key].append(line[1:])
    return kegg_dict


def kegg_smiles(kegg_id, server="http://rest.kegg.jp"):
    """Downloads a KEGG compound molecule object and converts it to SMILES."""
    # Ensure that the ID is a valid KEGG compound
    if not re.fullmatch("^C[0-9]{5}$", kegg_id):
        sError("\nWarning: '%s' is not a valid KEGG compound ID.\n" % str(kegg_id))
        return None
    # Set up the query
    rest_address = "/".join([server,"get","cpd:"+kegg_id,"mol"])
    # Set up con_attempts counter
    con_attempts = 0
    # Try connecting
    while con_attempts < 5:
        con_attempts += 1
        r = rget(rest_address)
        if r.status_code == 200:
            try:
                mol = Chem.MolFromMolBlock(r.text)
                smiles = Chem.MolToSmiles(mol)
                return smiles
            except:
                sError("\nWarning: SMILES could not be produced for KEGG ID '%s'.\n" % str(kegg_id))
                return None
        else:
            time.sleep(2)
    # The connection attempt limit was reached
    sError("\nWarning: Unable to download molecule data for '%s'.\n" % str(kegg_id))
    return None


class Reaction():
    """
    Basic KEGG reaction object constructed from a reaction dictionary.

    equation        The KEGG equation string
    compound_ids    A set of KEGG IDs representing reactants and products
    reactants       A set of KEGG IDs representing reactants only
    products        A set of KEGG IDs representing products only
    """
    def __init__(self, kegg_dict):
        # Reaction ID
        self.id = kegg_dict["ENTRY"][0][0]
        # Get the equation list
        eq_list = kegg_dict["EQUATION"][0]
        # Original equation string
        self.equation = " ".join(eq_list)
        # Create sets of compounds, reactants and products
        is_comp = re.compile("^C[0-9]{5}")
        self.compounds = set(filter(is_comp.match, eq_list))
        self.reactants = set(filter(is_comp.match, eq_list[0:eq_list.index("<=>")]))
        self.products = set(filter(is_comp.match, eq_list[eq_list.index("<=>")+1:]))


def parse_equation(equation):
    """Parses a KEGG reaction string into a tuple of lists of tuples."""
    # Split the equation into lists
    eq_list = re.split(" +", equation)
    reactants = eq_list[0:eq_list.index("<=>")]
    products = eq_list[eq_list.index("<=>")+1:]
    # Set up parser
    def stoichiometry_parse(stoichiometry_list):
        output = []
        s = 1
        for segment in stoichiometry_list:
            if re.match("^[0-9]+$", segment):
                s = int(segment)
                continue
            if segment == "+":
                continue
            else:
                output.append((s,segment))
                s = 1
        return output
    return (stoichiometry_parse(reactants), stoichiometry_parse(products))
