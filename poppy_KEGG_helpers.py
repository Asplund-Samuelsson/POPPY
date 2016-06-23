# Minet KEGG access functions

# Import modules
import re
import sys
import time
import queue
import threading
from requests import get as rget
from rdkit import Chem

# Import scripts
from poppy_helpers import *
from progress import Progress

# Define functions
def get_KEGG_text(kegg_id, krest="http://rest.kegg.jp"):
    """
    Downloads the raw text entry for the provided KEGG compound or reaction ID,
    via the KEGG rest API @ http://rest.kegg.jp/
    """

    # Check ID and construct query
    if re.fullmatch("^R[0-9]{5}$", kegg_id):
        # KEGG reactions
        krest = "/".join([krest,"get","rn:"+kegg_id])
    elif re.fullmatch("^C[0-9]{5}$", kegg_id):
        # KEGG compounds
        krest = "/".join([krest,"get","cpd:"+kegg_id])
    else:
        # Invalid ID
        s_err("Warning: '" + str(kegg_id) + \
        "' is not a valid KEGG reaction or compound ID.\n")
        return None

    n = 0

    while True:
        r = rget(krest)
        if r.status_code == 200:
            return r.text
        else:
            # Server not returning a result, try again
            n += 1
            if n >= 5:
                s_err("Warning: Unable to download KEGG data for '" + \
                str(kegg_id) + "'.\n")
                return None
            time.sleep(2)


def KEGG_rest_dict(kegg_text):
    """
    Parses a KEGG rest text record into a dictionary. Accepts a single record,
    as it stops after the first '///'.
    """

    kegg_dict = {}

    for line in kegg_text.split("\n"):
        if line == "///":
            # The first entry ends here
            break
        if not line.startswith(" "):
            line = line.split()
            key = line[0]
            line = line[1:]
            try:
                kegg_dict[key].extend(line)
            except KeyError:
                kegg_dict[key] = line
        else:
            try:
                kegg_dict[key].extend(line.split())
            except NameError:
                s_err("Warning: KEGG text line '%s' has no key." % line)

    return kegg_dict


def format_KEGG_reaction(kegg_text):
    """Formats a reaction KEGG rest text record in the MINE database format."""

    # Parse the text into a dictionary
    kegg_dict = KEGG_rest_dict(kegg_text)

    # Check that the needed keys are there
    # ENZYME is missing for known non-enzymatic reactions

    # Ensure that there is an ID
    if "ENTRY" not in kegg_dict.keys():
        s_err("\nWarning: The following KEGG reaction record lacks an" + \
        " ENTRY key: '%s'\n" % str(kegg_dict))
        return None

    # Get the ID
    _id = kegg_dict['ENTRY'][0]

    # Check that the record deals with a compound
    if not re.fullmatch("^R[0-9]{5}$", _id):
        s_err("\nWarning: '%s' is not a valid KEGG compound ID.\n" % str(_id))
        return None

    # Ensure that there is an EQUATION
    if "EQUATION" not in kegg_dict.keys():
        s_err("\nWarning: The following KEGG reaction record lacks an" + \
        " EQUATION key: '%s'\n" % str(kegg_dict))
        return None

    # Get the Operators (enzyme ECs)
    try:
        Operators = kegg_dict['ENZYME']
    except KeyError:
        Operators = ['NA']

    # Re-format the equation
    reactants = True
    n = 1

    Reactants = []
    Products = []

    for segment in kegg_dict['EQUATION']:
        if segment == "<=>":
            reactants = False
            n = 1
            continue
        if segment == "+":
            n = 1
            continue
        if re.fullmatch("^[0-9]+$", segment):
            n = int(segment)
            continue
        # If we make it here, we have a compound ID
        if reactants:
            Reactants.append([n,segment])
        else:
            Products.append([n,segment])

    # Parse the reactant pairs
    RPair = {}
    if 'RPAIR' in kegg_dict.keys():
        p_list = []
        for segment in kegg_dict['RPAIR']:
            if segment.startswith('['):
                continue
            p_list.append(segment)
            if len(p_list) == 3:
                RPair[p_list[0]] = (p_list[1], p_list[2])
                p_list = []
        # The p_list should be empty at the end of iteration
        if len(p_list):
            s_err("Warning: Unexpected format of RPair list '" + \
            str(kegg_dict['RPAIR']) + "'.")

    return {
        "_id":_id, "Operators":Operators,
        "Reactants":Reactants, "Products":Products,
        "RPair":RPair
    }


def get_KEGG_mol_smiles(kegg_id, krest="http://rest.kegg.jp"):
    """Downloads a KEGG compound molecule object and converts it to SMILES."""

    if not re.fullmatch("^C[0-9]{5}$", kegg_id):
        s_err("\nWarning: '" + str(kegg_id) + "' is not a valid KEGG compound ID.\n")
        return None

    # Set up the query
    krest = "/".join([krest,"get","cpd:"+kegg_id,"mol"])

    # Contact server (several times if necessary)
    n = 0
    while True:
        r = rget(krest)
        if r.status_code == 200:
            try:
                mol = Chem.MolFromMolBlock(r.text)
                smiles = Chem.MolToSmiles(mol)
                return smiles
            except:
                s_err("\nWarning: SMILES could not be produced for" + \
                " KEGG ID '%s'.\n" % str(kegg_id))
                return None
        else:
            # Server not returning a result, try again
            n += 1
            if n >= 5:
                s_err("\nWarning: Unable to download molecule data for" + \
                " '%s'.\n" % str(kegg_id))
                return None
            time.sleep(2)


def format_KEGG_compound(kegg_text):
    """Formats a compound KEGG rest text record in the MINE database format."""
    kegg_dict = KEGG_rest_dict(kegg_text)

    # Ensure that the kegg_dict is a dictionary
    if not type(kegg_dict) == dict:
        s_err("\nWarning: KEGG text record did not yield a dictionary: '" + \
        str(kegg_text) + "'")
        return None

    compound = {}

    # Ensure that there is an ID
    if "ENTRY" not in kegg_dict.keys():
        s_err("\nWarning: The following KEGG compound record lacks an" + \
        " ENTRY key: '%s'\n" % str(kegg_dict))
        return None

    # Add ID
    compound['_id'] = kegg_dict['ENTRY'][0]

    # Check that the record deals with a compound
    if not re.fullmatch("^C[0-9]{5}$", compound['_id']):
        s_err("\nWarning: '" + str(compound['_id']) + \
        "' is not a valid KEGG compound ID.\n")
        return None

    # Add DB_links
    compound['DB_links'] = {'KEGG':[compound['_id']]}

    # Add SMILES if possible
    smiles = get_KEGG_mol_smiles(compound['_id'])
    if smiles:
        compound['SMILES'] = smiles

    # Add Names if possible
    if 'NAME' in kegg_dict.keys():
        names = []
        name = []
        for segment in kegg_dict['NAME']:
            if segment.endswith(";"):
                name.append(segment.rstrip(";"))
                names.append(" ".join(name))
                name = []
            else:
                name.append(segment)
        if len(name):
            # Catch trailing name
            names.append(" ".join(name))
        compound['Names'] = names

    # Add Formula if possible
    if 'FORMULA' in kegg_dict.keys():
        compound['Formula'] = kegg_dict['FORMULA'][0]

    # Add reactions if possible
    if 'REACTION' in kegg_dict.keys():
        compound['Reactions'] = kegg_dict['REACTION']

    return compound


def get_KEGG_comps(comp_id_list, num_workers=128):
    """
    Threaded implementation of get_KEGG_text and format_KEGG_compound,
    taking a list of KEGG compound ids as input.
    """
    def worker():
        while True:
            comp_id = work.get()
            try:
                if comp_id is None:
                    work.task_done()
                    break
                kegg_text = get_KEGG_text(comp_id)
                if not kegg_text is None:
                    kegg_comp = format_KEGG_compound(kegg_text)
                    output.put(kegg_comp)
                else:
                    output.put(None)
                work.task_done()
            except:
                work.put(comp_id)
                work.task_done()

    work = queue.Queue()
    output = queue.Queue()

    threads = []

    for i in range(num_workers):
        t = threading.Thread(target=worker)
        t.start()
        threads.append(t)

    for comp_id in comp_id_list:
        work.put(comp_id)

    # Report progress
    M = len(comp_id_list)
    p = Progress(design='pbct', max_val=M)
    while M - output.qsize():
        n = output.qsize()
        p.write(n)
        time.sleep(1)
    n = output.qsize()
    p.write(n)
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

    return list(filter(None, comps))


def get_KEGG_rxns(rxn_id_list, num_workers=128):
    """
    Threaded implementation of get_KEGG_text and format_KEGG_reaction,
    taking a list of KEGG reaction ids as input.
    """
    def worker():
        while True:
            rxn_id = work.get()
            try:
                if rxn_id is None:
                    work.task_done()
                    break
                kegg_text = get_KEGG_text(rxn_id)
                if not kegg_text is None:
                    kegg_rxn = format_KEGG_reaction(kegg_text)
                    output.put(kegg_rxn)
                else:
                    output.put(None)
                work.task_done()
            except:
                work.put(rxn_id)
                work.task_done()

    work = queue.Queue()
    output = queue.Queue()

    threads = []

    for i in range(num_workers):
        t = threading.Thread(target=worker)
        t.start()
        threads.append(t)

    for rxn_id in rxn_id_list:
        work.put(rxn_id)

    # Report progress
    M = len(rxn_id_list)
    p = Progress(design='pbct', max_val=M)
    while M - output.qsize():
        n = output.qsize()
        p.write(n)
        time.sleep(1)
    n = output.qsize()
    p.write(n)
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

    return list(filter(None, rxns))
