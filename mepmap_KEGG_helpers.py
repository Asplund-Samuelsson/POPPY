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
from mepmap_helpers import *

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

def test_get_KEGG_text(capsys):
    kegg1 = "R01393"
    kegg2 = "C00001"
    kegg3 = "C99999"
    kegg4 = "X99999"
    exp1 = rget("http://rest.kegg.jp/get/rn:R01393").text
    exp2 = rget("http://rest.kegg.jp/get/cpd:C00001").text
    assert get_KEGG_text(kegg1) == exp1
    assert get_KEGG_text(kegg2) == exp2
    assert get_KEGG_text(kegg3) == None
    assert get_KEGG_text(kegg4) == None
    out, err = capsys.readouterr()
    assert err == "\n".join([
        "Warning: Unable to download KEGG data for 'C99999'.",
        "Warning: 'X99999' is not a valid KEGG reaction or compound ID.\n"
    ])


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

def test_KEGG_rest_dict():
    assert KEGG_rest_dict(get_KEGG_text("R05735"))['ENZYME'] == ['6.4.1.6']
    assert KEGG_rest_dict(get_KEGG_text("C18020"))['FORMULA'] == ['C10H18O2']
    assert set(KEGG_rest_dict(get_KEGG_text("C18020")).keys()) == set([
        "ENTRY","NAME","FORMULA",
        "EXACT_MASS","MOL_WEIGHT","REACTION",
        "ENZYME","DBLINKS","ATOM","BOND"
        ])
    names = ["Ferrocytochrome","b5"]
    assert KEGG_rest_dict(get_KEGG_text("C00999"))['NAME'] == names
    assert len(KEGG_rest_dict(get_KEGG_text("C00999"))['ENZYME']) == 21
    assert KEGG_rest_dict(get_KEGG_text("R06519"))['EQUATION'][9] == "<=>"


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

def test_format_KEGG_reaction():
    _id = "R01393"
    Operators = ["4.1.1.40"]
    Reactants = [[1,"C00168"]]
    Products = [[1,"C00266"],[1,"C00011"]]
    RPair = {
        "RP01475":("C00168_C00266","main"),
        "RP06553":("C00011_C00168","leave")
    }
    rxn1 = {
        "_id":_id, "Operators":Operators, "Reactants":Reactants,
        "Products":Products, "RPair":RPair
    }

    _id = "R05735"
    Operators = ["6.4.1.6"]
    Reactants = [[1,"C00207"],[1,"C00011"],[1,"C00002"],[2,"C00001"]]
    Products = [[1,"C00164"],[1,"C00020"],[2,"C00009"]]
    RPair = {
        "RP00010":("C00002_C00009","ligase"),
        "RP00274":("C00164_C00207","main"),
        "RP05676":("C00001_C00009","leave"),
        "RP05804":("C00011_C00164","leave"),
        "RP12346":("C00002_C00020","ligase"),
        "RP12391":("C00002_C00009","ligase")
    }
    rxn2 = {
        "_id":_id, "Operators":Operators, "Reactants":Reactants,
        "Products":Products,"RPair":RPair
    }

    _id = "R06519"
    Operators = ["1.14.19.17"]
    Reactants = [[1,"C12126"],[2,"C00999"],[1,"C00007"],[2,"C00080"]]
    Products = [[1,"C00195"],[2,"C00996"],[2,"C00001"]]
    RPair = {
        "RP00013":("C00001_C00007","cofac"),
        "RP05667":("C00195_C12126","main")
    }
    rxn3 = {
        "_id":_id, "Operators":Operators, "Reactants":Reactants,
        "Products":Products,"RPair":RPair
    }

    _id = "R00178"
    Operators = ["4.1.1.50"]
    Reactants = [[1,"C00019"],[1,"C00080"]]
    Products = [[1,"C01137"],[1,"C00011"]]
    RPair = {
        "RP03935":("C00019_C01137","main"),
        "RP08122":("C00011_C00019","leave")
    }
    rxn4 = {
        "_id":_id, "Operators":Operators, "Reactants":Reactants,
        "Products":Products, "RPair":RPair
    }

    assert format_KEGG_reaction(get_KEGG_text("R01393")) == rxn1
    assert format_KEGG_reaction(get_KEGG_text("R05735")) == rxn2
    assert format_KEGG_reaction(get_KEGG_text("R06519")) == rxn3
    assert format_KEGG_reaction(get_KEGG_text("R00178")) == rxn4


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

def test_get_KEGG_mol_smiles(capsys):
    assert get_KEGG_mol_smiles("C06099") == "C=C(C)C1CC=C(C)CC1"
    assert get_KEGG_mol_smiles("C09908") == "Cc1ccc(C(C)C)c(O)c1"
    assert get_KEGG_mol_smiles("C06142") == "CCCCO"
    assert get_KEGG_mol_smiles("C01412") == "CCCC=O"
    assert get_KEGG_mol_smiles("XYZ") == None
    assert get_KEGG_mol_smiles("C00999") == None
    assert get_KEGG_mol_smiles("C99999") == None

    out, err = capsys.readouterr()
    assert err == "".join([
        "\nWarning: 'XYZ' is not a valid KEGG compound ID.\n",
        "\nWarning: Unable to download molecule data for 'C00999'.\n",
        "\nWarning: Unable to download molecule data for 'C99999'.\n"
    ])

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

def test_format_KEGG_compound():
    comps = [
        {"_id" : "C06142",
        "SMILES" : "CCCCO",
        "Reactions" : ['R03544','R03545'],
        "Names" : ['1-Butanol','n-Butanol'],
        "DB_links" : {'KEGG':['C06142']},
        "Formula" : "C4H10O"},
        {"_id" : "C00999",
        "Reactions" : KEGG_rest_dict(get_KEGG_text("C00999"))['REACTION'],
        "Names" : ["Ferrocytochrome b5"],
        "DB_links" : {'KEGG':['C00999']}},
        {"_id" : "C00006",
        "SMILES" : 'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC' + \
        '3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
        "Reactions" : KEGG_rest_dict(get_KEGG_text("C00006"))['REACTION'],
        "Names" : [
            "NADP+", "NADP", "Nicotinamide adenine dinucleotide phosphate",
            "beta-Nicotinamide adenine dinucleotide phosphate", "TPN",
            "Triphosphopyridine nucleotide", "beta-NADP+"],
        "DB_links" : {'KEGG':['C00006']},
        "Formula" : "C21H29N7O17P3"}
    ]

    for comp in comps:
        assert comp == format_KEGG_compound(get_KEGG_text(comp['_id']))


def get_KEGG_comps(comp_id_list, num_workers=128):
    """
    Threaded implementation of get_KEGG_text and format_KEGG_compound,
    taking a list of KEGG compound ids as input.
    """
    def worker():
        while True:
            comp_id = work.get()
            if comp_id is None:
                work.task_done()
                break
            s_out("\rHandling compound query '%s'." % str(comp_id))
            kegg_text = get_KEGG_text(comp_id)
            if not kegg_text is None:
                kegg_comp = format_KEGG_compound(kegg_text)
                if not kegg_comp is None:
                    output.put(kegg_comp)
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

def test_get_KEGG_comps():
    comp_ids = ["C04625","C13929","C10269","C05119","C02419"]
    comps_1 = [format_KEGG_compound(get_KEGG_text(x)) for x in comp_ids]
    comps_2 = get_KEGG_comps(comp_ids)
    assert len(comps_1) == len(comps_2)
    for comp in comps_1:
        assert comp in comps_2
    for comp in comps_2:
        assert comp in comps_1
    assert get_KEGG_comps(["C99999"]) == []


def get_KEGG_rxns(rxn_id_list, num_workers=128):
    """
    Threaded implementation of get_KEGG_text and format_KEGG_reaction,
    taking a list of KEGG reaction ids as input.
    """
    def worker():
        while True:
            rxn_id = work.get()
            if rxn_id is None:
                work.task_done()
                break
            s_out("\rHandling reaction query '%s'." % str(rxn_id))
            kegg_text = get_KEGG_text(rxn_id)
            if not kegg_text is None:
                kegg_rxn = format_KEGG_reaction(kegg_text)
                if not kegg_rxn is None:
                    output.put(kegg_rxn)
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

def test_get_KEGG_rxns():
    rxn_ids = ["R10430","R07960","R04715","R07211","R10332"]
    rxns_1 = [format_KEGG_reaction(get_KEGG_text(x)) for x in rxn_ids]
    rxns_2 = get_KEGG_rxns(rxn_ids)
    assert len(rxns_1) == len(rxns_2)
    for rxn in rxns_1:
        assert rxn in rxns_2
    for rxn in rxns_2:
        assert rxn in rxns_1
    assert get_KEGG_rxns(["R99999"]) == []
