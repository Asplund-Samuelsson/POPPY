# Minet KEGG access functions

# Add repository root to the path
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

# Import the script to be tested
from poppy_KEGG_helpers import *

# Define functions
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
