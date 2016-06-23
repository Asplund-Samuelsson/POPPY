#!/usr/bin/env python3

# Add repository root to the path
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

# Import the script to be tested
from pykegg import *

# Define tests
def test_kegg_get(capsys):
    kegg1 = "R01393"
    kegg2 = "C00001"
    kegg3 = "C99999"
    kegg4 = "X99999"
    kegg5 = "M00010"
    assert kegg_get(kegg1) == rget("http://rest.kegg.jp/get/rn:R01393").text
    assert kegg_get(kegg2) == rget("http://rest.kegg.jp/get/cpd:C00001").text
    assert kegg_get(kegg3) == None
    assert kegg_get(kegg4) == None
    out, err = capsys.readouterr()
    assert err == "Warning: Unable to download KEGG data for 'C99999'.\nWarning: Unable to download KEGG data for 'X99999'.\n"
    assert kegg_get(kegg5) == rget("http://rest.kegg.jp/get/M00010").text


def test_threaded_kegg_get():
    keggs = ["R01393", "C00001", "C99999", "X99999", "M00010"]
    expected = [kegg_get(k) for k in keggs]
    t_results = threaded_kegg_get(keggs)
    assert expected == [t_results[k] for k in keggs]


def test_create_kegg_dict():
    assert create_kegg_dict(kegg_get("R05735"))['ENZYME'] == [['6.4.1.6']]
    assert create_kegg_dict(kegg_get("C18020"))['FORMULA'] == [['C10H18O2']]
    assert set(create_kegg_dict(kegg_get("C18020")).keys()) == set([
        "ENTRY","NAME","FORMULA",
        "EXACT_MASS","MOL_WEIGHT","REACTION",
        "ENZYME","DBLINKS","ATOM","BOND"
        ])
    assert create_kegg_dict(kegg_get("C00999"))['NAME'] == [["Ferrocytochrome","b5"]]
    assert len(create_kegg_dict(kegg_get("C00999"))['ENZYME']) == 6
    assert len([x for y in create_kegg_dict(kegg_get("C00999"))['ENZYME'] for x in y]) == 21
    assert create_kegg_dict(kegg_get("R06519"))['EQUATION'][0][9] == "<=>"
    assert " ".join(create_kegg_dict(kegg_get("M00009"))['NAME'][0][0:2]) == "Citrate cycle"


def test_kegg_smiles(capsys):
    assert kegg_smiles("C06099") == "C=C(C)C1CC=C(C)CC1"
    assert kegg_smiles("C09908") == "Cc1ccc(C(C)C)c(O)c1"
    assert kegg_smiles("C06142") == "CCCCO"
    assert kegg_smiles("C01412") == "CCCC=O"
    assert kegg_smiles("XYZ") == None
    assert kegg_smiles("C00999") == None
    assert kegg_smiles("C99999") == None

    out, err = capsys.readouterr()
    assert err == "".join([
        "\nWarning: 'XYZ' is not a valid KEGG compound ID.\n",
        "\nWarning: Unable to download molecule data for 'C00999'.\n",
        "\nWarning: Unable to download molecule data for 'C99999'.\n"
    ])


def test_Reaction():
    r1 = Reaction(create_kegg_dict(kegg_get("R05735")))
    r2 = Reaction(create_kegg_dict(kegg_get("R00259")))
    assert r1.id == "R05735"
    assert r1.equation == "C00207 + C00011 + C00002 + 2 C00001 <=> C00164 + C00020 + 2 C00009"
    assert r1.compounds == set(['C00207', 'C00011', 'C00002', 'C00001', 'C00164', 'C00020', 'C00009'])
    assert r1.reactants == set(['C00207', 'C00011', 'C00002', 'C00001'])
    assert r1.products == set(['C00164', 'C00020', 'C00009'])
    assert r2.id == "R00259"
    assert r2.equation == "C00024 + C00025 <=> C00010 + C00624"
    assert r2.compounds == set(['C00024', 'C00025', 'C00010', 'C00624'])
    assert r2.reactants == set(['C00024', 'C00025'])
    assert r2.products == set(['C00010', 'C00624'])


def test_parse_equation():
    eq = "A + B <=> 3 C + D"
    peq = ([(1,"A"),(1,"B")],[(3,"C"),(1,"D")])
    assert parse_equation(eq) == peq
    eq = "1 C00001 + 2 C00002 <=> 3 C00003 + 4 C00004"
    peq = ([(1,"C00001"),(2,"C00002")],[(3,"C00003"),(4,"C00004")])
    assert parse_equation(eq) == peq
    eq = "Comp-A + 2 Comp-B <=> 2 Comp-C"
    peq = ([(1,"Comp-A"),(2,"Comp-B")],[(2,"Comp-C")])
    assert parse_equation(eq) == peq
