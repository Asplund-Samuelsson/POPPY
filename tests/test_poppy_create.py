#!/usr/bin/env python3

# Add repository root to the path
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

# Import the script to be tested
from poppy_create import *

# Define tests
def test_allow_reaction_listing():
    # C1 is not a cofactor
    cpd = {"_id":"C1", "Reactions":['R1','R2'], "Formula":"C10"}
    rxn = {"_id":"R1", "RPair":{"RP1":("C1_C2","main"),"RP2":("C3_C4","cofac")}}
    assert allow_reaction_listing(cpd, rxn)

    # C1 is not in an RPair
    rxn = {"_id":"R1", "RPair":{"RP1":("C5_C2","main"),"RP2":("C3_C4","cofac")}}
    assert allow_reaction_listing(cpd, rxn)

    # ATP/ADP reaction pair should not be listed (RP00003)
    cpd = {"_id":"C00002", "Reactions":['R1','R2'], 'Formula':'C10'}
    rxn = {"_id":"R1",
           "Reactants":[[1,'C10000'],[1,'C00002']],
           "Products":[[1,'C00008'],[1,'C00100']]}
    assert not allow_reaction_listing(cpd, rxn)

    # A real ATP to AMP reaction
    cpd = format_KEGG_compound(get_KEGG_text("C00002"))
    rxn = format_KEGG_reaction(get_KEGG_text("R00085")) # ATP -> AMP
    assert not allow_reaction_listing(cpd, rxn)

    # Nucleotides might be involved in other reactions
    cpd = format_KEGG_compound(get_KEGG_text("C00008"))
    rxn = format_KEGG_reaction(get_KEGG_text("R00125")) # AppppA -> ADP
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

    # ...as well as CO2 and other inorganic carbon
    cpd = format_KEGG_compound(get_KEGG_text("C00011"))
    rxn = format_KEGG_reaction(get_KEGG_text(cpd['Reactions'][89]))
    assert not allow_reaction_listing(cpd, rxn)

    cpd = format_KEGG_compound(get_KEGG_text("C00237"))
    rxn = format_KEGG_reaction(get_KEGG_text(cpd['Reactions'][5]))
    assert not allow_reaction_listing(cpd, rxn)

    # Keep single-carbon reduced compounds though
    cpd = format_KEGG_compound(get_KEGG_text("C00132")) # Methanol
    rxn = format_KEGG_reaction(get_KEGG_text(cpd['Reactions'][23]))
    assert allow_reaction_listing(cpd, rxn)

    # Compounds without a formula are disallowed
    cpd = format_KEGG_compound(get_KEGG_text("C00139")) # Ox. ferredoxin
    rxn = format_KEGG_reaction(get_KEGG_text("R05742"))
    assert not allow_reaction_listing(cpd, rxn)

    # Nucleotides taking part in a reaction with nucleotides on both sides
    # are not allowed
    nucleotides = [
        'C00002','C00008','C00020', # ATP, ADP, AMP
        'C00075','C00015','C00105', # UTP, UDP, UMP
        'C00044','C00035','C00144', # GTP, GDP, GMP
        'C00063','C00112','C00055'  # CTP, CDP, CMP
    ]
    for n1 in nucleotides:
        for n2 in nucleotides:
            cpd1 = {'_id':n1, 'Formula':'C10'}
            cpd2 = {'_id':n2, 'Formula':'C10'}
            rxn0 = {
            '_id':'R_' + n1 + '_' + n2,
            'Reactants':[[1,n1],[1,'C15000']],
            'Products':[[1,n2],[1,'C15001']]
            }
            rxn1f = {
            '_id':'R_' + n1,
            'Reactants':[[1,n1],[1,'C15000']],
            'Products':[[1,'C15001']]
            }
            rxn1r = {
            '_id':'R_' + n1,
            'Reactants':[[1,'C15000']],
            'Products':[[1,n1],[1,'C15001']]
            }
            rxn2f = {
            '_id':'R_' + n2,
            'Reactants':[[1,n2],[1,'C15000']],
            'Products':[[1,'C15001']]
            }
            rxn2r = {
            '_id':'R_' + n2,
            'Reactants':[[1,'C15000']],
            'Products':[[1,n2],[1,'C15001']]
            }
            assert not allow_reaction_listing(cpd1, rxn0)
            assert not allow_reaction_listing(cpd2, rxn0)
            assert allow_reaction_listing(cpd1, rxn1f)
            assert allow_reaction_listing(cpd1, rxn1r)
            assert allow_reaction_listing(cpd2, rxn2f)
            assert allow_reaction_listing(cpd2, rxn2r)

    # Cofactor molecules partaking as cofactors are disallowed
    cpd1 = format_KEGG_compound(get_KEGG_text("C00004")) # NADH
    cpd2 = format_KEGG_compound(get_KEGG_text("C00003")) # NAD+
    rxn = format_KEGG_reaction(get_KEGG_text("R00605"))
    assert not allow_reaction_listing(cpd1, rxn)
    assert not allow_reaction_listing(cpd2, rxn)

    cpd1 = format_KEGG_compound(get_KEGG_text("C00006")) # NADPH
    cpd2 = format_KEGG_compound(get_KEGG_text("C00005")) # NADP+
    rxn = format_KEGG_reaction(get_KEGG_text("R01452"))
    assert not allow_reaction_listing(cpd1, rxn)
    assert not allow_reaction_listing(cpd2, rxn)

    cpd1 = format_KEGG_compound(get_KEGG_text("C01352")) # FADH2
    cpd2 = format_KEGG_compound(get_KEGG_text("C00016")) # FAD
    rxn = format_KEGG_reaction(get_KEGG_text("R07934"))
    assert not allow_reaction_listing(cpd1, rxn)
    assert not allow_reaction_listing(cpd2, rxn)


def test_sort_KEGG_reactions():
    # C1 is cofactor in one reaction, not in another
    # C2 is inorganic
    # C3 is a reactant in one reaction, product in another
    # C4 is a product of C3, and lists a reaction that doesn't exist
    # C5 lists a reaction in which it is not listed
    # C6 does not list reactions
    kegg_comp_dict = {
        "C00003":{"_id":"C00003","Reactions":["R1","R2"],"Formula":"C10H18O2"},
        "C2":{"_id":"C2","Reactions":["R3","R4"],"Formula":"XeF4"},
        "C3":{"_id":"C3","Reactions":["R5","R6"],"Formula":"C10H12O3"},
        "C4":{"_id":"C4","Reactions":["R5","RX"],"Formula":"C2H5O"},
        "C5":{"_id":"C5","Reactions":["R7"],"Formula":"CH3O"},
        "C6":{"_id":"C6","Formula":"C12"}
    }
    kegg_rxn_dict = {
        "R1":{"_id":"R1","Reactants":[[1,"C00003"]],"Products":[[1,"X1"]]},
        "R2":{"_id":"R2","Reactants":[[1,"C00003"],[1,"C100"]],
              "Products":[[1,"C101"],[2,"C00004"]]},
        "R3":{"_id":"R3","Reactants":[[1,"C2"]],"Products":[[1,"X2"]]},
        "R4":{"_id":"R3","Reactants":[[1,"Z2"]],"Products":[[1,"C2"]]},
        "R5":{"_id":"R5","Reactants":[[1,"C3"],[1,"Z9"]],"Products":[[1,"C4"]]},
        "R6":{"_id":"R6","Reactants":[[1,"C9"]],"Products":[[1,"C8"],[1,"C3"]]},
        "R7":{"_id":"R7","Reactants":[[1,"X4"]],"Products":[[1,"Z4"]]}
    }
    expected_comp_dict = {
        "C00003":{"_id":"C00003","Reactions":["R1","R2"],"Formula":"C10H18O2",
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


def test_quicksearch():

    # Set up connection
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    db = "KEGGexp2"

    assert quicksearch(con, db, 'C00022')[0]['Names'][0] == 'Pyruvate'
    assert quicksearch(con, db, 'random_query') == []


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


def test_extract_comp_reaction_ids():
    C1 = {'_id':'C1', 'Reactant_in':['R1']}
    C2 = {'_id':'C1', 'Reactant_in':['R2'], 'Product_of':['R3']}
    C3 = {'_id':'C1', 'Product_of':['R3', 'R4']}
    C4 = {'_id':'C4'}
    assert extract_comp_reaction_ids(C1) == ['R1']
    assert extract_comp_reaction_ids(C2) == ['R2','R3']
    assert extract_comp_reaction_ids(C3) == ['R3','R4']
    assert extract_comp_reaction_ids(C4) == []


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


def test_operators_identical():
    assert operators_identical('1.1.1.-', '1.1.1.-')
    assert operators_identical('2.42.1.a', '2.42.-1.a')
    assert not operators_identical('1.2.3.-', '1.2.4.-')


def test_extract_ints():
    L1 = ['1','2','23','c']
    L2 = ['4','42','-1','-']
    L3 = ['1','1','1','22']
    assert extract_ints(L1) == [1,2,23]
    assert extract_ints(L2) == [4,42,-1]
    assert extract_ints(L3) == [1,1,1,22]


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
        {'_id':'X4', 'DB_links':{'KEGG':['C40000','C40001','C99999']}},
    ]
    # Listing the KEGG compound IDs - note that C99999 is missing
    kegg_comp_ids = [
        'C10000','C10001','C20000','C30000',
        'C30001','C40000','C40001','C00002',
        'C00003'
    ]
    exp_rxns = [
        {'_id':'R1_0', 'Operators':['1.1.1.a'],
        'Reactants':[[1,'C10000'],[1,'C20000']],
        'Products':[[1,'C30000'],[1,'C40000']]
        },
        {'_id':'R1_1', 'Operators':['1.1.1.a'],
        'Reactants':[[1,'C10000'],[1,'C20000']],
        'Products':[[1,'C30001'],[1,'C40000']]
        },
        {'_id':'R1_2', 'Operators':['1.1.1.a'],
        'Reactants':[[1,'C10001'],[1,'C20000']],
        'Products':[[1,'C30000'],[1,'C40000']]
        },
        {'_id':'R1_3', 'Operators':['1.1.1.a'],
        'Reactants':[[1,'C10001'],[1,'C20000']],
        'Products':[[1,'C30001'],[1,'C40000']]
        },
        {'_id':'R5_0', 'Operators':['1.2.1.a'],
        'Reactants':[[1,'C00003'],[1,'C20000']],
        'Products':[[1,'C30000'],[1,'C00002']]
        },
        {'_id':'R5_1', 'Operators':['1.2.1.a'],
        'Reactants':[[1,'C00003'],[1,'C20000']],
        'Products':[[1,'C30001'],[1,'C00002']]
        },
        {'_id':'R6_0', 'Operators':['1.3.1.a'],
        'Reactants':[[2,'C10000'],[1,'C10000']],
        'Products':[[1,'C40000'],[1,'C20000']]
        },
        {'_id':'R6_1', 'Operators':['1.3.1.a'],
        'Reactants':[[2,'C10000'],[1,'C10000']],
        'Products':[[1,'C40001'],[1,'C20000']]
        },
        {'_id':'R6_2', 'Operators':['1.3.1.a'],
        'Reactants':[[2,'C10001'],[1,'C10000']],
        'Products':[[1,'C40000'],[1,'C20000']]
        },
        {'_id':'R6_3', 'Operators':['1.3.1.a'],
        'Reactants':[[2,'C10001'],[1,'C10000']],
        'Products':[[1,'C40001'],[1,'C20000']]
        },
        {'_id':'R9_0', 'Operators':['1.1.4.a', '3.4.2.-'],
        'Reactants':[[1,'C20000']],
        'Products':[[1,'C40000']]
        }
    ]
    new_rxns = KEGG_rxns_from_MINE_rxns(rxns, comps, kegg_comp_ids)
    assert len(exp_rxns) == len(new_rxns)
    assert exp_rxns == new_rxns


def test_add_MINE_rxns_to_KEGG_comps():
    comps = [
        {'_id':'C00016', 'Formula':'C2H6', 'Reactant_in':['R00001']},
        {'_id':'C00002', 'Formula':'C2H6',
            'Product_of':['R00001'], 'Reactant_in':['R00002']},
        {'_id':'C00003', 'Formula':'C2H6'},
        {'_id':'C01352', 'Formula':'C2H6', 'Product_of':['R00004']},
        {'_id':'C99999', 'Formula':'C2H6', 'Product_of':['R99999']}
    ]
    rxns = [
        {'_id':'R1_0',
        'Reactants':[[1,'C00002'],[1,'C01352']],
        'Products':[[1,'C00003'],[1,'C00016']],
        },
        {'_id':'R2_2',
        'Reactants':[[1,'C00003']],
        'Products':[[1,'C00016']]
        },
        {'_id':'R3_12',
        'Reactants':[[1,'C00016'],[1,'C00003'],[1,'C01352']],
        'Products':[[1,'C00002']],
        }
    ]
    exp_comps = [
        {'_id':'C00016', 'Formula':'C2H6',
            'Reactant_in':['R00001','R3_12'],
            'Product_of':['R2_2']},
        {'_id':'C00002', 'Formula':'C2H6',
            'Product_of':['R00001','R3_12'],
            'Reactant_in':['R00002','R1_0']},
        {'_id':'C00003', 'Formula':'C2H6',
            'Reactant_in':['R2_2','R3_12'],
            'Product_of':['R1_0']},
        {'_id':'C01352', 'Formula':'C2H6',
            'Product_of':['R00004'],
            'Reactant_in':['R3_12']},
        {'_id':'C99999', 'Formula':'C2H6',
            'Product_of':['R99999']}
    ]
    new_comps = add_MINE_rxns_to_KEGG_comps(comps, rxns)
    assert len(exp_comps) == len(new_comps)
    assert exp_comps == new_comps


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
    {'Products': [[1, 'C01412'], [1, 'C00004'], [1, 'C00080']],
    'Reactants': [[1, 'C06142'], [1, 'C00003']],
    'Operators': ['1.1.1.a'],
    '_id': 'R86e63dd5a75dedff25511e9535e77e2316e4c7af_0'
    },
    {'Products': [[1, 'C06142'], [1, 'C00011']],
    'Reactants': [[1, 'C02804']],
    'Operators': ['4.1.1.k'],
    '_id': 'R766e3bcdbdf841f470c146b7ad9b74ca35d5c3e6_0'}
    ]

    # The expected KEGG reaction dictionary
    exp_KEGG_rxn_dict = deepcopy(KEGG_rxn_dict)
    # It is expected that the second MINE reaction will be included under a new
    # ID, and the first MINE reaction will be merged with the already existing
    # KEGG reaction

    # Addition of the second MINE reaction
    exp_KEGG_rxn_dict['RM1'] = MINE_rxns[1]
    exp_KEGG_rxn_dict['RM1']['MINE_id'] = [MINE_rxns[1]['_id']]
    exp_KEGG_rxn_dict['RM1']['Operators'] = sorted(['M:' + o for o in MINE_rxns[1]['Operators']])
    exp_KEGG_rxn_dict['RM1']['_id'] = 'RM1'

    # Merger of the first MINE reaction
    exp_KEGG_rxn_dict['R03544']['Operators'].extend(['M:1.1.-1.a','M:1.1.1.a'])
    exp_KEGG_rxn_dict['R03544']['Operators'] = sorted(exp_KEGG_rxn_dict['R03544']['Operators'])

    # The expected KEGG compound dictionary
    exp_KEGG_comp_dict = deepcopy(KEGG_comp_dict)

    exp_KEGG_comp_dict['C06142']['Product_of'].append('RM1')
    exp_KEGG_comp_dict['C02804']['Reactant_in'] = ['RM1']

    # Produce enhanced KEGG dictionaries
    enh_KEGG_comp_dict, enh_KEGG_rxn_dict = enhance_KEGG_with_MINE(
        KEGG_comp_dict, KEGG_rxn_dict
    )

    # Assert equality
    assert exp_KEGG_comp_dict == enh_KEGG_comp_dict
    assert exp_KEGG_rxn_dict == enh_KEGG_rxn_dict


def test_KEGG_rxns_Equilibrator_filter():
    rxns = {
    'R1':{'_id':'R1', 'Reactants':[[1,'C00011']], 'Products':[[1,'C00469']]},
    'R2':{'_id':'R2', 'Reactants':[[1,'C00011']], 'Products':[[1,'C00229']]},
    'R3':{'_id':'R3', 'Reactants':[[1,'C01328']], 'Products':[[1,'C06142']]}
    }
    exp_rxns = {
    'R1':{'_id':'R1', 'Reactants':[[1,'C00011']], 'Products':[[1,'C00469']]}
    }
    KEGG_rxns_Equilibrator_filter(rxns)
    assert rxns == exp_rxns


def test_merge_MINE_KEGG_rxns():
    MINE_rxns = [
        {'_id':'R25c7d', 'Operators':['1.1.1.a'],
               'Products':[[1,'C10380'],[1,'C00183']],
               'Reactants':[[1,'C01212'],[1,'C00012']]},
        {'_id':'Rfa4df', 'Operators':['1.5.-1.-'],
               'Reactants':[[1,'C00120'],[1,'C00380']],
               'Products':[[1,'C01212'],[1,'C00012']]},
        {'_id':'R198a2', 'Operators':['5.4.3.c'],
               'Products':[[1,'C10380'],[1,'C00183']],
               'Reactants':[[1,'C00012'],[1,'C01212']]},
        {'_id':'R1faa3', 'Operators':['1.2.3.b','2.1.-3.a'],
               'Products':[[1,'C02010'],[2,'C00110']],
               'Reactants':[[1,'C10000']]},
        {'_id':'R410cc', 'Operators':['1.2.-3.b','2.1.3.a'],
               'Products':[[1,'C10000']],
               'Reactants':[[1,'C02010'],[2,'C00110']]},
        {'_id':'R198a2', 'Operators':['2.4.3.c'],
               'Products':[[1,'C10380']],
               'Reactants':[[1,'C00012']]},
        {'_id':'R1', # Reverse of R2
                'Operators':['2.7.-1.a'],
                'Products':[[1,'C1'],[1,'X1']],
                'Reactants':[[1,'X2'],[1,'C2']]},
        {'_id':'R2',
                'Operators':['2.7.1.a'],
                'Reactants':[[1,'C1'],[1,'X1']],
                'Products':[[1,'C2'],[1,'X2']]},
        {'_id':'R6',
                'Operators':['1.1.1.a', '2.2.2.-'],
                'Reactants':[[1,'C4']],
                'Products':[[2,'C5'],[1,'X3']]},
        {'_id':'R7', # Reverse of R6
                'Operators':['2.2.-2.-', '1.1.-1.a'],
                'Products':[[1,'C4']],
                'Reactants':[[2,'C5'],[1,'X3']]},
        {'_id':'R8',
                'Operators':['1.1.1.a', '2.2.2.-', '3.3.-3.g'],
                'Reactants':[[1,'C4']],
                'Products':[[2,'C5'],[1,'X4']]},
        {'_id':'R9', # Reverse of R8
                'Operators':['2.2.-2.-', '1.1.-1.a', '3.3.3.g'],
                'Products':[[1,'C4']],
                'Reactants':[[2,'C5'],[1,'X4']]},
        {'_id':'R10',
                'Operators':['1.2.3.a', '4.5.-6.-'],
                'Reactants':[[1,'C4']],
                'Products':[[2,'C9'],[1,'X3']]},
        {'_id':'R11', # Reverse of R10 and second encountered
                'Operators':['1.2.-3.a', '4.5.6.-'],
                'Products':[[1,'C4']],
                'Reactants':[[2,'C9'],[1,'X3']]},
        {'_id':'R12', # Reverse of R10 and third encountered
                'Operators':['1.2.-3.a', '4.5.6.-'],
                'Products':[[1,'C4']],
                'Reactants':[[2,'C9'],[1,'X3']]},
        {'_id':'R13', # Reverse of R14
                'Operators':['2.8.-1.a'],
                'Products':[[1,'C8'],[1,'X8']],
                'Reactants':[[1,'C9'],[1,'X9']]},
        {'_id':'R14',
                'Operators':['2.8.1.a'],
                'Reactants':[[1,'C8'],[1,'X8']],
                'Products':[[1,'C9'],[1,'X9']]},
    ]
    KEGG_rxns = [
        {'_id':'R10000', 'Operators':['1.1.1.2'],
              'Reactants':[[1,'C10380'],[1,'C00183']],
               'Products':[[1,'C01212'],[1,'C00012']]},
        {'_id':'R20000', 'Operators':['3.1.4.1'],
              'Reactants':[[1,'C00380'],[1,'C00120']],
               'Products':[[1,'C00012'],[1,'C01212']]},
        {'_id':'R30000', 'Operators':['2.2.2.2'],
              'Reactants':[[1,'C00100']],
               'Products':[[1,'C00200']]}
    ]
    exp_KR = [
        {'_id':'R10000', 'Operators':['1.1.1.2','M:1.1.1.a','M:5.4.3.c'],
              'Reactants':[[1,'C10380'],[1,'C00183']],
               'Products':[[1,'C01212'],[1,'C00012']]},
        {'_id':'R20000', 'Operators':['3.1.4.1','M:1.5.-1.-'],
              'Reactants':[[1,'C00380'],[1,'C00120']],
               'Products':[[1,'C00012'],[1,'C01212']]},
        {'_id':'R30000', 'Operators':['2.2.2.2'],
              'Reactants':[[1,'C00100']],
               'Products':[[1,'C00200']]}
    ]
    exp_MR = [
        {'_id':'RM1',
              'Operators':['M:1.2.-3.b','M:1.2.3.b','M:2.1.-3.a','M:2.1.3.a'],
              'MINE_id':['R1faa3','R410cc'],
               'Products':[[1,'C02010'],[2,'C00110']],
              'Reactants':[[1,'C10000']]},
        {'_id':'RM2',
              'MINE_id':['R198a2'], 'Operators':['M:2.4.3.c'],
               'Products':[[1,'C10380']],
              'Reactants':[[1,'C00012']]},
        {'_id':'RM3',
              'Operators':['M:2.7.-1.a','M:2.7.1.a'],
              'MINE_id':['R1','R2'],
               'Products':[[1,'C1'],[1,'X1']],
              'Reactants':[[1,'X2'],[1,'C2']]},
        {'_id':'RM4',
              'Operators':['M:1.1.-1.a','M:1.1.1.a','M:2.2.-2.-','M:2.2.2.-'],
              'MINE_id':['R6','R7'],
              'Reactants':[[1,'C4']],
               'Products':[[2,'C5'],[1,'X3']]},
        {'_id':'RM5',
              'Operators':['M:1.1.-1.a', 'M:1.1.1.a', 'M:2.2.-2.-',
                           'M:2.2.2.-', 'M:3.3.-3.g', 'M:3.3.3.g'],
              'MINE_id':['R8','R9'],
              'Reactants':[[1,'C4']],
               'Products':[[2,'C5'],[1,'X4']]},
        {'_id':'RM6',
              'Operators':['M:1.2.-3.a','M:1.2.3.a','M:4.5.-6.-','M:4.5.6.-'],
              'MINE_id':['R10','R11','R12'],
              'Reactants':[[1,'C4']],
               'Products':[[2,'C9'],[1,'X3']]},
        {'_id':'RM7',
              'Operators':['M:2.8.-1.a', 'M:2.8.1.a'],
              'MINE_id':['R13','R14'],
               'Products':[[1,'C8'],[1,'X8']],
              'Reactants':[[1,'C9'],[1,'X9']]}
    ]

    MINE_rxns, KEGG_rxns = merge_MINE_KEGG_rxns(MINE_rxns, KEGG_rxns)
    for KEGG_rxn in KEGG_rxns:
        assert KEGG_rxn in exp_KR
    for MINE_rxn in MINE_rxns:
        assert MINE_rxn in exp_MR
    assert MINE_rxns == exp_MR
    assert KEGG_rxns == exp_KR


def test_formula_to_dict():
    assert formula_to_dict("H5C3O2") == {'C':3.0,'O':2.0}
    assert formula_to_dict("H12C6O5Cl") == {'C':6.0,'O':5.0,'Cl':1.0}
    assert formula_to_dict("HCO3", H=True) == {'H':1.0,'C':1.0,'O':3.0}
    assert formula_to_dict("C25H42N7O17P3S", H=True) == {
        'C':25.0, 'H':42.0, 'N':7.0, 'O':17.0, 'P':3.0, 'S':1.0
    }
    assert formula_to_dict("ClCSH") == {'Cl':1.0,'C':1.0,'S':1.0}
    assert formula_to_dict("C14H18O4(C5H8)n", H=True) == {
        'C':19.0, 'H':26.0, 'O':4.0
    }


def test_is_balanced():
    # Get balanced reactions from KEGG
    balanced_rxns = get_KEGG_rxns([
        'R03544','R04429','R08939','R00006',
        'R00351','R01324','R00709','R08549',
        'R00405','R01082','R00342','R02163'
        ])

    # Get unbalanced reactions from KEGG
    unbalanced_rxns = get_KEGG_rxns([
        'R01725','R08609','R05539','R02129'
        ])

    # Get the compounds from KEGG
    comp_ids = set(['C00002', 'C03736', 'C02739', 'C01328'])
    for rxn in balanced_rxns + unbalanced_rxns:
        comp_ids = comp_ids.union(extract_reaction_comp_ids(rxn))
    comp_dict = dict([(c['_id'], c) for c in get_KEGG_comps(comp_ids)])

    # Check the balance status for all reactions
    for rxn in balanced_rxns:
        assert is_balanced(rxn, comp_dict)
    for rxn in unbalanced_rxns:
        assert not is_balanced(rxn, comp_dict)

    # Check output for reactions with float coefficients
    rxn_1 = {'_id': 'RM7', 'Operators': ['M:3.2.2a'],
             'Products': [[1.0, 'C00002'], [1.0, 'C03736']],
             'Reactants': [[1.0, 'C02739'], [1.0, 'C00001']],
             'MINE_id': ['R2d41']}
    rxn_2 = {'_id': 'RM7', 'Operators': ['M:3.2.2a'],
             'Products': [[1.5, 'C00002'], [1.5, 'C03736']],
             'Reactants': [[1.5, 'C02739'], [1.5, 'C00001']],
             'MINE_id': ['R2d41']}
    rxn_3 = {'_id': 'RM7', 'Operators': ['M:3.2.2a'],
             'Products': [[1, 'C00002'], [1, 'C03736']],
             'Reactants': [[1, 'C02739'], [1, 'C00001']],
             'MINE_id': ['R2d41']}
    rxn_4 = {'_id': 'RM7', 'Operators': ['M:3.2.2a'],
             'Products': [[1.5, 'C00002'], [1, 'C03736']],
             'Reactants': [[1, 'C02739'], [5, 'C00001']],
             'MINE_id': ['R2d41']}
    rxn_5 = {'_id': 'RM7', 'Operators': ['M:3.2.2a'],
             'Products': [[1.0, 'C00002'], [1.0, 'C03736']],
             'Reactants': [[1.0, 'C02739'], [1.0, 'C01328']],
             'MINE_id': ['R2d41']}
    assert  is_balanced(rxn_1, comp_dict)
    assert  is_balanced(rxn_2, comp_dict)
    assert  is_balanced(rxn_3, comp_dict)
    assert  not is_balanced(rxn_4, comp_dict)
    assert  not is_balanced(rxn_5, comp_dict)
