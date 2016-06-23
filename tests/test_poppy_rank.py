#!/usr/bin/env python3

# Add repository root to the path
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

# Import the script to be tested
from poppy_rank import *

# Define tests
def test_generate_pathway_hash():
    pathways = [
        "\n".join([
            "RM46785\tC00166 + C00025 <=> C02057 + C00026",
            "RM18439\tC02057 <=> C10438 + C00014",
            "RM64848\tC10438 + C00007 <=> C00048 + C00261"
        ]),
        "\n".join([
            "RM18423\tC00079 <=> C10438 + C00014",
            "RM64848\tC10438 + C00007 <=> C00048 + C00261"
        ]),
        "\n".join([
            "RM2803\tC00266 + C00003 <=> C14448 + C00004 + C00080",
            "RM50787\tC14448 + C00084 <=> C19377 + C00007",
            "RM5767\tC19377 + C00004 + C00080 <=> C01412 + C00003",
            "R03545\tC01412 + C00005 + C00080 <=> C06142 + C00006"
        ])
    ]
    exp_hashes = ['2543cc4384', '993c18ec9b', 'd8ed879133']
    assert exp_hashes == [generate_pathway_hash(P) for P in pathways]
    h1 = generate_pathway_hash('A\tC1 <=> C2\nB\tC1 <=> C3\nC\tC2 + C3 <=> C4')
    h2 = generate_pathway_hash('D\tC1 <=> C3\nE\tC1 <=> C2\nF\tC3 + C2 <=> C4')
    assert h1 == h2


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


def test_reaction_gibbs():
    equation = "(S)-Malate <=> Fumarate + H2O"
    dfG_dict = {"(S)-Malate":10, "Fumarate":20, "H2O":-10}
    assert reaction_gibbs(equation, dfG_dict) == 0
    # R00024
    equation = "C01182 + C00011 + C00001 <=> 2 C00197"
    dfG_dict = {
        "C01182":-2124.3,
        "C00011":-386.0,
        "C00001":-157.6,
        "C00197":-1348.1
        }
    assert reaction_gibbs(equation, dfG_dict) == -28.3
    assert reaction_gibbs("A <=> B", {"A":None,"B":1.2}) == None


def test_read_pathways_text():
    p_text = "\n//\n".join([
        "\n".join([
            ">088ac0e920 MDF 5.0 kJ/mol",
            "RM46785\tC00166 + C00025 <=> C02057 + C00026",
            "RM18439\tC02057 <=> C10438 + C00014",
            "RM64848\tC10438 + C00007 <=> C00048 + C00261"
        ]),
        "\n".join([
            "RM18423\tC00079 <=> C10438 + C00014",
            "RM64848\tC10438 + C00007 <=> C00048 + C00261"
        ]),
        "\n".join([
            "RM2803\tC00266 + C00003 <=> C14448 + C00004 + C00080",
            "RM50787\tC14448 + C00084 <=> C19377 + C00007",
            "RM5767\tC19377 + C00004 + C00080 <=> C01412 + C00003",
            "R03545\tC01412 + C00005 + C00080 <=> C06142 + C00006"
        ])
    ]) + "\n//\n"

    pathways = [
        "\n".join([
            "RM46785\tC00166 + C00025 <=> C02057 + C00026",
            "RM18439\tC02057 <=> C10438 + C00014",
            "RM64848\tC10438 + C00007 <=> C00048 + C00261"
        ]),
        "\n".join([
            "RM18423\tC00079 <=> C10438 + C00014",
            "RM64848\tC10438 + C00007 <=> C00048 + C00261"
        ]),
        "\n".join([
            "RM2803\tC00266 + C00003 <=> C14448 + C00004 + C00080",
            "RM50787\tC14448 + C00084 <=> C19377 + C00007",
            "RM5767\tC19377 + C00004 + C00080 <=> C01412 + C00003",
            "R03545\tC01412 + C00005 + C00080 <=> C06142 + C00006"
        ])
    ]

    assert read_pathways_text(p_text) == pathways


def test_load_dfG_dict():
    import tempfile

    pathways = [
        "\n".join([
            "R1\tC00024 + 2 C00025 <=> 2 C00010",
            "R2\tC00010 + C00624 <=> C00399 + C00042",
            "R3\tC00025 + C00390 <=> C00122 + C00047"
        ]),
        "\n".join([
            "R4\tC00007 + C00990 + C00011 <=> C00001 + C00207 + C00002",
            "R5\tC00207 + 2 C00002 + C00164 <=> C00020 + C00009"
        ])

    ]

    # Three pHs per compound (6.5, 7.0 and 7.5); Columns are cpd, pH, dG0, error
    # Sorted output
    exp_data_dfG = {
        6.5: {'C00011': -386.0, 'C00009': -1060.2, 'C00122': -525.3,
              'C00042': -530.2, 'C00002': -2331.2, 'C00624': -480.8,
              'C00047': 206.9, 'C00164': -300.0, 'C00207': 65.9,
              'C00399': 638.6, 'C00990': 303.0, 'C00024': -1952.9,
              'C00010': -1887.9, 'C00025': -391.2, 'C00390': 612.3,
              'C00007': 16.4, 'C00001': -163.3, 'C00020': -583.4}
    }

    # Write the expected output file
    dfG_json = tempfile.NamedTemporaryFile()
    dfG_json.write(str.encode(json.dumps(exp_data_dfG)))
    dfG_json.flush()

    # Check output
    loaded_dfG = load_dfG_dict(pH = 6.5, dfG_json = dfG_json.name)
    downloaded_dfG = load_dfG_dict(pathways, 6.5)

    assert json.load(open(dfG_json.name, 'r'))["6.5"] == exp_data_dfG[6.5]
    assert loaded_dfG == downloaded_dfG

    d1 = load_dfG_dict(pathways, 7.0)

    assert len(d1) == 18


def test_drGs_for_pathway():
    pathway = "\n".join(["R1\tC1 + 2 C2 <=> C3",
                         "R2\tC3 <=> C4",
                         "R3\tC4 + C5 <=> C6 + C7"])
    dfG_dict = {"C1":-1, "C2":2, "C3":-2, "C4":5, "C5":-3, "C6":10, "C7":-2}
    exp_drGs = {"R1":-5.0, "R2":7.0, "R3":6.0}
    assert drGs_for_pathway(pathway, dfG_dict) == exp_drGs


def test_pathways_to_mdf():
    pathways = [
        "\n".join(["R1\tC1 + C2 <=> C3 + C4",
                   "R2\tC3 <=> C5",
                   "R3\tC5 + C6 <=> C7 + C8"]),
        "\n".join(["R4\tX1 + X2 <=> X3 + X4",
                   "R5\tX3 <=> X5",
                   "R6\tX5 + X6 <=> X7 + X8"]),
        "\n".join(["R0\tZ1 <=> Z2"])
    ]

    # Set all dfGs to -1
    dfG_dict = dict(zip(['C' + str(i) for i in range(1,9)], [-1]*8))
    dfG_dict.update(dict(zip(['X' + str(i) for i in range(1,9)], [-1]*8)))
    dfG_dict.update(dict(zip(['Z' + str(i) for i in range(1,3)], [-1]*8)))

    # Load some simple constraints
    ineq_constraints = mdf.read_constraints("C1\t0.0001\t0.001")
    eq_constraints = mdf.read_ratio_constraints("X1\tX3\t1")

    # Perform MDF for the pathways
    pw_mdf_dict = {}

    # Pathway 1
    S = mdf.read_reactions(pathways[0])
    A = mdf.mdf_A(S)
    c = mdf.mdf_c(S)
    drGs_dict = drGs_for_pathway(pathways[0], dfG_dict)
    drGs_text = "\n".join(['{}\t{}'.format(k,v) for k,v in drGs_dict.items()])
    drGs = mdf.read_reaction_drGs(drGs_text)
    b = mdf.mdf_b(S, drGs, ineq_constraints)
    A_eq = None
    b_eq = None
    mdf_result = mdf.mdf(c, A, b, A_eq, b_eq)
    pw_mdf_dict[pathways[0]] = mdf_result.x[-1] * 298.15 * 8.31e-3

    # Pathway 2
    S = mdf.read_reactions(pathways[1])
    A = mdf.mdf_A(S)
    c = mdf.mdf_c(S)
    drGs_dict = drGs_for_pathway(pathways[1], dfG_dict)
    drGs_text = "\n".join(['{}\t{}'.format(k,v) for k,v in drGs_dict.items()])
    drGs = mdf.read_reaction_drGs(drGs_text)
    b = mdf.mdf_b(S, drGs, ineq_constraints)
    A_eq = mdf.mdf_A_eq(S, eq_constraints)
    b_eq = mdf.mdf_b_eq(eq_constraints)
    mdf_result = mdf.mdf(c, A, b, A_eq, b_eq)
    pw_mdf_dict[pathways[1]] = mdf_result.x[-1] * 298.15 * 8.31e-3

    # Pathway 3
    pw_mdf_dict[pathways[2]] = None # This pathway cannot be optimized

    # Additional constraints to make pathway 3 optimization fail
    ne_con = mdf.read_constraints("C1\t0.0001\t0.001\nZ1\t0.1\t0.1\nZ2\t0.2\t0.2")
    eq_con = mdf.read_ratio_constraints("X1\tX3\t1\nZ1\tZ2\t1")

    assert pathways_to_mdf(pathways, dfG_dict, ne_con, eq_con) == pw_mdf_dict


def test_format_output():

    P1 = "\n".join([
            "R1\tC00024 + 2 C00025 <=> 2 C00010",
            "R2\tC00010 + C00624 <=> C00399 + C00042",
            "R3\tC00025 + C00390 <=> C00122 + C00047"
        ])

    P2 = "\n".join([
            "R4\tC00007 + C00990 + C00011 <=> C00001 + C00207 + C00002",
            "R5\tC00207 + 2 C00002 + C00164 <=> C00020 + C00009"
        ])

    P3 = "\n".join([
            "R6\tC1 + C8 + C1 <=> 2 C2 + C11",
            "R7\tC11 + C4 <=> 2 C20",
            "R8\tC3 + C20 <=> C8 + C21",
            "R9\tC1 + C21 <=> C2 + Z80"
    ])

    P4 = "\n".join([
            "U1\t2 C00025 <=> 2 C00010",
            "U2\tC00624 <=> C00399 + C00042"
        ])

    mdf_dict = {
            P1 : 10.0011,
            P2 : None,
            P3 : 5.0,
            P4 : 10.0011
        }

    exp_output = "\n".join([
        ">" + 'e2729fe9a7' + " MDF 10.001 kJ/mol",
        P4,
        "//",
        ">" + '40e7bed38d' + " MDF 10.001 kJ/mol",
        P1,
        "//",
        ">" + '6f6fe001fc' + " MDF 5.000 kJ/mol",
        P3,
        "//",
        ">" + '65457c1467' + " MDF FAILED",
        P2,
        "//"
    ]) + "\n"

    assert format_output(mdf_dict) == exp_output
