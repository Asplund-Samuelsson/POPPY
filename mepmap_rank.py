#!/usr/bin/env python3

# Import modules
import sys
import argparse
import multiprocessing as mp
import time
from datetime import timedelta as delta
import threading
import re
from itertools import product
from decimal import Decimal
import json
import hashlib

# Import scripts
from progress import Progress
import mineclient3 as mc
from mepmap_origin_helpers import *
from mepmap_helpers import *
import mdf
from equilibrator_query import *

# Define functions
def generate_pathway_hash(pathway):
    """Returns the first 10 characters of an MD5 hash hexdigest"""
    # Remove reaction identifiers and sort to reduce pathway to its essence
    essence = []
    for line in pathway.split("\n"):
        reaction = line.split("\t")[1].split(" <=> ")
        reactants = " + ".join(sorted(reaction[0].split(" + ")))
        products = " + ".join(sorted(reaction[1].split(" + ")))
        essence.append(reactants + " <=> " + products)
    pathway = "\n".join(sorted(essence))
    return hashlib.md5(pathway.encode()).hexdigest()[:10]

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


def reaction_gibbs(equation, dfG_dict):
    """Calculate standard Gibbs reaction energy."""
    s = parse_equation(equation)
    if None in [dfG_dict[c[1]] for c in s[0] + s[1]]:
        # Cannot calculate drG when dfG is missing
        return None
    else:
        p = sum([c[0]*Decimal(str(dfG_dict[c[1]])) for c in s[1]])
        r = sum([c[0]*Decimal(str(dfG_dict[c[1]])) for c in s[0]])
        return float(p - r)

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


def read_pathways_text(pathways_text):
    """Read a pathways string and return it as a list"""

    pathways = []

    # Split by the double forward slashes
    for pathway in pathways_text.split("//"):

        # Remove newline at start and end
        pathway = pathway.strip()

        reactions = []

        # Remove lines beginning with ">"
        for line in pathway.split("\n"):
            if line.startswith(">"):
                continue
            else:
                reactions.append(line)

        # Skip empty pathway
        pathways.append("\n".join(reactions))

    return list(filter(None, pathways))

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


def load_dfG_dict(pathways = None, pH = 7.0, dfG_json = None):
    """Loads dfG json if provided, otherwise queries equilibrator"""

    if dfG_json:
        # Load the JSON file containing dfG values
        try:
            return json.load(open(dfG_json, 'r'))[str(pH)]
        except KeyError:
            pass

    # Identify all compounds in the pathways
    compounds = set()
    for pathway in pathways:
        for reaction in pathway.split("\n"):
            for side in parse_equation(reaction.rstrip().split("\t")[1]):
                for element in side:
                    compounds.add(element[1])

    # Set up equilibrator queries
    queries = list(zip(compounds, [pH]*len(compounds)))

    # Perform equilibrator query
    equilibrator_results = threaded_equilibrator_gibbf(queries)

    # Format dfG dictionary
    dfG_dict = {}
    for query in queries:
        dfG_dict[query[0]] = equilibrator_results[query][0]

    # Return finished dictionary
    return dfG_dict

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
    assert load_dfG_dict(pathways, 7.0) == load_dfG_dict(pH=7.0, dfG_json=dfG_json)


def drGs_for_pathway(pathway, dfG_dict):
    drGs = {}
    for reaction in pathway.split("\n"):
        r_id = reaction.split("\t")[0]
        r_eq = reaction.split("\t")[1]
        drGs[r_id] = reaction_gibbs(r_eq, dfG_dict)
    return drGs

def test_drGs_for_pathway():
    pathway = "\n".join(["R1\tC1 + 2 C2 <=> C3",
                         "R2\tC3 <=> C4",
                         "R3\tC4 + C5 <=> C6 + C7"])
    dfG_dict = {"C1":-1, "C2":2, "C3":-2, "C4":5, "C5":-3, "C6":10, "C7":-2}
    exp_drGs = {"R1":-5.0, "R2":7.0, "R3":6.0}
    assert drGs_for_pathway(pathway, dfG_dict) == exp_drGs


def pathways_to_mdf(pathways, dfGs, ne_con, eq_con):
    """Create a dictionary with pathways and their MDF values"""

    mdf_dict = {}

    # Perform MDF analysis for each pathway

    p = Progress(design='cp', max_val=len(pathways))
    n = 0

    for pathway in pathways:

        # Construct standard reaction Gibbs energy dataframe
        drGs_d = drGs_for_pathway(pathway, dfGs)
        drGs_t = "\n".join(['{}\t{}'.format(k,v) for k,v in drGs_d.items()])
        drGs = mdf.read_reaction_drGs(drGs_t)

        # Construct stoichiometric matrix
        S = mdf.read_reactions(pathway)

        # Construct c vector
        c = mdf.mdf_c(S)

        # Construct A matrix
        A = mdf.mdf_A(S)

        # Construct b vector
        b = mdf.mdf_b(S, drGs, ne_con)

        # Filter the ratio constraints to those relevant to S
        eq_con_f = eq_con[eq_con['cpd_id_num'].isin(S.index)]
        eq_con_f = eq_con_f[eq_con_f['cpd_id_den'].isin(S.index)]

        # Construct A_eq matrix and b_eq vector if equality constraints exist
        if not eq_con_f.empty:
            A_eq = mdf.mdf_A_eq(S, eq_con_f)
            b_eq = mdf.mdf_b_eq(eq_con_f)
        else:
            A_eq = None
            b_eq = None

        # Run MDF optimization
        mdf_result = mdf.mdf(c, A, b, A_eq, b_eq)

        # Add a result to the dictionary
        if mdf_result.success:
            mdf_dict[pathway] = mdf_result.x[-1]
        else:
            mdf_dict[pathway] = None

        # Report progress
        n += 1
        s_out("\rPerforming MDF analysis for pathways... %s" % p.to_string(n))

    print("")
    return mdf_dict

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
    pw_mdf_dict[pathways[0]] = mdf_result.x[-1]

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
    pw_mdf_dict[pathways[1]] = mdf_result.x[-1]

    # Pathway 3
    pw_mdf_dict[pathways[2]] = None # This pathway cannot be optimized

    # Additional constraints to make pathway 3 optimization fail
    ne_con = mdf.read_constraints("C1\t0.0001\t0.001\nZ1\t0.1\t0.1\nZ2\t0.2\t0.2")
    eq_con = mdf.read_ratio_constraints("X1\tX3\t1\nZ1\tZ2\t1")

    assert pathways_to_mdf(pathways, dfG_dict, ne_con, eq_con) == pw_mdf_dict


def format_output(mdf_dict):
    """Format pathway output based on the MDF results dictionary"""
    output = []

    # Group pathways by their MDF
    mdf_to_pathway = {}
    for pathway in mdf_dict:
        if mdf_dict[pathway] is None:
            try:
                mdf_to_pathway["NA"].add(pathway)
            except KeyError:
                mdf_to_pathway["NA"] = set([pathway])
        else:
            try:
                mdf_to_pathway[mdf_dict[pathway]].add(pathway)
            except KeyError:
                mdf_to_pathway[mdf_dict[pathway]] = set([pathway])

    # Sort by descending MDF with failures at the end
    for mdf in sorted(mdf_to_pathway, key=lambda x: (x != "NA", x),
        reverse=True):

        # Sort the pathways by step length
        pathways = sorted(
            list(mdf_to_pathway[mdf]), key=lambda x: x.count("\n")
        )

        # Create title for each pathway and add to output
        for pathway in pathways:
            H = generate_pathway_hash(pathway)
            if mdf != "NA":
                title = ">{0} MDF {1:.3f} kJ/mol".format(H, mdf)
            else:
                title = ">{0} MDF FAILED".format(H)
            output.append(title)
            output.append(pathway)
            output.append("//")

    # Join output into one string and return it
    return "\n".join(output) + "\n"


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


# Main code block
def main(pathway_file, outfile, dfG_json, pH, ne_con_file, eq_con_file):

    print("")

    # Load pathways
    pathways = read_pathways_text(open(pathway_file, 'r').read())

    # Load standard formation Gibbs energy dictionary
    dfG_dict = load_dfG_dict(pathways, pH, dfG_json)

    # Load inequality constraints
    if ne_con_file:
        ne_con_text = open(ne_con_file,'r').read()
    else:
        ne_con_text = ""

    ne_con = mdf.read_constraints(ne_con_text)

    # Load equality constraints
    if eq_con_file:
        eq_con_text = open(eq_con_file,'r').read()
    else:
        eq_con_text = ""

    eq_con = mdf.read_ratio_constraints(eq_con_text)

    # Perform MDF
    mdf_dict = pathways_to_mdf(pathways, dfG_dict, ne_con, eq_con)

    # Format output
    output = format_output(mdf_dict)

    # Write to outfile
    with open(outfile, 'w') as output_file:
        output_file.write(output)

    print("")

if __name__ == "__main__":
    # Read arguments from the commandline
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-g', '--gibbs', type=str,
        help='Read transformed Gibbs standard formation energies (JSON).'
    )
    parser.add_argument(
        '--pH', type=float, default=7.0,
        help='Specify the pH for the thermodynamics calculations.'
    )
    parser.add_argument(
        '-c', '--constraints', type=str,
        help='Read metabolite concentration bounds (inequality constraints).'
    )
    parser.add_argument(
        '-r', '--ratios', type=str,
        help='Read metabolite concentration ratios (equality constraints).'
    )
    parser.add_argument(
        'pathways', type=str,
        help='Read mepmap pathways text file.'
    )
    parser.add_argument(
        'outfile', type=str, default=False,
        help='Save modified pathways text file with MDF values.'
    )
    args = parser.parse_args()
    main(args.pathways, args.outfile, args.gibbs,
         args.pH, args.constraints, args.ratios)
