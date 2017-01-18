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
import pandas as pd

# Import scripts
from progress import Progress
import mineclient3 as mc
from poppy_origin_helpers import *
from poppy_helpers import *
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


def drGs_for_pathway(pathway, dfG_dict):
    drGs = {}
    for reaction in pathway.split("\n"):
        r_id = reaction.split("\t")[0]
        r_eq = reaction.split("\t")[1]
        # Remove compartment tags from compound IDs
        r_eq = re.sub("_\[[a-z]{3}\]", "", r_eq)
        drGs[r_id] = reaction_gibbs(r_eq, dfG_dict)
    return drGs


def pathways_to_mdf(pathways, dfGs, ne_con, eq_con, n_procs=4,
    T=298.15, R=8.31e-3, network_text="", x_max=0.1, x_min=0.0000001):
    """Create a dictionary with pathways and their MDF values"""

    # The cytosol is selected as the default compartment
    # This is where the pathway reactions take place
    # The "_[cyt]" tag will be removed
    #network_text = network_text.replace("_[cyt]", "")

    # Create a stoichiometric matrix of the background network
    S_net = mdf.read_reactions(network_text)

    # Define the worker
    def worker(dfGs, S_net):
        while True:
            pw_int_chunk = work.get()
            if pw_int_chunk is None:
                break
            mdf_results = []
            for pathway, pw_int in [(pathways[n], n) for n in pw_int_chunk]:

                # Construct standard reaction Gibbs energy dataframe
                if network_text:
                    pathway_and_network = pathway.strip() + "\n" + network_text.strip()
                else:
                    pathway_and_network = pathway
                drGs_d = drGs_for_pathway(pathway_and_network, dfGs)
                drGs_t = "\n".join(['{}\t{}'.format(k,v) for k,v in drGs_d.items()])
                drGs = mdf.read_reaction_drGs(drGs_t)

                # Construct stoichiometric matrix
                S_pat = mdf.read_reactions(pathway)
                S = pd.concat([S_pat, S_net], axis={0,1}).fillna(0)

                # Construct c vector
                c = mdf.mdf_c(S)

                # Construct A matrix
                A = mdf.mdf_A(S, list(S_net.columns))

                # Construct b vector
                b = mdf.mdf_b(S, drGs, ne_con, x_max, x_min)

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

                # Add a result to the list
                if mdf_result.success:
                    mdf_results.append((pw_int, mdf_result.x[-1] * T * R))
                else:
                    mdf_results.append((pw_int, None))
            output.extend(mdf_results)
            with lock:
                n_work_done.value += 1

    # Set up reporter thread
    def reporter():
        n_work = len(pathways)
        n_left = len(pathways)
        p = Progress(design='cp', max_val=n_work)
        while n_left > 0:
            # Check progress
            n_done = len(output)
            n_left = n_work - n_done
            s_out("\rPerforming pathway MDF analysis... %s" % p.to_string(n_done))
            time.sleep(1)

    with mp.Manager() as manager:
        # Initialize Work queue in manager
        work = manager.Queue()

        for pw_int_chunk in chunks(list(range(len(pathways))), 200):
            work.put(pw_int_chunk)

        # Place stop signals on queue
        for i in range(n_procs):
            work.put(None)

        n_work = work.qsize()

        # Initialize output list in manager
        output = manager.list()

        # Initialize number of tasks done counter
        n_work_done = mp.Value('i', 0)
        lock = mp.Lock()

        # Start reporter
        reporter = threading.Thread(target=reporter)
        reporter.start()

        # Start processes
        procs = []
        for i in range(n_procs):
            p = mp.Process(target=worker, args=(dfGs, S_net))
            procs.append(p)
            p.start()

        # Wait until all work is done
        while n_work_done.value != n_work - n_procs:
            time.sleep(1)

        # Terminate the processes
        for p in procs:
            if p.is_alive():
                p.terminate()

        # Stop reporter
        reporter.join()
        print("")

        # Construct MDF dictionary
        mdf_dict = {}
        for mdf_result in list(output):
            mdf_dict[pathways[mdf_result[0]]] = mdf_result[1]
        return mdf_dict


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


# Main code block
def main(pathway_file, outfile, dfG_json, pH, ne_con_file, eq_con_file,
         n_procs=1, T=298.15, R=8.31e-3):

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
    mdf_dict = pathways_to_mdf(
        pathways, dfG_dict, ne_con, eq_con, n_procs, T, R
    )

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
        '-p', '--processes', type=int, default=1,
        help='Number of parallel processes to run.'
    )
    parser.add_argument(
        '-T', type=float, default=298.15,
        help='Temperature (K).'
    )
    parser.add_argument(
        '-R', type=float, default=8.31e-3,
        help='Universal gas constant (kJ/(mol*K)).'
    )
    parser.add_argument(
        'pathways', type=str,
        help='Read poppy pathways text file.'
    )
    parser.add_argument(
        'outfile', type=str, default=False,
        help='Save modified pathways text file with MDF values.'
    )
    args = parser.parse_args()
    main(args.pathways, args.outfile, args.gibbs,
         args.pH, args.constraints, args.ratios,
         args.processes, args.T, args.R)
