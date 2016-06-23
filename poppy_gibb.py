#!/usr/bin/env python3

# Import modules
import networkx as nx
import sys
import argparse
import pickle
import multiprocessing as mp
import time
from datetime import timedelta as delta
import threading
import re
import json
from itertools import product

# Import scripts
from progress import Progress
import mineclient3 as mc
from poppy_origin_helpers import *
from poppy_helpers import *
from poppy_KEGG_helpers import *
from equilibrator_query import equilibrator_gibbf, threaded_equilibrator_gibbf

# Define functions


# Main code block
def main(compounds, out_dfG, pHs, ionic_strength=0.1):
    # Query Equilibrator
    equilibrator_queries = []
    for pH in pHs:
        for compound in compounds:
            equilibrator_queries.append((compound, pH, ionic_strength))
    equilibrator_results = set()
    dfG_dicts = {}
    for e in threaded_equilibrator_gibbf(equilibrator_queries).items():
        e = list(e)
        # Extract the dfG value
        if e[1] is None or not e[1]:
            dfG = None
        else:
            dfG = e[1][0]
        # Extract the pH value
        pH = e[0][1]
        # Extract the compound ID
        compound = e[0][0]
        # Add compound dfG to dictionary under each pH
        try:
            dfG_dicts[pH][compound] = dfG
        except KeyError:
            dfG_dicts[pH] = { compound : dfG }
    print("")

    # Write to outfile
    with open(out_dfG, 'w') as out:
        json.dump(dfG_dicts, out, sort_keys = True, indent = 4)


if __name__ == "__main__":
    # Read arguments from the commandline
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'outfile', type=str, default=False,
        help='Write KEGG compound transformed standard formation Gibbs energies to json file.'
    )
    args = parser.parse_args()

    # Download KEGG compound IDs
    s_out("Downloading KEGG compound list...")
    compounds = []
    r = rget("/".join(["http://rest.kegg.jp","list","compound"]))
    if r.status_code == 200:
        for line in r.text.split("\n"):
            if line == "": break # The end
            kegg_comp_id = line.split()[0].split(":")[1]
            compounds.append(kegg_comp_id)
    else:
        msg = "Error: Unable to download KEGG rest compound list.\n"
        sys.exit(msg)
    s_out(" Done.\n")

    pHs = [6.8, 7.4, 8.0]

    main(compounds, args.outfile, pHs)
