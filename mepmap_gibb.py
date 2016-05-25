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
from mepmap_origin_helpers import *
from mepmap_helpers import *
from mepmap_KEGG_helpers import *
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

def test_main():

    import tempfile

    # Three KEGG reactions and one manual
    compounds = [
        "C00024", "C00025", "C00010", "C00624",
        "C00399", "C00042", "C00390", "C00122",
        "C00047", "C00007", "C00990", "C00011",
        "C00001", "C00207", "C00002", "C00164",
        "C00020", "C00009"
    ]

    # Three pHs per compound (6.5, 7.0 and 7.5); Columns are cpd, pH, dG0, error
    # Sorted output
    exp_data_dfG = {
        7.0: {'C00009': -1056.2, 'C00624': -455.0, 'C00122': -519.6,
              'C00025': -368.3, 'C00207': 83.1, 'C00002': -2295.8,
              'C00047': 249.7, 'C00399': 712.7, 'C00164': -285.7,
              'C00042': -518.6, 'C00990': 340.1, 'C00024': -1855.8,
              'C00390': 692.1, 'C00010': -1796.5, 'C00011': -386.0,
              'C00007': 16.4, 'C00001': -157.6, 'C00020': -548.9},
        6.5: {'C00011': -386.0, 'C00009': -1060.2, 'C00122': -525.3,
              'C00042': -530.2, 'C00002': -2331.2, 'C00624': -480.8,
              'C00047': 206.9, 'C00164': -300.0, 'C00207': 65.9,
              'C00399': 638.6, 'C00990': 303.0, 'C00024': -1952.9,
              'C00010': -1887.9, 'C00025': -391.2, 'C00390': 612.3,
              'C00007': 16.4, 'C00001': -163.3, 'C00020': -583.4},
        7.5: {'C00009': -1052.9, 'C00624': -429.3, 'C00122': -513.9,
              'C00042': -507.2, 'C00207': 100.2, 'C00002': -2261.1,
              'C00011': -386.0, 'C00047': 292.5, 'C00164': -271.5,
              'C00025': -345.6, 'C00990': 377.2, 'C00399': 786.9,
              'C00024': -1758.8, 'C00010': -1705.3, 'C00390': 772.0,
              'C00007': 16.4, 'C00001': -151.9, 'C00020': -514.5}
    }

    # Write the expected output file
    exp_file_dfG = tempfile.NamedTemporaryFile()
    exp_file_dfG.write(str.encode(json.dumps(exp_data_dfG)))
    exp_file_dfG.flush()

    # Prepare output file
    out_dfG = tempfile.NamedTemporaryFile()
    out_dfG.flush()

    # Run main()
    main(compounds, out_dfG.name, pHs = [6.5,7.0,7.5])

    # Check output
    assert json.load(open(out_dfG.name, 'r')) == json.load(open(exp_file_dfG.name, 'r'))


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

    pHs = [7.0]

    main(compounds, args.outfile, pHs)
