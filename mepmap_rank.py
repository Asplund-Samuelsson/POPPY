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
from itertools import product

# Import scripts
from progress import Progress
import mineclient3 as mc
from mepmap_origin_helpers import *
from mepmap_helpers import *

# Define functions


# Main code block
def main():
    pass

if __name__ == "__main__":
    # Read arguments from the commandline
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'infile',
        help='Read mepmap network pickle.'
    )
    parser.add_argument(
        '-o', '--outfile', type=str, default=False,
        help='Save identified pathways in pickle.'
    )
    args = parser.parse_args()
    main()
