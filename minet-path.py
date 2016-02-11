#!/usr/bin/env python3

# Import modules
import networkx as nx
import MineClient3 as mc
import sys
import argparse
import pickle
import multiprocessing

# Define functions
def CountRxns(path):

def



def FindPaths(network, start_id, target_id, reaction_limit):
    """Find all simple paths from a starting compound to a target compound,
    limiting the total number of reactions."""
    start_node = ('c', start_id)
    target_node = ('c', target_id)
    paths = []
    try:
        for path in nx.shortest_simple_paths(network, start_node, target_node, weight=weight):
            reaction_number = 0
            for node in path:
                if node[0] != 'c':
                    reaction_number += 1
            if reaction_number > reaction_limit: break
            paths.append(path)
    except nx.NetworkXNoPath:
        pass
    return paths

def test_FindPaths():
    # Test here
    assert FindPaths() != None


def EnumeratePaths(network, start_id, target_id, reaction_limit):
    """Enumerate branched pathways from the starting compound to the target compound,
    limiting the total number of reactions."""
    # Code here
    return None

def test_EnumeratePaths():
    # Test here
    assert EnumeratePaths() != None


# Main code block
def main(infile_name, compound, reaction_limit, n_threads, outfile_name):
    # Code here


if __name__ == "__main__":
    # Read arguments from the commandline
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='Read minet network pickle.')
    parser.add_argument('outfile', help='')
    parser.add_argument('-c', type=str, help='')
    parser.add_argument('-r', type=int, default=10, help='')
    parser.add_argument('-N', type=int, default=1, help='')
    args = parser.parse_args()
    main()
