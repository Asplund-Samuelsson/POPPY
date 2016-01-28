#!/usr/bin/env python3

# Import modules
import networkx as nx
import MineClient3 as mc
import re
import sys

# Define functions
def ReadCompounds(filename):
    compounds = [line.rstrip() for line in open(filename, 'r')]
    for c in compounds:
        if re.fullmatch("^C[0-9]{5}$", c) == None:
            msg = "Warning: The supplied string '", c, "' is not a valid KEGG compound ID."
            sys.exit(msg)
    return compounds

def test_ReadCompounds():

    import pytest
    import tempfile

    t1 = str.encode("C10000\nC40055\nC13482\n")
    t2 = str.encode("C13854\nR10309\nC33190\n")
    t3 = str.encode("\naaaa\ndsofi0309¤¤¤\n309fu80")

    f1 = tempfile.NamedTemporaryFile()
    f1.write(t1)
    f1.flush()
    f2 = tempfile.NamedTemporaryFile()
    f2.write(t2)
    f2.flush()
    f3 = tempfile.NamedTemporaryFile()
    f3.write(t3)
    f3.flush()

    assert set(ReadCompounds(f1.name)) == set(['C10000','C40055','C13482'])
    with pytest.raises(SystemExit): ReadCompounds(f2.name)
    with pytest.raises(SystemExit): ReadCompounds(f3.name)

    f1.close()
    f2.close()
    f3.close()


def KeggToMineId(kegg_ids):
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    return [con.quick_search("KEGGexp2", kegg_id)[0]['_id'] for kegg_id in kegg_ids]

def test_KeggToMineId():
    assert KeggToMineId(['C15667', 'C16519', 'C00130']) == ['C023e725c27a385fd9057c8171ca1992a32a3e9a4',
    'Cfa3c20a209e12ac4a70f10530cf43bea2b4422fe',
    'Cae5be0a0a2bb6b6baef29e4d4c6b3f4c1a67ad19']


def SeedNetwork(nx_graph, compounds):
    # Code here
    return nx_graph

def test_SeedNetwork():
    server_url = "http://bio-data-1.mcs.anl.gov/services/mine-database"
    con = mc.mineDatabaseServices(server_url)
    G = nx.Graph()
    new_compounds = ["Cefbaa83ea06e7c31820f93c1a5535e1378aba42b"]
    assert ExtendNetwork(G, new_compounds) == []


def GrowNetwork(nx_graph):
    # Code here
    return nx_graph

def test_GrowNetwork():
    # Test code here
    G = nx.Graph()
    assert GrowNetwork(G) == []

# Main code block
def main():
    # Code here
    pass

def test_main():
    assert main() == ""

if __name__ == "__main__":
    main(options)
