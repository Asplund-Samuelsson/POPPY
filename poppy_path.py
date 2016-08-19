#!/usr/bin/env python3

# Import modules
import networkx as nx
import sys
import os
import subprocess as sub
import argparse
import pickle
import multiprocessing as mp
import time
from datetime import timedelta as delta
import threading
import re
import pandas as pd
import gzip
from itertools import product
from shutil import copyfile

# Import scripts
from progress import Progress
import mineclient3 as mc
from poppy_origin_helpers import *
from poppy_helpers import *
import poppy_rank as rank
from poppy_create import extract_reaction_comp_ids

# Specify path to repository
global repo_dir
repo_dir = os.path.dirname(__file__)

# Define functions
def count_reactions(network):
    """Count the number of unique reactions in a network."""
    rxns = set()
    for n in network.nodes():
        if network.node[n]['type'] != 'c':
            rxns.add(network.node[n]['mid'])
    return len(rxns)


def find_paths(network, reactant_node, compound_node, reaction_limit):
    """Find all simple paths from an origin reactant node to a target
    compound node, limiting the total number of reactions."""

    paths = []

    if nx.has_path(network, reactant_node, compound_node):
        for path in nx.shortest_simple_paths(network, reactant_node, compound_node):
            if len(path)/3 > reaction_limit:
                # Dividing path length by three yields reaction step number
                break
            if nx.is_directed_acyclic_graph(network.subgraph(path)):
                # Only paths that represent an acyclic sub-network are allowed
                # Paths may not fold back onto themselves
                paths.append(path)

    return paths


def generate_paths(network, target_node, reaction_limit, n_procs=1, quiet=False):
    """Generate a list of paths to a target node from origin reactant nodes."""

    if not quiet:
        s_out("Generating paths...\n")

    # Define the worker
    def worker():
        while True:
            origin_node = Work.get()
            if origin_node is None:
                break
            paths = find_paths(network, origin_node, target_node, reaction_limit)
            output.extend(paths)
            with lock:
                n_work_done.value += 1

    # Set up reporter thread
    def reporter():
        n_left = Work.qsize()
        time_p = Progress(design = 't', max_val = n_work)
        prog_p = Progress(design = 'p', max_val = n_work)
        status_format = "{0:<10} {1:<25} {2:<25}"
        while n_left > 0:
            # Check progress
            n_made = len(output)
            n_left = Work.qsize()
            n_done = n_work - n_left
            # Calculate speed and time left
            progress = prog_p.to_string(n_done)
            found = str(n_made) + " paths found."
            t_out = "Time left: " + time_p.to_string(n_done)
            status = status_format.format(progress, found, t_out)
            s_out('\r' + status)
            time.sleep(1)

    with mp.Manager() as manager:
        # Initialize Work queue in manager
        Work = manager.Queue()

        for origin_node in find_valid_reactant_nodes(network):
            Work.put(origin_node)

        # Place stop signals on queue
        for i in range(n_procs):
            Work.put(None)

        n_work = Work.qsize()

        # Initialize output list in manager
        output = manager.list()

        # Initialize number of tasks done counter
        n_work_done = mp.Value('i', 0)
        lock = mp.Lock()

        # Start reporter
        if not quiet:
            reporter = threading.Thread(target=reporter)
            reporter.start()

        # Start processes
        procs = []
        for i in range(n_procs):
            p = mp.Process(target=worker)
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
        if not quiet:
            reporter.join()
            s_out("\nDone.\n")

        return list(output)


def nodes_being_produced(network):
    """Create set of compound nodes produced by the reactions in the network."""
    produced_nodes = set()
    for node in network.nodes():
        if network.node[node]['type'] in {'pf','pr'}:
            produced_nodes = produced_nodes.union(network.node[node]['c'])
    return produced_nodes


def nodes_being_consumed(network):
    """Create set of compound nodes consumed by the reactions in the network."""
    consumed_nodes = set()
    for node in network.nodes():
        if network.node[node]['type'] in {'rf','rr'}:
            consumed_nodes = consumed_nodes.union(network.node[node]['c'])
    return consumed_nodes


def remove_incomplete_reactions(network):
    """
    Removes reactions that require reactants not provided as start compounds or
    as products in other reactions of the network.

    Removes compound nodes if they are connected only to the affected reaction
    and are not a starting compound.

    The process is iterative, as removal of one reaction may make
    additional reactions 'incomplete' in terms of reactant presence.
    """

    start_comp_nodes = find_start_comp_nodes(network)

    while True:

        # Stop iteration if no nodes were removed in the previous round
        try:
            if prev_node_count == len(network.nodes()):
                break
        except NameError:
            # First round yields a NameError
            pass

        # Determine what compounds are available
        available = start_comp_nodes.union(nodes_being_produced(network))

        # Determine what reactions should be removed based on
        # unfulfilled reactant requirements
        nodes_to_remove = set()
        for node in network.nodes():
            if network.node[node]['type'] in {'rf','rr'}:
                if not network.node[node]['c'].issubset(available):
                    nodes_to_remove.add(node)
                    nodes_to_remove.add(network.successors(node)[0])
                    # Reactant nodes have one successor, i.e. a product node
                    # Go through and remove compounds directly downstream of the
                    # reaction, if they are connected only to this reaction
                    downstream = network.successors(network.successors(node)[0])
                    for comp_node in downstream:
                        if network.node[comp_node]['type'] != 'c':
                            s_err("Warning: '" + str(comp_node) + \
                            "' is not a compound node as expected (successor" +\
                            " of " + str(network.successors(node)[0]) + ")")
                            continue
                        else:
                            dependent = network.predecessors(comp_node) == \
                            network.successors(node)
                            start = network.node[comp_node]['start']
                            if dependent and not start:
                                nodes_to_remove.add(comp_node)

        # Count the number of nodes before purging
        prev_node_count = len(network.nodes())

        # Remove the nodes
        for node in nodes_to_remove:
            network.remove_node(node)


def find_branch_nodes(network, severed=False):
    """Identifies reactant nodes that represent branch points in the pathway
    network.

    IMPORTANT: This function assumes that non-start predecessors leading to a
    severed branch are not present in the network, as would be the case for
    compounds not found in a direct path from origin to target."""

    branch_nodes = set()

    for node in network.nodes():

        # Branch nodes are reactant nodes
        if network.node[node]['type'] in {'rf','rr'}:
            S = 0

            # Iterate over predecessors, which are compound nodes
            for predecessor in network.predecessors(node):
                if network.node[predecessor]['start']:
                    S += 1

            # Start compound nodes should always be included in the network
            # Severed non-start compound nodes might not be included
            N = len(network.node[node]['c']) - S

            # Branch reactant nodes have more than one predecessor
            if len(network.node[node]['c']) > 1:
                if not severed:
                    branch_nodes.add(node)
                else:
                    # A severed branch node lacks 1 or more predecessors
                    # compared to what is listed in the 'c' set
                    pred_num = len(network.predecessors(node))
                    c_num = len(network.node[node]['c'])
                    if pred_num < c_num:
                        branch_nodes.add(node)

    return branch_nodes


def find_switch_nodes(network):
    """
    Identifies and returns a list of 'switch nodes' in the network.

    Switch nodes are compound nodes with multiple predecessors and thereby
    multiple paths of synthesis.
    """
    switch_nodes = set()
    for node in network.nodes():
        if network.node[node]['type'] == 'c':
            if len(network.predecessors(node)) > 1:
                switch_nodes.add(node)
    return switch_nodes


def generate_termini(network):
    """
    Generates termini and eliminates start-compound induced branching by
    cutting all start compound to reactant node edges.
    """
    for node in network.nodes():
        if network.node[node]['type'] == 'c':
            if network.node[node]['start']:
                network.remove_edges_from(network.out_edges(node))


def digraph_connected_component(network, target_node):
    """
    Performs a reverse depth-first search to find all nodes that have a path
    leading to the target node. Then returns that component of the network.
    """
    return set(nx.dfs_preorder_nodes(network.reverse(), target_node))


def subnetwork_from_paths(network, paths, target_node):

    s_out("\nConstructing sub-network from identified paths...")

    # Acquire basic data
    path_nodes = set([n for p in paths for n in p])
    start_comp_nodes = find_start_comp_nodes(network)

    # Construct a sub-network of all paths,
    # start compounds and produced compounds
    subnet = network.subgraph(path_nodes.union(start_comp_nodes))
    subnet = network.subgraph(
        set(subnet.nodes()).union(nodes_being_produced(subnet))
    )

    n_removed = 1
    while n_removed:
        n_before = len(subnet)
        # Identify the incomplete reactions and remove them
        remove_incomplete_reactions(subnet)
        # Reduce to the connected component
        try:
            subnet = subnet.subgraph(digraph_connected_component(subnet, target_node))
        except KeyError:
            sys.exit("\nError: Subnetwork cannot generate target compound.\n")
        n_removed = n_before - len(subnet)
        # Exit with an error message if the target node was removed
        if target_node not in subnet.nodes():
            sys.exit("\nError: Target node cannot be produced by the sub-network.\n")

    s_out(" Done.\n")

    return subnet


def format_graphml(network, subnet):
    # Need to label the origin reactant nodes
    origins = find_valid_reactant_nodes(network)

    for node in subnet.nodes():
        if subnet.node[node]['type'] in {'rf','rr'}:
            if node in origins:
                subnet.node[node]['origin'] = True
            else:
                subnet.node[node]['origin'] = False

    # Also need to create compound name labels
    for node in subnet.nodes():
        if subnet.node[node]['type'] == 'c':
            mid = network.node[node]['mid']
            try:
                common_name = network.graph['mine_data'][mid]['Names'][0]
            except KeyError:
                common_name = network.graph['mine_data'][mid]['Formula']
            subnet.node[node]['common_name'] = common_name


    # Furthermore, label reaction nodes with the reactants or products
    for node in subnet.nodes():
        if subnet.node[node]['type'] != 'c':
            c_names = []
            for c_node in subnet.node[node]['c']:
                mid = network.node[c_node]['mid']
                try:
                    common_name = network.graph['mine_data'][mid]['Names'][0]
                except KeyError:
                    common_name = network.graph['mine_data'][mid]['Formula']
                c_names.append(common_name)
            subnet.node[node]['common_name'] = " + ".join(c_names)

    # Copy and remove the incompatible stuff from nodes and graph
    outnet = subnet.copy()

    for node in outnet.nodes():
        if outnet.node[node]['type'] != 'c':
            del outnet.node[node]['c']

    data_keys = list(outnet.graph.keys())
    for key in data_keys:
        del outnet.graph[key]

    return outnet


def paths_to_pathways(network, paths, target_node, rxn_lim=10, shallow=False):
    """Enumerate complete branched pathways capable of producing the target"""

    # Function for determining whether a path is part of a network
    def path_in_network(path, network):
        for i in range(1, len(path)):
            n = path[i]
            if n in network.nodes():
                predecessors = network.predecessors(n)
            else:
                predecessors = []
            if path[i-1] not in predecessors:
                return False
        return True

    # Construct a subnetwork
    subnet = subnetwork_from_paths(network, paths, target_node)

    # Determine compounds nodes available in the subnetwork
    start_comp_nodes = find_start_comp_nodes(network) # Start compounds
    prod_comp_nodes = nodes_being_produced(subnet)
    tot_comp = start_comp_nodes.union(prod_comp_nodes)

    # Filter the supplied paths to those that are present in the subnetwork
    paths_filtered = []
    p = Progress(design = 'pt', max_val = len(paths))
    n = 0
    for path in paths:
        n += 1
        s_out("\rFiltering paths... %s" % p.to_string(n))
        if path_in_network(path, subnet):
            paths_filtered.append(path)

    print("")

    # Generate a dictionary with path segments producing the key node
    segments = {}
    p = Progress(max_val = len(paths_filtered))
    m = 0

    # Go through all filtered paths
    for path in paths_filtered:

        # Report progress
        m += 1
        s_out("\rGenerating path segments... %s" % p.to_string(m))

        # Iterate over the elements of the path
        for element in enumerate(path):
            i = element[0]
            n = element[1]

            # Check if product node
            if subnet.node[n]['type'] in {'pf','pr'}:

                # Construct segments
                path_segs = []
                offset = 0
                while True:
                    path_segs.append(path[i - (1 + offset):i + 1])
                    # If shallow enumeration is desired
                    if shallow:
                        # Stop construction after a single step
                        break
                    try:
                        # If the upstream compound node is not start
                        if not subnet.node[path[i - (2 + offset)]]['start']:
                            # Stop construction
                            break
                    # If there are no more nodes upstream
                    except IndexError:
                        # Stop construction
                        break
                    offset += 3

                # Iterate over products
                for segment in path_segs:
                    for c_node in subnet.node[n]['c']:

                        # Save the segment(s) under every produced node
                        try:
                            segments[c_node].add(tuple(segment + [c_node]))
                        except KeyError:
                            segments[c_node] = set([tuple(segment + [c_node])])

    print("")

    # Pathway enumeration
    print("\nEnumerating pathways...")

    # Storage container for finished and unfinished pathways
    finished_pathways = set()
    try:
        unfinished_pathways = set([frozenset(p) for p in segments[target_node]])
    except KeyError:
        sys.exit("\nError: Subnetwork cannot generate target compound.\n")

    # Keep track of (partial) pathways that have already been generated
    generated_pathways = set()

    # Progress setup
    max_length = 0
    min_length = count_reactions(subnet)
    p = Progress(design='s')
    F = '{0} Finished: {1:<12} Unfinished: {2:<10} Reactions (min/max):{3:^4}/{4:^4}'

    def report_progress():
        D = len(finished_pathways)
        L = len(unfinished_pathways)
        p_msg = "\r" + F.format(p.to_string(), D, L, min_length, max_length)
        s_out(p_msg)

    # Iterate through unfinished pathways
    while unfinished_pathways:

        # Report progress
        report_progress()

        # Pop a path off the set of unfinished pathways
        path = set(unfinished_pathways.pop())

        # Extract the path's sub-network
        pathnet = subnet.subgraph(path)

        # Discard the pathway if it contains cycles
        # Cycles are 1) wasteful and 2) indicate bootstrap compounds
        if not nx.is_directed_acyclic_graph(pathnet):
            continue

        # Check the number of reactions
        n_rxn = count_reactions(pathnet)
        if n_rxn > rxn_lim:
            # Discard if exceeding the limit
            continue

        # Determine what compounds are produced by the path
        cpd_prod = nodes_being_produced(pathnet)

        # Determine what compounds are consumed by the path
        cpd_cons = nodes_being_consumed(pathnet)

        # Check if the path is complete on its own
        if cpd_cons.issubset(cpd_prod.union(start_comp_nodes)):
            max_length = max(max_length, n_rxn)
            min_length = min(min_length, n_rxn)
            finished_pathways.add(frozenset(path))

        # If not, add all combinations of segments that might complement it
        else:

            # Identify the missing compound nodes
            missing = cpd_cons - cpd_prod - start_comp_nodes

            # Construct sets of complementary segments that will satisfy the
            # current missing compound nodes
            try:
                for complement in product(*[segments[i] for i in missing]):
                    new_pw = frozenset(set(path).union(*complement))
                    if new_pw not in generated_pathways:
                        if count_reactions(subnet.subgraph(new_pw)) <= rxn_lim:
                            unfinished_pathways.add(new_pw)
                            generated_pathways.add(new_pw)

            except KeyError:
                # If a complement cannot be found, discard the pathway
                pass

    # Report final progress
    report_progress()

    print("")

    return finished_pathways


def parse_compound(compound, network):
    """Determines the type of compound identifier and returns the node."""

    node = None

    mid_match = re.match('^[CX]{1}[0-9,a-f]{40}$', compound)
    kegg_match = re.match('^C{1}[0-9]{5}$', compound)

    if mid_match:
        try:
            node = network.graph['cmid2node'][compound]
            if not node in network.nodes():
                s_err("Error: MINE ID '" + compound + \
                "' appears to not be available in the network.\n")
                node = None
        except KeyError:
            s_err("Error: MINE ID '" + compound + \
            "' appears to not be available in the network.\n")
    elif kegg_match:
        try:
            nodes = network.graph['kegg2nodes'][compound]
            if len(nodes) > 1:
                nodes = "'\n'".join(
                    sorted([network.node[n]['mid'] for n in nodes])
                )
                s_err("Error: '" + compound + "' refers to multiple IDs." + \
                " Use --exact_comp_id and one of:\n'" + \
                nodes + "'\n")
            else:
                node = list(nodes)[0]
                if not node in network.nodes():
                    s_err("Error: KEGG ID '" + compound + \
                    "' appears to not be available in the network.\n")
                    node = None
        except KeyError:
            s_err("Error: KEGG ID '" + compound + \
            "' appears to not be available in the network.\n")
    else:
        try:
            nodes = network.graph['name2nodes'][compound]
            if len(nodes) > 1:
                nodes = "'\n'".join(
                    sorted([network.node[n]['mid'] for n in nodes])
                )
                s_err("Error: '" + compound + "' refers to multiple IDs." + \
                " Use --exact_comp_id and one of:\n'" + \
                nodes + "'\n")
            else:
                node = list(nodes)[0]
                if not node in network.nodes():
                    s_err("Error: Name '" + compound + \
                    "' appears to not be available in the network.\n")
                    node = None
        except KeyError:
            s_err("Error: Name '" + compound + \
            "' appears to not be available in the network.\n")

    return node


def update_start_compounds(network, start_comp_ids):
    """Update the start compound status for each compound node."""
    S = set(start_comp_ids)
    for n in network.nodes():
        if network.node[n]['type'] == 'c':
            if network.node[n]['mid'] in S:
                start = True
            else:
                start = False
            network.node[n]['start'] = start


def format_reaction_text(reaction, reverse=False):
    if not reverse:
        r = reaction['Reactants']
        p = reaction['Products']
    else:
        r = reaction['Products']
        p = reaction['Reactants']

    # Flatten lists and add plus characters
    r = [a for b in list(joinit(r, ['+'])) for a in b]
    p = [a for b in list(joinit(p, ['+'])) for a in b]

    # Combine and add reaction arrow
    rxn_elements = r + ['<=>'] + p

    # Remove "1" coefficients
    rxn_elements = list(filter(lambda x : x is not 1, rxn_elements))

    # Create text
    return " ".join([str(x) for x in rxn_elements])


def format_pathway_text(network, pathways, target_node, pw_sep=True):
    """Create pathways record in text format"""
    pathway_lines = []

    # Go through pathways in length order
    for pathway in sorted(list(pathways), key=len):

        # Create a subnetwork
        subnet = network.subgraph(pathway)

        # Identify the reactant nodes in the pathway
        reactant_nodes = [
            n[0] for n in subnet.nodes(data=True) if n[1]['type'] in {'rf','rr'}
        ]

        # Find the order of the reactions
        def distance(rxn_node):
            if nx.has_path(subnet, rxn_node, target_node):
                return len(nx.shortest_path(subnet, rxn_node, target_node))
            else:
                return 0

        rxn_nodes = sorted(reactant_nodes, key = distance, reverse = True)

        # Add reactions in the detected order
        for n in rxn_nodes:
            rxn_type = subnet.node[n]['type']
            rxn_id = subnet.node[n]['mid']
            rxn = network.graph['mine_data'][rxn_id]
            if rxn_type == 'rf':
                rxn_text = format_reaction_text(rxn)
            else:
                rxn_text = format_reaction_text(rxn, reverse=True)

            # Create text
            pathway_lines.append(rxn_id + "\t" + rxn_text)

        # Add pathway divider
        if pw_sep:
            pathway_lines.append("//")

    # Return formatted pathways text
    return "\n".join(pathway_lines) + "\n"


def format_mdf_summary(mdf_dict, network):

    # Set up column series
    pw_hash = []
    pw_nrxn = []
    pw_mdfs = []
    pw_rxns = []
    pw_txts = []

    # Iterate over pathways in the MDF dictionary
    for pw_txt in mdf_dict:

        # Append hash
        pw_hash.append(rank.generate_pathway_hash(pw_txt))

        # Append number of reactions
        pw_nrxn.append(len(set([
            x.split("\t")[0] for x in filter(None, pw_txt.split("\n"))
        ])))

        # Append MDF
        pw_mdfs.append(mdf_dict[pw_txt])

        # Generate list of reaction IDs with directionality
        rxn_ids = []
        for rxn_line in filter(None, pw_txt.split("\n")):
            rxn_id = rxn_line.split("\t")[0]
            rxn_txt = rxn_line.split("\t")[1]
            rxn = network.graph['mine_data'][rxn_id]

            # Check if the reaction is presented in reverse direction
            if rxn_txt != format_reaction_text(rxn):
                rxn_id = "-" + rxn_id

            # Append to list of reactions in correct order
            rxn_ids.append(rxn_id)

        pw_rxns.append(",".join(rxn_ids))

        # Append pathway text
        pw_txts.append(pw_txt)

    # Construct dataframe
    pw_df = pd.DataFrame({
        'pathway'   : pw_hash,
        'length'    : pw_nrxn,
        'MDF'       : pw_mdfs,
        'reactions' : pw_rxns,
        'pw_txt'   : pw_txts
    }, columns = ['pathway','MDF','length','reactions','pw_txt'])

    # Sort the dataframe
    pw_df = pw_df.sort_values(['MDF', 'length'], ascending=[False, True])

    # Turn dataframe into text
    text_output = pw_df.to_csv(
        na_rep='NA', index=False, float_format='%.3f', sep='\t',
        columns = ['pathway','MDF','length','reactions']
    )
    return (pw_df, text_output)


def format_pathway_html(pw_df, network, target_node, depth, rxn_lim, n_pw=200):

    html = []

    # Reduce dataframe to top pathways
    pw_df_top = pw_df.head(n_pw)

    # Gather general metrics and descriptive information
    number_of_pathways = pw_df.shape[0]

    timestamp = time.strftime('%Y-%m-%d %H:%M:%S')

    target_cpd_id = network.node[target_node]['mid']
    try:
        target_cpd_name = network.graph['mine_data'][target_cpd_id]['Names'][0]
    except (KeyError, IndexError) as e:
        target_cpd_name = target_cpd_id + '(no name)'

    n_feasible = sum(pw_df['MDF'] > 0)
    n_infeasible = number_of_pathways - n_feasible

    # Add the HTML header section
    html.extend([
        '<!DOCTYPE html>',
        '<html>',
        '<head>',
        '  <link rel="stylesheet" href="style.css">',
        '<base target="_blank">',
        '</head>',
        '<body>',
        ''
    ])

    # Add the summary section of the report
    html.extend([
        '<table border=0>',
        '<tr>',
        '  <td>',
        '  <h3><b>Report: %i %s pathways</b></h3>' % (number_of_pathways, target_cpd_name),
        '  </td>',
        '  <td>',
        '    <h5><b>MDF</b></h5>',
        '  </td>',
        '  <td>',
        '    <h5><b>Pathway length</b></h5>',
        '  </td>',
        '</tr>',
        '<tr>',
        '  <td>',
        '  <li><h5><b>Report generated on:</b> %s</h5>' % timestamp,
        '  <li><h5><b>Target compound:</b> %s (%s)</h5>' % (target_cpd_name, target_cpd_id),
        '  <li><h5><b>Reaction depth:</b> %i</h5>' % depth,
        '  <li><h5><b>Reaction limit:</b> %i</h5>' % rxn_lim,
        '  <li><h5><b>Number of pathways:</b> %i</h5>' % number_of_pathways,
        '  <li><h5><b>Feasible/infeasible:</b> %i/%i</h5>' % (n_feasible, n_infeasible),
        '  </td>',
        '  <td>',
        '  <img src="mdf_summary.png" style="width: 200px;"/>',
        '  </td>',
        '  <td>',
        '  <img src="length_summary.png" style="width: 200px;"/>',
        '  </td>',
        '</tr>',
        '</table>',
        ''
    ])

    # Define some formatting functions
    def kegg_url(kegg_id):
        return 'http://www.genome.jp/dbget-bin/www_bget?' + kegg_id

    def cpd_td(cpd_type):
        color = {'o':'#c7eae5','i':'#f6e8c3','t':'#c2a5cf'}[cpd_type]
        return '<td style="background-color:%s">' % color

    def enz_link(operator):
        if operator.startswith('M:'):
            operator = operator.lstrip('M:') + '*'
        op_list = [n.strip('-') for n in operator.split(".")]
        if not re.match('^[0-9]+$', op_list[-1]):
            op_list[-1] = '-'
        enz_url = "http://enzyme.expasy.org/EC/" + '.'.join(op_list)
        return '<a href="%s">%s</a>' % (enz_url, operator)

    def insert_img(n, cpd_id):
        if n > 1:
            img_str_l = str(n) + ' <img src="cpd_png/'
        else:
            img_str_l = '<img src="cpd_png/'
        return img_str_l + '%s.png" style="width: 75px;"/>' % cpd_id

    # Add report for each pathway
    for pw_row in pw_df_top.iterrows():
        # Add HR line
        html.extend(['<hr>',''])

        # Add title and MDF
        if pd.isnull(pw_row[1]['MDF']):
            m = 'MDF FAILED'
        else:
            m = '%.3f kJ/mol' % pw_row[1]['MDF']
        p = pw_row[1]['pathway']
        n = pw_row[1]['length']
        html.append('<h5>Pathway "%s": %s, %i reactions</h5>' % (p, m, n))
        html.append('')

        # Iterate over reactions, add one table for each
        for rxn_line in filter(None, pw_row[1]['pw_txt'].split("\n")):

            # Table start
            html.extend(['<table border=0>','<tr>'])

            # Parse elements of the reaction line
            rxn_id = rxn_line.split("\t")[0]
            rxn_eq = rank.parse_equation(rxn_line.split("\t")[1])

            # Add link to reaction
            if rxn_id.startswith('RM'):
                html.append('  <td>%s</td>' % rxn_id)
            else:
                url = kegg_url(rxn_id)
                html.append('  <td><a href="%s">%s</a></td>' % (url, rxn_id))

            # Iterate over compounds for the first table row
            name_cells = []
            for cpd_id in [e[1] for e in rxn_eq[0] + rxn_eq[1]]:

                # Identify the compound type (target, origin or intermediate)
                cpd_node = network.graph['cmid2node'][cpd_id]

                if cpd_node == target_node:
                    cpd_type = 't'
                elif network.node[cpd_node]['start']:
                    cpd_type = 'o'
                else:
                    cpd_type = 'i'

                # Generate KEGG URL and name
                url = kegg_url(cpd_id)
                try:
                    cpd_name = network.graph['mine_data'][cpd_id]['Names'][0]
                except (KeyError, IndexError) as e:
                    cpd_name = cpd_id + '(no name)'

                # Add table cells
                name_cells.append(
                    ''.join(['  ', cpd_td(cpd_type), '<a href="',
                             url, '">', cpd_name, '</a></td>'])
                )

            html.extend(joinit(name_cells, '  <td></td>'))

            html.append('</tr>')

            # Add second row of table (compound images)
            html.append('<tr>')

            # Add links to Expasy for the EC classes/operators
            operators = network.graph['mine_data'][rxn_id]['Operators']
            html.append('  <td><p>')
            html.extend(['    ' + enz_link(op) + '<br>' for op in operators])
            html.append('  </p></td>')

            # Iterate over compounds for the second table row
            img_cells = []
            for el in rxn_eq[0]:
                img_cells.append('  <td>' + insert_img(el[0], el[1]) + '</td>')
            html.extend(joinit(img_cells, '  <td>+</td>'))

            html.append('  <td><=></td>')

            img_cells = []
            for el in rxn_eq[1]:
                img_cells.append('  <td>' + insert_img(el[0], el[1]) + '</td>')
            html.extend(joinit(img_cells, '  <td>+</td>'))
            html.append('</tr>')

            # End table
            html.extend(['</table>',''])

    return "\n".join(html)


# Main code block
def main(infile_name, compound, start_comp_id_file, exact_comp_id,
    rxn_lim, depth, n_procs, sub_network_out, pathway_pickle, shallow,
    pathway_text, pathway_html, n_pw_out, bounds, ratios, dfG_json, net_file,
    pH, T, R):

    # Default results are empty
    results = {}

    # Load the network
    s_out("\nLoading network pickle...")
    network = pickle.load(open(infile_name, 'rb'))
    s_out(" Done.\n")

    # Update starting compounds
    if start_comp_id_file:
        S = [L.rstrip() for L in open(start_comp_id_file, 'r').readlines()]
        update_start_compounds(network, S)

    # Pathway enumeration
    if not exact_comp_id:
        target_node = parse_compound(compound, network)
    else:
        try:
            target_node = network.graph['cmid2node'][compound]
        except KeyError:
            target_node = None

    if target_node == None:
        sys.exit("Error: Target node was not found. Check compound '" + \
        compound + "'.\n")

    paths = generate_paths(network, target_node, depth, n_procs)

    # Save sub-network graphml
    if sub_network_out:
        subnet = subnetwork_from_paths(network, paths, target_node)
        s_out("\nWriting sub-network to graphml...")
        subnet = format_graphml(network, subnet)
        nx.write_graphml(subnet, sub_network_out)
        s_out(" Done.\n")

    # Enumerate pathways and save results
    if pathway_pickle or pathway_text or pathway_html:

        pathways = paths_to_pathways(
            network, paths, target_node, rxn_lim, shallow
        )

        # If there are no pathways, exit with an error message
        if not pathways:
            sys.exit("\nError: Could not construct any complete pathways.\n")

        if pathway_pickle:
            s_out("\nWriting pathways to pickle...")
            pickle.dump(pathways, open(pathway_pickle, 'wb'))
            s_out(" Done.\n")

        if pathway_text:
            s_out("\nWriting pathways to text...")
            with open(pathway_text, 'w') as txt:
                txt.write(format_pathway_text(network, pathways, target_node))
            s_out(" Done.\n")

        if pathway_html:

            # Format pathways so that ranking functions can read them
            pw_rank_txt = format_pathway_text(network, pathways, target_node)
            pw_rank = rank.read_pathways_text(pw_rank_txt)

            # Load standard formation Gibbs energy dictionary
            dfG_dict = rank.load_dfG_dict(pw_rank, pH, dfG_json)

            # Load metabolic network text
            if net_file:
                net_text = open(net_file, 'r').read()
            else:
                net_text = ""

            # Load inequality constraints
            if bounds:
                ne_con_text = open(bounds,'r').read()
            else:
                ne_con_text = ""

            ne_con = rank.mdf.read_constraints(ne_con_text)

            # Load equality constraints
            if ratios:
                eq_con_text = open(ratios,'r').read()
            else:
                eq_con_text = ""

            eq_con = rank.mdf.read_ratio_constraints(eq_con_text)

            # Perform MDF
            mdf_dict = rank.pathways_to_mdf(
                pw_rank, dfG_dict, ne_con, eq_con, n_procs, T, R, net_text
            )

            # Create output directory structure
            out_dir = os.path.abspath(pathway_html)
            rxn_file = os.path.join(out_dir, 'reactions.txt')
            pw_file = os.path.join(out_dir, 'pathways.txt')
            mdf_fig = os.path.join(out_dir, 'mdf_summary.png')
            len_fig = os.path.join(out_dir, 'length_summary.png')
            html_file = os.path.join(out_dir, 'report.html')
            css_file = os.path.join(out_dir, 'style.css')

            try:
                # Exit if the directory already exists
                os.stat(out_dir)
                sys.exit("Error: HTML output directory already exists.\n")
            except:
                os.mkdir(out_dir)

            os.mkdir(os.path.join(out_dir, 'cpd_png'))

            # Format and write reactions to text file
            pw_nodes = set([n for p in pathways for n in p])
            pw_rxn_ids = sorted(list(set([
                network.node[n]['mid'] for n in pw_nodes \
                if network.node[n]['type'] != 'c'
            ])))

            with open(rxn_file, 'w') as f:
                for rxn_id in pw_rxn_ids:
                    rxn = network.graph['mine_data'][rxn_id]
                    rxn_txt = format_reaction_text(rxn)
                    f.write(rxn_id + "\t" + rxn_txt + "\n")

            # Format and write pathways summary to text file
            pw_df, pw_sum = format_mdf_summary(mdf_dict, network)
            with open(pw_file, 'w') as f:
                f.write(pw_sum)

            # Create summary figures; call external R script
            sub.call([
                'Rscript', '--vanilla',
                os.path.join(repo_dir, 'poppy_sumfigs.R'),
                pw_file,
                mdf_fig,
                len_fig
            ])

            # Load image dictionary
            cpd_img_file = os.path.join(repo_dir, 'data/mol_png_images.pkl.gz')
            cpd_images = pickle.load(gzip.open(cpd_img_file))

            # Populate compound image directory
            compounds = set()
            for pw_row in pw_df.head(n_pw_out).iterrows():
                rxn_rows = filter(None, pw_row[1]['pw_txt'].split("\n"))
                for rxn_id in [r.split("\t")[0] for r in rxn_rows]:
                    rxn = network.graph['mine_data'][rxn_id]
                    cpd_ids = extract_reaction_comp_ids(rxn)
                    compounds = compounds.union(cpd_ids)

            for cpd in compounds:
                img_file = os.path.join(out_dir, 'cpd_png', cpd + '.png')
                with open(img_file, 'wb') as f:
                    try:
                        f.write(cpd_images[cpd])
                    except KeyError:
                        s_err("Warning: No image found for '%s'." % cpd)

            # Add CSS file to HTML directory
            css_loc = os.path.join(repo_dir, 'data/style.css')
            copyfile(css_loc, css_file)

            # Format and write pathways summary with top X to HTML
            html = format_pathway_html(
                        pw_df, network, target_node,
                        depth, rxn_lim, n_pw_out
                        )

            with open(html_file, 'w') as f:
                f.write(html)


if __name__ == "__main__":
    # Read arguments from the commandline
    parser = argparse.ArgumentParser()

    # Required input: Reaction network and target compound
    parser.add_argument(
        'infile',
        help='Read reaction network pickle.'
    )
    parser.add_argument(
        'compound', type=str, default=False,
        help='Target compound.'
    )

    # Number of processes to run for path-finding and MDF analysis
    parser.add_argument(
        '-p', '--processes', type=int, default=1,
        help='Number of parallel processes to run.'
    )

    # Pathway enumeration parameters
    parser.add_argument(
        '-S', '--start_comp_ids', type=str,
        help='Provide new starting compounds in a file.'
    )
    parser.add_argument(
        '-e', '--exact_comp_id', action='store_true',
        help='Look for exact compound ID.'
    )
    parser.add_argument(
        '-r', '--reactions', type=int, default=10,
        help='Maximum number of reactions.'
    )
    parser.add_argument(
        '-d', '--depth', type=int, default=5,
        help='Maximum direct path depth.'
    )
    parser.add_argument(
        '--shallow', action='store_true', default=False,
        help='Shallow enumeration; do not go through start compounds.'
    )

    # Output options
    parser.add_argument(
        '--pathway_pickle', type=str, default=False,
        help='Save identified pathways in pickle.'
    )
    parser.add_argument(
        '--pathway_text', type=str, default=False,
        help='Save identified pathways as a text file.'
    )
    parser.add_argument(
        '--pathway_html', type=str, default=False,
        help='Save ranked pathway html output to directory.'
    )
    parser.add_argument(
        '--sub_network', type=str, default=False,
        help='Save sub-network as graphml.'
    )
    parser.add_argument(
        '--n_html_pathways', type=int, default=200,
        help='The number of pathways to save in HTML output.'
    )

    # MDF analysis
    parser.add_argument(
        '--bounds', type=str,
        help='Read metabolite concentration bounds (inequality constraints).'
    )
    parser.add_argument(
        '--ratios', type=str,
        help='Read metabolite concentration ratios (equality constraints).'
    )
    parser.add_argument(
        '--gibbs', type=str,
        help='Read transformed Gibbs standard formation energies (JSON).'
    )
    parser.add_argument(
        '--model', type=str,
        help='Read reactions corresponding to a reduced network model.'
    )
    parser.add_argument(
        '--pH', type=float, default=7.0,
        help='Specify the pH for the thermodynamics calculations.'
    )
    parser.add_argument(
        '-T', type=float, default=298.15,
        help='Temperature (K).'
    )
    parser.add_argument(
        '-R', type=float, default=8.31e-3,
        help='Universal gas constant (kJ/(mol*K)).'
    )

    args = parser.parse_args()

    # Run main function
    main(
        args.infile, args.compound, args.start_comp_ids, args.exact_comp_id,
        args.reactions, args.depth, args.processes, args.sub_network,
        args.pathway_pickle, args.shallow, args.pathway_text, args.pathway_html,
        args.n_html_pathways, args.bounds, args.ratios, args.gibbs, args.model,
        args.pH, args.T, args.R
    )
