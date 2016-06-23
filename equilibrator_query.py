#!/usr/bin/env python3

# Import modules
import re
import threading
import queue
import sys
from requests import get as rget
from time import sleep

# Define functions
def equilibrator_gibbf(kegg_comp_id, pH=7.0, ionic_strength=0.1):
    """
    Performs an Equilibrator query with the supplied parameters, returning a
    tuple containing the Gibb's formation energy and the error.
    """
    address_root = "http://equilibrator.weizmann.ac.il/compound"
    c = "?compoundId=" + kegg_comp_id
    p = "&ph=" + str(pH)
    i = "&ionic_strength=" + str(ionic_strength)

    address = address_root + c + p + i

    n = 0

    while n < 5:
        n += 1
        r = rget(address)

        if r.status_code == 200:
            # The Delta G value might not be available
            if "Cannot estimate" in r.text:
                return None

            # Split into table rows
            eq_data = re.sub("[\n \t]*", "", r.text)
            eq_data = re.split("</tr>", eq_data)

            # Find rows containing "Delta"
            is_delta_line = re.compile("Delta")
            delta_lines = list(filter(is_delta_line.search, eq_data))


            # Extract standard formation Gibb's line
            try:
                delta_line_2 = delta_lines[1]
            except IndexError:
                # Delta G is not available
                return None

            # Clean up and split line
            delta_line_2 = list(filter(None,re.split("<.*?>|&[A-Za-z]+;|\[", delta_line_2)))

            # Extract floating point numbers
            is_float_string = re.compile("^-*\d+?\.*\d*?$")
            return tuple([float(x) for x in filter(is_float_string.match, delta_line_2)])
        elif r.status_code == 500:
            return None
        else:
            sleep(2)


def threaded_equilibrator_gibbf(queries):
    """Threaded implementation of equilibrator_gibbf."""

    def worker():
        while True:
            query = work.get()
            if query is None:
                break
            try:
                result = equilibrator_gibbf(*query)
                if result is None:
                    out = (query, None)
                else:
                    out = (query, result)
                output.put(out)
                work.task_done()
            except:
                # In case of exception, put query back on queue
                work.put(query)
                work.task_done()

    # Initialise queues
    work = queue.Queue()
    output = queue.Queue()

    for query in queries:
        work.put(query)

    # Start threads
    threads = []
    for i in range(16):
        t = threading.Thread(target=worker)
        t.start()
        threads.append(t)

    # Report on progress
    while True:
        progress = float(output.qsize() / len(queries) * 100)
        sys.stdout.write("\rDownloading standard formation energies... %0.1f%%" % progress)
        sys.stdout.flush()
        if output.qsize() == len(queries):
            print("")
            break
        sleep(0.5)

    # Join work queue
    work.join()

    # Stop workers
    for i in range(len(threads)):
        work.put(None)
    for t in threads:
        t.join()

    # Get results
    results = {}
    while not output.empty():
        result = output.get()
        results[result[0]] = result[1]

    return results
