#!/usr/bin/env python3

# Add repository root to the path
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

# Import the script to be tested
from equilibrator_query import *

# Define tests
def test_equilibrator_gibbf():
    # http://equilibrator.weizmann.ac.il/compound?compoundId=C09908&ph=6.75&ionic_strength=0.175
    assert equilibrator_gibbf("C09908", pH=6.75, ionic_strength=0.175) == (482.2, 7.1)
    assert equilibrator_gibbf("C06142") == (234.0, 4.8)
    assert equilibrator_gibbf("C06142", 4) == (62.9, 4.8)
    assert equilibrator_gibbf("C00229") is None # ACP (C00229) does not have a value
    assert equilibrator_gibbf("C99999") is None
    assert equilibrator_gibbf("C01102", 8.75, 0.125) == (-955.8, 3.0)


def test_threaded_equilibrator_gibbf():
    queries = [
        ("C09908", 6.75, 0.175),
        ("C06142",),
        ("C06142",4),
        ("C00229",),
        ("C99999",),
        ("C01102", 8.75, 0.125)
    ]
    results = {
        ("C09908", 6.75, 0.175):(482.2, 7.1),
        ("C06142",):(234.0, 4.8),
        ("C06142",4):(62.9, 4.8),
        ("C00229",):None,
        ("C99999",):None,
        ("C01102", 8.75, 0.125):(-955.8, 3.0)
    }
    assert threaded_equilibrator_gibbf(queries) == results
