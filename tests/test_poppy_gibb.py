#!/usr/bin/env python3

# Add repository root to the path
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

# Import the script to be tested
from poppy_gibb import *

# Define tests
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
