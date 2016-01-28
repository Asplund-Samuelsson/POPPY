## minet - Tool for creating a NetworkX graph from MINE data

#### NAME
minet

#### SYNOPSIS
**minet** [-s *reaction_steps*] *infile.txt* *outfile.pickle*

#### DESCRIPTION
**minet** reads a list of KEGG compound identifiers (e.g. C00469), then contacts the MINE database "KEGGexp2" and creates a NetworkX graph object from the connected reactions and compounds. The graph is saved as a Python Pickle. The maximum number of iterations (reaction steps) can be specified by the user. The tool is implemented in Python 3.5.1.

#### OPTIONS
-s *reaction_steps*
Number of iterations (reaction steps) to use when building the network.

#### DEPENDENCIES
Python NetworkX: https://github.com/networkx/networkx
Python MINE-API: https://github.com/JamesJeffryes/MINE-API

#### AUTHOR
Johannes Asplund-Samuelsson johannes.asplund.samuelsson@scilifelab.se
