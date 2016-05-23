### mepmap - Tools for creating and exploring metabolic reaction networks

#### Description
**mepmap_create.py** constructs a network of potential metabolic reactions using
resources supplied by KEGG (http://www.kegg.jp/) and MINE
(http://minedatabase.mcs.anl.gov/).

**mepmap_path.py** performs path-finding, sub-network extraction and enumeration
of putative biosynthetic pathways.

**mepmap_rank.py** ranks identified pathways in terms of thermodynamic
limitations using the MDF (Max-min Driving Force) approach.

**mepmap_gibb.py** downloads Equilibrator (http://equilibrator.weizmann.ac.il/)
transformed Gibbs standard formation energy values for a range of pHs and stores
them in a json file for use in the MDF analysis (mepmap_rank.py).

#### Dependencies
- Python ≥ 3.5.1
- Python NetworkX: https://github.com/networkx/networkx
- Python RDKit: https://github.com/rdkit/rdkit
- Python MINE-API: https://github.com/JamesJeffryes/MINE-API (included as **mineclient3.py**)

#### Author
Johannes Asplund-Samuelsson (<johannes.asplund.samuelsson@scilifelab.se>)
