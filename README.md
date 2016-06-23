![alt text](poppy.png "Prospecting Optimal Pathways with PYthon")
### Prospecting Optimal Pathways with PYthon
Tools for creating and exploring metabolic reaction networks.

#### Description
**poppy_create.py** constructs a network of potential metabolic reactions using
resources supplied by KEGG (http://www.kegg.jp/) and MINE
(http://minedatabase.mcs.anl.gov/).

**poppy_path.py** performs path-finding, sub-network extraction and enumeration
of putative biosynthetic pathways. Can write HTML reports for pathway
enumeration.

**poppy_rank.py** ranks identified pathways in terms of thermodynamic
limitations using the MDF (Max-min Driving Force) approach.

**poppy_gibb.py** downloads Equilibrator (http://equilibrator.weizmann.ac.il/)
transformed Gibbs standard formation energy values for a range of pHs and stores
them in a json file for use in the MDF analysis.

#### Dependencies
- Python ≥ 3.5.1
- Python NetworkX: https://github.com/networkx/networkx
- Python RDKit: https://github.com/rdkit/rdkit
- Python MINE-API: https://github.com/JamesJeffryes/MINE-API (included as **mineclient3.py**)
- R ≥ 3.0 with ggplot2 and RColorBrewer

#### Author
Johannes Asplund-Samuelsson (<johannes.asplund.samuelsson@scilifelab.se>)
