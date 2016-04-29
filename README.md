### mepmap - Tools for creating and exploring metabolic reaction networks

#### Description
**mepmap_create.py** constructs a network of potential metabolic reactions using
resources supplied by KEGG (http://www.kegg.jp/) and MINE
(http://minedatabase.mcs.anl.gov/).

**mepmap_path.py** performs path-finding, sub-network extraction and enumeration
of putative biosynthetic pathways.

#### Dependencies
- Python ≥ 3.5.1
- Python NetworkX: https://github.com/networkx/networkx
- Python MINE-API: https://github.com/JamesJeffryes/MINE-API (included as **mineclient3.py**)

#### Author
Johannes Asplund-Samuelsson (<johannes.asplund.samuelsson@scilifelab.se>)
