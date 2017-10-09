![alt text](poppy.png "Prospecting Optimal Pathways with PYthon")

# Prospecting Optimal Pathways with PYthon

Tools for creating and exploring metabolic reaction networks.

---

### 1. Construct reaction network

`poppy_create.py` constructs a network of potential metabolic reactions using
resources supplied by KEGG ([http://www.kegg.jp/](http://www.kegg.jp/)) and MINE
([http://minedatabase.mcs.anl.gov/](http://minedatabase.mcs.anl.gov/)).

##### _Example: Create a combined KEGG and MINE reaction network_

`./poppy_create.py --kegg --enhance --equilibrator_filter --infile examples/E_coli.origins.txt network.pkl`

---

### 2. Obtain thermodynamic data

`poppy_gibb.py` downloads Equilibrator ([http://equilibrator.weizmann.ac.il/](http://equilibrator.weizmann.ac.il/))
transformed Gibbs standard formation energy values for specified pHs and stores
them in a JSON file for use in the thermodynamic analyses.

##### _Example: Obtain formation delta G:s for_ E. coli _and_ Synechocystis _pH:s_

`./poppy_gibb.py 7.6,8.4 dfgs.json`

We modify the JSON file with values for ferredoxin and carboxysomal metabolites:

`perl -pi -e 's/null/19.9/ if /C00138/' dfgs.json`

`perl -pi -e 's/null/-19.9/ if /C00139/' dfgs.json`

`perl -pi -e 'my $I = $_; my $O = $I =~ s/\"\:/_a\"\:/gr; print $O if $I =~ /C00001|C00007|C00011|C00197|C00988|C01182/' dfgs.json`

---

### 3. Enumerate and evaluate pathways

`poppy_path.py` performs path-finding, sub-network extraction and enumeration
of putative biosynthetic pathways. Can write HTML reports for pathway
enumeration.

##### _Example: Enumerate pathways to 4-hydroxybutanoic acid in_ E. coli _and_ Synechocystis

`./poppy_path.py -p 4 -d 3 -r 5 --model examples/E_coli.model.tab --bounds examples/E_coli.concentrations.tab --ratios examples/E_coli.ratios.tab --gibbs dfgs.json --pH 7.6 --c_min 0.0000001 --c_max 0.1 --pathway_html E_coli_pathways network.pkl C00989`

`./poppy_path.py -p 4 -d 3 -r 5 -S examples/Synechocystis.origins.txt --model examples/Synechocystis.model.tab --bounds examples/Synechocystis.concentrations.tab --ratios examples/Synechocystis.ratios.tab --gibbs dfgs.json --pH 8.4 --c_min 0.0000001 --c_max 0.1 --pathway_html Synechocystis_pathways network.pkl C00989`

---

### 4. Calculate model network reaction Gibbs free energy changes

`poppy_rank.py` ranks identified pathways in terms of thermodynamic
driving forces. Also used for calculating transformed standard reaction Gibbs
free energy changes.

##### _Example: Calculate reaction delta G:s for the_ E. coli _and_ Synechocystis _models_

`./poppy_rank.py --write_gibbs --gibbs dfgs.json --pH 7.6 examples/E_coli.model.tab E_coli.model_drgs.tab`

`./poppy_rank.py --write_gibbs --gibbs dfgs.json --pH 8.4 examples/Synechocystis.model.tab Synechocystis.model_drgs.tab`

---

### 5. Standalone MDF and NEM analysis

`mdf.py` performs Max-min Driving Force (MDF; [Noor _et al._, 2014](http://doi.org/10.1371/journal.pcbi.1003483)) and Network-Embedded
MDF (NEM) analysis.

##### _Example: Perform NEM analysis on lysine biosynthesis in_ E. coli _and_ Synechocystis

`./mdf.py --min_conc 0.0000001 --max_conc 0.1 --constraints examples/E_coli.concentrations.tab --ratios examples/E_coli.Lys_opt_ratios.tab --pathway examples/E_coli.Lys_pathway.txt examples/E_coli.model.tab E_coli.model_drgs.tab E_coli_Lys_nem.csv`

`./mdf.py --min_conc 0.0000001 --max_conc 0.1 --constraints examples/Synechocystis.concentrations.tab --ratios examples/Synechocystis.Lys_opt_ratios.tab --pathway examples/Synechocystis.Lys_pathway.txt examples/Synechocystis.model.tab Synechocystis.model_drgs.tab Synechocystis_Lys_nem.csv`

---

## Dependencies
- Python ≥ 3.5.1
- Python SciPy and NumPy ([https://www.scipy.org/](https://www.scipy.org/))
- Python NetworkX ([https://github.com/networkx/networkx](https://github.com/networkx/networkx))
- Python RDKit ([https://github.com/rdkit/rdkit](https://github.com/rdkit/rdkit))
- Python MINE-API ([https://github.com/JamesJeffryes/MINE-API](https://github.com/JamesJeffryes/MINE-API); included as `mineclient3.py`)
- R ≥ 3.0 with ggplot2 and RColorBrewer

---

## Author
Johannes Asplund-Samuelsson (<johannes.asplund.samuelsson@scilifelab.se>)
