MinPath version 1.6 as of May 11, 2021
(older versions: MinPath1.5 released on Sep 21, 2020, MinPath1.1 released on July 31, 2010; MinPath1.0 released on July 3, 2009; MinPath1.2 on Oct 20, 2010, MinPath1.4 on Sep 11, 2018; older releases are available at https://omics.informatics.indiana.edu/MinPath/)

Developer: Yuzhen Ye <yye@indiana.edu>
Affiliation: School of Informatics, Computing and Engineering, Indiana University, Bloomington

0. What's new
  MinPath1.6: updated KEGG and MetaCyc Pathways as of May 10, 2021
  MinPath1.5: python 2 & 3 compatible, with updated KEGG and MetaCyc pathways
  MinPath1.4: fixed a bug that did not work with fig families
  MinPath1.2: works with any pathway system (see 3.1.3; and usages)
  MinPath1.1: outputs detailed pathway assignments 

1. Introduction
MinPath (Minimal set of Pathways) is a parsimony approach for biological pathway reconstructions 
using protein family predictions, achieving a more conservative, yet more faithful, estimation of 
the biological pathways for a query dataset.

MinPath is free software under the terms of the GNU General Public
License as published by the Free Software Foundation.

If you use MinPath, please cite,
Yuzhen Ye and Thomas G. Doak. A parsimony approach to biological pathway reconstruction/inference 
for genomes and metagenomes (PLoS Computational Biology, 5(8):e1000465, 2009)

2. Installation
Once you clone the repository, it is ready to be used. 
Check 3.2 if glpk doesn't work.

3. Requirements
 3.1 files of reference pathways with function to pathway (subsystem) mapping
  3.1.1 Files required for KEGG pathway reconstruction by MinPath (placed under data/)
    * KEGG-pathway.txt, KEGG-family.txt, KEGG-mapping.txt
    * These files are updated regularly

  3.1.2 File required for MetaCyc metabolic pathway reconstruction by MinPath (placed under data/)
    * MetaCyc-pathway.txt, MetaCyc-family.txt, MetaCyc-mapping.txt 
    * These files are updated regularly

  3.1.3 File required for SEED subsystem reconstruction by MinPath (placed under data/)
    * figfam_subsystem.dat (this file isn't updated since 2010)
    (note: current file is copied from figfam release6: 
     Release6/FigfamsData.Release.6/release_history/figfam_subsystem.dat; 
     please note it is an outdated version. )

  3.1.4 Any pathway system
   * what you need to prepare is a pathway-function mapping file similar to KEGG-mapping.txt or MetaCyc-mapping.txt
   * format of this file: one row for a mapping of a function to a pathway; pathway-id comes before functions, seperated by a tab;

 3.2 LP package for integer programming of the minimal pathway reconstruction
  package glpk-4.6 under glpk-4.6/ directory
  (you may need to re-install this package to use MinPath correctly, if you are NOT using it on a red hat linux machine;
   just go to glpk-4.6 directory, first call ./configure, then call make to re-install the package;
   make sure that you can call command glpsol under glpk-4.6/examples on your computer after you reinstall it)

4. Running MinPath
 4.1 Simply call MinPath.py to check its usage. Also check out the included examples as follows.

 4.2 Examples (under examples/ directory):
      under examples/
	#for KEGG pathway
	python ../MinPath.py -ko demo.ko -report demo.ko.minpath -details demo.ko.minpath.details
	#for MetaCyc pathway
	python ../MinPath.py -ec demo.ec -report demo.ec.minpath -details demo.ec.minpath.details
	#for SEED subsystem
	python ../MinPath.py -fig demo.fig -report demo.fig.minpath -details demo.fig.minpath.details
	#for any pathway
	python ../MinPath.py -any demo.any -map ../data/MetaCyc-mapping.txt -report demo.any.minpath -details demo.any.minpath.details

    Input of MinPath: ko (or fig, ec, or any) annotations (check demo.ko, demo.fig and demo.ec under examples/ directory)
    Format of the inputs: 
        one line per annotation;
        the first column is the protein/sequence id, and the second column is the annotation0;
        the two columns are seperated by a tab

    Outputs of MinPath: demo.ko.minpath, and demo.ko.minpath.details (and so on) 

 4.3 How to read the MinPath report file?
    e.g, demo.ko.minpath
    ...
    path 00030 kegg n/a  naive 1  minpath 1  fam0  42  fam-found  18  name  Pentose phosphate pathway
    path 00031 kegg n/a  naive 1  minpath 0  fam0  12  fam-found  2  name  Inositol metabolism
    ...
    1) path 00030, path 00031 are the KEGG pathway IDs (if your input file has fig families, the pathways will then be SEED subsystems)
    2) kegg n/a, indicates pathway reconstruction of your input dataset if not available (note: this information is only available for the genomes annotated 
	 in KEGG database); for input with fig families, KEGG is replaced by SEED 
    3) naive 1 or 0: the pathway is reconstructed, or not, by the naive mapping approach
    4) minpath 1 or 0: the pathway is kept, or removed by MinPath
    5) fam0: the total number of families involved in the corresponding pathway
    6) fam-found: the total number of involved families that are annotated
    7) name: the description of the corresponding pathway (subsystem)
  
 4.4 How to read the MinPath detailed report file?
    This report file lists all the pathways found MinPath, and a list of families each pathway includes
    e.g. demo.ko.minpath.details
    ...
    path 00010 fam0 56 fam-found 27 # Glycolysis / Gluconeogenesis
       K00001 hits 6 # E1.1.1.1, adh
       K00002 hits 1 # E1.1.1.4, adh
    ...
    1) path 00010 is the KEGG pathway ID, and fam0 and fam-found are the same as in 4.3
    2) this pathway includes families, K00001 and K00002, and so on 
       (here K numbers are used for KEGG families, and FIG ids for FIG families)
    3) family K00001 has 6 hits (i.e., 6 proteins/reads are annotated as this family)

5. MinPath on the web:
  http://omics.informatics.indiana.edu/MinPath
  (please check the project home page for updates and newer version of the MinPath)

6. Acknowledgments
Development was supported by NIH grant 1R01HG004908 to YY
