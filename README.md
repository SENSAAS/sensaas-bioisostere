# sensaas-bioisostere
A computational method for 3D shape-guided bioisosteric replacements and scaffold-hopping

[![badgepython](https://forthebadge.com/images/badges/made-with-python.svg)](https://www.python.org/downloads/release/python-370/)  [![forthebadge](https://forthebadge.com/images/badges/built-with-science.svg)](https://chemoinfo.ipmc.cnrs.fr/)

**SENSAAS** is a shape-based alignment method which allows to superimpose molecules in a rigid or flexible manner. It is based on the publications [SenSaaS: Shape-based Alignment by Registration of Colored Point-based Surfaces](https://onlinelibrary.wiley.com/doi/full/10.1002/minf.202000081) and [SENSAAS-Flex: a joint optimization approach for aligning 3D shapes and exploring the molecular conformation space](https://doi.org/10.1093/bioinformatics/btae105). 

**SENSAAS-Bioisostere** is an add-on of SENSAAS and SENSAAS-Flex that specifically replaces a query fragment by a new molecular fragment extracted from known bioactive molecules.


![example](/images/alignment.png)


**Documentation**: Full documentation is available at [https://github.com/SENSAAS/sensaas-bioisostere/docs/SENSAAS-Bioisostere-documentation.pdf](https://github.com/SENSAAS/sensaas-bioisostere/blob/main/docs/SENSAAS-Bioisostere-documentation.pdf)

**Website**: A web server to use SENSAAS, SENSAAS-Flex or SENSAAS-Bioisostere is available at https://chemoinfo.ipmc.cnrs.fr


## Main features

 - Allows bioisosteric replacements and scaffold-hopping

 - Bioisosteres are selected based on their 3D shape and pharmacophoric features similarities

 - Allows scaffold-hopping with fragments that have more than two substituents

 - Provides the superimposition of bioisosteres on the input molecular structure, which is kept as the reference coordinates

 - Provides similarity scores

 - Is free and open source


## Requirements

SENSAAS relies on the open-source library [Open3D](http://www.open3d.org). The current release of SENSAAS uses Open3D version 0.12.0 along with Python3.7 and the cheminformatics toolkit [RDKit](https://www.rdkit.org/docs/index.html).

Visit the following URL for using Python packages distributed via PyPI: [http://www.open3d.org/docs/release/getting_started.html](http://www.open3d.org/docs/release/getting_started.html) or conda: [https://anaconda.org/open3d-admin/open3d/files](https://anaconda.org/open3d-admin/open3d/files). For example, for windows-64, you can download *win-64/open3d-0.12.0-py37_0.tar.bz2*


## Virtual environment for python with conda (for Windows for example)

Install conda or Miniconda from [https://conda.io/miniconda.html](https://conda.io/miniconda.html)  
Launch Anaconda Prompt, then complete the installation:

	conda update conda
	conda create -n sensaas
	conda activate sensaas
	conda install python=3.7 numpy
	conda install perl

Once Open3D downloaded:
  
 	conda install open3d-0.12.0-py37_0.tar.bz2

Install RDKit (Open-Source Cheminformatics Software) that is compatible with Python 3.7 (eg: v2022.9.5 [look at older versions of RDKit](https://pypi.org/project/rdkit-pypi/))

	pip install rdkit-pypi
	

Retrieve and unzip sensaas-bioisostere repository in your desired folder. See below for running the programs **sensaas.py** or **meta-sensaas.py** for using the rigid version, **sensaasflex.py** or **meta-sensaasflex.py** for using the flexible version and **sensaas-bioisostere.pl** for scaffold-hopping. The directory containing executables is called sensaas-bioisostere-main.

## Linux

Scripts written in Perl are used in SENSAAS-Flex and SENSAAS-Bioisostere, thus, check for a Perl package on your system.

Install:

1. Python3.7 and numpy
2. Open3D version 0.12.0 (more information at [http://www.open3d.org/docs/release/getting_started.html](http://www.open3d.org/docs/release/getting_started.html))

The open-Source Cheminformatics Software RDKit must be installed (eg: version 2022.9.5). More information on RDKit can be found at [https://www.rdkit.org/docs/Install.html](https://www.rdkit.org/docs/Install.html) 

3. RDKit
  
Retrieve and unzip sensaas-bioisostere repository. The directory containing executables is called sensaas-bioisostere-main.

## MacOS

	Not tested


## Molecular viewer

We suggest you install a molecular viewer so you can visualize the molecular alignments.

PyMOL is a well-known software program for molecular visualisation available for Windows, Linux or macOS. More information at [PyMOL Wiki](https://pymolwiki.org) or [https://pymol.org/](https://pymol.org/)

Of note, we had trouble installing PyMOL in the current conda environment (python 3.7 and open3D 0.12.0). The installation method described in [sensaas-py](https://github.com/SENSAAS/sensaas-py/blob/main/) no longer works.


##  SENSAAS and SENSAAS-Flex

Readme, documentation and tutorials are available in the GitHub repository [sensaas-flex](https://github.com/SENSAAS/sensaas-flex/blob/main/) or in the documentation at [https://github.com/SENSAAS/sensaas-bioisostere/docs/SENSAAS-Bioisostere-documentation.pdf](https://github.com/SENSAAS/sensaas-bioisostere/blob/main/docs/SENSAAS-Bioisostere-documentation.pdf)


## SENSAAS-Bioisostere

**Example: Scaffold hopping of type 4 to optimize inhibitors of ATR kinase**

[Link to the publication](https://doi.org/10.1021/acs.jmedchem.4c00734) 

![example](/images/ATR-kinase.png)

Following files are in the folder examples:

The file 7jjg_B_VCD.sdf is the AZ20 molecule extracted from the publication. Coordinates of the ligand (VCD) were extracted from the co-structure PDB 7JJG [PDB entry](https://www.rcsb.org/structure/7JJG). It is the reference (or input) molecule.

The file list-atoms-7jjg_B_VCD contains the number of the atoms of the fragment to replace in molecule 7jjg_B_VCD.sdf.

The file AD1058.sdf is the synthesized bioisostere extracted from the publication and aligned on the AZ20 using SENSAAS-Flex.

**1. Create an SDF file containing the fragments that constitute your molecule**

**1a-** Firstly, use a molecular viewer to identify the number of atoms in the fragment that you want to replace. Then, write each number in a text file, separating them by a space. 
In our example, the file list-atoms-7jjg_B_VCD contains the following numbers: 

	2 3 4 5 22 23 18 20 21
	
Those numbers correspond to the highlighted atoms in the SDF file 7jjg_B_VCD.sdf as follows:

![example](/images/fragment-numbers.png)

**1b-** Secondly, use the following script to split the reference molecule into fragments. It tries to identify the fragment to replace by using the list of atoms from the file list-atoms-7jjg_B_VCD.

The syntax is:

	perl lea3d-MAKE_FGTS_LIB.pl <reference.sdf> <file with the list of atoms you want to replace> 

Example:

	perl lea3d-MAKE_FGTS_LIB.pl examples/7jjg_B_VCD.sdf examples/list-atoms-7jjg_B_VCD 

Here, the script creates 2 files : make_fgts_aggreg.sdf that contains 5 fragments and lib1.sdf that contains 4 fragments. 

Each SDF file contains the ensemble of fragments that can be combined to rebuild the reference molecule. This enables one fragment to be replaced by a bioisosteric fragment, after which the molecule can be rebuilt with all the other fragments remaining intact.

The output at your prompt may look like:

	lib1.sdf contains 4 fgts
	list of atoms to match: 2 3 4 5 22 23 18 20 21 in ../sensaas-bioisostere-main/examples/7jjg_B_VCD.sdf (fragmented in make_fgts_aggreg.sdf)
	2 fragments match the list of atoms
	Results:
	lib1.sdf fgt 2
	make_fgts_aggreg.sdf

The first line after the "Results:" tag indicates that the fragment number 2 in lib1.sdf seems to correspond to the list of selected atoms in text file list-atoms-7jjg_B_VCD. 
Visualize the file lib1.sdf by using a molecular viewer and check that the fragment number 2 is indeed the fragment you want to replace. For example, with PyMOL:

	pymol lib1.sdf

![example](/images/lib1-2.png)

Of note, you can choose another fragment to replace by selecting the lib name (lib1.sdf or make_fgts_aggreg.sdf) and by selecting the number of the fragment in the sdf file.

**2.Virtual screening of fragment libraries to identify bioisosteres for your selected fragment, followed by alignment of the bioisosteres on the input molecule**

The syntax is:

	perl sensaas-bioisostere.pl <reference.sdf> <file with fragments extracted from reference.sdf> <number of the fragment to replace> <library of fragment to screen>

In the db folder, 7 libraries of fragments are available. They were extracted from approved drugs or co-crystallized ligands from the PDBbind dataset. Combinations are substructures composed of 2 or 3 fragments. A test libary called test.sdf is also available in the folder db.

Table 1. Description of the 7 libraries of fragments extracted from bioactive molecules.

 * - **Library name**
   - **Number of fragments**
   - **Features of the fragments**
   - **Source**
 * - test.sdf
   - 9
   - 
   - 
 * - drug.sdf
   - 4509
   - single fragments
   - [e-Drug3D (2118 structures)](https://chemoinfo.ipmc.cnrs.fr/edrug3d)
 * - drug2.sdf
   - 2542
   - 2 combined fragments with MW ≤ 200
   - [e-Drug3D (2118 structures)](https://chemoinfo.ipmc.cnrs.fr/edrug3d)
 * - drug3.sdf
   - 2458
   - 3 combined fragments with MW ≤ 200
   - [e-Drug3D (2118 structures)](https://chemoinfo.ipmc.cnrs.fr/edrug3d)
 * - pdb.sdf
   - 7990
   - single fragments
   - PDBbind (12046 structures)
 * - pdb2.sdf
   - 8678
   - 2 combined fragments with MW ≤ 200
   - PDBbind (12046 structures)
 * - pdb3.sdf
   - 10333
   - 3 combined fragments with MW ≤ 200
   - PDBbind (12046 structures)

Example:

In our study, we will use a small library of 9 fragments called test (db/test.sdf):

	perl sensaas-bioisostere.pl examples/7jjg_B_VCD.sdf lib1.sdf 2 db/test.sdf

*Note 1: it takes about 1 hour to screen 500 fragments.*

*Note 2: it has been observed that the minimization with the standard process (RDKit, MMFF94 force field, 2000 Minimization cycles) is sometimes insufficient, resulting in a slight distortion of the structure. Finding a suitable trade-off was necessary here because this step is time-consuming.*

The output at your prompt may look like:

	perl scripts in ../sensaas-bioisostere-main
	SENSAAS-Bioisostere:
	Target ../sensaas-bioisostere-main/examples/7jjg_B_VCD.sdf
	4 fragments in lib1.sdf - referencex.sdf = libfgt_2.sdf contains the fragment to replace with 3 substituants
	Library: 6 fragments in ../sensaas-bioisostere-main/db/test.sdf having >= 3 substitution positions
	Screening using metasensaas (option -s mean -l 2)
	4 built bioisosteres
	Visualize the 4 aligned bioisosteres in ordered-bioisosteres-lib.sdf with PyMol (scores in ordered-bioisosteres-lib-score.txt) - Selected fragments from lib are in ordered-filtered-fragments-lib.sdf 		(scores in ordered-filtered-fragments-lib-score.txt)

The script creates 2 files: ordered-bioisosteres-lib.sdf that contains best aligned bioisosteres and ordered-bioisosteres-lib-score.txt that contains scores by decreasing order of similarity score.

**3. Visualization You can use any molecular viewer**

For instance, you can use PyMOL if installed.

	pymol 7jjg_B_VCD.sdf AD1058.sdf ordered-bioisosteres-lib.sdf

![example](/images/bioisostere.png)

Figure: one bioisostere possess a pyrrolo pyrimidine scaffold (in magenta). It is aligned on the input structure 7jjg_B_VCD.sdf (in green). The fitness score gfit+hfit = 1.54. It has the same molecular structure as AD1058.sdf.

Here, we show that SENSAAS-bioisostere found a pyrrolo pyrimidine fragment to replace the pyrimidine + cyclopropane fragment in AZ20 (7jjg_B_VCD.sdf) thus, reproducing the scaffold-hopping step described in the publication. When screened against larger libraries (see libraries in folder db), several other fragment bioisosteres can be identified. Those could also be interesting alternatives to the pyrimidine core.


## Licenses
SENSAAS code is released under [the 3-Clause BSD License](https://opensource.org/licenses/BSD-3-Clause)

## Copyright
Copyright (c) 2018-2021, CNRS, Inserm, Université Côte d'Azur, Dominique Douguet and Frédéric Payan, All rights reserved.

## Reference
[Biyuzan et al., SENSAAS-Flex: a joint optimization approach for aligning 3D shapes and exploring the molecular conformation space, *Bioinformatics*, **2024**, 40, 3](https://doi.org/10.1093/bioinformatics/btae105)

Bibtex format :

	@article{10.1093/bioinformatics/btae105,
    author = {Biyuzan, Hamza and Masrour, Mohamed-Akram and Grandmougin, Lucas and Payan, Frédéric and Douguet, Dominique},
    title = {SENSAAS-Flex: a joint optimization approach for aligning 3D shapes and exploring the molecular conformation space},
    journal = {Bioinformatics},
    volume = {40},
    number = {3},
    pages = {btae105},
    year = {2024},
    month = {02},
    doi = {10.1093/bioinformatics/btae105},
    url = {https://doi.org/10.1093/bioinformatics/btae105},
    eprint = {https://academic.oup.com/bioinformatics/article-pdf/40/3/btae105/56886864/btae105.pdf},
	}

[Douguet D. and Payan F., SenSaaS: Shape-based Alignment by Registration of Colored Point-based Surfaces, *Molecular Informatics*, **2020**, 8, 2000081](https://onlinelibrary.wiley.com/doi/full/10.1002/minf.202000081)
   
Bibtex format :

	@article{10.1002/minf.202000081,
	author 		= {Douguet, Dominique and Payan, Frédéric},
	title 		= {sensaas: Shape-based Alignment by Registration of Colored Point-based Surfaces},
	journal 	= {Molecular Informatics},
	volume 		= {39},
	number 		= {8},
	pages 		= {2000081},
	keywords 	= {Shape-based alignment, molecular surfaces, point clouds, registration, molecular similarity},
	doi 		= {https://doi.org/10.1002/minf.202000081},
	url 		= {https://onlinelibrary.wiley.com/doi/abs/10.1002/minf.202000081},
	eprint 		= {https://onlinelibrary.wiley.com/doi/pdf/10.1002/minf.202000081},
	year 		= {2020}
	}

