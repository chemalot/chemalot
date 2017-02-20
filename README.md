Overview of chemalot
===================================

Cheminformatics, like other scientific computing disciplines, requires a wide array of software tools, each with different data formats and program interfaces.  Frequently, we use these tools to compute properties and geometries of small molecules, pass the computed data to additional, multi-step, data-analyses, and finally build statistical models.  These chemistry data workflows commonly use upwards of twenty different programs in sequence and are run repeatedly and often automatically.  Thus software that connects and standardizes these disparate programs in a way that is amenable to large clusters or cloud computing is essential.  Here, we present chemalot, an open-source collection of command-line programs which both wraps and extends existing cheminformatics programs.  Each chemalot program reads and writes a single chemistry data format (SDF file) in a way that takes advantage of a standard, UNIX parallelization method called pipelines.  As a result, a cheminformatics scientist can quickly develop a chemistry workflow that performs a series of 2D and 3D computations on a library of thousands of molecules.  The workflow can take advantage of multiple cores and does not require any external "pipelining" software.

The package provides a comprehensive collection of tools for processing files with chemical structures and associated data. At Genentech these tools have been used for tasks such as:

* Chemical registration [[1](#1)]
* Post processing of docking results [[2](#2)]
* Clustering and prioritization of high throughput screening results [[3](#3)]
* Selection and analysis of compounds in HTS libraries [[4](#4)]
* Creation of QSAR models [[5](#5)]
* Analysis of strain in conformation of small molecules in crystal structures [[6](#6)]
* Compilation of project data files for SAR analysis [[7]](#7)

The package complements command line tools released in other packages, notably [Aestel](https://sourceforge.net/projects/aestel/) and [autocorrelator](https://github.com/chemalot/autocorrelator). The tools integrate well with command line tools from other open source projects such as [Open Babel](http://openbabel.org) and [RDKit](http://www.rdkit.org/) as well as with command line tools from commercial software such as [OpenEye](http://www.eyesopen.com/) and [CCG](https://www.chemcomp.com/).

For graphical programming and debugging of UNIX command tools that process SDF files on UNIX pipes the [chemalot_knime](https://github.com/chemalot/chemalot_knime/) package can be used.

This software is released under the "Apache License Version 2.0".

 Summary of dependencies:
-------------------------------

Most of the chemalot command line programs are implemented in java. Java 1.8 is thus required to run most programs.

chemalot uses command line programs and libraries from several Open Source and commercial 
software packages. For convenience the chemalot package includes the following libraries and executable:

  * Aestel_20160130.jar
  * autocorrelator.jar
  * commons-cli-1.2.jar
  * commons-pool-1.3.jar
  * cxf-2.2.9.jar
  * groovy-all-1.8.8.jar
  * httpcore-4.4.3.jar
  * jdom.jar
  * testng-6.2.jar
  * Open Babel

Not included are commercial libraries and executables for which licenses are required from the corresponding software vendor:

   * oejava-2016.Feb.1-Linux-x64.jar __(required)__
   * quacpac
     http://www.eyesopen.com/quacpac
   * szybki
     http://www.eyesopen.com/szybki
   * gaussian
     http://www.gaussian.com/
   * bmin (Macromodel from Schrodinger, LLC)
     https://www.schrodinger.com/
   * moebatch (MOE package from Chemical Computing Group)
     https://www.chemcomp.com/

Also not included is

   * R (The R Foundation for Statistical Computing)
     https://www.r-project.org/foundation/


 Licensing
-----------------------
Copies of the licenses for included packages can be found in the license subdirectory.

autocorrelator.jar
> The autocorrelator package is a separate open source package.
  The autocorrelator.jar is released under the
  [GNU General Public License v3](http://www.gnu.org/copyleft/gpl.html).
  The classes in autocorrelator.jar are not linked to any code in the
  chemalot package.
  More information including the autocorrelator is available at:
  [https://github.com/chemalot/autocorrelator](https://github.com/chemalot/autocorrelator)

OEChem Java library:
> *THE OPENEYE OECHEM JAR FILE IS NOT INCLUDED IN THIS PACKAGE*
  The oechem Java library is required for most of the programs in this package.
  To obtain the library and a license please contact:
  [OpenEye Scientific Software](http://www.eyesopen.com/contact-us)

Open Babel
> The Open Babel package is a separate open source package.
  Open Babel is licensed under the
  [GNU General Public License, version 2](http://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html).
  Only the babel executable is needed (and included) in the chemalot package.
  More information on Open Babel is available at:
  [http://openbabel.org](http://openbabel.org)

quacpac and szybki:
> *QUACPAC AND SZYBKI ARE NOT INCLUDED IN THIS PACKAGE*
  They are distributed by OpenEye Scientific Software Inc.
   For more information visit:
   [http://www.eyesopen.com/quacpac](http://www.eyesopen.com/quacpac)
   [http://www.eyesopen.com/szybki](http://www.eyesopen.com/szybki)

gaussian:
> *GAUSSIAN IS NOT INCLUDED IN THIS PACKAGE*
  Gaussain is a program distributed by Gaussian, Inc.
  For more information visit:
  [http://www.gaussian.com/](http://www.gaussian.com/)

moebatch (MOE):
> *MOE IS NOT INCLUDED IN THIS PACKAGE*
  MOE is a program distributed by Chemical Computing Group
  For more information visit:
  [https://www.chemcomp.com/](https://www.chemcomp.com/)

bmin (Macromodel):
> *MACROMODEL IS NOT INCLUDED IN THIS PACKAGE!*
  Macromodel is a commercial program distributed by Schrodinger, LLC
  For more information visit:
  [http://www.schrodinger.com/](http://www.schrodinger.com/)

R (The R Foundation for Statistical Computing):
> *R IS NOT INCLUDED IN THIS PACKAGE!*
  R package is a separate open source package and is licensed under the
  [GNU General Public License, version 2](http://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html).
  More information on R is available at:
  [https://www.r-project.org/foundation/](https://www.r-project.org/foundation/)

All other files
> are released under the
[Apache License Version 2.0](http://www.apache.org/licenses/LICENSE-2.0.html)


 Installation on LINUX
----------------------------
* Unpack the chemalot_XXX.tgz file

* Contact [OpenEye](http://www.eyesopen.com/contact-us)
  for a copy of their oechem Java toolkit.

   * copy the jar file from OpenEye into `chemalot/lib` directory.
     The examples were tested using `oejava-2014.Feb.3-Linux-x64.jar`

   * set the OE_LICENSE environment variable to point to your
     OpenEye license file. e.g.
     `setenv OE_LICENSE ~oechem/etc/oe_license.txt`

* Add the chemalot/bin path to your UNIX path, e.g.
  `set path=($path ~/chemalot/bin)`

* Define the AESTEL_DIR environment variable pointing to this directory e.g.
  `setenv AESTEL_DIR ~/chemalot`

* Define the OBABEL environment variable pointing to the Open Babel executable.
  This is necessary because a few of the command line wrappers use Open Babel
  to convert to and from files in gaussian format and we want to avoid confusion
  about which path to use.
  `setenv OBABEL ~/chemalot/bin/obabel`


 Compiling
------------------
Note: there is no need to compile the chemalot package as it already include a compiled
jar file. This description is for developers wjo whish to make changes to the code.

In order to compile this package you must have the **oechem Java
toolkit** installed as described under "Installation on LINUX".
In addition you must have [ant](http://ant.apache.org/)
installed.
Once these two conditions are met, simply type `ant` in the root directory
of the chemalot package and the chemalot jar file will be compiled.


 Creating JavaDoc
------------------
Issue `ant javaDoc`. The documentation will be created in the doc sub directory.


 Testing
-----------------
Requirements

* Python 2.6 or greater
* Libraries mentioned above.  Some tests will fail if you do not have all the dependencies above.

* Regression  tests:
   Regression tests can be executed by issuing `ant test`
   At this time important low level methods have regression tests but many
   higher level methods do not.

* Functional tests:
   A minimal functional test can be performed by issuing:
   `test/simple.csh`
   from the root directory of the package.
   This will execute the core command line programs and compare the resulting
   output to a reference file.

   Most command line programs have specific functional tests. To run all tests, issue:
   `cd test; runTests.py`


 Usage
--------------------
chemalot contains a set of command line wrapper scripts located in the bin subdirectory.
Each of the commands will print an explanation of the command and the available options
when executed without any command line options.

An example demonstrating the clustering process described in reference [[3](3)] can be
found in:
[examples/NovartisMalariaBox/=readme.md](examples/NovartisMalariaBox/=readme.md)

Also see the test directory which demonstrates simple executions of most commands.


 Other Operating Systems
------------------------
The Java code should run on other environments than LINUX as long as
the OpenEye library is available. However this has not been tested and
modified wrappers might be required.
The authors do not currently have experience with other environments
but would welcome reports describing the setup. For the time being, this
package should be considered LINUX only.


Brief Description of Command Line Programs
------------------------------------------

* **`g2sdf.pl`**
   Create an SDF file from a Gaussian output file.

* **`OEProps.csh`**
   Compute atom, bond, ring, and other count properties, TPSA, and other 2D properties. (see help text)

* **`QTorsionProfileGenerator.csh`**
   Perform a Gaussian torsion scan from an input molecule with a specified rotatable bond (requires Gaussian license)
   _Additional requirement: gaussian_

* **`sdf2DAlign.csh`**
   Transforms the 2D coordinates of input molecules according to the matching substructures specified in the template SDF file.

* **`sdf2g.pl`**
   Create Gaussian input files from molecules in an input SDF file.

* **`sdfAggregator.csh`**
   Given a set of input molecules with SDF tag data, group them by a specified tag value and then perform a grouping function (e.g. find the average "My Assay IC50" (grouping function) for each "Chemical Series" (the group-by SDF tag).

* **`sdfAlign.pl`**
   Transforms the 3D coordinates of input molecule conformers to a reference ligand and calculate the RMSD.

* **`sdfALogP.csh`**
   Calculates the ALogP and atom type counts.
   [Ghose AK, Crippen GM. J. Atomic Physicochemical Parameters for Three-Dimensional Structure-Directed Quantitative Structure-Activity Relations I. Partition Coefficients as a Measure of Hydrophobicity. Comput. Chem. 1986, 7 (4), 565-577; Viswanadhan VN, Reddy MR, Bacquet RJ, Erion MD. Assessment of methods used for predicting lipophilicity: Application to nucleosides and nucleoside bases. J. Comput. Chem. 1993, 14 (9), 1019-1026; Ghose AK, Viswanadhan VN, Wendoloski JJ. Prediction of Hydrophobic (Lipophilic) Properties of Small Organic Molecules Using Fragmental Methods: An Analysis of ALOGP and CLOGP Methods. J. Phys. Chem. A 1998, 102 (21), 3762-3772]

* **`sdfBinning.csh`**
   Groups numerical values into bins

* **`sdfCalcProps.csh`**
   Serves as the "properties" warehouse to which any calculator command line programs can be added, thus, enable a single point access to the properties.
   Feng JA, Aliagas I, Bergeron P, Blaney JM, Bradley EK, Koehler MFT, Lee M, Ortwine DF, Tsui V, Wu J, Gobbi A. An Integrated Suite of Modeling Tools That Empower Scientists in Structure- and Property-Based Drug Design. Journal of Computer-Aided Molecular Design 2015, 29 (6), 511-523.

* **`sdfCatsIndexer.csh`**
   Generate CATS fingerprints for input molecules.
   [Schneider G, Neidhart W, Giller T, Schmid G. "Scaffold-Hopping" by topological pharmacophore search: a contribution to virtual screening. Angew. Chem. Int. Ed. 1999, 38 (19), 2894-2896]

* **`sdfCFP.csh`**
   Generate circular fingerprints for input molecules.
   [Rogers D, Hahn M. Extended-Connectivity Fingerprints. J. Chem. Inf. Model., 2010, 50 (5), 742.754]

* **`sdfCluster.pl`**
   Clusters molecules using Atom-Atom-Path similarity and Sphere Exclusion algorithm.
   [Gobbi A, Giannetti AM, Chen H, Lee M. Atom-Atom-Path similarity and Sphere Exclusion clustering: tools for prioritizing fragment hits. J. Cheminform., 2015, 7:11]

* **`sdfConformerSampler.csh`**
   Generates conformers combinatorially modifying torsional angles as defined in a torsion file. Torsions of OH and HN2 groups are evaluated automatically rotated.

* **`sdfEnumerator.csh`**
   Enumerates a combinatorial library given the specified SMIRKS and corresponding reagent input files.

* **`sdfEStateCalculator.csh`**
   Compute the occurrence (count) of each E-state atom group in the input molecule and the corresponding sums as well as E-state indices of the individual atoms in a molecule.
   [Hall LH, Kier LB. Electrotopological State Indices for Atom Types: A Novel Combination of Electronic, Topological, and Valence State Information. J. Chem. Inf. Comput. Sci. 1995, 35 (6), 1039-1045]

* **`sdfFilter.csh`**
   Remove molecules from the SDF file based on the heavy atom count, number of components, invalid atoms, and max atomic number.

* **`sdfFingerprinter.csh`**
   Generate various types of fingerprints including, linear fingerprints and smarts based fragment fingerprints.

* **`sdfFPCluster.pl`**
   Use the specified fingerprints to cluster input molecules using the Sphere Exclusion clustering algorithm.  A radius of 0.5 is a good values for clustering HTS libraries.
   [Gobbi A, Lee M. DISE: Directed Sphere Exclusion. J. Chem. Inf. Comput. Sci. 2002, 43 (1), 317.323]

* **`sdfFPNNFinder.csh`**
   Use the specified fingerprints to identify most similar molecules (nearest neighbors) for each molecules in the input file based on their Tanimoto similarities. It can also be used to compute activity cliffs.

* **`sdfFPSphereExclusion.csh`**
   Use the specified fingerprint to compile a diverse sub set using the Sphere Exclusion algorithm.
   [Gobbi A., Lee M. DISE:. Directed Sphere Exclusion. J. Chem. Inf. Comput. Sci. 2002, 43 (1), 317.323]

* **`sdfGrep.pl`**
   Remove a molecule from the SDF file if the field of interest does not matching the specified requirement

* **`sdfLE.grvy`**
   Calculate various ligand efficiencies, i.e. LE, LLE
   [Hopkins AL, Groom CR, Alex A. Ligand efficiency: a useful metric for lead selection. 2004, 9 (10), 430-431; Leeson PD, Springthorpe B. The influence of drug-like concepts on decision-making in medicinal chemistry. Nat. Rev. Drug Disc. 2007, 6 (11), 881-890]

* **`sdfMACCSKeys.csh`**
   Generate MACCS keys or counts for input compounds [based on [ChemAxon](https://www.chemaxon.com/forum/ftopic8138.html) and [RDKit](http://rdkit.org/Python_Docs/rdkit.Chem.MACCSkeys-module.html) implementations]

* **`sdfMCSSNNFinder.csh`**
   User the Atom-Atom-Path similarity to identify the nearest neighbors for the input compounds. It can be used to compute activity cliffs.
   [Gobbi A, Giannetti AM, Chen H, Lee M. Atom-Atom-Path similarity and Sphere Exclusion clustering: tools for prioritizing fragment hits. J. Cheminform. (2015) 7:11]

* **`sdfMCSSSphereExclusion.csh`**
   Use MCSS or Atom-Atom-Path similarity to compile a diverse sub set using the Sphere Exclusion algorithm.
   [Gobbi A, Giannetti AM, Chen H, Lee M. Atom-Atom-Path similarity and Sphere Exclusion clustering: tools for prioritizing fragment hits. J. Cheminform. (2015) 7:11]

* **`sdfMDLSSSMatcher.csh`**
   Remove molecules from the SDF file that don't match any of substructures in MDL query file

* **`sdfModelCreateValidate.pl`**
   Use sdfRModelPredictor.pl to create a Machine Learning Model and validate at the same time using randomly selected training and test sets.
   Uses sdfR???ModelCreator.pl in the background.

* **`sdfMMConfAnalysis.pl`**
   Perform strain energy analysis of input conformers including generation and geometry optimization of a large number of conformations. Also evaluates the energy of the minimized input conformation with several restraint strengths.
   _Additional requirements: bmin, moebatch, szybki__

* **`sdfMMMinimize.csh`**
   Perform a geometry optimization using a molecular mechanics force field, with wrapped choices of Macromodel (Schrodinger), MOE (CCG) or SZYBKI (OpenEye). (License requirements)
   _Additional requirements: bmin, moebatch, szybki_

* **`sdfMolSeparator.csh`**
   Separate the disconnected molecules in a molfile (e.g. salt and compound) and output them in individual records.

* **`sdfMultiplexer.pl`**
   Parallelize the execution of command line scripts by executing multiple instances of a command line string and distributing the input molecules to the various instances. The output is combined back into a single file.

* **`sdfNormalizer.csh`**
   Normalize molecules according to Genentech's business rules. Unique tautomers are generated by Quacpac from OpenEye.
   [Gobbi A, Lee M. Handling of Tautomerism and Stereochemistry in Compound Registration. J. Chem. Inf. Model. 2011, 52 (2), 285.292]
   _Additional requirement: quacpac_

* **`sdfRExecutor.pl`**
   Apply costum R scripts to data in SD file and add computed fields to the output file.

* **`sdfRGroupCalcProps.pl`**
   Calculate the properties of molecule fragments with attachment points (e.g. [U+1], [U+2]); companion program to sdfRGroupExtractor.pl

* **`sdfRGroupExtractor.pl`**
   Fragment input molecules according to the specified transformations in SMIRKS or SMARTS format into R-groups with charged Uranium atoms representing attachment points

* **`sdfRingSystemExtraction.csh`**
   Fragments input molecules and outputs the largest (linked) ring system as well as the set of the basic rings as SMILES.

* **`sdfRModelPredictor.pl`**
   Compute the prediction according to the specified model created by sdfRRandomForestCreator.pl or sdfRSVMCreator.pl
   _Additional requirement: R_

* **`sdfRMSDNNFinder.csh`**
   Calculates RMSD values between conformers of molecules. This can align conformers by minimizing the RMSD

* **`sdfRMSDSphereExclusion.csh`**
   Applies Sphere Exclusion algorithm to find centroids of conformer clusters based on a given RMSD radius

* **`sdfRRandomForestCreator.pl`**
   Create models (R sessions) using Random Forest algorithm; companion program to sdfRModelPredictor.pl
   _Additional requirement: R_

* **`sdfRSVMCreator.pl`**
   Create models (R sessions) using Support Vector Machine algorithm; companion program to sdfRModelPredictor.pl
   _Additional requirement: R_

* **`sdfSelectivityCalculator.csh`**
   Computes selectivity (ratio) based on the specified numerator and denominator fields considering operator values

* **`sdfSliceByRe.pl`**
   Partition sdf files by ranges of rows

* **`sdfStructureTagger.csh`**
   Tag molecules with specified names based on the corresponding SMARTS or molfile (queries)

* **`sdfSubRMSD.csh`**
   Calculate the RMSD between a supplied fragment (e.g. core) and the matching part of the input molecule. Molecules need to be pre-aligned.

* **`sdfTopologicalIndexer.csh`**
   Compute topological indices, i.e. Balaban, Wiener, and Zagreb

* **`sdfTorsionScanner.csh`**
   Given one or more molecules (e.g. sdf), generate a set of conformers rotated around a single rotatable bond within the input molecules. Useful as pre-step to sdfMMMinimize.csh to calculate energy torsion scans. (Optional minimization requires a SZYBKI license).
   _Additional requirement: szybki_

* **`tabTagTool.pl`**
   Modify column header and filter tab-delimited files

References
--------------

<a name='1'>[1]</a> Gobbi A, Lee M. Handling of Tautomerism and Stereochemistry in Compound Registration. Journal of Chemical Information and Modeling 2011, 52 (2), 285-292.

<a name='2'>[2]</a> Feng JA, Aliagas I, Bergeron P, Blaney JM, Bradley EK, Koehler MFT, Lee M, Ortwine DF, Tsui V, Wu J, Gobbi A. An Integrated Suite of Modeling Tools That Empower Scientists in Structure- and Property-Based Drug Design. Journal of Computer-Aided Molecular Design 2015, 29 (6), 511-523.

<a name='3'>[3]</a> Gobbi A, Giannetti AM, Chen H, Lee M. Atom-Atom-Path Similarity and Sphere Exclusion Clustering: Tools for Prioritizing Fragment Hits. Journal of Cheminformatics 2015, 7:11.

<a name='4'>[4]</a> Beresini MH, Liu Y, Dawes TD, Clark KR, Orren L, Schmidt S, Turincio R, Jones SW, Rodriguez RA, Thana P, Hascall D, Gross DP, Skelton NJ. Small-Molecule Library Subset Screening as an Aid for Acceleration Lead Identification. Journal of Biomolecular Screening 2014, 19 (5), 758-770.

<a name='5'>[5]</a> Aliagas I, Gobbi A, Heffron T, Lee M, Ortwine DF, Zak M, Khojasteh SC. A Probabilistic Method to Report Predictions from a Human Liver Microsomes Stability QSAR Model: A Practical Tool for Drug Discovery. Journal of Computer-Aided Molecular Design 2015, 29 (4), 327-338.

<a name='6'>[6]</a> [50 years of the Cambridge Structural Database](https://www.ccdc.cam.ac.uk/News/csd50/CSD50_program_agenda_singles.pdf)

<a name='7'>[7]</a> Lee M, Aliagas IM, Feng JA, Gabriel T, O'Donnell TJ, Sellers B, Wiswedel B, Gobbi A. Command line programs as workflow tools for drug discovery. In preparation.
