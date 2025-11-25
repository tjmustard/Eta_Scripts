THIS CODE WAS INITIALLY CREATED BETWEEN 2010-2015 WHILE I WAS ATTENDING GRADUATE SCHOOL AT OREGON STATE UNIVERSITY. A VERSION OF THIS CODE CAN BE FOUND IN MY THESIS HERE: <https://ir.library.oregonstate.edu/concern/graduate_thesis_or_dissertations/xp68kn058>

Eta_Scripts
======

## Why Eta_Scripts?

These scripts were originally written to assist in creating and reading computational job files.

Credits:

These scripts were written mostly by myself, but I had some help and they deserve recognition:

* O. Maduka Ogba helped with numpy integration and Dist/Angle/Dihedral calculations.
* Prof. Paul Ha-Yeon Cheong was the inspiration for this project. While in grad school Paul wrote several scripts to assist his research. Most of this project is inspired by those scripts.
* Joshua J. Kincaid was instrumental in making the z-matrix to xyz conversion possible.
* The entire PHYC group for their inspiration and support. (Ryne C. Johston, O. Maduka Ogba, Wojtek Rajski, Ommi Pattawong, Lindsay Wills, Daniel Walden, Kevin Snyder, Jacob Buchanan, Ben Hanken, I-Ya Chang)

You may acknowledge the use of the scripts/programs in published material as:
Thomas J. L. Mustard, O. Maduka Ogba, Paul Ha-Yeon Cheong. unpublished results.

* [Eta_Scripts](http://github.com/tjmustard/Eta_Scripts)
* [PHYC Group](http://phyc.chem.oregonstate.edu/)

This repository started in ~2014, all the original code is in my [thesis](https://ir.library.oregonstate.edu/concern/graduate_thesis_or_dissertations/xp68kn058) and some updates were made over the years. The repository was deleted/taken down and then reuploaded with the same name. Recently an update to python 3 was made.

## How is everything organized?

Folder Setup:

    Eta_Scripts/
     |
      ——> EtaLib/      #Contains the majority of the code functions (core library modules)
          Psientific/  #Computational chemistry utility programs

## Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/tjmustard/Eta_Scripts.git
   cd Eta_Scripts
   ```

2. Install the package (recommended to use a virtual environment):

   ```bash
   pip install -e .
   ```

## Requirements

* Python 3.x
  * The codebase has been updated to Python 3
* NumPy
  * Required for computational operations
* See `requirements.txt` for the complete list of dependencies

# Psientific Code

This code was written starting in 2012 during my tenure as a graduate student at Oregon State University and continuing to present.

### nMAP.py

This will map one molecule onto another. This is very helpful for accelerated QM optimization, transition state search and ligand design. The idea is there are two set of molecules: Parents and Children. The Parents are the majority of the molecule and the children are the minor changes you are making. But it does not have to be this way, it effects mostly the output file naming.

To properly run a nMAP.py mapping you need:

* A Parent structure folder filled with Parent structures
* A Children structure folder filled with Children structures
* A list of atoms to align each set of structures.
  * One for both the Parents and Children
* A list of atoms to delete from each set of molecules
  * Again one for both the Parents and Children

To quickly create a template input file run:

    # nMAP.py -t
    #Input file
    parentFolder =
    #Output folder
    output =
    #Base name for files
    name =
    #Library folder to use
    childFolder =
    #Quatfit location
    quatfit = /opt/apps/Quatfit/quatfit
    #List of atoms to align to
    parent =
    child =
    #List of atoms to delete from the input (parent) file
    deleteParent =
    deleteChild =

A template folder is located in the Psientific folder to assist in understanding the program.

### nDist-Int-Splitter.py

This will split one molecule onto all the combinations possible. This is very helpful for generating input geometries for distortion-interaction calculations. The idea is there are several sections to a catalytic system. These are broken into seperate lines of atom numbers. Then the script will make all the combinations possible.

To properly run a nDist-Int-Splitter.py you need:

* A structure file
* An output folder location of interest (nDist-Int-Splitter.py will make one if it does not exist)
* A list of atoms for each section.
  * One line for each section
  * You can use "-" to denote all atoms from a to b

To quickly create a template input file run:

    # nDist-Int-Splitter.py -t
    #Input file name
    input = COH-Ph.xyz
    #Output folder and name
    output = Dist-Int
    #Molecule sections
    $SECTIONS       
    1-6
    7-11
    12 13 14
    $ENDSECTIONS

A template folder is located in the Psientific folder to assist in understanding the program.

### nDeprotonate.py

This will deprotonate one molecule onto all the combinations possible. This is very helpful for quickly generating input geometries for QM calculations. Then the script will make all the combinations possible for the number of deprotonations your interested in.

To properly run a nDeprotonate.py you need:

* A structure file
* An output folder location of interest (nDeprotonate.py will make one if it does not exist)
* The number of deprotonations you wish to run.
* A list of equivalent hydrogens.
  * One line for each set of hydrogens

To quickly create a template input file run:

    # nDeprotonate.py -t
    #Input file name
    input = COH-Ph.xyz
    #Output folder and name
    output = Deprotonate
    # Number of deprotonations
    deprotonations = 2
    #Protons to delete
    $PROTONS
    7 11
    8 10
    9
    14
    $ENDPROTONS

A template folder is located in the Psientific folder to assist in understanding the program.

### nPersistentOTS.py

This software wraps around your QM code to add functionality and accelerate OPTs, SCANs, and TS-SCANs. TS-SCANs are dynamic SCANs that will continue to search up a positionally restrained optimization scan to locate a input geometry for further transition state search.

To properly run nPersistentOTS you will need:

1. OPT acceleration:

* Run nPersistentOTS in OPT mode.
* Set your maximum number of steps to ~15-30. This allows nPersistentOTS to check periodically for "jitter" in the energy and displacement.

2. SCAN function:

* Run nPersistentOTS in SCAN mode.
* Have a correctly built and names Scan-.....txt file with SCAN turned on. (example below)
* If you want to pair the OPT acceleration with the SCAN function lower your maximum number of steps as seen in the OPT section above.

3. TS-SCAN function:

* Run nPersistentOTS in SCAN or TS-SCAN mode.
* Have a correctly built and names Scan-.....txt file with both SCAN and TS turned on. (example below)
* If you want to pair the OPT acceleration with the SCAN function lower your maximum number of steps as seen in the OPT section above.

##### (TS-)Scan Configuration File Example

####### Scan-COPE.txt
    #VARIABLE=bond 2 6
    #SCAN=yes
    #TS=dyn
    #SCANLENGTH=30
    #SCANSTEPNOW=0
    #STEPSIZE=-0.1
    #WAIT=0
    #MINTSSCANSIZE=none
    #MINORTS=one

**This file contains several pieces of data:**
**NOTE:** All variables can be either lower or upper case.

***
**\#VARIABLE**

The variable you will be SCANing through. It can be a "bond", "angle", or "dihedral" and the atoms associated with that variable.

Variable|Option|Definition
---------|:------------------|:---------------------------------------------------------------------

# VARIABLE=|bond 35 6|Will scan along the bond between atoms 35 and 6

|angle 3 5 2|Will scan along the angle between 3-5-2 with 5 being the middle atom.
|dihedral 6 7 3 4|Will scan along the dihedral 6-7-3-4.

***
**\#SCAN**

If you wish to do a SCAN

Variable|Option|Definition
---------|:------------------|:---------------------------------------------------------------------

# SCAN=|yes|SCAN function is turned on

|no (or anything else)|SCAN function is tuned off.

***
**\#TS**

If you wish to implement the TS-SCAN logic

Variable|Option|Definition
---------|:------------------|:---------------------------------------------------------------------

# TS=|yes|TS-SCAN logic is turned on

|dyn|TS-SCAN logic will be turned on when the energy starts to go uphill. Useful for starting scan from far out.
|no|TS-SCAN logic is off.

***
**\#SCANLENGTH**

How many maximum steps you wish to take.

Variable|Option|Definition
---------|:------------------|:---------------------------------------------------------------------

# SCANLENGTH=|30|The scan will run for a maximum of 30 steps

***
**\#SCANSTEPNOW**

What step numper you are currently on. You can edit this number to prematurely kill a SCAN.

Variable|Option|Definition
---------|:------------------|:---------------------------------------------------------------------

# SCANSTEPNOW=|0|This is the current step number you are on. This is editable

***
**\#STEPSIZE**

How much the variable will change between SCAN steps.

Variable|Option|Definition
---------|:------------------|:---------------------------------------------------------------------

# STEPSIZE=|-0.1|This number will be in Angstrom for bonds and Degrees for angle and dihedral

***
**\#WAIT**

If you want the TS-SCAN logic to wait a number of steps before initializing.

Variable|Option|Definition
---------|:------------------|:---------------------------------------------------------------------

# WAIT=|0|Wait X steps before turning on the TS-SCAN logic

***
**\#MINTSSCANSIZE**

What the minimum step size should be for a TS-SCAN. The minimum step size is determined from a simple calculation based on the mass of the atoms of interest. This is based on the fact that heavier atoms need a more accurate input geometry to accurately compute the TS.

Bond:
`#MINTSSCANSIZE = x/100.0/Mass of atom 1 + Mass of atom 2`

Angle:
`#MINTSSCANSIZE = x/Mass of atom 1 + Mass of atom 3`

Dihedral:
`#MINTSSCANSIZE = x/100.0/Mass of atom 1 + Mass of atom 4`

Variable|Option|Definition
---------|:------------------|:---------------------------------------------------------------------

# MINTSSCANSIZE=|none|x = 50 (C-C bond: 0.021Å)

|loose|x = 100 (C-C bond: 0.042Å)
|tight|x = 25 (C-C bond: 0.010Å)
|verytight|x = 5 (C-C bond: 0.002Å)

***
**\#MINORTS**

If we should explore minor TSs.

Variable|Option|Definition
---------|:------------------|:---------------------------------------------------------------------

# MINORTS=|one|Allow for one minor TS to be skipped. Useful if a geometric change is necessary for the TS to occur

|yes|Allow for an infinite number of minor TSs to be skipped.
|no|Do not allow a single minor TS to be skipped. Once the energy goes down nPersistenceOTS will go back to the last geometry and take a smaller step size.

I hope you like the scripts!

Thomas J L Mustard, Ph.D.
