###WARNING:
Several updates may have broken some scripts/libraries. Please use the 0.4a branch until this notice is changed.

Eta_Scripts
======

## Why Eta_Scripts?

These scripts were written to assist in creating and reading computational job files. Most, if not all, of the scripts are written with Python. This project is Gaussian heavy due my extensive use of it. But I do use other softwares and building upon the Gaussian scripts several other programs are supported. 

Credits:

These scripts were written mostly by myself, but I had some help and they deserve recognition:

* O. Maduka Ogba helped with numpy integration and Dist/Angle/Dihedral calculations.
* Prof. Paul Ha-Yeon Cheong was the inspiration for this project. While in grad school Paul wrote several scripts to assist his research. Most of this project is inspired by those scripts.
* Joshua J. Kincaid was instrumental in making the z-matrix to xyz conversion possible. I had, and still have, no idea how this math works, and Joshua explained the proccess perfectly.
* The entire PHYC group for their inspiration and support. (Ryne C. Johston, O. Maduka Ogba, Wojtek Rajski, Ommi Pattawong, Lindsay Wills, Daniel Walden, Kevin Snyder, Jacob Buchanan, Ben Hanken, I-Ya Chang)

You may acknowledge the use of the scripts/programs in published material as:
Thomas J. L. Mustard, O. Maduka Ogba, Paul Ha-Yeon Cheong. unpublished results.

* [Eta_Scripts](http://tjmustard.github.io/Eta_Scripts/) ("attracticve" one)
* [Eta_Scripts](http://github.com/tjmustard/Eta_Scripts) (here)
* [PHYC Group](http://phyc.chem.oregonstate.edu/)

This repository started in ~2014, all the original code is in my [thesis] (https://ir.library.oregonstate.edu/concern/graduate_thesis_or_dissertations/xp68kn058) and some updates were made over the years. The repository was deleted/taken down and then reuploaded with the same name. Recently an update to python 3 was made.


## How is everything organized?
Folder Setup:
    
    Eta_Scripts/
     |
      ——> EtaLib/      #Contains the majority of the code functions
          G0X/         #Gaussian related scripts
          Jaguar/      #Jaguar related scripts
          Turbomole/   #Turbomole related scripts
          Psientific/  #Programs written by the PHYC group and myself (and others when they happen)
          “PROGRAM”
             |
              ———> /scripts   #Contains the all the user accessible scripts
                   /snippets  #Contains files that augment the scripts
                   /aux       #Scripts and files that are called directly from other scripts and are
                               not meant to be user accessible

## Snippets:

Snippets are text files holding configuration information. In the case of input file generation, snippets hold the method and basis information, solvation, optimization settings, geometric restraints, etc. For submission scripts the snippet file holds particular scripting steps necessary for proper submission.



## Installation:

1. Unpack the "master.zip" file.
 
2. Add the below profile settings to your .zshrc or zshrc-global file. (For other shell environments you will need to edit the below as needed.)

3. Edit the first line to denote the absolute path to the master Eta_Scripts folder.

4. Reload your profile.

## Requirements:
* Python 2.7
  * All Python code was written in Python 2.7 and as such all the scripts ask for python2.7.
* Numpy 1.6+
  * Not sure if 1.6 is the minimum but I know it works with 1.6+
* Quatfit
  * Quatfit is needed for all scripts that utilize alignment of molecules. It is free and be downloaded from [here](http://www.ccl.net/cca/software/SOURCES/C/quaternion-mol-fit/).


##TODO

I have yet to write an install script. But for those of you who are linux/unix inclined you can figure it out for now.
I hope to have it done soon.

Incorporate the Quatfit python library into nMAP

## ZSH profile setup:
Place this text into your ZSH profile to

    ################################################################################
    #Eta_Scripts profile settings:
    ################################################################################
    export ETADIR=DIRECTORY_TO_SCRIPTS/Eta_Scripts
    export PATH=$ETADIR:$PATH
    export PYTHONPATH=$PYTHONPATH:$ETADIR/EtaLib
    
    #Gaussian snippet settings:
    export PATH=$ETADIR/G0X/scripts:$PATH
    alias 'nGT'='nGT.py'
    alias 'ngt'='nGT.py'
    alias 'nGT-all'='nGT-all.sh'
    alias 'ngt-all'='nGT-all.sh'
    alias 'nGList'='nGList.py'
    alias 'nGArchive'='nGArchive.py'
    compctl -g '*.log' nGT.py
    compctl -g '*.log' nGT
    compctl -g '*.log' ngt
    compctl -g "$ETADIR/G0X/snippets/*.g0x" nTranslate-XYZ-G09Input.py
    compctl -g "$ETADIR/G0X/snippets/*.g0x" nTranslate-G0XOutput-G0XInput.py
    
    #Turbomole snippet settings:
    export PATH=$ETADIR/Turbomole/scripts:$PATH
    alias 'nTMT'='nTMT.py'
    alias 'ntmt'='nTMT.py'
    alias 'nTMT-all'='nTMT-all.sh'
    alias 'ntmt-all'='nTMT-all.sh'
    alias 'nTMList'='nTMList.py'
    alias 'nTMArchive'='nTMArchive.py'
    compctl -g '*.Turbomole' nTMT.py
    compctl -g '*.Turbomole' nTMT
    compctl -g '*.Turbomole' ntmt
    compctl -g "$ETADIR/Turbomole/snippets/*.tm" nTranslate-XYZ-TMInput.py
    compctl -g "$ETADIR/Turbomole/snippets/*.tm" nTranslate-TMOutput-TMInput.py
    
    #Jaguar snippet settings:
    export PATH=$ETADIR/Jaguar/scripts:$PATH
    alias 'nJT'='nJT.py'
    alias 'njt'='nJT.py'
    alias 'nJT-all'='nJT-all.sh'
    alias 'njt-all'='nJT-all.sh'
    alias 'nJList'='nJList.py'
    alias 'nJArchive'='nJArchive.ppy'
    compctl -g '*.out' nJT.py
    compctl -g '*.out' njt
    compctl -g '*.out' nJT
    compctl -g "$ETADIR/Jaguar/snippets/*.jag" nTranslate-XYZ-JaguarInput.py
    compctl -g "$ETADIR/Jaguar/snippets/*.jag" nTranslate-JaguarOutput-JaguarInput.py
    
    #Terachem snippet settings:
    export PATH=$ETADIR/Terachem/scripts:$PATH
    compctl -g '$ETADIR/Terachem/snippets/*.tc' nTranslate-XYZ-TCInput.py
    compctl -g '$ETADIR/Terachem/snippets/*.tc' nTranslate-TCOutput-TCInput.py
    compctl -g '*.tc' Submit.One.TC*
    
    #VASP snippet settings:
    export PATH=$ETADIR/VASP/scripts:$PATH
    compctl -g '$ETADIR/VASP/snippets/*.vasp' nTranslate-XYZ-VASPInput.py
    compctl -g '$ETADIR/VASP/snippets/*.vasp' nTranslate-VASPOutput-VASPInput.py
    compctl -g '*.???' Submit.One.VASP*
    
    #Psientific settings
    export PATH=$ETADIR/Psientific:$PATH
    compctl -g '*.map' nMAP.py
    compctl -g '*.dep' nDeprotonate.py
    compctl -g '*.di' nDist-Int-Splitter.py
    
    #Babel settings
    export PATH=$ETADIR/Babel/scripts:$PATH
    
    ################################################################################


##Scripts Explained:
###nTranslate-... Scripts:
In general these scripts "translate" one file format to another.

####nTranslate-XYZ-G09Input.py
This script parses a XYZ geometry file and a snippet file to create a Gaussian input file. Auto completion is most
useful for snippets. You can implement this in the terminal using the above ZSH profile section

    nTranslate-XYZ-G09Input.py <snippet> -c <charge> -m <multiplicity> --mod "constraint" --ionic=reg -d
    -----------------------------------------------------------------------------------------------------------------
        -c    (0)              | Charge of system.
                          |
        -m    (1)           | Multiplicity of system.
                          |    Usually the .log or .out file.
                          |
        --mod="B 1 2 F"     | Adds a modredundant restraint.
                          |
        --ionic=reg          | Use an internal list for charge.
                          |
        -d                  | Print debug info
                          |
        -h                  | Print this help screen
                          |
    -----------------------------------------------------------------------------------------------------------------
    Example: nTranslate-XYZ-G09Input.py $ETA-SCRIPTS/G0X/snippets/SP_SOLV.PCM.TOLUENE=B3LYP.6-31Gd.g0x -c -1 -m 1
             --mod="B 1 2 F" --ionic=reg -d

###nXT.py:
These scripts parse an output file and print out a formatted information page.

    nGT.py <output.log>
    nJT.py <output.log>
    nTMT.py <output.log>
####nGT.py:
Example of use with Gaussian output file.

    # nGT.py H2O.log
    Analyzing Gaussian Output File: H2O.log
    Using  Gaussian 09, Revision B.01
    [#Processors=1  Memory=4GB  CheckPoint=H2O.chk]
    ==============================================================================
    # b3lyp/6-31+G** gfprint gfinput scf=(direct,tight,maxcycle=300,xqc) o
    pt=(maxcycle=250) freq=noraman
    ------------------------------------------------------------------------------
    Pointgroup=CS    Stoichiometry=H2O    CS[SG(H2O)]
    Charge = 0    Multiplicity = 1    Basis = 29
    ------------------------------------------------------------------------------
    SCF Energy = -76.4340476899    Predicted Change = -1.033990D-09
    ------------------------------------------------------------------------------
    ------------------------------------------------------------------------------
    Statistical Thermodynamic Analysis for H2O.log
    Temperature   298.150 Kelvin.  Pressure   1.00000 Atm.
    ==============================================================================
    SCF Energy=                                             -76.4340476899
    Zero-point correction (ZPE)=            0.021287        -76.412761
    Internal Energy (U)=                    0.024122        -76.409925
    Enthalpy (H)=                           0.025067        -76.408981
    Gibbs Free Energy (G)=                  0.002984        -76.431064
    Entropy (S) [Cal/Mol-Kelvin]=           46.478
    Entropy (S): electronic=                0.000
    Entropy (S): translational=             34.608
    Entropy (S): rotational=                11.862
    Entropy (S): vibrational=               0.008
    Frequencies: 2306.68  5480.63  5656.25
    ==============================================================================
    ##################              JOB COMPLETED               ##################
    ==============================================================================
   
###nXList.py:
These scripts parse an output file and print out a list of formatted information.

    nGList.py
    nJList.py
    nTMList.py
####nGList.py:
Example of use with Gaussian output file.

    # nGList.py
    --------------------------------------------------------------------------------------------------------------------------
         SCF           ZPE        Internal    Enthalpy     Gibbs     Freq1    Freq2  NBasis  PntGrp  Stoichiometry   File Name
    --------------------------------------------------------------------------------------------------------------------------
    -76.4340476899  -76.412761  -76.409925  -76.408981  -76.431064  2306.68  5480.63  29  CS[SG(H2O)]     H2O        H2O.log

###nXArchive.py:
These scripts read output files and create an ARCHIVE folder for easy storage of jobs.

    nGArchive.py
    nJArchive.py
    nTMArchive.py

Running this script will create the ARCHIVE folder and several sub folders. These will hold:
* Formatted .pub file for easy implementation into Supporting Information
* A XYZ file of the final geometry.
* Two comma separated table files of energy data:
  * An Equationed table that uses the RC format for copy and pasting into spreadsheets.
  * A Totals table that has the precomputed ZPE, U, H, S, and G.
* SOON TO COME: log file storage, input file storage, joblog storage, etc.
   
####Snippets Style:
Snippets are used to augmented the Eta_Scripts

    #CONFIG 'Starts the config section'
    'These lines hold the first lines in a Gaussian input file'
    #END 'Ends the config section'
    #ECP 'Starts ECP configuration'
    'This line holds any ECP information'
    #END 'Ends ECP section'
    #CLOSE 'This section starts the closing section'
    'These lines are for after the geometry and ecp'
    'They can be used to start subsequent jobs ie "link1"'
    #END 'Ends the close section'
Only the Gaussian nTranslate scripts (i.e. nTranslate-XYZ-G0XInput.py) implement the ECP and CLOSE sections. As new programs are supported this section will expand as necessary.


#Psientific Code:
This code was written starting in 2012 during my tenure as a graduate student at Oregon State University and continuing to present.

###nMAP.py:
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

###nDist-Int-Splitter.py:
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

###nDeprotonate.py:
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

###nPersistentOTS.py:
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
---------|:------------------|:---------------------------------------------------------------------|
#VARIABLE=|bond 35 6        |Will scan along the bond between atoms 35 and 6|
         |angle 3 5 2      |Will scan along the angle between 3-5-2 with 5 being the middle atom.|
         |dihedral 6 7 3 4 |Will scan along the dihedral 6-7-3-4.|

***
**\#SCAN**

If you wish to do a SCAN 

Variable|Option|Definition
---------|:------------------|:---------------------------------------------------------------------|
#SCAN=|yes |SCAN function is turned on
 |no (or anything else) |SCAN function is tuned off.

***
**\#TS**

If you wish to implement the TS-SCAN logic 

Variable|Option|Definition
---------|:------------------|:---------------------------------------------------------------------|
#TS=|yes |TS-SCAN logic is turned on.
 |dyn |TS-SCAN logic will be turned on when the energy starts to go uphill. Useful for starting scan from far out.
 |no |TS-SCAN logic is off.

***
**\#SCANLENGTH**

How many maximum steps you wish to take. 

Variable|Option|Definition
---------|:------------------|:---------------------------------------------------------------------|
#SCANLENGTH=|30 |The scan will run for a maximum of 30 steps.

***
**\#SCANSTEPNOW**

What step numper you are currently on. You can edit this number to prematurely kill a SCAN.

Variable|Option|Definition
---------|:------------------|:---------------------------------------------------------------------|
#SCANSTEPNOW=|0 |This is the current step number you are on. This is editable.

***
**\#STEPSIZE**

How much the variable will change between SCAN steps. 

Variable|Option|Definition
---------|:------------------|:---------------------------------------------------------------------|
#STEPSIZE=|-0.1 |This number will be in Angstrom for bonds and Degrees for angle and dihedral.

***
**\#WAIT**

If you want the TS-SCAN logic to wait a number of steps before initializing. 

Variable|Option|Definition
---------|:------------------|:---------------------------------------------------------------------|
#WAIT=|0 |Wait X steps before turning on the TS-SCAN logic.

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
---------|:------------------|:---------------------------------------------------------------------|
#MINTSSCANSIZE=|none | x = 50 (C-C bond: 0.021Å)
 |loose | x = 100 (C-C bond: 0.042Å)
 |tight | x = 25 (C-C bond: 0.010Å)
 |verytight | x = 5 (C-C bond: 0.002Å)

***
**\#MINORTS**

If we should explore minor TSs. 

Variable|Option|Definition
---------|:------------------|:---------------------------------------------------------------------|
#MINORTS|=one| Allow for one minor TS to be skipped. Useful if a geometric change is necessary for the TS to occur.
 |=yes| Allow for an infinite number of minor TSs to be skipped.
 |=no| Do not allow a single minor TS to be skipped. Once the energy goes down nPersistenceOTS will go back to the last geometry and take a smaller step size.

###nVisualize-Pymol-XYZ.py:
This will create a Pymol script for all the XYZ or PDB (nVisualize-Pymol-PDB.py) files in the current folder. With this script Pymol is launched using this script and assists in giving visually apealing pictures.

Several types of measurements can be augmented by renaming them with a specific suffix:
* Hydrogen bonds (*_h)
* Interactions (*_int)
* Dihedrals (*_dih)
* Angles (*_angle)
* Bond forming and braking TSs (*_ts)

Once renamed these settings can be activated by going to the 'File' menu clicking on 'Run' and selecting the PymolV1.3-Visualize.pml file. This reruns a sebset of the original script to augment the measurements.


I hope you like the scripts!

TJ Mustard

