Eta_Scripts
======

## Why Eta_Scripts?

These scripts were originally written to assist in creating and reading computational job files. Most, if not all, of the scripts are written with Python. Many of the scripts have been removed due to lack of compatibility.

Credits:

These scripts were written mostly by myself, but I had some help and they deserve recognition:

* Prof. Paul Ha-Yeon Cheong was the inspiration for this project. While in grad school Paul wrote several scripts to assist his research. Most of this project is inspired by those scripts.
* O. Maduka Ogba
* Joshua J. Kincaid
* The entire PHYC group for their inspiration and support. (Ryne C. Johston, O. Maduka Ogba, Wojtek Rajski, Ommi Pattawong, Lindsay Wills, Daniel Walden, Kevin Snyder, Jacob Buchanan, Ben Hanken, I-Ya Chang)

You may acknowledge the use of the scripts/programs in published material as:
Thomas J. L. Mustard, O. Maduka Ogba, Paul Ha-Yeon Cheong. unpublished results.

* [Eta_Scripts](http://tjmustard.github.io/Eta_Scripts/) ("attractive" one)
* [Eta_Scripts](http://github.com/tjmustard/Eta_Scripts) (here)
* [PHYC Group](http://phyc.chem.oregonstate.edu/)


## How is everything organized?
Folder Setup:
    
    Eta_Scripts/
     |
      ——> EtaLib/      #Contains the majority of the code functions
          Pymol/       #Pymol related scripts
          
          “PROGRAM”
             |
              ———> /scripts   #Contains the all the user accessible scripts
                   /snippets  #Contains files that augment the scripts
                   /aux       #Scripts and files that are called directly from other scripts and are
                               not meant to be user accessible

## Snippets:

Snippets are text files holding configuration information. In the case of input file generation, snippets hold the method and basis information, solvation, optimization settings, geometric restraints, etc. For submission scripts the snippet file holds particular scripting steps necessary for proper submission.



## Installation:

1. Unpack the "master.zip" file.
 
2. Add the below profile settings to your .bashrc, .zshrc, or similar file. (For other shell environments you will need to edit the below as needed.)

3. Edit the first line to denote the absolute path to the master Eta_Scripts folder.

4. Reload your profile.

## Requirements:
* Python 3
  * All Python code was written in Python 3 and as such all the scripts ask for python3.

## ZSH profile setup:
Place this text into your SH profile to

    ################################################################################
    #Eta_Scripts profile settings:
    ################################################################################
    export ETADIR=DIRECTORY_TO_SCRIPTS/Eta_Scripts
    export PATH=$ETADIR:$PATH
    export PYTHONPATH=$PYTHONPATH:$ETADIR/EtaLib
    
    #Pymol settings
    export PATH=$ETADIR/Pymol/scripts:$PATH
    
    ################################################################################


##Scripts Explained:
###nVisualize-... Scripts:
In general these scripts copy the XYZ or PDB files into a subdirectory named "Pymol-Picture" as well as some pymol scripts. You can then open the Visual-Script.pml from that directory or from within Pymol to load the structures and styles.


I hope you like the scripts!

TJ Mustard
