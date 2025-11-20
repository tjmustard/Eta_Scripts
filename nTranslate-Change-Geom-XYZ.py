#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright (c) 2014, Thomas J. L. Mustard, O. Maduka Ogba, Paul Ha-Yeon Cheong
#
#   PHYC Group
#   Oregon State University
#   College of Science
#   Department of Chemistry
#   153 Gilbert Hall
#   Corvallis OR, 97331
#   E-mail:  mustardt@onid.orst.edu
#   Ph.  (541)-737-2081
#   http://phyc.chem.oregonstate.edu/
#
#   All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are met:
#
#   * Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#
#   * Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#
#   * Neither the name of the {organization} nor the names of its
#     contributors may be used to endorse or promote products derived from
#     this software without specific prior written permission.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import getopt
import math
import os
import sys
from decimal import *
from sys import *

from numpy import *

from EtaLib import etaatom, etanumpy

### --- Arguments --- ###

program = "nTranslate-Change-Geom-XYZ.py"
debug = 0
ifile = ""
ofile = ""
printfile = 0
bond_angle_dihedral = []
newvalue = ""
allInDir = 0

### Read command line args

try:
    (myopts, args) = getopt.getopt(sys.argv[1:], "i:o:b:v:phd", ["bad=", "var="])
except getopt.GetoptError:
    print(
        program
        + " -i <inputfile.xyz> -o <outputfile.zmat> -b or --bad=<bond/angle/dihedral atom numbers>   -v or --newvalue=<new geometric value> "
    )
    sys.exit(2)

###############################
# o == option
# a == argument passed to the o
###############################

for (o, a) in myopts:
    if o == "-i":
        ifile = a
    elif o == "-o":
        ofile = a
    elif o == "--bad" or o == "-b":
        for item in a.split():
            bond_angle_dihedral.append(int(item))
    elif o == "--newvalue" or o == "-v":
        newvalue = a
    elif o == "-p":
        printfile = 1
    elif o == "-d":
        allInDir = 1
    elif o == "-D":
        debug += 1
    elif o == "-h":
        print(
            program
            + " -i <inputfile.xyz> -o <outputfile.zmat> --bad=<bond/angle/dihedral atom numbers>  --newvalue=<new geometric value> "
        )
        sys.exit(0)
    else:
        print(
            "Usage: %s  -i <inputfile.xyz> -o <outputfile.zmat> --bad=<bond/angle/dihedral atom numbers>  --newvalue=<new geometric value> "
            % sys.argv[0]
        )
        sys.exit(0)


def runChange(ifile, ofile, bond_angle_dihedral):

    # ## --- Open parent file --- ###

    ifilelol = etaatom.xyz_lol(ifile)

    # ## --- Renumber atoms to put the atoms of interest in the

    ifilelol_renumbered = etaatom.renumber_basic(ifilelol, bond_angle_dihedral)

    # #Change the bond/angle/dihedral of interest
    # if len(bond_angle_dihedral) == 2:
    #  list1 = etaatom.buildbondlist(bond_angle_dihedral[0], bond_angle_dihedral[1], ifilelol)
    #  list2 = etaatom.buildbondlist(bond_angle_dihedral[1], bond_angle_dihedral[0], ifilelol)
    # elif len(bond_angle_dihedral) == 3:
    #  list1 = etaatom.buildbondlist(bond_angle_dihedral[1], bond_angle_dihedral[2], ifilelol)
    #  list2 = etaatom.buildbondlist(bond_angle_dihedral[2], bond_angle_dihedral[1], ifilelol)
    # elif len(bond_angle_dihedral) == 4:
    #  list1 = etaatom.buildbondlist(bond_angle_dihedral[1], bond_angle_dihedral[2], ifilelol)
    #  list2 = etaatom.buildbondlist(bond_angle_dihedral[2], bond_angle_dihedral[1], ifilelol)

    # ## --- Generate a list of lists for the zmatrix of the structure --- ###

    zmatlol = []
    for i in range(0, len(ifilelol_renumbered)):
        linetemp = etanumpy.get_zmat(i, i - 1, i - 2, i - 3, ifilelol_renumbered)
        zmatlol.append(linetemp)

    # Change the bond/angle/dihedral of interest

    if len(bond_angle_dihedral) == 2:
        zmatlol[3][3] = newvalue
    elif len(bond_angle_dihedral) == 3:
        zmatlol[4][5] = newvalue
    elif len(bond_angle_dihedral) == 4:
        zmatlol[5][7] = newvalue

    # Grab the XYZ from the now changed zmat

    ifilelol_renumbered_fixed = etanumpy.get_xyz_from_zmat(zmatlol)

    # ## --- RE-Renumber atoms to put the atoms of interest in the
    # Iterate through the bond_angle_dihedral list and put them back to their locations

    for (i, atom) in enumerate(bond_angle_dihedral):
        ifilelol[atom + 1] = ifilelol_renumbered_fixed[i + 2]

    # Iterate throught the ifilelol and add the rest of the file to the renumbered

    usedBAD = 0
    for i in range(2, len(ifilelol)):
        if i - 1 not in bond_angle_dihedral:
            ifilelol[i] = ifilelol_renumbered_fixed[
                i + len(bond_angle_dihedral) - usedBAD
            ]
        else:
            usedBAD += 1

    # ## --- Output  format --- ###

    etaatom.output_xyz(ofile, ifilelol)

    # ## --- Print out the xyz from the zmatrix --- ###

    if printfile == 1:
        for i in range(len(ifilelol)):
            if i <= 1:
                print(ifilelol[i])
            else:
                print(
                    ifilelol[i].e
                    + "  "
                    + str("{:.6f}".format(ifilelol[i].x))
                    + "  "
                    + str("{:.6f}".format(ifilelol[i].y))
                    + "  "
                    + str("{:.6f}".format(ifilelol[i].z))
                )
    return


try:
    if allInDir == 1:
        for i in os.listdir(os.getcwd()):
            if i.endswith(".xyz"):
                ifile = i

                # ofile = "OUT-" + i

                ofile = etaatom.basename(ifile, ".xyz") + "-OUT.xyz"
                print(ifile + " ----> " + ofile)
                runChange(ifile, ofile, bond_angle_dihedral)
    else:
        print(ifile + " ----> " + ofile)
        runChange(ifile, ofile, bond_angle_dihedral)
except:
    print("Something failed. Make sure the/a xyz file(s) exists.")

    # #####################################################################
    # ## END OF SCRIPT
    # #####################################################################
