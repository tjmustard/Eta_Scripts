#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright (c) 2014, Thomas J. L. Mustard, O. Maduka Ogba, Paul Ha-Yeon Cheong
#
# PHYC Group
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

# Eta_Scripts default imports

import getopt
import math
import os
import shutil
import sys
import time
from decimal import *
from itertools import *
from sys import *

from EtaLib import etaatom, etanumpy

# Script specific imports


### --- Arguments --- ###

program = "nDist-Int-Splitter.py"
ifile = ""
ofile = ""
outputFolder = ""
parentFile = ""
verbose = 0
debug = 0
template = 0
start_time = time.time()
allInDir = 0

# Read command line args

try:
    (myopts, args) = getopt.getopt(sys.argv[1:], "i:p:hvtdD")
except getopt.GetoptError:
    print(program + " -i <inputfile.map> -d")
    sys.exit(2)

###############################
# o == option
# a == argument passed to the o
###############################

for (o, a) in myopts:
    if o == "-i":
        ifile = a
    elif o == "-p":
        deprotonations = int(a)
    elif o == "-h":
        print(program + "-i <inputfile.map>")
        sys.exit(0)
    elif o == "-v":
        verbose = 1
    elif o == "-d":
        allInDir = 1
    elif o == "-D":
        debug += 1
    elif o == "-t":
        template = 1
    else:
        print("Usage: %s -i inputfile.map" % sys.argv[0])
        sys.exit(0)

if template == 1:
    print("#Input file name")
    print("input = .xyz")
    print("#Output folder and name")
    print("output = ")
    print("#Molecule sections")
    print("$SECTIONS       <-- Each line is it's own section")
    print("1 2 3 4 5 6     <-- Atoms 4 and 6 are one section")
    print("7-10            <-- Atoms 7 through 10 are a different section")
    print("$ENDSECTIONS")
    sys.exit(0)

### --- Open input file and parse into variable --- ###

f = open(ifile)
inputList = f.readlines()
f.close()

sectionSwitch = "off"
sections = []
sectionsTemp = []
totalSections = 0
lineTemp = []
sectionLine = []

for var in inputList:
    line = var.split()
    if line[0] == "$SECTIONS":
        sectionSwitch = "on"
    if line[0] == "$ENDSECTIONS":
        sectionSwitch = "off"
    elif line[0] == "input":
        parentFile = line[2]
    elif line[0] == "output":
        outputFolder = line[2]
    if sectionSwitch == "on":
        if line[0][0] != "$" and line[0][0] != "#":
            for i in line:
                lineTemp.append(i)

            # print i

            totalSections += 1
            sectionsTemp.append(list(lineTemp))
            lineTemp[:] = []

### --- Expand upon the parentDeleteTemp list if
# there are certain characters to denote multiple atoms.
###

for i in range(len(sectionsTemp)):
    for j in range(len(sectionsTemp[i])):
        if "-" in str(
            sectionsTemp[i][j]
        ):  # Use the '-' symbol to denote all atoms between the former and latter
            line = str(sectionsTemp[i][j]).split("-")
            start = int(line[0])
            for k in range(int(line[0]), int(line[1]) + 1):

                # print k

                sectionLine.append(k)
        if "-" not in str(
            sectionsTemp[i][j]
        ):  # If not one of the special symbols just add it to the final list
            sectionLine.append(int(sectionsTemp[i][j]))
    sections.append(list(sectionLine))
    sectionLine[:] = []

# print sections

#### --- Script specific functions --- ###
###################################################################################

### --- Print out some general data --- ###

print(
    "Breaking up "
    + parentFile
    + " into "
    + str(len(sections))
    + " sections with atoms of:"
)
for i in range(len(sections)):
    line = ""
    for j in sections[i]:
        line += str(j) + " "
    print(line)
print("The output folder is: ./" + outputFolder)
print(
    "###################################################################################"
)

#### --- Make the output folder --- ###

etaatom.make_dir(outputFolder)

#### --- Open parent file --- ###

ifilelol = etaatom.xyz_lol(parentFile)

### --- Algorythm to find all the possible section combinations --- ###
### --- Generate all combinations WITHOUT replacements!
### --- For instance if there are 3 or 4 sections to your molecule you
### --- need to grab all the single, doubles and possibly triples.

permutationsList = []
for i in range(len(sections) - 1):
    permutationsListTemp = list(combinations(list(range(len(sections))), i + 1))
    for j in permutationsListTemp:
        permutationsList.append(j)

# print "permutationsList list:"
# print permutationsList

### --- Output the permutations in xyz files --- ###

print("Building and writing file:")
for i in range(len(permutationsList)):
    ofile = ""
    for j in range(len(permutationsList[i])):

        # print permutationsList[i][j]

        ofile += "S" + str(permutationsList[i][j] + 1) + "-"
    ofile += parentFile
    print(ofile)
    atomList = []
    for j in range(len(permutationsList[i])):

        # print sections[permutationsList[i][j]]

        for k in sections[permutationsList[i][j]]:
            atomList.append(k)

    # print atomList

    f = open(outputFolder + "/" + ofile, "w")
    f.write(str(len(atomList)) + "\n")
    f.write(ofile + "\n")
    for j in atomList:
        f.write(
            ifilelol[j + 1].e
            + "  "
            + str("{:.6f}".format(ifilelol[j + 1].x))
            + "  "
            + str("{:.6f}".format(ifilelol[j + 1].y))
            + "  "
            + str("{:.6f}".format(ifilelol[j + 1].z))
            + "\n"
        )
    f.close()

######################################################################
### END OF SCRIPT
######################################################################
