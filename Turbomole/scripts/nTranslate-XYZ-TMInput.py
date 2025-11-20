#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#   Copyright (c) 2014, Thomas J. L. Mustard, O. Maduka Ogba, Paul Ha-Yeon Cheong
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
from sys import *

from EtaLib import etaatom

### --- Arguments --- ###

program = "nTranslate-XYZ-TMInput.py"

# Grab the first argument from the command and use that as the snippet

try:
    snippet = sys.argv[1]
except IndexError:
    snippet = "-h"

# If help is wanted allow the skipping of a snippet

if snippet == "-h":
    argv.append("-h")
    argv.append("-h")

# Unless the charge is stated via the -s argument the charge will be assumed to be 0

charge = 0

# Unless the multiplicity is stated via the -m argument the multiplicity will be assumed to be 1

multi = 1
modred = []
writemod = 0
ionicCM = ""
ionicswitch = 0
debug = 0

# Constants

auToAngstrom = 1.889725988579

### Read command line args

try:
    (myopts, args) = getopt.getopt(sys.argv[2:], "c:m:dh", ["mod=", "ionic="])
except getopt.GetoptError:
    print(
        program
        + """ <snippet> -c <charge> -m <multiplicity> --mod \"constraint\" --ionic=reg/mno -d"""
    )
    sys.exit(2)

###############################
# o == option
# a == argument passed to the o
###############################

for (o, a) in myopts:
    if o == "-c":
        charge = a
    elif o == "-m":
        multi = a
    elif o == "-d":
        debug += 1
    elif o == "--mod":
        modred.append(int(a))
        writemod = 1
    elif o == "--ionic":
        ionicCM = a
        ionicswitch = 1
    elif o == "-h":
        print(
            program
            + """ <snippet> -c <charge> -m <multiplicity> --mod \"constraint\" --ionic=reg/mno -d"""
        )
        sys.exit(0)
    else:
        print(
            'Usage: %s  <snippet> -c <charge> -m <multiplicity> --mod "constraint" --ionic=reg/mno -d'
            % sys.argv[0]
        )
        sys.exit(0)

if debug >= 1:
    print("Charge: " + charge)
    print("Multiplicity: " + multi)
    print("Snippet: " + snippet)
    print("Modredundant Variable: ")
    print(modred)

# Grab the snippet information

(config, ecp, closelines, writeecp, writeclose) = etaatom.parse_snippet(snippet)

if debug >= 1:
    print("Configuration:")
    print(config)
    print("ECP Input: " + ecp)
    print("Additional Lines: ")
    print(closelines)

### --- Make NEWJOBS folder --- ###

if not os.path.exists("NEWJOBS.TM"):
    os.makedirs("NEWJOBS.TM")

### --- Open and parse the xyz files in the folder --- ###

print("Currently Translating:")
for i in os.listdir(os.getcwd()):
    if i.endswith(".xyz"):
        ifile = i
        print(i)

        # Parse XYZ file into a list of lists

        ifilelol = etaatom.xyz_lol(ifile)

        # ## --- Calculate charge from file --- ###

        if ionicswitch == 1:
            charge = etaatom.ionic(ionicswitch, ionicCM, ifilelol)

        # ## --- Out put the G0X input files --- ###

        etaatom.output_tm(
            ifile,
            config,
            charge,
            multi,
            writemod,
            modred,
            writeecp,
            [],
            writeclose,
            closelines,
            ifilelol,
        )

    # #####################################################################
    # ## END OF SCRIPT
    # #####################################################################
