#!/usr/bin/env python2.7
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

import os
import math
from sys import *
import sys
import getopt
import etaatom
import etanumpy
from decimal import *
from numpy import *

### --- Arguments --- ###

program = 'nTranslate-XYZ-Renumber.py'
mainfile = ''
debug = 0
autofind = 0
atomNum = -1

# Define what type of job we will be running

if len(argv) == 1:
    argv.append('-a')
    autofind = 1
elif argv[1] == '-h':
    print '======================================================================================================'
    print program \
        + '''<main file>.xyz -a <auto select> -h
======================================================================================================

This script will renumber all the xyz files within the current directory to a specified key structure.

If any of the structures within the directory are larger than the specified key file,
any unmatched atoms will be appended to the end of the file.

All renumbered files are saved in the folder "RENUMBERED" including the key structure.

If you do not specify a structure the smallest one will be chosen for you.'''
    sys.exit(0)
else:
    mainfile = argv[1]

### --- Make RENUMBERED folder --- ###

if not os.path.exists('RENUMBERED'):
    os.makedirs('RENUMBERED')

### --- If asked to do so iterate through the folder for other .xyz files and find the smallest one --- ###

if autofind == 1:
    for i in os.listdir(os.getcwd()):
        if i.endswith('.xyz'):
            f = open(i, 'r')
            for (j, line) in enumerate(f):
                if j >= 1:
                    break
                elif atomNum == -1:
                    atomNum = int(line)
                    mainfile = i
                elif int(line.strip()) < atomNum:
                    atomNum = int(line)
                    mainfile = i
            f.close()

### --- Open parent file --- ###

mainfileAtom = etaatom.xyz_lol(mainfile)

print 'Renumbering all structures to ' + mainfile + '.'

# Output the main file to the renumbered folder

etaatom.output_xyz('RENUMBERED/' + mainfile, mainfileAtom)

# Iterate through the folder for other .xyz files

for i in os.listdir(os.getcwd()):
    if i.endswith('.xyz') and i != mainfile:
        ifile = i

    # Print out the current job for the user

        print ifile

    # Get the XYZ lol/Atom data from the child file (ifile)

        childfileAtom = etaatom.xyz_lol(ifile)

    # Build a holder for the Atom data and copy the first two lines

        if len(childfileAtom) > len(mainfileAtom):
            ofileAtom = [0] * len(childfileAtom)
        else:
            ofileAtom = [0] * len(mainfileAtom)
        ofileAtom[0] = childfileAtom[0]
        ofileAtom[1] = childfileAtom[1]

    # Make a holder for matched and unmatched atom numbers to use later

        matchedatoms = []
        matchedatomsdist = []
        unmatchedatoms = []

    # Build a holder for the distance matrix

        distanceMatrix = [[0 for x in xrange(len(childfileAtom) - 2)]
                          for x in xrange(len(mainfileAtom) - 2)]

    # Build the distance matrix

        for j in range(2, len(mainfileAtom)):
            for k in range(2, len(childfileAtom)):
                distanceMatrix[j - 2][k - 2] = \
                    etanumpy.two_structure_get_distance(j, k,
                        mainfileAtom, childfileAtom)

    # Find the best nearest neighbor and if not already used

        for j in range(len(mainfileAtom) - 2):
            nearestNeighbor = 0
            nearestNeighborDist = distanceMatrix[j][0]
            for k in range(1, len(childfileAtom) - 2):
                if distanceMatrix[j][k] < nearestNeighborDist:
                    nearestNeighborDist = distanceMatrix[j][k]
                    nearestNeighbor = k

            if nearestNeighbor + 2 not in matchedatoms:

        # ofileAtom[j+2] = childfileAtom[nearestNeighbor + 2]

                matchedatoms.append(nearestNeighbor + 2)
                matchedatomsdist.append(nearestNeighborDist)
            elif nearestNeighbor + 2 in matchedatoms:

                for l in range(len(matchedatoms)):

          # print str(l) + "   " + str(matchedatoms[l]) + "  " + str(matchedatomsdist[l])

                    if nearestNeighbor + 2 == matchedatoms[l] \
                        and nearestNeighborDist < matchedatomsdist[l]:

            # ofileAtom[j+2] = childfileAtom[nearestNeighbor + 2]

                        matchedatoms[l] = nearestNeighbor + 2
                        matchedatomsdist[l] = nearestNeighborDist

    # Iterate throught the matched atom list and add them to the end on the ofile atom/lol

        for k in range(len(matchedatoms)):

      # print unmatchedatoms[k]

            ofileAtom[2 + k] = childfileAtom[matchedatoms[k]]

    # Iterate through the child file and if an atom is not matched add it to the unmatched atom list

        for k in range(2, len(childfileAtom)):
            if k not in matchedatoms:

        # print k

                unmatchedatoms.append(k)

    # Iterate throught the unmatched atom list and add them to the end on the ofile atom/lol

        for k in range(len(unmatchedatoms)):

      # print unmatchedatoms[k]

            ofileAtom[len(matchedatoms) + 2 + k] = \
                childfileAtom[unmatchedatoms[k]]

    # Output the newly renumbered file to the renumbered folder

        etaatom.output_xyz('RENUMBERED/' + ifile, ofileAtom)

######################################################################
### END OF SCRIPT
######################################################################
