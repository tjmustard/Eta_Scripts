#!/usr/bin/env python2.7
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

import os
import math
from sys import *
import sys
import getopt
from decimal import *
import shutil
import etaatom
import etanumpy

# Script specific imports

from itertools import *
import time

### --- Arguments --- ###

program = 'nDeprotonate.py'
ifile = ''
ofile = ''
verbose = 0
debug = 0
template = 0
start_time = time.time()

# Read command line args

try:
    (myopts, args) = getopt.getopt(sys.argv[1:], 'i:p:hvtdD')
except getopt.GetoptError:
    print program + ' -i <inputfile.map>'
    sys.exit(2)

###############################
# o == option
# a == argument passed to the o
###############################

for (o, a) in myopts:
    if o == '-i':
        ifile = a
    elif o == '-p':
        deprotonations = int(a)
    elif o == '-h':
        print program + ' -i <inputfile.map>'
        sys.exit(0)
    elif o == '-v':
        verbose = 1
    elif o == '-d':
        debug = 1
    elif o == '-D':
        debug = 2
    elif o == '-t':
        template = 1
    else:
        print 'Usage: %s -i inputfile.map' % sys.argv[0]
        sys.exit(0)

if template == 1:
    print '#Input file name'
    print 'input = .xyz'
    print '#Output folder and name'
    print 'output = '
    print '#Number of deprotonations'
    print 'deprotonations = 1'
    print '#Protons to delete'
    print '$PROTONS'
    print '4 6     <-- Hydrogens 4 and 6 are equivalent'
    print '8 9 10  <-- Hydrogens 8 9 & 10 are equivalent'
    print '$ENDPROTONS'
    sys.exit(0)

### --- Open input file and parse into variable --- ###

f = open(ifile)
inputList = f.readlines()
f.close()

protonSwitch = 'off'
protons = []
protonsTemp = []
symDifProtons = 0
totalProtons = 0
protonsList = []

for var in inputList:
    line = var.split()
    if line[0] == '$PROTONS':
        protonSwitch = 'on'
    if line[0] == '$ENDPROTONS':
        protonSwitch = 'off'
    elif line[0] == 'input':
        parentFile = line[2]
    elif line[0] == 'output':
        outputFolder = line[2]
    elif line[0] == 'deprotonations':
        deprotonations = int(line[2])
    if protonSwitch == 'on':
        if line[0][0] != '$' and line[0][0] != '#':
            for i in range(len(line)):
                protonsList.append(int(line[i]))
                protonsTemp.append(int(line[i]))
                totalProtons += 1
            symDifProtons += 1
            protons.append(list(protonsTemp))
            protonsTemp[:] = []


### --- Script specific functions --- ###
##################################################################################

### --- get the error of two different numbers --- ###

def getErrorDist(a, b):
    error = (a + b) / 2 / deprotonations * 0.10

  # print "DIST"
  # print str((a+b)/2) + " " + str(abs(a-b)) + " " + str(error)

    return error


def getErrorDih(a, b):
    error = (a + b) / 2 / deprotonations * 0.10

  # print "DIH"
  # print str((a+b)/2) + " " + str(abs(a-b)) + " " + str(error)

    return error


### --- if all the items in a list are the same return true --- ###

def all_same(items):
    return all(x == items[0] for x in items)


### --- "find" a target item in a list and return the index --- ###

def find(target, listoi):
    for (i, lst) in enumerate(listoi):
        for (j, item) in enumerate(lst):
            if item == target:
                return (i, j)
    return (None, None)


##################################################################################

### --- Make the output folder --- ###

etaatom.make_dir(outputFolder)

### --- Open parent file --- ###

ifilelol = etaatom.xyz_lol(parentFile)

### --- Get bonds --- ###

covBonds = []
covHBonds = []
covTMBonds = []
nearestNeighbor = []
neighborStart = [0, 1000000]

### --- Generate bond lists --- ###

for i in range(0, int(ifilelol[0])):
    nearestNeighbor.append(list(neighborStart))
    for j in range(0, int(ifilelol[0])):
        if i != j:
            distij = etanumpy.get_distance(i + 2, j + 2, ifilelol)
            if j > i:
                if distij <= 2.25 and ifilelol[i + 2].e != 1 \
                    and ifilelol[j + 2].e != 1 and (ifilelol[i
                        + 2].e.lower() not in etaatom.elementLarge
                        or ifilelol[j + 2].e.lower()
                        not in etaatom.elementLarge):
                    ifilelol[i + 2].neighbors.append(j + 1)
                    ifilelol[j + 2].neighbors.append(i + 1)
                    ifilelol[i + 2].neighborsdist.append(distij)
                    ifilelol[j + 2].neighborsdist.append(distij)
                elif distij <= 3 and ifilelol[i + 2].e.lower() \
                    in etaatom.elementLarge and ifilelol[j
                        + 2].e.lower() in etaatom.elementLarge:

                    ifilelol[i + 2].neighbors.append(j + 1)
                    ifilelol[j + 2].neighbors.append(i + 1)
                    ifilelol[i + 2].neighborsdist.append(distij)
                    ifilelol[j + 2].neighborsdist.append(distij)
                elif distij <= 1.3 and (ifilelol[i + 2].e == 1
                        or ifilelol[j + 2].e == 1):

                    ifilelol[i + 2].neighbors.append(j + 1)
                    ifilelol[j + 2].neighbors.append(i + 1)
                    ifilelol[i + 2].neighborsdist.append(distij)
                    ifilelol[j + 2].neighborsdist.append(distij)
    ifilelol[i + 2].nearest = ifilelol[i + 2].neighbors[0]
    for k in range(len(ifilelol[i + 2].neighbors) - 1):
        if ifilelol[i + 2].neighborsdist[k] > ifilelol[i
                + 2].neighborsdist[k + 1]:
            ifilelol[i + 2].nearest = ifilelol[i + 2].neighbors[k + 1]

### --- Algorythm to find all the possible and chemically relavent combinations --- ###
### --- Generate all combinations WITHOUT replacements --- ###

permutationsList = list(combinations(protonsList, deprotonations))

# print "permutationsList list:"
# print permutationsList

# Permutations time elapsed

elapsed_time = time.time() - start_time
print 'Permutations time elapsed: ' + str(elapsed_time)

### --- Setup the lists for holding all the data --- ###

permutationsListofIndex = []
permutationsListofIndexAndLast = []
permutationsDuplicates = []
permutationsDeduplicated = []

### --- Iterate through the list of permutations and build an index identifier list and and index plus last atom number list --- ###

a_list = []
for i in range(0, len(permutationsList)):
    if debug == 1 and float(i) * 100 % len(permutationsList) == 0:
        print 'Building permutations list ' + str(float(i)
                / len(permutationsList) * 100) + ' %'
    for j in range(0, len(permutationsList[i])):

    # Bring back the list index for the list of list number What index is 5 in?

        a_list.append(find(permutationsList[i][j], protons)[0])

  # Append the new list to the permutationsListofIndex list

    permutationsListofIndex.append(list(a_list))

  # Append the last atom number to the list

    a_list.append(permutationsList[i][len(permutationsList[i]) - 1])
    permutationsListofIndexAndLast.append(list(a_list))
    a_list[:] = []
    loiduplicate = 0

  # ## --- Iterate through the identifiers and remove chemically duplicate deprotonations --- ###

    for k in range(i):

    # ........### --- Check if the identifier lists are identical and remove one and set the "list of index" (loi) duplicate switch --- ###
    # ........if permutationsListofIndex[k] == permutationsListofIndex[i] and permutationsListofIndex[k][-1] == permutationsListofIndex[i][0] and loiduplicate == 0:
    # ............print str(permutationsList[k]) + " " + str(permutationsList[i]) + " " + str(permutationsListofIndex[k]) + " " + str(permutationsListofIndex[i])
    # ............loiduplicate = 1
    # ............#print "yes"
    # ............permutationsList[i] = [-1]
    # ............break
    # ........### --- Check if the identifier plus atom number lists are identical and if not reset the "list of index" (loi) duplicate switch --- ###
    # ........elif permutationsListofIndex[k] != permutationsListofIndex[i] :
    # ............loiduplicate = 0
    # ............#print "no"
    # ## --- Check if the identifier plus atom number lists are identical  and remove one--- ###

        if permutationsListofIndexAndLast[k] \
            == permutationsListofIndexAndLast[i]:
            permutationsList[i] = [-1]
            break

#### --- Iterate through the identifiers and remove chemically duplicate deprotonations --- ###
# for i in range(len(permutationsListofIndexAndLast)):
# ....if debug == 1 and float(i) * 100 % len(permutationsList) == 0:
# ........print "Parsing permutations list " + str(float(i)/len(permutationsListofIndexAndLast)*100) + " %"
# ....for j in range(i+1, len(permutationsListofIndexAndLast)):
# ........### --- Check if the identifier lists are identical and remove one and set the "list of index" (loi) duplicate switch --- ###
# ........if permutationsListofIndex[i] == permutationsListofIndex[j] and permutationsListofIndex[i][-1] == permutationsListofIndex[i][0] and loiduplicate == 0:
# ............print str(permutationsList[i]) + " " + str(permutationsList[j]) + " " + str(permutationsListofIndex[i]) + " " + str(permutationsListofIndex[j])
# ............loiduplicate = 1
# ............#print "yes"
# ............permutationsList[j] = [-1]
# ........### --- Check if the identifier plus atom number lists are identical  and remove one--- ###
# ........if permutationsListofIndexAndLast[i] == permutationsListofIndexAndLast[j]:
# ............permutationsList[j] = [-1]
# ............#print str(permutationsListofIndexAndLast[i]) + "\t" + str(permutationsListofIndexAndLast[j])
# ........### --- Check if the identifier plus atom number lists are identical and if not reset the "list of index" (loi) duplicate switch --- ###
# ........elif permutationsListofIndex[i] != permutationsListofIndex[j] :
# ............loiduplicate = 0
# ............#print "no"
# print permutationsListofIndex
# print permutationsListofIndexAndLast
# print permutationsList
# print permutationsList

# Math dedup time elapsed

elapsed_time = time.time() - start_time
print 'Math dedup time elapsed: ' + str(elapsed_time)

### --- make a new deduplicated list without the removed list entries (those which had [-1]) --- ###

permutationsDeduplicatedNoGeom = [x for x in permutationsList if x
                                  != [-1]]

#### --- Remove degenerate geometries by compairing distance and dihedrals for each --- ###

if deprotonations > 1:

  # ## --- Generate a matrix of distances and dihedrals --- ###

    protonList = []
    for i in protons:
        for j in i:
            protonList.append(int(j))

  # print protonList

    geomdistmat = [[0 for i in range(ifilelol[0])] for j in
                   range(ifilelol[0])]
    geomdihmat = [[0 for i in range(ifilelol[0])] for j in
                  range(ifilelol[0])]
    for i in range(len(protonList)):
        for j in range(i, len(protonList)):

      # print str(protonList[i]) + " " + str(protonList[j])
      # print etanumpy.get_distance(protonList[i], protonList[j], ifilelol)

            geomdistmat[protonList[i] - 1][protonList[j] - 1] = \
                etanumpy.get_distance(protonList[i], protonList[j],
                    ifilelol)
            geomdistmat[protonList[j] - 1][protonList[i] - 1] = \
                geomdistmat[i][j]
            if ifilelol[protonList[i] + 1].nearest \
                != ifilelol[protonList[j] + 1].nearest:

        # print etanumpy.get_dihedral(protonList[i], ifilelol[protonList[i]+1].nearest, ifilelol[protonList[j]+1].nearest, protonList[j], ifilelol)

                geomdihmat[protonList[i] - 1][protonList[j] - 1] = \
                    etanumpy.get_dihedral(protonList[i],
                        ifilelol[protonList[i] + 1].nearest,
                        ifilelol[protonList[j] + 1].nearest,
                        protonList[j], ifilelol)
                geomdihmat[protonList[j] - 1][protonList[i] - 1] = \
                    -geomdihmat[i][j]
            else:

        # print "0.0000000"

                geomdihmat[protonList[i] - 1][protonList[j] - 1] = 0.0
                geomdihmat[protonList[j] - 1][protonList[i] - 1] = 0.0

  # ##Print out these matrixes
  # for i in range(len(geomdistmat)):
  #  line = []
  #  for j in range(len(geomdistmat[i])):
  #    line.append(str(geomdistmat[i][j]) + " ")
  #  print line
  # for i in range(len(geomdihmat)):
  #  line = []
  #  for j in range(len(geomdihmat[i])):
  #    line.append(str(geomdihmat[i][j]) + " ")
  #  print line

  # Set up a list for compairing the distances and dihedrals

    geomdist = []
    geomdih = []

  # Sum up all the distances and dihedrals for future deduplication

    for i in range(len(permutationsDeduplicatedNoGeom)):
        dist = 0
        dih = 0

    # print permutationsDeduplicatedNoGeom[i]

        for j in range(1, len(permutationsDeduplicatedNoGeom[i])):

      # print str(permutationsDeduplicatedNoGeom[i][j-1]) + "  " + str(permutationsDeduplicatedNoGeom[i][j])

            dih += geomdihmat[permutationsDeduplicatedNoGeom[i][j - 1]
                              - 1][permutationsDeduplicatedNoGeom[i][j]
                                   - 1]
            dist += geomdistmat[permutationsDeduplicatedNoGeom[i][j
                                - 1]
                                - 1][permutationsDeduplicatedNoGeom[i][j]
                    - 1]
            if geomdihmat[permutationsDeduplicatedNoGeom[i][j - 1]
                          - 1][permutationsDeduplicatedNoGeom[i][j]
                               - 1] == 0:

        # print permutationsDeduplicatedNoGeom[i]

                permutationsDeduplicatedNoGeom[i] = [-1]
                break
        geomdist.append(dist)
        geomdih.append(dih)

  # ## --- make a new deduplicated list without the removed list entries (those which had [-1]) --- ###

    permutationsDeduplicatedSomeGeom = [x for x in
            permutationsDeduplicatedNoGeom if x != [-1]]

  # ## --- Look through the dist and dihedral lists and if they are within a cut-off they are degenerate. --- ###

    for i in range(len(permutationsDeduplicatedSomeGeom)):
        for k in range(i + 1, len(permutationsDeduplicatedSomeGeom)):

      # if abs(geomdist[i-1] - geomdist[i]) <= getError(geomdist[i-1], geomdist[i]) and abs(geomdih[i-1] - geomdih[i]) <= getError(geomdih[i-1], geomdih[i]):
      # if abs(geomdist[i-1] - geomdist[i]) <= 0.05*deprotonations and abs(geomdih[i-1] - geomdih[i]) <= 5*deprotonations:

            if abs(geomdist[i - 1] - geomdist[k - 1]) \
                <= getErrorDist(geomdist[i - 1], geomdist[k - 1]) \
                and abs(geomdih[i - 1] - geomdih[k - 1]) \
                <= getErrorDih(geomdih[i - 1], geomdih[k - 1]):
                permutationsDeduplicatedSomeGeom[k - 1] = [-1]
                geomdist[k - 1] = -1
                geomdih[k - 1] = -1
    geomdistDeduped = [x for x in geomdist if x != -1]
    geomdihDeduped = [x for x in geomdih if x != -1]
else:
    permutationsDeduplicatedSomeGeom = permutationsDeduplicatedNoGeom

### --- make a new deduplicated list without the removed list entries (those which had [-1]) --- ###

permutationsDeduplicated = [x for x in permutationsDeduplicatedSomeGeom
                            if x != [-1]]

#### --- List all the sums of distances and dihedrals for future removal of degenerate deprotonations --- ###
# if deprotonations > 1:
#  geomdist = []
#  geomdih = []
#  for i in range(len(permutationsDeduplicatedNoGeom)):
#    dist = 0
#    dih = 0
#    for j in range(len(permutationsDeduplicatedNoGeom[i])):
#      if ifilelol[permutationsDeduplicatedNoGeom[i][j - 1] + 1].nearest + 1 != ifilelol[
#              permutationsDeduplicatedNoGeom[i][j] + 1].nearest + 1:
#        dih += etanumpy.get_dihedral(permutationsDeduplicatedNoGeom[i][j - 1] + 1,
#                           ifilelol[permutationsDeduplicatedNoGeom[i][j - 1] + 1].nearest + 1,
#                           ifilelol[permutationsDeduplicatedNoGeom[i][j] + 1].nearest + 1,
#                           permutationsDeduplicatedNoGeom[i][j] + 1, ifilelol)
#      else:
#        permutationsDeduplicatedNoGeom[i] = [-1]
#        break
#      dist += etanumpy.get_distance(permutationsDeduplicatedNoGeom[i][j - 1] + 1, permutationsDeduplicatedNoGeom[i][j] + 1,
#                          ifilelol)
#    geomdist.append(dist)
#    geomdih.append(dih)
#
#  ### --- Look through the dist and dihedral lists and if they are within a cut-off they are degenerate. --- ###
#  for i in range(len(permutationsDeduplicatedNoGeom)):
#    for k in range(i + 1, len(permutationsDeduplicatedNoGeom)):
#      #if abs(geomdist[i-1] - geomdist[i]) <= getError(geomdist[i-1], geomdist[i]) and abs(geomdih[i-1] - geomdih[i]) <= getError(geomdih[i-1], geomdih[i]):
#      #if abs(geomdist[i-1] - geomdist[i]) <= 0.05*deprotonations and abs(geomdih[i-1] - geomdih[i]) <= 5*deprotonations:
#      if abs(geomdist[i - 1] - geomdist[k - 1]) <= getErrorDist(geomdist[i - 1], geomdist[k - 1]) and abs(
#                  geomdih[i - 1] - geomdih[k - 1]) <= getErrorDih(geomdih[i - 1], geomdih[k - 1]):
#        permutationsDeduplicatedNoGeom[k - 1] = [-1]
#        geomdist[k - 1] = -1
#        geomdih[k - 1] = -1
#  geomdistDeduped = [x for x in geomdist if x != -1]
#  geomdihDeduped = [x for x in geomdih if x != -1]

### --- Look to see that two hydrogens are not on the same oxegen --- ###
# ....for i in range(len(permutationsDeduplicatedNoGeom)):
# ........if permutationsDeduplicatedNoGeom[i][0] != -1:
# ............if ifilelol[permutationsDeduplicatedNoGeom[i][0]+1].nearest == ifilelol[permutationsDeduplicatedNoGeom[i][1]+1].nearest:
# ................permutationsDeduplicatedNoGeom[j] = [-1]

#### --- make a new deduplicated list without the removed list entries (those which had [-1]) --- ###
# permutationsDeduplicated = [x for x in permutationsDeduplicatedNoGeom if x != [-1]]

# Geom dedup time elapsed

elapsed_time = time.time() - start_time
print 'Geom dedup time elapsed: ' + str(elapsed_time)

### --- Print out data regarding the number of permutations (non degenerate states) --- ###

print 'Permutations (total: ' + str(len(permutationsDeduplicated)) \
    + ') for deprotonated geometries:'
if verbose == 1:
    print permutationsDeduplicated
print 'Output deprotonated files in ' + str(outputFolder) + '.'

### --- Output the permutations in xyz files --- ###

parentFileAtomNum = ifilelol[0]
newFileAtomNum = parentFileAtomNum - deprotonations
for i in range(0, len(permutationsDeduplicated)):
    newFileName = str(outputFolder) + '-' + str(deprotonations) + '-' \
        + str(i + 1).zfill(len(str(len(permutationsDeduplicated)))) \
        + '.xyz'
    f = open(outputFolder + '/' + newFileName, 'w+')
    f.write(str(newFileAtomNum) + '\n')
    f.write(str(newFileName) + '\n')
    if verbose == 1 or debug == 1:
        print 'Permutation: ' + str(i + 1)
        print 'Removing atom(s): ' + str(permutationsDeduplicated[i])
        if deprotonations > 1:
            print str(geomdistDeduped[i]) + ' ' + str(geomdihDeduped[i])
    linetemp = []
    for j in range(0, int(parentFileAtomNum)):
        if j + 1 not in permutationsDeduplicated[i]:
            linetemp.append(ifilelol[j + 2].e)
            linetemp.append(ifilelol[j + 2].x)
            linetemp.append(ifilelol[j + 2].y)
            linetemp.append(ifilelol[j + 2].z)
            line = '    '.join(str(e) for e in linetemp) + '\n'
            linetemp[:] = []
            f.write(line)
        else:
            if debug == 1:
                linetemp.append(ifilelol[j + 2].e)
                linetemp.append(ifilelol[j + 2].x)
                linetemp.append(ifilelol[j + 2].y)
                linetemp.append(ifilelol[j + 2].z)
                print '    '.join(str(e) for e in linetemp)
                linetemp[:] = []
    f.close()

# total time elapsed

elapsed_time = time.time() - start_time
print 'Total time elapsed: ' + str(elapsed_time)

######################################################################
### END OF SCRIPT
######################################################################
