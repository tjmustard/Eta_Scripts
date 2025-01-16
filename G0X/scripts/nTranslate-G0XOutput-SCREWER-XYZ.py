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
import time
import shutil

### --- Arguments --- ###

program = 'nTranslate-G0XOutput-SCREWER-XYZ.py'
ifile = ''
ofile = ''
allInDir = 0
trajectory = 0
geometry = 0
vibration = 1
negative = False
numAtoms = 0
finalenergy = 0
finalgeom = []
multiplyer = float(1.000)

# Grab the first argument from the command and use that as the input file

try:
    ifile = sys.argv[1]
except IndexError:
    ifile = '-h'

# If help is wanted allow the skipping of a input file

if ifile == '-h':
    argv.append('-h')
    argv.append('-h')
if ifile == '-d':
    argv.append('-d')
    argv.append('-d')

### Read command line args

try:
    (myopts, args) = getopt.getopt(sys.argv[2:], 'dnv:m:h')
except getopt.GetoptError:
    print program \
        + ' <input file (G0X log)> -v (vib. number) -d (all in dir) -n (negative vibration vectors) -m <multiplyer>'
    sys.exit(2)

###############################
# o == option
# a == argument passed to the o
###############################

for (o, a) in myopts:
    if o == '-v':
        vibration = int(a)
    elif o == '-d':
        allInDir = 1
    elif o == '-n':
        negative = True
    elif o == '-m':
        multiplyer = float(a)
    elif o == '-h':
        print program \
            + ' <input file (G0X log)> -v (vib. number) -d (all in dir) -n (negative vibration vectors) -m <multiplyer>'
        sys.exit(0)
    else:
        print 'Usage: %s  <input file (G0X log)> -v (vib. number) -d (all in dir) -n (negative vibration vectors) -m <multiplyer>' \
            % sys.argv[0]
        sys.exit(0)

# Make SCREWERed directory

etaatom.make_dir('SCREWERed')


def G0X_Screwer(ifile):

  # Print some information to the screen

    print 'A new XYZ file will be generated from ' + ifile
    print 'and will be augmented by the vibrational frequency ' \
        + str(vibration) + '.'

  # Extract the last geometry from the G0X log file

    (ifilelol, charge, multi) = etaatom.parse_output_g0x(ifile)

  # Grab the xyz changes for the SCREWERing and place them all into one big list

    freqswitch = 0
    freqwrite = 0
    writecount = -1
    freqlist = []
    the_file = open(ifile, 'r')
    for (idx, line) in enumerate(the_file):

    # print line

        if line.strip() != '':
            if line.split()[0] == 'NAtoms=':
                numAtoms = int(line.split()[1])
            if 'Low frequencies' in line:
                freqswitch = 1
            elif line.split()[0] == 'Frequencies' and freqswitch == 1:
                for i in range(numAtoms * 3):
                    freqlist.append(0.0000)
            elif line.split()[0] == 'Atom' and freqswitch == 1:
                freqwrite = 1
            elif freqwrite == 1:
                if len(line.split()) < 11:
                    freqwrite = 0
                    writecount += numAtoms * 3
                else:

                    freqlist[writecount + int(line.split()[0])
                             + numAtoms * 0] = line.split()[2] + ' ' \
                        + line.split()[3] + ' ' + line.split()[4]

                    freqlist[writecount + int(line.split()[0])
                             + numAtoms * 1] = line.split()[5] + ' ' \
                        + line.split()[6] + ' ' + line.split()[7]

                    freqlist[writecount + int(line.split()[0])
                             + numAtoms * 2] = line.split()[8] + ' ' \
                        + line.split()[9] + ' ' + line.split()[10]

  # Generate ofile name

    ofiletmp = ifile.split('.')
    ofile = ''
    for i in range(len(ofiletmp)):
        if i != len(ofiletmp) - 1:
            ofile += ofiletmp[i] + '.'
        elif i == len(ofiletmp) - 1:
            ofile += 'xyz'
    print ifile + ' ----> ' + ofile

  # Move the X Y and Z coordinates with regard to the vibrational mode of interest

    for i in range(2, len(ifilelol)):
        ifilelol[i].x += multiplyer * float(freqlist[(vibration - 1)
                * numAtoms + i - 2].split()[0])
        ifilelol[i].y += multiplyer * float(freqlist[(vibration - 1)
                * numAtoms + i - 2].split()[1])
        ifilelol[i].z += multiplyer * float(freqlist[(vibration - 1)
                * numAtoms + i - 2].split()[2])

  # Ouput the geometry(ies) to the ofile

    etaatom.output_xyz('SCREWERed/' + ofile, ifilelol)
    return


if negative:
    multiplyer *= -1
else:
    multiplyer = multiplyer

if allInDir == 1:
    for i in os.listdir(os.getcwd()):
        if i.endswith('.log'):
            ifile = i
            G0X_Screwer(ifile)
else:
    G0X_Screwer(ifile)

  # #####################################################################
  # ## END OF SCRIPT
  # #####################################################################
