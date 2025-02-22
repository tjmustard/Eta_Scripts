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

import os
import math
from sys import *
import sys
import getopt
import etaatom
import time
import shutil

### --- Arguments --- ###

program = 'nTranslate-G0XOutput-FREQ_List.py'
ifile = ''
ofile = ''
allindir = 0
numatoms = 0

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
    (myopts, args) = getopt.getopt(sys.argv[2:], 'dh')
except getopt.GetoptError:
    print program + ' <input file (G0X log)> -d (all in dir)'
    sys.exit(2)

###############################
# o == option
# a == argument passed to the o
###############################

for (o, a) in myopts:
    if o == '-d':
        allindir = 1
    elif o == '-h':
        print program + '  <input file (G0X log)> -d (all in dir)'
        sys.exit(0)
    else:
        print 'Usage: %s <input file (G0X log)> -d (all in dir)' \
            % sys.argv[0]
        sys.exit(0)


def g0x_vibrations(ifile):

  # Print some information to the screen

    print ifile + ' will be parsed to find all the vibrations.'

  # print "and will be augmented by the vibrational frequency " + str(vibration) + "."

  # Grab the xyz changes for the SCREWERing and place them all into one big list

    freqswitch = False
    freqlist = []
    the_file = open(ifile, 'r')
    for (idx, line) in enumerate(the_file):

    # print line

        if line.strip() != '':
            if line.split()[0] == 'NAtoms=':
                numatoms = int(line.split()[1])
            if 'Low frequencies' in line:
                freqswitch = True
            elif line.split()[0] == 'Frequencies' and freqswitch:

        # Store the frequencies

                freqlist.append(float(line.split()[2]))
                if len(line.split()) >= 4:
                    freqlist.append(float(line.split()[3]))
                if len(line.split()) >= 5:
                    freqlist.append(float(line.split()[4]))

    for vib in freqlist:
        print vib


if allindir == 1:
    for i in os.listdir(os.getcwd()):
        if i.endswith('.log'):
            ifile = i
            g0x_vibrations(ifile)
else:
    g0x_vibrations(ifile)

  # #####################################################################
  # ## END OF SCRIPT
  # #####################################################################
