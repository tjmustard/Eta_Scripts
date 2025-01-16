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

### --- Arguments --- ###

program = 'nTranslate-TMOut-XYZ.py'
ifile = ''
ofile = ''
allInDir = 0
trajectory = 0

numAtoms = 0
finalenergy = 0
finalgeom = []

# Grab the first argument from the command and use that as the input file

try:
    ifile = sys.argv[1]
except IndexError:
    ifile = '-h'

# If help is wanted allow the skipping of a input file

if ifile == '-h':
    argv.append('-h')
    argv.append('-h')
elif ifile == '-d':
    argv.append('-d')
elif ifile == '-t':
    argv.append('-t')

### Read command line args

try:
    (myopts, args) = getopt.getopt(sys.argv[2:], 'dtg:h')
except getopt.GetoptError:
    print program \
        + ' <input file (TM log)> -d (all in current directory) -t (output entire trajectory) -g ## (export geom number ##)'
    sys.exit(2)

###############################
# o == option
# a == argument passed to the o
###############################

for (o, a) in myopts:
    if o == '-d':
        allInDir = 1
    elif o == '-t':
        trajectory = 1
    elif o == '-g':
        geometry = int(a)
        trajectory = 2
    elif o == '-h':
        print program \
            + ' <input file (TM log)> -d (all in current directory) -t (output entire trajectory) -g ## (export geom number ##)'
        sys.exit(0)
    else:
        print 'Usage: %s  <input file (TM log)> -d (all in current directory) -t (output entire trajectory) -g ## (export geom number ##)' \
            % sys.argv[0]
        sys.exit(0)

if allInDir == 1:
    for i in os.listdir(os.getcwd()):
        if i.endswith('.log'):
            ifile = i
            ofile = etaatom.basename(ifile, '.log') + '.xyz'
            print ifile + ' ----> ' + ofile
            if trajectory == 1:

        # Extract the entire trajectory from the Turbomole log file

                (ifilelol, charge, multi) = \
                    etaatom.parse_output_tm_traj(ifile)
            elif trajectory == 2:

        # Extract the entire trajectory from the G0X log file

                (trajectory, charge, multi) = \
                    etaatom.parse_output_tm_traj(ifile)
                ifilelol = etaatom.get_geom_n(trajectory, geometry)
            else:

        # Extract the last geometry from the Turbomole log file

                (ifilelol, charge, multi) = \
                    etaatom.parse_output_tm(ifile)

      # Ouput the geometry(ies) to the ofile

            etaatom.output_xyz(ofile, ifilelol)
else:
    if trajectory == 1:

    # Extract the last geometry from the Turbomole log file

        (ifilelol, charge, multi) = etaatom.parse_output_tm_traj(ifile)
    elif trajectory == 2:

    # Extract the entire trajectory from the G0X log file

        (trajectory, charge, multi) = \
            etaatom.parse_output_tm_traj(ifile)
        ifilelol = etaatom.get_geom_n(trajectory, geometry)
    else:

    # Extract the last geometry from the Turbomole log file

        (ifilelol, charge, multi) = etaatom.parse_output_tm(ifile)
    ofiletmp = ifile.split('.')
    ofile = ''
    for i in range(len(ofiletmp)):
        if i != len(ofiletmp) - 1:
            ofile += ofiletmp[i] + '.'
        elif i == len(ofiletmp) - 1:
            ofile += 'xyz'
    print ifile + ' ----> ' + ofile

  # Ouput the geometry(ies) to the ofile

    etaatom.output_xyz(ofile, ifilelol)

  # #####################################################################
  # ## END OF SCRIPT
  # #####################################################################
