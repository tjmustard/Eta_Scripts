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
from decimal import *

### --- Arguments --- ###

program = 'nTranslate-XYZ-Remove_Atom.py'
ifile = ''
ofile = ''
deleteAtoms = []
allindir = 0

### Read command line args

try:
    (myopts, args) = getopt.getopt(sys.argv[1:], 'i:o:a:dh')
except getopt.GetoptError:
    print program \
        + ' -i <inputfile.xyz> -o <outputfile.xyz> -a <delete_atom> -d (all in dir)'
    sys.exit(2)

###############################
# o == option
# a == argument passed to the o
###############################

for (o, a) in myopts:
    if o == '-i':
        ifile = a
    elif o == '-o':
        ofile = a
    elif o == '-a':
        deleteAtoms.append(int(a))
    elif o == '-d':
        allindir = 1
    elif o == '-h':
        print program \
            + ' -i <inputfile.xyz> -o <outputfile.xyz> -d <delete_atom>'
        sys.exit(0)
    else:
        print 'Usage: %s -i <inputfile.xyz> -o <outputfile.xyz> -d <delete_atom>' \
            % sys.argv[0]
        sys.exit(0)

if allindir == 0:

  # ## --- Open parent file --- ###

    ifilelol = etaatom.xyz_lol(ifile)

  # ## --- Output file in Z-matrix format --- ###

    f = open(ofile, 'w+')
    f.write(str(ifilelol[0] - len(deleteAtoms)) + '\n')
    f.write(ifilelol[1] + '\n')
    for i in range(2, len(ifilelol)):
        if i - 1 not in deleteAtoms:
            line = ifilelol[i].e + '  ' \
                + str('{:.6f}'.format(ifilelol[i].x)) + '  ' \
                + str('{:.6f}'.format(ifilelol[i].y)) + '  ' \
                + str('{:.6f}'.format(ifilelol[i].z))
            f.write(line + '\n')
        else:
            print 'Removing atom: ' + str(i - 1)
    f.close()
elif allindir == 1:

    etaatom.make_dir('REMOVED')
    for i in os.listdir(os.getcwd()):
        if i.endswith('.xyz'):
            ifile = i
            ofile = 'REMOVED/' + ifile
            print i

      # ## --- Open parent file --- ###

            ifilelol = etaatom.xyz_lol(ifile)

      # ## --- Output file in Z-matrix format --- ###

            f = open(ofile, 'w+')
            f.write(str(ifilelol[0] - len(deleteAtoms)) + '\n')
            f.write(ifilelol[1] + '\n')
            for j in range(2, len(ifilelol)):
                if j - 1 not in deleteAtoms:
                    line = ifilelol[j].e + '  ' \
                        + str('{:.6f}'.format(ifilelol[j].x)) + '  ' \
                        + str('{:.6f}'.format(ifilelol[j].y)) + '  ' \
                        + str('{:.6f}'.format(ifilelol[j].z))
                    f.write(line + '\n')
                else:
                    print 'Removing atom: ' + str(j - 1)
            f.close()

######################################################################
### END OF SCRIPT
######################################################################
