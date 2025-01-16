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

program = 'nTranslate-XYZ-RMSD.py'
ifile = ''
ofile = ''
deleteAtoms = []

### Read command line args

try:
    (myopts, args) = getopt.getopt(sys.argv[1:], 'i:o:d:h')
except getopt.GetoptError:
    print program \
        + ' -i <inputfile.xyz> -o <outputfile.xyz> -d <delete_atom>'
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
    elif o == '-d':
        deleteAtoms.append(int(a))
    elif o == '-h':
        print program \
            + ' -i <inputfile.xyz> -o <outputfile.xyz> -d <delete_atom>'
        sys.exit(0)
    else:
        print 'Usage: %s -i <inputfile.xyz> -o <outputfile.xyz> -d <delete_atom>' \
            % sys.argv[0]
        sys.exit(0)

### --- Open parent file --- ###

ifilelol = etaatom.xyz_lol(ifile)

### --- Iterating through a folder of files --- ###

for i in os.listdir(os.getcwd()):
    if i.endswith('.xyz'):
        childFile = i

######################################################################
### END OF SCRIPT
######################################################################
