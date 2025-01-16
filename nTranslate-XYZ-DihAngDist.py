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

if argv[1] == '-h':
    print '<input file>.xyz atom1 atom2 atom3 atom4'
    sys.exit(0)

### --- Arguments --- ###

ifile = argv[1]
atomTemp = ' '.join(sys.argv[2:] + ['-1', '-1', '-1', '-1'])
atoms = atomTemp.split()
atom1 = int(atoms[0])
atom2 = int(atoms[1])
atom3 = int(atoms[2])
atom4 = int(atoms[3])

### Display input and output file name passed as the args
# print ("Input file : %s and output file: %s" % (ifile,ofile) )

### --- Open parent file --- ###

ifilelol = etaatom.xyz_lol(ifile)

### --- Print out the distance, angle, or dihedral --- ###

if atom3 == -1:
    print ifile + ': ' + str(etanumpy.get_distance(atom1 + 1, atom2
                             + 1, ifilelol))
elif atom4 == -1:
    print ifile + ': ' + str(etanumpy.get_angle(atom1 + 1, atom2 + 1,
                             atom3 + 1, ifilelol))
elif atom4 >= 0:
    print ifile + ': ' + str(etanumpy.get_dihedral(atom1 + 1, atom2
                             + 1, atom3 + 1, atom4 + 1, ifilelol))

######################################################################
### END OF SCRIPT
######################################################################
