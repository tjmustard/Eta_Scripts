#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright (c) 2025, Thomas J. L. Mustard, O. Maduka Ogba, Paul Ha-Yeon Cheong
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
from decimal import *
import shutil
import etaatom
import etamap

### --- Arguments --- ###

program = 'nMAP.py'
ifile = ''
ofile = ''
mapinput = etamap.MAPArguments()
mapinput.debug = 0
writetemplatemap = False
triangle = False

# Read command line args

try:
    (myopts, args) = getopt.getopt(sys.argv[1:], 'i:tThD')
except getopt.GetoptError:
    print(program \
        + ' -i <inputfile.map> -t <print out template map file> -T <only build the triangle of the posible combinations>')
    sys.exit(2)

###############################
# o == option
# a == argument passed to the o
###############################

for (o, a) in myopts:
    if o == '-i':
        ifile = a
    elif o == '-D':
        mapinput.debug += 1
    elif o == '-t':
        writetemplatemap = True
    elif o == '-T':
        triangle = True
    elif o == '-h':
        print(program \
            + ' -i <inputfile.map> -t <print out template map file> -T <only build the triangle of the posible combinations>')
        sys.exit(0)
    else:
        print('Usage: %s -i inputfile.map -t <print out template map file> -T <only build the triangle of the posible combinations>' \
            % sys.argv[0])
        sys.exit(0)

if writetemplatemap == 1:
    etamap.write_template()
    sys.exit(0)

mapinput = etamap.parse_input_file(ifile)

### --- Expand upon the parentDeleteTemp list if
# there are certain characters to denote multiple atoms.
###

mapinput = etamap.expand_temp_lists(mapinput)

### --- Make the output folder --- ###

etaatom.make_dir(mapinput.outputFolder)
etaatom.make_dir(mapinput.outputFolder + '/temp')

### --- Make the align files for quatfit. --- ###

etamap.make_align_file(mapinput)

### --- Align the parent files before aligning the substrates to them --- ###
### Purely for aesthetics and does not change the final outcome

etamap.align_parents(mapinput)

etamap.align_children(mapinput, triangle)

######################################################################
### END OF SCRIPT
######################################################################
