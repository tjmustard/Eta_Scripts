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

import os
import math
from sys import *
import sys
import getopt

### --- Arguments --- ###

program = 'nTM-TS.py'

# Grab the first argument from the command and use that as the snippet

try:
    controlfile = sys.argv[1]
except IndexError:
    controlfile = '-h'

# If help is wanted allow the skipping of a snippet

if controlfile == '-h':
    print program + ''' <control file> <ts = 1 opt = 0> -h (this)'''
    sys.exit(0)

if sys.argv[2] == str(0) or sys.argv[2] == str(1):
    yep = 0
else:
    print 'You must specify 1 for TS or 0 for OPT.'
    sys.exit(0)

f = open(controlfile, 'r')
controllines = f.readlines()
f.close()
f = open(controlfile, 'w')
for (i, line) in enumerate(controllines):
    if 'itrvec' in line:
        f.write('   itrvec ' + sys.argv[2] + '\n')
    else:
        f.write(line)
f.close()

  # #####################################################################
  # ## END OF SCRIPT
  # #####################################################################
