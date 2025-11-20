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

import getopt
import math
import os
import sys
from sys import *

### --- Arguments --- ###

program = "nTM-MEMORY.py"

# Grab the first argument from the command and use that as the snippet

try:
    controlfile = sys.argv[1]
    memory = int(sys.argv[2])
    cpus = int(sys.argv[3])
    ratio = sys.argv[4]
except IndexError:
    controlfile = "-h"
    memory = "-h"

# If help is wanted allow the skipping of a snippet

if memory == "-h":
    print(
        program
        + """ <control file> <memory> <cpus> <ratio= genral or freq> -h (this)"""
    )
    sys.exit(0)

if ratio == "general":
    ratiori = 1
    ratiomax = 0.9
    ratiorpa = 0.1
elif ratio == "freq":
    ratiori = 1
    ratiomax = 0.9
    ratiorpa = 0.1
else:
    print(program + """ <control file> <memory> <cpus> <ratio> -h (this)""")
    sys.exit(0)

f = open(controlfile, "r")
controllines = f.readlines()
f.close()
f = open(controlfile, "w")
for (i, line) in enumerate(controllines):
    if "$ricore" in line:
        f.write("$ricore " + str(memory / cpus * 1) + "\n")
    elif "$maxcor" in line:
        f.write("$maxcor " + str(memory * 0.9) + "\n")
    elif "$rpacor" in line:
        f.write("$rpacor " + str(memory * 0.1) + "\n")
    else:
        f.write(line)
f.close()

# #####################################################################
# ## END OF SCRIPT
# #####################################################################
