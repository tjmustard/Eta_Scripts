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

from EtaLib import etaatom

### --- Arguments --- ###

program = "nJT.py"

# Grab the first argument from the command and use that as the snippet

try:
    logfile = sys.argv[1]
except IndexError:
    logfile = "-h"

# If help is wanted allow the skipping of a snippet

if logfile == "-h":
    argv.append("-h")
    argv.append("-h")
if logfile == "-a":
    argv.append("-a")
    argv.append("-a")
debug = 0
printType = "none"
allInDir = 0

### Read command line args

try:
    (myopts, args) = getopt.getopt(sys.argv[2:], "dhela")
except getopt.GetoptError:
    print(program + " <log file> -d -e -l -a")
    sys.exit(2)

###############################
# o == option
# a == argument passed to the o
###############################

for (o, a) in myopts:
    if o == "-d":
        debug += 1
    elif o == "-l":
        printType = "long"
    elif o == "-e":
        printType = "extra"
    elif o == "-a":
        allInDir = 1
    elif o == "-h":
        print(program + """ <log file> -d -l -e -a""")
        sys.exit(0)
    else:
        print("Usage: %s  <log file> -d -l -e -a" % sys.argv[0])
        sys.exit(0)

if debug >= 1:
    print(logfile)


def parse():
    if allInDir == 1:
        for i in os.listdir(os.getcwd()):
            if i.endswith(".out"):
                logfile = i
                if debug >= 1:
                    print(logfile)
                print(
                    "##############################################################################"
                )
                printLines = etaatom.translate_jaguar_output(
                    logfile, "translate", printType
                )
                for line in printLines:
                    print(line)
    else:
        printLines = etaatom.translate_jaguar_output(logfile, "translate", printType)
        for line in printLines:
            print(line)
    return


try:
    parse()
except:
    if debug >= 1:
        parse()
    else:
        print("Something failed. Make sure the/a log file exists.")

# #####################################################################
# ## END OF SCRIPT
# #####################################################################
