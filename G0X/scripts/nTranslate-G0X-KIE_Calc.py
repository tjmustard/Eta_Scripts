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

import getopt
import math
import os
import shutil
import sys

import numpy

from EtaLib import etaatom, etanumpy

### --- Arguments --- ###

program = "nTranslate-G0X-KIE_Calc.py"
tsnatfile = ""
tsisofile = ""
subnatfile = ""
subisofile = ""
numatoms = 0
debug = 0
temp = 298.15
scaling = 1.0000
method = ""

### Read command line args

try:
    (myopts, args) = getopt.getopt(
        sys.argv[1:], "t:s:c:S:dh", ["tsnat=", "tsiso=", "subnat=", "subiso="]
    )
except getopt.GetoptError:
    print(
        program
        + " --tsnat=<natural isotope TS log> --tsiso=<unnatural isotope TS log> \n"
        + "--subnat=<natural isotope sub log> --subiso=<unnatural isotope sub log>\n"
        + "-t <temp in K> -c <temp in C> -s <scalling> -S <method B3LYP, M06, etc>"
    )
    sys.exit(2)

###############################
# o == option
# a == argument passed to the o
###############################

for (o, a) in myopts:
    if o == "-d":
        debug += 1
    elif o == "-t":
        temp = float(a)
    elif o == "-c":
        temp = float(a) + 273.15
    elif o == "-s":
        scaling = float(a)
    elif o == "-S":
        method = a
    elif o == "--tsnat":
        tsnatfile = a
    elif o == "--tsiso":
        tsisofile = a
    elif o == "--subnat":
        subnatfile = a
    elif o == "--subiso":
        subisofile = a
    elif o == "-h":
        print(
            program
            + " --tsnat=<natural isotope TS log> --tsiso=<unnatural isotope TS log> \n"
            + "--subnat=<natural isotope sub log> --subiso=<unnatural isotope sub log>\n"
            + "-t <temp in K> -c <temp in C> -s <scalling> -S <method B3LYP, M06, etc>"
        )
        sys.exit(0)
    else:
        print(
            "Usage: %s --tsnat=<natural isotope TS log> --tsiso=<unnatural isotope TS log> \n"
            + "--subnat=<natural isotope sub log> --subiso=<unnatural isotope sub log>\n"
            + "-t <temp in K> -c <temp in C> -s <scalling> -S <method B3LYP, M06, etc>"
            % sys.argv[0]
        )
        sys.exit(0)

# Print out some information

print("\nComputing the KIE and quantum KIE from four G0X output files:")
print("- the TS w/ natural isotopes")
print("- the TS w/ an unnatural isotope")
print("- the Substrate w/ natural isotopes")
print("- the Substrate w/ an unnatural isotope")
print("")
print(
    "Make sure that the isotope replacement is identical for both the TS and Substrate.\n"
)

# If method is set, set the scaling factor

if method.upper() == "M062X":
    print("Using the scaling factor 0.961 for method M062X.\n")
    scaling = 0.961
elif method.upper() == "B3LYP":
    print("Using the scaling factor 0.961 for method B3LYP.\n")
    scaling = 0.963

# Run the KIE computation

(kie, qkie) = etanumpy.compute_kie(
    etaatom.g0x_vibrations(tsnatfile),
    etaatom.g0x_vibrations(tsisofile),
    etaatom.g0x_vibrations(subnatfile),
    etaatom.g0x_vibrations(subisofile),
    temp,
    scaling,
)
print("KIE:               KIE w/ Quantum Effects:")
print(str(kie) + "      " + str(qkie))

# #####################################################################
# ## END OF SCRIPT
# #####################################################################
