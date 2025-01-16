#!/usr/bin/env python2.7
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
import etaatom

### --- Make NEWJOBS folder --- ###

if not os.path.exists('ARCHIVE'):
    os.makedirs('ARCHIVE')
if not os.path.exists('ARCHIVE/PUB'):
    os.makedirs('ARCHIVE/PUB')
if not os.path.exists('ARCHIVE/TABLE'):
    os.makedirs('ARCHIVE/TABLE')
f = open('ARCHIVE/TABLE/Table-Energy-Equationed.txt', 'w')
f.write('File_Name,SCF,ZPE,T*S,H,G,dG(Hart),dG(kcal/mol)\n')
f.close()
f = open('ARCHIVE/TABLE/Table-Energy-Totals.txt', 'w')
f.write('File_Name,SCF,ZPE,T*S,H,G,dG(Hart),dG(kcal/mol)\n')
f.close()
if not os.path.exists('ARCHIVE/XYZ'):
    os.makedirs('ARCHIVE/XYZ')
if not os.path.exists('ARCHIVE/JOBLOG'):
    os.makedirs('ARCHIVE/JOBLOG')
if not os.path.exists('ARCHIVE/LOG'):
    os.makedirs('ARCHIVE/LOG')
if not os.path.exists('ARCHIVE/INPUT'):
    os.makedirs('ARCHIVE/INPUT')

pwd = os.getcwd()

# Archive each log file

for i in os.listdir(pwd):
    if i.endswith('.out'):
        basefile = i.replace('.out', '', 1)
        logfile = i
        print 'Archiving ' + basefile

    # Grab the final geometry and the translated output for the pub file

        (finalgeom, charge, multi) = etaatom.parse_output_g0x(logfile)
        (printLines, tableLines) = \
            etaatom.translate_g0x_output(logfile, 'archive', 'none')

    # Output the pub files

        f = open('ARCHIVE/PUB/' + basefile + '.pub', 'w')
        for line in printLines:
            if 'GEOM DATA!' in line:
                for j in range(len(finalgeom)):

          # if i <= 1:
          #  print finalgeom[i]

                    if j >= 2:
                        f.write(etaatom.get_element_name(etaatom.get_element_num(str(finalgeom[j].e)))
                                + '      ' + str(finalgeom[j].x)
                                + '      ' + str(finalgeom[j].y)
                                + '      ' + str(finalgeom[j].z) + '\n')
            else:
                f.write(line + '\n')
        f.close()

    # Output the Table files

        for (j, line) in enumerate(tableLines):
            if j == 0:
                f = open('ARCHIVE/TABLE/Table-Energy-Equationed.txt',
                         'a')
                f.write(line + '\n')
                f.close()
            elif j == 1:
                f = open('ARCHIVE/TABLE/Table-Energy-Totals.txt', 'a')
                f.write(line + '\n')
                f.close()
            elif j == 2:
                f = open('ARCHIVE/TABLE/Table-Energy-SCF.txt', 'a')
                f.write(line + '\n')
                f.close()

    # Output the final geometry in xyz format

        etaatom.output_xyz('ARCHIVE/XYZ/' + basefile + '.xyz',
                           finalgeom)

  # #####################################################################
  # ## END OF SCRIPT
  # #####################################################################
