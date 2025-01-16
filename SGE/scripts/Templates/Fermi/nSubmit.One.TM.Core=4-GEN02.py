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

### --- Arguments --- ###

program = 'Submit.One.TM.Core=2-GEN02.py'

# Grab the first argument from the command and use that as the snippet

try:
    inputfile = sys.argv[1]
except IndexError:
    inputfile = '-h'

# If help is wanted allow the skipping of a snippet

if inputfile == '-h':
    argv.append('-h')
    argv.append('-h')
elif inputfile == '-d':
    argv.append('-d')
    argv.append('-d')
elif inputfile == '--hold':
    argv.append('--hold')
    argv.append('--hold')
debug = 0
hold = 0

### Read command line args

try:
    (myopts, args) = getopt.getopt(sys.argv[2:], 'dh', ['hold'])
except getopt.GetoptError:
    print program + ''' <TM input file> -d --hold'''
    sys.exit(2)

###############################
# o == option
# a == argument passed to the o
###############################

for (o, a) in myopts:
    if o == '-d':
        debug += 1
    elif o == '--hold':
        hold = 1
    elif o == '-h':
        print program + ''' <TM input file> -d --hold'''
        sys.exit(0)
    else:
        print 'Usage: %s  <TM input file> -d --hold' % sys.argv[0]
        sys.exit(0)

if debug >= 1:
    print inputfile

### --- Get the 'basename' of the file --- ###
# Assuming the input file ends with '.in'
# basename = inputfile[::-1].replace('ni.', '')[::-1]

basename = etaatom.basename(inputfile, '.Turbomole')

#########################################################
### --- Job settings --- ###

program = 'TM'
nodeprocs = 8  # Number of procs on node
nodemem = 64000  # Total memory on node in MB
nproc = 4  # Number of procs to be used
memory = str(int(math.floor(nodemem / nodeprocs * nproc * 0.90)))
queue = 'gen02'  # Queue name
EtaDir = os.environ['ETADIR']

header = EtaDir + '/SGE/headers/Fermi.SMP.sget'
snippet = sys.argv[1]  # EtaDir + "/SGE/snippets/hidden/TM.sget"

#########################################################

headerLines = etaatom.return_modified_snippet(
    basename,
    program,
    queue,
    nproc,
    memory,
    header,
    )
snippetLines = etaatom.return_modified_snippet(
    basename,
    program,
    queue,
    nproc,
    memory,
    snippet,
    )

f = open(basename + '.sge', 'w')
for line in headerLines:
    f.write(line)
for line in snippetLines:
    f.write(line)
f.close()

if hold == 0:
    os.system('qsub ' + basename + '.sge |tee ' + basename + '.joblog')
else:
    print 'The job ' + basename + '.sge will be made but not submitted.'

######################################################################
### END OF SCRIPT
###################################################################
