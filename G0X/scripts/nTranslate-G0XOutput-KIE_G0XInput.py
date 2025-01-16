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
from sys import *
import sys
import getopt
import etaatom

### --- Arguments --- ###

program = 'nTranslate-G0XOutput-KIE_g0xinput.py'

# Grab the first argument from the command and use that as the snippet

try:
    snippet = sys.argv[1]
except IndexError:
    snippet = '-h'

# If help is wanted allow the skipping of a snippet

if snippet == '-h':
    argv.append('-h')
    argv.append('-h')

# If interactive mode is wanted allow the skipping of a snippet

if snippet == '-i':
    argv.append('-i')
    argv.append('-i')

# Unless the charge is stated via the -c argument the charge will be assumed to be 0

charge = 0
chargeset = False

# Unless the multiplicity is stated via the -m argument the multiplicity will be assumed to be 1

multi = 1
temp = 298.15
pressure = 1.0
multiset = False
modred = []
writemod = 0
ionicCM = ''
ionicswitch = 0
debug = 0
interactive = 0

### Read command line args

try:
    (myopts, args) = getopt.getopt(sys.argv[2:], 'c:m:t:p:dih')
except getopt.GetoptError:
    print program \
        + ''' <snippet> -c <charge> -m <multiplicity> -t <temp in kelvin> -p <pressure> -d'''
    sys.exit(2)

###############################
# o == option
# a == argument passed to the o
###############################

for (o, a) in myopts:
    if o == '-c':
        charge = a
        chargeset = True
    elif o == '-m':
        multi = a
        multiset = True
    elif o == '-d':
        debug += 1
    elif o == '-t':
        temp = float(a)
    elif o == '-p':
        pressure = float(a)
    elif o == '-h':
        print program \
            + ''' <snippet> -c <charge> -m <multiplicity> -t <temp in kelvin> -p <pressure> -d'''
        sys.exit(0)
    else:
        print 'Usage: %s <snippet> -c <charge> -m <multiplicity> -t <temp in kelvin> -p <pressure> -d' \
            % sys.argv[0]
        sys.exit(0)

if debug >= 1:
    print 'Charge: ' + str(charge)
    print 'Multiplicity: ' + str(multi)
    print 'Snippet: ' + snippet
    print 'Temp: ' + temp
    print 'Pressure: ' + pressure

### --- OK Now setup the input variables for the first time --- ###

g0xinput = etaatom.InputArguments()

# if interactive == 1:
#  g0xinput = etaatom.interactive_config
# else:
#  #Grab the snippet information
#  g0xinput = etaatom.parse_snippet_g0x(snippet, g0xinput)

### --- Make NEWJOBS folder --- ###

if not os.path.exists('NEWJOBS.G0X'):
    os.makedirs('NEWJOBS.G0X')

### --- Open and parse the xyz files in the folder --- ###

print 'Currently Translating:'
for i in os.listdir(os.getcwd()):
    if i.endswith('.log'):
        ifile = i
        print i

        basename = etaatom.basename(ifile, '.log')

        if debug >= 2:
            print snippet

    # Parse XYZ file into a list of lists

        (ifilelol, g0xinput.charge, g0xinput.multi) = \
            etaatom.parse_output_g0x(ifile)

    # Change the charge if it has been set by the user

        if chargeset:
            g0xinput.charge = charge
        if multiset:
            g0xinput.multi = multi

    # ## --- Generate the configurationa and ECP list for files if needed --- ###

        g0xinput = etaatom.parse_snippet_g0x(snippet, g0xinput)

        if g0xinput.writeecp == 1:
            (g0xinput.config, g0xinput.ecplines) = \
                etaatom.g0x_ecp(g0xinput.ecp, g0xinput.config, ifilelol)

        if g0xinput.writecloseecp == 1:
            (g0xinput.closelines, g0xinput.closeecplines) = \
                etaatom.g0x_ecp(g0xinput.closeecp, g0xinput.closelines,
                                ifilelol)

    # ## --- Make a KIE output for each C and H in the system --- ###

        g0xinput.writemod = 1

        changeiso = []
        for i in range(2, len(ifilelol)):
            if ifilelol[i].en == 1 or ifilelol[i].en == 6:
                changeiso.append([i - 1, ifilelol[i].en])

        for i in range(len(changeiso)):
            g0xinput.modred = []
            g0xinput.modred.append(str(temp) + ' ' + str(pressure))
            for j in range(2, len(ifilelol)):
                if j == changeiso[i][0] + 1:
                    g0xinput.modred.append(str(etaatom.elementIsotope[ifilelol[j].en
                            - 1][1]))
                else:
                    g0xinput.modred.append(str(etaatom.elementIsotope[ifilelol[j].en
                            - 1][0]))

      # ## --- Out put the G0X input files --- ###

            etaatom.output_g0x(etaatom.get_element_name(changeiso[i][1])
                               + str(changeiso[i][0]) + '-' + ifile,
                               g0xinput, ifilelol)

            etaatom.copy_file(basename + '.chk', 'NEWJOBS.G0X/'
                              + etaatom.get_element_name(changeiso[i][1])
                              + str(changeiso[i][0]) + '-' + basename
                              + '.chk')

    # ## --- RESET THE VARIABLES --- ###

        g0xinput = etaatom.reset_input_variables(g0xinput)

######################################################################
### END OF SCRIPT
######################################################################

