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

program = 'nTranslate-XYZ-G09Input.py'

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
multiset = False
modred = []
writemod = 0
ionicCM = ''
ionicswitch = 0
debug = 0
interactive = 0

### Read command line args

try:
    (myopts, args) = getopt.getopt(sys.argv[2:], 'c:m:dih', ['mod=',
                                   'ionic='])
except getopt.GetoptError:
    print program \
        + ''' <snippet> -c <charge> -m <multiplicity> --mod=\"constraint\" --ionic=reg/mno -d'''
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
    elif o == '--mod':
        modred.append(a)
        writemod = 1
    elif o == '--ionic':
        ionicCM = a
        ionicswitch = 1
    elif o == '-i':
        interactive = 1
    elif o == '-h':
        print program \
            + ''' <snippet> -c <charge> -m <multiplicity> --mod=\"constraint\" --ionic=reg/mno -d'''
        sys.exit(0)
    else:
        print 'Usage: %s  <snippet> -c <charge> -m <multiplicity> --mod="constraint" --ionic=reg/mno -d' \
            % sys.argv[0]
        sys.exit(0)

if debug >= 1:
    print 'Charge: ' + str(charge)
    print 'Multiplicity: ' + str(multi)
    print 'Snippet: ' + snippet
    print 'Modredundant Variable: '
    print modred

### --- OK Now setup the input variables for the first time --- ###

G0XInput = etaatom.InputArguments()

# if interactive == 1:
#  G0XInput = etaatom.interactive_config
# else:
#  #Grab the snippet information
#  G0XInput = etaatom.parse_snippet_g0x(snippet, G0XInput)

### --- Make NEWJOBS folder --- ###

if not os.path.exists('NEWJOBS.G0X'):
    os.makedirs('NEWJOBS.G0X')

### --- Open and parse the xyz files in the folder --- ###

print 'Currently Translating:'
for i in os.listdir(os.getcwd()):
    if i.endswith('.xyz'):
        ifile = i
        print i

        if debug >= 2:
            print snippet

    # Parse XYZ file into a list of lists

        ifilelol = etaatom.xyz_lol(ifile)

    # ## Copy the charge, multi and modredundant lines to the G0XInput

        if len(modred) > 0:
            G0XInput.writemod = 1
        for line in modred:
            G0XInput.modred.append(line)
        if chargeset:
            print 'YES'
            G0XInput.charge = charge
        if multiset:
            G0XInput.multi = multi

    # ## --- Calculate charge from file --- ###

        if ionicswitch == 1:
            G0XInput.charge = etaatom.ionic(ionicswitch, ionicCM,
                    ifilelol)

    # ## --- Generate the ECP list for files if needed --- ###

        if interactive == 0:
            G0XInput = etaatom.parse_snippet_g0x(snippet, G0XInput)

        if G0XInput.writeecp == 1:
            (G0XInput.config, G0XInput.ecplines) = \
                etaatom.g0x_ecp(G0XInput.ecp, G0XInput.config, ifilelol)

        if G0XInput.writecloseecp == 1:
            (G0XInput.closelines, G0XInput.closeecplines) = \
                etaatom.g0x_ecp(G0XInput.closeecp, G0XInput.closelines,
                                ifilelol)

    # ## --- Out put the G0X input files --- ###

        etaatom.output_g0x(ifile, G0XInput, ifilelol)

    # ## --- RESET THE VARIABLES --- ###

        G0XInput = etaatom.reset_input_variables(G0XInput)

    # #####################################################################
    # ## END OF SCRIPT
    # #####################################################################
