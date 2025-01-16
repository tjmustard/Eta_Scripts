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
import math
import sys
from decimal import *
import shutil
import etaatom
import quatfit


###==================================================================================================================###
### --- Input Arguments class --- ###

class MAPArguments(object):

    def __init__(self):

    # Hold the string for the PARENT folder

        self.parentFolder = ''

    # Hold the string for the CHILDREN folder

        self.childFolder = ''

    # Hold the string for the OUTPUT folder

        self.outputFolder = ''

    # Hold the list of align parent atoms

        self.alignParent = []

    # Hold the list of align children atoms

        self.alignChild = []

    # Hold the list of align weights

        self.alignWeights = []

    # Hold the list of delete parent atoms

        self.parentDelete = []

    # Hold the list of temp delete parent atoms

        self.parentDeleteTemp = []

    # Hold the list of delete children atoms

        self.childDelete = []

    # Hold the list of temp delete children atoms

        self.childDeleteTemp = []

    # Hold the location of quatfit

        self.quatfit = 'quatfit'

    # Hold the string for the basename

        self.baseName = ''

    # Hold the string for the endname

        self.endName = ''

    # Hold the base state for debug

        self.debug = 0


### --- Reset Input Arguments class --- ###

def reset_map_variables(mapinput):

  # Hold the string for the PARENT folder

    mapinput.parentFolder = ''

  # Hold the string for the CHILDREN folder

    mapinput.childFolder = ''

  # Hold the string for the OUTPUT folder

    mapinput.outputFolder = ''

  # Hold the list of align parent atoms

    mapinput.alignParent = []

  # Hold the list of align children atoms

    mapinput.alignChild = []

  # Hold the list of align weights

    mapinput.alignWeights = []

  # Hold the list of delete parent atoms

    mapinput.parentDelete = []

  # Hold the list of temp delete parent atoms

    mapinput.parentDeleteTemp = []

  # Hold the list of delete children atoms

    mapinput.childDelete = []

  # Hold the list of temp delete children atoms

    mapinput.childDeleteTemp = []

  # Hold the location of quatfit

    mapinput.quatfit = 'quatfit'

  # Hold the string for the basename

    mapinput.baseName = ''

  # Hold the string for the endname

    mapinput.endName = ''

  # Hold the base state for debug

    mapinput.debug = 0
    return mapinput


# Constants

auToAngstrom = 1.889725988579
hartreeTokcalmol = 627.509469
hartreeToeV = 27.21138386
kcalmolTokjoulemol = 4.184


###==================================================================================================================###

### --- Functions for Submission --- ###
###==================================================================================================================###

def write_template():
    print '#Input file'
    print 'parentFolder = '
    print '#Library folder to use'
    print 'childFolder = '
    print '#Output folder'
    print 'outputFolder = '
    print '#Base name for files'
    print 'baseName = '
    print '#End name for files'
    print 'endName = '
    print '#Quatfit location'
    print 'quatfit = /opt/apps/Quatfit/quatfit'
    print '#List of atoms to align to'
    print 'parent = '
    print 'child = '
    print 'weights = '
    print '#List of atoms to delete from the input (parent) file'
    print 'deleteParent = '
    print 'deleteChild = '
    return


def parse_input_file(ifile):

  # ## --- Open input settings file and parse into variable --- ###

    f = open(ifile)
    inputlist = f.readlines()
    f.close()

  # ## --- Set up variables for the settings file --- ###

    mapinput = MAPArguments()
    mapinput.parentFolder = ''
    mapinput.childFolder = ''
    mapinput.outputFolder = ''
    mapinput.alignParent = []
    mapinput.alignChild = []
    mapinput.alignWeights = []
    mapinput.parentDelete = []
    mapinput.parentDeleteTemp = []
    mapinput.childDelete = []
    mapinput.childDeleteTemp = []
    mapinput.quatfit = 'quatfit'
    mapinput.baseName = ''
    mapinput.endName = ''

  # ## --- Parse line by line the settings file --- ###

    for var in inputlist:
        line = var.split()
        if line[0] == 'parentFolder' and str(line[-1]) != '=':
            mapinput.parentFolder = str(line[-1])
        elif line[0] == 'childFolder' and str(line[-1]) != '=':
            mapinput.childFolder = str(line[-1])
        elif line[0] == 'outputFolder' and str(line[-1]) != '=':
            mapinput.outputFolder = str(line[-1])
        elif line[0] == 'baseName' and str(line[-1]) != '=':
            mapinput.baseName = str(line[-1])
        elif line[0] == 'endName' and str(line[-1]) != '=':
            mapinput.endName = str(line[-1])
        elif line[0] == 'quatfit' and str(line[-1]) != '=':
            mapinput.quatfit = str(line[-1])
        elif line[0] == 'parent' and str(line[-1]) != '=':
            for i in range(len(line)):  # For multiple entries add them to a list
                if i >= 2:
                    mapinput.alignParent.append(line[i])
        elif line[0] == 'child' and str(line[-1]) != '=':
            for i in range(len(line)):
                if i >= 2:
                    mapinput.alignChild.append(line[i])
        elif line[0] == 'weights' and str(line[-1]) != '=':
            for i in range(len(line)):
                if i >= 2:
                    mapinput.alignWeights.append(line[i])
        elif line[0] == 'deleteParent' and str(line[-1]) != '=':
            for i in range(len(line)):
                if i >= 2:
                    mapinput.parentDeleteTemp.append(line[i])
        elif line[0] == 'deleteChild' and str(line[-1]) != '=':
            for i in range(len(line)):
                if i >= 2:
                    mapinput.childDeleteTemp.append(line[i])
    return mapinput


### --- Expand upon the parentDeleteTemp list if
# there are certain characters to denote multiple atoms.
###

def expand_temp_lists(mapinput):
    for i in range(len(mapinput.parentDeleteTemp)):
        if '-' in str(mapinput.parentDeleteTemp[i]):  # Use the '-' symbol to denote all atoms between the former and latter
            line = str(mapinput.parentDeleteTemp[i]).split('-')
            start = int(line[0])
            for j in xrange(int(line[0]), int(line[1]) + 1):

        # print j

                mapinput.parentDelete.append(j)
        if '-' not in str(mapinput.parentDeleteTemp[i]):  # If not one of the special symbols just add it to the final list
            mapinput.parentDelete.append(int(mapinput.parentDeleteTemp[i]))
    for i in range(len(mapinput.childDeleteTemp)):
        if '-' in str(mapinput.childDeleteTemp[i]):  # Use the '-' symbol to denote all atoms between the former and latter
            line = str(mapinput.childDeleteTemp[i]).split('-')
            start = int(line[0])
            for j in xrange(int(line[0]), int(line[1]) + 1):

        # print j

                mapinput.childDelete.append(j)
        if '-' not in str(mapinput.childDeleteTemp[i]):  # If not one of the special symbols just add it to the final list
            mapinput.childDelete.append(int(mapinput.childDeleteTemp[i]))
    return mapinput


### --- Make the output folder --- ###
# etaatom.make_dir(outputFolder)
# etaatom.make_dir(outputFolder + "/temp")

### --- Make the align files for quatfit. --- ###

def make_align_file(mapinput):
    if len(mapinput.alignWeights) < len(mapinput.alignParent):
        for i in range(len(mapinput.alignParent)):
            mapinput.alignWeights.append(1)
    f = open(mapinput.outputFolder + '/align.txt', 'w+')
    f.write(str(len(mapinput.alignParent)) + '\n')
    for i in range(len(mapinput.alignParent)):
        linetemp = str(mapinput.alignParent[i]) + ' ' \
            + str(mapinput.alignChild[i]) + ' ' \
            + str(mapinput.alignWeights[i]) + '\n'  # ie. 3 4 1 (atoms 3 and 4 will be aligned with a weight of 1)
        f.write(linetemp)
    f.close()
    f = open(mapinput.outputFolder + '/alignparent.txt', 'w+')
    f.write(str(len(mapinput.alignParent)) + '\n')
    for i in range(len(mapinput.alignParent)):
        linetemp = str(mapinput.alignParent[i]) + ' ' \
            + str(mapinput.alignParent[i]) + ' ' \
            + str(mapinput.alignWeights[i]) + '\n'  # ie. 3 4 1 (atoms 3 and 4 will be aligned with a weight of 1)
        f.write(linetemp)
    f.close()
    return


### --- Align the parent files before aligning the substrates to them --- ###
### Purely for aesthetics and does not change the final outcome

def align_parents(mapinput):
    mainparent = ''
    print 'Aligning parent strutures:'
    for aF in os.listdir(mapinput.parentFolder):
        if aF.endswith('.xyz'):
            parentfile = mapinput.parentFolder + '/' + aF
            if mainparent == '':
                mainparent = parentfile
            else:
                print 'Aligning ' + parentfile + ' to ' + mainparent
                os.system(mapinput.quatfit + ' -r' + mainparent + ' -f '
                           + parentfile + ' -p '
                          + mapinput.outputFolder
                          + '/alignparent.txt -o ' + parentfile
                          + '.tmp >/dev/null')
                shutil.move(parentfile + '.tmp', parentfile)
    return


### --- Align and output the parents with the children --- ###

def align_children(mapinput, triangle):

  # ## --- If the parents and children are the same or triangle is enforced only run the triangle of matches --- ###

    if mapinput.parentFolder == mapinput.childFolder or triangle:
        print 'Running only the triangle TRIANGLE!'

    # ## --- Iterating through the parent files in the parent folder and children files in the children folder making lists --- ###

        parentfilelist = []
        childfilelist = []
        for pf in os.listdir(mapinput.parentFolder):
            if pf.endswith('.xyz'):
                parentfilelist.append(mapinput.parentFolder + '/' + pf)
        for cf in os.listdir(mapinput.childFolder):
            if cf.endswith('.xyz'):
                childfilelist.append(cf)

    # ## --- Iterate through the lists just made and align and merge them

        for i in range(len(parentfilelist)):
            parentfile = parentfilelist[i]
            parentfilelol = etaatom.xyz_lol(parentfile)
            for j in range(i, len(childfilelist)):
                childfile = childfilelist[j]
                childfiletemp = mapinput.outputFolder + '/temp/' \
                    + childfile

                print 'Aligning ' + parentfile + ' and ' + childfile

        # ## Running program outside the script. ###

                if mapinput.debug == 1:
                    os.system(mapinput.quatfit + ' -r' + parentfile
                              + ' -f ' + mapinput.childFolder + '/'
                              + childfile + ' -p '
                              + mapinput.outputFolder + '/align.txt -o '
                               + childfiletemp)
                    print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                else:
                    os.system(mapinput.quatfit + ' -r' + parentfile
                              + ' -f ' + mapinput.childFolder + '/'
                              + childfile + ' -p '
                              + mapinput.outputFolder + '/align.txt -o '
                               + childfiletemp + ' >/dev/null')

        # ## --- Open child library file and input each line into a list --- ###

                libfilelol = etaatom.xyz_lol(childfiletemp)

        # ## --- Get new file length --- ###

                newfileatomnum = int(parentfilelol[0]) \
                    + int(libfilelol[0]) - len(mapinput.parentDelete) \
                    - len(mapinput.childDelete)
                if mapinput.baseName == '':
                    if mapinput.endName == '':
                        newfilename = ''.join(parentfilelol[1]) + '-' \
                            + ''.join(libfilelol[1]) + '.xyz'
                    else:
                        newfilename = ''.join(parentfilelol[1]) + '-' \
                            + ''.join(libfilelol[1]) + '-' \
                            + mapinput.endName + '.xyz'
                else:
                    if mapinput.endName == '':
                        newfilename = mapinput.baseName + '-' \
                            + ''.join(parentfilelol[1]) + '-' \
                            + ''.join(libfilelol[1]) + '.xyz'
                    else:
                        newfilename = mapinput.baseName + '-' \
                            + ''.join(parentfilelol[1]) + '-' \
                            + ''.join(libfilelol[1]) + '-' \
                            + mapinput.endName + '.xyz'
                print 'Exporting ' + newfilename

        # ## --- Output file WITHOUT parentDelete atoms and with the childfile atoms --- ###

                f = open(mapinput.outputFolder + '/' + newfilename, 'w+'
                         )
                f.write(str(newfileatomnum) + '\n')
                f.write(str(newfilename)[:-4] + '\n')
                for k in range(0, int(parentfilelol[0])):
                    if k + 1 not in mapinput.parentDelete:
                        linetemp = parentfilelol[k + 2].e + '  ' \
                            + str('{:.6f}'.format(parentfilelol[k
                                  + 2].x)) + '  ' \
                            + str('{:.6f}'.format(parentfilelol[k
                                  + 2].y)) + '  ' \
                            + str('{:.6f}'.format(parentfilelol[k
                                  + 2].z)) + '\n'

            # linetemp = '     '.join(str(e) for e in parentfilelol[k + 2]) + "\n"

                        f.write(linetemp)
                for k in range(0, int(libfilelol[0])):
                    if k + 1 not in mapinput.childDelete:
                        linetemp = libfilelol[k + 2].e + '  ' \
                            + str('{:.6f}'.format(libfilelol[k + 2].x)) \
                            + '  ' + str('{:.6f}'.format(libfilelol[k
                                + 2].y)) + '  ' \
                            + str('{:.6f}'.format(libfilelol[k + 2].z)) \
                            + '\n'

            # linetemp = '     '.join(str(e) for e in libfilelol[k + 2]) + "\n"

                        f.write(linetemp)
                f.close()
    else:

  # ##############################################

        for pf in os.listdir(mapinput.parentFolder):
            if pf.endswith('.xyz'):
                parentfile = mapinput.parentFolder + '/' + pf

        # ## --- Iterating through a folder of files --- ###

                for i in os.listdir(mapinput.childFolder):
                    if i.endswith('.xyz'):
                        childfile = i
                        print 'Aligning ' + pf + ' and ' + i

            # ## Running program outside the script. ###

                        if mapinput.debug >= 1:
                            os.system(mapinput.quatfit + ' -r'
                                    + parentfile + ' -f '
                                    + mapinput.childFolder + '/'
                                    + childfile + ' -p '
                                    + mapinput.outputFolder
                                    + '/align.txt -o '
                                    + mapinput.outputFolder + '/temp/'
                                    + childfile)
                            print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                        else:
                            os.system(mapinput.quatfit + ' -r'
                                    + parentfile + ' -f '
                                    + mapinput.childFolder + '/'
                                    + childfile + ' -p '
                                    + mapinput.outputFolder
                                    + '/align.txt -o '
                                    + mapinput.outputFolder + '/temp/'
                                    + childfile + ' >/dev/null')

        # ## --- Open parent file and input each line into a list --- ###

                parentfilelol = etaatom.xyz_lol(parentfile)

        # ## --- Iterating through a folder of files --- ###

                for i in os.listdir(mapinput.outputFolder + '/temp'):
                    if i.endswith('.xyz'):
                        childfile = mapinput.outputFolder + '/temp/' + i

            # print childfile
            # ## --- Open child library file and input each line into a list --- ###

                        libfilelol = etaatom.xyz_lol(childfile)

            # ## --- Get new file length --- ###

                        newfileatomnum = int(parentfilelol[0]) \
                            + int(libfilelol[0]) \
                            - len(mapinput.parentDelete) \
                            - len(mapinput.childDelete)
                        if mapinput.baseName == '':
                            if mapinput.endName == '':
                                newfilename = ''.join(parentfilelol[1]) \
                                    + '-' + ''.join(libfilelol[1]) \
                                    + '.xyz'
                            else:
                                newfilename = ''.join(parentfilelol[1]) \
                                    + '-' + ''.join(libfilelol[1]) \
                                    + '-' + mapinput.endName + '.xyz'
                        else:
                            if mapinput.endName == '':
                                newfilename = mapinput.baseName + '-' \
                                    + ''.join(parentfilelol[1]) + '-' \
                                    + ''.join(libfilelol[1]) + '.xyz'
                            else:
                                newfilename = mapinput.baseName + '-' \
                                    + ''.join(parentfilelol[1]) + '-' \
                                    + ''.join(libfilelol[1]) + '-' \
                                    + mapinput.endName + '.xyz'
                        print 'Exporting ' + newfilename

            # ## --- Output file WITHOUT parentDelete atoms and with the childfile atoms --- ###

                        f = open(mapinput.outputFolder + '/'
                                 + newfilename, 'w+')
                        f.write(str(newfileatomnum) + '\n')
                        f.write(str(newfilename)[:-4] + '\n')
                        for j in range(0, int(parentfilelol[0])):
                            if j + 1 not in mapinput.parentDelete:
                                linetemp = parentfilelol[j + 2].e \
                                    + '  ' \
                                    + str('{:.6f}'.format(parentfilelol[j
                                        + 2].x)) + '  ' \
                                    + str('{:.6f}'.format(parentfilelol[j
                                        + 2].y)) + '  ' \
                                    + str('{:.6f}'.format(parentfilelol[j
                                        + 2].z)) + '\n'

                # linetemp = '     '.join(str(e) for e in parentfilelol[j + 2]) + "\n"

                                f.write(linetemp)
                        for j in range(0, int(libfilelol[0])):
                            if j + 1 not in mapinput.childDelete:
                                linetemp = libfilelol[j + 2].e + '  ' \
                                    + str('{:.6f}'.format(libfilelol[j
                                        + 2].x)) + '  ' \
                                    + str('{:.6f}'.format(libfilelol[j
                                        + 2].y)) + '  ' \
                                    + str('{:.6f}'.format(libfilelol[j
                                        + 2].z)) + '\n'

                # linetemp = '     '.join(str(e) for e in libfilelol[j + 2]) + "\n"

                                f.write(linetemp)
                        f.close()
    return


