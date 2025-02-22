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
from sys import *
import shutil


###==================================================================================================================###
### --- Atom class --- ###

class Atom(object):

    def __init__(self):

    # Hold the string for the element

        self.e = ''

    # Hold the int for the element

        self.en = int(0)

    # Hold the float for the x coordinate

        self.x = float(0)

    # Hold the float for the y coordinate

        self.y = float(0)

    # Hold the float for the z coordinate

        self.z = float(0)

    # Hold a list of the neighbors

        self.neighbors = []

    # Hold a list of the distances for the neighbors

        self.neighborsdist = []

    # Hold the int for the atom nearest to this one

        self.nearest = int(0)

    # Hold the int for the atom nearest to this one

        self.nearestdist = float(100000000)

    # Hold a string for hybridization for this atom

        self.hybridization = ''

    # Hold an int/float for the charge of this atom

        self.charge = float(0)

    # Hold an int

        self.int = 0

    # Hold a string

        self.string = ''

    # Hold a string for the atomtype of this atom

        self.atomtype = ''


class AtomZMAT(object):

    def __init__(self):

    # Hold the string for the element

        self.e = ''

    # Hold the int for the element

        self.en = int(0)

    # Hold the float for the distance

        self.dist = float(0)

    # Hold the float for the angle

        self.angle = float(0)

    # Hold the float for the dihedral

        self.dihedral = float(0)

    # Hold the float for the distance

        self.atomdist = int(0)

    # Hold the float for the angle

        self.atomangle = int(0)

    # Hold the float for the dihedral

        self.atomdihedral = int(0)


### --- Input Arguments class --- ###

class InputArguments(object):

    def __init__(self):

    # ifile, config, charge, multi, writemod, modred, writeecp, ecplines, writeclose, closelines
    # Hold the string for the input file

        self.ifile = ''

    # Hold the int for the charge

        self.charge = 0

    # Hold the int for the multiplicity

        self.multi = 1

    # Hold the list of config lines

        self.config = []

    # Hold the int/logic for the modredundant

        self.writemod = 0

    # Hold the list of modredundant lines

        self.modred = []

    # Hold the ecp string

        self.ecp = ''

    # Hold the int/logic for the ecp lines

        self.writeecp = 0

    # Hold the list of ecp lines

        self.ecplines = []

    # Hold the int/logic for the close lines

        self.writeclose = 0

    # Hold the list of closing lines

        self.closelines = []

    # Hold the close ecp string

        self.closeecp = ''

    # Hold the int/logic for the close ecp lines

        self.writecloseecp = 0

    # Hold the list of closing ecp lines

        self.closeecplines = []


### --- Reset Input Arguments class --- ###

def reset_input_variables(g0xinput):

  # Hold the string for the input file

    g0xinput.ifile = ''

  # Hold the int for the charge

    g0xinput.charge = 0

  # Hold the int for the multiplicity

    g0xinput.multi = 1

  # Hold the list of config lines

    g0xinput.config = []

  # Hold the int/logic for the modredundant

    g0xinput.writemod = 0

  # Hold the list of modredundant lines

    g0xinput.modred = []

  # Hold the ecp string

    g0xinput.ecp = ''

  # Hold the int/logic for the ecp lines

    g0xinput.writeecp = 0

  # Hold the list of ecp lines

    g0xinput.ecplines = []

  # Hold the int/logic for the close lines

    g0xinput.writeclose = 0

  # Hold the list of closing lines

    g0xinput.closelines = []

  # Hold the close ecp string

    g0xinput.closeecp = ''

  # Hold the int/logic for the close ecp lines

    g0xinput.writecloseecp = 0

  # Hold the list of closing ecp lines

    g0xinput.closeecplines = []
    return g0xinput


# Constants

auToAngstrom = 1.889725988579
hartreeTokcalmol = 627.509469
hartreeToeV = 27.21138386
kcalmolTokjoulemol = 4.184


###==================================================================================================================###

### --- Functions for Submission --- ###
###==================================================================================================================###
### --- Modify the text and return it --- ###

def return_modified_snippet(
    basename,
    program,
    queue,
    nproc,
    memory,
    sfile,
    ):

  # Setup some variables for snippets

    try:
        etadir = os.environ['ETADIR']
    except os.error:
        etadir = '$ETADIR'
    try:
        appdir = os.environ['APPDIR']
    except os.error:
        appdir = '$APPDIR'
    try:
        scriptdir = os.environ['SCRIPTDIR']
    except os.error:
        scriptdir = '$SCRIPTDIR'

  # Start building up the modified text

    printlines = []
    f = open(sfile, 'r')
    for line in f:
        if 'BASENAME!' in line:
            line = line.replace('BASENAME!', str(basename))
        if 'PROGRAM!' in line:
            line = line.replace('PROGRAM!', str(program))
        if 'QUE!' in line:
            line = line.replace('QUE!', str(queue))
        if 'PROCS!' in line:
            line = line.replace('PROCS!', str(nproc))
        if 'MEM!' in line:
            line = line.replace('MEM!', str(memory))
        if 'ETADIR!' in line:
            line = line.replace('ETADIR!', str(etadir))
        if 'SCRIPTDIR!' in line:
            line = line.replace('SCRIPTDIR!', str(scriptdir))
        if 'APPDIR!' in line:
            line = line.replace('APPDIR!', str(appdir))
        printlines.append(line)
    f.close()
    return printlines


def return_modified_string_list(itext, otext, flist):

  # Start building up the modified text

    printlines = []
    for line in flist:
        if itext in line:
            print('Found instance of: ' + itext)
            line = line.replace(itext, otext)
        printlines.append(line)
    return printlines


###==================================================================================================================###

### --- Functions for files and directories --- ###
###==================================================================================================================###
### --- Make directories if they do not exist --- ###

def make_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return


### --- Delete file if it exists --- ###

def del_file(sfile):
    if os.path.isfile(sfile):
        os.remove(sfile)
    return


### --- Copy file if it exists --- ###

def copy_file(ifile, ofile):
    if os.path.isfile(ifile):
        shutil.copyfile(ifile, ofile)
    return


### --- Make directories if they do not exist --- ###

def basename(name, extension):
    basename = name[::-1].replace(extension[::-1], '')[::-1]
    return basename


###==================================================================================================================###

### --- Functions to get and give element numbers and names --- ###
###==================================================================================================================###
# A lowercase list of the elements

elementList = [
    'h',
    'he',
    'li',
    'be',
    'b',
    'c',
    'n',
    'o',
    'f',
    'ne',
    'na',
    'mg',
    'al',
    'si',
    'p',
    's',
    'cl',
    'ar',
    'k',
    'ca',
    'sc',
    'ti',
    'v',
    'cr',
    'mn',
    'fe',
    'co',
    'ni',
    'cu',
    'zn',
    'ga',
    'ge',
    'as',
    'se',
    'br',
    'kr',
    'rb',
    'sr',
    'y',
    'zr',
    'nb',
    'mo',
    'tc',
    'ru',
    'rh',
    'pd',
    'ag',
    'cd',
    'in',
    'sn',
    'sb',
    'te',
    'i',
    'xe',
    'cs',
    'ba',
    'la',
    'ce',
    'pr',
    'nd',
    'pm',
    'sm',
    'eu',
    'gd',
    'tb',
    'dy',
    'ho',
    'er',
    'tm',
    'yb',
    'lu',
    'hf',
    'ta',
    'w',
    're',
    'os',
    'ir',
    'pt',
    'au',
    'hg',
    'tl',
    'pb',
    'bi',
    'po',
    'at',
    'rn',
    'fr',
    'ra',
    'ac',
    'th',
    'pa',
    'u',
    'np',
    'pu',
    'am',
    'cm',
    'bk',
    'cf',
    'es',
    'fm',
    'md',
    'no',
    'lr',
    'rf',
    'db',
    'sg',
    'bh',
    'hs',
    'mt',
    'ds',
    'rg',
    'cn',
    'uut',
    'fl',
    'uup',
    'lv',
    'uus',
    'uuo',
    ]

# A list of the elements with there standard capitalization.

elementNames = [
    'H',
    'He',
    'Li',
    'Be',
    'B',
    'C',
    'N',
    'O',
    'F',
    'Ne',
    'Na',
    'Mg',
    'Al',
    'Si',
    'P',
    'S',
    'Cl',
    'Ar',
    'K',
    'Ca',
    'Sc',
    'Ti',
    'V',
    'Cr',
    'Mn',
    'Fe',
    'Co',
    'Ni',
    'Cu',
    'Zn',
    'Ga',
    'Ge',
    'As',
    'Se',
    'Br',
    'Kr',
    'Rb',
    'Sr',
    'Y',
    'Zr',
    'Nb',
    'Mo',
    'Tc',
    'Ru',
    'Rh',
    'Pd',
    'Ag',
    'Cd',
    'In',
    'Sn',
    'Sb',
    'Te',
    'I',
    'Xe',
    'Cs',
    'Ba',
    'La',
    'Ce',
    'Pr',
    'Nd',
    'Pm',
    'Sm',
    'Eu',
    'Gd',
    'Tb',
    'Dy',
    'Ho',
    'Er',
    'Tm',
    'Yb',
    'Lu',
    'Hf',
    'Ta',
    'W',
    'Re',
    'Os',
    'Ir',
    'Pt',
    'Au',
    'Hg',
    'Tl',
    'Pb',
    'Bi',
    'Po',
    'At',
    'Rn',
    'Fr',
    'Ra',
    'Ac',
    'Th',
    'Pa',
    'U',
    'Np',
    'Pu',
    'Am',
    'Cm',
    'Bk',
    'Cf',
    'Es',
    'Fm',
    'Md',
    'No',
    'Lr',
    'Rf',
    'Db',
    'Sg',
    'Bh',
    'Hs',
    'Mt',
    'Ds',
    'Rg',
    'Cn',
    'Uut',
    'Fl',
    'Uup',
    'Lv',
    'Uus',
    'Uuo',
    ]

# A list of the masses for each element (NOTE: I cannot attest to the scientific accuracy for these masses)

elementMass = [
    1.0079,
    4.0026,
    6.941,
    9.0122,
    10.811,
    12.0107,
    14.0067,
    15.9994,
    18.9984,
    20.1797,
    22.9897,
    24.305,
    26.9815,
    28.0855,
    30.9738,
    32.065,
    35.453,
    39.0983,
    39.948,
    40.078,
    44.9559,
    47.867,
    50.9415,
    51.9961,
    54.938,
    55.845,
    58.6934,
    58.9332,
    63.546,
    65.39,
    69.723,
    72.64,
    74.9216,
    78.96,
    79.904,
    83.8,
    85.4678,
    87.62,
    88.9059,
    91.224,
    92.9064,
    95.94,
    98,
    101.07,
    102.9055,
    106.42,
    107.8682,
    112.411,
    114.818,
    118.71,
    121.76,
    126.9045,
    127.6,
    131.293,
    132.9055,
    137.327,
    138.9055,
    140.116,
    140.9077,
    144.24,
    145,
    150.36,
    151.964,
    157.25,
    158.9253,
    162.5,
    164.9303,
    167.259,
    168.9342,
    173.04,
    174.967,
    178.49,
    180.9479,
    183.84,
    186.207,
    190.23,
    192.217,
    195.078,
    196.9665,
    200.59,
    204.3833,
    207.2,
    208.9804,
    209,
    210,
    222,
    223,
    226,
    227,
    231.0359,
    232.0381,
    237,
    238.0289,
    243,
    244,
    247,
    247,
    251,
    252,
    257,
    258,
    259,
    261,
    262,
    262,
    264,
    266,
    268,
    272,
    277,
    ]

# A list of the masses for each element (NOTE: I cannot attest to the scientific accuracy for these masses)

elementIsotope = [
    [1, 2, 3],
    [4, 3],
    [7, 6],
    [9],
    [11, 10],
    [12, 13],
    [14, 15],
    [16, 18, 17],
    [19],
    [20, 22, 21],
    [23],
    [24, 26, 25],
    [27],
    [28, 29, 30],
    [32],
    [32, 34, 33, 36],
    [35, 37],
    [40, 36, 38],
    [39, 41, 40],
    [
        40,
        44,
        42,
        48,
        43,
        46,
        ],
    [45],
    [48, 46, 47, 49, 50],
    [51, 50],
    [52, 53, 50, 54],
    [55],
    [56, 54, 57, 58],
    [59],
    [58, 60, 62, 61, 64],
    [63, 65],
    [64, 66, 68, 67, 70],
    [69, 71],
    [74, 72, 70, 73, 76],
    [75],
    [
        80,
        78,
        76,
        82,
        77,
        74,
        ],
    [79, 81],
    [
        84,
        86,
        82,
        83,
        80,
        78,
        ],
    [85, 87],
    [88, 86, 87, 84],
    [89],
    [90, 94, 92, 91, 96],
    [93],
    [
        98,
        96,
        95,
        92,
        97,
        94,
        ],
    [99],
    [
        102,
        104,
        101,
        99,
        100,
        96,
        98,
        ],
    [103],
    [
        106,
        108,
        105,
        110,
        104,
        102,
        ],
    [107, 109],
    [
        114,
        112,
        111,
        110,
        113,
        116,
        106,
        108,
        ],
    [115, 113],
    [
        120,
        118,
        116,
        119,
        117,
        124,
        122,
        112,
        114,
        116,
        ],
    [121, 123],
    [
        130,
        128,
        126,
        125,
        124,
        122,
        123,
        120,
        ],
    [127],
    [
        132,
        129,
        131,
        134,
        136,
        130,
        128,
        126,
        124,
        ],
    [133],
    [
        138,
        137,
        136,
        135,
        134,
        130,
        132,
        ],
    [139, 138],
    [140, 142, 138, 136],
    [141],
    [
        142,
        144,
        146,
        143,
        145,
        148,
        150,
        ],
    [147],
    [
        152,
        154,
        147,
        149,
        148,
        150,
        144,
        ],
    [153, 151],
    [
        158,
        160,
        156,
        157,
        155,
        154,
        152,
        ],
    [159],
    [
        164,
        162,
        163,
        161,
        160,
        158,
        156,
        ],
    [165],
    [
        166,
        168,
        167,
        170,
        164,
        162,
        ],
    [169],
    [
        174,
        172,
        173,
        171,
        176,
        170,
        168,
        ],
    [175, 176],
    [
        180,
        178,
        177,
        179,
        176,
        174,
        ],
    [181],
    [184, 186, 182, 183, 180],
    [187, 185],
    [
        192,
        190,
        189,
        188,
        187,
        186,
        184,
        ],
    [193, 191],
    [
        195,
        194,
        196,
        198,
        192,
        190,
        ],
    [197],
    [
        202,
        200,
        199,
        201,
        198,
        204,
        196,
        ],
    [205, 203],
    [208, 206, 207, 204],
    [209],
    [210],
    [210],
    [222],
    [223],
    [226],
    [227],
    [232],
    [231],
    [238],
    [237],
    [244],
    [243],
    [247],
    [247],
    [251],
    [252],
    [257],
    [258],
    [259],
    [260],
    ]

# A list of the atoms that are 3rd row or below (aka large)

elementLarge = [
    'na',
    'mg',
    'al',
    'si',
    'p',
    's',
    'cl',
    'ar',
    'k',
    'ca',
    'sc',
    'ti',
    'v',
    'cr',
    'mn',
    'fe',
    'co',
    'ni',
    'cu',
    'zn',
    'ga',
    'ge',
    'as',
    'se',
    'br',
    'kr',
    'rb',
    'sr',
    'y',
    'zr',
    'nb',
    'mo',
    'tc',
    'ru',
    'rh',
    'pd',
    'ag',
    'cd',
    'in',
    'sn',
    'sb',
    'te',
    'i',
    'xe',
    'cs',
    'ba',
    'la',
    'ce',
    'pr',
    'nd',
    'pm',
    'sm',
    'eu',
    'gd',
    'tb',
    'dy',
    'ho',
    'er',
    'tm',
    'yb',
    'lu',
    'hf',
    'ta',
    'w',
    're',
    'os',
    'ir',
    'pt',
    'au',
    'hg',
    'tl',
    'pb',
    'bi',
    'po',
    'at',
    'rn',
    'fr',
    'ra',
    'ac',
    'th',
    'pa',
    'u',
    'np',
    'pu',
    'am',
    'cm',
    'bk',
    'cf',
    'es',
    'fm',
    'md',
    'no',
    'lr',
    'rf',
    'db',
    'sg',
    'bh',
    'hs',
    'mt',
    'ds',
    'rg',
    'cn',
    'uut',
    'fl',
    'uup',
    'lv',
    'uus',
    'uuo',
    ]


# Check the elementList for the lowercase string and return the index for that number + 1 (the atomic number)

def get_element_num(at1):
    element = elementList.index(at1.lower())
    return element + 1


# Check the elementName list for index (at1 - 1) and return the standard capitalization

def get_element_name(at1):
    element = elementNames[at1 - 1]
    return element


# Check the elementMass list for index (at1 - 1) and return the mass

def get_element_mass(at1):
    atomicmass = elementMass[at1 - 1]
    return float(atomicmass)


###==================================================================================================================###

### --- Lists and function to get ECP data and molecular charge --- ###
###==================================================================================================================###
# A list of the element in the LANL2DZ ECP basis set

elementLANL2DZ_K = [
    'k',
    'ca',
    'sc',
    'ti',
    'v',
    'cr',
    'mn',
    'fe',
    'co',
    'ni',
    'cu',
    'zn',
    'ga',
    'ge',
    'as',
    'se',
    'br',
    'kr',
    'rb',
    'sr',
    'y',
    'zr',
    'nb',
    'mo',
    'tc',
    'ru',
    'rh',
    'pd',
    'ag',
    'cd',
    'in',
    'sn',
    'sb',
    'te',
    'i',
    'xe',
    'cs',
    'ba',
    'la',
    'hf',
    'ta',
    'w',
    're',
    'os',
    'ir',
    'pt',
    'au',
    'hg',
    'tl',
    'pb',
    'bi',
    'u',
    'np',
    'pu',
    ]

# A list of the element in the LANL08 ECP basis set

elementLANL08 = [
    'sc',
    'ti',
    'v',
    'cr',
    'mn',
    'fe',
    'co',
    'ni',
    'cu',
    'zn',
    'rb',
    'sr',
    'y',
    'zr',
    'nb',
    'mo',
    'tc',
    'ru',
    'rh',
    'pd',
    'ag',
    'cd',
    'in',
    'sn',
    'sb',
    'te',
    'i',
    'xe',
    'cs',
    'ba',
    'la',
    'hf',
    'ta',
    'w',
    're',
    'os',
    'ir',
    'pt',
    'au',
    'hg',
    'tl',
    'pb',
    'bi',
    ]

# A custom charge list for automatic charge calculation for ionic systems. (Update as needed)

ionicReg = [
    1,
    0,
    1,
    2,
    3,
    -4,
    -3,
    -2,
    -1,
    0,
    1,
    2,
    3,
    'si',
    3,
    's',
    -1,
    'ar',
    1,
    2,
    'sc',
    'ti',
    'v',
    'cr',
    'mn',
    'fe',
    'co',
    '0',
    '1',
    '0',
    3,
    'ge',
    'as',
    'se',
    -1,
    'kr',
    'rb',
    'sr',
    'y',
    'zr',
    5,
    'mo',
    'tc',
    'ru',
    'rh',
    0,
    1,
    'cd',
    3,
    'sn',
    'sb',
    'te',
    'i',
    'xe',
    'cs',
    'ba',
    'la',
    'ce',
    'pr',
    'nd',
    'pm',
    'sm',
    'eu',
    'gd',
    'tb',
    'dy',
    'ho',
    'er',
    'tm',
    'yb',
    'lu',
    3,
    'ta',
    'w',
    're',
    'os',
    'ir',
    0,
    1,
    'hg',
    'tl',
    'pb',
    'bi',
    'po',
    'at',
    'rn',
    'fr',
    'ra',
    'ac',
    'th',
    'pa',
    'u',
    'np',
    'pu',
    'am',
    'cm',
    'bk',
    'cf',
    'es',
    'fm',
    'md',
    'no',
    'lr',
    'rf',
    'db',
    'sg',
    'bh',
    'hs',
    'mt',
    'ds',
    'rg',
    'cn',
    'uut',
    'fl',
    'uup',
    'lv',
    'uus',
    'uuo',
    ]


# Check the chargeList for the element at1 and return the charge associated with it.

def get_element_charge(at1, chargelist):
    element = chargelist[at1]
    return element


### --- Calculate charge from file --- ###

def ionic(ionicswitch, ioniccm, ifilelol):

  # Start with charge zero

    charge = 0

  # Logic to use the ionicReg charge list

    if ioniccm.lower() == 'reg':

    # Iterate through the list of atoms (skipping the number of atoms and the second informational line) in the ifilelol

        for i in range(2, len(ifilelol)):

      # Sum the charge returned from the get_element_charge function, but first cleaning the input with get_element_num

            charge += \
                int(get_element_charge(get_element_num(ifilelol[i].e)
                    - 1, ionicReg))
    return charge


### --- General structure functions --- ###
###==================================================================================================================###
### --- Parse the input XYZ file and convert it into a Atom Class list of lists --- ###

def xyz_lol(ifile):

  # ## --- Open parent file, read it in and close it --- ###

    f = open(ifile, 'r')
    ifilelist = f.readlines()
    f.close()

  # ## --- Create some variables needed for parsing the input file --- ###
  # ifilelength is used to gauge what line we are on and also allows us to determine the number of lines in the file
  # ifilelol is a list of lists (lol) that will hold all the data of the input file
  # alist is a temp list that will be appended to the ifilelol
  # ###############################
  # Find out the final length of the ifilelol

    ifilelength = len(ifilelist)

  # Fill the ifilelol with 0's so that there is a place to put the atomic data

    ifilelol = [0] * ifilelength

  # ### --- Input/parse the input file into a list of lists called ifilelol --- ###

    for i in range(ifilelength):

    # If i == 0 then put the number of atoms in as an int

        if i == 0:
            ifilelol[i] = int(ifilelist[i].rstrip())
        elif i == 1:

    # if i == 1 then put some information in as a string (here it is the second line of the XYZ file)

            ifilelol[i] = str(ifilelist[i].rstrip())
        elif i >= 2:

    # if i >= 2 then make the instance of ifilelol an Atom class then fill it with element data, and coordinate data

            line = ifilelist[i].rstrip().split()
            ifilelol[i] = Atom()
            ifilelol[i].e = line[0]
            ifilelol[i].en = get_element_num(line[0])
            ifilelol[i].x = float(line[1])
            ifilelol[i].y = float(line[2])
            ifilelol[i].z = float(line[3])

  # If the XYZ file is longer or shorter than the list of data should be, warn the user.

    if int(ifilelol[0]) + 2 != ifilelength:
        print('Your file is the wrong length!')
        print('Make sure you have the right number of atoms.')
        exit(0)
    return ifilelol


### --- Output the structural data in an XYZ format

def output_xyz(ofile, ifilelol):

  # Open the file for writing

    f = open(ofile, 'w+')

  # iterate through the whole ifilelol and write each line in XYZ format to the file

    for i in range(len(ifilelol)):

    # if the instance is in Atom form print the 'element X-Coord Y-Coord Z-Coord'

        if isinstance(ifilelol[i], Atom):
            line = ifilelol[i].e + '  ' \
                + str('{:.6f}'.format(ifilelol[i].x)) + '  ' \
                + str('{:.6f}'.format(ifilelol[i].y)) + '  ' \
                + str('{:.6f}'.format(ifilelol[i].z))
        else:

    # else print the string (number of atoms or comment line)

            line = str(ifilelol[i])
        f.write(line + '\n')
    f.close()
    return


### --- Output the trajectory data in an XYZ format and append it to the file

def output_xyz_traj(ofile, ifilelol):

  # Open the file for writing

    f = open(ofile, 'a')

  # iterate through the whole ifilelol and write each line in XYZ format to the file

    for i in range(len(ifilelol)):

    # if the instance is in Atom form print the 'element X-Coord Y-Coord Z-Coord'

        if isinstance(ifilelol[i], Atom):
            line = ifilelol[i].e + '  ' \
                + str('{:.6f}'.format(ifilelol[i].x)) + '  ' \
                + str('{:.6f}'.format(ifilelol[i].y)) + '  ' \
                + str('{:.6f}'.format(ifilelol[i].z))
        else:

    # else print the string (number of atoms or comment line)

            line = str(ifilelol[i])
        f.write(line + '\n')
    f.close()
    return


### --- Renumber atoms to put the atoms of interest at the front in the correct order --- ###

def renumber_basic(ifilelol, badlist):

  # Create the real ifilelol and add the first two lines

    ifilelol_renumbered = ['0'] * len(ifilelol)
    ifilelol_renumbered[0] = ifilelol[0]
    ifilelol_renumbered[1] = ifilelol[1]

  # Iterate through the bond_angle_dihedral list and put them at the top of the list

    for (i, atom) in enumerate(badlist):
        ifilelol_renumbered[i + 2] = ifilelol[atom + 1]

  # Iterate throught the ifilelol and add the rest of the file to the renumbered

    usedbad = 0
    for i in range(2, len(ifilelol)):
        if i - 1 not in badlist:
            ifilelol_renumbered[i + len(badlist) - usedbad] = \
                ifilelol[i]
        else:
            usedbad += 1
    return ifilelol_renumbered


### --- Import a large trajectory and export only the Xth geometry --- ###

def get_geom_n(trajectory, geometry):
    numatoms = int(trajectory[0])
    geometrylength = numatoms + 2
    startline = geometrylength * (geometry - 1)
    endline = startline + geometrylength
    ifilelol = [0] * geometrylength
    count = 0
    for i in range(startline, endline):
        if isinstance(trajectory[i], Atom):
            ifilelol[count] = Atom()
            ifilelol[count] = trajectory[i]
        else:
            ifilelol[count] = trajectory[i]
        count += 1
    return ifilelol


### --- Function to easily print out a geom using the atom class for debuging --- ###

def print_geom(geom):
    for i in range(len(geom)):
        if isinstance(geom[i], Atom):
            line = geom[i].e + '  ' + str('{:.6f}'.format(geom[i].x)) \
                + '  ' + str('{:.6f}'.format(geom[i].y)) + '  ' \
                + str('{:.6f}'.format(geom[i].z))
        else:

    # else print the string (number of atoms or comment line)

            line = str(geom[i])
        print(line)
    return


### --- Snippet functions --- ###
###==================================================================================================================###
### --- Open and parse the snippet file --- ###

def parse_snippet(snippet):

  # open the file for reading

    f = open(snippet)
    sniplist = f.readlines()
    f.close()

  # set up some blank/zeroed variables

    configswitch = 0
    config = []
    ecpswitch = 0
    ecp = ''
    closeswitch = 0
    closelines = []
    writeecp = 0
    writeclose = 0

  # iterate through the length of the snippet file line by line

    for (i, var) in enumerate(sniplist):

    # create a split line

        line = var.strip()

    # if the end of a section is set stop recording

        if line == '#END':
            ecpswitch = 0
            closeswitch = 0
            configswitch = 0
        elif configswitch == 1:

    # For the config section start collecting the lines

            config.append(line)
        elif ecpswitch == 1:

    # For the ecp section start collecting the lines

            ecp = line
        elif closeswitch == 1:

    # For the close section start collecting the lines

            closelines.append(line)
        elif line == '#CONFIG':

    # use the # lines to turn on collection of sections

            configswitch = 1
        elif line == '#ECP':
            writeecp = 1
            ecpswitch = 1
        elif line == '#CLOSE':
            writeclose = 1
            closeswitch = 1
    return (config, ecp, closelines, writeecp, writeclose)


### --- Interactive mode for CONFIG --- ###

def interactive_config(SOMETHING):

  # set up some blank/zeroed variables

    configswitch = 0
    config = []
    ecpswitch = 0
    ecp = ''
    closeswitch = 0
    closelines = []
    writeecp = 0
    writeclose = 0

  # Figure out the method

    var = eval(input('Method:'))
    if var == 1:
        config.append('#B3LYP/')
    elif var == 2:
        config.append('#M06/')
    elif var == 3:
        config.append('#M06-2X/')
    elif var == 4:
        config.append('#M06-2X/')

  # Figure out the basis set

    var = eval(input('Will there be ECPs? (y/n)'))
    if var.lower() == 'y':
        var = eval(input('ECPs Type:'))
        if var == 1:
            config.append('gen ')
            ecp.append('LANL2DZ')
        elif var == 2:
            config.append('gen ')
            ecp.append('SOMETHINGELSE')
    var = eval(input('Basis set:'))
    if var == 1:
        config.append('6-31G* ')
    elif var == 2:
        config.append('6-311++G** ')
    elif var == 2:
        config.append('6-311++G** ')

  # var = input("Job Type:")
  # var = input("Solvation:")
  # var = input("Extras:")

    if SOMETHING == 1:

    # if the end of a section is set stop recording

        if line == '#END':
            ecpswitch = 0
            closeswitch = 0
            configswitch = 0
        elif configswitch == 1:

    # For the config section start collecting the lines

            config.append(line)
        elif ecpswitch == 1:

    # For the ecp section start collecting the lines

            ecp = line
        elif closeswitch == 1:

    # For the close section start collecting the lines

            closelines.append(line)
        elif line == '#CONFIG':

    # use the # lines to turn on collection of sections

            configswitch = 1
        elif line == '#ECP':
            writeecp = 1
            ecpswitch = 1
        elif line == '#CLOSE':
            writeclose = 1
            closeswitch = 1
    return (config, ecp, closelines, writeecp, writeclose)


### --- Gaussian Specific functions --- ###
###==================================================================================================================###
### --- Open and parse the snippet file --- ###

def parse_snippet_g0x(snippet, g0xinput):

  # open the file for reading

    f = open(snippet)
    sniplist = f.readlines()
    f.close()

  # set up some blank/zeroed variables

    configswitch = 0
    ecpswitch = 0
    closeswitch = 0
    closeecpswitch = 0

  # g0xinput.config = []
  # g0xinput.ecp = ""
  # g0xinput.closeecp = ""
  # g0xinput.closelines = []
  # g0xinput.writeecp = 0
  # g0xinput.writeclose = 0
  # g0xinput.writecloseecp = 0
  # Iterate through the length of the snippet file line by line

    for (i, var) in enumerate(sniplist):

    # create a split line

        line = var.strip()

    # if the end of a section is set stop recording

        if line == '#END':
            ecpswitch = 0
            closeswitch = 0
            configswitch = 0
            closeecpswitch = 0
        elif configswitch == 1:

    # For the config section start collecting the lines

            g0xinput.config.append(line)
        elif ecpswitch == 1:

    # For the ecp section start collecting the lines

            g0xinput.ecp = line
        elif closeswitch == 1:

    # For the close section start collecting the lines

            g0xinput.closelines.append(line)
        elif closeecpswitch == 1:

    # For the close ecp section start collecting the lines

            g0xinput.closeecp = line
        elif line == '#CONFIG':

    # use the # lines to turn on collection of sections

            configswitch = 1
        elif line == '#ECP':
            g0xinput.writeecp = 1
            ecpswitch = 1
        elif line == '#ECPCLOSE':
            g0xinput.writecloseecp = 1
            closeecpswitch = 1
        elif line == '#CLOSE':
            g0xinput.writeclose = 1
            closeswitch = 1
    return g0xinput


### --- Parses the input G0X file --- ###

def parse_input_g0x(ifile):
    config = []
    comment = ''
    charge = 0
    multi = 1
    initialgeom = []
    endlines = []
    inputnumatoms = 0
    emptyline = 0
    endline = 'no'
    the_file = open(ifile, 'r')
    for line in the_file:
        if line.strip() == '' and endline == 'no':
            emptyline += 1
        if line.strip() != '' and emptyline == 0:
            config.append(line.strip())
        if line.strip() != '' and emptyline == 1:
            comment = line.strip()
        if line.strip() != '' and emptyline == 3:
            initialgeom.append(line.strip())
            inputnumatoms += 1
        if line.strip() != '' and emptyline == 2:
            charge = int(line.split()[0])
            multi = int(line.split()[1])
            emptyline += 1
        if line.strip() != '' and emptyline == 4:
            endline = 'yes'
            emptyline += 1
        if emptyline >= 5:
            endlines.append(line.strip())
    ifilelol = [0] * (inputnumatoms + 2)
    ifilelol[0] = inputnumatoms
    ifilelol[1] = comment
    for (i, var) in enumerate(initialgeom):
        line = var.split()
        ifilelol[i + 2] = Atom()
        ifilelol[i + 2].e = line[0]
        ifilelol[i + 2].en = get_element_num(line[0])
        ifilelol[i + 2].x = float(line[1])
        ifilelol[i + 2].y = float(line[2])
        ifilelol[i + 2].z = float(line[3])
    return (config, endlines, charge, multi, ifilelol)


### --- Parses the geometry and energy created by G0X and gives back the last geometry in a list of lists --- ###

def parse_output_g0x(ifile):
    alist = []
    geom = 0
    numatoms = 0
    charge = 0
    multi = 1
    finalenergy = 0.000000
    the_file = open(ifile, 'r')
    for (idx, line) in enumerate(the_file):
        if line.strip() != '':
            if line.split()[0] == 'SCF':
                finalenergy = float(line.split()[4])
            if line.split()[0] == 'NAtoms=':
                numatoms = int(line.split()[1])
            if '-----------------------------------------------' \
                in line and geom >= 6:
                geom = 0
            if geom >= 1:
                geom += 1
            if line.strip() == 'Standard orientation:':
                geom = 1
            if geom >= 5:
                alist.append(line)
            elif 'Multiplicity' in line:

      # Get charge and multiplicity

                chargemulti = line.split()
                charge = chargemulti[2]
                multi = chargemulti[-1]
    the_file.close()
    alistdedup = [x for x in alist if x
                  != ' ---------------------------------------------------------------------\n'
                  ]
    finalgeom = [0] * (numatoms + 2)
    finalgeom[0] = str(numatoms)
    finalgeom[1] = 'Energy: ' + str(finalenergy)

  # print alistdedup

    for i in range(len(alistdedup) - numatoms, len(alistdedup)):
        line = alistdedup[i].split()
        j = i - len(alistdedup) + numatoms + 2
        finalgeom[j] = Atom()
        finalgeom[j].e = get_element_name(int(line[1]))
        finalgeom[j].en = int(line[1])
        finalgeom[j].x = float(line[3])
        finalgeom[j].y = float(line[4])
        finalgeom[j].z = float(line[5])
    return (finalgeom, charge, multi)


### --- Parses the geometry and energy created by G0X and gives back the trajectory in a list of lists --- ###

def parse_output_g0x_traj(ifile):
    alist = []
    geom = 0
    numatoms = 0
    allenergy = []
    trajectorygeom = []
    the_file = open(ifile, 'r')
    for (idx, line) in enumerate(the_file):
        if line.strip() != '':
            if line.split()[0] == 'SCF':
                allenergy.append(float(line.split()[4]))
            if line.split()[0] == 'NAtoms=':
                numatoms = int(line.split()[1])
            if '-----------------------------------------------' \
                in line and geom >= 6:
                geom = 0
            if geom >= 1:
                geom += 1
            if line.strip() == 'Standard orientation:':
                geom = 1
            if geom >= 5:
                alist.append(line)
            elif 'Multiplicity' in line:

      # Get charge and multiplicity

                chargemulti = line.split()
                charge = chargemulti[2]
                multi = chargemulti[-1]
    the_file.close()
    alistdedup = [x for x in alist
                  if '------------------------------------------------------'
                   not in x]
    trajectorygeom = [0] * (len(alistdedup) / numatoms * (numatoms + 2))
    skip = 0
    for i in range(len(alistdedup)):
        if i % numatoms == 0 and i / numatoms != len(allenergy):
            trajectorygeom[i + skip] = int(numatoms)
            trajectorygeom[i + skip + 1] = 'Energy: ' + str(allenergy[i
                    / numatoms])
            skip += 2
        elif i % numatoms == 0 and i / numatoms >= len(allenergy):
            trajectorygeom[i + skip] = int(numatoms)
            trajectorygeom[i + skip + 1] = 'Energy: ' \
                + str(allenergy[-1])
            skip += 2
        line = alistdedup[i].split()
        trajectorygeom[i + skip] = Atom()
        trajectorygeom[i + skip].e = get_element_name(int(line[1]))
        trajectorygeom[i + skip].en = int(line[0])
        trajectorygeom[i + skip].x = float(line[3])
        trajectorygeom[i + skip].y = float(line[4])
        trajectorygeom[i + skip].z = float(line[5])
    return (trajectorygeom, charge, multi)


### --- Logic for the ECP data for G0X --- ###

def g0x_ecp(ecp, config, ifilelol):
    ecplines = []
    ecplist = []
    genlist = []
    if 'LANL2DZ' in ecp:
        for i in range(2, len(ifilelol)):
            if ifilelol[i].e.lower() in elementLANL2DZ_K \
                and ifilelol[i].e not in ecplist:
                ecplist.append(ifilelol[i].e)
            if ifilelol[i].e.lower() not in elementLANL2DZ_K \
                and ifilelol[i].e not in genlist:
                genlist.append(ifilelol[i].e)
        ecplines.append(' '.join(str(e) for e in genlist) + ' 0')
        ecplines.append(ecp.split('/')[0])
        ecplines.append('****')
        if len(ecplist) > 0:
            ecplines.append(' '.join(str(e) for e in ecplist) + ' 0')
            ecplines.append('LANL2DZ')
            ecplines.append('****')
            ecplines.append('')
            ecplines.append(' '.join(str(e) for e in ecplist) + ' 0')
            ecplines.append('LANL2DZ')
        elif len(ecplist) == 0:
            for i in range(len(config)):
                config[i] = config[i].replace('pseudo=read ', '')
    elif 'LANL08' in ecp:
        for i in range(2, len(ifilelol)):
            if ifilelol[i].e.lower() in elementLANL08 and ifilelol[i].e \
                not in ecplist:
                ecplist.append(ifilelol[i].e)
            if ifilelol[i].e.lower() not in elementLANL08 \
                and ifilelol[i].e not in genlist:
                genlist.append(ifilelol[i].e)
        ecplines.append(' '.join(str(e) for e in genlist) + ' 0')
        ecplines.append(ecp.split('/')[0])
        ecplines.append('****')
        if len(ecplist) > 0:

      # ecplines.append(' '.join(str(e) for e in ecplist) + " 0")

            ecplines.append('@/share/apps/g09/basis/lanl08-bs.gbs')

      # ecplines.append("****")

            ecplines.append('')

      # ecplines.append(' '.join(str(e) for e in ecplist) + " 0")

            ecplines.append('@/share/apps/g09/basis/lanl08-ecp.gbs')
        elif len(ecplist) == 0:
            for i in range(len(config)):
                config[i] = config[i].replace('pseudo=read ', '')
    else:
        for i in range(2, len(ifilelol)):
            if ifilelol[i].e.lower() in elementLANL2DZ_K \
                and ifilelol[i].e not in ecplist:
                ecplist.append(ifilelol[i].e)
            if ifilelol[i].e.lower() not in elementLANL2DZ_K \
                and ifilelol[i].e not in genlist:
                genlist.append(ifilelol[i].e)
        ecplines.append(' '.join(str(e) for e in genlist) + ' 0')
        ecplines.append(ecp.split('/')[0])
        ecplines.append('****')
        if len(ecplist) > 0:
            ecplines.append(' ')
            ecplines.append(' '.join(str(e) for e in ecplist) + ' 0')
            ecplines.append(ecp.split('/')[1])
            ecplines.append('****')
            ecplines.append('')
            ecplines.append(' '.join(str(e) for e in ecplist) + ' 0')
            ecplines.append(ecp.split('/')[1])
        elif len(ecplist) == 0:
            for i in range(len(config)):
                config[i] = config[i].replace('pseudo=read ', '')
    return (config, ecplines)


### --- Out put the G0X input files --- ###

def output_g0x(ifile, g0xinput, ifilelol):

  # config, charge, multi, writemod, modred, writeecp, ecplines, writeclose, closelines

    ofiletmp = ifile.split('.')
    ofile = ''
    for i in range(len(ofiletmp)):
        if i != len(ofiletmp) - 1:
            ofile += ofiletmp[i] + '.'
        elif i == len(ofiletmp) - 1:
            ofile += 'com'

  # print ofile

    f = open('NEWJOBS.G0X/' + ofile, 'w+')
    for (i, var) in enumerate(g0xinput.config):
        f.write(var + '\n')
    f.write('\n')
    f.write(ifile + '\n')
    f.write('\n')
    f.write(str(g0xinput.charge) + '  ' + str(g0xinput.multi) + '\n')
    for i in range(2, len(ifilelol)):
        line = ifilelol[i].e + '  ' \
            + str('{:.6f}'.format(ifilelol[i].x)) + '  ' \
            + str('{:.6f}'.format(ifilelol[i].y)) + '  ' \
            + str('{:.6f}'.format(ifilelol[i].z))
        f.write(line + '\n')
    f.write('\n')
    if g0xinput.writemod == 1:
        if len(g0xinput.modred) > 0 and g0xinput.modred[0] != '':
            for (i, var) in enumerate(g0xinput.modred):
                f.write(var + '\n')
        else:
            f.write('MODRED!\n')
        f.write('\n')
    if g0xinput.writeecp == 1:
        for (i, var) in enumerate(g0xinput.ecplines):
            f.write(var + '\n')
        f.write('\n')
    if g0xinput.writeclose == 1:
        for (i, var) in enumerate(g0xinput.closelines):
            f.write(var + '\n')
        f.write('\n')
    if g0xinput.writecloseecp == 1:
        for (i, var) in enumerate(g0xinput.closeecplines):
            f.write(var + '\n')
        f.write('\n')
    f.close()
    return


# Grab the vibrations from the G0X log file

def g0x_vibrations(ifile):
    freqswitch = False
    freqlist = []
    the_file = open(ifile, 'r')
    for (idx, line) in enumerate(the_file):
        if line.strip() != '':
            if line.split()[0] == 'NAtoms=':
                numatoms = int(line.split()[1])
            if 'Low frequencies' in line:
                freqswitch = True
            elif line.split()[0] == 'Frequencies' and freqswitch:

        # Store the frequencies

                freqlist.append(float(line.split()[2]))
                if len(line.split()) >= 4:
                    freqlist.append(float(line.split()[3]))
                if len(line.split()) >= 5:
                    freqlist.append(float(line.split()[4]))
    return freqlist


# Parse the Gaussian output file and print data to the screen

def translate_g0x_output(logfile, trans, printtype):
    initlines = \
        ['\033[1;34mPython Gaussian 0X Output File Parsing Script:   nGT.py\033[1;0m'
         ,
         '\033[1;34mScripted by Thomas J. L. Mustard, O. Maduka Ogba, Paul Ha-Yeon Cheong\033[1;0m'
         ,
         '------------------------------------------------------------------------------'
         ]
    if trans == 'translate':
        for line in initlines:
            print(line)

  # Global Variables

    config = ''
    configswitch = 0
    printlines = []
    programversion = ''
    chk = ''
    nproc = ''
    mem = ''
    configlines = []
    configlinesall = ''
    modswitch = 0
    modlines = []
    pointstoic = ''
    point = ''
    stoic = ''
    chargemulti = ''
    charge = ''
    multi = ''
    basisfunc = ''
    scfenergy = ''
    predechange = ''
    optswitch = 0
    optcomplete = 0
    optstep = 0
    if printtype == 'long' or printtype == 'extra':
        maxforce = []
        rmsdisp = []
        rmsforce = []
        maxdisp = []
    else:
        maxforce = ''
        rmsforce = ''
        maxdisp = ''
        rmsdisp = ''
    freqswitch = 0
    freqtemp = ''
    if trans == 'list':
        frequencies = [
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            ]
        negfreq = 0
        zeropointe = '0.000000'
        therme = '0.000000'
        enthalpye = '0.000000'
        gibbse = '0.000000'
        zeropointesum = '0.000000'
        thermesum = '0.000000'
        enthalpyesum = '0.000000'
        gibbsesum = '0.000000'
        stotal = '0.000000'
        selec = '0.000000'
        srot = '0.000000'
        strans = '0.000000'
        svib = '0.000000'
    else:
        frequencies = [
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            ]
        negfreq = 0
        zeropointe = ''
        therme = ''
        enthalpye = ''
        gibbse = ''
        zeropointesum = ''
        thermesum = ''
        enthalpyesum = ''
        gibbsesum = ''
        stotal = ''
        selec = ''
        srot = ''
        strans = ''
        svib = ''
    jobdone = 0
    memswitch = 0
    f = open(logfile)
    filelist = f.readlines()
    f.close()
    for (i, var) in enumerate(filelist):
        if i == 0 and 'Entering Gaussian System' not in var:
            print('Error! Not a Gaussian log file.')
            exit(0)
        elif 'Cite this work as:' in var:

    # Find the Gaussian Version

            programversion = filelist[i + 1].rstrip().rstrip(',')
        elif '%chk=' in var:

    # Find the chk file

            chk = var.split('=')[1].rstrip()
        elif '%nproc=' in var:

    # Find the nproc used

            nproc = var.split('=')[1].rstrip()
        elif '%mem' in var:

    # Find the mem used

            mem = var.split('=')[1].rstrip()
            memswitch = 1
        elif memswitch == 1 and '------------------' in var:

    # Find the job configuration used

            for j in range(1, 5):
                if '-------' not in filelist[i + j]:
                    configlines.append(filelist[i + j].strip())
                    configlinesall += filelist[i + j].strip()
                else:
                    memswitch = 0
                    break

      # print configlinesall

            if ' opt=' in configlinesall:
                optswitch = 1
            if 'modredundant' in configlinesall:
                modswitch = 1
            if 'freq' in configlinesall:
                freqswitch = 1
        elif modswitch == 1:

    # Check to see if modredundant is in the config

            if 'The following ModRedundant input section has been read:' \
                in var:
                for j in range(1, 20):
                    if len(filelist[i + j]) > 2:
                        modlines.append(filelist[i + j].strip())
                    else:
                        break

    # Get pointgroup and stoic

        if 'Framework group' in var:
            pointstoic = var.split()[-1]
            point = pointstoic.split('[')[0]
        elif 'Stoichiometry' in var:
            stoic = var.split()[-1]
        elif 'Multiplicity' in var:

    # Get charge and multiplicity

            chargemulti = var.split()
            charge = chargemulti[2]
            multi = chargemulti[-1]
        elif 'basis functions,' in var:

    # Find basis functions
    # elif 'Two-electron integral symmetry is turned on.' in var:
    #  basisfunc = filelist[i+1].strip().split()[0]

            basisfunc = var.strip().split()[0]
        elif 'SCF Done:' in var:

    # Get SCF Energy and predicted change

            scfenergy = var.strip().split()[4]
        elif 'Predicted change in Energy=' in var:
            predechange = var.strip().split('=')[-1]
        elif optswitch == 1:

    # Optimization parse

            if 'Optimization completed.' in var:
                optcomplete = 1
            elif 'Optimization completed on the basis of negligible forces.' \
                in var:
                optcomplete = 2
            if printtype == 'none':
                if 'Maximum Force' in var:
                    maxforce = var.rstrip().split()
                    if maxforce[-1] == 'NO':
                        maxforce[-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif maxforce[-1] == 'YES':
                        maxforce[-1] = '\033[1;32mYES\033[1;0m'
                elif 'RMS     Force' in var:
                    rmsforce = var.rstrip().split()
                    if rmsforce[-1] == 'NO':
                        rmsforce[-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif rmsforce[-1] == 'YES':
                        rmsforce[-1] = '\033[1;32mYES\033[1;0m'
                elif 'Maximum Displacement' in var:
                    maxdisp = var.rstrip().split()
                    if maxdisp[-1] == 'NO':
                        maxdisp[-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif maxdisp[-1] == 'YES':
                        maxdisp[-1] = '\033[1;32mYES\033[1;0m'
                elif 'RMS     Displacement' in var:
                    rmsdisp = var.rstrip().split()
                    if rmsdisp[-1] == 'NO':
                        rmsdisp[-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif rmsdisp[-1] == 'YES':
                        rmsdisp[-1] = '\033[1;32mYES\033[1;0m'
                elif 'Step number' in var and 'out of a maximum' in var:
                    optstep = var.split()[2].strip()
            elif printtype == 'long' or printtype == 'extra':
                if 'Step number' in var and 'out of a maximum' in var:
                    optstep = var.split()[2].strip()
                elif 'Maximum Force' in var:
                    maxforce.append(list(var.rstrip().split()))
                    if maxforce[-1][-1] == 'NO':
                        maxforce[-1][-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif maxforce[-1][-1] == 'YES':
                        maxforce[-1][-1] = '\033[1;32mYES\033[1;0m'
                elif 'RMS     Force' in var:
                    rmsforce.append(list(var.rstrip().split()))
                    if rmsforce[-1][-1] == 'NO':
                        rmsforce[-1][-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif rmsforce[-1][-1] == 'YES':
                        rmsforce[-1][-1] = '\033[1;32mYES\033[1;0m'
                elif 'Maximum Displacement' in var:
                    maxdisp.append(list(var.rstrip().split()))
                    if maxdisp[-1][-1] == 'NO':
                        maxdisp[-1][-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif maxdisp[-1][-1] == 'YES':
                        maxdisp[-1][-1] = '\033[1;32mYES\033[1;0m'
                elif 'RMS     Displacement' in var:
                    rmsdisp.append(list(var.rstrip().split()))
                    if rmsdisp[-1][-1] == 'NO':
                        rmsdisp[-1][-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif rmsdisp[-1][-1] == 'YES':
                        rmsdisp[-1][-1] = '\033[1;32mYES\033[1;0m'

    # Frequency and

        if freqswitch == 1:
            if 'Temperature' in var:
                freqtemp = var.strip()
            elif 'Frequencies --' in var:

      # elif 'Vibrational temperatures:' in var:
      #  frequencies = var.strip().split(':')[-1].strip()

                if frequencies[0] == '0.00':
                    frequencies = []
                if '---' in var:
                    frequencies.append(var.split()[1].strip().replace('---'
                            , '-'))
                    frequencies.append(var.split()[2].strip())
                    frequencies.append(var.split()[3].strip())
                elif len(var.split()) == 3:
                    frequencies.append(var.split()[2].strip())
                    frequencies.append('0.000')
                    frequencies.append('0.000')
                    frequencies.append('0.000')
                    frequencies.append('0.000')
                    frequencies.append('0.000')
                elif len(var.split()) == 4:
                    frequencies.append(var.split()[2].strip())
                    frequencies.append(var.split()[3].strip())
                    frequencies.append('0.000')
                    frequencies.append('0.000')
                    frequencies.append('0.000')
                    frequencies.append('0.000')
                else:
                    frequencies.append(var.split()[2].strip())
                    frequencies.append(var.split()[3].strip())
                    frequencies.append(var.split()[4].strip())
            elif 'imaginary frequencies' in var:
                if var.strip().split()[0] == '******':
                    negfreq = int(var.strip().split()[1])
                else:
                    negfreq = int(var.strip().split()[0])
            elif 'Zero-point correction' in var:
                zeropointe = str(var.strip().split()[-2])
            elif 'Thermal correction to Energy' in var:
                therme = str(var.strip().split()[-1])
            elif 'Thermal correction to Enthalpy' in var:
                enthalpye = str(var.strip().split()[-1])
            elif 'Thermal correction to Gibbs Free Energy' in var:
                gibbse = str(var.strip().split()[-1])
            elif 'Sum of electronic and zero-point Energies' in var:
                zeropointesum = str(var.strip().split()[-1])
            elif 'Sum of electronic and thermal Energies' in var:
                thermesum = str(var.strip().split()[-1])
            elif 'Sum of electronic and thermal Enthalpies' in var:
                enthalpyesum = str(var.strip().split()[-1])
            elif 'Sum of electronic and thermal Free Energies' in var:
                gibbsesum = str(var.strip().split()[-1])
                stotal = str(filelist[i + 4].strip().split()[-1])
                selec = str(filelist[i + 5].strip().split()[-1])
                strans = str(filelist[i + 6].strip().split()[-1])
                srot = str(filelist[i + 7].strip().split()[-1])
                svib = str(filelist[i + 8].strip().split()[-1])
        if 'Normal termination of Gaussian' in var and len(filelist) \
            - i < 10:
            jobdone = 1
        elif 'Error termination' in var and i - len(filelist) < 10:
            jobdone = 2

    if trans == 'translate':

    # Printed output list
    # printlines = ['\033[1;44m                                                                              \033[1;0m']

        printlines = ['Analyzing Gaussian Output File: ' + logfile]
        printlines.append('Using ' + str(programversion))
        printlines.append('[#Processors=' + nproc + '  Memory=' + mem
                          + '  CheckPoint=' + chk + ']')
        printlines.append('=============================================================================='
                          )
    elif trans == 'archive':

    # Initial information

        printlines.append('Supporting Information: ' + logfile)
        printlines.append('------------------------------------------------------------------------------'
                          )
        printlines.append('Using: ' + programversion)
        printlines.append('=============================================================================='
                          )

    if trans == 'archive' or trans == 'translate':

    # Print config

        for line in configlines:
            printlines.append(line)

    # Modredundant output

        if modswitch == 1:
            for line in modlines:
                printlines.append(line)

    # Print pointgroup and stoichiometry

        printlines.append('------------------------------------------------------------------------------'
                          )
        printlines.append('Pointgroup=' + point + '\tStoichiometry='
                          + stoic + '\t' + pointstoic)

    # Print Charge Multi and BasisFunc

        printlines.append('Charge = ' + charge + '\tMultiplicity = '
                          + multi + '\tBasis = ' + basisfunc)
        printlines.append('------------------------------------------------------------------------------'
                          )

    # Print Energy

        printlines.append('SCF Energy = ' + scfenergy
                          + '\tPredicted Change = ' + predechange)
        printlines.append('------------------------------------------------------------------------------'
                          )

    # Optimization output

        if optswitch == 1:
            if optcomplete >= 1 and optstep != 0:
                printlines.append('\033[1;32mOptimization completed.\033[0;0m  ----------------------------  Optimization step: '
                                   + optstep)
            elif optcomplete >= 1 and optstep == 0:
                printlines.append('\033[1;32mOptimization completed.\033[0;0m'
                                  )
            elif optcomplete == 0 and optstep != 0:
                printlines.append('\033[1;31mOptimization incomplete.\033[0;0m  ----------------------------  Optimization step: '
                                   + optstep)
            elif optcomplete == 0 and optstep == 0:
                printlines.append('\033[1;31mOptimization incomplete.\033[0;0m'
                                  )
            if optcomplete == 2:
                printlines.append('\033[1;33mOptimization completed on the basis of negligible forces.\033[0;0m'
                                  )
            if len(maxforce) > 0 and printtype == 'none':
                printlines.append('Item      Max Val.    Criteria    Pass?      RMS Val.    Criteria    Pass?'
                                  )
                printlines.append('Force     ' + maxforce[-3] + ' || '
                                  + maxforce[-2] + '    '
                                  + maxforce[-1] + '        '
                                  + rmsforce[-3] + ' || '
                                  + rmsforce[-2] + '    '
                                  + rmsforce[-1])
                printlines.append('Disp      ' + maxdisp[-3] + ' || '
                                  + maxdisp[-2] + '    ' + maxdisp[-1]
                                  + '        ' + rmsdisp[-3] + ' || '
                                  + rmsdisp[-2] + '    ' + rmsdisp[-1])
                printlines.append('------------------------------------------------------------------------------'
                                  )
            elif len(maxforce) > 0 and printtype == 'long':
                printlines.append('Opt Step:      Force [Max || RMS]     Displacement [Max || RMS]'
                                  )
                for i in range(len(maxforce)):
                    printlines.append(' ' + str(i + 1) + '  ['
                            + maxforce[i][-3] + ' (' + maxforce[i][-1]
                            + ') || ' + rmsforce[i][-3] + ' ('
                            + rmsforce[i][-1] + ')]  ['
                            + maxdisp[i][-3] + ' (' + maxdisp[i][-1]
                            + ') || ' + rmsdisp[i][-3] + ' ('
                            + rmsdisp[i][-1] + ')]')

          # printlines.append('Force     ' + maxforce[i][-3] + ' || ' + maxforce[i][-2] + '    ' + maxforce[i][-1] + '        ' + rmsforce[i][-3] + ' || ' + rmsforce[i][-2] + '    ' + rmsforce[i][-1])
          # printlines.append('Disp      ' + maxdisp[i][-3] + ' || ' + maxdisp[i][-2] + '    ' + maxdisp[i][-1] + '        ' + rmsdisp[i][-3] + ' || ' + rmsdisp[i][-2] + '    ' + rmsdisp[i][-1])

                printlines.append('------------------------------------------------------------------------------'
                                  )
            elif len(maxforce) > 0 and printtype == 'extra':
                printlines.append('Opt Step:      Force [Max || RMS]     Displacement [Max || RMS]'
                                  )
                if len(maxforce) > 5:
                    lastlines = len(maxforce) - 5
                else:
                    lastlines = 0
                for i in range(lastlines, len(maxforce)):
                    printlines.append(' ' + str(i + 1) + '  ['
                            + maxforce[i][-3] + ' (' + maxforce[i][-1]
                            + ') || ' + rmsforce[i][-3] + ' ('
                            + rmsforce[i][-1] + ')]  ['
                            + maxdisp[i][-3] + ' (' + maxdisp[i][-1]
                            + ') || ' + rmsdisp[i][-3] + ' ('
                            + rmsdisp[i][-1] + ')]')

          # printlines.append('Force     ' + maxforce[i][-3] + ' || ' + maxforce[i][-2] + '    ' + maxforce[i][-1] + '        ' + rmsforce[i][-3] + ' || ' + rmsforce[i][-2] + '    ' + rmsforce[i][-1])
          # printlines.append('Disp      ' + maxdisp[i][-3] + ' || ' + maxdisp[i][-2] + '    ' + maxdisp[i][-1] + '        ' + rmsdisp[i][-3] + ' || ' + rmsdisp[i][-2] + '    ' + rmsdisp[i][-1])

                printlines.append('------------------------------------------------------------------------------'
                                  )

    if trans == 'archive':

        printlines.append('    Atomic        Coordinates (Angstroms)')
        printlines.append('    Type:            X       Y       Z')
        printlines.append('------------------------------------------------------------------------------'
                          )
        printlines.append('GEOM DATA!')

    if trans == 'translate' or trans == 'archive':

    # Frequency output

        if freqswitch == 1 and len(enthalpye) > 0:
            printlines.append('Statistical Thermodynamic Analysis for '
                              + logfile)
            printlines.append(freqtemp)
            printlines.append('=============================================================================='
                              )
            printlines.append('SCF Energy=                                             '
                               + scfenergy)
            printlines.append('Zero-point correction (ZPE)=            '
                               + zeropointe + '        '
                              + zeropointesum)
            printlines.append('Internal Energy (U)=                    '
                               + therme + '        ' + thermesum)
            printlines.append('Enthalpy (H)=                           '
                               + enthalpye + '        ' + enthalpyesum)
            printlines.append('Gibbs Free Energy (G)=                  '
                               + gibbse + '        ' + gibbsesum)
            printlines.append('Entropy (S) [Cal/Mol-Kelvin]=           '
                               + stotal)
            printlines.append('Entropy (S): electronic=                '
                               + selec)
            printlines.append('Entropy (S): translational=             '
                               + strans)
            printlines.append('Entropy (S): rotational=                '
                               + srot)
            printlines.append('Entropy (S): vibrational=               '
                               + svib)
            printlines.append('Frequencies: ' + frequencies[0] + '  '
                              + frequencies[1] + '  ' + frequencies[2]
                              + '  ' + frequencies[3] + '  '
                              + frequencies[4])
            if negfreq == 1:
                printlines.append('\033[1;33m              !!!Warning there were '
                                   + str(negfreq)
                                  + ' negative frequencies!!!\033[0;m')
            elif negfreq >= 2:
                printlines.append('\033[1;31m              !!!Warning there were '
                                   + str(negfreq)
                                  + ' negative frequencies!!!\033[0;m')
            printlines.append('------------------------------------------------------------------------------'
                              )

    if trans == 'translate':

    # Print end of file for long and extra

        if printtype == 'extra' or printtype == 'long':
            for i in range(len(filelist) - 10, len(filelist)):
                printlines.append(filelist[i].strip())

    # Completed job

        if jobdone == 1:
            printlines.append('=============================================================================='
                              )
            printlines.append('\033[1;32m##################              JOB COMPLETED               ##################\033[0;0m'
                              )
            printlines.append('=============================================================================='
                              )
        elif jobdone == 2:
            printlines.append('=============================================================================='
                              )
            printlines.append('\033[1;31m!!!!!!!!!!!!!!!!!!            !!!JOB FAILED!!!              !!!!!!!!!!!!!!!!!!\033[0;0m'
                              )
            printlines.append('=============================================================================='
                              )
    elif trans == 'list':

        line = scfenergy + '  '
        line += zeropointesum + '  ' + thermesum + '  ' + enthalpyesum \
            + '  ' + gibbsesum + '  '
        if int(negfreq) == 0:
            line += '\033[1;32m' + frequencies[0] + '  ' \
                + frequencies[1] + '\033[1;0m' + '  '
        elif int(negfreq) == 1:
            line += '\033[1;33m' + frequencies[0] + '  ' \
                + frequencies[1] + '\033[1;0m' + '  '
        elif int(negfreq) >= 2:
            line += '\033[1;31m' + frequencies[0] + '  ' \
                + frequencies[1] + '\033[1;0m' + '  '
        line += basisfunc + '  ' + pointstoic + '  ' + stoic + '  '
        if jobdone == 0:
            line += logfile.split('/')[-1]
        elif jobdone == 1:
            filename = logfile.split('/')[-1]
            line += '\033[1;32m' + filename + '\033[0;0m'
        elif jobdone == 2:
            filename = logfile.split('/')[-1]
            line += '\033[1;31m' + filename + '\033[0;0m'
        printlines.append(line)

    if trans != 'archive':
        return printlines
    if trans == 'archive':
        tablelines = []
        tablelines.append(logfile.split('/')[-1] + ',' + scfenergy + ','
                           + zeropointe + ',=' + therme + '+RC[-2],='
                          + enthalpye + '+RC[-3],=' + gibbse
                          + '+RC[-4],=RC[-1]-MIN(C[-1]),=RC[-1]*627.5095'
                          )
        tablelines.append(logfile.split('/')[-1] + ',' + scfenergy + ','
                           + zeropointe + ',' + thermesum + ','
                          + enthalpyesum + ',' + gibbsesum
                          + ',=RC[-1]-MIN(C[-1]),=RC[-1]*627.5095')
        tablelines.append(logfile.split('/')[-1] + ',' + scfenergy)
        return (printlines, tablelines)


### --- Jaguar Specific Functions --- ###
###==================================================================================================================###
### --- Parse the Jaguar OUT file for geometry information --- ###

def parse_output_jaguar(ifile):

  # ## --- Open parent file --- ###

    f = open(ifile, 'r')
    ifilelist = f.readlines()
    f.close()

  # ## --- Create some variables needed for parsing the input file --- ###
  # ifilelength is used to gauge what line we are on and also
  # allows us to determine the number of lines in the file
  # ifilelol is a list of lists (lol) that will hold all the data of the input file
  # alist is a temp list that will be appended to the ifilelol
  # ###############################

    ifilelength = 0
    ifilelol = []
    alist = []

  # ### --- Input/parse the input file into a list of lists called ifilelol --- ###

    collectgeom = 0
    for (i, var) in enumerate(ifilelist):

    # Grab charge and multiplicity

        if 'net molecular charge:' in var:
            charge = var.rstrip().split()[-1]
        elif 'multiplicity' in var:
            multi = var.rstrip().split()[-1]
        elif 'SCFE:' in var:
            finalenergy = var.strip().split('hartrees')[0].split()[-1]
        elif 'new geometry:' in var or 'final geometry:' in var:
            collectgeom = 1
            numatoms = 0
            ifilelol = []
        elif len(var.strip()) == 0:

      # alist = []

            collectgeom = 0
        elif collectgeom == 1 and 'angstroms' not in var and 'atom' \
            not in var:
            alist.append(var)
            numatoms += 1
    finalgeom = [0] * (numatoms + 2)
    finalgeom[0] = str(numatoms)
    finalgeom[1] = 'Energy : ' + str(finalenergy)
    for i in range(len(alist) - numatoms, len(alist)):
        line = alist[i].split()
        j = i - len(alist) + numatoms + 2
        finalgeom[j] = Atom()
        finalgeom[j].e = ''.join([k for k in line[0]
                                 if not k.isdigit()])
        finalgeom[j].en = get_element_num(finalgeom[j].e)
        finalgeom[j].x = float(line[1])
        finalgeom[j].y = float(line[2])
        finalgeom[j].z = float(line[3])
    return (finalgeom, charge, multi)


### --- Parse the Jaguar OUT file for geometry information --- ###

def parse_output_jaguar_traj(ifile):

  # ## --- Open parent file --- ###

    f = open(ifile, 'r')
    ifilelist = f.readlines()
    f.close()

  # ## --- Create some variables needed for parsing the input file --- ###
  # ifilelength is used to gauge what line we are on and also
  # allows us to determine the number of lines in the file
  # ifilelol is a list of lists (lol) that will hold all the data of the input file
  # alist is a temp list that will be appended to the ifilelol
  # ###############################

    alist = []

  # ### --- Input/parse the input file into a list of lists called ifilelol --- ###

    collectgeom = 0
    allenergy = []
    for (i, var) in enumerate(ifilelist):

    # Grab charge and multiplicity

        if 'net molecular charge:' in var:
            charge = var.rstrip().split()[-1]
        elif 'multiplicity' in var:
            multi = var.rstrip().split()[-1]
        elif 'SCFE:' in var:
            allenergy.append(var.strip().split('hartrees'
                             )[0].split()[-1])
        elif 'new geometry:' in var or 'final geometry:' in var:
            collectgeom = 1
            numatoms = 0
            ifilelol = []
        elif len(var.strip()) == 0:

      # alist = []

            collectgeom = 0
        elif collectgeom == 1 and 'angstroms' not in var and 'atom' \
            not in var:
            alist.append(var)
            numatoms += 1
    trajectorygeom = [0] * (len(alist) / numatoms * (numatoms + 2))
    skip = 0
    for i in range(len(alist)):
        if i % numatoms == 0 and i / numatoms != len(allenergy):
            trajectorygeom[i + skip] = str(numatoms)
            trajectorygeom[i + skip + 1] = 'Energy: ' + str(allenergy[i
                    / numatoms])
            skip += 2
        elif i % numatoms == 0 and i / numatoms >= len(allenergy):
            trajectorygeom[i + skip] = str(numatoms)
            trajectorygeom[i + skip + 1] = 'Energy: ' \
                + str(allenergy[-1])
            skip += 2
        line = alist[i].split()
        trajectorygeom[i + skip] = Atom()
        trajectorygeom[i + skip].e = ''.join([j for j in line[0]
                if not j.isdigit()])
        trajectorygeom[i + skip].en = get_element_num(trajectorygeom[i
                + skip].e)
        trajectorygeom[i + skip].x = float(line[1])
        trajectorygeom[i + skip].y = float(line[2])
        trajectorygeom[i + skip].z = float(line[3])
    return (trajectorygeom, charge, multi)


### --- Out put the Jaguar input files --- ###

def output_jaguar(
    ifile,
    config,
    charge,
    multi,
    writemod,
    modred,
    writeecp,
    ecplines,
    writeclose,
    closelines,
    ifilelol,
    ):

    ofiletmp = ifile.split('.')
    ofile = ''
    for i in range(len(ofiletmp)):
        if i != len(ofiletmp) - 1:
            ofile += ofiletmp[i] + '.'
        elif i == len(ofiletmp) - 1:
            ofile += 'in'

  # print ofile

    f = open('NEWJOBS.JAG/' + ofile, 'w+')
    f.write('&gen\n')
    for (i, var) in enumerate(config):
        f.write(var + '\n')
    f.write('molchg=' + str(charge) + '\n')
    f.write('multip=' + str(multi) + '\n')
    f.write('&\n')
    f.write('entry_name: ' + ifile + '\n')
    f.write('&zmat\n')
    for i in range(2, len(ifilelol)):
        line = str(ifilelol[i].e) + str(i - 1) + ' ' \
            + str(ifilelol[i].x) + ' ' + str(ifilelol[i].y) + ' ' \
            + str(ifilelol[i].z)
        f.write(line + '\n')
    f.write('&\n')
    if writemod == 1:
        f.write('&coord\n')
        if len(modred) > 0 and modred[0] != '':
            for (i, var) in enumerate(modred):
                f.write(var + '\n')
        else:
            f.write('MODRED!\n')
        f.write('&\n')
    f.close()
    return


# Parse the Jaguar output file and print data to the screen

def translate_jaguar_output(logfile, trans, printtype):
    initlines = \
        ['\033[1;34mPython Jaguar Output File Parsing Script:   nJT.py\033[1;0m'
         ,
         '\033[1;34mScripted by Thomas J. L. Mustard, O. Maduka Ogba, Paul Ha-Yeon Cheong\033[1;0m'
         ,
         '------------------------------------------------------------------------------'
         ]
    if trans == 'translate':
        for line in initlines:
            print(line)

  # Global Variables

    config = ''
    configswitch = 0
    printlines = []
    programversion = ''
    procs = ''
    tmpfiles = ''
    configlines = []
    configlinesall = ''
    modswitch = 0
    modlines = []
    pointstoic = ''
    point = ''
    stoic = ''
    chargemulti = ''
    charge = ''
    multi = ''
    basisfunc = ''
    scfenergy = ''
    predechange = ''
    freqTemp = ''
    optswitch = 0
    optcomplete = 0
    optstep = 0
    if printtype == 'none':
        maxforce = ''
        rmsforce = ''
        maxdisp = ''
        rmsdisp = ''
        energycon = [' ---N/A--- ', 'N/A', '(', '0.000000', ')']
    else:
        maxforce = []
        rmsforce = []
        maxdisp = []
        rmsdisp = []
        energycon = [' ---N/A--- ', 'N/A', '(', '0.000000', ')']
    freqswitch = 0
    freqtemp = ''
    freqpres = ''
    if trans == 'list':
        frequencies = ['0.00', '0.00']
        negfreq = 0
        zeropointe = '0.000000'
        therme = '0.000000'
        enthalpye = '0.000000'
        gibbse = '0.000000'
        zeropointesum = '0.000000'
        thermesum = '0.000000'
        enthalpyesum = '0.000000'
        gibbsesum = '0.000000'
        stotal = '0.000000'
        selec = '0.000000'
        srot = '0.000000'
        strans = '0.000000'
        svib = '0.000000'
    else:
        frequencies = []
        negfreq = 0
        zeropointe = ''
        therme = ''
        enthalpye = ''
        gibbse = ''
        zeropointesum = ''
        thermesum = ''
        enthalpyesum = ''
        gibbsesum = ''
        stotal = ''
        selec = ''
        srot = ''
        strans = ''
        svib = ''
    jobdone = 0
    memswitch = 0
    f = open(logfile)
    filelist = f.readlines()
    f.close()
    for (i, var) in enumerate(filelist):

    # if i < 30 and 'Jaguar version' not in var and jaguartext == 0:
    #  print "Error! Not a Jaguar out file."
    # exit(0)
    # Find the Gaussian Version

        if 'Jaguar version' in var:
            programversion = var.rstrip().split(',')[0]
        elif 'Temporary files :' in var:

    # Find the tmp file

            tmpfiles = var.split(':')[1].rstrip()
        elif 'Using up to ' in var:

    # Find the processors/threads used

            procs = var.split('to')[1].rstrip()

    # Get pointgroup and stoic

        if 'Point Group used:' in var:
            point = var.split()[-1]
        elif 'Stoichiometry:' in var:
            stoic = var.split()[-1]
        elif 'net molecular charge:' in var:

    # Get charge and multiplicity

            charge = var.rstrip().split()[-1]
        elif 'multiplicity' in var:
            multi = var.rstrip().split()[-1]
        elif 'basis set:' in var:

    # Find basis functions

            basisfunc = var.strip().split()[-1]

    # Find the job configuration used

        if 'Non-default options chosen:' in var:
            for j in range(1, 30):
                if 'Temporary integer options' not in filelist[i + j]:
                    configlines.append(filelist[i + j].strip())
                    if 'optimized' in filelist[i + j]:
                        optswitch = 1
                    if 'modredundant' in filelist[i + j]:
                        modswitch = 1
                    if 'frequencies' in filelist[i + j]:
                        freqswitch = 1
                else:
                    break
        elif modswitch == 1:

    # Check to see if modredundant is in the config

            if 'The following ModRedundant input section has been read:' \
                in var:
                for j in range(1, 20):
                    if len(filelist[i + j]) > 2:
                        modlines.append(filelist[i + j].strip())
                    else:
                        break
        elif 'SCFE:' in var:

    # Get SCF Energy and predicted change

            scfenergy = var.strip().split('hartrees')[0].split()[-1]
        elif 'predicted energy change:' in var:
            predechange = var.strip().split(':')[-1]
        elif optswitch == 1:

    # Optimization parse

            if 'Geometry optimization was OK.' in var:
                optcomplete = 1
            elif 'geometry optimization step' in var:
                optstep = var.split()[-1]
            if printtype == 'none':
                if 'gradient maximum:' in var:
                    maxforce = var.rstrip().split(':')[-1].split()
                    if maxforce[1] != '.':
                        maxforce[1] = '\033[1;32mYES\033[0;m'
                    else:
                        maxforce[1] = '\033[1;31mNO!\033[0;m'
                elif 'gradient rms:' in var:
                    rmsforce = var.rstrip().split(':')[-1].split()
                    if rmsforce[1] != '.':
                        rmsforce[1] = '\033[1;32mYES\033[0;m'
                    else:
                        rmsforce[1] = '\033[1;31mNO!\033[0;m'
                elif 'displacement maximum:' in var:
                    maxdisp = var.rstrip().split(':')[-1].split()
                    if maxdisp[1] != '.':
                        maxdisp[1] = '\033[1;32mYES\033[0;m'
                    else:
                        maxdisp[1] = '\033[1;31mNO!\033[0;m'
                elif 'displacement rms:' in var:
                    rmsdisp = var.rstrip().split(':')[-1].split()
                    if rmsdisp[1] != '.':
                        rmsdisp[1] = '\033[1;32mYES\033[0;m'
                    else:
                        rmsdisp[1] = '\033[1;31mNO!\033[0;m'
                elif 'energy change:' in var and 'predicted' not in var:
                    energycon = var.rstrip().split(':')[-1].split()
                    if energycon[1] != '.':
                        energycon[1] = '\033[1;32mYES\033[0;m'
                    else:
                        energycon[1] = '\033[1;31mNO!\033[0;m'
            elif printtype == 'extra' or printtype == 'long':
                if 'gradient maximum:' in var:

          # Fix the fact that sometimes the energy convergence data is not printed out for the first step

                    if 'energy change:' not in filelist[i - 1]:
                        energycon.append([' ---N/A--- ', 'N/A', '(',
                                '0.000000', ')'])
                    maxforce.append(var.rstrip().split(':')[-1].split())
                    if maxforce[-1][1] != '.':
                        maxforce[-1][1] = '\033[1;32mYES\033[0;m'
                    else:
                        maxforce[-1][1] = '\033[1;31mNO!\033[0;m'
                elif 'gradient rms:' in var:
                    rmsforce.append(var.rstrip().split(':')[-1].split())
                    if rmsforce[-1][1] != '.':
                        rmsforce[-1][1] = '\033[1;32mYES\033[0;m'
                    else:
                        rmsforce[-1][1] = '\033[1;31mNO!\033[0;m'
                elif 'displacement maximum:' in var:
                    maxdisp.append(var.rstrip().split(':')[-1].split())
                    if maxdisp[-1][1] != '.':
                        maxdisp[-1][1] = '\033[1;32mYES\033[0;m'
                    else:
                        maxdisp[-1][1] = '\033[1;31mNO!\033[0;m'
                elif 'displacement rms:' in var:
                    rmsdisp.append(var.rstrip().split(':')[-1].split())
                    if rmsdisp[-1][1] != '.':
                        rmsdisp[-1][1] = '\033[1;32mYES\033[0;m'
                    else:
                        rmsdisp[-1][1] = '\033[1;31mNO!\033[0;m'
                elif 'energy change:' in var and 'predicted' not in var:
                    energycon.append(var.rstrip().split(':'
                            )[-1].split())
                    if energycon[-1][1] != '.':
                        energycon[-1][1] = '\033[1;32mYES\033[0;m'
                    else:
                        energycon[-1][1] = '\033[1;31mNO!\033[0;m'

    # Frequency and

        if freqswitch == 1:
            if 'Thermochemical properties at' in var:
                freqpres = var.strip().split('s at')[-1].strip()
            elif 'T =' in var:
                freqtemp = var.strip().split('=')[-1].strip()
            elif 'frequencies' in var and 'frequencies' \
                == var.split()[0]:
                if len(frequencies) != 0 and frequencies[0] == '0.00':
                    frequencies = []
                if len(frequencies) < 5:
                    for i in range(1,len(var.split())):
                        frequencies.append(var.split()[i])
                    #frequencies.append(var.split()[2])
                    #frequencies.append(var.split()[3])
                    #frequencies.append(var.split()[4])
                    #frequencies.append(var.split()[5])
                if len(frequencies) <= 1:
                    frequencies.append('0.00')
                    frequencies.append('0.00')
                    frequencies.append('0.00')
                    frequencies.append('0.00')

            elif 'Number of imaginary frequencies:' in var:

      # elif 'vibrational temperatures:' in var:
      #  frequencies = filelist[i+2].strip().split(':')[-1].strip()

                negfreq = int(var.strip().split(':')[-1])
            elif 'The zero point energy (ZPE):' in var:
                zeropointe = str(var.strip().split(':')[-1])
            elif 'Total internal energy, Utot (SCFE + ZPE + U):' in var:
                thermesum = str(var.strip().split(':')[-1])
                stotal = str(filelist[i - 2].strip().split()[3])
                selec = str(filelist[i - 3].strip().split()[3])
                strans = str(filelist[i - 6].strip().split()[3])
                srot = str(filelist[i - 5].strip().split()[3])
                svib = str(filelist[i - 4].strip().split()[3])
                therme = str(filelist[i - 2].strip().split()[1]).strip()
                enthalpye = str(filelist[i
                                - 2].strip().split()[4]).strip()
                gibbse = str(filelist[i - 2].strip().split()[5]).strip()
            elif 'Total enthalpy, Htot (Utot + pV):' in var:
                enthalpyesum = str(var.strip().split(':')[-1])
            elif 'Total Gibbs free energy, Gtot (Htot - T*S):' in var:
                gibbsesum = str(var.strip().split(':')[-1])
        if 'Job' in var and 'completed on' in var and i - len(filelist) \
            < 10:
            jobdone = 1

    if trans == 'translate':

    # Printed output list

        printlines = ['Analyzing Jaguar Output File: ' + logfile]
        printlines.append('Using ' + str(programversion))
        printlines.append('[#Processors=' + procs + '  Temp Files='
                          + tmpfiles + ']')
        printlines.append('=============================================================================='
                          )
    elif trans == 'archive':

    # Initial information

        printlines.append('Supporting Information: ' + logfile)
        printlines.append('------------------------------------------------------------------------------'
                          )
        printlines.append('Using: ' + programversion)
        printlines.append('=============================================================================='
                          )

    if trans == 'translate' or trans == 'archive':

    # Print config

        configlinestmp = ' '.join(str(e) for e in configlines)
        segments = list()
        while len(configlinestmp) > 0:
            segments.append(configlinestmp[:78])
            configlinestmp = configlinestmp[78:]
        for line in segments:
            printlines.append(line)

    # printlines.append(' '.join(str(e) for e in configlines))

    # Modredundant output

        if modswitch == 1:
            for line in modlines:
                printlines.append(line)

    # Print pointgroup and stoichiometry

        printlines.append('------------------------------------------------------------------------------'
                          )
        printlines.append('Pointgroup = ' + point + '\tStoichiometry = '
                           + stoic)

    # Print Charge Multi and BasisFunc

        printlines.append('Charge = ' + charge + '\tMultiplicity = '
                          + multi + '\tBasis = ' + basisfunc)
        printlines.append('------------------------------------------------------------------------------'
                          )

    # Print Energy

        printlines.append('SCF Energy = ' + scfenergy
                          + '\tPredicted Change = ' + predechange)
        printlines.append('------------------------------------------------------------------------------'
                          )

    # Optimization output

        if optswitch == 1:
            if optcomplete == 1:
                printlines.append('\033[1;32mOptimization completed.\033[0;m  ----------------------------  Optimization step: '
                                   + str(optstep))
            else:
                printlines.append('\033[1;31mOptimization incomplete.\033[0;m  ----------------------------  Optimization step: '
                                   + str(optstep))
            if len(maxforce) > 0 and printtype == 'none':
                printlines.append('Item      Max Val.    Criteria    Pass?      RMS Val.    Criteria    Pass?'
                                  )
                printlines.append('Energy    ' + energycon[-5] + ' || '
                                  + energycon[-2] + '   '
                                  + energycon[-4])
                printlines.append('Force     ' + maxforce[-5] + ' || '
                                  + maxforce[-2] + '   ' + maxforce[-4]
                                  + '       ' + rmsforce[-5] + ' || '
                                  + rmsforce[-2] + '   ' + rmsforce[-4])
                printlines.append('Disp      ' + maxdisp[-5] + ' || '
                                  + maxdisp[-2] + '   ' + maxdisp[-4]
                                  + '       ' + rmsdisp[-5] + ' || '
                                  + rmsdisp[-2] + '   ' + rmsdisp[-4])
                printlines.append('------------------------------------------------------------------------------'
                                  )
            elif len(maxforce) > 0 and printtype == 'long':
                printlines.append('Opt Step:   Energy       Force [Max || RMS]          Displacement [Max || RMS]'
                                  )
                for i in range(len(maxforce)):
                    printlines.append(' ' + str(i + 1) + '  '
                            + energycon[i][-5] + ' ('
                            + energycon[i][-4] + ')  ['
                            + maxforce[i][-5] + ' (' + maxforce[i][-4]
                            + ') || ' + rmsforce[i][-5] + ' ('
                            + rmsforce[i][-4] + ')]  ['
                            + maxdisp[i][-5] + ' (' + maxdisp[i][-4]
                            + ') || ' + rmsdisp[i][-5] + ' ('
                            + rmsdisp[i][-4] + ')]')
                printlines.append('------------------------------------------------------------------------------'
                                  )
            elif len(maxforce) > 0 and printtype == 'extra':
                printlines.append('Opt Step:   Energy       Force [Max || RMS]          Displacement [Max || RMS]'
                                  )
                if len(maxforce) > 5:
                    lastlines = len(maxforce) - 5
                else:
                    lastlines = 0
                for i in range(lastlines, len(maxforce)):
                    printlines.append(' ' + str(i + 1) + '  '
                            + energycon[i][-5] + ' ('
                            + energycon[i][-4] + ')  ['
                            + maxforce[i][-5] + ' (' + maxforce[i][-4]
                            + ') || ' + rmsforce[i][-5] + ' ('
                            + rmsforce[i][-4] + ')]  ['
                            + maxdisp[i][-5] + ' (' + maxdisp[i][-4]
                            + ') || ' + rmsdisp[i][-5] + ' ('
                            + rmsdisp[i][-4] + ')]')
                printlines.append('------------------------------------------------------------------------------'
                                  )

    if trans == 'archive':
        printlines.append('    Atomic        Coordinates (Angstroms)')
        printlines.append('    Type:            X       Y       Z')
        printlines.append('------------------------------------------------------------------------------'
                          )
        printlines.append('GEOM DATA!')

    if trans == 'translate' or trans == 'archive':

    # Frequency output

        if freqswitch == 1 and len(frequencies) > 0:
            printlines.append('Statistical Thermodynamic Analysis for '
                              + logfile)
            printlines.append('Computed at:  ' + freqtemp + '  and  '
                              + freqpres)
            printlines.append('=============================================================================='
                              )
            printlines.append('SCF Energy=                                             '
                               + scfenergy)
            printlines.append('Zero-point correction (ZPE)=                         '
                               + zeropointe)
            printlines.append('Internal Energy (U)=                    '
                               + therme + '       ' + thermesum)
            printlines.append('Enthalpy (H)=                           '
                               + enthalpye + '       ' + enthalpyesum)
            printlines.append('Gibbs Free Energy (G)=                  '
                               + gibbse + '       ' + gibbsesum)
            printlines.append('Entropy (S) [Cal/Mol-Kelvin]=                            '
                               + stotal)
            printlines.append('Entropy (S): electronic=                                 '
                               + selec)
            printlines.append('Entropy (S): translational=                              '
                               + strans)
            printlines.append('Entropy (S): rotational=                                 '
                               + srot)
            printlines.append('Entropy (S): vibrational=                                '
                               + svib)

            for i, var in enumerate(frequencies):
                freqTemp += var + '   '
            printlines.append('Frequencies:  ' + freqTemp)
                              #frequencies[0] + '   '
                              #+ frequencies[1] + '   ' + frequencies[2]
                              #+ '   ' + frequencies[3] + '   '
                              #+ frequencies[4])
            if negfreq == 1:
                printlines.append('\033[1;33m              !!!Warning there were '
                                   + str(negfreq)
                                  + ' negative frequencies!!!\033[0;m')
            elif negfreq >= 2:
                printlines.append('\033[1;31m              !!!Warning there were '
                                   + str(negfreq)
                                  + ' negative frequencies!!!\033[0;m')

    if trans == 'translate':

    # Print end of file for long and extra

        if printtype == 'extra' or printtype == 'long':
            printlines.append('------------------------------------------------------------------------------'
                              )
            for i in range(len(filelist) - 10, len(filelist)):
                printlines.append(filelist[i].strip())

    # Completed job

        if jobdone == 1:
            printlines.append('=============================================================================='
                              )
            printlines.append('\033[1;32m##################              JOB COMPLETED               ##################\033[0;m'
                              )
            printlines.append('=============================================================================='
                              )

    if trans == 'list':
        line = scfenergy + '  '
        line += zeropointe.split()[0] + '  ' + thermesum.split()[0] \
            + '  ' + enthalpyesum.split()[0] + '  ' \
            + gibbsesum.split()[0] + '  '
        if negfreq == 0:
            line += '\033[1;32m' + frequencies[0] + '  ' \
                + frequencies[1] + '\033[1;0m' + '  '
        elif negfreq == 1:
            line += '\033[1;33m' + frequencies[0] + '  ' \
                + frequencies[1] + '\033[1;0m' + '  '
        elif negfreq >= 2:
            line += '\033[1;31m' + frequencies[0] + '  ' \
                + frequencies[1] + '\033[1;0m' + '  '
        line += basisfunc + '  ' + pointstoic + '  ' + stoic + '  '
        if jobdone == 0:
            line += logfile.split('/')[-1]
        elif jobdone == 1:
            filename = logfile.split('/')[-1]
            line += '\033[1;32m' + filename + '\033[0;0m'
        elif jobdone == 2:
            filename = logfile.split('/')[-1]
            line += '\033[1;31m' + filename + '\033[0;0m'
        printlines.append(line)

    # if negfreq > 0 and trans != "archive":
    #    printlines.append('\033[1;31m              !!!Warning there were ' + negfreq + ' negative frequencies!!!\033[0;m')

    if trans != 'archive':
        return printlines
    if trans == 'archive':
        tablelines = []
        tablelines.append(logfile.split('/')[-1] + ',' + scfenergy + ','
                           + zeropointe + ',=' + therme + '+RC[-2],='
                          + enthalpye + '+RC[-3],=' + gibbse
                          + '+RC[-4],=RC[-1]-MIN(C[-1]),=RC[-1]*627.5095'
                          )
        tablelines.append(logfile.split('/')[-1] + ',' + scfenergy + ','
                           + zeropointe + ',' + thermesum + ','
                          + enthalpyesum + ',' + gibbsesum
                          + ',=RC[-1]-MIN(C[-1]),=RC[-1]*627.5095')
        tablelines.append(logfile.split('/')[-1] + ',' + scfenergy)
        return (printlines, tablelines)


### --- Turbomole Specific Functions --- ###
###==================================================================================================================###
### --- Parse the Turbomole LOG file for geometry information --- ###

def parse_output_tm(ifile):

  # ## --- Open parent file --- ###

    f = open(ifile, 'r')
    ifilelist = f.readlines()
    f.close()

  # ## --- Create some variables needed for parsing the input file --- ###
  # ifilelength is used to gauge what line we are on and also
  # allows us to determine the number of lines in the file
  # ifilelol is a list of lists (lol) that will hold all the data of the input file
  # alist is a temp list that will be appended to the ifilelol
  # ###############################

    charge = 0
    multi = -1
    alist = []
    geom = 0
    numatoms = 0
    finalenergy = 0
    finalgeom = []
    for (i, line) in enumerate(ifilelist):
        if line.strip() != '':
            if line.split()[0] == '|':
                if line.split()[1] == 'total':
                    if line.split()[2] == 'energy':
                        finalenergy = float(line.split()[4])
            if line.split()[0] == 'total:':
                numatoms = int(line.split()[1])
            if line.split()[0] == 'center':
                geom = 0
            if geom == 1:
                geom += 1
            if line.split()[0] == 'atomic':
                if line.split()[1] == 'coordinates':
                    geom = 1
            if geom >= 2:
                alist.append(line)
            elif 'charge' in line \
                and '------------------------------------------------------------------------------' \
                in ifilelist[i + 1]:
                charge = ifilelist[i + 2].strip().split()[-1]
            elif 'molecular charge' in line:
                charge = line.rstrip().split()[-1]
            elif 'real representations' in line:
                multi = line.strip().split()[2]
    alistdedup = [x for x in alist if x != '']
    finalgeom = [0] * (numatoms + 2)
    finalgeom[0] = numatoms
    finalgeom[1] = 'Energy: ' + str(finalenergy)
    for i in range(len(alistdedup) - numatoms, len(alistdedup)):
        line = alistdedup[i].split()
        j = i - len(alistdedup) + numatoms + 2
        finalgeom[j] = Atom()
        finalgeom[j].e = line[3]
        finalgeom[j].en = get_element_num(line[3])
        finalgeom[j].x = float(line[0]) / auToAngstrom
        finalgeom[j].y = float(line[1]) / auToAngstrom
        finalgeom[j].z = float(line[2]) / auToAngstrom
    return (finalgeom, charge, multi)


### --- Parse the Turbomole LOG file for geometry information --- ###

def parse_output_tm_traj(ifile):

  # ## --- Open parent file --- ###

    f = open(ifile, 'r')
    ifilelist = f.readlines()
    f.close()

  # ## --- Create some variables needed for parsing the input file --- ###
  # ifilelength is used to gauge what line we are on and also
  # allows us to determine the number of lines in the file
  # ifilelol is a list of lists (lol) that will hold all the data of the input file
  # alist is a temp list that will be appended to the ifilelol
  # ###############################

    charge = 0
    multi = -1
    alist = []
    geom = 0
    skipgeom = True
    numatoms = 0
    allenergy = []
    trajectorygeom = []
    for (i, line) in enumerate(ifilelist):
        if line.strip() != '':
            if line.split()[0] == 'ENERGY' and line.split()[1] == '=':
                allenergy.append(float(line.split()[2]))
            elif line.split()[0] == 'total:':
                numatoms = int(line.split()[1])
            elif line.split()[0] == 'center' or '************' in line:
                geom = 0
            elif 'S T A T P T' in line:
                skipgeom = False
            elif line.split()[0] == 'ATOM' and line.split()[1] \
                == 'CARTESIAN' and line.split()[2] == 'COORDINATES':

        # if line.split()[1] == "coordinates":

                if not skipgeom:
                    geom = 1
                skipgeom = True
            elif geom == 1:
                alist.append(line)
            elif 'charge' in line \
                and '------------------------------------------------------------------------------' \
                in str(ifilelist[i + 1]):
                charge = ifilelist[i + 2].strip().split()[-1]
            elif 'molecular charge' in line:
                charge = line.rstrip().split()[-1]
            elif 'real representations' in line:
                multi = line.strip().split()[2]
    alistdedup = [x for x in alist if x != '']
    trajectorygeom = [0] * (len(alistdedup) / numatoms * (numatoms + 2))
    skip = 0
    for i in range(len(alistdedup)):
        if i % numatoms == 0 and i / numatoms != len(allenergy):
            trajectorygeom[i + skip] = str(numatoms)
            trajectorygeom[i + skip + 1] = 'Energy: ' + str(allenergy[i
                    / numatoms])
            skip += 2
        elif i % numatoms == 0 and i / numatoms >= len(allenergy):
            trajectorygeom[i + skip] = str(numatoms)
            trajectorygeom[i + skip + 1] = 'Energy: ' \
                + str(allenergy[-1])
            skip += 2
        line = alistdedup[i].split()
        trajectorygeom[i + skip] = Atom()
        trajectorygeom[i + skip].en = get_element_num(line[1])
        trajectorygeom[i + skip].e = \
            get_element_name(get_element_num(line[1]))
        trajectorygeom[i + skip].x = float(line[2]) / auToAngstrom
        trajectorygeom[i + skip].y = float(line[3]) / auToAngstrom
        trajectorygeom[i + skip].z = float(line[4]) / auToAngstrom
    return (trajectorygeom, charge, multi)


### --- Out put the TM input files --- ###

def output_tm(
    ifile,
    config,
    charge,
    multi,
    writemod,
    modred,
    writeecp,
    ecplines,
    writeclose,
    closelines,
    ifilelol,
    ):

    ofiletmp = ifile.split('.')
    ofile = ''
    currentdir = os.getcwd()
    for i in range(len(ofiletmp)):
        if i != len(ofiletmp) - 1:
            ofile += ofiletmp[i] + '.'
            ifilebase = ofiletmp[i]
            outdir = currentdir + '/NEWJOBS.TM/' + ofiletmp[i] \
                + '.Turbomole'
        elif i == len(ofiletmp) - 1:
            ofile += 'tmol'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    f = open(outdir + '/' + ofile, 'w+')
    f.write('$coord\n')
    for i in range(2, len(ifilelol)):
        if writemod == 1 and len(modred) > 0 and modred[0] != '':
            for var in modred:
                if var == i - 1:
                    line = str(ifilelol[i].x * auToAngstrom)[0:11] \
                        + '   ' + str(ifilelol[i].y
                            * auToAngstrom)[0:11] + '   ' \
                        + str(ifilelol[i].z * auToAngstrom)[0:11] \
                        + '   ' + str(ifilelol[i].e).lower() + '  f'
        else:
            line = str(ifilelol[i].x * auToAngstrom)[0:11] + '   ' \
                + str(ifilelol[i].y * auToAngstrom)[0:11] + '   ' \
                + str(ifilelol[i].z * auToAngstrom)[0:11] + '   ' \
                + str(ifilelol[i].e).lower()
        f.write(line + '\n')
    f.close()

  # Build define.in file for the define program

    f = open(outdir + '/define.in', 'w+')
    for (i, var) in enumerate(config):
        if 'CHARGE!' in var:
            f.write(str(charge) + '\n')
        elif 'TMOL!' in var:
            f.write('a coord\n')
        elif 'TITLE!' in var:

      # f.write("a " + ofile + "\n")

            f.write(ifilebase + '\n')
        else:
            f.write(var + '\n')
    f.close()

  # Run the define program on define.in and output define.out

    os.chdir(outdir)
    os.system('ln -s ' + ofile + ' coord')
    os.system('define < define.in > define.out')
    f = open('control', 'r')
    control = f.readlines()
    f.close()
    f = open('control', 'w')
    for line in control:
        if '$end' in line:
            f.write('$statpt\n')
            f.write('   itrvec 0\n')
            f.write('   hssfreq 0\n')
            f.write('$maxcor    1000\n')
            f.write('$rpacor    1000\n')
            f.write(line)
        elif '$rij' in line:
            f.write('$ricore_slave    1\n')
            f.write(line)
        else:
            f.write(line)
    f.close()
    os.chdir(currentdir)
    return


# Parse the Turbomole output file and print data to the screen

def translate_tm_output(logfile, trans, printtype):
    initlines = \
        ['\033[1;34mPython Turbomole Output File Parsing Script:   nTMT.py\033[1;0m'
         ,
         '\033[1;34mScripted by Thomas J. L. Mustard, O. Maduka Ogba, Paul Ha-Yeon Cheong\033[1;0m'
         ,
         '------------------------------------------------------------------------------'
         ]
    if trans == 'translate':
        for line in initlines:
            print(line)

  # Global Variables

    config = ''
    configswitch = 0
    printlines = []
    programversion = ''
    procs = ''
    tmpfiles = ''
    configlines = []
    configlinesall = ''
    modswitch = 0
    modlines = []
    pointstoic = ''
    point = ''
    stoic = ''
    numatoms = ''
    chargemulti = ''
    charge = ''
    multi = ''
    basisfunc = ''
    scfenergy = ''
    predechange = ''
    optswitch = 0
    optcomplete = 0
    optstep = 0
    if printtype == 'none':
        energycon = ''
        maxforce = ''
        rmsforce = ''
        maxdisp = ''
        rmsdisp = ''
    else:
        energycon = []
        maxforce = []
        rmsforce = []
        maxdisp = []
        rmsdisp = []
    freqswitch = 0
    if trans == 'list':
        frequencies = []
        zeropointesum = '0.00000'
        thermesum = '0.00000'
        enthalpyesum = '0.00000'
        gibbsesum = '0.00000'
    else:
        frequencies = []
        zeropointesum = ''
        thermesum = ''
        enthalpyesum = ''
        gibbsesum = ''
    freqtemp = ''
    freqpres = ''
    negfreq = 0
    zeropointe = ''
    therme = ''
    enthalpye = ''
    gibbse = ''
    stotal = ''
    entropyesum = ''
    selec = ''
    srot = ''
    strans = ''
    svib = ''
    freehswitch = 0
    jobdone = 0
    riccswitch = 0
    rhfenergy = ''
    correnergy = ''
    mp2energy = ''
    findbasis = 0
    findtheory = 0

    f = open(logfile)
    filelist = f.readlines()
    f.close()
    for (i, var) in enumerate(filelist):

    # if i < 30 and 'Jaguar version' not in var and jaguartext == 0:
    #  print "Error! Not a Jaguar out file."
    # exit(0)
    # Find the Gaussian Version

        if 'TURBOMOLE V' in var:
            programversion = var.split(':')[1].split()[0] + ' ' \
                + var.split(':')[1].split()[1]
        elif 'TURBOTMPDIR environment variable set to' in var:

    # Find the tmp file

            tmpfiles = var.split('"')[1].strip()
        elif 'PROCESSORS!' in var:

    # Find the processors/threads used

            procs = var.split('ON')[1].strip().split()[0]

    # Get pointgroup

        if 'symmetry group' in var:
            point = var.split()[-1]
        elif 'type   atoms' in var:

    # Get stoic and number of atoms

            stoic = ''
            for j in range(2, 100):
                if '---------------------------------------------------------------------------' \
                    not in filelist[i + j]:
                    stoic += filelist[i + j].split()[0].upper() \
                        + filelist[i + j].split()[1]
                else:

                    break
            if 'total' in filelist[i + 7]:
                numatoms = filelist[i + 7].split()[1]
        elif 'charge' in var \
            and '------------------------------------------------------------------------------' \
            in filelist[i + 1]:

    # Get charge and multiplicity

            charge = filelist[i + 2].strip().split()[-1]
        elif 'molecular charge' in var:
            charge = var.rstrip().split()[-1]
        elif 'real representations' in var:
            multi = var.strip().split()[2]
        elif 'number of basis functions' in var:

    # Find basis functions

            basisfunc = var.strip().split(':')[-1].strip()

    # Find the job configuration used

        if 'Non-default options chosen:' in var:
            for j in range(1, 30):
                if 'Temporary integer options' not in filelist[i + j]:
                    configlines.append(filelist[i + j].strip())
                    if 'modredundant' in filelist[i + j]:
                        modswitch = 1
                else:
                    break
        elif 'type   atoms  prim   cont   basis' in var and findbasis \
            == 0:
            configlines.append('basis=' + filelist[i + 2].split()[-2])
            findbasis += 1
        elif 'density functional' in var and '--------' in filelist[i
                + 1] and findtheory == 0:

            configlines.append('theory=' + filelist[i + 2].split()[-1])
            findtheory += 1
        elif 'jobex' in var and optswitch == 0:
            optswitch = 1
            configlines.append('optimization')
        elif 'a o f o r c e' in var or 'N u m F o r c e' in var:
            freqswitch = 1
            configlines.append('freq:aoforce')
        elif 'f r e e   e n t h a l p y' in var or 'freeh' in var:
            freehswitch = 1
            configlines.append('freq:numforce')
        elif modswitch == 1:

    # Check to see if modredundant is in the config

            if 'The following ModRedundant input section has been read:' \
                in var:
                for j in range(1, 20):
                    if len(filelist[i + j]) > 2:
                        modlines.append(filelist[i + j].strip())
                    else:
                        break
        elif 'total energy      =' in var:

    # Get SCF Energy and predicted change

            scfenergy = var.strip().split('=')[1].split()[0]
        elif 'ricc2' in var:

    # Check for ricc2 calculation

            riccswitch = 1
            if 'RHF  energy' in var:
                rhfenergy = var.split(':')[1].split()[0].strip()
            elif 'Final MP2 energy' in var:
                mp2energy = var.split(':')[1].split()[0].strip()
            elif 'correlation energy' in var:
                correnergy = var.split(':')[1].split()[0].strip()
        elif 'Predicted energy change:' in var:
            predechange = var.strip().split(':')[-1].strip()
        elif optswitch == 1:

    # Optimization parse

            if 'CONVERGENCY CRITERIA FULFILLED IN CYCLE' in var:
                optcomplete = 1
            elif 'OPTIMIZATION CYCLE' in var:
                optstep = int(var.split()[-1])
            if printtype == 'none':
                if 'MAX gradient' in var:
                    maxforce = var.strip().split()
                    if maxforce[-3] == 'yes':
                        maxforce[-3] = '\033[1;32mYES\033[0;m'
                    else:
                        maxforce[-3] = '\033[1;31mNO!\033[0;m'
                elif 'RMS of gradient' in var:
                    rmsforce = var.strip().split()
                    if rmsforce[-3] == 'yes':
                        rmsforce[-3] = '\033[1;32mYES\033[0;m'
                    else:
                        rmsforce[-3] = '\033[1;31mNO!\033[0;m'
                elif 'MAX displacement' in var:
                    maxdisp = var.strip().split()
                    if maxdisp[-3] == 'yes':
                        maxdisp[-3] = '\033[1;32mYES\033[0;m'
                    else:
                        maxdisp[-3] = '\033[1;31mNO!\033[0;m'
                elif 'RMS of displacement' in var:
                    rmsdisp = var.strip().split()
                    if rmsdisp[-3] == 'yes':
                        rmsdisp[-3] = '\033[1;32mYES\033[0;m'
                    else:
                        rmsdisp[-3] = '\033[1;31mNO!\033[0;m'
                elif 'Energy change' in var:
                    energycon = var.strip().split()
                    if energycon[-3] == 'yes':
                        energycon[-3] = '\033[1;32mYES\033[0;m'
                    else:
                        energycon[-3] = '\033[1;31mNO!\033[0;m'
            if printtype == 'long' or printtype == 'extra':
                if 'MAX gradient' in var and 'Threshold' not in var:
                    maxforce.append(var.strip().split())
                    if maxforce[-1][-3] == 'yes':
                        maxforce[-1][-3] = '\033[1;32mYES\033[0;m'
                    else:
                        maxforce[-1][-3] = '\033[1;31mNO!\033[0;m'
                elif 'RMS of gradient' in var and 'Threshold' \
                    not in var:
                    rmsforce.append(var.strip().split())
                    if rmsforce[-1][-3] == 'yes':
                        rmsforce[-1][-3] = '\033[1;32mYES\033[0;m'
                    else:
                        rmsforce[-1][-3] = '\033[1;31mNO!\033[0;m'
                elif 'MAX displacement' in var and 'Threshold' \
                    not in var:
                    maxdisp.append(var.strip().split())
                    if maxdisp[-1][-3] == 'yes':
                        maxdisp[-1][-3] = '\033[1;32mYES\033[0;m'
                    else:
                        maxdisp[-1][-3] = '\033[1;31mNO!\033[0;m'
                elif 'RMS of displacement' in var and 'Threshold' \
                    not in var:
                    rmsdisp.append(var.strip().split())
                    if rmsdisp[-1][-3] == 'yes':
                        rmsdisp[-1][-3] = '\033[1;32mYES\033[0;m'
                    else:
                        rmsdisp[-1][-3] = '\033[1;31mNO!\033[0;m'
                elif 'Energy change' in var and 'ratio' not in var:
                    energycon.append(var.strip().split())
                    if energycon[-1][-3] == 'yes':
                        energycon[-1][-3] = '\033[1;32mYES\033[0;m'
                    else:
                        energycon[-1][-3] = '\033[1;31mNO!\033[0;m'

    # Frequency logic

        if freqswitch == 1:
            if 'frequency' in var and not 'zero' in var:
                for j in range(1, len(var.split())):
                    frequencies.append(var.split()[j])
            if 'zero point vibrational energy' in var:
                zeropointe = str(float(filelist[i + 2].split()[1])
                                 / 4.184 / 2625.49962)
            elif 'T        p       ln(qtrans) ln(qrot) ln(qvib) chem.pot.   energy    entropy' \
                in var:
                freqtemp = filelist[i + 3].split()[0]
                freqpres = filelist[i + 3].split()[1]

        # ##therme = str((float(filelist[i + 3].split()[1]) + (298.15 * 8.3144621 / 1000)) / 2625.49962)

                enthalpye = str((float(filelist[i + 7].split()[5])
                                + 298.15 * 8.3144621 / 1000)
                                / 2625.49962)
                gibbse = str(float(filelist[i + 3].split()[4])
                             / 2625.49962)
                stotal = str(float(filelist[i + 3].split()[5]) * 298.15
                             / 2625.49962)
                thermesum = str(float(therme) + float(scfenergy))
                enthalpyesum = str(float(enthalpye) + float(scfenergy))
                gibbsesum = str(float(gibbse) + float(scfenergy))
                entropyesum = str(float(stotal) + float(scfenergy))

    # Job complete logic

        if i - len(filelist) < 10 and 'total wall-time' in var:  # i > len(filelist) - 30 and
            jobdone = 1
        elif i - len(filelist) < 10 and ('next step' in var
                or 'OPTIMIZATION CYCLE' in var or 'RUNNING PROGRAM'
                in var or 'STARTING ' in var):

                                                                                                                                            # i > len(filelist)-30 and

            jobdone = 0

    if trans == 'translate':

    # Printed output list

        printlines = ['Analyzing Turbomole Output File: ' + logfile]
        printlines.append('Using: ' + str(programversion))
        printlines.append('[#Processors=' + procs + '  Temp Files='
                          + tmpfiles + ']')
        printlines.append('=============================================================================='
                          )
    elif trans == 'archive':

    # Initial information

        printlines.append('Supporting Information: ' + logfile)
        printlines.append('------------------------------------------------------------------------------'
                          )
        printlines.append('Using: ' + programversion)
        printlines.append('=============================================================================='
                          )

    if trans == 'translate' or trans == 'archive':

        # Print config

        configlinestmp = ' '.join(str(e) for e in configlines)
        segments = list()
        while len(configlinestmp) > 0:
            segments.append(configlinestmp[:78])
            configlinestmp = configlinestmp[78:]
        for line in segments:
            printlines.append(line)

    # printlines.append(' '.join(str(e) for e in configlines))

    # Modredundant output

        if modswitch == 1:
            for line in modlines:
                printlines.append(line)

    # Print pointgroup and stoichiometry

        printlines.append('------------------------------------------------------------------------------'
                          )
        printlines.append('Pointgroup = ' + point
                          + '      Stoichiometry = ' + stoic
                          + '       # of Atoms = ' + numatoms)

    # Print Charge Multi and BasisFunc

        printlines.append('Charge = ' + charge[0:5]
                          + '       Multiplicity = ' + multi
                          + '                  Basis = ' + basisfunc)
        printlines.append('------------------------------------------------------------------------------'
                          )

    # Print Energy

        printlines.append('SCF Energy = ' + scfenergy
                          + '\tPredicted Change = ' + predechange)
        printlines.append('------------------------------------------------------------------------------'
                          )

    # Print MP2 Energy

        if riccswitch == 1:
            printlines.append('RHF Energy: ' + rhfenergy
                              + ' |  MP2 Energy: ' + mp2energy
                              + ' | Correlation Energy: ' + correnergy)
            printlines.append('------------------------------------------------------------------------------'
                              )

    # Optimization output

        if optswitch == 1:
            if optcomplete == 1 and optstep == 0:
                printlines.append('\033[1;32mOptimization completed.\033[0;m'
                                  )
            elif optcomplete == 1 and optstep >= 0:
                printlines.append('\033[1;32mOptimization completed.\033[0;m  ----------------------------  Optimization step: '
                                   + str(optstep))
            elif optcomplete == 0 and optstep == 0:
                printlines.append('\033[1;31mOptimization incomplete.\033[0;m'
                                  )
            elif optcomplete == 0 and optstep >= 0:
                printlines.append('\033[1;31mOptimization incomplete.\033[0;m  ----------------------------  Optimization step: '
                                   + str(optstep))
            if len(maxforce) > 0 and printtype == 'none':
                printlines.append('Item      Max Val.    Criteria    Pass?      RMS Val.    Criteria    Pass?'
                                  )
                printlines.append('Energy    ' + energycon[-2] + ' || '
                                  + energycon[-1] + '    '
                                  + energycon[-3])
                printlines.append('Force     ' + maxforce[-2] + ' || '
                                  + maxforce[-1] + '    '
                                  + maxforce[-3] + '        '
                                  + rmsforce[-2] + ' || '
                                  + rmsforce[-1] + '    '
                                  + rmsforce[-3])
                printlines.append('Disp      ' + maxdisp[-2] + ' || '
                                  + maxdisp[-1] + '    ' + maxdisp[-3]
                                  + '        ' + rmsdisp[-2] + ' || '
                                  + rmsdisp[-1] + '    ' + rmsdisp[-3])
                printlines.append('------------------------------------------------------------------------------'
                                  )
            elif len(maxforce) > 0 and printtype == 'long':
                printlines.append('Opt Step:   Energy       Force [Max || RMS]          Displacement [Max || RMS]'
                                  )
                for i in range(len(maxforce)):
                    printlines.append(' ' + str(i + 1) + '  '
                            + energycon[i][-2] + ' ('
                            + energycon[i][-3] + ')  ['
                            + maxforce[i][-2] + ' (' + maxforce[i][-3]
                            + ') || ' + rmsforce[i][-2] + ' ('
                            + rmsforce[i][-3] + ')]  ['
                            + maxdisp[i][-2] + ' (' + maxdisp[i][-3]
                            + ') || ' + rmsdisp[i][-2] + ' ('
                            + rmsdisp[i][-3] + ')]')
                printlines.append('------------------------------------------------------------------------------'
                                  )
            elif len(maxforce) > 0 and printtype == 'extra':
                printlines.append('Opt Step:   Energy       Force [Max || RMS]          Displacement [Max || RMS]'
                                  )
                if len(maxforce) > 5:
                    lastlines = len(maxforce) - 5
                else:
                    lastlines = 0
                for i in range(lastlines, len(maxforce)):
                    printlines.append(' ' + str(i + 1) + '  '
                            + energycon[i][-2] + ' ('
                            + energycon[i][-3] + ')  ['
                            + maxforce[i][-2] + ' (' + maxforce[i][-3]
                            + ') || ' + rmsforce[i][-2] + ' ('
                            + rmsforce[i][-3] + ')]  ['
                            + maxdisp[i][-2] + ' (' + maxdisp[i][-3]
                            + ') || ' + rmsdisp[i][-2] + ' ('
                            + rmsdisp[i][-3] + ')]')
                printlines.append('------------------------------------------------------------------------------'
                                  )

    if trans == 'archive':
        printlines.append('------------------------------------------------------------------------------'
                          )
        printlines.append('    Atomic        Coordinates (Angstroms)')
        printlines.append('    Type:            X       Y       Z')
        printlines.append('------------------------------------------------------------------------------'
                          )
        printlines.append('GEOM DATA!')

    if trans == 'translate' or trans == 'archive':

    # Frequency output

        if freqswitch == 1 and len(frequencies) > 0:
            printlines.append('Statistical Thermodynamic Analysis for '
                              + logfile)
        if freqswitch == 1 and freehswitch == 1:
            printlines.append('Computed at:  ' + freqtemp + '  and  '
                              + freqpres)
            printlines.append('=============================================================================='
                              )
            printlines.append('SCF Energy=                                             '
                               + scfenergy)
            printlines.append('Zero-point correction (ZPE)=                         '
                               + zeropointe)
            printlines.append('Internal Energy (U)=                    '
                               + therme + '       ' + thermesum)
            printlines.append('Enthalpy (H)=                           '
                               + enthalpye + '       ' + enthalpyesum)
            printlines.append('Gibbs Free Energy (G)=                  '
                               + gibbse + '       ' + gibbsesum)
            printlines.append('Entropy (S) [Cal/Mol-Kelvin]=                            '
                               + stotal)
            printlines.append('Entropy (S): electronic=                                 '
                               + selec)
            printlines.append('Entropy (S): translational=                              '
                               + strans)
            printlines.append('Entropy (S): rotational=                                 '
                               + srot)
            printlines.append('Entropy (S): vibrational=                                '
                               + svib)
        if freqswitch == 1 and len(frequencies) > 0:
            frequenciesnozero = [x for x in frequencies if x != '0.00']
            printlines.append('Frequencies: ' + frequenciesnozero[0]
                              + '  ' + frequenciesnozero[1] + '  '
                              + frequenciesnozero[2] + '  '
                              + frequenciesnozero[3] + '  '
                              + frequenciesnozero[4] + '  '
                              + frequenciesnozero[5])
            for var in frequencies:
                if 'i' in var:
                    negfreq += 1
            if negfreq == 1:
                printlines.append('\033[1;33m              !!!Warning there were '
                                   + str(negfreq)
                                  + ' negative frequencies!!!\033[0;m')
            elif negfreq >= 2:
                printlines.append('\033[1;31m              !!!Warning there were '
                                   + str(negfreq)
                                  + ' negative frequencies!!!\033[0;m')

    if trans == 'translate':

    # Print end of file for long and extra

        if printtype == 'extra' or printtype == 'long':
            for i in range(len(filelist) - 10, len(filelist)):
                printlines.append(filelist[i].strip())

    # Completed job

        if jobdone == 1:
            printlines.append('=============================================================================='
                              )
            printlines.append('\033[1;32m##################              JOB COMPLETED               ##################\033[0;m'
                              )
            printlines.append('=============================================================================='
                              )

    if trans == 'list':
        if len(frequencies) == 2:
            frequenciesnozero = frequencies
        else:
            frequenciesnozero = [x for x in frequencies if x != '0.00']
        line = scfenergy + '  '
        line += zeropointesum + '  ' + thermesum + '  ' + enthalpyesum \
            + '  ' + gibbsesum + '  '
        frequenciesnozero = [x for x in frequencies if x != '0.00']
        if len(frequenciesnozero) == 0:
            frequenciesnozero = ['0.00000', '0.00000']
        for var in frequencies:
            if 'i' in var:
                negfreq += 1
        if negfreq == 0:
            line += '\033[1;32m' + frequenciesnozero[0] + '  ' \
                + frequenciesnozero[1] + '\033[1;0m' + '  '
        elif negfreq == 1:
            line += '\033[1;33m' + frequenciesnozero[0] + '  ' \
                + frequenciesnozero[1] + '\033[1;0m' + '  '
        elif negfreq >= 2:
            line += '\033[1;31m' + frequenciesnozero[0] + '  ' \
                + frequenciesnozero[1] + '\033[1;0m' + '  '
        line += basisfunc + '  ' + pointstoic + '  ' + stoic + '  '
        if jobdone == 0:
            line += logfile.split('/')[-1]
        elif jobdone == 1:
            filename = logfile.split('/')[-1]
            line += '\033[1;32m' + filename + '\033[0;0m'
        elif jobdone == 2:
            filename = logfile.split('/')[-1]
            line += '\033[1;31m' + filename + '\033[0;0m'
        printlines.append(line)

    if trans != 'archive':
        return printlines
    if trans == 'archive':
        tablelines = []
        tablelines.append(logfile.split('/')[-1] + ',' + scfenergy + ','
                           + zeropointe + ',=' + therme + '+RC[-2],='
                          + enthalpye + '+RC[-3],=' + gibbse
                          + '+RC[-4],=RC[-1]-MIN(C[-1]),=RC[-1]*627.5095'
                          )
        tablelines.append(logfile.split('/')[-1] + ',' + scfenergy + ','
                           + zeropointe + ',' + thermesum + ','
                          + enthalpyesum + ',' + gibbsesum
                          + ',=RC[-1]-MIN(C[-1]),=RC[-1]*627.5095')
        tablelines.append(logfile.split('/')[-1] + ',' + scfenergy)
        return (printlines, tablelines)


### --- Terachem Specific Functions --- ###
###==================================================================================================================###
### --- Parse the Terachem optim file for geometry information --- ###

def parse_optim_terachem(optimfile):
    f = open(optimfile)
    geomlist = f.readlines()
    f.close()
    numatoms = int(geomlist[0])
    allenergy = []
    for line in geomlist:
        if 'frame' in line:
            allenergy.append(float(line.split()[0]))
    finalgeom = [0] * (numatoms + 2)
    finalgeom[0] = str(numatoms)
    for i in range(len(geomlist) - numatoms - 1, len(geomlist)):
        j = i - len(geomlist) + numatoms + 2
        if i == len(geomlist) - numatoms - 1:
            line = geomlist[i].split()
            finalgeom[1] = 'Energy: ' + str(allenergy[-1])
        if i > len(geomlist) - numatoms - 1:
            line = geomlist[i].split()
            finalgeom[j] = Atom()
            finalgeom[j].e = line[0]
            finalgeom[j].en = get_element_num(line[0])
            finalgeom[j].x = float(line[1])
            finalgeom[j].y = float(line[2])
            finalgeom[j].z = float(line[3])
    return finalgeom


def parse_optim_terachem_traj(optimfile):
    f = open(optimfile)
    geomlist = f.readlines()
    f.close()
    numatoms = int(geomlist[0])
    allenergy = []
    for line in geomlist:
        if 'frame' in line:
            allenergy.append(float(line.split()[0]))
    trajectorygeom = [0] * len(geomlist)
    skip = 0
    for i in range(len(geomlist)):
        if i % (numatoms + 2) == 0 and i / (numatoms + 2) \
            != len(allenergy):
            trajectorygeom[i] = geomlist[i].strip()
        elif 'frame' in geomlist[i]:
            trajectorygeom[i] = 'Energy: ' \
                + str(float(geomlist[i].split()[0]))
        else:
            line = geomlist[i].split()
            trajectorygeom[i] = Atom()
            trajectorygeom[i].e = line[0]
            trajectorygeom[i].en = get_element_num(line[0])
            trajectorygeom[i].x = float(line[1]) / auToAngstrom
            trajectorygeom[i].y = float(line[2]) / auToAngstrom
            trajectorygeom[i].z = float(line[3]) / auToAngstrom
    return trajectorygeom


### --- Out put the G0X input files --- ###

def output_tc(
    ifile,
    config,
    charge,
    multi,
    writemod,
    modred,
    writeecp,
    ecplines,
    writeclose,
    closelines,
    ifilelol,
    ):

    ofiletmp = ifile.split('.')
    ofile = ''
    for i in range(len(ofiletmp)):
        if i != len(ofiletmp) - 1:
            ofile += ofiletmp[i] + '.'
        elif i == len(ofiletmp) - 1:
            ofile += 'tc'
    os.system('cp ' + ifile + ' NEWJOBS.TC/')
    f = open('NEWJOBS.TC/' + ofile, 'w+')
    for (i, var) in enumerate(config):
        if 'BASE!' in var:
            linetmp = var.replace('BASE!', ifile)
            f.write(linetmp + '\n')
        elif 'SCRDIR!' in var:
            linetmp = var.replace('SCRDIR!', './scr-' + ifile.lower())
            f.write(linetmp + '\n')
        elif 'CHARGE!' in var:
            linetmp = var.replace('CHARGE!', str(charge))
            f.write(linetmp + '\n')
        elif 'MULTI!' in var:
            linetmp = var.replace('MULTI!', str(multi))
            f.write(linetmp + '\n')
        else:
            f.write(var + '\n')
    if writemod == 1:
        f.write('$constraints\n')
        if len(modred) > 0 and modred[0] != '':
            for (i, var) in enumerate(modred):
                f.write(var + '\n')
        else:
            f.write('MODRED!\n')
        f.write('$end\n')
    if writeclose == 1:
        for (i, var) in enumerate(closelines):
            f.write(var + '\n')
        f.write('\n')
    f.close()
    return


# Parse the Gaussian output file and print data to the screen

def translate_tc_output(logfile, trans, printtype):
    initlines = \
        ['\033[1;34mPython Terachem Output File Parsing Script:   nTCT.py\033[1;0m'
         ,
         '\033[1;34mScripted by Thomas J. L. Mustard, O. Maduka Ogba, Paul Ha-Yeon Cheong\033[1;0m'
         ,
         '------------------------------------------------------------------------------'
         ]
    if trans == 'translate':
        for line in initlines:
            print(line)

  # Global Variables

    config = ''
    configswitch = 0
    printlines = []
    programversion = ''
    chk = ''
    nproc = ''
    mem = ''
    configlines = []
    configlinesall = ''
    modswitch = 0
    modlines = []
    pointstoic = ''
    point = ''
    stoic = ''
    chargemulti = ''
    charge = ''
    multi = ''
    basisfunc = ''
    scfenergy = ''
    predechange = ''
    optswitch = 0
    optcomplete = 0
    optstep = 0
    if printtype == 'long' or printtype == 'extra':
        maxforce = []
        rmsdisp = []
        rmsforce = []
        maxdisp = []
        energycon = []
    else:
        maxforce = ''
        rmsforce = ''
        maxdisp = ''
        rmsdisp = ''
        energycon = ''
    freqswitch = 0
    freqtemp = ''
    if trans == 'list':
        frequencies = [
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            ]
        negfreq = 0
        zeropointe = '0.000000'
        therme = '0.000000'
        enthalpye = '0.000000'
        gibbse = '0.000000'
        zeropointesum = '0.000000'
        thermesum = '0.000000'
        enthalpyesum = '0.000000'
        gibbsesum = '0.000000'
        stotal = '0.000000'
        selec = '0.000000'
        srot = '0.000000'
        strans = '0.000000'
        svib = '0.000000'
    else:
        frequencies = [
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            ]
        negfreq = 0
        zeropointe = 0.000000
        therme = 0.000000
        enthalpye = 0.000000
        gibbse = 0.000000
        zeropointesum = 0.000000
        thermesum = 0.000000
        enthalpyesum = 0.000000
        gibbsesum = 0.000000
        stotal = 0.000000
        selec = 0.000000
        srot = 0.000000
        strans = 0.000000
        svib = 0.000000
    jobdone = 0
    memswitch = 0
    f = open(logfile)
    filelist = f.readlines()
    f.close()
    for (i, var) in enumerate(filelist):

    # Find the Terachem Version

        if 'TeraChem v' in var:
            programversion = 'TeraChem ' + var.split()[2]
        elif 'Scratch directory:' in var:

    # Find the chk file

            chk = var.split(':')[-1].strip()
        elif 'CUDA devices' in var:

    # Find the nproc used

            nproc = var.split()[1]
        elif 'GPU Memory Available:' in var:

    # Find the mem used

            mem = var.split(':')[1].strip()
        elif 'Method:' in var:

    # Find the job configuration used

            configlines.append('method=' + var.split()[-1])
        elif 'Using dynamic DFT grids.' in var:
            configlines.append('dynamicgrid=yes')
        elif 'Basis set:' in var:
            configlines.append('basis=' + var.split()[-1])
        elif 'Maximum number of SCF iterations:' in var:
            configlines.append('maxit=' + var.split()[-1])
        elif 'Max number of optimization cycles:' in var:
            configlines.append('nstep=' + var.split()[-1])
        elif 'Optimization will be performed' in var:
            configlines.append('min_coordinates=' + var.split()[-2])
        elif 'RUNNING GEOMETRY OPTMIZATION' in var:
            optswitch = 1
        elif 'modredundant' in filelist[i]:
            modswitch = 1
        elif 'RUNNING THERMOCHEMISTRY ANALYSIS' in filelist[i]:
            freqswitch = 1
        elif modswitch == 1:

    # Check to see if modredundant is in the config

            if 'The following ModRedundant input section has been read:' \
                in var:
                for j in range(1, 20):
                    if len(filelist[i + j]) > 2:
                        modlines.append(filelist[i + j].strip())
                    else:
                        break

    # Get pointgroup and stoic

        if 'Framework group' in var:
            pointstoic = var.split()[-1]
            point = pointstoic.split('[')[0]
        elif 'Stoichiometry' in var:
            stoic = var.split()[-1]
        elif 'Spin multiplicity:' in var:

    # Get charge and multiplicity

            multi = var.split()[-1]
        elif 'Total charge:' in var:
            charge = var.split()[-1]
        elif 'Total orbitals:' in var:

    # Find basis functions

            basisfunc = filelist[i].split()[-1]
        elif 'FINAL ENERGY:' in var:

    # Get SCF Energy and predicted change

            scfenergy = var.split()[-2]
        elif 'Predicted change in Energy=' in var:
            predechange = var.strip().split('=')[-1]
        elif optswitch == 1:

    # Optimization parse

            if 'Converged!' in var and 'converged' in filelist[i + 1]:
                optcomplete = 1
            if printtype == 'none':
                if 'Testing convergence  in cycle' in var:
                    optstep = var.split()[-1]

        # elif 'Energy change' in var and 'ratio' not in var:

                    energycon = filelist[i + 1].strip().split()
                    if energycon[-1] == 'no':
                        energycon[-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif energycon[-1] == 'yes':
                        energycon[-1] = '\033[1;32mYES\033[1;0m'
                elif 'Max grad' in var:
                    maxforce = var.rstrip().split()
                    if maxforce[-1] == 'no':
                        maxforce[-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif maxforce[-1] == 'yes':
                        maxforce[-1] = '\033[1;32mYES\033[1;0m'
                elif 'RMS grad' in var:
                    rmsforce = var.rstrip().split()
                    if rmsforce[-1] == 'no':
                        rmsforce[-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif rmsforce[-1] == 'yes':
                        rmsforce[-1] = '\033[1;32mYES\033[1;0m'
                elif 'Max step' in var:
                    maxdisp = var.rstrip().split()
                    if maxdisp[-1] == 'no':
                        maxdisp[-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif maxdisp[-1] == 'yes':
                        maxdisp[-1] = '\033[1;32mYES\033[1;0m'
                elif 'RMS step' in var:
                    rmsdisp = var.rstrip().split()
                    if rmsdisp[-1] == 'no':
                        rmsdisp[-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif rmsdisp[-1] == 'yes':
                        rmsdisp[-1] = '\033[1;32mYES\033[1;0m'
            elif printtype == 'long' or printtype == 'extra':
                if 'Testing convergence  in cycle' in var:
                    optstep = var.split()[2].strip()
                    if 'Energy' not in filelist[i + 1]:
                        energycon.append([
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '\x1b[1;31mNO!\x1b[1;0m',
                            ])

            # energycon.append(["0.0000", "0.0000", "0.0000", "0.0000", "0.0000", "0.0000", "0.0000", "\033[1;31mNO!\033[1;0m"])

                        maxdisp.append([
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '\x1b[1;31mNO!\x1b[1;0m',
                            ])
                        rmsdisp.append([
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '\x1b[1;31mNO!\x1b[1;0m',
                            ])
                    else:
                        energycon.append(list(filelist[i
                                + 1].rstrip().split()))
                    if energycon[-1][-1] == 'no':
                        energycon[-1][-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif energycon[-1][-1] == 'yes':
                        energycon[-1][-1] = '\033[1;32mYES\033[1;0m'
                elif 'Max grad' in var and 'component' not in var:
                    maxforce.append(list(var.rstrip().split()))
                    if maxforce[-1][-1] == 'no':
                        maxforce[-1][-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif maxforce[-1][-1] == 'yes':
                        maxforce[-1][-1] = '\033[1;32mYES\033[1;0m'
                elif 'RMS grad' in var:
                    rmsforce.append(list(var.rstrip().split()))
                    if rmsforce[-1][-1] == 'no':
                        rmsforce[-1][-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif rmsforce[-1][-1] == 'yes':
                        rmsforce[-1][-1] = '\033[1;32mYES\033[1;0m'
                elif 'Max step' in var:
                    maxdisp.append(list(var.rstrip().split()))
                    if maxdisp[-1][-1] == 'no':
                        maxdisp[-1][-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif maxdisp[-1][-1] == 'yes':
                        maxdisp[-1][-1] = '\033[1;32mYES\033[1;0m'
                elif 'RMS step' in var:
                    rmsdisp.append(list(var.rstrip().split()))
                    if rmsdisp[-1][-1] == 'no':
                        rmsdisp[-1][-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif rmsdisp[-1][-1] == 'yes':
                        rmsdisp[-1][-1] = '\033[1;32mYES\033[1;0m'

    # Frequency and

        if freqswitch == 1:
            if 'Temperature' in var:
                freqtemp = var.strip()
            elif 'Thermochemical analysis' in var:

      # elif 'Vibrational temperatures:' in var:
      #  frequencies = var.strip().split(':')[-1].strip()

                if frequencies[0] == '0.00':
                    frequencies = []
                frequencies.append(str(filelist[i
                                   + 3].split()[2].strip()))
                frequencies.append(str(filelist[i
                                   + 4].split()[2].strip()))
                frequencies.append(str(filelist[i
                                   + 5].split()[2].strip()))
                frequencies.append(str(filelist[i
                                   + 6].split()[2].strip()))
                frequencies.append(str(filelist[i
                                   + 7].split()[2].strip()))
                for j in frequencies:
                    line = float(j)
                    if line < 0:
                        negfreq += 1
            elif 'total ZPE' in var:
                zeropointe = float(var.strip().split()[-2])
            elif 'total E vib' in var:
                enthalpye = float(var.strip().split()[-2])
            elif 'total S vib' in var:

        # therme = str(var.strip().split()[-2])

                svib = float(var.strip().split()[-2])

      #  enthalpye = str(var.strip().split()[-1])
      # elif 'Thermal correction to Gibbs Free Energy' in var:
      #  gibbse = str(var.strip().split()[-1])
      # elif 'Sum of electronic and zero-point Energies' in var:
      #  zeropointesum = str(var.strip().split()[-1])
      # elif 'Sum of electronic and thermal Energies' in var:
      #  thermesum = str(var.strip().split()[-1])
      # elif 'Sum of electronic and thermal Enthalpies' in var:
      #  enthalpyesum = str(var.strip().split()[-1])
      # elif 'Sum of electronic and thermal Free Energies' in var:
      #  gibbsesum = str(var.strip().split()[-1])
      #  stotal = str(filelist[i+4].strip().split()[-1])
      #  selec = str(filelist[i+5].strip().split()[-1])
      #  strans = str(filelist[i+6].strip().split()[-1])
      #  srot = str(filelist[i+7].strip().split()[-1])
      #  svib = str(filelist[i+8].strip().split()[-1])

        if 'Job finished:' in var and len(filelist) - i < 10:
            jobdone = 1
        elif 'Job finished:' in var and i - len(filelist) < 10:
            jobdone = 2

    if trans == 'translate':

    # Printed output list
    # printlines = ['\033[1;44m                                                                              \033[1;0m']

        printlines = ['Analyzing Terachem Output File: ' + logfile]
        printlines.append('Using ' + str(programversion))
        printlines.append('[#GPUs=' + nproc + '  Memory=' + mem
                          + '  Scratch Dir=' + chk + ']')
        printlines.append('=============================================================================='
                          )
    elif trans == 'archive':

    # Initial information

        printlines.append('Supporting Information: ' + logfile)
        printlines.append('------------------------------------------------------------------------------'
                          )
        printlines.append('Using: ' + programversion)
        printlines.append('=============================================================================='
                          )

    if trans == 'archive' or trans == 'translate':

    # Print config

        configlinestmp = ' '.join(str(e) for e in configlines)
        segments = list()
        while len(configlinestmp) > 0:
            segments.append(configlinestmp[:78])
            configlinestmp = configlinestmp[78:]
        for line in segments:
            printlines.append(line)

    # printlines.append(' '.join(str(e) for e in configlines))

    # Modredundant output

        if modswitch == 1:
            for line in modlines:
                printlines.append(line)

    # Print pointgroup and stoichiometry

        printlines.append('------------------------------------------------------------------------------'
                          )

    # printlines.append('Pointgroup=' + point + '\tStoichiometry=' + stoic + '\t' + pointstoic)

    # Print Charge Multi and BasisFunc

        printlines.append('Charge = ' + charge + '\tMultiplicity = '
                          + multi + '\tBasis = ' + basisfunc)
        printlines.append('------------------------------------------------------------------------------'
                          )

    # Print Energy

        if freqswitch == 0:
            printlines.append('SCF Energy = ' + scfenergy)  # + '\tPredicted Change = ' + predechange)
            printlines.append('------------------------------------------------------------------------------'
                              )

    # Optimization output

        if optswitch == 1:
            if optcomplete == 1 and optstep != 0:
                printlines.append('\033[1;32mOptimization completed.\033[0;0m  ----------------------------  Optimization step: '
                                   + optstep)
            elif optcomplete == 1 and optstep == 0:
                printlines.append('\033[1;32mOptimization completed.\033[0;0m'
                                  )
            elif optcomplete == 0 and optstep != 0:
                printlines.append('\033[1;31mOptimization incomplete.\033[0;0m  ----------------------------  Optimization step: '
                                   + optstep)
            elif optcomplete == 0 and optstep == 0:
                printlines.append('\033[1;31mOptimization incomplete.\033[0;0m'
                                  )
            if len(maxforce) > 0 and printtype == 'none':
                printlines.append('Item      Max Val.    Criteria    Pass?      RMS Val.    Criteria    Pass?'
                                  )
                printlines.append('Energy    ' + energycon[1] + ' || '
                                  + energycon[3] + '    '
                                  + energycon[-1])
                printlines.append('Force     ' + maxforce[2] + ' || '
                                  + maxforce[4] + '    ' + maxforce[-1]
                                  + '      ' + rmsforce[2] + ' || '
                                  + rmsforce[4] + '    ' + rmsforce[-1])
                printlines.append('Disp      ' + maxdisp[2] + ' || '
                                  + maxdisp[4] + '    ' + maxdisp[-1]
                                  + '      ' + rmsdisp[2] + ' || '
                                  + rmsdisp[4] + '    ' + rmsdisp[-1])
                printlines.append('------------------------------------------------------------------------------'
                                  )
            elif len(maxforce) > 0 and printtype == 'long':
                printlines.append('Opt Step:   Energy       Force [Max || RMS]          Displacement [Max || RMS]'
                                  )
                for i in range(len(maxforce)):
                    printlines.append(' ' + str(i + 1) + ' '
                            + energycon[i][1] + '(' + energycon[i][-1]
                            + ') [' + maxforce[i][2] + '('
                            + maxforce[i][-1] + ')||' + rmsforce[i][2]
                            + '(' + rmsforce[i][-1] + ')] ['
                            + maxdisp[i][2] + '(' + maxdisp[i][-1]
                            + ')||' + rmsdisp[i][2] + '('
                            + rmsdisp[i][-1] + ')]')
                printlines.append('------------------------------------------------------------------------------'
                                  )
            elif len(maxforce) > 0 and printtype == 'extra':
                printlines.append('Opt Step:   Energy       Force [Max || RMS]          Displacement [Max || RMS]'
                                  )
                if len(maxforce) > 5:
                    lastlines = len(maxforce) - 5
                else:
                    lastlines = 0
                for i in range(lastlines, len(maxforce)):
                    printlines.append(' ' + str(i + 1) + ' '
                            + energycon[i][1] + '(' + energycon[i][-1]
                            + ') [' + maxforce[i][2] + '('
                            + maxforce[i][-1] + ')||' + rmsforce[i][2]
                            + '(' + rmsforce[i][-1] + ')] ['
                            + maxdisp[i][2] + '(' + maxdisp[i][-1]
                            + ')||' + rmsdisp[i][2] + '('
                            + rmsdisp[i][-1] + ')]')
                printlines.append('------------------------------------------------------------------------------'
                                  )

    if trans == 'archive':

        printlines.append('    Atomic        Coordinates (Angstroms)')
        printlines.append('    Type:            X       Y       Z')
        printlines.append('------------------------------------------------------------------------------'
                          )
        printlines.append('GEOM DATA!')

    if trans == 'translate' or trans == 'archive':

    # Frequency output

        if freqswitch == 1 and enthalpye != '':
            printlines.append('Statistical Thermodynamic Analysis for '
                              + logfile)
            printlines.append(freqtemp)
            printlines.append('=============================================================================='
                              )

      # printlines.append('SCF Energy=                                             ' + scfenergy)

            printlines.append('Zero-point correction (ZPE)=            '
                               + str(zeropointe) + '  J/mol')  # + zeropointesum)

      # printlines.append('Internal Energy (U)=                    ' + therme + '  J/mol' + thermesum)

            printlines.append('Enthalpy (H)=                           '
                               + str(enthalpye) + '        ')  # + enthalpyesum)

      # printlines.append('Gibbs Free Energy (G)=                  ' + gibbse + '        ' + gibbsesum)
      # printlines.append('Entropy (S) [Cal/Mol-Kelvin]=           ' + stotal)
      # printlines.append('Entropy (S): electronic=                ' + selec)
      # printlines.append('Entropy (S): translational=             ' + strans)
      # printlines.append('Entropy (S): rotational=                ' + srot)

            printlines.append('Entropy (S): vibrational=               '
                               + str(svib) + '  J/mol/K')
            printlines.append('Frequencies: ' + frequencies[0] + '  '
                              + frequencies[1] + '  ' + frequencies[2]
                              + '  ' + frequencies[3] + '  '
                              + frequencies[4])
            if negfreq == 1:
                printlines.append('\033[1;33m              !!!Warning there were '
                                   + str(negfreq)
                                  + ' negative frequencies!!!\033[0;m')
            elif negfreq >= 2:
                printlines.append('\033[1;31m              !!!Warning there were '
                                   + str(negfreq)
                                  + ' negative frequencies!!!\033[0;m')
            printlines.append('------------------------------------------------------------------------------'
                              )

    if trans == 'translate':

    # Print end of file for long and extra

        if printtype == 'extra' or printtype == 'long':
            for i in range(len(filelist) - 10, len(filelist)):
                printlines.append(filelist[i].strip())

    # Completed job

        if jobdone == 1:
            printlines.append('=============================================================================='
                              )
            printlines.append('\033[1;32m##################              JOB COMPLETED               ##################\033[0;0m'
                              )
            printlines.append('=============================================================================='
                              )
        elif jobdone == 2:
            printlines.append('=============================================================================='
                              )
            printlines.append('\033[1;31m!!!!!!!!!!!!!!!!!!            !!!JOB FAILED!!!              !!!!!!!!!!!!!!!!!!\033[0;0m'
                              )
            printlines.append('=============================================================================='
                              )
    elif trans == 'list':

        line = scfenergy + '  '
        line += zeropointesum + '  ' + thermesum + '  ' + enthalpyesum \
            + '  ' + gibbsesum + '  '
        if negfreq == 0:
            line += '\033[1;32m' + frequencies[0] + '  ' \
                + frequencies[1] + '\033[1;0m' + '  '
        elif negfreq == 1:
            line += '\033[1;33m' + frequencies[0] + '  ' \
                + frequencies[1] + '\033[1;0m' + '  '
        elif negfreq >= 2:
            line += '\033[1;31m' + frequencies[0] + '  ' \
                + frequencies[1] + '\033[1;0m' + '  '
        line += basisfunc + '  ' + pointstoic + '  ' + stoic + '  '
        if jobdone == 0:
            line += logfile.split('/')[-1]
        elif jobdone == 1:
            filename = logfile.split('/')[-1]
            line += '\033[1;32m' + filename + '\033[0;0m'
        elif jobdone == 2:
            filename = logfile.split('/')[-1]
            line += '\033[1;31m' + filename + '\033[0;0m'
        printlines.append(line)

    if trans != 'archive':
        return printlines
    if trans == 'archive':
        entropyesum = str(float(svib) / 1000 / 2625.5)
        enthalpyesum = str(float(enthalpye) / 1000 / 2625.5)
        entropyehart = float(svib) / 1000 / 2625.5
        enthalpyehart = enthalpye / 1000 / 2625.5
        if len(freqtemp) > 0:
            gibbse = float(svib) / 1000 / 2625.5 \
                * float(freqtemp.split()[1]) + enthalpye / 1000 / 2625.5
            gibbsesum = gibbse
        else:
            gibbse = 0
            gibbsesum = 0
        tablelines = []
        tablelines.append(logfile.split('/')[-1] + ',' + scfenergy + ','
                           + str(zeropointe) + ',=' + str(entropyehart)
                          + '+RC[-2],=' + str(enthalpyehart)
                          + '+RC[-3],=' + str(gibbse)
                          + '+RC[-4],=RC[-1]-MIN(C[-1]),=RC[-1]*627.5095'
                          )
        tablelines.append(logfile.split('/')[-1] + ',' + scfenergy + ','
                           + str(zeropointe) + ',' + str(thermesum)
                          + ',' + str(enthalpyesum) + ','
                          + str(gibbsesum)
                          + ',=RC[-1]-MIN(C[-1]),=RC[-1]*627.5095')
        tablelines.append(logfile.split('/')[-1] + ',' + scfenergy)
        return (printlines, tablelines)


### --- GULP Specific Functions --- ###
###==================================================================================================================###
### --- Parse the GULP out file for geometry information --- ###

def parse_output_gulp(ifile):
    alist = []
    geom = 0
    numatoms = 0
    the_file = open(ifile, 'r')
    for (idx, line) in enumerate(the_file):
        if line.strip() != '':
            if 'number of k points' in line and geom >= 5:
                geom = 0
            if geom >= 1:
                geom += 1
            if line.strip() == 'Cartesian axes':
                geom = 1
            if geom >= 3:
                alist.append(line)
    print(alist)
    the_file.close()
    alistdedup = [x for x in alist if 'number of k points' not in x]
    numatoms = int(alistdedup[-1].split()[0])
    print(numatoms)
    finalgeom = [0] * (numatoms + 2)
    finalgeom[0] = str(numatoms)
    finalgeom[1] = 'Energy: '

  # print alistdedup

    for i in range(len(alistdedup) - numatoms, len(alistdedup)):
        line = alistdedup[i].split()
        j = i - len(alistdedup) + numatoms + 2
        finalgeom[j] = Atom()
        finalgeom[j].e = get_element_name(get_element_num(line[1]))
        finalgeom[j].en = get_element_num(line[1])
        finalgeom[j].x = float(line[6])
        finalgeom[j].y = float(line[7])
        finalgeom[j].z = float(line[8])
    return finalgeom


# Parse the GULP output file and print data to the screen

def translate_gulp_output(logfile, trans, printtype):

  # initlines = ['Python GULP Output File Parsing Script:   nGulpT.py\033[1;0m',
  #           'Scripted by Thomas J. L. Mustard, O. Maduka Ogba, Paul Ha-Yeon Cheong\033[1;0m',
  #           '-------------------------------------------------------------------------------------------------------']
  # if trans == "translate":
  #  for line in initlines:
  #    print line

  # Global Variables

    config = ''
    configswitch = 0
    printlines = []
    programversion = ''
    chk = ''
    nproc = ''
    mem = ''
    configlines = []
    configlinesall = ''
    modswitch = 0
    modlines = []
    pointstoic = ''
    point = ''
    stoic = ''
    chargemulti = ''
    charge = ''
    electrons = ''
    multi = ''
    basisfunc = ''
    scfenergy = ''
    predechange = ''
    optswitch = 0
    optcomplete = 0
    optstep = 0
    if printtype == 'long' or printtype == 'extra':
        maxforce = []
        rmsdisp = []
        rmsforce = []
        maxdisp = []
        energycon = []
    else:
        maxforce = ''
        rmsforce = ''
        maxdisp = ''
        rmsdisp = ''
        energycon = ''
    freqswitch = 0
    freqtemp = ''
    if trans == 'list':
        frequencies = [
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            ]
        negfreq = 0
        zeropointe = '0.000000'
        therme = '0.000000'
        enthalpye = '0.000000'
        gibbse = '0.000000'
        zeropointesum = '0.000000'
        thermesum = '0.000000'
        enthalpyesum = '0.000000'
        gibbsesum = '0.000000'
        stotal = '0.000000'
        selec = '0.000000'
        srot = '0.000000'
        strans = '0.000000'
        svib = '0.000000'
    else:
        frequencies = [
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            '0.00',
            ]
        negfreq = 0
        zeropointe = 0.000000
        therme = 0.000000
        enthalpye = 0.000000
        gibbse = 0.000000
        zeropointesum = 0.000000
        thermesum = 0.000000
        enthalpyesum = 0.000000
        gibbsesum = 0.000000
        stotal = 0.000000
        selec = 0.000000
        srot = 0.000000
        strans = 0.000000
        svib = 0.000000
    jobdone = 0
    memswitch = 0
    f = open(logfile)
    filelist = f.readlines()
    f.close()
    for (i, var) in enumerate(filelist):

    # Find the Terachem Version

        if 'Program ' in var:
            programversion = var.split()[1] + ' ' + var.split()[2]
        elif 'Scratch directory:' in var:

    # Find the chk file

            chk = var.split(':')[-1].strip()
        elif 'processors' in var:

    # Find the nproc used

            nproc = var.split()[-2]
        elif 'GPU Memory Available:' in var:

    # Find the mem used

            mem = var.split(':')[1].strip()
        elif 'Exchange-correlation' in var:

    # Find the job configuration used

            configlines.append('method=' + var.split('=')[1].split()[0]
                               + '\n')
        elif 'PseudoPot.' in var:
            configlines.append('basis=' + filelist[i + 1].strip() + '\n'
                               )
        elif 'simplified LDA+U calculation' in var:
            configlines.append('DFTU:\n')
            for j in range(1, 5):
                if filelist[i + j].split()[0] != '-----':
                    configlines.append(filelist[i + j].strip() + '\n')
                elif filelist[i + j].split()[0] == '-----':
                    break
        elif 'RUNNING GEOMETRY OPTMIZATION' in var:

            optswitch = 1
        elif 'modredundant' in filelist[i]:
            modswitch = 1
        elif 'RUNNING THERMOCHEMISTRY ANALYSIS' in filelist[i]:
            freqswitch = 1
        elif modswitch == 1:

    # Check to see if modredundant is in the config

            if 'The following ModRedundant input section has been read:' \
                in var:
                for j in range(1, 20):
                    if len(filelist[i + j]) > 2:
                        modlines.append(filelist[i + j].strip())
                    else:
                        break

    # Get pointgroup and stoic

        if 'Framework group' in var:
            pointstoic = var.split()[-1]
            point = pointstoic.split('[')[0]
        elif 'Stoichiometry' in var:
            stoic = var.split()[-1]
        elif 'Spin multiplicity:' in var:

    # Get charge and multiplicity

            multi = var.split()[-1]
        elif 'Reading input from' in var:
            charge = var.split('_')[-1].split('.qe')[0]
            if charge == '0':
                charge = '  0.0000'
            elif float(charge) > 0:
                charge = ' ' + charge
        elif 'starting charge' in var:
            if var.split()[-1].endswith('.00000'):
                electrons = var.split()[-1]
        elif 'Total orbitals:' in var:

    # Find basis functions

            basisfunc = filelist[i].split()[-1]
        elif 'total energy              =' in var:

    # Get SCF Energy and predicted change

            scfenergy = var.split('=')[-1].strip()
        elif 'Predicted change in Energy=' in var:
            predechange = var.strip().split('=')[-1]
        elif optswitch == 1:

    # Optimization parse

            if 'Converged!' in var and 'converged' in filelist[i + 1]:
                optcomplete = 1
            if printtype == 'none':
                if 'Testing convergence  in cycle' in var:
                    optstep = var.split()[-1]

        # elif 'Energy change' in var and 'ratio' not in var:

                    energycon = filelist[i + 1].strip().split()
                    if energycon[-1] == 'no':
                        energycon[-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif energycon[-1] == 'yes':
                        energycon[-1] = '\033[1;32mYES\033[1;0m'
                elif 'Max grad' in var:
                    maxforce = var.rstrip().split()
                    if maxforce[-1] == 'no':
                        maxforce[-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif maxforce[-1] == 'yes':
                        maxforce[-1] = '\033[1;32mYES\033[1;0m'
                elif 'RMS grad' in var:
                    rmsforce = var.rstrip().split()
                    if rmsforce[-1] == 'no':
                        rmsforce[-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif rmsforce[-1] == 'yes':
                        rmsforce[-1] = '\033[1;32mYES\033[1;0m'
                elif 'Max step' in var:
                    maxdisp = var.rstrip().split()
                    if maxdisp[-1] == 'no':
                        maxdisp[-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif maxdisp[-1] == 'yes':
                        maxdisp[-1] = '\033[1;32mYES\033[1;0m'
                elif 'RMS step' in var:
                    rmsdisp = var.rstrip().split()
                    if rmsdisp[-1] == 'no':
                        rmsdisp[-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif rmsdisp[-1] == 'yes':
                        rmsdisp[-1] = '\033[1;32mYES\033[1;0m'
            elif printtype == 'long' or printtype == 'extra':
                if 'Testing convergence  in cycle' in var:
                    optstep = var.split()[2].strip()
                    if 'Energy' not in filelist[i + 1]:
                        energycon.append([
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '\x1b[1;31mNO!\x1b[1;0m',
                            ])

            # energycon.append(["0.0000", "0.0000", "0.0000", "0.0000", "0.0000", "0.0000", "0.0000", "\033[1;31mNO!\033[1;0m"])

                        maxdisp.append([
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '\x1b[1;31mNO!\x1b[1;0m',
                            ])
                        rmsdisp.append([
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '0.0000',
                            '\x1b[1;31mNO!\x1b[1;0m',
                            ])
                    else:
                        energycon.append(list(filelist[i
                                + 1].rstrip().split()))
                    if energycon[-1][-1] == 'no':
                        energycon[-1][-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif energycon[-1][-1] == 'yes':
                        energycon[-1][-1] = '\033[1;32mYES\033[1;0m'
                elif 'Max grad' in var and 'component' not in var:
                    maxforce.append(list(var.rstrip().split()))
                    if maxforce[-1][-1] == 'no':
                        maxforce[-1][-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif maxforce[-1][-1] == 'yes':
                        maxforce[-1][-1] = '\033[1;32mYES\033[1;0m'
                elif 'RMS grad' in var:
                    rmsforce.append(list(var.rstrip().split()))
                    if rmsforce[-1][-1] == 'no':
                        rmsforce[-1][-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif rmsforce[-1][-1] == 'yes':
                        rmsforce[-1][-1] = '\033[1;32mYES\033[1;0m'
                elif 'Max step' in var:
                    maxdisp.append(list(var.rstrip().split()))
                    if maxdisp[-1][-1] == 'no':
                        maxdisp[-1][-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif maxdisp[-1][-1] == 'yes':
                        maxdisp[-1][-1] = '\033[1;32mYES\033[1;0m'
                elif 'RMS step' in var:
                    rmsdisp.append(list(var.rstrip().split()))
                    if rmsdisp[-1][-1] == 'no':
                        rmsdisp[-1][-1] = '\x1b[1;31mNO!\x1b[1;0m'
                    elif rmsdisp[-1][-1] == 'yes':
                        rmsdisp[-1][-1] = '\033[1;32mYES\033[1;0m'

    # Frequency and

        if freqswitch == 1:
            if 'Temperature' in var:
                freqtemp = var.strip()
            elif 'Thermochemical analysis' in var:

      # elif 'Vibrational temperatures:' in var:
      #  frequencies = var.strip().split(':')[-1].strip()

                if frequencies[0] == '0.00':
                    frequencies = []
                frequencies.append(str(filelist[i
                                   + 3].split()[2].strip()))
                frequencies.append(str(filelist[i
                                   + 4].split()[2].strip()))
                frequencies.append(str(filelist[i
                                   + 5].split()[2].strip()))
                frequencies.append(str(filelist[i
                                   + 6].split()[2].strip()))
                frequencies.append(str(filelist[i
                                   + 7].split()[2].strip()))
                for j in frequencies:
                    line = float(j)
                    if line < 0:
                        negfreq += 1
            elif 'total ZPE' in var:
                zeropointe = float(var.strip().split()[-2])
            elif 'total E vib' in var:
                enthalpye = float(var.strip().split()[-2])
            elif 'total S vib' in var:

        # therme = str(var.strip().split()[-2])

                svib = float(var.strip().split()[-2])

      #  enthalpye = str(var.strip().split()[-1])
      # elif 'Thermal correction to Gibbs Free Energy' in var:
      #  gibbse = str(var.strip().split()[-1])
      # elif 'Sum of electronic and zero-point Energies' in var:
      #  zeropointesum = str(var.strip().split()[-1])
      # elif 'Sum of electronic and thermal Energies' in var:
      #  thermesum = str(var.strip().split()[-1])
      # elif 'Sum of electronic and thermal Enthalpies' in var:
      #  enthalpyesum = str(var.strip().split()[-1])
      # elif 'Sum of electronic and thermal Free Energies' in var:
      #  gibbsesum = str(var.strip().split()[-1])
      #  stotal = str(filelist[i+4].strip().split()[-1])
      #  selec = str(filelist[i+5].strip().split()[-1])
      #  strans = str(filelist[i+6].strip().split()[-1])
      #  srot = str(filelist[i+7].strip().split()[-1])
      #  svib = str(filelist[i+8].strip().split()[-1])

        if 'Job finished:' in var and len(filelist) - i < 10:
            jobdone = 1
        elif 'Job finished:' in var and i - len(filelist) < 10:
            jobdone = 2

    if trans == 'translate':

    # Printed output list
    # printlines = ['\033[1;44m                                                                              \033[1;0m']
    # printlines = ['Analyzing GULP Output File for Supporting Information: ' + logfile]
    # printlines.append('Using ' + str(programversion) + " on #Procs=" + nproc)
    # printlines.append('[#Procs=' + nproc + '  Memory=' + mem + '  Scratch Dir=' + chk + ']')
    # printlines.append('===========================================================')

        printlines.append('Charge = ' + charge + '\tSCF Energy = '
                          + scfenergy)
        if electrons != '':
            printlines.append('Total Valence Electrons = ' + electrons)
    elif trans == 'archive':

    # printlines.append('-------------------------------------------------------------------------------------------------------')
    # printlines.append('Input File:')
    # printlines.append('===========================================================')

    # Initial information

        printlines.append('Supporting Information: ' + logfile)
        printlines.append('-------------------------------------------------------------------------------------------------------'
                          )
        printlines.append('Using: ' + programversion)
        printlines.append('==========================================================='
                          )

  # if trans == "archive" or trans == "translate":
  #  #Print config
  #  printlines.append(''.join(str(e) for e in configlines))
  #
  #  #Modredundant output
  #  if modswitch == 1:
  #    for line in modlines:
  #      printlines.append(line)
#
  #  #Print pointgroup and stoichiometry
  #  printlines.append('------------------------------------------------------------------------------')
  #  #printlines.append('Pointgroup=' + point + '\tStoichiometry=' + stoic + '\t' + pointstoic)
#
  #  #Print Charge Multi and BasisFunc
  #  printlines.append('Charge = ' + charge + '\tSCF Energy = ' + scfenergy)
  #  printlines.append('------------------------------------------------------------------------------')
#
  #  #Print Energy
  #  #if freqswitch == 0:
  #  #  printlines.append('SCF Energy = ' + scfenergy ) #+ '\tPredicted Change = ' + predechange)
  #  #  printlines.append('------------------------------------------------------------------------------')
  #  #
  #  #
  #  ##Optimization output
  #  #if optswitch == 1:
  #  #  if optcomplete == 1 and optstep != 0:
  #  #    printlines.append('\033[1;32mOptimization completed.\033[0;0m  ----------------------------  Optimization step: ' + optstep)
  #  #  elif optcomplete == 1 and optstep == 0:
  #  #    printlines.append('\033[1;32mOptimization completed.\033[0;0m')
  #  #  elif optcomplete == 0 and optstep != 0:
  #  #    printlines.append('\033[1;31mOptimization incomplete.\033[0;0m  ----------------------------  Optimization step: ' + optstep)
  #  #  elif optcomplete == 0 and optstep == 0:
  #  #    printlines.append('\033[1;31mOptimization incomplete.\033[0;0m')
  #  #  if len(maxforce) > 0 and printtype == "none":
  #  #    printlines.append('Item      Max Val.    Criteria    Pass?      RMS Val.    Criteria    Pass?')
  #  #    printlines.append('Energy    ' + energycon[1] + ' || ' + energycon[3] + '    ' + energycon[-1])
  #  #    printlines.append('Force     ' + maxforce[2] + ' || ' + maxforce[4] + '    ' + maxforce[-1] + '      ' + rmsforce[2] + ' || ' + rmsforce[4] + '    ' + rmsforce[-1])
  #  #    printlines.append('Disp      ' + maxdisp[2] + ' || ' + maxdisp[4] + '    ' + maxdisp[-1] + '      ' + rmsdisp[2] + ' || ' + rmsdisp[4] + '    ' + rmsdisp[-1])
  #  #    printlines.append('------------------------------------------------------------------------------')
  #  #  elif len(maxforce) > 0 and printtype == "long":
  #  #    printlines.append('Opt Step:   Energy       Force [Max || RMS]          Displacement [Max || RMS]')
  #  #    for i in range(len(maxforce)):
  #  #      printlines.append(" " + str(i+1) + " " + energycon[i][1] + "(" + energycon[i][-1] + ") [" + maxforce[i][2] + "(" + maxforce[i][-1] + ")||" + rmsforce[i][2] + "(" + rmsforce[i][-1] + ")] [" + maxdisp[i][2] + "(" + maxdisp[i][-1] + ")||" + rmsdisp[i][2] + "(" + rmsdisp[i][-1] + ")]")
  #  #    printlines.append('------------------------------------------------------------------------------')
  #  #  elif len(maxforce) > 0 and printtype == "extra":
  #  #    printlines.append('Opt Step:   Energy       Force [Max || RMS]          Displacement [Max || RMS]')
  #  #    if len(maxforce) > 5:
  #  #      lastlines = len(maxforce) - 5
  #  #    else:
  #  #      lastlines = 0
  #  #    for i in range(lastlines,len(maxforce)):
  #  #      printlines.append(" " + str(i+1) + " " + energycon[i][1] + "(" + energycon[i][-1] + ") [" + maxforce[i][2] + "(" + maxforce[i][-1] + ")||" + rmsforce[i][2] + "(" + rmsforce[i][-1] + ")] [" + maxdisp[i][2] + "(" + maxdisp[i][-1] + ")||" + rmsdisp[i][2] + "(" + rmsdisp[i][-1] + ")]")
  #  #    printlines.append('------------------------------------------------------------------------------')

    if trans == 'archive':

        printlines.append('    Atomic        Coordinates (Angstroms)')
        printlines.append('    Type:            X       Y       Z')
        printlines.append('------------------------------------------------------------------------------'
                          )
        printlines.append('GEOM DATA!')

    if trans == 'translate' or trans == 'archive':

    # Frequency output

        if freqswitch == 1 and enthalpye != '':
            printlines.append('Statistical Thermodynamic Analysis for '
                              + logfile)
            printlines.append(freqtemp)
            printlines.append('=============================================================================='
                              )

      # printlines.append('SCF Energy=                                             ' + scfenergy)

            printlines.append('Zero-point correction (ZPE)=            '
                               + str(zeropointe) + '  J/mol')  # + zeropointesum)

      # printlines.append('Internal Energy (U)=                    ' + therme + '  J/mol' + thermesum)

            printlines.append('Enthalpy (H)=                           '
                               + str(enthalpye) + '        ')  # + enthalpyesum)

      # printlines.append('Gibbs Free Energy (G)=                  ' + gibbse + '        ' + gibbsesum)
      # printlines.append('Entropy (S) [Cal/Mol-Kelvin]=           ' + stotal)
      # printlines.append('Entropy (S): electronic=                ' + selec)
      # printlines.append('Entropy (S): translational=             ' + strans)
      # printlines.append('Entropy (S): rotational=                ' + srot)

            printlines.append('Entropy (S): vibrational=               '
                               + str(svib) + '  J/mol/K')
            printlines.append('Frequencies: ' + frequencies[0] + '  '
                              + frequencies[1] + '  ' + frequencies[2]
                              + '  ' + frequencies[3] + '  '
                              + frequencies[4])
            if negfreq == 1:
                printlines.append('\033[1;33m              !!!Warning there were '
                                   + str(negfreq)
                                  + ' negative frequencies!!!\033[0;m')
            elif negfreq >= 2:
                printlines.append('\033[1;31m              !!!Warning there were '
                                   + str(negfreq)
                                  + ' negative frequencies!!!\033[0;m')
            printlines.append('------------------------------------------------------------------------------'
                              )

    if trans == 'translate':

    # Print end of file for long and extra

        if printtype == 'extra' or printtype == 'long':
            for i in range(len(filelist) - 10, len(filelist)):
                printlines.append(filelist[i].strip())

    # Completed job

        if jobdone == 1:
            printlines.append('=============================================================================='
                              )
            printlines.append('\033[1;32m##################              JOB COMPLETED               ##################\033[0;0m'
                              )
            printlines.append('=============================================================================='
                              )
        elif jobdone == 2:
            printlines.append('=============================================================================='
                              )
            printlines.append('\033[1;31m!!!!!!!!!!!!!!!!!!            !!!JOB FAILED!!!              !!!!!!!!!!!!!!!!!!\033[0;0m'
                              )
            printlines.append('=============================================================================='
                              )
    elif trans == 'list':

        line = scfenergy + '  '
        line += zeropointesum + '  ' + thermesum + '  ' + enthalpyesum \
            + '  ' + gibbsesum + '  '
        if negfreq == 0:
            line += '\033[1;32m' + frequencies[0] + '  ' \
                + frequencies[1] + '\033[1;0m' + '  '
        elif negfreq == 1:
            line += '\033[1;33m' + frequencies[0] + '  ' \
                + frequencies[1] + '\033[1;0m' + '  '
        elif negfreq >= 2:
            line += '\033[1;31m' + frequencies[0] + '  ' \
                + frequencies[1] + '\033[1;0m' + '  '
        line += basisfunc + '  ' + pointstoic + '  ' + stoic + '  '
        if jobdone == 0:
            line += logfile.split('/')[-1]
        elif jobdone == 1:
            filename = logfile.split('/')[-1]
            line += '\033[1;32m' + filename + '\033[0;0m'
        elif jobdone == 2:
            filename = logfile.split('/')[-1]
            line += '\033[1;31m' + filename + '\033[0;0m'
        printlines.append(line)

    if trans != 'archive':
        return printlines
    if trans == 'archive':
        entropyesum = str(float(svib) / 1000 / 2625.5)
        enthalpyesum = str(float(enthalpye) / 1000 / 2625.5)
        entropyehart = float(svib) / 1000 / 2625.5
        enthalpyehart = enthalpye / 1000 / 2625.5
        if len(freqtemp) > 0:
            gibbse = float(svib) / 1000 / 2625.5 \
                * float(freqtemp.split()[1]) + enthalpye / 1000 / 2625.5
            gibbsesum = gibbse
        else:
            gibbse = 0
            gibbsesum = 0
        tablelines = []
        tablelines.append(logfile.split('/')[-1] + ',' + scfenergy + ','
                           + str(zeropointe) + ',=' + str(entropyehart)
                          + '+RC[-2],=' + str(enthalpyehart)
                          + '+RC[-3],=' + str(gibbse)
                          + '+RC[-4],=RC[-1]-MIN(C[-1]),=RC[-1]*627.5095'
                          )
        tablelines.append(logfile.split('/')[-1] + ',' + scfenergy + ','
                           + str(zeropointe) + ',' + str(thermesum)
                          + ',' + str(enthalpyesum) + ','
                          + str(gibbsesum)
                          + ',=RC[-1]-MIN(C[-1]),=RC[-1]*627.5095')
        tablelines.append(logfile.split('/')[-1] + ',' + scfenergy)
        return (printlines, tablelines)


