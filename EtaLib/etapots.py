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

import os
import shutil
import time
from sys import *

from numpy import *

from . import etaatom, etanumpy

###==================================================================================================================###
### --- Atom class --- ###


class PersistentOTS(object):
    def __init__(self):

        # ## --- OPT Variables --- ###

        self.stepNumber = 300
        self.maxopt = 25
        self.breaktype = ""
        self.breakpoint = 0.0
        self.breakatoms = []
        self.optimized = False
        self.kill = False
        self.inputnumatoms = 0
        self.shortenStepSize = 0

        # ## --- Scan File Variables --- ###

        self.scanType = ""
        self.atomList = []
        self.scan = 0
        self.ts = 0
        self.minorTS = 0
        self.scanLength = 0
        self.scanStepNow = 0
        self.scanSize = 0
        self.scanSizeNOW = 0
        self.waitSize = 0
        self.breakpoint = 0.0
        self.breakatoms = []
        self.minTSScanVar = "none"

        # ## --- (TS-)SCAN Variables --- ###

        self.minTSScanSize = 0.0
        self.scanEnergy = []
        self.tsscanEnergy = []

        # ## --- Molecule Variables --- ###

        self.finalenergy = 0
        self.finalForce = 0
        self.finalDisp = 0
        self.allEnergy = []
        self.allForce = []
        self.numatoms = 0
        self.charge = 0
        self.multi = 1
        self.inputnumatoms = 0

        # ## --- Input File Variables (mostly for G0X) --- ###

        self.QMconfig = []
        self.QMconfigORIG = []
        self.QMcomment = ""
        self.QMcommentORIG = ""
        self.QMendlines = []
        self.QMendlinesORIG = []
        self.TCscratchDir = ""

        # ## --- Geometry holders (Atom class) --- ###

        self.inputGeomAtom = []
        self.finalGeomAtom = []
        self.trajectoryGeomAtom = []
        self.highestGeomAtom = []
        self.ZMATAtom = []

        # ## --- PersistentOTS variables --- ###

        self.baseinput = ""
        self.baselog = ""
        self.joblog = ""
        self.scratchdir = ""
        self.trajectoryfile = ""
        self.scanfile = ""
        self.outScanFile = ""
        self.linearScanFile = ""
        self.program = "nPersistentOTS.py"
        self.inputfile = ""
        self.outfile = ""
        self.jobname = ""
        self.pwdirectory = ""
        self.jobtype = ""
        self.QMtype = ""
        self.path2QM = ""
        self.debug = 0
        self.test = 0
        self.verbose = 0


# Constants

auToAngstrom = 1.889725988579
hartreeTokcalmol = 627.509469
hartreeToeV = 27.21138386
kcalmolTokjoulemol = 4.184


### --- Setup and reset all the variables --- ###


def reset_pots_variables(pots):

    # ## --- OPT Variables --- ###
    # pots.stepNumber = 300
    # pots.maxopt = 25

    pots.breaktype = ""
    pots.breakpoint = 0.0
    pots.breakatoms = []
    pots.optimized = False

    # pots.kill = False

    pots.inputnumatoms = 0
    pots.shortenStepSize = 0

    # ## --- Scan File Variables --- ###

    pots.scanType = ""
    pots.atomList = []
    pots.scan = 0
    pots.ts = 0
    pots.minorTS = 0
    pots.scanLength = 0
    pots.scanStepNow = 0
    pots.scanSize = 0
    pots.scanSizeNOW = 0
    pots.waitSize = 0
    pots.breakpoint = 0.0
    pots.breakatoms = []
    pots.minTSScanVar = "none"

    # ## --- (TS-)SCAN Variables --- ###

    pots.minTSScanSize = 0.0

    # pots.scanEnergy = []
    # pots.tsscanEnergy = []

    # ## --- Molecule Variables --- ###

    pots.finalenergy = 0
    pots.finalForce = 0
    pots.finalDisp = 0
    pots.allEnergy = []
    pots.allForce = []

    # pots.numatoms = 0
    # pots.charge = 0
    # pots.multi = 1

    pots.inputnumatoms = 0

    # ## --- Input File Variables (mostly for G0X) --- ###

    pots.QMconfig = []

    # pots.QMconfigORIG = []

    pots.QMcomment = ""

    # pots.QMcommentORIG = ''

    pots.QMendlines = []

    # pots.QMendlinesORIG = []

    pots.TCscratchDir = ""

    # ## --- Geometry holders (Atom class) --- ###

    pots.inputGeomAtom = []
    pots.finalGeomAtom = []
    pots.trajectoryGeomAtom = []

    # pots.highestGeomAtom = []

    pots.ZMATAtom = []

    # ## --- PersistentOTS variables --- ###
    # pots.baseinput = ''
    # pots.baselog = ''
    # pots.joblog = ''
    # pots.scratchdir = ''
    # pots.trajectoryfile = ''
    # pots.scanfile = ''
    # pots.outScanFile = ''
    # pots.linearScanFile = ''
    # pots.program = "nPersistentOTS.py"
    # pots.inputfile = ''
    # pots.outfile = ''
    # pots.jobname = ''
    # pots.pwdirectory = ''
    # pots.jobtype = ''
    # pots.QMtype = ''
    # pots.path2QM = ''
    # pots.debug = 0
    # pots.test = 0

    return pots


### --- FUNCTIONS --- ###
########################################################################################################################
### --- Parse the scan file  --- ###


def parse_scan_file(pots):
    if pots.debug >= 1:
        print("parse_scan_file")
    f = open(pots.scanfile, "r")
    scanlist = f.readlines()
    f.close()
    for var in scanlist:
        line = var.split("=")
        if line[0].strip() == "#STEPSIZE":
            pots.scanSize = float(line[1].strip())
        elif line[0].strip() == "#SCAN":
            if line[1].split()[0].upper() == "YES":
                pots.scan = 1
        elif line[0].strip() == "#TS":
            if line[1].split()[0].upper() == "YES":
                pots.ts = 1
            if line[1].split()[0].upper() == "DYN":
                pots.ts = 2
        elif line[0].strip() == "#SCANLENGTH":
            pots.scanLength = int(line[1].split()[0])
        elif line[0].strip() == "#SCANSTEPNOW":
            pots.scanStepNow = int(line[1].split()[0])
        elif line[0].strip() == "#MINTSSCANSIZE":
            pots.minTSScanVar = str(line[1].split()[0].lower())
        elif line[0].strip() == "#WAIT":
            pots.waitSize = int(line[1].split()[0])
        elif line[0].strip() == "#VARIABLE":
            line2 = line[1].split()
            pots.scanType = line2[0].lower()
            for i in range(1, len(line2)):
                pots.atomList.append(int(line2[i]))
        elif line[0].strip() == "#BREAKPOINT":
            line2 = line[1].split()
            pots.breakpoint = float(line2[0])
            for i in range(1, len(line2)):
                pots.breakatoms.append(int(line2[i]))
        elif line[0].strip() == "#MINORTS":
            if line[1].split()[0].upper() == "ONE":
                pots.minorTS = 0
            if line[1].split()[0].upper() == "YES":
                pots.minorTS = 2
            if line[1].split()[0].upper() == "NO":
                pots.minorTS = 3
        if pots.QMtype == "Terachem":
            if var.strip() != "":
                line2 = var.strip().split()
                if line2[0].strip() == "#TCSCRDIR":
                    pots.TCscratchDir = (line2[1])[1:]
    if pots.scan == 0:
        print_and_write(
            [
                "The scan feature is turned off in the scan text file.",
                "Please turn this on (#SCAN=yes) if you want this feature.",
            ],
            pots.verbose,
            pots.joblog,
        )
    if pots.debug >= 1:
        print((pots.scanSize))
        print((pots.scan))
        print((pots.ts))
        print((pots.scanLength))
        print((pots.scanStepNow))
        print((pots.minTSScanVar))
        print((pots.waitSize))
        print((pots.scanType))
        print((pots.atomList))
        print((pots.breakpoint))
        print((pots.breakatoms))
        print((pots.minorTS))
        print((pots.TCscratchDir))
    return pots


### --- Parse the scan file  --- ###


def make_temp_scan_file(pots):
    if pots.debug >= 1:
        print("make_temp_scan_file")
    if os.path.isfile("Scan.txt"):
        if pots.QMtype == "G09" or pots.QMtype == "G03":
            for i in os.listdir(os.getcwd()):
                if i.endswith(".com"):
                    ofile = "Scan-" + etaatom.basename(i, ".com") + ".txt"
                    shutil.copyfile("Scan.txt", ofile)
                    print(ofile)
        elif pots.QMtype == "Turbomole":
            for i in os.listdir(os.getcwd()):
                if i.endswith(".Turbomole"):
                    ofile = (
                        etaatom.basename(i, ".Turbomole")
                        + "/Scan-"
                        + etaatom.basename(i, ".Turbomole")
                        + ".txt"
                    )
                    shutil.copyfile("Scan.txt", ofile)
                    print(ofile)
            for i in os.listdir(os.getcwd()):
                if i.endswith(".tmol"):
                    ofile = "Scan-" + etaatom.basename(i, ".tmol") + ".txt"
                    shutil.copyfile("Scan.txt", ofile)
                    print(ofile)
        elif pots.QMtype == "Terachem":
            for i in os.listdir(os.getcwd()):
                if i.endswith(".tc"):
                    ofile = "Scan-" + etaatom.basename(i, ".tc") + ".txt"
                    shutil.copyfile("Scan.txt", ofile)
                    print(ofile)
        elif pots.QMtype.upper() == "JAGUAR":
            for i in os.listdir(os.getcwd()):
                if i.endswith(".in"):
                    ofile = "Scan-" + etaatom.basename(i, ".in") + ".txt"
                    shutil.copyfile("Scan.txt", ofile)
                    print(ofile)
        shutil.move("Scan.txt", "Scan.distributed")
    else:
        f = open("Scan.txt", "w")
        print("#SCAN = YES/NO")
        print("#TS = YES/DYN/NO")
        print("#VARIABLE = bond/angle/dihedral 1 2 3 4")
        print("#SCANLENGTH = 30")
        print("#SCANSTEPNOW = 0")
        print("#STEPSIZE = -0.05")
        print("#MINTSSCANSIZE = loose/none/tight/verytight")
        print("#WAIT = 0")
        print("#BREAKPOINT: 2.3 2 4")
        print("#MINORTS = YES/ONE/NO")
        print("#TCSCRDIR = ./scr-h2o")
        print(
            "============================================================================"
        )
        print(
            "Edit the new Scan.txt file and rerun nPersistantOTS.py -T and the Scan.txt "
        )
        print("will be distributed for all input files.")

        # Write the temp scan file

        f.write("#SCAN = YES\n")
        f.write("#TS = YES\n")
        f.write("#VARIABLE = bond 1 2\n")
        f.write("#SCANLENGTH = 30\n")
        f.write("#SCANSTEPNOW = 0\n")
        f.write("#STEPSIZE = -0.05\n")
        f.write("#MINTSSCANSIZE = none\n")
        f.write("#WAIT = 0\n")
        f.write("#BREAKPOINT: 2.3 2 4\n")
        f.write("#MINORTS = ONE\n")
        f.write("#TCSCRDIR = \n")
        f.close()
    return pots


### --- Write the final geom to the SCAN xyz file. This can be visualized by molden or others--- ###


def print_and_write(alist, verbose, sfile):
    f = open(sfile, "a")
    for i in alist:
        if verbose == 1:
            print(i)
        elif verbose == 3:
            print((i.center(80, " ")))
        f.write(str(i) + "\n")
    f.close()
    return


#### --- Write the trajectory of the last optimization to the trajectory file--- ###
# def write_trajectory(pots):
#  if pots.debug >= 1:
#    print_and_write(["write_trajectory"], pots.verbose, pots.joblog)
#  #List of ALL
#  f = open(pots.trajectoryfile, 'a')
#  for i in range(len(pots.trajectoryGeomAtom)):
#    #if the instance is in Atom form print the 'element X-Coord Y-Coord Z-Coord'
#    if isinstance(ifilelol[i], etaatom.Atom):
#      line = ifilelol[i].e + "  " + str("{:.6f}".format(ifilelol[i].x)) + \
#             "  " + str("{:.6f}".format(ifilelol[i].y)) + "  " + str("{:.6f}".format(ifilelol[i].z))
#    #else print the string (number of atoms or comment line)
#    else:
#      line = str(ifilelol[i])
#    f.write(line + "\n")
#  f.close()
#  return

### --- Write the final geom to the SCAN xyz file. This can be visualized by molden or others--- ###


def write_final_geom_to_scan(pots):
    if pots.debug >= 1:
        print_and_write(["write_final_geom_to_scan"], pots.verbose, pots.joblog)

    # List of ALL

    scangeom = []

    # List of Bond Angle Dihedral info

    scanbad = []
    if os.path.isfile(pots.outScanFile):

        # Grab the scan geometries Bond Angle Dihedral (BAD) numbers

        f = open(pots.outScanFile, "r")
        for line in f:
            scangeom.append(line)
            if line.rstrip().split()[0] == "Energy:":
                scanbad.append(float(line.rstrip().split()[3]))
        f.close()

        # Output the Geometries

        usedgeom = 0
        f = open(pots.outScanFile, "w+")
        for i in range(len(scanbad)):
            if (
                float(pots.finalGeomAtom[1].rstrip().split()[3]) < float(scanbad[i])
                and usedgeom == 0
            ):
                usedgeom = 1

                # Write finalgeom

                for j in range(len(pots.finalGeomAtom)):

                    # if the instance is in Atom form print the 'element X-Coord Y-Coord Z-Coord'

                    if isinstance(pots.finalGeomAtom[j], etaatom.Atom):
                        line = (
                            pots.finalGeomAtom[j].e
                            + "  "
                            + str("{:.6f}".format(pots.finalGeomAtom[j].x))
                            + "  "
                            + str("{:.6f}".format(pots.finalGeomAtom[j].y))
                            + "  "
                            + str("{:.6f}".format(pots.finalGeomAtom[j].z))
                        )
                    else:

                        # else print the string (number of atoms or comment line)

                        line = str(pots.finalGeomAtom[j])
                    f.write(line + "\n")

                # Write scangeom

                for j in range(len(scangeom)):
                    if (
                        len(pots.finalGeomAtom) * i
                        <= j
                        < len(pots.finalGeomAtom) * (i + 1)
                    ):
                        f.write(scangeom[j])
            else:

                # Write scangeom

                for j in range(len(scangeom)):
                    if (
                        len(pots.finalGeomAtom) * i
                        <= j
                        < len(pots.finalGeomAtom) * (i + 1)
                    ):
                        f.write(scangeom[j])
        if usedgeom == 0:

            # Write finalgeom

            for j in range(len(pots.finalGeomAtom)):

                # if the instance is in Atom form print the 'element X-Coord Y-Coord Z-Coord'

                if isinstance(pots.finalGeomAtom[j], etaatom.Atom):
                    line = (
                        pots.finalGeomAtom[j].e
                        + "  "
                        + str("{:.6f}".format(pots.finalGeomAtom[j].x))
                        + "  "
                        + str("{:.6f}".format(pots.finalGeomAtom[j].y))
                        + "  "
                        + str("{:.6f}".format(pots.finalGeomAtom[j].z))
                    )
                else:

                    # else print the string (number of atoms or comment line)

                    line = str(pots.finalGeomAtom[j])
                f.write(line + "\n")
    else:
        f = open(pots.outScanFile, "w+")

        # Write finalgeom

        for j in range(len(pots.finalGeomAtom)):

            # if the instance is in Atom form print the 'element X-Coord Y-Coord Z-Coord'

            if isinstance(pots.finalGeomAtom[j], etaatom.Atom):
                line = (
                    pots.finalGeomAtom[j].e
                    + "  "
                    + str("{:.6f}".format(pots.finalGeomAtom[j].x))
                    + "  "
                    + str("{:.6f}".format(pots.finalGeomAtom[j].y))
                    + "  "
                    + str("{:.6f}".format(pots.finalGeomAtom[j].z))
                )
            else:

                # else print the string (number of atoms or comment line)

                line = str(pots.finalGeomAtom[j])
            f.write(line + "\n")
    f.close()
    return


### --- Write the final geom to the SCAN xyz file. This can be visualized by molden or others--- ###


def write_final_geom_to_scan_linear(pots):
    if pots.debug >= 1:
        print_and_write(["write_final_geom_to_scan_linear"], pots.verbose, pots.joblog)
    if os.path.isfile(pots.linearScanFile):
        f = open(pots.linearScanFile, "a")
    else:
        f = open(pots.linearScanFile, "w+")

    # Write finalgeom

    for i in range(len(pots.finalGeomAtom)):

        # if the instance is in Atom form print the 'element X-Coord Y-Coord Z-Coord'

        if isinstance(pots.finalGeomAtom[i], etaatom.Atom):
            line = (
                pots.finalGeomAtom[i].e
                + "  "
                + str("{:.6f}".format(pots.finalGeomAtom[i].x))
                + "  "
                + str("{:.6f}".format(pots.finalGeomAtom[i].y))
                + "  "
                + str("{:.6f}".format(pots.finalGeomAtom[i].z))
            )
        else:

            # else print the string (number of atoms or comment line)

            line = str(pots.finalGeomAtom[i])
        f.write(line + "\n")
    f.close()
    return


### --- Write the final geom to the SCAN xyz file. This can be visualized by molden or others--- ###


def write_high_geom(pots):
    if pots.debug >= 1:
        print_and_write(["write_high_geom"], pots.verbose, pots.joblog)
    if os.path.isfile(pots.scratchdir + "/HIGHEST-" + pots.jobname + ".xyz"):
        shutil.move(
            pots.scratchdir + "/HIGHEST-" + pots.jobname + ".xyz",
            pots.scratchdir + "/SECOND-" + pots.jobname + ".xyz",
        )
    f = open(pots.scratchdir + "/HIGHEST-" + pots.jobname + ".xyz", "w+")

    # Write finalgeom

    for i in range(len(pots.finalGeomAtom)):

        # if the instance is in Atom form print the 'element X-Coord Y-Coord Z-Coord'

        if isinstance(pots.finalGeomAtom[i], etaatom.Atom):
            line = (
                pots.finalGeomAtom[i].e
                + "  "
                + str("{:.6f}".format(pots.finalGeomAtom[i].x))
                + "  "
                + str("{:.6f}".format(pots.finalGeomAtom[i].y))
                + "  "
                + str("{:.6f}".format(pots.finalGeomAtom[i].z))
            )
        else:

            # else print the string (number of atoms or comment line)

            line = str(pots.finalGeomAtom[i])
        f.write(line + "\n")
    f.close()
    return


### --- Grabs the highest energy geometry held in the SCRATCH folder  --- ###


def save_highest_geom(pots):
    if pots.debug >= 1:
        print_and_write(
            ["save_highest_geom", "Saving this geometry as the highest."],
            pots.verbose,
            pots.joblog,
        )
    pots.highestGeomAtom = [0] * len(pots.finalGeomAtom)
    for i in range(len(pots.finalGeomAtom)):
        if i <= 1:
            pots.highestGeomAtom[i] = pots.finalGeomAtom[i]
        elif i >= 2:
            pots.highestGeomAtom[i] = etaatom.Atom()
            pots.highestGeomAtom[i].e = pots.finalGeomAtom[i].e
            pots.highestGeomAtom[i].en = pots.finalGeomAtom[i].en
            pots.highestGeomAtom[i].x = pots.finalGeomAtom[i].x
            pots.highestGeomAtom[i].y = pots.finalGeomAtom[i].y
            pots.highestGeomAtom[i].z = pots.finalGeomAtom[i].z
            pots.highestGeomAtom[i].int = pots.finalGeomAtom[i].int
            pots.highestGeomAtom[i].en = pots.finalGeomAtom[i].en
            pots.highestGeomAtom[i].neighbors = pots.finalGeomAtom[i].neighbors
            pots.highestGeomAtom[i].neighborsdist = pots.finalGeomAtom[i].neighborsdist
            pots.highestGeomAtom[i].nearest = pots.finalGeomAtom[i].nearest
            pots.highestGeomAtom[i].nearestdist = pots.finalGeomAtom[i].nearestdist
            pots.highestGeomAtom[i].hybridization = pots.finalGeomAtom[i].hybridization
            pots.highestGeomAtom[i].charge = pots.finalGeomAtom[i].charge
            pots.highestGeomAtom[i].string = pots.finalGeomAtom[i].string
            pots.highestGeomAtom[i].atomtype = pots.finalGeomAtom[i].atomtype
    return pots


### --- Grabs the highest energy geometry held in the SCRATCH folder  --- ###


def grab_highest_geom(pots):
    if pots.debug >= 1:
        print_and_write(["grab_highest_geom"], pots.verbose, pots.joblog)
    for i in range(len(pots.highestGeomAtom)):
        if i <= 1:
            pots.highestGeomAtom[i] = pots.finalGeomAtom[i]
        elif i >= 2:
            pots.finalGeomAtom[i] = etaatom.Atom()
            pots.finalGeomAtom[i].e = pots.highestGeomAtom[i].e
            pots.finalGeomAtom[i].en = pots.highestGeomAtom[i].en
            pots.finalGeomAtom[i].x = pots.highestGeomAtom[i].x
            pots.finalGeomAtom[i].y = pots.highestGeomAtom[i].y
            pots.finalGeomAtom[i].z = pots.highestGeomAtom[i].z
            pots.finalGeomAtom[i].int = pots.highestGeomAtom[i].int
            pots.finalGeomAtom[i].en = pots.highestGeomAtom[i].en
            pots.finalGeomAtom[i].neighbors = pots.highestGeomAtom[i].neighbors
            pots.finalGeomAtom[i].neighborsdist = pots.highestGeomAtom[i].neighborsdist
            pots.finalGeomAtom[i].nearest = pots.highestGeomAtom[i].nearest
            pots.finalGeomAtom[i].nearestdist = pots.highestGeomAtom[i].nearestdist
            pots.finalGeomAtom[i].hybridization = pots.highestGeomAtom[i].hybridization
            pots.finalGeomAtom[i].charge = pots.highestGeomAtom[i].charge
            pots.finalGeomAtom[i].string = pots.highestGeomAtom[i].string
            pots.finalGeomAtom[i].atomtype = pots.highestGeomAtom[i].atomtype
    return pots


#  if os.path.isfile(pots.scratchdir + "/HIGHEST-" + pots.jobname + ".xyz"):
#    alist = []
#    the_file = open(pots.scratchdir + "/HIGHEST-" + pots.jobname + ".xyz",'r')
#    for line in the_file:
#      alist.append(line.rstrip())
#    the_file.close()
#    for i in range(len(alist)):
#      if i <= 1:
#        pots.finalGeomAtom[i] = alist[i]
#      elif i >= 2:
#        line = alist[i].split()
#        pots.finalGeomAtom[i] = etaatom.Atom()
#        pots.finalGeomAtom[i].e = line[0]
#        #pots.finalGeomAtom[i].en = etaatom.get_element_num(alist[j][0])
#        pots.finalGeomAtom[i].x = float(line[1])
#        pots.finalGeomAtom[i].y = float(line[2])
#        pots.finalGeomAtom[i].z = float(line[3])
#  else:
#    print_and_write(["Highest geom file does not exist."], pots.verbose, pots.joblog)
#  return pots

### --- Grabs a previous nth geometry from a xyz file. geomNum  --- ###


def grab_previous_geom(pots, geomnum):
    if pots.debug >= 1:
        print_and_write(["grab_previous_geom"], pots.verbose, pots.joblog)
    alist = []
    the_file = open(pots.linearScanFile, "r")
    for line in the_file:
        alist.append(line.rstrip())
    the_file.close()

    # Extract the final geometry from the trajectory

    pots.finalGeomAtom = [0] * (pots.trajectoryGeomAtom[0] + 2)
    for i in range(
        len(alist) - geomnum * len(pots.finalGeomAtom),
        len(alist) - (geomnum - 1) * len(pots.finalGeomAtom),
    ):
        j = i - (len(alist) - geomnum * len(pots.finalGeomAtom))
        line = alist[j].split()
        if j <= 1:
            pots.finalGeomAtom[j] = alist[i]
        elif j >= 2:
            pots.finalGeomAtom[j] = etaatom.Atom()
            pots.finalGeomAtom[j].e = line[0]

            # pots.finalGeomAtom[i].en = etaatom.get_element_num(alist[j][0])

            pots.finalGeomAtom[j].x = float(line[1])
            pots.finalGeomAtom[j].y = float(line[2])
            pots.finalGeomAtom[j].z = float(line[3])
    return pots


### --- Make a backup copy of the original geom in Atom form --- ###


def make_geom_backup(pots):
    backupgeom = [0] * len(pots.finalGeomAtom)
    backupgeom[0] = pots.finalGeomAtom[0]
    backupgeom[1] = pots.finalGeomAtom[1]
    for i in range(2, len(pots.finalGeomAtom)):
        backupgeom[i] = etaatom.Atom()

        # print_and_write([pots.finalGeomAtom[i].e], pots.verbose, pots.joblog)

        backupgeom[i].e = pots.finalGeomAtom[i].e
        backupgeom[i].en = pots.finalGeomAtom[i].en
        backupgeom[i].x = pots.finalGeomAtom[i].x
        backupgeom[i].y = pots.finalGeomAtom[i].y
        backupgeom[i].z = pots.finalGeomAtom[i].z
        backupgeom[i].int = pots.finalGeomAtom[i].int
        backupgeom[i].en = pots.finalGeomAtom[i].en
        backupgeom[i].neighbors = pots.finalGeomAtom[i].neighbors
        backupgeom[i].neighborsdist = pots.finalGeomAtom[i].neighborsdist
        backupgeom[i].nearest = pots.finalGeomAtom[i].nearest
        backupgeom[i].nearestdist = pots.finalGeomAtom[i].nearestdist
        backupgeom[i].hybridization = pots.finalGeomAtom[i].hybridization
        backupgeom[i].charge = pots.finalGeomAtom[i].charge
        backupgeom[i].string = pots.finalGeomAtom[i].string
        backupgeom[i].atomtype = pots.finalGeomAtom[i].atomtype

    return backupgeom


### --- Given a number and a list, how many entries in the list are smaller than the given number? --- ###


def number_smaller(number, alist):
    numsmall = 0
    for i in alist:
        if number < i:
            numsmall += 1
    return numsmall


### --- Given a number and a list, how many entries in the list are larger than the given number? --- ###


def number_larger(number, alist):
    numlarge = 0
    for i in alist:
        if number > i:
            numlarge += 1
    return numlarge


### --- Move the atoms of interest to the front of the file for subsequent changes --- ###


def renumber_fix(pots):
    if pots.debug >= 1:
        print_and_write(["renumber_fix"], pots.verbose, pots.joblog)

    # Place the original locations into each atom as an int for unfixing the geometry

    for i in range(2, len(pots.finalGeomAtom)):
        pots.finalGeomAtom[i].int = i

    # Make a backup geom to work from

    oldfinalgeomatom = make_geom_backup(pots)

    # oldfinalgeomatom = pots.finalGeomAtom
    # print_geom(oldfinalgeomatom)
    # Move the atoms of interest to the tom of the molecule file in the order they are listed
    # (needed for accurate angle and diherdral SCANs)

    for i in range(len(pots.atomList)):
        for j in range(2, len(pots.finalGeomAtom)):
            if j - 1 == pots.atomList[i]:
                pots.finalGeomAtom[i + 2].e = oldfinalgeomatom[j].e
                pots.finalGeomAtom[i + 2].en = oldfinalgeomatom[j].en
                pots.finalGeomAtom[i + 2].x = oldfinalgeomatom[j].x
                pots.finalGeomAtom[i + 2].y = oldfinalgeomatom[j].y
                pots.finalGeomAtom[i + 2].z = oldfinalgeomatom[j].z
                pots.finalGeomAtom[i + 2].int = oldfinalgeomatom[j].int
                pots.finalGeomAtom[i + 2].en = oldfinalgeomatom[j].en
                pots.finalGeomAtom[i + 2].neighbors = oldfinalgeomatom[j].neighbors
                pots.finalGeomAtom[i + 2].neighborsdist = oldfinalgeomatom[
                    j
                ].neighborsdist
                pots.finalGeomAtom[i + 2].nearest = oldfinalgeomatom[j].nearest
                pots.finalGeomAtom[i + 2].nearestdist = oldfinalgeomatom[j].nearestdist
                pots.finalGeomAtom[i + 2].hybridization = oldfinalgeomatom[
                    j
                ].hybridization
                pots.finalGeomAtom[i + 2].charge = oldfinalgeomatom[j].charge
                pots.finalGeomAtom[i + 2].string = oldfinalgeomatom[j].string
                pots.finalGeomAtom[i + 2].atomtype = oldfinalgeomatom[j].atomtype

    # print_geom(pots.finalGeomAtom)
    # For all other atoms move them down the appropriate number of places

    for i in range(2, len(pots.finalGeomAtom)):
        if i - 1 not in pots.atomList:
            j = number_smaller(i - 1, pots.atomList)

            # k = number_larger(i - 1, pots.atomList)

            pots.finalGeomAtom[i + j].e = oldfinalgeomatom[i].e
            pots.finalGeomAtom[i + j].en = oldfinalgeomatom[i].en
            pots.finalGeomAtom[i + j].x = oldfinalgeomatom[i].x
            pots.finalGeomAtom[i + j].y = oldfinalgeomatom[i].y
            pots.finalGeomAtom[i + j].z = oldfinalgeomatom[i].z
            pots.finalGeomAtom[i + j].int = oldfinalgeomatom[i].int
            pots.finalGeomAtom[i + j].en = oldfinalgeomatom[i].en
            pots.finalGeomAtom[i + j].neighbors = oldfinalgeomatom[i].neighbors
            pots.finalGeomAtom[i + j].neighborsdist = oldfinalgeomatom[i].neighborsdist
            pots.finalGeomAtom[i + j].nearest = oldfinalgeomatom[i].nearest
            pots.finalGeomAtom[i + j].nearestdist = oldfinalgeomatom[i].nearestdist
            pots.finalGeomAtom[i + j].hybridization = oldfinalgeomatom[i].hybridization
            pots.finalGeomAtom[i + j].charge = oldfinalgeomatom[i].charge
            pots.finalGeomAtom[i + j].string = oldfinalgeomatom[i].string
            pots.finalGeomAtom[i + j].atomtype = oldfinalgeomatom[i].atomtype

    # print_geom(pots.finalGeomAtom)

    return pots


### --- Move the atoms of interest back to the correct locations--- ###


def renumber_unfix(pots):
    if pots.debug >= 1:
        print_and_write(["renumber_unfix"], pots.verbose, pots.joblog)
    oldfinalgeomatom = make_geom_backup(pots)

    # print_geom(oldfinalgeomatom)
    # Using the original location saved in oldfinalgeomatom[i].int place the atoms back in their correct locations

    for i in range(2, len(pots.finalGeomAtom)):
        j = oldfinalgeomatom[i].int
        pots.finalGeomAtom[j].e = oldfinalgeomatom[i].e
        pots.finalGeomAtom[j].en = oldfinalgeomatom[i].en
        pots.finalGeomAtom[j].x = oldfinalgeomatom[i].x
        pots.finalGeomAtom[j].y = oldfinalgeomatom[i].y
        pots.finalGeomAtom[j].z = oldfinalgeomatom[i].z
        pots.finalGeomAtom[j].int = oldfinalgeomatom[i].int
        pots.finalGeomAtom[j].en = oldfinalgeomatom[i].en
        pots.finalGeomAtom[j].neighbors = oldfinalgeomatom[i].neighbors
        pots.finalGeomAtom[j].neighborsdist = oldfinalgeomatom[i].neighborsdist
        pots.finalGeomAtom[j].nearest = oldfinalgeomatom[i].nearest
        pots.finalGeomAtom[j].nearestdist = oldfinalgeomatom[i].nearestdist
        pots.finalGeomAtom[j].hybridization = oldfinalgeomatom[i].hybridization
        pots.finalGeomAtom[j].charge = oldfinalgeomatom[i].charge
        pots.finalGeomAtom[j].string = oldfinalgeomatom[i].string
        pots.finalGeomAtom[j].atomtype = oldfinalgeomatom[i].atomtype

    # print_geom(pots.finalGeomAtom)

    return pots


###


def generate_min_ts_scan_size(pots):

    # Set the minimum step size allowed if non was set manually

    if pots.minTSScanVar.lower().strip() == "loose":
        x = 100
        if pots.scanType == "bond":
            pots.minTSScanSize = round(
                x
                / 100.0
                / (
                    float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                    + float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                ),
                4,
            )
        elif pots.scanType == "angle":
            pots.minTSScanSize = round(
                x
                / (
                    float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                    + float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                ),
                4,
            )
        elif pots.scanType == "dihedral":
            pots.minTSScanSize = round(
                x
                / (
                    float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                    + float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                ),
                4,
            )
    elif pots.minTSScanVar.lower().strip() == "none":
        x = 50
        if pots.scanType == "bond":
            pots.minTSScanSize = round(
                x
                / 100.0
                / (
                    float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                    + float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                ),
                4,
            )
        elif pots.scanType == "angle":
            pots.minTSScanSize = round(
                x
                / (
                    float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                    + float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                ),
                4,
            )
        elif pots.scanType == "dihedral":
            pots.minTSScanSize = round(
                x
                / (
                    float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                    + float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                ),
                4,
            )
    elif pots.minTSScanVar.lower().strip() == "tight":
        x = 25
        if pots.scanType == "bond":
            pots.minTSScanSize = round(
                x
                / 100.0
                / (
                    float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                    + float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                ),
                4,
            )
        elif pots.scanType == "angle":
            pots.minTSScanSize = round(
                x
                / (
                    float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                    + float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                ),
                4,
            )
        elif pots.scanType == "dihedral":
            pots.minTSScanSize = round(
                x
                / (
                    float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                    + float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                ),
                4,
            )
    elif pots.minTSScanVar.lower().strip() == "verytight":
        x = 5
        if pots.scanType == "bond":
            pots.minTSScanSize = round(
                x
                / 100.0
                / (
                    float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                    + float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                ),
                4,
            )
        elif pots.scanType == "angle":
            pots.minTSScanSize = round(
                x
                / (
                    float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                    + float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                ),
                4,
            )
        elif pots.scanType == "dihedral":
            pots.minTSScanSize = round(
                x
                / (
                    float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                    + float(
                        etaatom.get_element_mass(
                            etaatom.get_element_num(
                                pots.finalGeomAtom[pots.atomList[0] + 1].e
                            )
                        )
                    )
                ),
                4,
            )

    # elif pots.minTSScanVar.lower() == 'var':
    #  if pots.scanType == "bond":
    #    pots.minTSScanSize = 0.001
    #  elif pots.scanType == "angle":
    #    pots.minTSScanSize = 0.1
    #  elif pots.scanType == "dihedral":
    #    pots.minTSScanSize = 0.1
    # Safe gaurd in case the user mispelled or mistyped the minTSScanVar keyword

    if pots.minTSScanSize == 0.0:
        print_and_write(
            [
                "You must have mispelled the #MINTSSCANSIZE keyword",
                "Your input was: " + pots.minTSScanVar,
            ],
            pots.verbose,
            pots.joblog,
        )

        #                     "The options are:",
        #                     "loose = ",
        #                     "none = ",
        #                     "tight = ",
        #                     "verytight = ",
        #                     "We will now be using: superveryubertight"], pots.verbose, pots.joblog)

        if pots.scanType == "bond":
            pots.minTSScanSize = 0.001
        elif pots.scanType == "angle":
            pots.minTSScanSize = 0.1
        elif pots.scanType == "dihedral":
            pots.minTSScanSize = 0.1
    return pots


### --- Logic for how the step size and geometry to be used should be. --- ###


def scan_ts_logic(pots):
    if pots.debug >= 1:
        print_and_write(["scan_ts_logic"], pots.verbose, pots.joblog)

    # energyList = []
    # ##Change the step number in the scan text file

    pots.scanStepNow += 1
    if pots.scanStepNow >= pots.scanLength:
        print_and_write(
            [
                "The scan has completed "
                + str(pots.scanStepNow)
                + " of "
                + str(pots.scanLength)
                + " steps.",
                "The job will stop now.",
            ],
            pots.verbose,
            pots.joblog,
        )
        pots.scan = 0
        pots.kill = True
    elif pots.scan == 1:

        # Start Logic for SCAN
        # Check for the scan energy file and open if it exists

        if os.path.isfile(pots.scratchdir + "/scan-" + pots.jobname + "-energy.txt"):
            the_file = open(
                pots.scratchdir + "/scan-" + pots.jobname + "-energy.txt", "a"
            )
            the_file.write(str(float(pots.finalenergy)) + "\n")
            the_file.close()
            if pots.test == 1:
                the_file = open(
                    pots.scratchdir + "/scan-" + pots.jobname + "-energy.txt", "r"
                )
                for line in the_file:
                    pots.scanEnergy.append(float(line.strip()))
                the_file.close()
            else:

                # print pots.scanEnergy
                # print float(pots.finalenergy) - pots.scanEnergy[-2]
                # Save the final energy to the scan energies list

                pots.scanEnergy.append(float(pots.finalenergy))

            # Print the energy

            print_and_write(
                [
                    "Last Energy: "
                    + str(pots.scanEnergy[-2])
                    + " Difference: "
                    + str(float(pots.finalenergy) - float(pots.scanEnergy[-2]))
                ],
                pots.verbose,
                pots.joblog,
            )

            # #Grab the highest energy from the list
            # if os.path.isfile(pots.scratchdir + "/HIGHEST-" + pots.jobname + ".xyz"):
            #  the_file = open(pots.scratchdir + "/HIGHEST-" + pots.jobname + ".xyz", 'r')
            #  for line in the_file:
            #    if "Energy" in line:
            #      var = line.strip().split()
            #      higheste = var[-1]
            # else:
            #  higheste = pots.scanEnergy[-1]

            # ##Start logic for TS-SCAN

            if pots.ts >= 1:
                pots.tsscanEnergy.append(float(pots.finalenergy))
            if pots.ts == 1:  # and pots.waitSize < pots.scanStepNow:

                # Check for a previous minor TS

                for i in range(len(pots.scanEnergy) - 2):
                    if (
                        pots.scanEnergy[i] > pots.scanEnergy[i + 1]
                        and pots.scanEnergy[i + 1] < pots.scanEnergy[i + 2]
                        and pots.minorTS == 0
                    ):
                        pots.minorTS = 1
                print_and_write(["Running TS-SCAN job."], pots.verbose, pots.joblog)

                # if pots.breakpoint != 0:
                #  print "YEP NEED TO FIX THIS!!"

                # Start logic for TS-SCAN

                if (
                    pots.ts == 1
                    and pots.waitSize < pots.scanStepNow
                    and pots.breakpoint == 0
                ):

                    # Did this energy equal the last one?

                    if pots.tsscanEnergy[-1] == pots.tsscanEnergy[-2]:
                        print_and_write(
                            [
                                "The energy did not change.",
                                "The job will now be terminated.",
                            ],
                            pots.verbose,
                            pots.joblog,
                        )
                        kill = True
                        scan = 0
                    elif len(pots.scanEnergy) == 2:

                        # ## --- For jobs that have run 2 steps
                        # Step ONE in the TS-SCAN.
                        # Did we go down in energy?

                        if (
                            pots.tsscanEnergy[-1] < pots.tsscanEnergy[-2]
                        ):  # Did the first step go the correct direction?
                            print_and_write(
                                [
                                    "The energy decreased and we assume we are going the wrong way.",
                                    "We will now turn around from the first geometry.",
                                ],
                                pots.verbose,
                                pots.joblog,
                            )

                            # pots = grab_previous_geom(pots, 1)  #Grab the first geometry (or second back from current)

                            pots = grab_highest_geom(pots)
                            pots.scanSize = -pots.scanSize
                        elif pots.tsscanEnergy[-1] > pots.tsscanEnergy[-2]:

                            # Did we go up in energy?

                            print_and_write(
                                [
                                    "The energy is increasing and we assume we are approaching the TS geometry.",
                                    "No changes to the step size will be made. Carry on.",
                                ],
                                pots.verbose,
                                pots.joblog,
                            )
                            write_high_geom(pots)
                            pots = save_highest_geom(pots)
                    elif len(pots.scanEnergy) == 3:

                        # ## --- For jobs that have run 3 steps
                        # Step TWO in the TS-SCAN.
                        # Did we go down twice in a row?

                        if (
                            pots.tsscanEnergy[-1]
                            < pots.tsscanEnergy[-2]
                            < pots.tsscanEnergy[-3]
                        ):
                            print_and_write(
                                [
                                    "We turned around and that didn't work.",
                                    "Most likely we were close but took to far of an initial step size.",
                                    "We will now turn back around and take a smaller step size.",
                                ],
                                pots.verbose,
                                pots.joblog,
                            )

                            # pots = grab_previous_geom(pots, 2)  #Grab the first geom (or third one back from the current)

                            pots = grab_highest_geom(pots)
                            pots.scanSize = -pots.scanSize / 3
                        elif (
                            pots.tsscanEnergy[-1] < pots.tsscanEnergy[-2]
                            and pots.tsscanEnergy[-2] > pots.tsscanEnergy[-3]
                            and (pots.minorTS == 0 or pots.minorTS == 2)
                        ):

                            # Did we go down after going up and we are allowing minor TSs?

                            print_and_write(
                                [
                                    "The energy decreased from the last step.",
                                    "One more step to see if we found a minor TS.",
                                ],
                                pots.verbose,
                                pots.joblog,
                            )
                        elif (
                            pots.tsscanEnergy[-1] < pots.tsscanEnergy[-2]
                            and pots.tsscanEnergy[-2] > pots.tsscanEnergy[-3]
                            and (pots.minorTS == 1 or pots.minorTS == 3)
                        ):

                            # Did we go down after going up and we don't want to allow minor TSs?

                            print_and_write(
                                [
                                    "The energy decreased from the last step.",
                                    "Going back to a more reasonable geometry, taking a smaller step size, and turning around.",
                                ],
                                pots.verbose,
                                pots.joblog,
                            )
                            pots = grab_highest_geom(pots)
                            pots.scanSize = -pots.scanSize / 3
                            del pots.tsscanEnergy[-1]
                        elif pots.tsscanEnergy[-1] > pots.tsscanEnergy[-2]:

                            # Did we go up at all?

                            print_and_write(
                                [
                                    "The energy is increasing and we assume we are approaching the TS geometry.",
                                    "No changes to the step size will be made. Carry on.",
                                ],
                                pots.verbose,
                                pots.joblog,
                            )
                            write_high_geom(pots)
                            pots = save_highest_geom(pots)
                    elif len(pots.scanEnergy) > 3:

                        # ## WE WILL NOW OUTPUT THE GEOMETRY TO THE SCAN-SCRATCH FOLDER BELOW
                        # ## --- For jobs that have run more than 3 steps
                        # Did we turn around twice and it was all down hill from the initial geom? If so turn around and take a smaller step size

                        if (
                            pots.tsscanEnergy[-1]
                            < pots.tsscanEnergy[-2]
                            < pots.tsscanEnergy[-3]
                            < pots.tsscanEnergy[-4]
                        ):
                            print_and_write(
                                [
                                    "Wow that was a bad idea.",
                                    "Going back to a more reasonable geometry, taking a smaller step size, and turning around.",
                                ],
                                pots.verbose,
                                pots.joblog,
                            )

                            # pots = grab_previous_geom(pots, 3)

                            pots = grab_highest_geom(pots)
                            pots.scanSize = -pots.scanSize / 3
                            del pots.tsscanEnergy[-1]
                            del pots.tsscanEnergy[-1]
                            del pots.tsscanEnergy[-1]
                        elif (
                            pots.tsscanEnergy[-1]
                            < pots.tsscanEnergy[-2]
                            < pots.tsscanEnergy[-3]
                            > pots.tsscanEnergy[-4]
                        ):

                            # Did we drop in energy twice in a row and were we comming from a high energy state?
                            # Better go back to the high energy state

                            print_and_write(
                                [
                                    "That extra step forward to see if we found a minor TS was a bad idea.",
                                    "Going back to a more reasonable geometry, taking a smaller step size, and turning around.",
                                ],
                                pots.verbose,
                                pots.joblog,
                            )

                            # pots = grab_previous_geom(pots, 2)

                            pots = grab_highest_geom(pots)
                            pots.scanSize = -pots.scanSize / 3
                            del pots.tsscanEnergy[-1]
                            del pots.tsscanEnergy[-1]
                        elif pots.tsscanEnergy[-1] < pots.tsscanEnergy[
                            -2
                        ] > pots.tsscanEnergy[-3] and (
                            pots.minorTS == 0 or pots.minorTS == 2
                        ):

                            # Did we drop in energy but it is the first time in recent history? If so check to see if it is a minor TS (if we havent found one yet already)

                            print_and_write(
                                [
                                    "The energy decreased from the last step.",
                                    "One more step to see if we found a minor TS.",
                                ],
                                pots.verbose,
                                pots.joblog,
                            )
                        elif (
                            pots.tsscanEnergy[-1]
                            < pots.tsscanEnergy[-2]
                            > pots.tsscanEnergy[-3]
                            and pots.minorTS == 3
                        ):

                            # Did we drop in energy but it is the first time in recent history? If so check to see if it is a minor TS (if we havent found one yet already)

                            print_and_write(
                                [
                                    "You have specified that you wish not to check for minor TSs.",
                                    "Going back to a more reasonable geometry, taking a smaller step size, and turning around.",
                                ],
                                pots.verbose,
                                pots.joblog,
                            )

                            # pots = grab_previous_geom(pots, 1)

                            pots = grab_highest_geom(pots)
                            pots.scanSize = -pots.scanSize / 3
                            del pots.tsscanEnergy[-1]
                        elif (
                            pots.tsscanEnergy[-1]
                            < pots.tsscanEnergy[-2]
                            > pots.tsscanEnergy[-3]
                            and pots.minorTS == 1
                        ):

                            # Did we drop in energy but it is the first time in recent history? AND did we already find a minor TS? Better go back and try again

                            print_and_write(
                                [
                                    "The energy is decreasing and we have already looked for a minor TS.",
                                    "Going back to a more reasonable geometry, taking a smaller step size, and turning around.",
                                ],
                                pots.verbose,
                                pots.joblog,
                            )

                            # pots = grab_previous_geom(pots, 1)

                            pots = grab_highest_geom(pots)
                            pots.scanSize = -pots.scanSize / 3
                            del pots.tsscanEnergy[-1]
                        elif (
                            pots.tsscanEnergy[-1]
                            > pots.tsscanEnergy[-2]
                            < pots.tsscanEnergy[-3]
                            < pots.tsscanEnergy[-4]
                        ):

                            # Did we increase in energy after going back down twice? Must be comming of a geometry change
                            # If this is the higher that the last high then save as high otherwise grab the highest geom

                            if (
                                pots.tsscanEnergy[-1]
                                > pots.highestGeomAtom[1].split()[1]
                            ):
                                print_and_write(
                                    [
                                        "The energy is increasing and we assume we are approaching the TS geometry.",
                                        "No changes to the step size will be made. Carry on.",
                                    ],
                                    pots.verbose,
                                    pots.joblog,
                                )
                                write_high_geom(pots)
                                pots = save_highest_geom(pots)
                                del pots.tsscanEnergy[-2]
                                del pots.tsscanEnergy[-2]
                            else:
                                print_and_write(
                                    [
                                        "The energy increased but it isn't the best geometry.",
                                        "Going back to a more reasonable geometry, taking a smaller step size, and turning around.",
                                    ],
                                    pots.verbose,
                                    pots.joblog,
                                )
                                pots = grab_highest_geom(pots)
                                pots.scanSize = -pots.scanSize / 3
                                del pots.tsscanEnergy[-1]
                                del pots.tsscanEnergy[-1]
                                del pots.tsscanEnergy[-1]
                        elif (
                            pots.tsscanEnergy[-1] > pots.tsscanEnergy[-2]
                            and pots.tsscanEnergy[-2] < pots.tsscanEnergy[-3]
                            and pots.tsscanEnergy[-3] > pots.tsscanEnergy[-4]
                            and (pots.minorTS == 0 or pots.minorTS == 2)
                        ):

                            # Did we increase in energy after going back down

                            print_and_write(
                                [
                                    "We found a minor TS and are saving this geometry as the highest energy.",
                                    "No changes to the step size will be made. Carry on.",
                                ],
                                pots.verbose,
                                pots.joblog,
                            )

                            # pots = grab_previous_geom(pots, 2)

                            write_high_geom(pots)
                            pots = save_highest_geom(pots)
                            del pots.tsscanEnergy[-2]
                        elif pots.tsscanEnergy[-1] > pots.tsscanEnergy[
                            -2
                        ] < pots.tsscanEnergy[-3] > pots.tsscanEnergy[-4] and (
                            pots.minorTS == 1 or pots.minorTS == 3
                        ):

                            if (
                                pots.tsscanEnergy[-1]
                                > pots.highestGeomAtom[1].split()[1]
                                and pots.finalGeomAtom[1].split()[-1]
                                - pots.highestGeomAtom[1].split()[-1]
                                - 2 * pots.scanSize
                                < 0.0001
                            ):
                                print_and_write(
                                    [
                                        "The energy is increasing and we found a better geometry.",
                                        "No changes to the step size will be made. Carry on.",
                                    ],
                                    pots.verbose,
                                    pots.joblog,
                                )
                                write_high_geom(pots)
                                pots = save_highest_geom(pots)
                                del pots.tsscanEnergy[-2]
                            else:
                                print_and_write(
                                    [
                                        "The energy increased but it isn't the best geometry.",
                                        "Going back to a more reasonable geometry, taking a smaller step size, and turning around.",
                                    ],
                                    pots.verbose,
                                    pots.joblog,
                                )
                                pots = grab_highest_geom(pots)
                                pots.scanSize = -pots.scanSize / 3
                                del pots.tsscanEnergy[-1]
                                del pots.tsscanEnergy[-1]
                        elif (
                            pots.tsscanEnergy[-1]
                            > pots.tsscanEnergy[-2]
                            > pots.tsscanEnergy[-3]
                            > pots.tsscanEnergy[-4]
                        ):

                            # Did we increase in energy twice in a row? Carry on.

                            print_and_write(
                                [
                                    "The energy is increasing and we assume we are approaching the TS geometry.",
                                    "No changes to the step size will be made. Carry on.",
                                ],
                                pots.verbose,
                                pots.joblog,
                            )
                            write_high_geom(pots)
                            pots = save_highest_geom(pots)

                        # ## WE WILL NOW OUTPUT THE GEOMETRY TO THE SCAN-SCRATCH FOLDER BELOW
                        # If the step size is below the minimum but above half of the minimum then make it the minimum and continue
                        # if pots.scanEnergy[-1] < pots.scanEnergy[-2]:

                        if (
                            pots.minTSScanSize
                            >= abs(pots.scanSize)
                            > pots.minTSScanSize / 2
                        ):
                            if pots.scanSize >= 0:
                                pots.scanSize = pots.minTSScanSize
                            else:
                                pots.scanSize = -pots.minTSScanSize
                        elif abs(pots.scanSize) < pots.minTSScanSize:

                            # If the step size is below the minimum and below or equal to half of the minimum then quit

                            print_and_write(
                                [
                                    "Scan step size has dropped below the threshold.",
                                    "This scan will now be terminated.",
                                ],
                                pots.verbose,
                                pots.joblog,
                            )
                            pots.scan = 0
                            pots.kill = True
                elif pots.ts == 1 and pots.waitSize > pots.scanStepNow:

                    # If the wait number is greater than the step now then wait

                    print_and_write(
                        [
                            "Waiting "
                            + str(pots.waitSize - pots.scanStepNow)
                            + " steps till TS-SCAN starts."
                        ],
                        pots.verbose,
                        pots.joblog,
                    )
                elif pots.ts == 1 and pots.waitSize == pots.scanStepNow:

                    # If the wait number is equal to the step now then stop waiting next time

                    print_and_write(
                        ["TS-SCAN will start next step."], pots.verbose, pots.joblog
                    )
            elif pots.ts == 2:

                # If dynamic is on then check if the energy is increasing and if so turn of dynamic and turn on ts

                print_and_write(
                    [
                        "Dynamic TS is turned on. The direction and step size will not change till the energy increases."
                    ],
                    pots.verbose,
                    pots.joblog,
                )
                if pots.tsscanEnergy[-1] > pots.tsscanEnergy[-2]:
                    print_and_write(
                        [
                            "The energy is increasing and we assume we are approaching the TS geometry.",
                            "Dynamic TS will now be turned off.",
                            "No changes to the step size will be made. Carry on.",
                        ],
                        pots.verbose,
                        pots.joblog,
                    )
                    pots.ts = 1

            # If this is step is higher in energy than the last then export it to the HIGHEST- file

            if (
                pots.scanEnergy[-1] > pots.scanEnergy[-2]
                and pots.breakpoint == 0
                and pots.ts == 0
            ):
                higheste = pots.scanEnergy[-1]
                if pots.highestGeomAtom:
                    pots = save_highest_geom(pots)
                elif higheste > pots.highestGeomAtom[1][-1]:
                    print_and_write(
                        ["This is the highest energy geometry as of yet."],
                        pots.verbose,
                        pots.joblog,
                    )
                    write_high_geom(pots)
                    pots = save_highest_geom(pots)
        else:

            # notHighest = 0
            # for line in pots.scanEnergy:
            #  if higheste < line:
            #    notHighest = 1
            # if notHighest == 0:
            #  print_and_write(["This is the highest energy geometry as of yet."], pots.verbose, pots.joblog)
            #  write_high_geom(pots)
            #  save_highest_geom(pots)
            # If said scan energy file does not exist make one and start adding to it

            print_and_write(
                [
                    "There is no known energy data about previous steps.",
                    "Creating this file now.",
                ],
                pots.verbose,
                pots.joblog,
            )
            the_file = open(
                pots.scratchdir + "/scan-" + pots.jobname + "-energy.txt", "w+"
            )
            the_file.write(str(pots.finalenergy) + "\n")
            the_file.close()

            # Save the final energy to the scan energies list

            pots.scanEnergy.append(float(pots.finalenergy))
            print_and_write(
                ["Energy :" + str(pots.finalenergy)], pots.verbose, pots.joblog
            )

            # Export the first optimized geometry

            pots = save_highest_geom(pots)
            if pots.ts >= 1:
                write_high_geom(pots)
                pots.tsscanEnergy.append(float(pots.finalenergy))

        # ################################################################################
        # ##Change the geometry for the next step

        if pots.debug >= 1:
            print_and_write(["bond_angle_dihedral"], pots.verbose, pots.joblog)

        # Move the atoms of interest around

        pots = renumber_fix(pots)

        # Using internal convert the xyz into a z-matrix

        pots = get_zmat(pots)
        if pots.verbose == 1:
            print_zmat(pots.ZMATAtom)
        if pots.scanType == "bond":

            # Grap the bond length now

            bondnow = pots.ZMATAtom[1].dist  # zmatlol[1][2]
            print_and_write(
                [
                    "Bond length was: "
                    + str(bondnow)
                    + ", it will be changed by "
                    + str(pots.scanSize)
                ],
                pots.verbose,
                pots.joblog,
            )

            # Change the bond length acording to the .tc file then replace the old bond length for this new one

            bondmath = bondnow + pots.scanSize
            print_and_write(
                ["Bond length now is: " + str(bondmath)], pots.verbose, pots.joblog
            )
            pots.ZMATAtom[1].dist = bondmath  # zmatlol[1][2]
            if bondmath < 0.5:
                pots.scan = 0
                print_and_write(
                    [
                        "The bond of interest crossed 0.5 angstroms.",
                        "This is unreasonable and the job will be stopped.",
                        "Are you doing nuclear chemistry?",
                    ],
                    pots.verbose,
                    pots.joblog,
                )
        elif pots.scanType == "angle":

            # Grap the bond length now

            anglenow = pots.ZMATAtom[2].angle  # zmatlol[2][4]
            print_and_write(
                [
                    "Angle was: "
                    + str(anglenow)
                    + ", it will be changed by "
                    + str(pots.scanSize)
                ],
                pots.verbose,
                pots.joblog,
            )

            # Change the bond length acording to the .tc file then replace the old bond length for this new one

            anglemath = anglenow + pots.scanSize
            print_and_write(
                ["Angle now is: " + str(anglemath)], pots.verbose, pots.joblog
            )
            pots.ZMATAtom[2].angle = anglemath  # zmatlol[2][4]

            # If the new angle will exceed 180 degrees then the scan size will be reversed for the next step.

            if (
                anglenow <= 180 <= anglemath
            ):  # or (anglenow >= -180 and anglemath <= -180):
                pots.scanSize = -pots.scanSize
                print_and_write(
                    [
                        "The angle of interest crossed 180 degrees and the scan step size has been reversed.",
                        "This will allow for the accurate continuation of this scan.",
                        "Original step size:\t" + str(-pots.scanSize),
                        "New step size:\t" + str(pots.scanSize),
                    ],
                    pots.verbose,
                    pots.joblog,
                )

            # If the new angle will exceed 0 degrees then the scan size will be reversed for the next step.

            if anglenow >= 0 >= anglemath:
                pots.scanSize = -pots.scanSize
                print_and_write(
                    [
                        "The angle of interest crossed 0 degrees and the scan step size has been reversed.",
                        "This will allow for the accurate continuation of this scan.",
                        "Original step size:\t" + str(-pots.scanSize),
                        "New step size:\t" + str(pots.scanSize),
                    ],
                    pots.verbose,
                    pots.joblog,
                )
        elif pots.scanType == "dihedral":

            # Grap the bond length now

            dihedralnow = pots.ZMATAtom[3].dihedral  # zmatlol[3][6]
            print_and_write(
                [
                    "Dihedral was: "
                    + str(dihedralnow)
                    + ", it will be changed by "
                    + str(pots.scanSize)
                ],
                pots.verbose,
                pots.joblog,
            )

            # Change the bond length acording to the .tc file then replace the old bond length for this new one

            dihedralmath = dihedralnow + pots.scanSize
            print_and_write(
                ["Dihedral now is: " + str(dihedralmath)], pots.verbose, pots.joblog
            )
            pots.ZMATAtom[3].dihedral = dihedralmath  # zmatlol[3][6]
        if pots.verbose == 1:
            print_zmat(pots.ZMATAtom)

        # Convert the augmented z-matrix to a xyz

        pots = get_xyz_from_zmat(pots)

        # Move the atoms back to their original locations

        pots = renumber_unfix(pots)

        # if pots.debug >= 2:
        # ....outputFixUnfixGeom(unfixgeom)

        # Open the scanfile and read it into a list

        the_file = open(pots.scanfile, "r")
        scanfilelist = the_file.readlines()
        the_file.close()

        # Make sure that the variable names and values are separated correctly

        for (i, line) in enumerate(scanfilelist):
            if "=" in line:
                line2 = line.replace("=", " = ")
                scanfilelist[i] = line2

        # Open the scanfile for writing it back in with changes

        f = open(pots.scanfile, "w+")
        for line in scanfilelist:
            if line.rstrip().split()[0] == "#SCAN" and pots.scan == 0:
                f.write("#SCAN = NO - DONE\n")
            elif line.rstrip().split()[0] == "#TS" and pots.ts == 0:
                f.write("#TS = NO\n")
            elif line.rstrip().split()[0] == "#SCANLENGTH" and pots.scan == 0:
                f.write("#SCANLENGTH = 0\n")
            elif line.rstrip().split()[0] == "#STEPSIZE":
                f.write("#STEPSIZE = " + str(pots.scanSize) + "\n")
            elif line.rstrip().split()[0] == "#SCANSTEPNOW":
                f.write("#SCANSTEPNOW = " + str(pots.scanStepNow) + "\n")
            else:
                f.write(line.replace("  =  ", " = "))
        print_and_write(
            [
                "Steps:\t\t" + str(pots.scanLength),
                "Step number:\t" + str(pots.scanStepNow),
                "Wait number:\t" + str(pots.waitSize),
                "Step size:\t" + str(pots.scanSize),
            ],
            pots.verbose,
            pots.joblog,
        )
        if pots.ts >= 1:
            print_and_write(
                ["Minimum step size:\t+/- " + str(pots.minTSScanSize)],
                pots.verbose,
                pots.joblog,
            )
        if pots.scan == 0:
            pots.kill = True
    return pots


###Depreciated function. This is now a part of scan_ts_logic
### --- Changes the geometries by converting them to zmatrix making changes to the file them converting them back to xyz --- ###


def bond_angle_dihedral(pots):
    if pots.debug >= 1:
        print_and_write(["bond_angle_dihedral"], pots.verbose, pots.joblog)

    # Move the atoms of interest around

    pots = renumber_fix(pots)

    # Using internal convert the xyz into a z-matrix

    pots = get_zmat(pots)
    if pots.verbose == 1:
        print_zmat(pots.ZMATAtom)
    if pots.scanType == "bond":

        # Grap the bond length now

        bondnow = pots.ZMATAtom[1].dist  # zmatlol[1][2]
        print_and_write(
            [
                "Bond length was: " + str(bondnow),
                "It will be changed by " + str(pots.scanSize),
            ],
            pots.verbose,
            pots.joblog,
        )

        # Change the bond length acording to the .tc file then replace the old bond length for this new one

        bondmath = bondnow + pots.scanSize
        print_and_write(
            ["Bond length now is: " + str(bondmath)], pots.verbose, pots.joblog
        )
        pots.ZMATAtom[1].dist = bondmath  # zmatlol[1][2]
        if bondmath < 0:
            pots.scan = 0
            print_and_write(
                [
                    "The bond of interest crossed 0 angstroms.",
                    "This is unreasonable and the job will be stopped.",
                ],
                pots.verbose,
                pots.joblog,
            )
    elif pots.scanType == "angle":

        # Grap the bond length now

        anglenow = pots.ZMATAtom[2].angle  # zmatlol[2][4]
        print_and_write(
            [
                "Angle was: " + str(anglenow),
                "It will be changed by " + str(pots.scanSize),
            ],
            pots.verbose,
            pots.joblog,
        )

        # Change the bond length acording to the .tc file then replace the old bond length for this new one

        anglemath = anglenow + pots.scanSize
        print_and_write(["Angle now is: " + str(anglemath)], pots.verbose, pots.joblog)
        pots.ZMATAtom[2].angle = anglemath  # zmatlol[2][4]

        # If the new angle will exceed 180 degrees then the scan size will be reversed for the next step.

        if anglenow <= 180 <= anglemath:  # or (anglenow >= -180 and anglemath <= -180):
            pots.scanSize = -pots.scanSize
            print_and_write(
                [
                    "The angle of interest crossed 180 degrees and the scan step size has been reversed.",
                    "This will allow for the accurate continuation of this scan.",
                ],
                pots.verbose,
                pots.joblog,
            )

        # If the new angle will exceed 0 degrees then the scan size will be reversed for the next step.

        if anglenow >= 0 >= anglemath:
            pots.scanSize = -pots.scanSize
            print_and_write(
                [
                    "The angle of interest crossed 0 degrees and the scan step size has been reversed.",
                    "This will allow for the accurate continuation of this scan.",
                ],
                pots.verbose,
                pots.joblog,
            )
    elif pots.scanType == "dihedral":

        # Grap the bond length now

        dihedralnow = pots.ZMATAtom[3].dihedral  # zmatlol[3][6]
        print_and_write(
            [
                "Dihedral was: " + str(dihedralnow),
                "It will be changed by " + str(pots.scanSize),
            ],
            pots.verbose,
            pots.joblog,
        )

        # Change the bond length acording to the .tc file then replace the old bond length for this new one

        dihedralmath = dihedralnow + pots.scanSize
        print_and_write(
            ["Dihedral now is: " + str(dihedralmath)], pots.verbose, pots.joblog
        )
        pots.ZMATAtom[3].dihedral = dihedralmath  # zmatlol[3][6]
    if pots.verbose == 1:
        print_zmat(pots.ZMATAtom)

    # Using internal convert the z-matrix to a xyz

    pots = get_xyz_from_zmat(pots)

    # Move the atoms back to their original locations

    pots = renumber_unfix(pots)

    # if pots.debug >= 2:
    # ....outputFixUnfixGeom(unfixgeom)

    return pots


### --- Convert coordinates to distance/angle/dihedral --- ###


def get_zmat(pots):
    if pots.debug >= 1:
        print_and_write(["get_zmat"], pots.verbose, pots.joblog)

    # Make the ZMATAtom only as long as the atoms in the finalGeomAtom list

    pots.ZMATAtom = [0] * (len(pots.finalGeomAtom) - 2)

    # Iterate only throught the atoms and not the first two lines

    for i in range(2, len(pots.finalGeomAtom)):

        # Make each ZMATAtom entry the class ZMAT_Atom and give it the elemental information

        pots.ZMATAtom[i - 2] = etaatom.AtomZMAT()
        pots.ZMATAtom[i - 2].e = pots.finalGeomAtom[i].e

        # If it is after the first atom give it the distance from the last atom

        if i > 2:
            pots.ZMATAtom[i - 2].atomdist = i - 3
            pots.ZMATAtom[i - 2].dist = etanumpy.get_distance(
                i - 1, i, pots.finalGeomAtom
            )

        # If it is at least the third atom give it the angle from the last two atoms

        if i > 3:
            pots.ZMATAtom[i - 2].atomangle = i - 4
            pots.ZMATAtom[i - 2].angle = etanumpy.get_angle(
                i - 2, i - 1, i, pots.finalGeomAtom
            )

        # If it is at least the fourth atom give it the dihedral from the last three atoms

        if i > 4:
            pots.ZMATAtom[i - 2].atomdihedral = i - 5
            pots.ZMATAtom[i - 2].dihedral = etanumpy.get_dihedral(
                i - 3, i - 2, i - 1, i, pots.finalGeomAtom
            )
    return pots


### --- Get the XYZ coordinates from distance, angle and dihedral data --- ###


def get_xyz_from_zmat(pots):
    if pots.debug >= 1:
        print_and_write(["get_xyz_from_zmat"], pots.verbose, pots.joblog)
    line = []
    x = 0
    y = 0
    z = 0
    xyzgeom = []
    for i in range(len(pots.ZMATAtom)):

        # ## Set up the variables to be used for the function
        # at4 = i

        if i > 0:
            dist = float(pots.ZMATAtom[i].dist)
            at3 = int(pots.ZMATAtom[i].atomdist)
        if i > 1:
            angle = float(pots.ZMATAtom[i].angle)
            anglerad = radians(angle)  # * math.pi / 180
            at2 = int(pots.ZMATAtom[i].atomangle)
        if i > 2:
            dihedral = float(pots.ZMATAtom[i].dihedral)
            dihedralrad = radians(dihedral)  # * math.pi / 180
            at1 = int(pots.ZMATAtom[i].atomdihedral)

        # ## Start to place the atoms in their locations
        # First atom get put at the origin

        if i == 0:
            x = 0.0
            y = 0.0
            z = 0.0
        elif i == 1:

            # Second atom gets put on the x axis

            x = dist
            y = 0.0
            z = 0.0
        elif i == 2:

            # Third atom is put on the xz plane with some trig

            a = xyzgeom[at3][1]
            b = dist
            x = a + dist * cos(math.pi - anglerad)
            y = 0.0
            z = -dist * sin(math.pi - anglerad)
        elif i >= 3:

            # Fourth and later atoms are placed by using vector math, of which I understand only slightly
            # ###The at4 x,y,z coordinates from spherical coord

            sx = dist * sin(anglerad) * cos(dihedralrad)
            sy = (
                -dist * sin(anglerad) * sin(dihedralrad)
            )  # For some reason I need to have a negative here to get the correct sign in the output..... weird
            sz = dist * cos(anglerad)
            at4l = [sx, sy, sz]

            # ##Finding the angle theta
            # Make the list of lists for the three point (z-axis, origin, and translated atom 2) needed for an angle calculation

            z32 = [
                [0, 0, 0, 1],
                [0, 0, 0, 0],
                [
                    0,
                    xyzgeom[at2][1] - xyzgeom[at3][1],
                    xyzgeom[at2][2] - xyzgeom[at3][2],
                    xyzgeom[at2][3] - xyzgeom[at3][3],
                ],
            ]

            # Get theta using the getangle function

            theta = radians(etanumpy.get_angle(0, 1, 2, z32))

            # ##Rodrigues' rotation formula
            # Create the vectprs needed to calculate k

            vector3 = array([xyzgeom[at3][1], xyzgeom[at3][2], xyzgeom[at3][3]])
            vector2 = array([xyzgeom[at2][1], xyzgeom[at2][2], xyzgeom[at2][3]])
            vector0 = array([0, 0, 1])

            # Calculate k for the Rodrigues rotation formula

            k = cross(vector2 - vector3, vector0) / linalg.norm(
                cross(vector2 - vector3, vector0)
            )

            # Generate an array for translated 1

            t1 = [
                xyzgeom[at1][1] - xyzgeom[at3][1],
                xyzgeom[at1][2] - xyzgeom[at3][2],
                xyzgeom[at1][3] - xyzgeom[at3][3],
            ]

            # Calculate the Rodrigues rotation matrix

            rr23t1 = (
                dot(t1, cos(theta))
                + dot(cross(k, t1), sin(theta))
                + dot(dot(k, dot(k, t1)), 1 - cos(theta))
            )

            # Make the list of lists for the four points (x-axis, z-axis, origin, and rotated translated 1) needed for a dihedral calculation

            xz31 = [
                [0, 1, 0, 0],
                [0, 0, 0, 1],
                [0, 0, 0, 0],
                [0, rr23t1[0], rr23t1[1], rr23t1[2]],
            ]

            # Get phi using the getdihedral function

            phi = radians(etanumpy.get_dihedral(0, 1, 2, 3, xz31))

            # ##Rotation matrix
            # Create the array for the rotation matrix including dihedral phi

            rm = array([[cos(phi), sin(phi), 0], [-sin(phi), cos(phi), 0], [0, 0, 1]])

            # Calculate the dot product of the rotation matrix and the coordinates for 4 (from spherical)

            rm4 = dot(rm, at4l)

            # Calculate the rotated coordinates of the rotated coordinates of atom 4

            rrn23rm4 = (
                dot(rm4, cos(-theta))
                + dot(cross(k, rm4), sin(-theta))
                + dot(dot(k, dot(k, rm4)), 1 - cos(-theta))
            )

            # Final coordinates that are rotated, rotated and translated

            x = rrn23rm4[0] + xyzgeom[at3][1]
            y = rrn23rm4[1] + xyzgeom[at3][2]
            z = rrn23rm4[2] + xyzgeom[at3][3]

        # Putting everything into a list to send back

        line.append(pots.ZMATAtom[i].e)
        line.append(x)
        line.append(y)
        line.append(z)
        xyzgeom.append(list(line))
        line[:] = []

    # Add the first two lines to the newgeom lol

    for i in range(len(xyzgeom)):

        # print_and_write([xyzgeom[i]], pots.verbose, pots.joblog)

        pots.finalGeomAtom[i + 2].e = str(xyzgeom[i][0])
        pots.finalGeomAtom[i + 2].x = float(xyzgeom[i][1])
        pots.finalGeomAtom[i + 2].y = float(xyzgeom[i][2])
        pots.finalGeomAtom[i + 2].z = float(xyzgeom[i][3])
    return pots


### --- Function to easily print out a geom using the atom class for debuging --- ###


def print_geom(geom):
    for i in range(len(geom)):
        if isinstance(geom[i], etaatom.Atom):
            line = (
                geom[i].e
                + "  "
                + str("{:.6f}".format(geom[i].x))
                + "  "
                + str("{:.6f}".format(geom[i].y))
                + "  "
                + str("{:.6f}".format(geom[i].z))
            )
        else:

            # else print the string (number of atoms or comment line)

            line = str(geom[i])
        print_and_write([line], pots.verbose, pots.joblog)
    return


### --- Function to easily print out a geom using the atom class for debuging --- ###


def print_zmat(zmat):
    for i in range(len(zmat)):
        if isinstance(zmat[i], etaatom.AtomZMAT):
            line = (
                zmat[i].e
                + "  "
                + str(zmat[i].atomdist)
                + "  "
                + str("{:.6f}".format(zmat[i].dist))
                + "  "
                + str(zmat[i].atomangle)
                + "  "
                + str("{:.6f}".format(zmat[i].angle))
                + "  "
                + str(zmat[i].atomdihedral)
                + "  "
                + str("{:.6f}".format(zmat[i].dihedral))
            )
        else:

            # else print the string (number of atoms or comment line)

            line = str(geom[i])
        print_and_write([line], pots.verbose, pots.joblog)
    return


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

### --- G0X FUNCTIONS --- ###
########################################################################################################################
### --- Parses the input G0X file --- ###


def parse_input_g0x(pots):
    if pots.debug >= 1:
        print_and_write(["parse_input_g0x"], pots.verbose, pots.joblog)
    pots.QMconfig = []
    pots.QMcomment = ""
    pots.charge = 0
    pots.multi = 1
    pots.initialGeomAtom = []
    initialgeom = []
    pots.QMendlines = []
    inputnumatoms = 0
    emptyline = 0
    endline = "no"
    the_file = open(pots.baseinput, "r")
    for line in the_file:
        if line.strip() == "" and endline == "no":
            emptyline += 1
        if line.strip() != "" and emptyline == 0:
            pots.QMconfig.append(line.strip())
        if line.strip() != "" and emptyline == 1:
            pots.QMcomment = line.strip()
        if line.strip() != "" and emptyline == 3:
            initialgeom.append(line.strip())
            inputnumatoms += 1
        if line.strip() != "" and emptyline == 2:
            pots.charge = int(line.split()[0])
            pots.multi = int(line.split()[1])
            emptyline += 1
        if line.strip() != "" and emptyline == 4:
            endline = "yes"
            emptyline += 1
        if emptyline >= 5:
            pots.QMendlines.append(line.strip())
    pots.numatoms = inputnumatoms
    pots.initialGeomAtom = [0] * (pots.numatoms + 2)
    pots.initialGeomAtom[0] = inputnumatoms
    pots.initialGeomAtom[1] = pots.QMcomment
    for (i, var) in enumerate(initialgeom):
        line = var.split()
        pots.initialGeomAtom[i + 2] = etaatom.Atom()
        pots.initialGeomAtom[i + 2].e = line[0]
        pots.initialGeomAtom[i + 2].en = etaatom.get_element_num(line[0])
        pots.initialGeomAtom[i + 2].x = float(line[1])
        pots.initialGeomAtom[i + 2].y = float(line[2])
        pots.initialGeomAtom[i + 2].z = float(line[3])

    # TODO make sure no other "opt"s get caught in this...
    # Set the maxcycles for opt to 25 if larger than 100

    rewrite = False
    for (i, line) in enumerate(pots.QMconfig):
        sub = line.split()
        for subline in sub:
            if "opt=(" in subline and "maxcycle" in subline:
                subsubline = subline.split(",")
                for var in subsubline:
                    if "maxcycle=" in var:
                        pots.maxopt = int(var.split("=")[-1].replace(")", ""))
                        if pots.maxopt > 100:
                            print_and_write(
                                [
                                    "The max number of optimization cycles was too low.",
                                    "It has now been set to 25.",
                                ],
                                pots.verbose,
                                pots.joblog,
                            )
                            pots.QMconfig[i] = line.replace(
                                "maxcycle=" + str(pots.maxopt), "maxcycle=25"
                            )
                            rewrite = True
            elif "opt=(" in subline and "maxcycle" not in subline:
                pots.QMconfig[i] = line.replace("opt=(", "opt=(maxcycle=25,")
                rewrite = True
            elif "opt" in subline:
                pots.QMconfig[i] = line.replace("opt", "opt=(maxcycle=25) ")
                rewrite = True
    if rewrite:
        pots.maxopt = 25
        pots.finalGeomAtom = [0] * len(pots.initialGeomAtom)
        for i in range(len(pots.initialGeomAtom)):
            if i <= 1:
                pots.finalGeomAtom[i] = pots.initialGeomAtom[i]
            elif i >= 2:
                pots.finalGeomAtom[i] = etaatom.Atom()
                pots.finalGeomAtom[i].e = pots.initialGeomAtom[i].e
                pots.finalGeomAtom[i].en = pots.initialGeomAtom[i].en
                pots.finalGeomAtom[i].x = pots.initialGeomAtom[i].x
                pots.finalGeomAtom[i].y = pots.initialGeomAtom[i].y
                pots.finalGeomAtom[i].z = pots.initialGeomAtom[i].z
        rebuild_g0x_com(pots)
        pots.finalGeomAtom = []
    return pots


### --- Parses the geometry and energy created by G0X and gives back the trajectory in a list of lists --- ###


def parse_output_g0x(pots):
    if pots.debug >= 1:
        print_and_write(["parse_output_g0x"], pots.verbose, pots.joblog)
    alist = []
    f_list = []
    geom = 0
    force = 0
    numatoms = 0
    pots.allEnergy = []
    pots.trajectoryGeomAtom = []
    the_file = open(pots.baselog, "r")
    for (idx, line) in enumerate(the_file):
        if line.strip() != "":
            if line.split()[0] == "SCF":
                pots.allEnergy.append(float(line.split()[4]))
            if line.split()[0] == "Maximum" and line.split()[1] == "Displacement":
                pots.finalDisp = float(line.split()[2])
            if line.split()[0] == "NAtoms=":
                pots.numatoms = int(line.split()[1])
            if "-----------------------------------------------" in line and geom >= 6:
                geom = 0
            if geom >= 1:
                geom += 1
            if line.strip() == "Standard orientation:":
                geom = 1
            if geom >= 5:
                alist.append(line)

            # Get forces

            if "-----------------------------------------------" in line and force >= 5:
                force = 0
            if force >= 1:
                force += 1
            if "Forces (Hartrees/Bohr)" in line:
                force = 1
            if force >= 4:
                f_list.append(line)
            elif "Multiplicity" in line:

                # Get charge and multiplicity

                chargemulti = line.split()
                pots.charge = chargemulti[2]
                pots.multi = chargemulti[-1]
    the_file.close()

    # Get the final energy from the energy list

    pots.finalenergy = float(pots.allEnergy[-1])
    alistdedup = [
        x
        for x in alist
        if "------------------------------------------------------" not in x
    ]

    # #####################################################################
    # Create the trajectory

    pots.trajectoryGeomAtom = [0] * (
        len(alistdedup) / pots.numatoms * (pots.numatoms + 2)
    )
    skip = 0
    for i in range(len(alistdedup)):
        if i % pots.numatoms == 0 and i / pots.numatoms != len(pots.allEnergy):
            pots.trajectoryGeomAtom[i + skip] = int(pots.numatoms)
            pots.trajectoryGeomAtom[i + skip + 1] = "Energy: " + str(
                pots.allEnergy[i / pots.numatoms]
            )
            skip += 2
        elif i % pots.numatoms == 0 and i / pots.numatoms >= len(pots.allEnergy):
            pots.trajectoryGeomAtom[i + skip] = int(pots.numatoms)
            pots.trajectoryGeomAtom[i + skip + 1] = "Energy: " + str(pots.allEnergy[-1])
            skip += 2
        line = alistdedup[i].split()
        pots.trajectoryGeomAtom[i + skip] = etaatom.Atom()
        pots.trajectoryGeomAtom[i + skip].e = etaatom.get_element_name(int(line[1]))
        pots.trajectoryGeomAtom[i + skip].en = int(line[0])
        pots.trajectoryGeomAtom[i + skip].x = float(line[3])
        pots.trajectoryGeomAtom[i + skip].y = float(line[4])
        pots.trajectoryGeomAtom[i + skip].z = float(line[5])

    # #####################################################################
    # Extract the final geometry from the trajectory

    pots.finalGeomAtom = [0] * (pots.trajectoryGeomAtom[0] + 2)
    for i in range(len(pots.finalGeomAtom)):
        j = len(pots.trajectoryGeomAtom) - len(pots.finalGeomAtom) + i
        if i <= 1:
            pots.finalGeomAtom[i] = pots.trajectoryGeomAtom[j]
        elif i >= 2:
            pots.finalGeomAtom[i] = etaatom.Atom()
            pots.finalGeomAtom[i].e = pots.trajectoryGeomAtom[j].e
            pots.finalGeomAtom[i].en = pots.trajectoryGeomAtom[j].en
            pots.finalGeomAtom[i].x = pots.trajectoryGeomAtom[j].x
            pots.finalGeomAtom[i].y = pots.trajectoryGeomAtom[j].y
            pots.finalGeomAtom[i].z = pots.trajectoryGeomAtom[j].z

    # #####################################################################
    # Extract the final forces

    pots.allForce = [0] * pots.numatoms
    for i in range(len(pots.allForce)):
        j = len(f_list) - len(pots.allForce) + i
        line = f_list[j].strip().split()
        pots.allForce[i] = line[2] + "  " + line[3] + "  " + line[4]
    return pots


### --- G0X Optimization Logic --- ###


def opt_logic_g0x(pots):
    if pots.debug >= 1:
        print_and_write(["opt_logic_g0x"], pots.verbose, pots.joblog)

    # print_and_write(["BEGIN OPT LOGIC!"
    # ## --- Logic for how the step size and geometry to be used should be. --- ###
    # Check to see if smaller step sizes are needed.

    print_and_write(["Checking for jitter in the energies."], pots.verbose, pots.joblog)
    jitter = 0
    f = open(pots.scratchdir + "/optenergies-" + pots.jobname + ".txt", "a")
    for (i, line) in enumerate(pots.allEnergy):
        f.write(str(line) + "\n")
        if i > len(pots.allEnergy) / 2:
            if (
                pots.allEnergy[i] < pots.allEnergy[i - 1]
                and pots.allEnergy[i - 1] > pots.allEnergy[i - 2]
                or pots.allEnergy[i] > pots.allEnergy[i - 1]
                and pots.allEnergy[i - 1] < pots.allEnergy[i - 2]
            ):
                jitter += 1
    if jitter > len(pots.allEnergy) / 4:
        print_and_write(
            ["Jitters in last half of OPT: " + str(jitter)], pots.verbose, pots.joblog
        )

        # + " > " + str(len(allEnergy)/4) + "\n")

        print_and_write(
            [
                "Due to the jitter in energies from the last OPT iteration,",
                "the step size of the optimization algorithm will be shortened.",
            ],
            pots.verbose,
            pots.joblog,
        )
        pots.shortenStepSize = 1

    # Check to see if the job ran too long.

    f = open(pots.scratchdir + "/optenergies-" + pots.jobname + ".txt", "r")
    fileenergylist = f.readlines()
    if len(fileenergylist) > pots.stepNumber:
        pots.kill = True
        print_and_write(
            ["This job will end due to: exhausted all optimization steps."],
            pots.verbose,
            pots.joblog,
        )
    f.close()

    # Check to see if the job finished optimization.

    count = 0
    f = open(pots.baselog, "r")
    for line in reversed(open(pots.baselog).readlines()):
        if "Normal termination" in line or "Optimization completed" in line:
            if pots.jobtype == "SCAN" or pots.jobtype == "TS-SCAN":
                pots.optimized = True
                print_and_write(
                    [
                        "This job will continue due to: 'Normal termination' and/or 'Optimization Complete'."
                    ],
                    pots.verbose,
                    pots.joblog,
                )
            else:
                pots.kill = True
                print_and_write(
                    [
                        "This job will end due to: 'Normal termination' and/or 'Optimization Complete'."
                    ],
                    pots.verbose,
                    pots.joblog,
                )
            break
        count += 1
        if count > 1000:
            break

    # Check if a 'break' has occured

    breaklogic = 0
    line = pots.QMcomment.split()
    for var in line:
        if pots.breaktype == "bond" and breaklogic == 4:
            pots.breakpoint = float(var)
            breaklogic = 0
        elif pots.breaktype == "angle" and breaklogic == 5:
            pots.breakpoint = float(var)
            breaklogic = 0
        elif pots.breaktype == "dihedral" and breaklogic == 6:
            pots.breakpoint = float(var)
            breaklogic = 0
        elif breaklogic >= 2:
            pots.breakatoms.append(int(var))
            breaklogic += 1
        elif breaklogic == 1:
            pots.breaktype = var
            breaklogic += 1
        if var == "Break:":
            breaklogic = 1
    if pots.breaktype != "":
        print_and_write(
            [
                "Break point logic implemented.",
                "Checking if "
                + str(pots.breaktype)
                + " "
                + str(pots.breakatoms)
                + " are within "
                + str(pots.breakpoint)
                + ".",
            ],
            pots.verbose,
            pots.joblog,
        )

    # Break point logic implementation

    if pots.breaktype == "bond":
        if pots.breakpoint < etanumpy.get_distance(
            pots.breakatoms[0] - 1, pots.breakatoms[1] - 1, pots.finalGeomAtom
        ):
            pots.kill = True
            print_and_write(
                ["This job will end due to: a break point failure."],
                pots.verbose,
                pots.joblog,
            )
    elif pots.breaktype == "angle":
        if pots.breakpoint < etanumpy.get_angle(
            pots.breakatoms[0] - 1,
            pots.breakatoms[1] - 1,
            pots.breakatoms[2] - 1,
            pots.finalGeomAtom,
        ):
            pots.kill = True
            print_and_write(
                ["This job will end due to: a break point failure."],
                pots.verbose,
                pots.joblog,
            )
    elif pots.breaktype == "dihedral":
        if pots.breakpoint < etanumpy.get_dihedral(
            pots.breakatoms[0] - 1,
            pots.breakatoms[1] - 1,
            pots.breakatoms[2] - 1,
            pots.breakatoms[3] - 1,
            pots.finalGeomAtom,
        ):
            pots.kill = True
            print_and_write(
                ["This job will end due to: a break point failure."],
                pots.verbose,
                pots.joblog,
            )

    # Stepsize logic implementation

    iop = 0
    newsub = 0
    if pots.shortenStepSize == 1 and not pots.kill and not pots.optimized:
        for (i, line) in enumerate(pots.QMconfig):
            if "iop(1/8=" in line:
                var = line.split()
                iop = 1
                for sub in var:
                    if "iop(1/8=" in sub:
                        print_and_write(
                            ["Your step size was: " + str(sub)],
                            pots.verbose,
                            pots.joblog,
                        )
                        newsub = (
                            "iop(1/8="
                            + str(
                                int(
                                    round(float(sub.split("=")[1].split(")")[0]) * 0.75)
                                )
                            )
                            + ")"
                        )
                        print_and_write(
                            ["It will now be: " + str(newsub)],
                            pots.verbose,
                            pots.joblog,
                        )
                        pots.QMconfig[i] = line.replace(sub, newsub)
        if iop == 0:
            for (i, line) in enumerate(pots.QMconfig):
                if "opt=(" in line:
                    print_and_write(
                        ["We are now implementing a smaller step size of iop(1/8=18)"],
                        pots.verbose,
                        pots.joblog,
                    )
                    pots.QMconfig[i] = line + " iop(1/8=18)"

    # Implement speedups for G09/G03

    if not pots.kill and not pots.optimized:
        for (i, line) in enumerate(pots.QMconfig):
            if "opt=(" in line and "guess=read" not in line:
                print_and_write(
                    ["We are now implementing guess=read for the next iteration."],
                    pots.verbose,
                    pots.joblog,
                )
                pots.QMconfig[i] = line + " guess=read"
        for (i, line) in enumerate(pots.QMconfig):
            if "opt=(" in line and "readfc" not in line:
                print_and_write(
                    ["We are now implementing readfc for the next iteration."],
                    pots.verbose,
                    pots.joblog,
                )
                var = line.split()
                for sub in var:
                    if "opt=(" in sub:
                        newsub = sub[:-1] + ",readfc)"
                        pots.QMconfig[i] = line.replace(sub, newsub)

    # Kill logic implementation
    # if pots.kill == True:
    #  pots.QMcomment = "!KILL! " + pots.QMcomment

    return pots


### --- Output the G09 com file for subsequent computation --- ###


def rebuild_g0x_com(pots):
    if pots.debug >= 1:
        print_and_write(["rebuild_g0x_com"], pots.verbose, pots.joblog)

    # Output the new geometry in G09 input format

    datetime = time.strftime("%Y.%m.%d-%H:%M:%S")

    # Backup .com file

    if os.path.exists(pots.baseinput):
        if os.path.isfile(pots.baseinput):
            shutil.copyfile(
                pots.baseinput,
                pots.scratchdir
                + "/Previous_Files/"
                + datetime
                + "-"
                + pots.jobname
                + ".com",
            )

    # Backup .chk file

    if os.path.exists(pots.jobname + ".chk"):
        if os.path.isfile(pots.jobname + ".chk"):
            shutil.copyfile(
                pots.jobname + ".chk",
                pots.scratchdir
                + "/Previous_Files/"
                + datetime
                + "-"
                + pots.jobname
                + ".chk",
            )

    # Backup .log file

    if os.path.exists(pots.baselog):
        if os.path.isfile(pots.baselog):
            shutil.copyfile(
                pots.baselog,
                pots.scratchdir
                + "/Previous_Files/"
                + datetime
                + "-"
                + pots.jobname
                + ".log",
            )
    f = open(pots.baseinput, "w+")

    # for line in beginText:

    for line in pots.QMconfig:
        f.write(line + "\n")
    f.write("\n")
    f.write(pots.QMcomment + "\n")
    f.write("\n")
    f.write(str(pots.charge) + "  " + str(pots.multi) + "\n")
    for i in range(len(pots.finalGeomAtom)):
        if i >= 2:
            f.write(
                pots.finalGeomAtom[i].e
                + "  "
                + str("{:.6f}".format(pots.finalGeomAtom[i].x))
                + "  "
                + str("{:.6f}".format(pots.finalGeomAtom[i].y))
                + "  "
                + str("{:.6f}".format(pots.finalGeomAtom[i].z))
                + "\n"
            )

    # for line in endText:

    f.write("\n")
    for line in pots.QMendlines:
        f.write(line + "\n")
    f.close()
    return


### --- JAGUAR FUNCTIONS --- ###
########################################################################################################################
### --- Parses the input G0X file --- ###


def parse_input_jaguar(pots):
    if pots.debug >= 1:
        print_and_write(["parse_input_jaguar"], pots.verbose, pots.joblog)
    pots.QMconfig = []
    pots.QMcomment = ""
    pots.charge = 0
    pots.multi = 1
    pots.initialGeomAtom = []
    pots.numatoms = 0
    initialgeom = []
    pots.QMendlines = []
    xyzline = 0
    the_file = open(pots.baseinput, "r")
    for line in the_file:
        if xyzline == 0 and "&zmat" in line:
            xyzline = 1
        elif xyzline == 0:
            pots.QMconfig.append(line)
        elif xyzline == 1 and "&" in line:
            xyzline = 2
        elif xyzline == 1:
            initialgeom.append(line)
            pots.numatoms += 1
        elif xyzline == 2:
            pots.QMendlines.append(line)
    pots.initialGeomAtom = [0] * (pots.numatoms + 2)
    pots.initialGeomAtom[0] = pots.numatoms
    pots.initialGeomAtom[1] = pots.QMcomment
    for (i, var) in enumerate(initialgeom):
        line = var.split()
        pots.initialGeomAtom[i + 2] = etaatom.Atom()
        pots.initialGeomAtom[i + 2].e = "".join([j for j in line[0] if not j.isdigit()])
        pots.initialGeomAtom[i + 2].en = etaatom.get_element_num(
            pots.initialGeomAtom[i + 2].e
        )
        pots.initialGeomAtom[i + 2].x = float(line[1])
        pots.initialGeomAtom[i + 2].y = float(line[2])
        pots.initialGeomAtom[i + 2].z = float(line[3])

    # Set the maxcycles for opt to 25 if larger than 100
    # for i, line in enumerate(pots.QMconfig):
    #  sub = line.split()
    #  for subline in sub:
    #    if "opt=" in subline:
    #      subsubline = subline.split(",")
    #      for var in subsubline:
    #        if "maxcycle=" in var:
    #          pots.maxopt = int(var.split('=')[-1].replace(')',''))
    #          if pots.maxopt > 100:
    #            print_and_write(["The max number of optimization cycles was too low.",\
    #                                "It has now been set to 25."], pots.verbose, pots.joblog)
    #            pots.QMconfig[i] = line.replace("maxcycle=" + str(pots.maxopt), "maxcycle=25")
    #            pots.finalGeomAtom = [0] * len(pots.initialGeomAtom)
    #            pots.maxopt = 25
    #            for i in range(len(pots.initialGeomAtom)):
    #              if i <= 1:
    #                pots.finalGeomAtom[i] = pots.initialGeomAtom[i]
    #              elif i >= 2:
    #                pots.finalGeomAtom[i] = etaatom.Atom()
    #                pots.finalGeomAtom[i].e             = pots.initialGeomAtom[i].e
    #                pots.finalGeomAtom[i].en            = pots.initialGeomAtom[i].en
    #                pots.finalGeomAtom[i].x             = pots.initialGeomAtom[i].x
    #                pots.finalGeomAtom[i].y             = pots.initialGeomAtom[i].y
    #                pots.finalGeomAtom[i].z             = pots.initialGeomAtom[i].z
    #            rebuild_g0x_com(pots)
    #            pots.finalGeomAtom = []

    return pots


### --- Parses the geometry and energy created by G0X and gives back the trajectory in a list of lists --- ###


def parse_output_jaguar(pots):
    if pots.debug >= 1:
        print_and_write(["parse_output_jaguar"], pots.verbose, pots.joblog)
    alist = []
    f_list = []
    force = 0

    # numatoms = 0

    collectgeom = 0
    pots.allEnergy = []
    pots.trajectoryGeomAtom = []
    the_file = open(pots.baselog, "r")
    for (idx, line) in enumerate(the_file):

        # Grab charge and multiplicity

        if "net molecular charge:" in line:
            pots.charge = line.rstrip().split()[-1]
        elif "multiplicity" in line:
            pots.multi = line.rstrip().split()[-1]
        elif "SCFE:" in line:
            pots.allEnergy.append(line.strip().split("hartrees")[0].split()[-1])
        elif "new geometry:" in line or "final geometry:" in line:
            collectgeom = 1
            pots.numatoms = 0
        elif len(line.strip()) == 0:
            collectgeom = 0
        elif collectgeom == 1 and "angstroms" not in line and "atom" not in line:
            alist.append(line)
            pots.numatoms += 1
        elif line.split()[0] == "displacement" and line.split()[1] == "maximum:":

            # Get the max displacement

            pots.finalDisp = float(line.split()[2])

        # Get forces

        if "-------------" in line and force >= 6:
            force = 0
        if force >= 1:
            force += 1
        if "forces (hartrees/bohr)" in line:
            force = 1
        if force >= 5:
            f_list.append(line)

    the_file.close()

    # Get the final energy from the energy list

    pots.finalenergy = float(pots.allEnergy[-1])
    alistdedup = [
        x
        for x in alist
        if "------------------------------------------------------" not in x
    ]

    # #####################################################################
    # Create the trajectory

    pots.trajectoryGeomAtom = [0] * (
        len(alistdedup) / pots.numatoms * (pots.numatoms + 2)
    )
    skip = 0
    for i in range(len(alistdedup)):
        if i % pots.numatoms == 0 and i / pots.numatoms != len(pots.allEnergy):
            pots.trajectoryGeomAtom[i + skip] = int(pots.numatoms)
            pots.trajectoryGeomAtom[i + skip + 1] = "Energy: " + str(
                pots.allEnergy[i / pots.numatoms]
            )
            skip += 2
        elif i % pots.numatoms == 0 and i / pots.numatoms >= len(pots.allEnergy):
            pots.trajectoryGeomAtom[i + skip] = int(pots.numatoms)
            pots.trajectoryGeomAtom[i + skip + 1] = "Energy: " + str(pots.allEnergy[-1])
            skip += 2
        line = alistdedup[i].split()
        pots.trajectoryGeomAtom[i + skip] = etaatom.Atom()
        pots.trajectoryGeomAtom[i + skip].e = "".join(
            [j for j in line[0] if not j.isdigit()]
        )
        pots.trajectoryGeomAtom[i + skip].en = etaatom.get_element_num(
            pots.trajectoryGeomAtom[i + skip].e
        )
        pots.trajectoryGeomAtom[i + skip].x = float(line[1])
        pots.trajectoryGeomAtom[i + skip].y = float(line[2])
        pots.trajectoryGeomAtom[i + skip].z = float(line[3])

    # #####################################################################
    # Extract the final geometry from the trajectory

    pots.finalGeomAtom = [0] * (pots.trajectoryGeomAtom[0] + 2)
    for i in range(len(pots.finalGeomAtom)):
        j = len(pots.trajectoryGeomAtom) - len(pots.finalGeomAtom) + i
        if i <= 1:
            pots.finalGeomAtom[i] = pots.trajectoryGeomAtom[j]
        elif i >= 2:
            pots.finalGeomAtom[i] = etaatom.Atom()
            pots.finalGeomAtom[i].e = pots.trajectoryGeomAtom[j].e
            pots.finalGeomAtom[i].en = pots.trajectoryGeomAtom[j].en
            pots.finalGeomAtom[i].x = pots.trajectoryGeomAtom[j].x
            pots.finalGeomAtom[i].y = pots.trajectoryGeomAtom[j].y
            pots.finalGeomAtom[i].z = pots.trajectoryGeomAtom[j].z

    # #####################################################################
    # Extract the final forces

    pots.allForce = [0] * pots.numatoms
    for i in range(len(pots.allForce)):
        j = len(f_list) - len(pots.allForce) + i
        line = f_list[j].strip().split()
        pots.allForce[i] = line[2] + "  " + line[3] + "  " + line[4]
    return pots


### --- G0X Optimization Logic --- ###


def opt_logic_jaguar(pots):
    if pots.debug >= 1:
        print_and_write(["OPTLogicJaguar"], pots.verbose, pots.joblog)

    # Check to see if the job finished optimization.

    count = 0
    f = open(pots.baselog, "r")
    for line in reversed(open(pots.baselog).readlines()):
        if "**           Geometry optimization complete           **" in line:
            if pots.jobtype == "SCAN" or pots.jobtype == "TS-SCAN":
                pots.optimized = True
                print_and_write(
                    [
                        "This job will continue due to: 'Normal termination' and/or 'Optimization Complete'."
                    ],
                    pots.verbose,
                    pots.joblog,
                )
            else:
                pots.optimized = True
                pots.kill = True
                print_and_write(
                    [
                        "This job will end due to: 'Normal termination' and/or 'Optimization Complete'."
                    ],
                    pots.verbose,
                    pots.joblog,
                )
            break
        count += 1
        if count > 100000:
            break

    # ## Check for jitter in the energy

    print_and_write(["Checking for jitter in the energies."], pots.verbose, pots.joblog)
    jitter = 0
    f = open(pots.scratchdir + "/optenergies-" + pots.jobname + ".txt", "a")
    for (i, line) in enumerate(pots.allEnergy):
        f.write(str(line) + "\n")
        if i > len(pots.allEnergy) / 2:
            if (
                pots.allEnergy[i] < pots.allEnergy[i - 1]
                and pots.allEnergy[i - 1] > pots.allEnergy[i - 2]
                or pots.allEnergy[i] > pots.allEnergy[i - 1]
                and pots.allEnergy[i - 1] < pots.allEnergy[i - 2]
            ):
                jitter += 1
    if jitter > len(pots.allEnergy) / 4 and not pots.optimized:
        print_and_write(
            ["Jitters in last half of OPT: " + str(jitter)], pots.verbose, pots.joblog
        )

        # + " > " + str(len(allEnergy)/4) + "\n")

        print_and_write(
            [
                "Due to the jitter in energies from the last OPT iteration,",
                "the step size of the optimization algorithm will be shortened.",
            ],
            pots.verbose,
            pots.joblog,
        )
        pots.shortenStepSize = 1

    # Check to see if the job ran too long.

    f = open(pots.scratchdir + "/optenergies-" + pots.jobname + ".txt", "r")
    fileenergylist = f.readlines()
    if len(fileenergylist) > pots.stepNumber:
        pots.kill = True
        print_and_write(
            ["This job will end due to: exhausted all optimization steps."],
            pots.verbose,
            pots.joblog,
        )
    f.close()

    # Check if a 'break' has occured

    breaklogic = 0
    line = pots.QMcomment.split()
    for var in line:
        if pots.breaktype == "bond" and breaklogic == 4:
            pots.breakpoint = float(var)
            breaklogic = 0
        elif pots.breaktype == "angle" and breaklogic == 5:
            pots.breakpoint = float(var)
            breaklogic = 0
        elif pots.breaktype == "dihedral" and breaklogic == 6:
            pots.breakpoint = float(var)
            breaklogic = 0
        elif breaklogic >= 2:
            pots.breakatoms.append(int(var))
            breaklogic += 1
        elif breaklogic == 1:
            pots.breaktype = var
            breaklogic += 1
        if var == "Break:":
            breaklogic = 1
    if pots.breaktype != "":
        print_and_write(
            [
                "Break point logic implemented.",
                "Checking if "
                + str(pots.breaktype)
                + " "
                + str(pots.breakatoms)
                + " are within "
                + str(pots.breakpoint)
                + ".",
            ],
            pots.verbose,
            pots.joblog,
        )

    # Break point logic implementation

    if pots.breaktype == "bond":
        if pots.breakpoint < etanumpy.get_distance(
            pots.breakatoms[0] - 1, pots.breakatoms[1] - 1, pots.finalGeomAtom
        ):
            pots.kill = True
            print_and_write(
                ["This job will end due to: a break point failure."],
                pots.verbose,
                pots.joblog,
            )
    elif pots.breaktype == "angle":
        if pots.breakpoint < etanumpy.get_angle(
            pots.breakatoms[0] - 1,
            pots.breakatoms[1] - 1,
            pots.breakatoms[2] - 1,
            pots.finalGeomAtom,
        ):
            pots.kill = True
            print_and_write(
                ["This job will end due to: a break point failure."],
                pots.verbose,
                pots.joblog,
            )
    elif pots.breaktype == "dihedral":
        if pots.breakpoint < etanumpy.get_dihedral(
            pots.breakatoms[0] - 1,
            pots.breakatoms[1] - 1,
            pots.breakatoms[2] - 1,
            pots.breakatoms[3] - 1,
            pots.finalGeomAtom,
        ):
            pots.kill = True
            print_and_write(
                ["This job will end due to: a break point failure."],
                pots.verbose,
                pots.joblog,
            )

    # ### TODO grab the rerun file for the speedup
    # ## Implement speedups for Jaguar
    # Read in the restart file and get all the goods

    if not pots.kill and not pots.optimized:
        pots.QMconfig = []
        pots.QMendlines = []
        pots.QMcomment = ""
        pots.numatoms = 0
        initialgeom = []
        xyzline = 0
        the_file = open(pots.pwdirectory + "/restart.in", "r")
        for line in the_file:
            if xyzline == 2:
                pots.QMendlines.append(line)
            elif xyzline == 1 and "&" in line:
                xyzline = 2
            elif xyzline == 0 and "&zmat" in line:
                xyzline = 1
            elif "MAEFILE:" in line:
                pots.QMcomment = line.strip()
            elif xyzline == 0:
                pots.QMconfig.append(line)
        the_file.close()

        # Stepsize logic implementation

        iop = 0
        newsub = 0
        if pots.shortenStepSize == 1:
            for (i, line) in enumerate(pots.QMconfig):
                if "trust=" in line:
                    var = line.split()
                    iop = 1
                    for sub in var:
                        if "trust=" in sub:
                            print_and_write(
                                ["Your step size was: " + str(sub)],
                                pots.verbose,
                                pots.joblog,
                            )
                            newsub = "trust=" + str(
                                int(
                                    round(float(sub.split("=")[1].split(")")[0]) * 0.75)
                                )
                            )
                            print_and_write(
                                ["It will now be: " + str(newsub)],
                                pots.verbose,
                                pots.joblog,
                            )
                            pots.QMconfig[i] = line.replace(sub, newsub)
            if iop == 0:
                for (i, line) in enumerate(pots.QMconfig):
                    if "&gen" in line:
                        print_and_write(
                            [
                                "We are now implementing a smaller step size of trust=0.2"
                            ],
                            pots.verbose,
                            pots.joblog,
                        )
                        pots.QMconfig[i] = line + "trust=0.2\n"

    # for i, line in enumerate(pots.QMconfig):
    #  if "&gen" in line and "readfc" not in line:
    #    print_and_write(["We are now implementing readfc for the next iteration."], pots.verbose, pots.joblog)
    #    var = line.split()
    #    for sub in var:
    #      if "opt=(" in sub:
    #        newsub = sub[:-1] + ",readfc)"
    #        pots.QMconfig[i] = line.replace(sub, newsub)

    return pots


### --- Output the G09 com file for subsequent computation --- ###


def rebuild_jaguar_in(pots):
    if pots.debug >= 1:
        print_and_write(["rebuild_jaguar_in"], pots.verbose, pots.joblog)

    # Output the new geometry in Jaguar input format

    datetime = time.strftime("%Y.%m.%d-%H:%M:%S")

    # Backup .in file

    if os.path.exists(pots.baseinput):
        if os.path.isfile(pots.baseinput):
            shutil.copyfile(
                pots.baseinput,
                pots.scratchdir + "/Previous_Files/" + datetime + "-" + pots.inputfile,
            )

    # Backup .log file

    if os.path.exists(pots.jobname + ".log"):
        if os.path.isfile(pots.jobname + ".log"):
            shutil.copyfile(
                pots.jobname + ".log",
                pots.scratchdir
                + "/Previous_Files/"
                + datetime
                + "-"
                + pots.jobname
                + ".log",
            )

    # Backup .out file

    if os.path.exists(pots.baselog):
        if os.path.isfile(pots.baselog):
            shutil.copyfile(
                pots.baselog,
                pots.scratchdir
                + "/Previous_Files/"
                + datetime
                + "-"
                + pots.jobname
                + ".out",
            )

    # print ofile

    f = open(pots.baseinput, "w+")
    for (i, var) in enumerate(pots.QMconfig):
        f.write(var)
    f.write("&zmat\n")
    for i in range(2, len(pots.finalGeomAtom)):
        line = (
            str(pots.finalGeomAtom[i].e)
            + str(i - 1)
            + " "
            + str(pots.finalGeomAtom[i].x)
            + " "
            + str(pots.finalGeomAtom[i].y)
            + " "
            + str(pots.finalGeomAtom[i].z)
        )
        f.write(line + "\n")
    f.write("&\n")
    for (i, var) in enumerate(pots.QMendlines):
        f.write(var)
    f.close()
    return


######################################################################
### END OF LIBRARY
######################################################################
