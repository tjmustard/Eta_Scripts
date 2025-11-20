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
import os
import sys
from sys import *

from EtaLib import etaatom

### --- Arguments --- ###

program = "nVisualize-Pymol-PDB.py"
dist = False
EtaDir = os.environ["ETADIR"]
build = "Pymolv1.3-Build.pml"
visual = "Pymolv1.3-Visualize.pml"

### Read command line args

try:
    (myopts, args) = getopt.getopt(sys.argv[1:], "b:v:dh")
except getopt.GetoptError:
    print(program + " -d (generate distances) -h")
    sys.exit(2)

###############################
# o == option
# a == argument passed to the o
###############################

for (o, a) in myopts:
    if o == "-d":
        dist = True
    elif o == "-v":
        visual = a
    elif o == "-b":
        build = a
    elif o == "-h":
        print(program + " -h")
        sys.exit(0)
    else:
        print("Usage: %s -h" % sys.argv[0])
        sys.exit(0)

### --- Make the Pymol folder --- ###

etaatom.make_dir("Pymol-Picture")

### --- Grab the build and visualize configuration file --- ###

f = open(EtaDir + "/Pymol/snippets/" + build, "r")
buildlines = f.readlines()
f.close()
f = open(EtaDir + "/Pymol/snippets/" + visual, "r")
vislines = f.readlines()
f.close()
os.system("cp " + EtaDir + "/Pymol/snippets/Pymolv1.3-Visualize.pml Pymol-Picture/")

### --- Creat the Visual-Script.pml file --- ###

f = open("Pymol-Picture/Visual-Script.pml", "w")
f.write(
    """#==========================================
#Start
"""
)

### --- Iterating through a folder of files --- ###

for i in os.listdir(os.getcwd()):
    if i.endswith(".pdb"):
        ifile = i
        basename = etaatom.basename(ifile, ".pdb")
        os.system("cp " + basename + ".pdb Pymol-Picture/" + basename + "_orig.pdb")
        os.system("cp " + basename + ".pdb Pymol-Picture/" + basename + "_high.pdb")
        os.system("cp " + basename + ".pdb Pymol-Picture/" + basename + "_med.pdb")
        os.system("cp " + basename + ".pdb Pymol-Picture/" + basename + "_low.pdb")
        os.system("cp " + basename + ".pdb Pymol-Picture/" + basename + "_dat.pdb")
        f.write("load " + basename + "_orig.pdb\n")
        f.write("load " + basename + "_high.pdb\n")
        f.write("load " + basename + "_med.pdb\n")
        f.write("load " + basename + "_low.pdb\n")
        f.write("load " + basename + "_dat.pdb\n")
        f.write("#grouping objects\n")
        f.write(
            "#------------------------------------------------------------------------------------------------------------------------------------------\n"
        )
        f.write("group " + basename + ", " + basename + "_orig\n")
        f.write("group " + basename + ", " + basename + "_high\n")
        f.write("group " + basename + ", " + basename + "_med\n")
        f.write("group " + basename + ", " + basename + "_low\n")
        f.write("group " + basename + ", " + basename + "_dat\n")
        f.write("\n")

        if dist:
            f.write(
                "#------------------------------------------------------------------------------------------------------------------------------------------\n"
            )
            f.write("#DISTANCE, DISTANCES\n")
            f.write("#Defining distances (dist)\n")
            f.write(
                "distance Distance_HX_"
                + basename
                + "_DIST = (elem H) and "
                + basename
                + "_high,  (neighbor elem H) and "
                + basename
                + "_high, 2.0;color Grey, Distance_HX\n"
            )
            f.write(
                "distance Distance_CX_"
                + basename
                + "_DIST = (elem C) and "
                + basename
                + "_high, (neighbor elem C) and "
                + basename
                + "_high, 2.0;color Grey, Distance_CX\n"
            )
            f.write(
                "distance Distance_HB_"
                + basename
                + "_DIST = (elem B) and "
                + basename
                + "_high, (elem H and neighbor elem B) and "
                + basename
                + "_high, 2.0;color purple, Distance_HB\n"
            )
            f.write(
                "distance Distance_HC_"
                + basename
                + "_DIST = (elem C) and "
                + basename
                + "_high, (elem H and neighbor elem C) and "
                + basename
                + "_high, 2.0;color Black, Distance_HC\n"
            )
            f.write(
                "distance Distance_HN_"
                + basename
                + "_DIST = (elem N) and "
                + basename
                + "_high, (elem H and neighbor elem N) and "
                + basename
                + "_high, 2.0;color Blue, Distance_HN\n"
            )
            f.write(
                "distance Distance_HO_"
                + basename
                + "_DIST = (elem O) and "
                + basename
                + "_high, (elem H and neighbor elem O) and "
                + basename
                + "_high, 2.0;color Red, Distance_HO\n"
            )
            f.write(
                "distance Distance_HSi_"
                + basename
                + "_DIST = (name Si) and "
                + basename
                + "_high, (elem H and neighbor elem Si) and "
                + basename
                + "_high, 2.0;color Orange, Distance_HSi\n"
            )
            f.write(
                "distance Distance_HP_"
                + basename
                + "_DIST = (elem P) and "
                + basename
                + "_high, (elem H and neighbor elem P) and "
                + basename
                + "_high, 2.0;color Orange, Distance_HP\n"
            )
            f.write(
                "distance Distance_HS_"
                + basename
                + "_DIST = (elem S) and "
                + basename
                + "_high, (elem H and neighbor elem S) and "
                + basename
                + "_high, 2.0;color Orange, Distance_HS\n"
            )
            f.write(
                "distance Distance_BC_"
                + basename
                + "_DIST = (elem B) and "
                + basename
                + "_high, (elem C and neighbor elem B) and "
                + basename
                + "_high, 2.0;color Black, Distance_BC\n"
            )
            f.write(
                "distance Distance_BN_"
                + basename
                + "_DIST = (elem N) and "
                + basename
                + "_high, (elem B and neighbor elem N) and "
                + basename
                + "_high, 2.0;color Blue, Distance_BN\n"
            )
            f.write(
                "distance Distance_BO_"
                + basename
                + "_DIST = (elem O) and "
                + basename
                + "_high, (elem B and neighbor elem O) and "
                + basename
                + "_high, 2.0;color Red, Distance_BO\n"
            )
            f.write(
                "distance Distance_CC_"
                + basename
                + "_DIST = (elem C) and "
                + basename
                + "_high, (elem C) and "
                + basename
                + "_high, 1.6;color Black, Distance_CC\n"
            )
            f.write(
                "distance Distance_CN_"
                + basename
                + "_DIST = (elem N) and "
                + basename
                + "_high, (elem C and neighbor elem N) and "
                + basename
                + "_high, 2.0;color Blue, Distance_CN\n"
            )
            f.write(
                "distance Distance_CO_"
                + basename
                + "_DIST = (elem O) and "
                + basename
                + "_high, (elem C and neighbor elem O) and "
                + basename
                + "_high, 2.0;color Red, Distance_CO\n"
            )
            f.write(
                "distance Distance_CSi_"
                + basename
                + "_DIST = (name Si) and "
                + basename
                + "_high, (elem C and neighbor elem Si) and "
                + basename
                + "_high, 2.0;color Orange, Distance_CSi\n"
            )
            f.write(
                "distance Distance_CP_"
                + basename
                + "_DIST = (elem P) and "
                + basename
                + "_high, (elem C and neighbor elem P) and "
                + basename
                + "_high, 2.0;color Orange, Distance_CP\n"
            )
            f.write(
                "distance Distance_CS_"
                + basename
                + "_DIST = (elem S) and "
                + basename
                + "_high, (elem C and neighbor elem S) and "
                + basename
                + "_high, 2.0;color Orange, Distance_CS\n"
            )
            f.write(
                "distance Distance_NN_"
                + basename
                + "_DIST = (elem N) and "
                + basename
                + "_high, (elem N and neighbor elem N) and "
                + basename
                + "_high, 2.0;color Blue, Distance_NN\n"
            )
            f.write(
                "distance Distance_NO_"
                + basename
                + "_DIST = (elem O) and "
                + basename
                + "_high, (elem N neighbor elem O) and "
                + basename
                + "_high, 2.0;color Blue, Distance_NO\n"
            )
            f.write(
                "distance Distance_NSi_"
                + basename
                + "_DIST = (name Si) and "
                + basename
                + "_high, (elem N and neighbor elem Si) and "
                + basename
                + "_high, 2.0;color Orange, Distance_NSi\n"
            )
            f.write(
                "distance Distance_NP_"
                + basename
                + "_DIST = (elem P) and "
                + basename
                + "_high, (elem N and neighbor elem P) and "
                + basename
                + "_high, 2.0;color Orange, Distance_NP\n"
            )
            f.write(
                "distance Distance_NS_"
                + basename
                + "_DIST = (elem S) and "
                + basename
                + "_high, (elem N and neighbor elem S) and "
                + basename
                + "_high, 2.0;color Orange, Distance_NS\n"
            )
            f.write(
                "distance Distance_OO_"
                + basename
                + "_DIST = (elem O) and "
                + basename
                + "_high, (elem O and neighbor elem O) and "
                + basename
                + "_high, 2.0;color Red, Distance_OO\n"
            )
            f.write(
                "distance Distance_OSi_"
                + basename
                + "_DIST = (name Si) and "
                + basename
                + "_high, (elem O and neighbor elem Si) and "
                + basename
                + "_high, 2.0;color Orange, Distance_OSi\n"
            )
            f.write(
                "distance Distance_OP_"
                + basename
                + "_DIST = (elem P) and "
                + basename
                + "_high, (elem O and neighbor elem P) and "
                + basename
                + "_high, 2.0;color Orange, Distance_OP\n"
            )
            f.write(
                "distance Distance_OS_"
                + basename
                + "_DIST = (elem S) and "
                + basename
                + "_high, (elem O and neighbor elem S) and "
                + basename
                + "_high, 2.0;color Orange, Distance_OS\n"
            )
            f.write(
                "distance Distance_PRh_"
                + basename
                + "_DIST = (elem P) and "
                + basename
                + "_high, (name Rh and neighbor elem P) and "
                + basename
                + "_high, 2.5;color Orange, Distance_PRh\n"
            )
            f.write("\n")
            f.write(
                "#------------------------------------------------------------------------------------------------------------------------------------------\n"
            )
            f.write("#HYDROGEN-BONDING, H-BOND, H-BONDING\n")
            f.write("#Defining hydrogen bonding distances (h_bond, hbond)\n")
            f.write("set h_bond_cutoff_center, 3.6\n")
            f.write("set h_bond_cutoff_edge, 3.2\n")
            f.write("\n")
            f.write(
                "#------------------------------------------------------------------------------------------------------------------------------------------\n"
            )
            f.write(
                "#TRANSITION STATE BONDS (TS, TSS, TS-BOND, TS-BONDS, TSBOND, TSBONDS)\n"
            )
            f.write("#Defining ts bonds (ts,ts_bond)\n")
            f.write(
                "distance TS_C-C_"
                + basename
                + "_DIST = (/////C) and "
                + basename
                + "_high, (/////C) and "
                + basename
                + "_high, 3.0\n"
            )
            f.write("\n")
            f.write(
                "#------------------------------------------------------------------------------------------------------------------------------------------\n"
            )
            f.write("#STERIC (STERIC, STERICS)\n")
            f.write(
                "#Defining steric interactions between hydrogens (hh, steric, sterics)\n"
            )
            f.write(
                "distance Steric_HH_21_"
                + basename
                + "_DIST = (elem H) and "
                + basename
                + "_high, (elem H) and "
                + basename
                + "_high, 2.1; color Grey, Steric_HH_21\n"
            )
            f.write(
                "distance Steric_HH_22_"
                + basename
                + "_DIST = (elem H) and "
                + basename
                + "_high, (elem H) and "
                + basename
                + "_high, 2.2; color Grey, Steric_HH_22\n"
            )
            f.write(
                "distance Steric_HH_23_"
                + basename
                + "_DIST = (elem H) and "
                + basename
                + "_high, (elem H) and "
                + basename
                + "_high, 2.3; color Grey, Steric_HH_23\n"
            )
            f.write(
                "distance Steric_HH_24_"
                + basename
                + "_DIST = (elem H) and "
                + basename
                + "_high, (elem H) and "
                + basename
                + "_high, 2.4; color Grey, Steric_HH_24\n"
            )
            f.write(
                "#------------------------------------------------------------------------------------------------------------------------------------------\n"
            )
            f.write("#ELECTROSTATIC (ELECTROSTATIC, ELECTROSTATICS, ESP, ESPS)\n")
            f.write(
                "#Defining electrostatic interactions and contacts (electrostatic, esp, esc)\n"
            )
            f.write(
                "select XH_32, /////H and (neighbor /////N or neighbor /////O or neighbor /////S)\n"
            )
            f.write(
                "distance Classic_ESP_NH_32_"
                + basename
                + "_DIST = /////N, XH_32, 3.2; color Navy_Blue, Classic_ESP_NH_*\n"
            )
            f.write(
                "distance Classic_ESP_OH_32_"
                + basename
                + "_DIST = /////O, XH_32, 3.2; color Red, Classic_ESP_OH_*\n"
            )
            f.write(
                "distance Classic_ESP_HX_32_"
                + basename
                + "_DIST = XH_32 and "
                + basename
                + "_high, (/////N or /////O or /////S) and "
                + basename
                + "_high, 3.2; color Grey, Classic_ESP_HX_*; delete XH_32\n"
            )
            f.write(
                "distance ESP_CH_32_"
                + basename
                + "_DIST = (elem C) and "
                + basename
                + "_high, ((elem H) and not (neighbor (elem H))) and "
                + basename
                + "_high, 3.2; color Grey, ESP_CH_*\n"
            )
            f.write(
                "distance ESP_NH_32_"
                + basename
                + "_DIST = (elem N) and "
                + basename
                + "_high, ((elem H) and not (neighbor (elem H))) and "
                + basename
                + "_high, 3.2; color Navy_Blue, ESP_NH_*\n"
            )
            f.write(
                "distance ESP_OH_32_"
                + basename
                + "_DIST = (elem O) and "
                + basename
                + "_high, ((elem H) and not (neighbor (elem H))) and "
                + basename
                + "_high, 3.2; color Pink, ESP_OH_*\n"
            )
            f.write(
                "distance ESP_SH_32_"
                + basename
                + "_DIST = (elem S) and "
                + basename
                + "_high, ((elem H) and not (neighbor (elem H))) and "
                + basename
                + "_high, 3.2; color Yellow, ESP_SH_*\n"
            )
            f.write(
                "distance ESP_BN_32_"
                + basename
                + "_DIST = (elem B) and "
                + basename
                + "_high, (elem N) and "
                + basename
                + "_high, 3.2; color Pink, ESP_BN_*\n"
            )
            f.write(
                "distance ESP_OO_32_"
                + basename
                + "_DIST = (elem O) and "
                + basename
                + "_high, (elem O) and "
                + basename
                + "_high, 3.2; color Red, ESP_OO_*\n"
            )
            f.write(
                "distance ESP_ON_32_"
                + basename
                + "_DIST = (elem O) and "
                + basename
                + "_high, (elem N) and "
                + basename
                + "_high, 3.2; color purple, ESP_ON_*\n"
            )
            f.write(
                "distance ESP_SS_32_"
                + basename
                + "_DIST = (elem S) and "
                + basename
                + "_high, (elem S) and "
                + basename
                + "_high, 3.2; color sulfur, ESP_SS_*\n"
            )
            f.write(
                "distance ESP_SO_32_"
                + basename
                + "_DIST = (elem S) and "
                + basename
                + "_high, (elem O) and "
                + basename
                + "_high, 3.2; color violet, ESP_SO_*\n"
            )
            f.write(
                "distance ESP_SN_32_"
                + basename
                + "_DIST = (elem S) and "
                + basename
                + "_high, (elem N) and "
                + basename
                + "_high, 3.2; color purple, ESP_SN_*\n"
            )
            f.write(
                "distance ESP_ClH_32_"
                + basename
                + "_DIST = (name Cl) and "
                + basename
                + "_high, (elem H) and "
                + basename
                + "_high, 3.2; color Green, ESP_ClH_*\n"
            )
            f.write("\n")
            f.write(
                "#------------------------------------------------------------------------------------------------------------------------------------------\n"
            )
            f.write(
                "distance TMBond_SiO_"
                + basename
                + "_DIST = (elem O) and "
                + basename
                + "_high, (elem S) and "
                + basename
                + "_high, 2.5;  Bond_SiO\n"
            )
            f.write(
                "distance TMBond_SnO_"
                + basename
                + "_DIST = (elem O) and "
                + basename
                + "_high, (elem S) and "
                + basename
                + "_high, 2.5;  Bond_SnO\n"
            )
            f.write(
                "distance TMBond_ZnO_"
                + basename
                + "_DIST = (elem O) and "
                + basename
                + "_high, (elem Zn) and "
                + basename
                + "_high, 2.5;  Bond_ZnO\n"
            )
            f.write(
                "distance TMBond_InO_"
                + basename
                + "_DIST = (elem O) and "
                + basename
                + "_high, (elem In) and "
                + basename
                + "_high, 2.5;  Bond_InO\n"
            )
            f.write(
                "distance TMBond_GaO_"
                + basename
                + "_DIST = (elem O) and "
                + basename
                + "_high, (elem Ga) and "
                + basename
                + "_high, 2.5;  Bond_GaO\n"
            )
            f.write(
                "distance TMBond_AlO_"
                + basename
                + "_DIST = (elem O) and "
                + basename
                + "_high, (elem Al) and "
                + basename
                + "_high, 2.5;  Bond_AlO\n"
            )
            f.write("\n")
            f.write("group " + basename + ", *_" + basename + "_DIST\n")
    if dist:
        f.write(
            "#------------------------------------------------------------------------------------------------------------------------------------------\n"
        )
        f.write("#Defining all distances\n")
        f.write(
            "distance ALL_MainGroup_Distance = *, *, 1.75, 1; color Black, ALL_MainGroup_Distance\n"
        )
        f.write(
            "distance ALL_TMGroup_Distance = *, *, 2.6, 1; color Black, ALL_TMGroup_Distance\n"
        )
        f.write(
            "distance ALL_Distance= *, neighbor, 3.0, 1; color Black, ALL_Distance\n"
        )

for line in buildlines:
    f.write(line)
for line in vislines:
    f.write(line)

f.close()

os.system("open Pymol-Picture/Visual-Script.pml")

######################################################################
### END OF SCRIPT
######################################################################
