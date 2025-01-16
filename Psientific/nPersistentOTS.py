#!/usr/bin/env python2.7
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
import math
from sys import *
import sys
import getopt
from decimal import *
import shutil
import etaatom
import etanumpy
import etapots

### --- Arguments --- ###

program = 'nPersistentOTS.py'
inputfile = ''
outfile = ''
jobname = ''
pwdirectory = ''
jobtype = ''
QMtype = ''
path2QM = ''
cpus = 1
debug = 0
test = 0
maketemplate = False

### --- OK Now setup the variables for the first time --- ###
# etapots.reset_pots_variables()

pots = etapots.PersistentOTS()

### Read command line args

try:
    (myopts, args) = getopt.getopt(sys.argv[1:], 'i:o:n:d:t:j:p:c:hDTv'
                                   , ['template'])
except getopt.GetoptError:
    print program \
        + ' -i <input file> -o <output file> -n <job name> -d <directory> -t <QM code> -j <job type> -p <path to QM code> -h (print help)'
    sys.exit(2)

###############################
# o == option
# a == argument passed to the o
###############################

for (o, a) in myopts:
    if o == '-i':
        pots.inputfile = a
    elif o == '-o':
        pots.outfile = a
    elif o == '-n':
        pots.jobname = a
    elif o == '-d':
        pots.pwdirectory = a
    elif o == '-t':
        pots.QMtype = a
    elif o == '-j':
        pots.jobtype = a
    elif o == '-p':
        pots.path2QM = a
    elif o == '-c':
        cpus = a
    elif o == '-D':
        pots.debug += 1
    elif o == '-T':
        pots.test = 1
    elif o == '-v':
        pots.verbose = 1
    elif o == '--template':
        maketemplate = True
    elif o == '-h':
        print program + ' -idjtjphD'
        print '-----------------------------------------------------------------------------------------------------------------'
        print "	-i	(H2O.com)		| Input file."
        print "					|    For Turbomole jobs you will specify the define input file."
        print "					|    For Gaussian jobs you will specify the .com/.inp file."
        print "					|    For Terachem jobs you will specify the configuration file."
        print "					| "
        print "	-o	(H2O.log)		| Output file."
        print "					|    Usually the .log or .out file."
        print "					| "
        print "	-n	(H2O)			| Job name. "
        print "					|    Needs to be unique from other jobs in the same directory."
        print "					|    All outputs will have this string in their names."
        print "					|"
        print "	-d	(/home/user/H2O_OPT)	| Directory where files will be saved/updated."
        print "					|    If you are running the job on a cluster, this can be the originating "
        print "					|    folder of the input and output files. This way updates can be "
        print "					|    propogated from the node's local drive every iteration of the job."
        print "					| "
        print "	-t	(G09)			| QM code used (Gaussian (G03/G09), Turbomole (TM), Terachem (TC))"
        print "					|    Indicating the type of QM code used allows for correct parsing of "
        print "					|    output and input files."
        print "					|"
        print "	-j	(OPT)			| Job type (OPT, SCAN, TS-SCAN)"
        print "					|    OPT: "
        print "					|    Run subsequent iterations of the input file while updating "
        print "					|    the geometry. Useful for very long or difficult optimizations."
        print "					|    SCAN: "
        print "					|    In conjunction with a 'Scan' configuration file, this will"
        print "					|    run subsequent positionally restrained optimizations (POPT)"
        print "					|    until the number of steps has been fulfilled."
        print "					|    TS-SCAN: "
        print "					|    In conjunction with a 'Scan' configuration file, this will"
        print "					|    run subsequent positionally restrained optimizations (POPT)"
        print "					|    with added logic to take smaller scanning step sizes and "
        print "					|    the ability to reverse the direction of the scan to find the "
        print "					|    optimal input geometry for a TS search."
        print "					|    "
        print "	-p	(/share/apps/Gaussian)	| Path to the QM code"
        print "					|    "
        print "	-h				| Print this help screen"
        print "					| "
        print '-----------------------------------------------------------------------------------------------------------------'
        print 'Example: ' + program \
            + ' -i H2O.com -o H2O.log -n H20 -d /home/user/H2O_OPT -t G09 -j OPT -p /share/apps/Gaussian'
        sys.exit(0)
    else:
        print 'Usage: %s -i <input file> -o <output file> -n <job name> -d <directory> -t <QM code> -j <job type> -p <path to QM code> -h (print help)' \
            % sys.argv[0]
        sys.exit(0)

if maketemplate:
    etapots.make_temp_scan_file(pots)
    sys.exit(0)

######################################################################
### START OF SCRIPT
######################################################################
###Start off

pots.kill = False

###Set up variables

pots.baseinput = pots.pwdirectory + '/' + pots.inputfile
pots.baselog = pots.pwdirectory + '/' + pots.outfile
pots.joblog = pots.pwdirectory + '/' + pots.jobname + '.joblog'
pots.scratchdir = pots.pwdirectory + '/SCRATCH-' + pots.jobname
pots.trajectoryfile = pots.scratchdir + '/TRAJECTORY-' + pots.jobname \
    + '.xyz'
pots.scanfile = pots.pwdirectory + '/Scan-' + pots.jobname + '.txt'
pots.outScanFile = pots.pwdirectory + '/SCAN-' + pots.jobname + '.xyz'
pots.linearScanFile = pots.scratchdir + '/LINEAR-SCAN-' + pots.jobname \
    + '.xyz'

etapots.print_and_write([
    '################################################################################'
        ,
    '',
    '                 *********************************************',
    '                               nPersistentOTS.py',
    '                 *********************************************',
    '',
    ' Copyright (c) 2014, Thomas J. L. Mustard, O. Maduka Ogba, Paul Ha-Yeon Cheong'
        ,
    '',
    '                                   PHYC Group',
    '                            Oregon State University',
    '                               College of Science',
    '                            Department of Chemistry',
    '                                153 Gilbert Hall',
    '                              Corvallis OR, 97331',
    '                        E-mail:  mustardt@onid.orst.edu',
    '                              Ph.  (541)-737-2081',
    '                       http://phyc.chem.oregonstate.edu/',
    '',
    '################################################################################'
        ,
    '',
    '  All rights reserved.',
    '',
    '  Redistribution and use in source and binary forms, with or without'
        ,
    '  modification, are permitted provided that the following conditions are met:'
        ,
    '',
    '  * Redistributions of source code must retain the above copyright notice, this'
        ,
    '    list of conditions and the following disclaimer.',
    '',
    '  * Redistributions in binary form must reproduce the above copyright notice,'
        ,
    '    this list of conditions and the following disclaimer in the documentation'
        ,
    '    and/or other materials provided with the distribution.',
    '',
    '  * Neither the name of the {organization} nor the names of its',
    '    contributors may be used to endorse or promote products derived from'
        ,
    '    this software without specific prior written permission.',
    '',
    '  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"'
        ,
    '  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE'
        ,
    '  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE'
        ,
    '  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE'
        ,
    '  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL'
        ,
    '  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR'
        ,
    '  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER'
        ,
    '  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,'
        ,
    '  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE'
        ,
    '  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.'
        ,
    '',
    '################################################################################'
        ,
    ], pots.verbose, pots.joblog)

###Start the loop

while not pots.kill:
    etapots.print_and_write(['', 'QM used for job: '
                            + str(pots.QMtype), 'Type of job: '
                            + str(pots.jobtype), ''], pots.verbose,
                            pots.joblog)
    if debug >= 1:
        etapots.print_and_write([
            pots.baseinput,
            pots.baselog,
            pots.scratchdir,
            pots.trajectoryfile,
            pots.scanfile,
            pots.outScanFile,
            ], pots.verbose, pots.joblog)

  # #####################################################################
  # ## --- #Make required directories if needed --- ###

    if pots.debug >= 1:
        etapots.print_and_write(['makeDirectories'], pots.verbose,
                                pots.joblog)
    etaatom.make_dir(pots.scratchdir)
    etaatom.make_dir(pots.scratchdir + '/Previous_Files')

  # #####################################################################
  # ## --- Grab the last configuration from the input file --- ###

    if pots.QMtype == 'G09' or pots.QMtype == 'G03':
        pots = etapots.parse_input_g0x(pots)
    elif pots.QMtype == 'Turbomole':
        etapots.print_and_write(['Not sure what to due here.',
                                pots.verbose, pots.joblog])
    elif pots.QMtype == 'Terachem':
        etapots.print_and_write(['Not sure what to due here.',
                                pots.verbose, pots.joblog])
    elif pots.QMtype.upper() == 'JAGUAR':
        pots = etapots.parse_input_jaguar(pots)

  # Save the original configurations of the file

    if pots.QMconfigORIG == []:
        for i in pots.QMconfig:
            pots.QMconfigORIG.append(i)
    if pots.QMcommentORIG == '':
        pots.QMcommentORIG = pots.QMcomment
    if pots.QMendlinesORIG == []:
        for i in pots.QMendlines:
            pots.QMendlinesORIG.append(i)

    if pots.test != 1:

    # ## Run the QM code ###
    # ## G09

        if pots.QMtype == 'G09':
            etapots.print_and_write(['Running G09'], pots.verbose,
                                    pots.joblog)
            try:
                os.system('g09 < ' + pots.baseinput + ' > '
                          + pots.baselog)
            except getopt.GetoptError:
                etapots.print_and_write(['Gaussian09 Failed!!!!'],
                        pots.verbose, pots.joblog)
        elif pots.QMtype == 'G03':

    # ## G03

            etapots.print_and_write(['Running G03'], pots.verbose,
                                    pots.joblog)
            try:
                os.system('g03 < ' + pots.baseinput + ' > '
                          + pots.baselog)
            except getopt.GetoptError:
                etapots.print_and_write(['Gaussian03 Failed!!!!'],
                        pots.verbose, pots.joblog)
        elif pots.QMtype == 'Turbomole':

    # ## Turbomole

            etapots.print_and_write(['Running Turbomole Jobex'],
                                    pots.verbose, pots.joblog)
            try:
                os.system('jobex -ri -c 300 -v -outfile '
                          + pots.baselog)
            except getopt.GetoptError:
                etapots.print_and_write(['Turbomole Failed!!!!'],
                        pots.verbose, pots.joblog)
        elif pots.QMtype == 'Terachem':

    # ## Terachem

            etapots.print_and_write(['Running Terachem'], pots.verbose,
                                    pots.joblog)
            try:
                os.system('/share/psw/TeraChem/hidden/RUN-TC.sh '
                          + pots.baseinput + ' >> ' + pots.baselog)
            except getopt.GetoptError:
                etapots.print_and_write(['TeraChem Failed!!!!'],
                        pots.verbose, pots.joblog)
        elif pots.QMtype.upper() == 'JAGUAR':

    # ## Jaguar

            etapots.print_and_write(['Running Jaguar'], pots.verbose,
                                    pots.joblog)
            try:
                os.system('/export/apps/schrodinger2014-2/jaguar run -TPP '
                           + cpus + ' -WAIT -HOST localhost '
                          + pots.baseinput + ' >> ' + pots.baselog)
            except getopt.GetoptError:
                etapots.print_and_write(['Jaguar Failed!!!!'],
                        pots.verbose, pots.joblog)

  # #####################################################################
  # ## --- Grab the last geom from the output file --- ###

    if pots.QMtype == 'G09' or pots.QMtype == 'G03':

    # pots.trajectoryGeomAtom, pots.charge, pots.multi = etaatom.parse_output_g0x_traj(baselog)  #Parse log file for final geometry
    # pots.QMconfig, pots.QMendlines, pots.charge, pots.multi, pots.inputGeomAtom = etaatom.parse_input_g0x(baseinput)

        pots = etapots.parse_output_g0x(pots)  # Parse log file for final geometry
        pots = etapots.parse_input_g0x(pots)
    elif pots.QMtype == 'Turbomole':

    # pots.trajectoryGeomAtom, pots.charge, pots.multi = etaatom.parse_output_tm_traj(baselog)  #Parse log file for final geometry

        pots = etaatom.parse_output_tm_traj(pots)  # Parse log file for final geometry
    elif pots.QMtype == 'Terachem':

    # pots.trajectoryGeomAtom = etaatom.parse_optim_terachem_traj(pwdirectory + TCscratchDir + "optim.xyz")  #Parse optim.xyz file for final geometry

        pots = etapots.parse_optim_terachem_traj(pots)  # Parse optim.xyz file for final geometry
    elif pots.QMtype.upper() == 'JAGUAR':
        pots = etapots.parse_output_jaguar(pots)
        pots = etapots.parse_input_jaguar(pots)

  # #####################################################################
  # Write the trajectory out

    etaatom.output_xyz_traj(pots.trajectoryfile,
                            pots.trajectoryGeomAtom)

  # #####################################################################
  # ##For OPT jobs
  # if pots.jobtype == "OPT":
  # #####################################################################
  # ## --- Output a new input file for the next step in the scan --- ###

    if pots.QMtype == 'G09' or pots.QMtype == 'G03':

    # Decide what the next step should be

        pots = etapots.opt_logic_g0x(pots)
    elif pots.QMtype == 'Turbomole':

    # Move/convert files for the next computation
    # etapots.rebuild_g0x_com(pots)
    # Decide what the next step should be

        etapots.print_and_write(['OPT LOGIC HERE!'], pots.verbose,
                                pots.joblog)
        pots = etapots.opt_logic_tm(pots)
    elif pots.QMtype == 'Terachem':

    # output the fixed geometry
    # outputFixUnfixGeom(pots.finalGeomAtom)
    # Move/convert files for the next computation
    # etapots.rebuild_turbomole_tmol(pots)
    # Decide what the next step should be

        etapots.print_and_write(['OPT LOGIC HERE!'], pots.verbose,
                                pots.joblog)
        pots = etapots.opt_logic_tc(pots)
    elif pots.QMtype.upper() == 'JAGUAR':

    # Move/convert files for the next computation
    # etapots.rebuild_terachem_xyz(pots)
    # ##For SCAN and TS-SCAN jobs
    # Decide what the next step should be

        pots = etapots.opt_logic_jaguar(pots)

    # Move/convert files for the next computation
    # etapots.rebuild_jaguar_in(pots)
  # #####################################################################

    if (pots.jobtype == 'SCAN' or pots.jobtype == 'TS-SCAN') \
        and pots.optimized:
        if not os.path.isfile(pots.scanfile):
            etapots.print_and_write(['The scan file nessecary for both the SCAN and TS-SCAN function does not exist.'
                                    ,
                                    'Please include this file and submit again.'
                                    ], pots.verbose, pots.joblog)
            sys.exit(0)

    # Parse the scan file and print out important information

        pots = etapots.parse_scan_file(pots)
        pots = etapots.generate_min_ts_scan_size(pots)
        etapots.print_and_write(['Atoms to be moved: '
                                + str(pots.scanType) + '  '
                                + str(pots.atomList)], pots.verbose,
                                pots.joblog)

    # Update the second line for the SCAN file

        if pots.scanType == 'bond':
            pots.finalGeomAtom[1] = pots.finalGeomAtom[1] + '   Bond: ' \
                + str(etanumpy.get_distance(pots.atomList[0] + 1,
                      pots.atomList[1] + 1, pots.finalGeomAtom))
        if pots.scanType == 'angle':
            pots.finalGeomAtom[1] = pots.finalGeomAtom[1] \
                + '   Angle: ' \
                + str(etanumpy.get_angle(pots.atomList[0] + 1,
                      pots.atomList[1] + 1, pots.atomList[2] + 1,
                      pots.finalGeomAtom))
        if pots.scanType == 'dihedral':
            pots.finalGeomAtom[1] = pots.finalGeomAtom[1] \
                + '   Dihedral: ' \
                + str(etanumpy.get_dihedral(pots.atomList[0] + 1,
                      pots.atomList[1] + 1, pots.atomList[2] + 1,
                      pots.atomList[3] + 1, pots.finalGeomAtom))

    # Write last finalgeom to the SCAN xyz file

        etapots.write_final_geom_to_scan(pots)
        etapots.write_final_geom_to_scan_linear(pots)

    # Decide what the next step should be

        pots = etapots.scan_ts_logic(pots)

    # ## Implemented in scanTSLogic so that the scan file is updated
    # #Change the geometry for the next step
    # pots = etapots.bond_angle_dihedral(pots)
    # Delete the optimize files if they exist

        etaatom.del_file(pots.scratchdir + '/optenergies-'
                         + pots.jobname + '.txt')

    # Reset the configuration for the next scan step. This will hasten optimization if a smaller step size was taken

        pots.QMconfig = []
        pots.QMendlines = []
        for line in pots.QMconfigORIG:
            pots.QMconfig.append(line)
        pots.QMcomment = pots.QMcommentORIG
        for line in pots.QMendlinesORIG:
            pots.QMendlines.append(line)

  # #####################################################################
  # ## --- Output a new input file for the next step in the OPT/(TS-)SCAN --- ###

    if pots.QMtype == 'G09' or pots.QMtype == 'G03':

    # Move/convert files for the next computation

        etapots.rebuild_g0x_com(pots)
    elif pots.QMtype == 'Turbomole':

    # output the fixed geometry
    # outputFixUnfixGeom(unfixgeom)
    # Move/convert files for the next computation

        etapots.rebuild_turbomole_tmol(pots)
    elif pots.QMtype == 'Terachem':

    # Move/convert files for the next computation

        etapots.rebuild_terachem_xyz(pots)
    elif pots.QMtype.upper() == 'JAGUAR':

    # Move/convert files for the next computation

        etapots.rebuild_jaguar_in(pots)
    etapots.print_and_write(['################################################################################'
                            ], pots.verbose, pots.joblog)

  # ## --- RESET THE VARIABLES --- ###

    pots = etapots.reset_pots_variables(pots)
    if pots.test == 1:
        pots.kill = True
etapots.print_and_write(['====== END ======', ''], pots.verbose,
                        pots.joblog)

  # #####################################################################
  # ## END OF SCRIPT
  # #####################################################################
