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


### --- Build Gaussian .sge submit script for SGE --- ###

def submit_g09(basename, queue, nproc):

  # Start building up the SGE script

    printlines = [
        '#!/bin/bash\n',
        '\n',
        '#$ -S /bin/sh\n',
        '\n',
        '# Give a jobname. This option is used to  append a string to any jobname starting with\n'
            ,
        '# numbers instead of alphabets because SGE cannot accept jobname starting with\n'
            ,
        '# numbers.\n',
        '#$ -N G09-' + basename + '\n',
        '#\n',
        '#  The following items pertain to this script\n',
        '#  Use current working directory\n',
        '#$ -cwd\n',
        '#  input           = /dev/null\n',
        '#  output          = ' + basename + '.joblog\n',
        '#$ -o ' + basename + '.joblog\n',
        '#  error           = Merged with joblog\n',
        '#$ -j y\n',
        '\n',
        '# Resources\n',
        '#$ -pe smp ' + str(nproc) + '\n',
        '#$ -q ' + str(queue) + '\n',
        '#$ -l mem_free=' + memory + '\n',
        '#$ -l h_vmem=' + memory + '\n',
        '\n',
        '#  Email address to notify\n',
        '#$ -M ' + adminemail + '\n',
        '\n',
        '#  Notify at beginning and end of job\n',
        '#$ -m bea\n',
        '\n',
        '#  Job is not rerunable\n',
        '#$ -r n\n',
        '\n',
        '# Initialization for Gaussian09 serial execution\n',
        '  g09root=/export/apps\n',
        '  mkdir /state/partition1/$USER-' + basename + '-$JOB_ID\n',
        '  GAUSS_SCRDIR="/state/partition1/$USER-' + basename
            + '-$JOB_ID"\n',
        '  export g09root GAUSS_SCRDIR\n',
        '  . $g09root/g09/bsd/g09.profile\n',
        '\n',
        '  echo ""\n',
        '  echo "____________________________________________________________________________________"\n'
            ,
        '  echo "Gaussian09 SGE SUBMISSION SCRIPT"\n',
        '  echo "  Started by USER:      $USER"\n',
        '  echo "  JOB ID:               $JOB_ID"\n',
        '  echo "  APP directory:         $g09root/g09"\n',
        '  echo "  PWD directory:        $PWD"\n',
        '  echo "  SCRATCH directory:    $GAUSS_SCRDIR"\n',
        '\n',
        '# Run the Gaussian09 program\n',
        '  echo "  Started on:           " `hostname -s`\n',
        '  echo "  Started at:           " `date`\n',
        '\n',
        'SNIPPET!\n',
        '\n',
        '  echo "  Finished at:          " `date`\n',
        '  rm -rvf /state/partition1/$USER-' + basename + '-$JOB_ID\n',
        '\n',
        '\n',
        '',
        ]

  # ###############################################################################################
  # Gausssian setup

  # ###############################################################################################
  # Print out general data

  # ###############################################################################################
  # Program information below.
  # if snippet == "none":
  # printlines.append("  g09 " + basename + ".com\n")
  # ###############################################################################################
  # Finish the script

    return (printlines, memory)


### --- Build Gaussian .sge submit script for SGE --- ###

def submit_jaguar(basename, queue, nproc):

  # Start building up the SGE script

    printlines = []
    printlines.append('#!/bin/bash\n')
    printlines.append('\n')
    printlines.append('#$ -S /bin/sh\n')
    printlines.append('\n')
    printlines.append('# Give a jobname. This option is used to  append a string to any jobname starting with\n'
                      )
    printlines.append('# numbers instead of alphabets because SGE cannot accept jobname starting with\n'
                      )
    printlines.append('# numbers.\n')
    printlines.append('#$ -N Jaguar-' + basename + '\n')
    printlines.append('#\n')
    printlines.append('#  The following items pertain to this script\n')
    printlines.append('#  Use current working directory\n')
    printlines.append('#$ -cwd\n')
    printlines.append('#  input           = /dev/null\n')
    printlines.append('#  output          = ' + basename + '.joblog\n')
    printlines.append('#$ -o ' + basename + '.joblog\n')
    printlines.append('#  error           = Merged with joblog\n')
    printlines.append('#$ -j y\n')
    printlines.append('\n')
    printlines.append('# Resources\n')
    printlines.append('#$ -pe smp ' + str(nproc) + '\n')
    printlines.append('#$ -q ' + str(queue) + '\n')
    printlines.append('#$ -l mem_free=' + memory + '\n')
    printlines.append('#$ -l h_vmem=' + memory + '\n')
    printlines.append('\n')
    printlines.append('#  Email address to notify\n')
    printlines.append('#$ -M ' + adminemail + '\n')
    printlines.append('\n')
    printlines.append('#  Notify at beginning and end of job\n')
    printlines.append('#$ -m bea\n')
    printlines.append('\n')
    printlines.append('#  Job is not rerunable\n')
    printlines.append('#$ -r n\n')
    printlines.append('\n')

  # ###############################################################################################
  # Jaguar setup
  # Initialization for Schrodinger Variables

    printlines.append('export SCHRODINGER="/export/apps/schrodinger2014-2"'
                      )
    printlines.append('path=(')
    printlines.append("	/export/apps/schrodinger2014-2")
    printlines.append("	/bin")
    printlines.append("	/sbin")
    printlines.append("	/usr/bin")
    printlines.append("	/usr/sbin")
    printlines.append("	/share/apps/openbabel-2.3.0/bin")
    printlines.append("	/opt/gridengine/bin/lx26-amd64")
    printlines.append("	$path")
    printlines.append(')')
    printlines.append('export path')
    printlines.append('export JAGUAR_SCRATCH=/state/partition1/$USER-'
                      + basename + '-$JOB_ID')
    printlines.append('mkdir $JAGUAR_SCRATCH')
    printlines.append('basedir=`pwd`')
    printlines.append('')

  # ###############################################################################################

    printlines.append('')
    printlines.append('echo ""')
    printlines.append('echo "____________________________________________________________________________________"'
                      )
    printlines.append('echo "JAGUAR 2014-2 SGE SUBMISSION SCRIPT"')
    printlines.append('echo "Started by USER:      $USER"')
    printlines.append('echo "JOB ID:               $JOB_ID')
    printlines.append('echo "PWD directory:        $PWD"')
    printlines.append('')
    printlines.append('# Run the JAGUAR program')
    printlines.append('echo "Started on:           " `hostname -s`')
    printlines.append('echo "Started at:           " `date`')
    printlines.append('')

  # ###############################################################################################
  # Setting the memory

    printlines.append('ulimit -s unlimited')
    printlines.append('ulimit -c unlimited')
    printlines.append('ulimit -a')
    printlines.append('')

  # ###############################################################################################
  # Make local directories on the node and copy all files to that directory

    printlines.append('echo ""')
    printlines.append('echo "Making temp directory and moving the files to the node."'
                      )
    printlines.append('mkdir /state/partition1/$USER-' + basename
                      + '-$JOB_ID')
    printlines.append('cp -rfv * /state/partition1/$USER-' + basename
                      + '-$JOB_ID')
    printlines.append('cd /state/partition1/$USER-' + basename
                      + '-$JOB_ID')
    printlines.append('rm ' + basename + '.log')
    printlines.append('#rm ' + basename + '.joblog')
    printlines.append('')

  # ###############################################################################################
  # This script checks for hanges in the log/out file and copies changes to the head node job directory

    printlines.append('touch lockfile')
    printlines.append('function slave2 {')
    printlines.append('while test -f lockfile')
    printlines.append('do')
    printlines.append("	sleep 15")
    printlines.append("	baselog=`wc -l " + basename + '.tmp`')
    printlines.append("	nodelog=`wc -l " + basename + '.out`')
    printlines.append("	if [[ $baselog != $nodelog ]]")
    printlines.append("	then")
    printlines.append("		cat BASE!.log > $basedir/" + basename + '.log')
    printlines.append("		cat BASE!.out > $basedir/" + basename + '.out')
    printlines.append("		cat BASE!.out > " + basename + '.tmp')
    printlines.append("	fi")
    printlines.append("	sleep 15")
    printlines.append('done')
    printlines.append('}')
    printlines.append('slave2 &')
    printlines.append('')

  # ###############################################################################################
  # Jaguar execution

    printlines.append('echo "Running JAGUAR!"')

  # printlines.append("/export/apps/schrodinger2014-2/jaguar run -TPP 6 -WAIT -HOST localhost " + basename + ".in")

    printlines.append('SNIPPET!\n')

  # ###############################################################################################

  # ###############################################################################################
  # Copy finished files back to the head node

    printlines.append('')
    printlines.append('echo ""')
    printlines.append('rm -rfv lockfile')
    printlines.append('echo "Copying everything back and removing the files from the node."'
                      )
    printlines.append('tar -czvf SCRATCH-' + basename
                      + '.tar.gz /state/partition1/$USER-' + basename
                      + '-$JOB_ID/*')
    printlines.append('mv -fv /state/partition1/$USER-' + basename
                      + '-$JOB_ID/SCRATCH-' + basename
                      + '.tar.gz $basedir')
    printlines.append('cat ' + basename + '.out > $basedir/' + basename
                      + '.out')
    printlines.append('cat ' + basename + '.log > $basedir/' + basename
                      + '.log')
    printlines.append('rm -rfv /state/partition1/$USER-' + basename
                      + '-$JOB_ID')
    printlines.append('cd $basedir')

    return (printlines, memory)


### --- Build Turbomole .sge submit script for SGE --- ###

def submit_tm(basename, queue, nproc):
    printlines = []

  # Start building up the SGE script

    printlines.append('#!/bin/bash\n')
    printlines.append('\n')
    printlines.append('#$ -S /bin/sh\n')
    printlines.append('\n')
    printlines.append('# Give a jobname. This option is used to  append a string to any jobname starting with\n'
                      )
    printlines.append('# numbers instead of alphabets because SGE cannot accept jobname starting with\n'
                      )
    printlines.append('# numbers.\n')
    printlines.append('#$ -N TM-' + basename + '\n')
    printlines.append('#\n')
    printlines.append('#  The following items pertain to this script\n')
    printlines.append('#  Use current working directory\n')
    printlines.append('#$ -cwd\n')
    printlines.append('#  input           = /dev/null\n')
    printlines.append('#  output          = ' + basename + '.joblog\n')
    printlines.append('#$ -o ' + basename + '.joblog\n')
    printlines.append('#  error           = Merged with joblog\n')
    printlines.append('#$ -j y\n')
    printlines.append('\n')
    printlines.append('# Resources\n')
    printlines.append('#$ -pe smp ' + str(nproc) + '\n')
    printlines.append('#$ -q ' + str(queue) + '\n')
    printlines.append('#$ -l mem_free=' + memory + '\n')
    printlines.append('#$ -l h_vmem=' + memory + '\n')
    printlines.append('\n')
    printlines.append('#  Email address to notify\n')
    printlines.append('#$ -M ' + adminemail + '\n')
    printlines.append('\n')
    printlines.append('#  Notify at beginning and end of job\n')
    printlines.append('#$ -m bea\n')
    printlines.append('\n')
    printlines.append('#  Job is not rerunable\n')
    printlines.append('#$ -r n\n')
    printlines.append('\n')

  # ###############################################################################################
  # Turbomole setup

    printlines.append('# Initialization for Turbomole execution\n')
    printlines.append('export PARA_ARCH=SMP\n')
    printlines.append('export PARNODES=' + nproc + '\n')
    printlines.append('export TURBOTMPDIR=/state/partition1/$USER-'
                      + basename + '-$JOB_ID-TEMP\n')
    printlines.append('export TURBODIR=/share/apps/COSMOlogic11/TURBOMOLE\n'
                      )
    printlines.append('export TURBOWINDIR=/share/apps/COSMOlogic11/TURBOMOLE/bin/x86_64-unknown-linux-gnu\n'
                      )
    printlines.append('export PATH=/bin:/sbin:/usr/bin:/usr/sbin:/share/apps/openbabel-2.3.0/bin:/opt/gridengine/bin/lx26-amd64\n'
                      )
    printlines.append('export PATH=$TURBODIR/scripts:/share/apps/psw/local/Turbomole/6.XX/aux:$PATH\n'
                      )
    printlines.append('export PATH=$TURBODIR/bin/`sysname`:$PATH\n')
    printlines.append('mkdir $TURBOTMPDIR\n')
    printlines.append('\n')
    printlines.append('basedir=`pwd`\n')

  # ###############################################################################################

    printlines.append('\n')
    printlines.append('  echo ""\n')
    printlines.append('  echo "____________________________________________________________________________________"\n'
                      )
    printlines.append('  echo "Turbomole SGE SUBMISSION SCRIPT"\n')
    printlines.append('  echo "  Started by USER:      $USER"\n')
    printlines.append('  echo "  JOB ID:               $JOB_ID"\n')
    printlines.append('  echo "  APP directory:         $g09root/g09"\n'
                      )
    printlines.append('  echo "  PWD directory:        $PWD"\n')
    printlines.append('  echo "  SCRATCH directory:    $TURBOTMPDIR"\n')
    printlines.append('\n')
    printlines.append('# Run the Turbomole program\n')
    printlines.append('  echo "  Started on:           " `hostname -s`\n'
                      )
    printlines.append('  echo "  Started at:           " `date`\n')
    printlines.append('\n')

  # ###############################################################################################
  # Program information below.
  # snippetfile = open(snippet, "r")
  # for line in snippetfile:
  #  f.write(line)
  # snippetfile.close()

    printlines.append('SNIPPET!\n')

  # ###############################################################################################

    printlines.append('\n')
    printlines.append('  echo "  Finished at:          " `date`\n')
    printlines.append('  rm -rvf /state/partition1/$USER-' + basename
                      + '-$JOB_ID\n')
    printlines.append('\n')
    printlines.append('\n')
    printlines.append('')
    return (printlines, memory)


### --- Build Terachem .sge submit script for SGE --- ###

def submit_tc(
    basename,
    queue,
    nproc,
    snippet,
    ):

  # Start building up the SGE script

    printlines = []
    printlines.append('#!/bin/bash\n')
    printlines.append('\n')
    printlines.append('#$ -S /bin/sh\n')
    printlines.append('\n')
    printlines.append('# Give a jobname. This option is used to  append a string to any jobname starting with\n'
                      )
    printlines.append('# numbers instead of alphabets because SGE cannot accept jobname starting with\n'
                      )
    printlines.append('# numbers.\n')
    printlines.append('#$ -N TC-' + basename + '\n')
    printlines.append('#\n')
    printlines.append('#  The following items pertain to this script\n')
    printlines.append('#  Use current working directory\n')
    printlines.append('#$ -cwd\n')
    printlines.append('#  input           = /dev/null\n')
    printlines.append('#  output          = ' + basename + '.joblog\n')
    printlines.append('#$ -o ' + basename + '.joblog\n')
    printlines.append('#  error           = Merged with joblog\n')
    printlines.append('#$ -j y\n')
    printlines.append('\n')
    printlines.append('# Resources\n')
    printlines.append('#$ -pe smp ' + str(nproc) + '\n')
    printlines.append('#$ -q ' + str(queue) + '\n')
    printlines.append('#$ -l mem_free=' + memory + '\n')
    printlines.append('#$ -l h_vmem=' + memory + '\n')
    printlines.append('\n')
    printlines.append('#  Email address to notify\n')
    printlines.append('#$ -M ' + adminemail + '\n')
    printlines.append('\n')
    printlines.append('#  Notify at beginning and end of job\n')
    printlines.append('#$ -m bea\n')
    printlines.append('\n')
    printlines.append('#  Job is not rerunable\n')
    printlines.append('#$ -r n\n')
    printlines.append('\n')

  # ###############################################################################################
  # Terachem setup

    printlines.append('# Put your program run command below this line.\n'
                      )
    printlines.append('export LD_LIBRARY_PATH=/usr/local/cuda/lib:/usr/local/cuda/lib64:/share/apps/TeraChem/stdlib:/opt/gridengine/lib/lx26-amd64:/opt/openmpi/lib\n'
                      )
    printlines.append('export TeraChem=/share/apps/TeraChem\n')
    printlines.append('export NBOEXE=/share/apps/TeraChem/nbo6.exe\n')
    printlines.append('export OMP_NUM_THREADS=4\n')

  # ###############################################################################################

    printlines.append('echo "____________________________________________________________________________________"\n'
                      )
    printlines.append('echo "TERACHEM 1.5K SGE SUBMISSION SCRIPT"\n')
    printlines.append('echo "Started by USER:      $USER"\n')
    printlines.append('echo "JOB ID: $JOB_ID"\n')
    printlines.append('echo "PWD directory:        $PWD"\n')
    printlines.append('echo ""\n')
    printlines.append('echo "Started on:           " `hostname -s`\n')
    printlines.append('echo "Started at:           " `date`\n')
    printlines.append('echo "=============================================================================================================================="\n'
                      )
    printlines.append('hostname\n')
    printlines.append('date \n')

  # ###############################################################################################
  # Program information bel\now.

    printlines.append('hostname\n')

  # printlines.append('/share/psw/TeraChem/hidden/RUN-TC.sh $base.tc >> $base.joblog\n')'

    printlines.append('SNIPPET!\n')

  # ###############################################################################################

    printlines.append('rm core.*\n')
    printlines.append('echo "=============================================================================================================================="\n'
                      )
    printlines.append('date \n')

    return (printlines, memory)


