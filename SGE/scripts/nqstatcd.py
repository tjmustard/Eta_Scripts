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

jobnum = str(sys.argv[1])
joblist = []
print ''
if jobnum == 'all':

    qstat = os.popen('qstat').read().splitlines()
    for i in range(2, len(qstat)):
        joblist.append(qstat[i].split()[0])
    for i in joblist:
        SGEWD = os.popen('qstat -explain c -j ' + i
                         + " |grep sge_o_workdir |sed -e 's/sge_o_workdir:              //'"
                         ).read()
        SGEFN = os.popen('qstat -explain c -j ' + i
                         + " |grep stdout_path_list |sed -e 's/stdout_path_list:           //'|sed -e 's/.joblog//'"
                         ).read()
        print 'Job Name:      ' + str(SGEFN).split(':')[-1].strip()
        print 'Job Directory: ' + str(SGEWD).strip()
elif jobnum == '':
    print 'This script changes to the location of the job ID you have supplied.'
    print 'Please supply a job ID.'
    print 'For example: qstatw 12324'
else:
    print 'Parsing the physical location of JOB ID: ' + jobnum
    SGEWD = os.popen('qstat -explain c -j ' + jobnum
                     + " |grep sge_o_workdir |sed -e 's/sge_o_workdir:              //'"
                     ).read()
    SGEFN = os.popen('qstat -explain c -j ' + jobnum
                     + " |grep stdout_path_list |sed -e 's/stdout_path_list:           //'|sed -e 's/.joblog//'"
                     ).read()
    print 'Job Name:      ' + str(SGEFN).split(':')[-1].strip()
    print 'Job Directory: ' + str(SGEWD).strip()

