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
from decimal import *
import shutil
import etaatom
import etamap
from decimal import *
from Tkinter import *
import tkFileDialog

# Create the window

main = Tk()

# Modify the window and make it 100x100 pixels

main.title('nMAP GUI')
main.geometry('450x400')


def get_rel_path():
    filePath = tkFileDialog.askdirectory()
    rel_path = os.path.relpath(filePath)
    return rel_path


def set_file_path_pb():
    rel_path = get_rel_path()
    pe.delete(0, 'end')
    pe.insert(0, rel_path)
    return


def set_file_path_cb():
    rel_path = get_rel_path()
    ce.delete(0, 'end')
    ce.insert(0, rel_path)
    return


def set_file_path_ob():
    rel_path = get_rel_path()
    oe.delete(0, 'end')
    oe.insert(0, rel_path)
    return


#  mapinput.parentFolder = ''

pl = Label(main, text='Parent Folder:')
pl.grid(row=1, column=0, sticky=W)
pe = Entry(main)
pe.grid(row=1, column=1)
pb = Button(main, text='Browse', command=set_file_path_pb)
pb.grid(row=1, column=2)

#  mapinput.childFolder = ''

cl = Label(main, text='Children Folder:')
cl.grid(row=2, column=0, sticky=W)
ce = Entry(main)
ce.grid(row=2, column=1)
cb = Button(main, text='Browse', command=set_file_path_cb)
cb.grid(row=2, column=2)

#  mapinput.outputFolder = ''

ol = Label(main, text='Output Folder:')
ol.grid(row=3, column=0, sticky=W)
oe = Entry(main)
oe.grid(row=3, column=1)
ob = Button(main, text='Browse', command=set_file_path_ob)
ob.grid(row=3, column=2)

#  mapinput.alignParent = []

pal = Label(main, text='Parent Alignment Atoms:')
pal.grid(row=4, column=0, sticky=W)
pae = Entry(main)
pae.grid(row=4, column=1)

#  mapinput.alignChild = []

cal = Label(main, text='Children Alignment Atoms:')
cal.grid(row=5, column=0, sticky=W)
cae = Entry(main)
cae.grid(row=5, column=1)

#  mapinput.alignWeights = []

awl = Label(main, text='Alignment Weights:')
awl.grid(row=6, column=0, sticky=W)
awe = Entry(main)
awe.grid(row=6, column=1)

#  mapinput.parentDelete = []

pdl = Label(main, text='Parent Delete Atoms:')
pdl.grid(row=7, column=0, sticky=W)
pde = Entry(main)
pde.grid(row=7, column=1)

#  mapinput.childDelete = []

cdl = Label(main, text='Children Delete Atoms:')
cdl.grid(row=8, column=0, sticky=W)
cde = Entry(main)
cde.grid(row=8, column=1)

#  mapinput.quatfit = "quatfit"

ql = Label(main, text='Quatfit:')
ql.grid(row=9, column=0, sticky=W)
qe = Entry(main)
qe.grid(row=9, column=1)
qe.insert(0, 'quatfit')

#  mapinput.baseName = ''

bl = Label(main, text='Base Name:')
bl.grid(row=10, column=0, sticky=W)
be = Entry(main)
be.grid(row=10, column=1)

#  mapinput.endName = ""

el = Label(main, text='End Name:')
el.grid(row=11, column=0, sticky=W)
ee = Entry(main)
ee.grid(row=11, column=1)


##  triangle
# triangle = IntVar()
# tc = Checkbutton(main, text="Triangle", variable = triangle, onvalue = 1, offvalue = 0)
# tc.grid(row=12, column=1, sticky=W)

def grab_arguments():

  # ## --- Arguments --- ###

    program = 'nMAP.py'
    ifile = ''
    ofile = ''
    mapinput = etamap.MAPArguments()
    mapinput.debug = 0
    writeTemplateMap = False
    triangle = False

    if writeTemplateMap == 1:
        etamap.write_template()
        sys.exit(0)

    mapinput.parentFolder = pe.get()
    mapinput.childFolder = ce.get()
    mapinput.outputFolder = oe.get()
    alignParentTemp = pae.get()
    alignChildTemp = cae.get()
    alignWeightsTemp = awe.get()
    for i in alignParentTemp.split():
        mapinput.alignParent.append(i)
    for i in alignChildTemp.split():
        mapinput.alignChild.append(i)
    for i in alignWeightsTemp.split():
        mapinput.alignWeights.append(i)

  # mapinput.parentDelete = pde.get()

    parentDeleteTemp = pde.get()
    for i in parentDeleteTemp.split():
        mapinput.parentDeleteTemp.append(i)

  # mapinput.childDelete = []

    childDeleteTemp = cde.get()
    for i in childDeleteTemp.split():
        mapinput.childDeleteTemp.append(i)
    mapinput.quatfit = qe.get()
    mapinput.baseName = be.get()
    mapinput.endName = ee.get()
    return mapinput


def runnMAP():

    mapinput = grab_arguments()

  # mapinput = etamap.parse_input_file(ifile)

  # ## --- Expand upon the parentDeleteTemp list if
  # there are certain characters to denote multiple atoms.
  # ##

    mapinput = etamap.expand_temp_lists(mapinput)

  # ## --- Make the output folder --- ###

    etaatom.make_dir(mapinput.outputFolder)
    etaatom.make_dir(mapinput.outputFolder + '/temp')

  # ## --- Make the align files for quatfit. --- ###

    etamap.make_align_file(mapinput)

  # ## --- Align the parent files before aligning the substrates to them --- ###
  # ## Purely for aesthetics and does not change the final outcome

    etamap.align_parents(mapinput)

    triangle = False
    etamap.align_children(mapinput, triangle)
    return


def saveConfig():

    mapinput = grab_arguments()
    ofile = tkFileDialog.asksaveasfilename()
    f = open(ofile, 'w')
    f.write('#Input file' + '\n')
    f.write('parentFolder = ' + mapinput.parentFolder + '\n')
    f.write('#Library folder to use' + '\n')
    f.write('childFolder = ' + mapinput.childFolder + '\n')
    f.write('#Output folder' + '\n')
    f.write('outputFolder = ' + mapinput.outputFolder + '\n')
    f.write('#Base name for files' + '\n')
    f.write('baseName = ' + mapinput.baseName + '\n')
    f.write('#End name for files' + '\n')
    f.write('endName = ' + mapinput.endName + '\n')
    f.write('#Quatfit location' + '\n')
    f.write('quatfit = ' + mapinput.quatfit + '\n')
    f.write('#List of atoms to align to' + '\n')
    alignptemp = ''
    alignctemp = ''
    alignwtemp = ''
    alignpdtemp = ''
    aligncdtemp = ''
    for i in mapinput.alignParent:
        alignptemp += i + ' '
    for i in mapinput.alignChild:
        alignctemp += i + ' '
    for i in mapinput.alignWeights:
        alignwtemp += i + ' '
    for i in mapinput.parentDeleteTemp:
        alignpdtemp += i + ' '
    for i in mapinput.childDeleteTemp:
        aligncdtemp += i + ' '
    f.write('parent = ' + alignptemp + '\n')
    f.write('child = ' + alignctemp + '\n')
    f.write('weights = ' + alignwtemp + '\n')
    f.write('#List of atoms to delete from the input (parent) file'
            + '\n')
    f.write('deleteParent = ' + alignpdtemp + '\n')
    f.write('deleteChild = ' + aligncdtemp + '\n')
    f.close()
    return


def loadConfig():
    ofile = tkFileDialog.askopenfilename(defaultextension='map',
            initialdir='./')
    mapinput = etamap.MAPArguments()
    mapinput = etamap.parse_input_file(ofile)
    pe.delete(0, 'end')
    ce.delete(0, 'end')
    oe.delete(0, 'end')
    pae.delete(0, 'end')
    cae.delete(0, 'end')
    awe.delete(0, 'end')
    pde.delete(0, 'end')
    cde.delete(0, 'end')
    qe.delete(0, 'end')
    be.delete(0, 'end')
    ee.delete(0, 'end')
    pe.insert(0, mapinput.parentFolder)
    ce.insert(0, mapinput.childFolder)
    oe.insert(0, mapinput.outputFolder)
    alignptemp = ''
    alignctemp = ''
    alignwtemp = ''
    alignpdtemp = ''
    aligncdtemp = ''
    for i in mapinput.alignParent:
        alignptemp += i + ' '
    for i in mapinput.alignChild:
        alignctemp += i + ' '
    for i in mapinput.alignWeights:
        alignwtemp += i + ' '
    for i in mapinput.parentDeleteTemp:
        alignpdtemp += i + ' '
    for i in mapinput.childDeleteTemp:
        aligncdtemp += i + ' '
    pae.insert(0, alignptemp)
    cae.insert(0, alignctemp)
    awe.insert(0, alignwtemp)
    pde.insert(0, alignpdtemp)
    cde.insert(0, aligncdtemp)
    qe.insert(0, mapinput.quatfit)
    be.insert(0, mapinput.baseName)
    ee.insert(0, mapinput.endName)

    return


Run = Button(main, text='Run nMAP', command=runnMAP)
Run.grid(row=14, column=0)

Save = Button(main, text='Save Configuration', command=saveConfig)
Save.grid(row=13, column=0)

Load = Button(main, text='Load Configuration', command=loadConfig)
Load.grid(row=13, column=1)

# Start off the gui event loop

main.mainloop()

#######################################################################
#### END OF SCRIPT
########################################################################
