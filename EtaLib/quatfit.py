#! /usr/bin/env python3

##############################################################################
# The program to superimpose atoms of two molecules by quaternion method
#
# David J. Heisterberg
# The Ohio Supercomputer Center
# 1224 Kinnear Rd.
# Columbus, OH  43212-1163
# (614)292-6036
# djh@ccl.net    djh@ohstpy.bitnet    ohstpy::djh
#
# Translated to C from fitest.f program and interfaced with Xmol program
# by Jan Labanowski,  jkl@ccl.net   jkl@ohstpy.bitnet   ohstpy::jkl
#
# Translated to python from quatfit.c
# by Thomas J. L. Mustard, mustardt@onid.orst.edu
#
# Copyright: Ohio Supercomputer Center, David J. Heisterberg, 1990.
# The program can be copied and distributed freely, provided that
# this copyright in not removed. You may acknowledge the use of the
# program in published material as:
# David J. Heisterberg, 1990, unpublished results.
#
#
# #include <iostream>
# #include <fstream>
# #include <cstdlib>
# #include <vector>
# #include <string>
# #include <cassert>
# #include <cmath>
# #include <algorithm>
# #include <sstream>
# #using namespace std
#
# #define MAXPOINTS     400
# #define MAXLINELEN    250
# #define OPTIONS       "r:f:p:o:s:"
#
#  options
#  -r refmol    reference molecule xmol file. If this option is not given, the
#               information is read from standard input.
#
#  -f fitmol    fitted molecule input xmol file If this option is not given, the
#               information is read from standard input.
#
#  -p pairs     file with the list of fitted atom pairs and weights. If this
#               option is not specified, pairs and weights are taken from
#               stdin. If file name "none" is used (i.e. -p none), atoms of the
#               fitted molecule are fitted to atoms of the reference
#               molecule with the same number, and all weights are assumed 1.0.
#               If molecules do not have the same number of atoms, the
#               smaller number of atoms is fitted.
#
#  -o outmol    output file for transformed fitted molecule. If this option
#               is not given, the information is written to standard output.
#
#  -s statfile  file with fit statistics (distances between fitted atoms, etc).
#               If this option is not given, the information is written to
#               standard output.
#
#  If any files are read from stdin, the order is: refmol, fitmol, pairs.
#  If any files are written to stdout, the order is: outmol, statfile.
#
#The file formats are:
# The refmol, fitmol and outmol files are in the XYZ format used with xmol
# program from Minnesota Supercomputer Institute. The format is:
#   1st line: number of atoms
#   2nd line: title
#   3rd and next lines have the format depending on the kind of information:
#       T  X  Y  Z                (total of 4 columns)
#       T  X  Y  Z  C             (total of 5 columns)
#       T  X  Y  Z  Mx My Mz      (total of 7 columns)
#       T  X  Y  Z  C Mx My Mz    (total of 8 columns)
#     where T is atom type (usually element symbol), X, Y, Z are cartesian
#     coordinates of an atom in Angstroms, C is atomic charge, and Mx, My, Mz
#     are normal modes.
#
# The pairs file format is:
#   1st line: number of pairs
#   2nd and next lines:
#       Ar   Af    W
#     where Ar is the atom number of the reference molecule, Af is the atom
#     number of fitted molecule, w is the statistical weight. Weights W are
#     related to the square of expected deviation "sigma" between the reference
#     and fitted molecule atoms and allow to make fit of some atom pairs more
#     tight. W is proportional to 1/sigma^2. The larger the weight, the more
#     tight will be the resulting fit for the given pair.
#
# The statfile lists results of the fit with explanation.
#
# there was a typo in the formula for RMS. In the part calculating the
# RMS, there was a line:
#    wnorm += s
# while it should be:
#    wnorm += s*s
# Martin Lema mlema[-at-]unq.edu.ar was kind to find it. THANKS!!! jkl.
# It did not affect fitting too much. But the RMS value was not right when
# weights were different than 1.
#
####################################################################################

import os
from sys import *
import math
import sys


###==================================================================================================================###
### --- Atom class --- ###
class Atom(object):
  def __init__(self):
    #Hold the string for the element
    self.e = ''
    #Hold the int for the element
    self.en = int(0)
    #Hold the float for the x coordinate
    self.x = float(0)
    #Hold the float for the y coordinate
    self.y = float(0)
    #Hold the float for the z coordinate
    self.z = float(0)
    #Hold an int/float for the charge of this atom
    self.charge = ''
    #Hold the float for the x coordinate
    self.mx = ''
    #Hold the float for the y coordinate
    self.my = ''
    #Hold the float for the z coordinate
    self.mz = ''

### --- Parse the input XYZ file and convert it into a Atom Class list of lists --- ###
def parseXYZ(ifile):
  ### --- Open parent file, read it in and close it --- ###
  f = open(ifile, "r")
  ifilelist = f.readlines()
  f.close()
  #Find out the final length of the ifilelol
  ifilelength = len(ifilelist)
  #Fill the ifilelol with 0's so that there is a place to put the atomic data
  ifilelol = [0] * (ifilelength - 2)
  #### --- Input/parse the input file into a list of lists called ifilelol --- ###
  for i in range(2, ifilelength):
    line = ifilelist[i].rstrip().split()
    ifilelol[i-2] = Atom()
    ifilelol[i-2].e = line[0]
    ifilelol[i-2].x = float(line[1])
    ifilelol[i-2].y = float(line[2])
    ifilelol[i-2].z = float(line[3])
    if len(line) >= 5:
      ifilelol[i-2].charge = float(line[4])
    if len(line) == 8:
      ifilelol[i-2].mx = float(line[5])
      ifilelol[i-2].my = float(line[6])
      ifilelol[i-2].mz = float(line[7])
  #If the XYZ file is longer or shorter than the list of data should be, warn the user.
  if int(ifilelist[0]) != len(ifilelol):
    print("Your file is the wrong length!")
    print("Make sure you have the right number of atoms.")
    exit(0)
  return ifilelol

### --- Output the structural data in an XYZ format
def outputXYZ(ofile, ofilelol, printlevel):
  #Open the file for writing
  f = open(ofile, "w+")
  f.write(str(len(ofilelol)) + "\n")
  f.write(ofile + "\n")
  #iterate through the whole ifilelol and write each line in XYZ format to the file
  for i in range(len(ofilelol)):
    #if the instance is in Atom form print the 'element X-Coord Y-Coord Z-Coord'
    if isinstance(ofilelol[i], Atom):
      if printlevel == 1:
        line = ofilelol[i].e + "  " + str("{:.6f}".format(ofilelol[i].x)) + "  " + str("{:.6f}".format(ofilelol[i].y)) + "  " + str("{:.6f}".format(ofilelol[i].z))
      elif printlevel == 2:
        line = ofilelol[i].e + "  " + str("{:.6f}".format(ofilelol[i].x)) + "  " + str("{:.6f}".format(ofilelol[i].y)) + "  " + str("{:.6f}".format(ofilelol[i].z)) + "  " + str("{:.6f}".format(ofilelol[i].charge))
      elif printlevel == 3:
        line = ofilelol[i].e + "  " + str("{:.6f}".format(ofilelol[i].x)) + "  " + str("{:.6f}".format(ofilelol[i].y)) + "  " + str("{:.6f}".format(ofilelol[i].z)) + "  " + str("{:.6f}".format(ofilelol[i].charge)) + "  " + str("{:.6f}".format(ofilelol[i].mx)) + "  " + str("{:.6f}".format(ofilelol[i].my)) + "  " + str("{:.6f}".format(ofilelol[i].mz))
    f.write(line + "\n")
  f.close()
  return

### --- Parse the pairs file and return both the pairs and the weights --- ###
def parsePairs(pairsfile):
  ### --- Open parent file, read it in and close it --- ###
  f = open(pairsfile, "r")
  pairsfilelist = f.readlines()
  f.close()
  #Find out the final length of the ifilelol
  if int(pairsfilelist[0]) != len(pairsfilelist) - 1:
    print("The pairs file is the incorrect length.\nIt should be " \
          + str(int(pairsfilelist[0])) + " pairs long and it is " + str(len(pairsfilelist) - 1) + " pairs long.")
  else:
    pairs = []
    weights=[]
    for i in range(len(pairsfilelist)):
      if i == 0:
        pairnum = int(pairsfilelist[i])
      elif i >= 1:
        line = pairsfilelist[i].strip().split()
        pairs.append([int(line[0]), int(line[1])])
        weights.append(int(line[2]))
  return pairs, weights

### --- Create the reference XYZ files for future alignment --- ###
def createRefGeom(reflol, fitlol, pairs):
  #Create ref_xyz and fit_xyz
  ref_xyz = [0] * len(pairs)
  fit_xyz = [0] * len(pairs)
  #Copy the pair atoms from the original molecule into the reference xyz file
  for i in range(len(pairs)):
    #Copy the ref pairs
    ref_xyz[i] = Atom()
    ref_xyz[i].e = reflol[pairs[i][0]-1].e
    ref_xyz[i].x = reflol[pairs[i][0]-1].x
    ref_xyz[i].y = reflol[pairs[i][0]-1].y
    ref_xyz[i].z = reflol[pairs[i][0]-1].z
    ref_xyz[i].charge = reflol[pairs[i][0]-1].charge
    ref_xyz[i].mx = reflol[pairs[i][0]-1].mx
    ref_xyz[i].my = reflol[pairs[i][0]-1].my
    ref_xyz[i].mz = reflol[pairs[i][0]-1].mz
    #Copy the fit pairs
    fit_xyz[i] = Atom()
    fit_xyz[i].e = fitlol[pairs[i][1]-1].e
    fit_xyz[i].x = fitlol[pairs[i][1]-1].x
    fit_xyz[i].y = fitlol[pairs[i][1]-1].y
    fit_xyz[i].z = fitlol[pairs[i][1]-1].z
    fit_xyz[i].charge = fitlol[pairs[i][1]-1].charge
    fit_xyz[i].mx = fitlol[pairs[i][1]-1].mx
    fit_xyz[i].my = fitlol[pairs[i][1]-1].my
    fit_xyz[i].mz = fitlol[pairs[i][1]-1].mz
  return ref_xyz, fit_xyz

### --- Print the molecule class object --- ###
def printAtom(ofilelol):
  for i in range(len(ofilelol)):
    print(ofilelol[i].e + "  " + str("{:.6f}".format(ofilelol[i].x)) + "  " + str("{:.6f}".format(ofilelol[i].y)) + "  " + str("{:.6f}".format(ofilelol[i].z)) + "  " + str("{:.6f}".format(ofilelol[i].charge)) + "  " + str("{:.6f}".format(ofilelol[i].mx)) + "  " + str("{:.6f}".format(ofilelol[i].my)) + "  " + str("{:.6f}".format(ofilelol[i].mz)))
  return

### --- center the coordinates, or translate them to some xyz --- ###
def center(filelol, weights, centerswitch, centerxyz):
  """====================================================================
  CENTER
   center or translate a molecule.
   atomnum (n) - number of atoms
   filelol (x) - on input  - original xyz coordinates of a molecule
       on output - moved xyz coordinates (see io for modes).

   weights (w) - if centerswitch=1, weights of atoms
       if centerswitch=2 or 3, unused

   centerswitch (io) - 1 weighted geometric center of the molecule will be at (0,0,0)
        2 molecule will be moved by a vector -center (i.e., components of a vector center
          will be subtracted from atom coordinates).
        3 molecule will be moved by a vector +center (i.e., components of a vector center
          will be added atom coordinates).

   centerxyz (o) - if centerswitch=1, output, center of original coordinates
       if centerswitch=2, input, vector center will be subtracted from atomic coordinates
       if centerswitch=3, input, vector center will be added to atomic coordinates

  ====================================================================="""
  wnorm = float(0.00)
  modif = float(0.00)
  #int i
  if centerswitch == 2:
    modif = -1.0
  elif centerswitch == 3:
    modif = 1.0
  else:
    modif = -1.0
    centerxyz[0] = float(0.0)
    centerxyz[1] = float(0.0)
    centerxyz[2] = float(0.0)
    wnorm = 0.0
    for i in range(len(filelol)):
      centerxyz[0] += filelol[i].x * math.sqrt(weights[i])
      centerxyz[1] += filelol[i].y * math.sqrt(weights[i])
      centerxyz[2] += filelol[i].z * math.sqrt(weights[i])
      wnorm += math.sqrt(weights[i])
    centerxyz[0] = centerxyz[0] / wnorm
    centerxyz[1] = centerxyz[1] / wnorm
    centerxyz[2] = centerxyz[2] / wnorm
  for i in range(len(filelol)):
    filelol[i].x += modif*centerxyz[0]
    filelol[i].y += modif*centerxyz[1]
    filelol[i].z += modif*centerxyz[2]
  return centerxyz, filelol


### --- ROTMOL --- ###
def rotmol(filelol, rotmat):
  """
  ROTMOL
  rotate a molecule
  n - number of atoms
  filelol (x) - input coordinates
  filelol (y) - rotated coordinates y = u * x
  rotmat (u) - left rotation matrix
  """
  yx = float(0.0)
  yy = float(0.0)
  yz = float(0.0)

  for i in range(len(filelol)):
    yx = rotmat[0][0] * filelol[i].x + rotmat[1][0] * filelol[i].y + rotmat[2][0] * filelol[i].z
    yy = rotmat[0][1] * filelol[i].x + rotmat[1][1] * filelol[i].y + rotmat[2][1] * filelol[i].z
    yz = rotmat[0][2] * filelol[i].x + rotmat[1][2] * filelol[i].y + rotmat[2][2] * filelol[i].z
    filelol[i].x = yx  #x
    filelol[i].y = yy  #y
    filelol[i].z = yz  #z
  return filelol


### --- JACOBI --- ###
def jacobi(matrix, maxsweeps):
  """
  JACOBI
  Jacobi diagonalizer with sorted output. It is only good for 4x4 matrices.
  (was too lazy to do pointers...)
  matrix (a) - input: matrix to diagonalize
  eigenvect (v) - output: eigenvectors
  eigenval (d) - output: eigenvalues
  maxsweeps (nrot) - input: maximum number of sweeps
  """

  eigenvect = [[float(0.0) for x in range(4)] for x in range(4)]
  eigenval = [float(0.0) for x in range(4)]
  onorm = float(0.0)
  dnorm = float(0.0)
  b = float(0.0)
  dma = float(0.0)
  q = float(0.0)
  t = float(0.0)
  c = float(0.0)
  s = float(0.0)
  atemp = float(0.0)
  vtemp = float(0.0)
  dtemp = float(0.0)

  for j in range(4):
    #for i in range(4):
    #  eigenvect[i][j] = 0.0
    eigenvect[j][j] = 1.0
    eigenval[j] = matrix[j][j]

  for m in range(maxsweeps):
    dnorm = 0.0
    onorm = 0.0
    for j in range(4):
      dnorm += math.fabs(eigenval[j])
      for i in range(j):
        onorm += math.fabs(matrix[i][j])
    if onorm/dnorm <= 1.0e-12: break  #goto Exit_now;
    for j in range(1, 4):
      for i in range(j):
        b = matrix[i][j]
        if math.fabs(b) > 0.0:
          dma = eigenval[j] - eigenval[i]
          if (math.fabs(dma) + math.fabs(b)) <=  math.fabs(dma):
            t = b / dma
          else:
            q = 0.5 * dma / b
            t = 1.0/(math.fabs(q) + math.sqrt(1.0+q*q))
            if q < 0.0:
              t = -t
          c = 1.0/math.sqrt(t * t + 1.0)
          s = t * c
          matrix[i][j] = 0.0
          for k in range(i):
            atemp = c * matrix[k][i] - s * matrix[k][j]
            matrix[k][j] = s * matrix[k][i] + c * matrix[k][j]
            matrix[k][i] = atemp
          for k in range(i+1, j):
            atemp = c * matrix[i][k] - s * matrix[k][j]
            matrix[k][j] = s * matrix[i][k] + c * matrix[k][j]
            matrix[i][k] = atemp
          for k in range(j+1, 4):
            atemp = c * matrix[i][k] - s * matrix[j][k]
            matrix[j][k] = s * matrix[i][k] + c * matrix[j][k]
            matrix[i][k] = atemp
          for k in range(4):
            vtemp = c * eigenvect[k][i] - s * eigenvect[k][j]
            eigenvect[k][j] = s * eigenvect[k][i] + c * eigenvect[k][j]
            eigenvect[k][i] = vtemp
          dtemp = c*c*eigenval[i] + s*s*eigenval[j] - 2.0*c*s*b
          eigenval[j] = s*s*eigenval[i] + c*c*eigenval[j] +  2.0*c*s*b
          eigenval[i] = dtemp

  maxsweeps = m

  for j in range(3):
    k = j
    dtemp = eigenval[k]
    for i in range(j+1, 4):
      if eigenval[i] < dtemp:
        k = i
        dtemp = eigenval[k]

    if k > j:
      eigenval[k] = eigenval[j]
      eigenval[j] = dtemp
      for i in range(4):
        dtemp = eigenvect[i][k]
        eigenvect[i][k] = eigenvect[i][j]
        eigenvect[i][j] = dtemp

  return eigenvect, eigenval, maxsweeps



### --- Q2MAT --- ###
def q2mat(quaternion):
  """
  Q2MAT
  Generate a left rotation matrix from a normalized quaternion

  INPUT
    quaternion (q)      - normalized quaternion

  OUTPUT
    rotmat (u)      - the rotation matrix
  """
  rotmat = [[float(0.0) for x in range(3)] for x in range(3)]
  rotmat[0][0] = quaternion[0]*quaternion[0] + quaternion[1]*quaternion[1] - quaternion[2]*quaternion[2] - quaternion[3]*quaternion[3]
  rotmat[1][0] = 2.0 * (quaternion[1] * quaternion[2] - quaternion[0] * quaternion[3])
  rotmat[2][0] = 2.0 * (quaternion[1] * quaternion[3] + quaternion[0] * quaternion[2])
  rotmat[0][1] = 2.0 * (quaternion[2] * quaternion[1] + quaternion[0] * quaternion[3])
  rotmat[1][1] = quaternion[0]*quaternion[0] - quaternion[1]*quaternion[1] + quaternion[2]*quaternion[2] - quaternion[3]*quaternion[3]
  rotmat[2][1] = 2.0 * (quaternion[2] * quaternion[3] - quaternion[0] * quaternion[1])
  rotmat[0][2] = 2.0 *(quaternion[3] * quaternion[1] - quaternion[0] * quaternion[2])
  rotmat[1][2] = 2.0 * (quaternion[3] * quaternion[2] + quaternion[0] * quaternion[1])
  rotmat[2][2] = quaternion[0]*quaternion[0] - quaternion[1]*quaternion[1] - quaternion[2]*quaternion[2] + quaternion[3]*quaternion[3]
  return rotmat


### --- QTRFIT --- ###
def qtrfit(fit_xyz, ref_xyz, weights, maxsweeps):
  """
   QTRFIT
   Find the quaternion, q,[and left rotation matrix, u] that minimizes

     |qTXq - Y| ^ 2  [|uX - Y| ^ 2]

   This is equivalent to maximizing Re (qTXTqY).

   This is equivalent to finding the largest eigenvalue and corresponding
   eigenvector of the matrix

   [A2   AUx  AUy  AUz ]
   [AUx  Ux2  UxUy UzUx]
   [AUy  UxUy Uy2  UyUz]
   [AUz  UzUx UyUz Uz2 ]

   where

     A2   = Xx Yx + Xy Yy + Xz Yz
     Ux2  = Xx Yx - Xy Yy - Xz Yz
     Uy2  = Xy Yy - Xz Yz - Xx Yx
     Uz2  = Xz Yz - Xx Yx - Xy Yy
     AUx  = Xz Yy - Xy Yz
     AUy  = Xx Yz - Xz Yx
     AUz  = Xy Yx - Xx Yy
     UxUy = Xx Yy + Xy Yx
     UyUz = Xy Yz + Xz Yy
     UzUx = Xz Yx + Xx Yz

   The left rotation matrix, u, is obtained from q by

     u = qT1q

   INPUT
     n      - number of points
     fit_xyz (x)      - fitted molecule coordinates
     ref_xyz (y)      - reference molecule coordinates
     weights (w)      - weights

   OUTPUT
     quaternion (q)      - the best-fit quaternion
     rotmat (u)      - the best-fit left rotation matrix
     maxsweeps (nr)     - max number of jacobi sweeps

  """


  #Create variables/lists/matrixes
  matrix = [[float(0.0) for x in range(4)] for x in range(4)]  #double c[4][4]
  quaternion = [float(0.0) for x in range(4)]

  # generate the upper triangle of the quadratic form matrix
  xxyx = float(0.0)
  xxyy = float(0.0)
  xxyz = float(0.0)
  xyyx = float(0.0)
  xyyy = float(0.0)
  xyyz = float(0.0)
  xzyx = float(0.0)
  xzyy = float(0.0)
  xzyz = float(0.0)

  for i in range(len(fit_xyz)):
    xxyx = xxyx + fit_xyz[i].x * ref_xyz[i].x * weights[i]
    xxyy = xxyy + fit_xyz[i].x * ref_xyz[i].y * weights[i]
    xxyz = xxyz + fit_xyz[i].x * ref_xyz[i].z * weights[i]
    xyyx = xyyx + fit_xyz[i].y * ref_xyz[i].x * weights[i]
    xyyy = xyyy + fit_xyz[i].y * ref_xyz[i].y * weights[i]
    xyyz = xyyz + fit_xyz[i].y * ref_xyz[i].z * weights[i]
    xzyx = xzyx + fit_xyz[i].z * ref_xyz[i].x * weights[i]
    xzyy = xzyy + fit_xyz[i].z * ref_xyz[i].y * weights[i]
    xzyz = xzyz + fit_xyz[i].z * ref_xyz[i].z * weights[i]

  matrix[0][0] = xxyx + xyyy + xzyz
  matrix[0][1] = xzyy - xyyz
  matrix[1][1] = xxyx - xyyy - xzyz
  matrix[0][2] = xxyz - xzyx
  matrix[1][2] = xxyy + xyyx
  matrix[2][2] = xyyy - xzyz - xxyx
  matrix[0][3] = xyyx - xxyy
  matrix[1][3] = xzyx + xxyz
  matrix[2][3] = xyyz + xzyy
  matrix[3][3] = xzyz - xxyx - xyyy

  # diagonalize c
  eigenvect, eigenval, maxsweeps = jacobi(matrix, maxsweeps)

  # extract the desired quaternion
  quaternion[0] = eigenvect[0][3]
  quaternion[1] = eigenvect[1][3]
  quaternion[2] = eigenvect[2][3]
  quaternion[3] = eigenvect[3][3]

  # generate the rotation matrix
  rotmat = q2mat(quaternion)

  return quaternion, rotmat, maxsweeps

### --- Run the classic quatfit --- ###
def quatfitClassic(reffile, fitfile, pairsfile, ofile, statfile):
  """
  This runs quatfit in it's classic form you must give it:
  reffile - reference coordinate file (xyz)
  fitfile - fit coordinate file (xyz) (The one you want to change)
  pairsfile - the pairs file
  ofile - the output file name
  statfile - the stat file name if wanted
  """

  #Check to make sure all the files are accounted for
  if ofile == '':
    print("You must specify an output file.")
    sys.exit(0)
  elif reffile == '':
    print("You must specify an reference file.")
    sys.exit(0)
  elif fitfile == '':
    print("You must specify an fit file.")
    sys.exit(0)

  #Create the xyz coord in Atom class
  reffilelol = parseXYZ(reffile)
  fitfilelol = parseXYZ(fitfile)

  #Parse the pairs file
  pairs, weights = parsePairs(pairsfile)

  #Create the reference coords based on the pairs
  ref_xyz, fit_xyz = createRefGeom(reffilelol, fitfilelol, pairs)

  #Center the reference coords around 0,0,0
  refcenter, ref_xyz = center(ref_xyz, weights, 1, [float(0), float(0), float(0)])
  fitcenter, fit_xyz = center(fit_xyz, weights, 1, [float(0), float(0), float(0)])

  #fit the specified atom coords of the fit to reference
  quaternion, rotmat, maxsweeps = qtrfit(fit_xyz, ref_xyz, weights, 30)

  #subtract coordinates of the center of fitted atoms of the fitted molecule
  #from all atom coordinates of the fitted molecule (note that weight is
  #a dummy parameter)
  fitcenter, fitfilelol = center(fitfilelol, weights, 2, fitcenter)

  # rotate the fitted molecule by the rotation matrix u
  fitfilelol = rotmol(fitfilelol, rotmat)

  # same with set of fitted atoms of the fitted molecule
  fit_xyz = rotmol(fit_xyz, rotmat)

  ### Haven't yet fully translated this section
  # if modes given in fitted molecule, rotate the modes too
  #if n_fields_f > 4:
  #  rotmol(nat_f, modes_f, modes_f, u)
  #  # calculate dot product of reference and fitted molecule modes
  #  if n_fields_r > 4:
  #    dotm = 0.0
  #    #for (i = 1; i <= npairs; i++) {
  #    for i in range(1,len(pairs)):
  #      #for (j = 1; j <= 3; j++) {
  #      for j in range(1,3):
  #        dotm += modes_r[j][atoms_r[i]]*modes_f[j][atoms_f[i]]

  # translate atoms of the fitted molecule to the center
  # of fitted atoms of the reference molecule
  refcenter, fitfilelol = center(fitfilelol, weights, 3, refcenter)

  # same with set of fitted atoms of the fitted molecule
  refcenter, fit_xyz = center(fit_xyz, weights, 3, refcenter)

  # translate fitted atoms of reference molecule to their orig. location
  refcenter, ref_xyz = center(ref_xyz, weights, 3, refcenter)

  # write modified XYZ file for fitted molecule
  outputXYZ(ofile, fitfilelol, 1)

  # find distances between fitted and reference atoms and print them in out file
  printlines = ["", "Distances and weighted distances between fitted atoms",
                "Ref.At. Fit.At.  Distance  Dist*sqrt(weight)  weight"]

  rms = 0.0
  wnorm = 0.0
  for i in range(len(pairs)):
    d = 0.0
    for j in range(3):
      if j == 0:
        s = ref_xyz[i].x - fit_xyz[i].x
      elif j == 1:
        s = ref_xyz[i].y - fit_xyz[i].y
      elif j == 2:
        s = ref_xyz[i].z - fit_xyz[i].z
      d += s*s
    rms += d
    printlines.append("  " + str("{0: <4}".format(ref_xyz[i].e)) + "    " + str("{0: <5}".format(fit_xyz[i].e)) + "  "
                      + str("{:.6f}".format(d)) + "      " + str("{:.6f}".format(weights[i]*d)) + "      "
                      + str("{:.6f}".format(weights[i])))
  rms = math.sqrt(rms/len(pairs))

  printlines.append("\nWeighted root mean square = " + str("{:.6f}".format(rms)))

  printlines.append("\nCenter of reference molecule fitted atoms")
  printlines.append("Xc = " + str("{:.6f}".format(refcenter[0])) + "  Yc = " + str("{:.6f}".format(refcenter[1])) +
                    "  Zc = " + str("{:.6f}".format(refcenter[2])))

  printlines.append("\nCenter of fitted molecule fitted atoms")
  printlines.append("Xc = " + str("{:.6f}".format(fitcenter[0])) + "  Yc = " + str("{:.6f}".format(fitcenter[1])) +
                    "  Zc = " + str("{:.6f}".format(fitcenter[2])))

  printlines.append("\nLeft rotation matrix")
  for i in range(3):
    printlines.append(" " + str("{:10.6f}".format(rotmat[i][0])) + "  " + str("{:10.6f}".format(rotmat[i][1])) +
                      "  " + str("{:10.6f}".format(rotmat[i][2])))

  #Write out the stats if there is a statfile else print to screen
  if statfile == '':
    for line in printlines:
      print(line)
  else:
    f = open(statfile, 'w+')
    for line in printlines:
      f.write(line + "\n")
    f.close()

  #/* {
  # double rmsd = 0.0;
  # double x,y,z;
  # for(int i = 1; i <= npairs; i++)
  # {
  #   x = ref_xyz[1][i] - fit_xyz[1][i];
  #   y = ref_xyz[2][i] - fit_xyz[2][i];
  #   z = ref_xyz[3][i] - fit_xyz[3][i];
  #
  #   rmsd += x*x;
  #   rmsd += y*y;
  #   rmsd += z*z;
  #   //rmsd += (pow(ref_xyz[1][i] - fit_xyz[1][i], 2) + pow(ref_xyz[2][i] - fit_xyz[2][i], 2) + pow(ref_xyz[3][i] - fit_xyz[3][i], 2));
  # }
  # //rmsd /= static_cast<double>(npairs);
  # rms = sqrt(rmsd/npairs);
  # }*/
  #if((n_fields_f > 4) && (n_fields_r > 4)) {
  #  fprintf(statfile,
  #  "\nDot product of normal modes on fitted atom pairs =%11.6f\n", dotm);
  #  }

  return

### --- Run the classic quatfit --- ###
def quatfitGetMolecule(reffilelol, fitfilelol, pairs, weights):
  """
  This runs quatfit in a updated form
  reffilelol - reference coordinate in Atom format
  fitfilelol - fit coordinate in molecule format
  pairs - the pairs in a list of lists (i.e. [[2, 3], [3, 13], ...])
  weights - the weights in a list (i.e. [1, 1, ...])
  """

  #Create the reference coords based on the pairs
  ref_xyz, fit_xyz = createRefGeom(reffilelol, fitfilelol, pairs)

  #Center the reference coords around 0,0,0
  refcenter, ref_xyz = center(ref_xyz, weights, 1, [float(0), float(0), float(0)])
  fitcenter, fit_xyz = center(fit_xyz, weights, 1, [float(0), float(0), float(0)])

  #fit the specified atom coords of the fit to reference
  quaternion, rotmat, maxsweeps = qtrfit(fit_xyz, ref_xyz, weights, 30)

  #subtract coordinates of the center of fitted atoms of the fitted molecule
  #from all atom coordinates of the fitted molecule (note that weight is
  #a dummy parameter)
  fitcenter, fitfilelol = center(fitfilelol, weights, 2, fitcenter)

  # rotate the fitted molecule by the rotation matrix u
  fitfilelol = rotmol(fitfilelol, rotmat)

  # same with set of fitted atoms of the fitted molecule
  fit_xyz = rotmol(fit_xyz, rotmat)

  # translate atoms of the fitted molecule to the center
  # of fitted atoms of the reference molecule
  refcenter, fitfilelol = center(fitfilelol, weights, 3, refcenter)

  # same with set of fitted atoms of the fitted molecule
  refcenter, fit_xyz = center(fit_xyz, weights, 3, refcenter)

  # translate fitted atoms of reference molecule to their orig. location
  refcenter, ref_xyz = center(ref_xyz, weights, 3, refcenter)

  rms = 0.0
  wnorm = 0.0
  for i in range(len(pairs)):
    d = 0.0
    for j in range(3):
      if j == 0:
        s = ref_xyz[i].x - fit_xyz[i].x
      elif j == 1:
        s = ref_xyz[i].y - fit_xyz[i].y
      elif j == 2:
        s = ref_xyz[i].z - fit_xyz[i].z
      d += s*s
    rms += d

  rms = math.sqrt(rms/len(pairs))

  return fitfilelol, rms
