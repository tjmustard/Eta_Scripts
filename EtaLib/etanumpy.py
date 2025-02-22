#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#   Copyright (c) 2014, Thomas J. L. Mustard, Joshua J. Kincaid, O. Maduka Ogba, Paul Ha-Yeon Cheong
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

### --- Special Acknowledgments --- ###
# Joshua J. Kincaid made the z-matrix to xyz conversion function possible.
# Without his help several scripts/programs would be useless.
# Thank you for your assistance Joshua.

import os
from sys import *
import math
from . import etaatom
from numpy import *
from decimal import *


### --- Functions to get geometric data --- ###
###==================================================================================================================###

### --- Distance function --- ###

def get_distance(at1, at2, lol):
    if isinstance(lol[at1], etaatom.Atom):
        atom1 = array([lol[at1].x, lol[at1].y, lol[at1].z])
        atom2 = array([lol[at2].x, lol[at2].y, lol[at2].z])
    else:
        atom1 = array([float(lol[at1][1]), float(lol[at1][2]),
                      float(lol[at1][3])])
        atom2 = array([float(lol[at2][1]), float(lol[at2][2]),
                      float(lol[at2][3])])
    vector1 = atom2 - atom1
    dist = linalg.norm(vector1)

  # print dist

    return dist


### --- Angle function --- ###

def get_angle(
    at1,
    at2,
    at3,
    lol,
    ):

    if isinstance(lol[at1], etaatom.Atom):

    # put positions in array

        atom1 = array([lol[at1].x, lol[at1].y, lol[at1].z])
        atom2 = array([lol[at2].x, lol[at2].y, lol[at2].z])
        atom3 = array([lol[at3].x, lol[at3].y, lol[at3].z])
    else:
        atom1 = array([float(lol[at1][1]), float(lol[at1][2]),
                      float(lol[at1][3])])
        atom2 = array([float(lol[at2][1]), float(lol[at2][2]),
                      float(lol[at2][3])])
        atom3 = array([float(lol[at3][1]), float(lol[at3][2]),
                      float(lol[at3][3])])

  # making appropriate vectors and normals

    vector1 = atom1 - atom2
    vector2 = atom3 - atom2
    angle = arccos(dot(vector1, vector2) / (linalg.norm(vector1)
                   * linalg.norm(vector2)))

  # print degrees(angle)

    return degrees(angle)


### --- Dihedral angle function --- ###

def get_dihedral(
    at1,
    at2,
    at3,
    at4,
    lol,
    ):

    if isinstance(lol[at1], etaatom.Atom):

    # put positions in array

        atom1 = array([lol[at1].x, lol[at1].y, lol[at1].z])
        atom2 = array([lol[at2].x, lol[at2].y, lol[at2].z])
        atom3 = array([lol[at3].x, lol[at3].y, lol[at3].z])
        atom4 = array([lol[at4].x, lol[at4].y, lol[at4].z])
    else:
        atom1 = array([float(lol[at1][1]), float(lol[at1][2]),
                      float(lol[at1][3])])
        atom2 = array([float(lol[at2][1]), float(lol[at2][2]),
                      float(lol[at2][3])])
        atom3 = array([float(lol[at3][1]), float(lol[at3][2]),
                      float(lol[at3][3])])
        atom4 = array([float(lol[at4][1]), float(lol[at4][2]),
                      float(lol[at4][3])])

  # making appropriate vectors and normals

    vector1 = atom2 - atom1
    vector2 = atom3 - atom2
    plane1 = cross(vector1, vector2)
    vector3 = atom2 - atom3
    vector4 = atom4 - atom3
    plane2 = cross(vector3, vector4)

  # finding dihedral angle

    dihedral = arccos(-dot(plane1, plane2) / (linalg.norm(plane1)
                      * linalg.norm(plane2)))

  # checking the sign of the dihedral then displaying result

    if dot(plane1, vector4) > 0:

    # print degrees(dihedral)

        return degrees(dihedral)
    else:

    # print -degrees(dihedral)

        return -degrees(dihedral)


### WHY DO I HAVE TWO OF THESE!?
#### --- TWO STRUCTURE Distance function --- ###
# def two_structure_get_distance(at1, at2, lol1, lol2):
#  #if isinstance(lol1[at1], etaatom.Atom):
#  atom1 = array([lol1[at1].x, lol1[at1].y, lol1[at1].z])
#  atom2 = array([lol2[at2].x, lol2[at2].y, lol2[at2].z])
#  #else:
#  #  atom1 = array([float(lol1[at1][1]), float(lol1[at1][2]), float(lol1[at1][3])])
#  #  atom2 = array([float(lol2[at2][1]), float(lol2[at2][2]), float(lol2[at2][3])])
#  vector1 = atom2-atom1
#  dist = linalg.norm(vector1)
#  #print dist
#  return dist

### --- Distance function TWOSTRUCTURE --- ###

def two_structure_get_distance(
    at1,
    at2,
    lol1,
    lol2,
    ):

    if isinstance(lol1[at1], etaatom.Atom):
        atom1 = array([lol1[at1].x, lol1[at1].y, lol1[at1].z])
        atom2 = array([lol2[at2].x, lol2[at2].y, lol2[at2].z])
    else:
        atom1 = array([float(lol1[at1][1]), float(lol1[at1][2]),
                      float(lol1[at1][3])])
        atom2 = array([float(lol2[at2][1]), float(lol2[at2][2]),
                      float(lol2[at2][3])])
    vector1 = atom2 - atom1
    dist = linalg.norm(vector1)

  # print dist

    return dist


###==================================================================================================================###

### --- Bond elucidation section --- ###
###==================================================================================================================###
### --- Get bonds --- ###

def get_bonds(fileatom):
    ifileatomnum = fileatom[0]
    ifilename = fileatom[1]
    covbonds = []
    covhbonds = []
    covtmbonds = []
    nearestneighbor = []
    neighborstart = [0, 1000000]

  # Iterate through the file

    for (i, linei) in enumerate(fileatom):

    # Set a placeholder for the nearest neighbor

        nearestneighbor.append(list(neighborstart))

    # Iterate throught the file again

        for (j, linej) in enumerate(fileatom):
            if i != j:
                distij = getdistance(i + 2, j + 2, fileatom)
                if j > i:
                    if distij <= 2.25 and fileatom[i + 2][0] != 1 \
                        and fileatom[j + 2][0] != 1 \
                        and etaatom.get_element_name(fileatom[i
                            + 2][0]).lower() \
                        not in etaatom.elementLarge \
                        and etaatom.get_element_name(fileatom[j
                            + 2][0]).lower() \
                        not in etaatom.elementLarge:
                        distlist = [i + 1, j + 1, distij]
                        covbonds.append(distlist)
                    elif distij <= 1.3 and (fileatom[i + 2][0] == 1
                            or fileatom[j + 2][0] == 1):

            # print str(i + 2) + "\t" + str(j + 2) + "\t" + str(distij)

                        distlist = [i + 1, j + 1, distij]
                        covhbonds.append(distlist)
                    elif distij <= 3 \
                        and etaatom.get_element_name(fileatom[i
                            + 2][0]).lower() in etaatom.elementLarge \
                        and etaatom.get_element_name(fileatom[j
                            + 2][0]).lower() in etaatom.elementLarge:

                        distlist = [i + 1, j + 1, distij]
                        covtmbonds.append(distlist)
                if distij < nearestneighbor[i][1]:
                    nearestneighbor[i][0] = j + 1
                    nearestneighbor[i][1] = distij

  # ## --- Remove hydrogen bonds from bond list --- ###

    for i in range(0, len(covhbonds)):
        if covhbonds[i][0] != nearestneighbor[covhbonds[i][0]][0] \
            and covhbonds[i][1] != nearestneighbor[covhbonds[i][0]][0]:
            del covhbonds[i]

  # print "Covalent bonds to Hydrogen:"
  # print covhbonds
  # print "Covalent bonds between \"heavy\" atoms."
  # print covbonds
  # print "Covalent bonds between TM atoms."
  # print covtmbonds
  # print nearestneighbor

    return fileatom


### --- Build a list that contains all atoms starting with at1
#  on one side of the molecule that is opposite the other atom atm2

def build_bond_list(at1, at2, filelol):
    bondedlist = [1, 2, 3, 4]
    return bondedlist


###==================================================================================================================###

### --- ZMAT to XYZ to ZMAT conversion section --- ###
###==================================================================================================================###

### --- get the distance, angle and dihedrals for a Z-matrix --- ###
# This function is iteratable through a file, allowing for:
# 1. grabbing of a single atoms's distance, angle and dihedral
# 2. double ended zmatrix files
# 3. complex zmatrix files dependent on bonding
# You must give the (up to) 4 atoms nessecary to get the zmat for the atom

def get_zmat(
    i,
    j,
    k,
    l,
    ifilelol,
    ):

  # If the line i is in Atom format find the relevant distance, angle, and dihedral

    if isinstance(ifilelol[i], etaatom.Atom):
        line = []
        line.append(i - 1)
        line.append(ifilelol[i].e)

    # If the atom is not the first in the file get the distance to the last atom

        if i > 2:
            line.append(j - 1)
            dist = getdistance(j, i, ifilelol)
            line.append(dist)

    # If the atom is not the first or second in the file get the angle to the last two atoms

        if i > 3:
            line.append(k - 1)
            angle = getangle(k, j, i, ifilelol)
            line.append(angle)

    # If the atom is not the first, second, or third in the file get the dihedral to the last three atoms

        if i > 4:
            line.append(l - 1)
            dihedral = getdihedral(l, k, j, i, ifilelol)
            line.append(dihedral)
    else:

  # If not an Atom (ie first two lines) print them back out

        line = [ifilelol[i]]
    line.append(-1)
    line.append(-1)
    line.append(-1)
    line.append(-1)
    line.append(-1)
    line.append(-1)
    line.append('\n')
    return line


### --- Get the XYZ coordinates from distance, angle and dihedral data --- ###

def get_xyz_from_zmat(zmatlol):
    xyzlol = ['0'] * len(zmatlol)
    xyzlol[0] = zmatlol[0][0]
    xyzlol[1] = zmatlol[1][0]
    for i in range(2, len(zmatlol)):

    # ## Set up the variables to be used for the function

        dist = float(zmatlol[i][3])
        angle = float(zmatlol[i][5])
        dihedral = float(zmatlol[i][7])
        anglerad = radians(angle)  # * math.pi / 180
        dihedralrad = radians(dihedral)  # * math.pi / 180
        at1 = int(zmatlol[i][6]) + 1
        at2 = int(zmatlol[i][4]) + 1
        at3 = int(zmatlol[i][2]) + 1
        at4 = int(zmatlol[i][0]) + 1
        x = 0
        y = 0
        z = 0

    # ## Start to place the atoms in their locations

        if at4 == 2:
            x = 0.00000
            y = 0.00000
            z = 0.00000
        elif at4 == 3:

            x = dist
            y = 0.00000
            z = 0.00000
        elif at4 == 4:

            a = xyzlol[at3].x
            b = dist
            x = a + dist * cos(math.pi - anglerad)
            y = 0.00000
            z = -dist * sin(math.pi - anglerad)
        elif at4 >= 5:

      # ###The at4 x,y,z coordinates from spherical coord

            sx = dist * sin(anglerad) * cos(dihedralrad)
            sy = -dist * sin(anglerad) * sin(dihedralrad)  # For some reason I need to have a negative here to get the correct sign in the output..... weird
            sz = dist * cos(anglerad)
            at4l = [sx, sy, sz]

      # ##Finding the angle theta
      # Make the list of lists for the three point (z-axis, origin, and translated atom 2) needed for an angle calculation

            z32 = [[0, 0, 0, 1], [0, 0, 0, 0], [0, xyzlol[at2].x
                   - xyzlol[at3].x, xyzlol[at2].y - xyzlol[at3].y,
                   xyzlol[at2].z - xyzlol[at3].z]]

      # Get theta using the getangle function

            theta = radians(getangle(0, 1, 2, z32))

      # print "theta: " +str(theta)

      # ##Rodrigues' rotation formula
      # Create the vectors needed to calculate k

            vector3 = array([xyzlol[at3].x, xyzlol[at3].y,
                            xyzlol[at3].z])
            vector2 = array([xyzlol[at2].x, xyzlol[at2].y,
                            xyzlol[at2].z])
            vector0 = array([0, 0, 1])

      # Calculate k for the Rodrigues rotation formula

            k = cross(vector2 - vector3, vector0) \
                / linalg.norm(cross(vector2 - vector3, vector0))

      # Generate an array for translated 1

            t1 = [xyzlol[at1].x - xyzlol[at3].x, xyzlol[at1].y
                  - xyzlol[at3].y, xyzlol[at1].z - xyzlol[at3].z]

      # Calculate the Rodrigues rotation matrix

            rr23t1 = dot(t1, cos(theta)) + dot(cross(k, t1),
                    sin(theta)) + dot(dot(k, dot(k, t1)), 1
                    - cos(theta))

      # Make the list of lists for the four points (x-axis, z-axis, origin, and rotated translated 1) needed for a dihedral calculation

            xz31 = [[0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0], [0,
                    rr23t1[0], rr23t1[1], rr23t1[2]]]

      # Get phi using the getdihedral function

            phi = radians(getdihedral(0, 1, 2, 3, xz31))

      # ##Rotation matrix
      # Create the array for the rotation matrix including dihedral phi

            rm = array([[cos(phi), sin(phi), 0], [-sin(phi), cos(phi),
                       0], [0, 0, 1]])

      # Calculate the dot product of the rotation matrix and the coordinates for 4 (from spherical)

            rm4 = dot(rm, at4l)

      # Calculate the rotated coordinates of the rotated coordinates of atom 4

            rrn23rm4 = dot(rm4, cos(-theta)) + dot(cross(k, rm4),
                    sin(-theta)) + dot(dot(k, dot(k, rm4)), 1
                    - cos(-theta))

      # Final coordinates that are rotated, rotated and translated

            x = rrn23rm4[0] + xyzlol[at3].x
            y = rrn23rm4[1] + xyzlol[at3].y
            z = rrn23rm4[2] + xyzlol[at3].z

    # Putting everything into a Atom list to send back

        xyzlol[i] = etaatom.Atom()
        xyzlol[i].e = zmatlol[i][1]
        xyzlol[i].en = etaatom.get_element_num(zmatlol[i][1])
        xyzlol[i].x = x
        xyzlol[i].y = y
        xyzlol[i].z = z
    return xyzlol


### --- Generate a list of lists for the zmatrix of the structure --- ###

def gen_zmat_basic(ifilelol):
    zmatlol = []
    for i in range(0, len(ifilelol)):
        linetemp = getzmat(i, i - 1, i - 2, i - 3, ifilelol)
        zmatlol.append(linetemp)
    return zmatlol


def gen_zmat_double_sided(ifilelol, list1, list2):
    zmatlol = []
    for i in range(0, len(ifilelol)):
        linetemp = getzmat(i, i - 1, i - 2, i - 3, ifilelol)
        zmatlol.append(linetemp)
    return zmatlol


###==================================================================================================================###

# Compute the f and qt for the vibrational lists. Used for KIE

def compute_f_qt(vib, hkt, scaling):
    u = [0.00000] * len(vib)
    e = [0.00000] * len(vib)
    n = [0.00000] * len(vib)
    vibprod = [0.00000] * len(vib)
    for i in range(len(vib)):
        u[i] = vib[i] * hkt * scaling
    for i in range(len(vib)):
        if vib[i] > 0:
            e[i] = math.exp(u[i] / 2)
    for i in range(len(vib)):
        if vib[i] > 0:
            n[i] = 1 - math.exp(-u[i])
    for i in range(len(vib)):
        if vib[i] > 0:
            vibprod[i] = e[i] * n[i] / u[i]
    if vib[0] < 0:
        qt = u[0] / 2 / math.sin(u[0] / 2)
        del vibprod[0]
        f = prod(vibprod) / u[0]
    else:
        qt = u[0] / 2 / math.sin(u[0] / 2)
        f = prod(vibprod)
    return (f, qt)  # u, e, n, prod,


# Compute the KIE

def compute_kie(
    tsnatvib,
    tsisovib,
    subnatvib,
    subisovib,
    temp,
    scaling,
    ):

    # Calculate hkt

    hkt = 6.626E-34 / (1.381E-23 * temp) * 30000000000

    # Calculate the u and qt for each vibration list

    (subnatu, subnatqt) = compute_f_qt(subnatvib, hkt, scaling)
    (subisou, subisoqt) = compute_f_qt(subisovib, hkt, scaling)
    (tsnatu, tsnatqt) = compute_f_qt(tsnatvib, hkt, scaling)
    (tsisou, tsisoqt) = compute_f_qt(tsisovib, hkt, scaling)

    # Calculate the kie and quantum kie

    kie = tsisou * subnatu / (tsnatu * subisou)
    qkie = tsisou * subnatu / (tsnatu * subisou) * tsnatqt / tsisoqt
    return (kie, qkie)


