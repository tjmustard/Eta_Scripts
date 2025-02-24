# coding=utf-8

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
# This program was heavily modified by Daniel Kratzert
#
# Then *very* slightly modified again by Thomas J L Mustard
# Thank you for making this useful again Daniel!

import copy
from math import fabs, sqrt

def matrix_minus_vect(m, v):
    """
    """
    result = []
    for coords in m:
      result.append([coords[0] - v[0], coords[1] - v[1], coords[2] - v[2]])
    return result

def matrix_plus_vect(m, v):
  """
  """
  result = []
  for coords in m:
    result.append([coords[0] + v[0], coords[1] + v[1], coords[2] + v[2]])
  return result

def transpose(a):
    """
    transposes a matrix
    """
    return list(zip(*a))

def rotmol(frag_atoms, rotmat):
    """
    ROTMOL
    rotate a molecule
    y = u * x

    frag_atoms (x) - input coordinates
    rotmat (u) - left rotation matrix
    output (y) - rotated coordinates
    """
    atoms = copy.deepcopy(frag_atoms)
    for i in range(len(atoms)):
        yx = rotmat[0][0] * atoms[i][0] + rotmat[1][0] * atoms[i][1] + rotmat[2][0] * atoms[i][2]
        yy = rotmat[0][1] * atoms[i][0] + rotmat[1][1] * atoms[i][1] + rotmat[2][1] * atoms[i][2]
        yz = rotmat[0][2] * atoms[i][0] + rotmat[1][2] * atoms[i][1] + rotmat[2][2] * atoms[i][2]
        atoms[i][0] = yx  # x
        atoms[i][1] = yy  # y
        atoms[i][2] = yz  # z
    return atoms


def jacobi(matrix, maxsweeps):
    """
    JACOBI
    Jacobi diagonalizer with sorted output. It is only good for 4x4 matrices.
    (was too lazy to do pointers...)

    INPUTS
    matrix (a) - input: matrix to diagonalize
    maxsweeps (nrot) - input: maximum number of sweeps

    OUTPUTS
    eigenvect (v) - output: eigenvectors
    eigenval (d) - output: eigenvalues
    """
    eigenvect = [[float(0.0) for _ in range(4)] for _ in range(4)]
    eigenval = [float(0.0) for _ in range(4)]
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
    m = 0
    for j in range(4):
        eigenvect[j][j] = 1.0
        eigenval[j] = matrix[j][j]
    for m in range(maxsweeps):
        dnorm = 0.0
        onorm = 0.0
        for j in range(4):
            dnorm += fabs(eigenval[j])
            for i in range(j):
                onorm += fabs(matrix[i][j])
        if onorm / dnorm <= 1.0e-12:
            # Solution converged
            break
        for j in range(1, 4):
            for i in range(j):
                b = matrix[i][j]
                if fabs(b) > 0.0:
                    dma = eigenval[j] - eigenval[i]
                    if (fabs(dma) + fabs(b)) <= fabs(dma):
                        t = b / dma
                    else:
                        q = 0.5 * dma / b
                        t = 1.0 / (fabs(q) + sqrt(1.0 + q * q))
                        if q < 0.0:
                            t = -t
                    c = 1.0 / sqrt(t * t + 1.0)
                    s = t * c
                    matrix[i][j] = 0.0
                    for k in range(i):
                        atemp = c * matrix[k][i] - s * matrix[k][j]
                        matrix[k][j] = s * matrix[k][i] + c * matrix[k][j]
                        matrix[k][i] = atemp
                    for k in range(i + 1, j):
                        atemp = c * matrix[i][k] - s * matrix[k][j]
                        matrix[k][j] = s * matrix[i][k] + c * matrix[k][j]
                        matrix[i][k] = atemp
                    for k in range(j + 1, 4):
                        atemp = c * matrix[i][k] - s * matrix[j][k]
                        matrix[j][k] = s * matrix[i][k] + c * matrix[j][k]
                        matrix[i][k] = atemp
                    for k in range(4):
                        vtemp = c * eigenvect[k][i] - s * eigenvect[k][j]
                        eigenvect[k][j] = s * eigenvect[k][i] + c * eigenvect[k][j]
                        eigenvect[k][i] = vtemp
                    dtemp = c * c * eigenval[i] + s * s * eigenval[j] - 2.0 * c * s * b
                    eigenval[j] = s * s * eigenval[i] + c * c * eigenval[j] + 2.0 * c * s * b
                    eigenval[i] = dtemp
    # print(m)
    for j in range(3):
        k = j
        dtemp = eigenval[k]
        for i in range(j + 1, 4):
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
    return eigenvect, eigenval


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
    rotmat[0][0] = quaternion[0] * quaternion[0] + quaternion[1] * quaternion[1] - quaternion[2] * quaternion[2] - \
                   quaternion[3] * quaternion[3]
    rotmat[1][0] = 2.0 * (quaternion[1] * quaternion[2] - quaternion[0] * quaternion[3])
    rotmat[2][0] = 2.0 * (quaternion[1] * quaternion[3] + quaternion[0] * quaternion[2])
    rotmat[0][1] = 2.0 * (quaternion[2] * quaternion[1] + quaternion[0] * quaternion[3])
    rotmat[1][1] = quaternion[0] * quaternion[0] - quaternion[1] * quaternion[1] + quaternion[2] * quaternion[2] - \
                   quaternion[3] * quaternion[3]
    rotmat[2][1] = 2.0 * (quaternion[2] * quaternion[3] - quaternion[0] * quaternion[1])
    rotmat[0][2] = 2.0 * (quaternion[3] * quaternion[1] - quaternion[0] * quaternion[2])
    rotmat[1][2] = 2.0 * (quaternion[3] * quaternion[2] + quaternion[0] * quaternion[1])
    rotmat[2][2] = quaternion[0] * quaternion[0] - quaternion[1] * quaternion[1] - quaternion[2] * quaternion[2] + \
                   quaternion[3] * quaternion[3]
    return rotmat



def qtrfit(parent_coords, child_coords, maxsweeps=30):
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
       child_coords (x)      - fitted molecule coordinates
       parent_coords (y)      - reference molecule coordinates
       maxsweeps (nr)     - max number of jacobi sweeps

       #n      - number of points
       #weights (w)      - weights

     OUTPUT
       quaternion (q)     - the best-fit quaternion
       rotmat (u)         - the best-fit left rotation matrix

    """

    # Create variables/lists/matrixes
    matrix = [[float(0.0) for x in range(4)] for x in range(4)]
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
    for i in range(len(parent_coords)):
        xxyx += float(parent_coords[i][0]) * float(child_coords[i][0])
        xxyy += float(parent_coords[i][0]) * float(child_coords[i][1])
        xxyz += float(parent_coords[i][0]) * float(child_coords[i][2])
        xyyx += float(parent_coords[i][1]) * float(child_coords[i][0])
        xyyy += float(parent_coords[i][1]) * float(child_coords[i][1])
        xyyz += float(parent_coords[i][1]) * float(child_coords[i][2])
        xzyx += float(parent_coords[i][2]) * float(child_coords[i][0])
        xzyy += float(parent_coords[i][2]) * float(child_coords[i][1])
        xzyz += float(parent_coords[i][2]) * float(child_coords[i][2])

    matrix[0][0] = xxyx + xyyy + xzyz
    matrix[0][1] = xzyy - xyyz
    matrix[0][2] = xxyz - xzyx
    matrix[0][3] = xyyx - xxyy
    matrix[1][1] = xxyx - xyyy - xzyz
    matrix[1][2] = xxyy + xyyx
    matrix[1][3] = xzyx + xxyz
    matrix[2][2] = xyyy - xzyz - xxyx
    matrix[2][3] = xyyz + xzyy
    matrix[3][3] = xzyz - xxyx - xyyy

    # diagonalize c
    eigenvect, eigenval = jacobi(matrix, maxsweeps)

    #print(eigenvect)
    #print(eigenval)

    # extract the desired quaternion
    quaternion[0] = eigenvect[0][3]
    quaternion[1] = eigenvect[1][3]
    quaternion[2] = eigenvect[2][3]
    quaternion[3] = eigenvect[3][3]

    # generate the rotation matrix
    rotmat = q2mat(quaternion)

    return quaternion, rotmat


def centroid(X):
    """
    Calculate the centroid from a vectorset X.

    https://en.wikipedia.org/wiki/Centroid
    Centroid is the mean position of all the points in all of the coordinate
    directions.

    C = sum(X)/len(X)

    Parameters
    ----------
    X : np.array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    C : float
        centeroid
    """
    s = [0.0, 0.0, 0.0]
    num = 0
    for line in X:
        num += 1
        for n, v in enumerate(line):
            s[n] += v
    return (s[0] / num, s[1] / num, s[2] / num)


def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.

    Parameters
    ----------
    V : np.array
        (N,D) matrix, where N is points and D is dimension.
    W : np.array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    fit : float
        Root-mean-square deviation

    """
    d = 0
    for v, w in zip(V, W):
        d += sum([(h - k) ** 2 for h, k in zip(v, w)])
    return sqrt(d / len(V))

def fit_fragment(child_atoms, child_align_atoms, parent_align_atoms):
    """
    Takes a list of fragment atoms (child_atoms) and fits them to the position of target atoms.
    The child_allign_atoms are a fraction of the child_atoms fragment to be fitted on the parent_align_atoms.

    Parameters
    ----------
    child_atoms: list
        complete set of atoms of a fragment
    child_allign_atoms: list
        subsection of fragment atoms
    parent_align_atoms: list
        target position for child_allign_atoms

    Returns
    -------
    rotated_fragment: list
        list of coordinates from the fitted fragment
    rmsd: float
        RMSD (root mean square deviation)
    """
    #Find the centroids of the two alignment coordinate lists
    Ccentroid = centroid(child_align_atoms)
    Pcentroid = centroid(parent_align_atoms)

    # Move child_allign_atoms and parent_align_atoms to origin:
    centered_child_allign_atoms = matrix_minus_vect(child_align_atoms, Ccentroid)
    centered_parent_align_atoms = matrix_minus_vect(parent_align_atoms, Pcentroid)

    # get the Kabsch rotation matrix:
    quaternion, U = qtrfit(centered_parent_align_atoms, centered_child_allign_atoms, maxsweeps=30)

    #Move the child_atoms to the center before rotation
    centered_child_atoms = matrix_minus_vect(child_atoms, Ccentroid)

    # rotate all child_atoms (instead of child_allign_atoms):
    aligned_child = rotmol(centered_child_atoms, U)

    # rotate also source atoms for rmsd calculation:
    rotated_child = rotmol(centered_child_allign_atoms, U)

    # move fragment back from zero (be aware that the translation is still wrong!):
    aligned_child = matrix_plus_vect(aligned_child, Pcentroid)

    rms = rmsd(centered_parent_align_atoms, rotated_child)
    return aligned_child, rms

