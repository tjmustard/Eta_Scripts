#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright (c) 2025, Thomas J. L. Mustard, O. Maduka Ogba, Paul Ha-Yeon Cheong
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
import etamap
import etastructure
import quatfit

import argparse
import os


def parse_args():
    """
    Parse command-line arguments using ArgParse.

    Returns:
        Namespace: An object containing parsed argument values.
    """

    # Create a parser object
    parser = argparse.ArgumentParser(description='Example ArgParse Script')
    parser.add_argument('-p', '--parent_file', required=True, help='Path to parent file')
    parser.add_argument('--parent_align', nargs='+', type=int, help='List of atoms to align with from parent structure.')
    parser.add_argument('--parent_delete', nargs='+', type=int, help='List of atoms to delete from parent structure.')
    parser.add_argument('-c', '--child_file', required=True, help='Path to child file')
    parser.add_argument('--child_align', nargs='+', type=int, help='List of atoms to align with from child structure.')
    parser.add_argument('--child_delete', nargs='+', type=int, help='List of atoms to delete from child structure.')
    parser.add_argument('-m', '--map_file', required=False, help='Path to map file')

    # Add boolean flag for verbose output
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help='Print extra information')

    return parser.parse_args()


def read_input_file(file_name):
    """
    Read an input file with key-value pairs and store them in variables.

    Parameters:
        file_name (str): The name of the input file

    Returns:
        A dictionary containing the parsed key-value pairs
    """

    # Initialize an empty dictionary to store the data
    data = {}

    try:
        # Open the input file for reading
        with open(file_name, 'r') as f:

            # Iterate over each line in the file
            for line in f.readlines():
                # Remove leading/trailing whitespace and split into key-value pair
                line = line.strip()
                if '=' in line:
                    key, value = [x.strip() for x in line.split('=')]

                    # Store the key-value pair in the dictionary
                    data[key] = value

    except FileNotFoundError:
        print(f"Error: Input file '{file_name}' not found.")

    return data

def pull_xyz_coords(struct, align, all=False):
    coords = []
    if all:
        for atom in struct.atoms:
            #print("{} {} {}".format(atom['x'], atom['y'], atom['z']))
            coords.append([atom['x'], atom['y'], atom['z']])
    else:
        for atom in align:
            print("{} {} {} {}".format(struct.atoms[atom-1]['element'], struct.atoms[atom-1]['x'], struct.atoms[atom-1]['y'], struct.atoms[atom-1]['z']))
            coords.append([struct.atoms[atom-1]['x'], struct.atoms[atom-1]['y'], struct.atoms[atom-1]['z']])
    return coords

def main():
    """
    The main function where everything gets executed.

    Returns:
        None: This is the typical signature of a Python `main` method that doesn't explicitly indicate any value being returned. If needed in your code, consider returning something (if applicable).
    """

    # Get parsed arguments
    args = parse_args()

    if len(args.parent_align) != len(args.child_align):
        print("Error: Number of Parent and Child alignment atoms do not match.")
        exit(1)

    # Get the path to the input file
    pst = etastructure.ChemicalStructure.load_xyz(args.parent_file)
    cst = etastructure.ChemicalStructure.load_xyz(args.child_file)

    print("Child align atom indexes:")
    print(args.child_align)
    print("Parent align atom indexes:")
    print(args.parent_align)

    pcoords = pull_xyz_coords(pst, args.parent_align)

    print("Align Parent coords:")
    print("____________________________________")
    print(pcoords)

    ccoords = pull_xyz_coords(cst, args.child_align)

    print("Align Child coords:")
    print("____________________________________")
    print(ccoords)

    accoords = pull_xyz_coords(cst, args.child_align, all=True)

    print("ALL Child coords:")
    print("____________________________________")
    print(len(accoords))
    print(accoords)

    aligned_child, rms = quatfit.fit_fragment(accoords, ccoords, pcoords)

    #aligned_child = quatfit.find_coordinates(len(pcoords), pcoords, ccoords, accoords)

    #print("Aligned coords:")
    #print("____________________________________")
    #print(len(aligned_child))
    #print(aligned_child)
    print("RMS:")
    print("____________________________________")
    print(rms)

    for atomindex, coord in enumerate(aligned_child):
        cst.atoms[atomindex]['x'] = coord[0]
        cst.atoms[atomindex]['y'] = coord[1]
        cst.atoms[atomindex]['z'] = coord[2]

    etastructure.ChemicalStructure.export_xyz(cst, "aligned.xyz")

    # Check if verbose output is enabled
    if args.verbose:
        print("Verbose mode: ON")

    if args.map_file:
        data = read_input_file(args.map_file)
        print(data)


if __name__ == '__main__':
    main()
