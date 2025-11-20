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

import argparse
from sys import *

from EtaLib import etastructure, quatfit


def parse_args():
    """
    Parse command-line arguments using ArgParse.

    Returns:
        Namespace: An object containing parsed argument values.
    """

    # Create a parser object
    parser = argparse.ArgumentParser(description="Example ArgParse Script")
    parser.add_argument(
        "-p", "--parent_file", required=True, help="Path to parent file"
    )
    parser.add_argument(
        "--parent_align",
        nargs="+",
        type=int,
        help="List of atoms to align with from parent structure.",
    )
    parser.add_argument(
        "--parent_delete",
        nargs="+",
        type=str,
        help="List of atoms to delete from parent structure.",
    )
    parser.add_argument("-c", "--child_file", required=True, help="Path to child file")
    parser.add_argument(
        "--child_align",
        nargs="+",
        type=int,
        help="List of atoms to align with from child structure.",
    )
    parser.add_argument(
        "--child_delete",
        nargs="+",
        type=str,
        help="List of atoms to delete from child structure.",
    )
    parser.add_argument("-m", "--map_file", required=False, help="Path to map file")

    # Add boolean flag for verbose output
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Print extra information",
    )

    return parser.parse_args()


def read_input_file(file_name):
    """
    Read an input file with key-value pairs and store them in variables.

    Args:
        file_name (str): The name of the input file

    Returns:
        A dictionary containing the parsed key-value pairs
    """

    # Initialize an empty dictionary to store the data
    data = {}

    try:
        # Open the input file for reading
        with open(file_name, "r") as f:

            # Iterate over each line in the file
            for line in f.readlines():
                # Remove leading/trailing whitespace and split into key-value pair
                line = line.strip()
                if "=" in line:
                    key, value = [x.strip() for x in line.split("=")]

                    # Store the key-value pair in the dictionary
                    data[key] = value

    except FileNotFoundError:
        print(f"Error: Input file '{file_name}' not found.")

    return data


def pull_xyz_coords(struct, align, all=False):
    """
    Extract just the X, Y, and Z coordinates from the ChemicalStructure class and return as a list of tuples.
    If an alignment list is provided return only those atoms in the order of the atom indexes provided
    If all is True return the entire list of atom coordinates.

    Args:
        struct: ChemicalStructure class
        allign: List of atom indexes to select
        all: Boolean to return all atoms or only a list of atoms.

    Returns:
        A list of coordinate tuples
    """

    coords = []
    if all:
        for atom in struct.atoms:
            # print("{} {} {}".format(atom['x'], atom['y'], atom['z']))
            coords.append([atom["x"], atom["y"], atom["z"]])
    else:
        for atom in align:
            # print("{} {} {} {}".format(struct.atoms[atom-1]['element'], struct.atoms[atom-1]['x'], struct.atoms[atom-1]['y'], struct.atoms[atom-1]['z']))
            coords.append(
                [
                    struct.atoms[atom - 1]["x"],
                    struct.atoms[atom - 1]["y"],
                    struct.atoms[atom - 1]["z"],
                ]
            )
    return coords


def expand_temp_lists(atom_list):
    """
    Expand the list of atoms to delete. If a "-" is used to denote all atom indexes between
    the numbers split it out to alist ionlcuding all those numbers.

    Args:
        atom_list: List of strings provided by the user.

    Return:
        A list of ints
    """

    expanded_atom_list = []
    for i in range(len(atom_list)):
        if "-" in str(
            atom_list[i]
        ):  # Use the '-' symbol to denote all atoms between the former and latter
            line = str(atom_list[i]).split("-")
            start = int(line[0])
            for j in range(int(line[0]), int(line[1]) + 1):
                expanded_atom_list.append(j)
        if "-" not in str(
            atom_list[i]
        ):  # If not one of the special symbols just add it to the final list
            expanded_atom_list.append(int(atom_list[i]))
    return expanded_atom_list


def merge_lists(child_remove_atoms, child_atoms, parent_remove_atoms, parent_atoms):
    """
    Removes specific indexes from two lists of tuples and combines them into a single list.

    Args:
        - child_remove_atoms (list): A list of indexes to be removed from the child list.
        - child_atoms (list): The initial child list of tuples.
        - parent_remove_atoms (list): A list of indexes or values to be removed from the parent list based on their position in this list.
        - parent_atoms (list): The initial parent list of tuples.

    Returns:
        - result_list (list): A combined and filtered version of both lists as a single output.
    """
    # Filter out specific values and indexes from child_atoms
    remaining_child_atoms = [
        child_tuple
        for idx, child_tuple in enumerate(child_atoms)
        if idx + 1 not in child_remove_atoms
    ]

    # Remove items based on their position specified by parent_remove_atoms
    remaining_parent_atoms = [
        parent_tuple
        for idx, parent_tuple in enumerate(parent_atoms)
        if idx + 1 not in parent_remove_atoms
    ]

    return remaining_child_atoms + remaining_parent_atoms


def main():
    """
    A tool to align two molecules to each other based on a list of atoms provided for both. (Must be equal in length.)
    Also removes atoms provided by the user and creates an aligned and mapped output.
    """

    # Get parsed arguments
    args = parse_args()

    if len(args.parent_align) != len(args.child_align):
        print("Error: Number of Parent and Child alignment atoms do not match.")
        exit(1)

    # Get the path to the input file
    pst = etastructure.ChemicalStructure.load_xyz(args.parent_file)
    cst = etastructure.ChemicalStructure.load_xyz(args.child_file)

    # Pull coordinates for allignemnt
    pcoords = pull_xyz_coords(pst, args.parent_align)
    ccoords = pull_xyz_coords(cst, args.child_align)
    accoords = pull_xyz_coords(cst, args.child_align, all=True)

    # Run alignment
    aligned_child, rms = quatfit.fit_fragment(accoords, ccoords, pcoords)

    print(("Alignment RMSD: {}Å".format(rms)))

    # Save aligned coordinates back to the child molecule
    for atomindex, coord in enumerate(aligned_child):
        cst.atoms[atomindex]["x"] = coord[0]
        cst.atoms[atomindex]["y"] = coord[1]
        cst.atoms[atomindex]["z"] = coord[2]

    # Export the aligned child before atom deletion
    etastructure.ChemicalStructure.export_xyz(cst, "aligned.xyz")

    # Create a new ChemicalStructure
    newst = etastructure.ChemicalStructure()

    # Remove delete atoms from child and parent and merge them to one list.
    newatoms = merge_lists(
        expand_temp_lists(args.child_delete),
        cst.atoms,
        expand_temp_lists(args.parent_delete),
        pst.atoms,
    )

    # Add merged list of atoms to the new ChemicalStructure
    for atom in newatoms:
        newst.add_atom(atom["element"], atom["x"], atom["y"], atom["z"])

    # Export the aligned and merged structure
    etastructure.ChemicalStructure.export_xyz(newst, "replaced.xyz")

    # Check if verbose output is enabled
    if args.verbose:
        print("Verbose mode: ON")

    # Potential to use a map file in the future
    if args.map_file:
        data = read_input_file(args.map_file)
        print(data)


if __name__ == "__main__":
    main()
