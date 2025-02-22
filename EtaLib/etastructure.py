#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#   Copyright (c) 2025, Thomas J. L. Mustard, O. Maduka Ogba, Paul Ha-Yeon Cheong
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
import shutil


class ChemicalStructure:
    def __init__(self):
        self.title = None  # title of the molecule
        self.type = None # type of chemical structure (MOLecule, BIOlogical, CRYstal, HELM, etc)
        self.properties = {
            'charge': None,  # charge on the molecule (in elementary charges)
            'multiplicity': None  # multiplicity of the molecule (1 for singlet, etc.)
        }
        self.atoms = {}  # dictionary to store atom data

    def add_    def __init__(self):
        self.title = None  # title of the molecule
        self.type = None # type of chemical structure (MOLecule, BIOlogical, CRYstal, HELM, etc)
        self.properties = {
            'charge': None,  # charge on the molecule (in elementary charges)
            'multiplicity': None  # multiplicity of the molecule (1 for singlet, etc.)
        }
        self.atoms = []  # list of dictionaries


    def add_atom(self, element: str, x: float, y: float, z: float):
        """Adds an atom at the specified position with a specific element."""
        new_atom = {
            'element': str(element)
            'x': float(x),
            'y': float(y),
            'z': float(z)
            }
        self.atoms.append(new_atom)

    def get_atoms(self):
        """Return a list of all atoms"""
        return self.atoms

    def set_property(self, name, value=None):
        """Set the title or other properties of the molecule (including user-defined ones)"""
        if not isinstance(name, str):
            raise ValueError("Property names must be strings")

        self.properties[name] = value

    def set_atom_property(self, index, name, value=None):
        """Set the title or other properties of the molecule (including user-defined ones)"""
        if not isinstance(name, str):
            raise ValueError("Property names must be strings")

        self.atom[index][name] = value


    @classmethod
    def load_xyz(cls, filename):
        """Load XYZ file into a ChemicalStructure instance."""
        with open(filename, 'r') as f:
            lines = [line.strip() for line in f.readlines()]

        if len(lines) < 3 or int(lines[0]) + 2 != len(lines):
            raise ValueError(f"Invalid XYZ file: {filename}")

        cls_title = lines[1]
        atoms_info = lines[3:]

        cs_instance = cls()
        cs_instance.title = cls_title

        for line in atoms_info:
            element, x, y, z = line.split()
            cs_instance.add_atom(element, x, y, z)

        return cs_instance

    def export_xyz(cls, filename, title="Molecule", atoms=None):
        """Export molecule data to an XYZ file."""
        with open(filename, 'w') as f:
            if atoms is None:
                f.write(len(cls.atoms) + "\n")
                if title is None:
                    f.write(cls.title + "\n")
                else:
                    f.write(title + "\n")
                for atom in cls.get_atoms().values():
                    element = list(atom.keys())[0]
                    x, y, z = atom[element].values()
                    line = f"{element} {x:.5f} {y:.5f} {z:.5f}\n"
                    f.write(line)
            else:
                f.write(len(atoms) + "\n")
                if title is None:
                    f.write(cls.title + "\n")
                else:
                    f.write(title + "\n")
                for atom in atoms.values():
                    element = list(atom.keys())[0]
                    x, y, z = atom[element].values()
                    line = f"{element} {x:.5f} {y:.5f} {z:.5f}\n"
                    f.write(line)

        with open(filename, 'a') as f:
            f.write(f"\n{title}")
atom(self, element, x, y, z):
        """Add an atom to the structure"""
        if element in self.atoms:
            raise ValueError(f"Duplicate element {element} found")
        self.atoms[index] = {
            'element': str(element)
            'x': float(x),
            'y': float(y),
            'z': float(z)
        }

    def get_atoms(self):
        """Return a dictionary of all atoms"""
        return dict(self.atoms)

    def set_property(self, name, value=None):
        """Set the title or other properties of the molecule (including user-defined ones)"""
        if not isinstance(name, str):
            raise ValueError("Property names must be strings")

        self.properties[name] = value

    @classmethod
    def load_xyz(cls, filename):
        """Load XYZ file into a ChemicalStructure instance."""
        with open(filename, 'r') as f:
            lines = [line.strip() for line in f.readlines()]

        if len(lines) < 3 or int(lines[0]) + 2 != len(lines):
            raise ValueError(f"Invalid XYZ file: {filename}")

        cls_title = lines[1]
        atoms_info = lines[3:]

        cs_instance = cls()
        cs_instance.title = cls_title

        for line in atoms_info:
            element, x, y, z = line.split()
            cs_instance.add_atom(element, x, y, z)

        return cs_instance

    def export_xyz(cls, filename, title="Molecule", atoms=None):
        """Export molecule data to an XYZ file."""
        with open(filename, 'w') as f:
            if atoms is None:
                f.write(len(cls.atoms) + "\n")
                if title is None:
                    f.write(cls.title + "\n")
                else:
                    f.write(title + "\n")
                for atom in cls.get_atoms().values():
                    element = list(atom.keys())[0]
                    x, y, z = atom[element].values()
                    line = f"{element} {x:.5f} {y:.5f} {z:.5f}\n"
                    f.write(line)
            else:
                f.write(len(atoms) + "\n")
                if title is None:
                    f.write(cls.title + "\n")
                else:
                    f.write(title + "\n")
                for atom in atoms.values():
                    element = list(atom.keys())[0]
                    x, y, z = atom[element].values()
                    line = f"{element} {x:.5f} {y:.5f} {z:.5f}\n"
                    f.write(line)

        with open(filename, 'a') as f:
            f.write(f"\n{title}")

