import pytest

from EtaLib.etaatom import Atom, AtomZMAT, InputArguments


def test_atom_initialization():
    atom = Atom()
    assert atom.e == ""
    assert atom.en == 0
    assert atom.x == 0.0
    assert atom.neighbors == []

def test_atom_zmat_initialization():
    zmat = AtomZMAT()
    assert zmat.e == ""
    assert zmat.dist == 0.0

def test_input_arguments_initialization():
    args = InputArguments()
    assert args.ifile == ""
    assert args.charge == 0
    assert args.multi == 1
