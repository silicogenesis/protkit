#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class to calculate bond angles `BondAngles` to calculate bond angles
in a residue, chain or protein.
"""

from typing import List, Dict, Tuple

from protkit.structure.atom import Atom
from protkit.structure.residue import Residue
from protkit.structure.chain import Chain
from protkit.structure.protein import Protein

from protkit.geometry.math import Math

class BondAngles:
    HEAVY_ATOM_BOND_ANGLES = {
        "ALA": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB")
        },
        "ARG": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "CG"),
            ("CB", "CG", "CD"),
            ("CG", "CD", "NE"),
            ("CD", "NE", "CZ"),
            ("NE", "CZ", "NH1"),
            ("NE", "CZ", "NH2"),
            ("NH1", "CZ", "NH2")
        },
        "ASN": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "CG"),
            ("CB", "CG", "OD1"),
            ("CB", "CG", "ND2"),
            ("ND2", "CG", "OD1")
        },
        "ASP": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "CG"),
            ("CB", "CG", "OD1"),
            ("CB", "CG", "OD2"),
            ("OD1", "CG", "OD2")
        },
        "CYS": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "SG")
        },
        "GLU": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "CG"),
            ("CB", "CG", "CD"),
            ("CG", "CD", "OE1"),
            ("CG", "CD", "OE2"),
            ("OE1", "CD", "OE2")
        },
        "GLN": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "CG"),
            ("CB", "CG", "CD"),
            ("CG", "CD", "OE1"),
            ("CG", "CD", "NE2"),
            ("OE1", "CD", "NE2")
        },
        "GLY": {
            ("N", "CA", "C"),
            ("CA", "C", "O")
        },
        "HIS": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "CG"),
            ("CB", "CG", "ND1"),
            ("CB", "CG", "CD2"),
            ("CG", "ND1", "CE1"),
            ("ND1", "CE1", "NE2"),
            ("CE1", "NE2", "CD2"),
            ("NE2", "CD2", "CG"),
            ("CD2", "CG", "ND1")
        },
        "ILE": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "CG1"),
            ("CB", "CG1", "CD1"),
            ("CA", "CB", "CG2"),
            ("CG1", "CB", "CG2")
        },
        "LEU": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "CG"),
            ("CB", "CG", "CD1"),
            ("CB", "CG", "CD2"),
            ("CD1", "CG", "CD2")
        },
        "LYS": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "CG"),
            ("CB", "CG", "CD"),
            ("CG", "CD", "CE"),
            ("CD", "CE", "NZ")
        },
        "MET": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "CG"),
            ("CB", "CG", "SD"),
            ("CG", "SD", "CE")
        },
        "PHE": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "CG"),
            ("CB", "CG", "CD1"),
            ("CB", "CG", "CD2"),
            ("CD1", "CG", "CD2"),
            ("CG", "CD1", "CE1"),
            ("CD1", "CE1", "CZ"),
            ("CE1", "CZ", "CE2"),
            ("CZ", "CE2", "CD2"),
            ("CE2", "CD2", "CG")
        },
        "PRO": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "CG"),
            ("CB", "CG", "CD"),
            ("CG", "CD", "N"),
            ("CA", "N", "CD")
        },
        "SER": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "OG")
        },
        "THR": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "OG1"),
            ("CA", "CB", "CG2"),
            ("OG1", "CB", "CG2")
        },
        "TRP": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "CG"),
            ("CB", "CG", "CD1"),
            ("CB", "CG", "CD2"),
            ("CD1", "CG", "CD2"),
            ("CG", "CD1", "NE1"),
            ("CD1", "NE1", "CE2"),
            ("NE1", "CE2", "CD2"),
            ("CE2", "CD2", "CG"),
            ("CG", "CD2", "CE3"),
            ("NE1", "CE2", "CZ2"),
            ("CE3", "CD2", "CE2"),
            ("CD2", "CE2", "CZ2"),
            ("CE2", "CZ2", "CH2"),
            ("CZ2", "CH2", "CZ3"),
            ("CH2", "CZ3", "CE3"),
            ("CZ3", "CE3", "CD2")
        },
        "TYR": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "CG"),
            ("CB", "CG", "CD1"),
            ("CB", "CG", "CD2"),
            ("CD1", "CG", "CD2"),
            ("CG", "CD1", "CE1"),
            ("CG", "CD2", "CE2"),
            ("CD1", "CE1", "CZ"),
            ("CD2", "CE2", "CZ"),
            ("CE1", "CZ", "CE2"),
            ("CE1", "CZ", "OH"),
            ("CE2", "CZ", "OH")
        },
        "VAL": {
            ("N", "CA", "C"),
            ("CA", "C", "O"),
            ("N", "CA", "CB"),
            ("C", "CA", "CB"),
            ("CA", "CB", "CG1"),
            ("CA", "CB", "CG2"),
            ("CG1", "CB", "CG2")
        }
    }

    @staticmethod
    def angle(a1: Atom, a2: Atom, a3: Atom) -> float:
        """
        Calculates the angle between three atoms.

        Args:
            a1 (Atom): The first atom.
            a2 (Atom): The second atom.
            a3 (Atom): The third atom.

        Returns:
            float: The angle between the three atoms.
        """
        return Math.angle(a1.x, a1.y, a1.z, a2.x, a2.y, a2.z, a3.x, a3.y, a3.z)

    @staticmethod
    def bond_angles_of_residue(residue: Residue,
                               assign_attribute: bool = False,
                               key: str = "bond_angles") -> Dict[Tuple[str, str, str], float]:
        """
        Returns the bond angles of a residue

        Args:
            residue (Residue): The residue for which to return the bond angles
            assign_attribute (bool): Whether to assign the bond angles to the residue
            key (str): The key to use to access the bond angles

        Returns:
            Dict[Tuple[str, str, str], float]: A dictionary of bond angles
        """
        angles = {}
        defined_angles = BondAngles.HEAVY_ATOM_BOND_ANGLES.get(residue.residue_type, set())

        for (atom1, atom2, atom3) in defined_angles:
            a1 = residue.get_atom(atom1)
            a2 = residue.get_atom(atom2)
            a3 = residue.get_atom(atom3)
            if a1 and a2 and a3:
                angles[(atom1, atom2, atom3)] = BondAngles.angle(a1, a2, a3)
            else:
                angles[(atom1, atom2, atom3)] = None

        if assign_attribute:
            residue.set_attribute(key, angles)

        return angles

    @staticmethod
    def bond_angles_of_chain(chain: Chain,
                             assign_attribute: bool = False,
                             key: str = "bond_angles") -> Dict[str, Dict]:
        """
        Returns the bond angles of a chain

        Args:
            chain (Chain): The chain for which to return the bond angles
            assign_attribute (bool): Whether to assign the bond angles to the chain
            key (str): The key to use to access the bond angles

        Returns:
            Dict[str, Dict]: A dictionary of bond angles.
        """
        angles = {}

        for residue in chain.residues:
            angles[residue.residue_code] = BondAngles.bond_angles_of_residue(residue, assign_attribute=assign_attribute, key=key)

        return angles

    @staticmethod
    def bond_angles_of_protein(protein: Protein,
                               assign_attribute: bool = False,
                               key: str = "bond_angles") -> Dict[str, Dict]:
        """
        Returns the bond angles of a protein

        Args:
            protein (Protein): The protein for which to return the bond angles
            assign_attribute (bool): Whether to assign the bond angles to the protein
            key (str): The key to use to access the bond angles

        Returns:
            Dict[str, Dict]: A dictionary of bond angles.
        """
        angles = {}

        for chain in protein.chains:
            angles[chain.chain_id] = BondAngles.bond_angles_of_chain(chain, assign_attribute=assign_attribute, key=key)

        return angles
