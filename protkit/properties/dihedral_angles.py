#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

from typing import Dict

from protkit.structure.atom import Atom
from protkit.structure.residue import Residue
from protkit.structure.chain import Chain
from protkit.structure.protein import Protein

from protkit.geometry.math import Math

class DihedralAngles:
    # Definitions taken from: http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html
    # See also: https://www.cryst.bbk.ac.uk/PPS95/course/3_geometry/conform.html
    # See also: https://swissmodel.expasy.org/course/text/chapter3.htm#:~:text=Due%20to%20almost%20constant%20bond,cluster%20around%20energetically%20preferred%20conformations.
    # https://leimao.github.io/blog/Dihedral-Angles/
    # Code for dihedral angle calculation from:
    # https://charmm-gui.org/?doc=lecture&module=scientific&lesson=9
    # https://github.com/PDB-REDO/dssp/blob/trunk/src/dssp.cpp
    # https://gist.github.com/hypnopump/30d6bfdb8358d3a57d010c9a501fda56
    DIHEDRAL_ANGLES = {
        "PHI": ("C", "N", "CA", "C"),
        "PSI": ("N", "CA", "C", "N"),
        "OMEGA": ("CA", "C", "N", "CA"),
        # The peptide bond formed by the residues I and I + 1 is assigned to the residue I + 1. The same applies to the omega angle. For that reason no omega angle is assigned to the first residue.

        "ALA": {},
        "ARG": {
            "CHI1": ("N", "CA", "CB", "CG"),
            "CHI2": ("CA", "CB", "CG", "CD"),
            "CHI3": ("CB", "CG", "CD", "NE"),
            "CHI4": ("CG", "CD", "NE", "CZ"),
            "CHI5": ("CD", "NE", "CZ", "NH1")
        },
        "ASN": {
            "CHI1": ("N", "CA", "CB", "CG"),
            "CHI2": ("CA", "CB", "CG", "OD1")
        },
        "ASP": {
            "CHI1": ("N", "CA", "CB", "CG"),
            "CHI2": ("CA", "CB", "CG", "OD1")
        },
        "CYS": {
            "CHI1": ("N", "CA", "CB", "SG")
        },
        "GLU": {
            "CHI1": ("N", "CA", "CB", "CG"),
            "CHI2": ("CA", "CB", "CG", "CD"),
            "CHI3": ("CB", "CG", "CD", "OE1")
        },
        "GLN": {
            "CHI1": ("N", "CA", "CB", "CG"),
            "CHI2": ("CA", "CB", "CG", "CD"),
            "CHI3": ("CB", "CG", "CD", "OE1")
        },
        "GLY": {},
        "HIS": {
            "CHI1": ("N", "CA", "CB", "CG"),
            "CHI2": ("CA", "CB", "CG", "ND1")
        },
        "ILE": {
            "CHI1": ("N", "CA", "CB", "CG1"),
            "CHI2": ("CA", "CB", "CG1", "CD1")  # Website ref calls this CD
        },
        "LEU": {
            "CHI1": ("N", "CA", "CB", "CG"),
            "CHI2": ("CA", "CB", "CG", "CD1")
        },
        "LYS": {
            "CHI1": ("N", "CA", "CB", "CG"),
            "CHI2": ("CA", "CB", "CG", "CD"),
            "CHI3": ("CB", "CG", "CD", "CE"),
            "CHI4": ("CG", "CD", "CE", "NZ")
        },
        "MET": {
            "CHI1": ("N", "CA", "CB", "CG"),
            "CHI2": ("CA", "CB", "CG", "SD"),
            "CHI3": ("CB", "CG", "SD", "CE")
        },
        "PHE": {
            "CHI1": ("N", "CA", "CB", "CG"),
            "CHI2": ("CA", "CB", "CG", "CD1"),
        },
        "PRO": {
            "CHI1": ("N", "CA", "CB", "CG"),
            "CHI2": ("CA", "CB", "CG", "CD"),
        },
        "SER": {
            "CHI1": ("N", "CA", "CB", "OG")
        },
        "THR": {
            "CHI1": ("N", "CA", "CB", "OG1"),
        },
        "TRP": {
            "CHI1": ("N", "CA", "CB", "CG"),
            "CHI2": ("CA", "CB", "CG", "CD1"),
        },
        "TYR": {
            "CHI1": ("N", "CA", "CB", "CG"),
            "CHI2": ("CA", "CB", "CG", "CD1"),
        },
        "VAL": {
            "CHI1": ("N", "CA", "CB", "CG1"),
        }
    }

    @staticmethod
    def dihedral_angle(a1: Atom, a2: Atom, a3: Atom, a4: Atom) -> float:
        """
        Calculates the dihedral angle between four atoms.

        Args:
            a1 (Atom): The first atom.
            a2 (Atom): The second atom.
            a3 (Atom): The third atom.
            a4 (Atom): The fourth atom.

        Returns:
            float: The dihedral angle between the four atoms.
        """
        return Math.dihedral_angle(a1.x, a1.y, a1.z, a2.x, a2.y, a2.z, a3.x, a3.y, a3.z, a4.x, a4.y, a4.z)

    @staticmethod
    def dihedral_angles_of_residue(residue: Residue,
                                   previous_residue: Residue = None,
                                   assign_attribute: bool = False,
                                   key: str = "dihedral_angles") -> Dict[str, float]:
        """
        Returns the dihedral angles of a residue.

        Args:
            residue (Residue): The residue object
            previous_residue (Residue): The previous residue object
            assign_attribute (bool): Whether to assign the dihedral angles as attributes to the residue object
            key (str): The key to use when assigning the dihedral angles as attributes to the residue object

        Returns:
            dict: A dictionary containing the dihedral angles of the residue side chain. If the previous residue
                is defined, the phi, psi and omega angles are also calculated.
        """
        angles = {}
        defined_angles = DihedralAngles.DIHEDRAL_ANGLES.get(residue.residue_type, set())

        if previous_residue:
            n1 = previous_residue.get_atom("N")
            ca1 = previous_residue.get_atom("CA")
            c1 = previous_residue.get_atom("C")
            n2 = residue.get_atom("N")
            ca2 = residue.get_atom("CA")
            c2 = residue.get_atom("C")
            angles["PHI"] = DihedralAngles.dihedral_angle(c1, n2, ca2, c2)
            angles["PSI"] = DihedralAngles.dihedral_angle(n1, ca1, c1, n2)
            angles["OMEGA"] = DihedralAngles.dihedral_angle(ca1, c1, n2, ca2)

        for angle_name, (atom1, atom2, atom3, atom4) in defined_angles.items():
            a1 = residue.get_atom(atom1)
            a2 = residue.get_atom(atom2)
            a3 = residue.get_atom(atom3)
            a4 = residue.get_atom(atom4)
            if a1 and a2 and a3 and a4:
                angles[angle_name] = DihedralAngles.dihedral_angle(a1, a2, a3, a4)
            else:
                angles[angle_name] = None

        if assign_attribute:
            residue.set_attribute(key, angles)

        return angles

    @staticmethod
    def dihedral_angles_of_chain(chain: Chain,
                                 assign_attribute: bool = False,
                                 key: str = "dihedral_angles"):
        """
        Returns the dihedral angles of a chain.

        Args:
            chain (Chain): The chain object
            assign_attribute (bool): Whether to assign the dihedral angles as attributes to the chain object
            key (str): The key to use when assigning the dihedral angles as attributes to the chain object

        Returns:
            dict: A dictionary containing the dihedral angles of the chain
        """
        angles = {}

        previous_residue = None

        for residue in chain.residues:
            angles[residue.residue_code] = DihedralAngles.dihedral_angles_of_residue(
                residue,
                previous_residue=previous_residue,
                assign_attribute=assign_attribute, key=key)
            previous_residue = residue

        return angles

    @staticmethod
    def dihedral_angles_of_protein(protein: Protein,
                                   assign_attribute: bool = False,
                                   key: str = "dihedral_angles"):
        """
        Returns the dihedral angles of a protein.

        Args:
            protein (Protein): The protein object
            assign_attribute (bool): Whether to assign the dihedral angles as attributes to the protein object
            key (str): The key to use when assigning the dihedral angles as attributes to the protein object
        """

        angles = {}

        for chain in protein.chains:
            angles[chain.chain_id] = DihedralAngles.dihedral_angles_of_chain(chain, assign_attribute=assign_attribute, key=key)

        return angles