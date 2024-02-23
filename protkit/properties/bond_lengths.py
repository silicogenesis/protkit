#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `BondLengths` to calculate the bond lengths in a residue,
chain or protein.
"""

from typing import Dict, List, Tuple, Optional

from protkit.structure.atom import Atom
from protkit.structure.residue import Residue
from protkit.structure.chain import Chain
from protkit.structure.protein import Protein

from protkit.geometry.math import Math


class BondLengths:
    HEAVY_ATOM_BONDS = {
        "ALA": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB")
        ],
        "ASN": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "OD1"),
            ("CG", "ND2")
        ],
        "ASP": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "OD1"),
            ("CG", "OD2")
        ],
        "ARG": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "CD"),
            ("CD", "NE"),
            ("NE", "CZ"),
            ("CZ", "NH1"),
            ("CZ", "NH2")
        ],
        "CYS": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "SG")
        ],
        "GLN": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "CD"),
            ("CD", "OE1"),
            ("CD", "NE2")
        ],
        "GLY": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O")
        ],
        "GLU": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "CD"),
            ("CD", "OE1"),
            ("CD", "OE2")
        ],
        "HIS": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "ND1"),
            ("CG", "CD2"),
            ("ND1", "CE1"),
            ("CD2", "NE2"),
            ("CE1", "NE2")
        ],
        "ILE": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "CG1"),
            ("CB", "CG2"),
            ("CG1", "CD1")
        ],
        "LEU": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "CD1"),
            ("CG", "CD2")
        ],
        "LYS": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "CD"),
            ("CD", "CE"),
            ("CE", "NZ")
        ],
        "MET": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "SD"),
            ("SD", "CE")
        ],
        "PHE": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "CD1"),
            ("CG", "CD2"),
            ("CD1", "CE1"),
            ("CD2", "CE2"),
            ("CE1", "CZ"),
            ("CE2", "CZ")
        ],
        "PRO": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "CD"),
            ("N", "CD"),
        ],
        "SER": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "OG")
        ],
        "THR": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "OG1"),
            ("CB", "CG2")
        ],
        "TRP": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "CD1"),
            ("CG", "CD2"),
            ("CD1", "NE1"),
            ("CD2", "CE2"),
            ("CD2", "CE3"),
            ("NE1", "CE2"),
            ("CE2", "CZ2"),
            ("CE3", "CZ3"),
            ("CZ2", "CH2"),
            ("CZ3", "CH2"),
        ],
        "TYR": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "CG"),
            ("CG", "CD1"),
            ("CG", "CD2"),
            ("CD1", "CE1"),
            ("CD2", "CE2"),
            ("CE1", "CZ"),
            ("CE2", "CZ"),
            ("CZ", "OH")
        ],
        "VAL": [
            ("N", "CA"),
            ("CA", "C"),
            ("C", "O"),
            ("CA", "CB"),
            ("CB", "CG1"),
            ("CB", "CG2")
        ]
    }

    @staticmethod
    def atom_distance(atom1: Atom, atom2: Atom) -> float:
        """
        Calculates the Euclidean distance between two atoms.

        Args:
            atom1 (Atom): The first atom.
            atom2 (Atom): The second atom.

        Returns:
            float: The distance between the two atoms.

        Raises:
            None
        """
        return Math.euclidean_distance(atom1.x, atom1.y, atom1.z, atom2.x, atom2.y, atom2.z)


    @staticmethod
    def bond_lengths_of_residue(residue: Residue,
                                assign_attribute: bool = False,
                                key: str = "bond_lengths") -> Dict[Tuple[str, str], float]:
        """
        Returns the bond lengths of the residue.

        Args:
            residue (Residue): The residue for which to return the bond lengths.
            assign_attribute (bool): Whether to assign the bond lengths to the residue.
            key (str): The key to use for the attribute.

        Returns:
            Dict[Tuple[str, str], float]: A dictionary containing the bond lengths of the residue.
        """

        bonds = {}
        defined_bonds = BondLengths.HEAVY_ATOM_BONDS.get(residue.residue_type, [])

        for (atom1, atom2) in defined_bonds:
            a1 = residue.get_atom(atom1)
            a2 = residue.get_atom(atom2)
            if a1 and a2:
                bonds[(atom1, atom2)] = BondLengths.atom_distance(a1, a2)
            else:
                bonds[(atom1, atom2)] = None

        if assign_attribute:
            residue.set_attribute(key, bonds)

        return bonds

    @staticmethod
    def bond_lengths_of_chain(chain: Chain,
                              assign_attribute: bool = False,
                              key: str = "bond_lengths") -> Dict[str, Dict]:
        """
        Returns the bond lengths of the chain.

        Args:
            chain (Chain): The chain for which to return the bond lengths.
            assign_attribute (bool): Whether to assign the bond lengths to the chain.
            key (str): The key to use for the attribute.

        Returns:
            Dict[str, Dict]: A dictionary containing the bond lengths of the chain.
        """

        bonds = {}
        for residue in chain.residues:
            bonds[residue.residue_code] = BondLengths.bond_lengths_of_residue(residue, assign_attribute=assign_attribute, key=key)

        return bonds

    @staticmethod
    def bond_lengths_of_protein(protein: Protein,
                                assign_attribute: bool = False,
                                key: str = "bond_lengths") -> Dict[str, Dict]:
        """
        Returns the bond lengths of the protein.

        Args:
            protein (Protein): The protein for which to return the bond lengths.
            assign_attribute (bool): Whether to assign the bond lengths to the protein.
            key (str): The key to use for the attribute.

        Returns:
            Dict[str, Dict]: A dictionary containing the bond lengths of the protein.
        """

        bonds = {}
        for chain in protein.chains:
            bonds[chain.chain_id] = BondLengths.bond_lengths_of_chain(chain, assign_attribute=assign_attribute, key=key)

        return bonds

    @staticmethod
    def peptide_bond_length(residue1: Residue, residue2: Residue) -> Optional[float]:
        """
        Returns the peptide bond length between two residues.  Specifically, the bond length
        between the C atom of the first residue and the N atom of the second residue.

        Args:
            residue1 (Residue): The first residue.
            residue2 (Residue): The second residue.

        Returns:
            float: The bond length between the two residues.  If the atoms are not found, None is returned.
        """

        a1 = residue1.get_atom("C")
        a2 = residue2.get_atom("N")
        if a1 and a2:
            return BondLengths.atom_distance(a1, a2)
        else:
            return None

    @staticmethod
    def peptide_bond_lengths_of_chain(chain: Chain,
                                      assign_attribute: bool = False,
                                      key: str = "peptide_bond_lengths") -> List[float]:
        """
        Returns the peptide bond lengths of the chain.

        Args:
            chain (Chain): The chain for which to return the peptide bond lengths.
            assign_attribute (bool): Whether to assign the peptide bond lengths to the chain.
            key (str): The key to use for the attribute.

        Returns:
            dict: A dictionary containing the peptide bond lengths of the chain.
        """

        bonds = []
        for i in range(chain.num_residues - 1):
            bonds.append(BondLengths.peptide_bond_length(chain.get_residue(i), chain.get_residue(i + 1)))

        if assign_attribute:
            chain.set_attribute(key, bonds)

        return bonds

    @staticmethod
    def peptide_bond_lengths_of_protein(protein: Protein,
                                        assign_attribute: bool = False,
                                        key: str = "peptide_bond_lengths") -> Dict[str, List[float]]:
        """
        Returns the peptide bond lengths of the protein.

        Args:
            protein (Protein): The protein for which to return the peptide bond lengths.
            assign_attribute (bool): Whether to assign the peptide bond lengths to the protein.
            key (str): The key to use for the attribute.

        Returns:
            dict: A dictionary containing the peptide bond lengths of the protein.
        """

        bonds = {}
        for chain in protein.chains:
            bonds[chain.chain_id] = BondLengths.peptide_bond_lengths_of_chain(chain, assign_attribute=assign_attribute, key=key)

        return bonds
