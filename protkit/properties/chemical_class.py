#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `ChemicalClass` to represent the chemical class of a residue.

IMGT defines the following chemical classes for amino acids:
    1. Aliphatic (A, G, I, L, P, V)
    2. Aromatic (F, W, Y)
    3. Sulfur (C, M)
    4. Hydroxyl (S, T)
    5. Basic (R, H, K)
    6. Acidic (D, E)
    7. Amide (N, Q)

See: https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
"""

from typing import List

from protkit.structure.residue import Residue
from protkit.structure.chain import Chain
from protkit.structure.protein import Protein


class ChemicalClass:
    UNDEFINED = 0
    ALIPHATIC = 1
    AROMATIC = 2
    SULFUR = 3
    HYDROXYL = 4
    BASIC = 5
    ACIDIC = 6
    AMIDE = 7

    CHEMICAL_CLASS_STRING = [
        "Undefined",
        "Aliphatic",
        "Aromatic",
        "Sulfur",
        "Hydroxyl",
        "Basic",
        "Acidic",
        "Amide"
    ]

    CHEMICAL_CLASS = {
        "ALA": ALIPHATIC,
        "ARG": BASIC,
        "ASN": AMIDE,
        "ASP": ACIDIC,
        "CYS": SULFUR,
        "GLN": AMIDE,
        "GLU": ACIDIC,
        "GLY": ALIPHATIC,
        "HIS": BASIC,
        "ILE": ALIPHATIC,
        "LEU": ALIPHATIC,
        "LYS": BASIC,
        "MET": SULFUR,
        "PHE": AROMATIC,
        "PRO": ALIPHATIC,
        "SER": HYDROXYL,
        "THR": HYDROXYL,
        "TRP": AROMATIC,
        "TYR": AROMATIC,
        "VAL": ALIPHATIC,

        "A": ALIPHATIC,
        "R": BASIC,
        "N": AMIDE,
        "D": ACIDIC,
        "C": SULFUR,
        "Q": AMIDE,
        "E": ACIDIC,
        "G": ALIPHATIC,
        "H": BASIC,
        "I": ALIPHATIC,
        "L": ALIPHATIC,
        "K": BASIC,
        "M": SULFUR,
        "F": AROMATIC,
        "P": ALIPHATIC,
        "S": HYDROXYL,
        "T": HYDROXYL,
        "W": AROMATIC,
        "Y": AROMATIC,
        "V": ALIPHATIC
    }

    @staticmethod
    def chemical_class_of_residue(residue: Residue,
                                  assign_attribute: bool = False,
                                  key: str = "chemical_class") -> int:
        """
        Returns the chemical class of the residue.

        Args:
            residue (Residue): The residue for which to determine the chemical class.
            assign_attribute (bool): Whether to assign the chemical class to the residue.
            key (str): The key to use for the attribute.

        Returns:
            int: The chemical class of the residue.
        """
        chemical_class = ChemicalClass.CHEMICAL_CLASS.get(residue.residue_type, ChemicalClass.UNDEFINED)
        if assign_attribute:
            residue.set_attribute(key, chemical_class)
        return chemical_class

    @staticmethod
    def chemical_classes_of_chain(chain: Chain,
                                  assign_attribute: bool = False,
                                  key: str = "chemical_class") -> List[int]:
        """
        Returns the chemical classes of the residues in the chain.

        Args:
            chain (Chain): The chain for which to determine the chemical classes.
            assign_attribute (bool): Whether to assign the chemical classes to the chain.
            key (str): The key to use for the attribute.

        Returns:
            List[int]: The chemical classes of the residues in the chain.
        """
        chemical_classes = [ChemicalClass.chemical_class_of_residue(residue, assign_attribute=assign_attribute, key=key) for residue in chain.residues]
        if assign_attribute:
            chain.set_attribute(key, chemical_classes)
        return chemical_classes

    @staticmethod
    def chemical_classes_of_protein(protein: Protein,
                                    assign_attribute: bool = False,
                                    key: str = "chemical_class") -> List[List[int]]:
        """
        Returns the chemical classes of the residues in the protein.

        Args:
            protein (Protein): The protein for which to determine the chemical classes.
            assign_attribute (bool): Whether to assign the chemical classes to the protein.
            key (str): The key to use for the attribute.

        Returns:
            List[List[int]]: The chemical classes of the residues in the protein.
        """
        chemical_classes = [ChemicalClass.chemical_classes_of_chain(chain, assign_attribute=assign_attribute, key=key) for chain in protein.chains]
        if assign_attribute:
            protein.set_attribute(key, chemical_classes)
        return chemical_classes
