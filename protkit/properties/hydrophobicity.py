#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `Hydrophobicity` to calculate the hydrophobicity and
hydrophobicity class of a residue, chain, protein or sequence.

Hydrophobicity is calculated based on the Kyte-Doolittle scale.

For chains, the hydrophobicity is
the sum of the hydrophobicity values of the residues in the chain.  For
proteins, the hydrophobicity is the sum of the hydrophobicity values of
the chains in the protein.

Calculated values can be added as attributes to the respective objects.

For more information, see:
https://en.wikipedia.org/wiki/Hydrophobicity_scales#Kyte%E2%80%93Doolittle_scale

The Kyte-Doolittle scale was defined in:
Kyte, J. and Doolittle, R.F. (1982) A simple method for displaying the
hydropathic character of a protein. J. Mol. Biol. 157, 105-132.

Also see:
https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
"""

from typing import List

from protkit.structure.residue import Residue
from protkit.structure.chain import Chain
from protkit.structure.protein import Protein
from protkit.seq.sequence import Sequence


class Hydrophobicity:
    # Hydrophobicity values on the Kyte-Doolittle scale.
    HYDROPHOBICITY = {
        "ALA": 1.8,
        "ARG": -4.5,
        "ASN": -3.5,
        "ASP": -3.5,
        "CYS": 2.5,
        "GLN": -3.5,
        "GLU": -3.5,
        "GLY": -0.4,
        "HIS": -3.2,
        "ILE": 4.5,
        "LEU": 3.8,
        "LYS": -3.9,
        "MET": 1.9,
        "PHE": 2.8,
        "PRO": -1.6,
        "SER": -0.8,
        "THR": -0.7,
        "TRP": -0.9,
        "TYR": -1.3,
        "VAL": 4.2,

        "A": 1.8,
        "R": -4.5,
        "N": -3.5,
        "D": -3.5,
        "C": 2.5,
        "Q": -3.5,
        "E": -3.5,
        "G": -0.4,
        "H": -3.2,
        "I": 4.5,
        "L": 3.8,
        "K": -3.9,
        "M": 1.9,
        "F": 2.8,
        "P": -1.6,
        "S": -0.8,
        "T": -0.7,
        "W": -0.9,
        "Y": -1.3,
        "V": 4.2
    }

    # Hydrophobicity classes as assigned by IMGT. See
    # https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html?
    # The classes are:
    # Hydrophobic: I, V, L, F, C, M, A, W
    # Neutral: G, T, S, Y, P, H
    # Hydrophilic: N, D, Q, E, K, R
    # Note that IMGT assigns Tryptophan (-0.9 on the Kyte-Doolittle scale) to the
    # hydrophobic class, while it is often considered to be neutral.
    # Also note that Histidine (-3.2 on the Kyte-Doolittle scale) is considered to be
    # neutral by IMGT, while it is often considered to be hydrophilic.

    UNDEFINED = 0
    HYDROPHOBIC = 1
    NEUTRAL = 2
    HYDROPHILIC = 3

    HYDROPHOBICITY_CLASS_STRING = [
        "Undefined",
        "Hydrophobic",
        "Neutral",
        "Hydrophilic"
    ]

    HYDROPHOBICITY_CLASS = {
        "ALA": HYDROPHOBIC,
        "ARG": HYDROPHILIC,
        "ASN": HYDROPHILIC,
        "ASP": HYDROPHILIC,
        "CYS": HYDROPHOBIC,
        "GLN": HYDROPHILIC,
        "GLU": HYDROPHILIC,
        "GLY": NEUTRAL,
        "HIS": NEUTRAL,
        "ILE": HYDROPHOBIC,
        "LEU": HYDROPHOBIC,
        "LYS": HYDROPHILIC,
        "MET": HYDROPHOBIC,
        "PHE": HYDROPHOBIC,
        "PRO": NEUTRAL,
        "SER": NEUTRAL,
        "THR": NEUTRAL,
        "TRP": HYDROPHOBIC,
        "TYR": NEUTRAL,
        "VAL": HYDROPHOBIC,

        "A": HYDROPHOBIC,
        "R": HYDROPHILIC,
        "N": HYDROPHILIC,
        "D": HYDROPHILIC,
        "C": HYDROPHOBIC,
        "Q": HYDROPHILIC,
        "E": HYDROPHILIC,
        "G": NEUTRAL,
        "H": NEUTRAL,
        "I": HYDROPHOBIC,
        "L": HYDROPHOBIC,
        "K": HYDROPHILIC,
        "M": HYDROPHOBIC,
        "F": HYDROPHOBIC,
        "P": NEUTRAL,
        "S": NEUTRAL,
        "T": NEUTRAL,
        "W": HYDROPHOBIC,
        "Y": NEUTRAL,
        "V": HYDROPHOBIC
    }

    @staticmethod
    def hydrophobicity_of_residue(residue: Residue,
                                  assign_attribute: bool = False,
                                  key: str = "hydrophobicity") -> float:
        """
        Returns the hydrophobicity value of the residue.

        Args:
            residue (Residue): The residue for which the hydrophobicity value will be returned.
            assign_attribute (bool): If True, the hydrophobicity value will be added as an attribute to the residue.
            key (str): The name of the attribute that will be added to the residue.

        Returns:
            float: The hydrophobicity value of the residue.
        """
        hydrophobicity = Hydrophobicity.HYDROPHOBICITY.get(residue.residue_type, 0.0)
        if assign_attribute:
            residue.set_attribute(key, hydrophobicity)
        return hydrophobicity

    @staticmethod
    def hydrophobicity_class_of_residue(residue: Residue,
                                        assign_attribute: bool = False,
                                        key: str = "hydrophobicity_class") -> int:
        """
        Returns the hydrophobicity class of the residue.

        Args:
            residue (Residue): The residue for which the hydrophobicity class will be returned.
            assign_attribute (bool): If True, the hydrophobicity class will be added as an attribute to the residue.
            key (str): The name of the attribute that will be added to the residue.

        Returns:
            int: The hydrophobicity class of the residue.
        """
        hydrophobicity_class = Hydrophobicity.HYDROPHOBICITY_CLASS.get(residue.residue_type, Hydrophobicity.UNDEFINED)
        if assign_attribute:
            residue.set_attribute(key, hydrophobicity_class)
        return hydrophobicity_class

    @staticmethod
    def hydrophobicity_of_chain(chain: Chain,
                                assign_attribute: bool = False,
                                key: str = "hydrophobicity") -> float:
        """
        Returns the hydrophobicity value of the chain.

        Args:
            chain (Chain): The chain for which the hydrophobicity value will be returned.
            assign_attribute (bool): If True, the hydrophobicity value will be added as an attribute to the chain.
            key (str): The name of the attribute that will be added to the chain.

        Returns:
            float: The hydrophobicity value of the chain.
        """
        hydrophobicity = sum(Hydrophobicity.hydrophobicity_of_residue(residue, assign_attribute=assign_attribute, key=key) for residue in chain.residues)
        if assign_attribute:
            chain.set_attribute(key, hydrophobicity)
        return hydrophobicity

    @staticmethod
    def hydrophobicity_classes_of_chain(chain: Chain,
                                        assign_attribute: bool = False,
                                        key: str = "hydrophobicity_class") -> List[int]:
        """
        Returns the hydrophobicity class of the chain.

        Args:
            chain (Chain): The chain for which the hydrophobicity class will be returned.
            assign_attribute (bool): If True, the hydrophobicity class will be added as an attribute to the chain.
            key (str): The name of the attribute that will be added to the chain.

        Returns:
            list: The hydrophobicity classes of the residues of the chain.
        """
        hydrophobicity_classes = [Hydrophobicity.hydrophobicity_class_of_residue(residue, assign_attribute=assign_attribute, key=key) for residue in chain.residues]
        if assign_attribute:
            chain.set_attribute(key, hydrophobicity_classes)
        return hydrophobicity_classes

    @staticmethod
    def hydrophobicity_of_protein(protein: Protein,
                                  assign_attribute: bool = False,
                                  key: str = "hydrophobicity") -> float:
        """
        Returns the hydrophobicity value of the protein.

        Args:
            protein (Protein): The protein for which the hydrophobicity value will be returned.
            assign_attribute (bool): If True, the hydrophobicity value will be added as an attribute to the protein.
            key (str): The name of the attribute that will be added to the protein.

        Returns
            float: The hydrophobicity value of the protein.
        """
        hydrophobicity = sum(Hydrophobicity.hydrophobicity_of_chain(chain, assign_attribute=assign_attribute, key=key) for chain in protein.chains)
        if assign_attribute:
            protein.set_attribute(key, hydrophobicity)
        return hydrophobicity

    @staticmethod
    def hydrophobicity_classes_of_protein(protein: Protein,
                                          assign_attribute: bool = False,
                                          key: str = "hydrophobicity_class") -> List[List[int]]:
        """
        Returns the hydrophobicity classes of the residues of the protein.

        Args:
            protein (Protein): The protein for which the hydrophobicity class will be returned.
            assign_attribute (bool): If True, the hydrophobicity class will be added as an attribute to the protein.
            key (str): The name of the attribute that will be added to the protein.

        Returns:
            list: The hydrophobicity classes of the residues of the protein.
        """
        hydrophobicity_classes = [Hydrophobicity.hydrophobicity_classes_of_chain(chain, assign_attribute=assign_attribute, key=key) for chain in protein.chains]
        if assign_attribute:
            protein.set_attribute(key, hydrophobicity_classes)
        return hydrophobicity_classes

    @staticmethod
    def hydrophobicity_of_sequence(sequence: Sequence,
                                   assign_attribute: bool = False,
                                   key: str = "hydrophobicity") -> float:
        """
        Returns the hydrophobicity value of the sequence.

        Args:
            sequence (Sequence): The sequence for which the hydrophobicity value will be returned.
            assign_attribute (bool): If True, the hydrophobicity value will be added as an attribute to the sequence.
            key (str): The name of the attribute that will be added to the sequence.

        Returns:
            float: The hydrophobicity value of the sequence.
        """
        hydrophobicity = sum(Hydrophobicity.HYDROPHOBICITY.get(residue, 0.0) for residue in sequence)
        if assign_attribute:
            sequence.set_attribute(key, hydrophobicity)
        return hydrophobicity

    @staticmethod
    def hydrophobicity_classes_of_sequence(sequence: Sequence,
                                           assign_attribute: bool = False,
                                           key: str = "hydrophobicity_class") -> List[int]:
        """
        Returns the hydrophobicity classes of the residues of the sequence.

        Args:
            sequence (Sequence): The sequence for which the hydrophobicity class will be returned.
            assign_attribute (bool): If True, the hydrophobicity class will be added as an attribute to the sequence.
            key (str): The name of the attribute that will be added to the sequence.

        Returns:
            list: The hydrophobicity classes of the residues of the sequence.
        """
        hydrophobicity_classes = [Hydrophobicity.HYDROPHOBICITY_CLASS.get(residue, Hydrophobicity.UNDEFINED) for residue in sequence]
        if assign_attribute:
            sequence.set_attribute(key, hydrophobicity_classes)
        return hydrophobicity_classes
