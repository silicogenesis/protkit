#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `Polarity` to calculate the polarity of a residue,
chain, protein or sequence.

Polarity values related to residues are defined at:
https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html#charge
"""

from typing import List

from protkit.structure.residue import Residue
from protkit.structure.chain import Chain
from protkit.structure.protein import Protein
from protkit.seq.sequence import Sequence


class Polarity:
    UNDEFINED = 0
    POLAR = 1
    NON_POLAR = 2

    # Residues are classified as either polar or non-polar.
    # See: https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
    # Polar residues: R, N, D, Q, E, H, K, S, T, Y
    # Non-polar residues: A, C, F, G, I, L, M, P, V, W

    POLARITY_STRING = [
        "Undefined",
        "Polar",
        "Non-polar"
    ]

    POLARITY = {
        "ALA": NON_POLAR,
        "ARG": POLAR,
        "ASN": POLAR,
        "ASP": POLAR,
        "CYS": NON_POLAR,
        "GLN": POLAR,
        "GLU": POLAR,
        "GLY": NON_POLAR,
        "HIS": POLAR,
        "ILE": NON_POLAR,
        "LEU": NON_POLAR,
        "LYS": POLAR,
        "MET": NON_POLAR,
        "PHE": NON_POLAR,
        "PRO": NON_POLAR,
        "SER": POLAR,
        "THR": POLAR,
        "TRP": NON_POLAR,
        "TYR": POLAR,
        "VAL": NON_POLAR,

        "A": NON_POLAR,
        "R": POLAR,
        "N": POLAR,
        "D": POLAR,
        "C": NON_POLAR,
        "E": POLAR,
        "Q": POLAR,
        "G": NON_POLAR,
        "H": POLAR,
        "I": NON_POLAR,
        "L": NON_POLAR,
        "K": POLAR,
        "M": NON_POLAR,
        "F": NON_POLAR,
        "P": NON_POLAR,
        "S": POLAR,
        "T": POLAR,
        "W": NON_POLAR,
        "Y": POLAR,
        "V": NON_POLAR
    }

    @staticmethod
    def polarity_of_residue(residue: Residue,
                            assign_attribute: bool = False,
                            key: str = "polarity") -> int:
        """
        Returns the polarity of the residue.

        Args:
            residue (Residue): The residue for which to determine the polarity.
            assign_attribute (bool): Whether to assign the polarity to the residue.
            key (str): The key to use for the attribute.

        Returns:
            int: The polarity of the residue.
        """
        polarity = Polarity.POLARITY.get(residue.residue_type, Polarity.UNDEFINED)
        if assign_attribute:
            residue.set_attribute(key, polarity)
        return polarity

    @staticmethod
    def polarities_of_chain(chain: Chain,
                            assign_attribute: bool = False,
                            key: str = "polarity") -> List[int]:
        """
        Returns the polarities of the residues of a chain.

        Args:
            chain (Chain): The chain for which to determine the polarity.
            assign_attribute (bool): Whether to assign the polarity to the chain.
            key (str): The key to use for the attribute.

        Returns:
            int: The polarity of the chain.
        """
        polarities = [Polarity.polarity_of_residue(residue, assign_attribute=assign_attribute, key=key) for residue in chain.residues]
        if assign_attribute:
            chain.set_attribute(key, polarities)
        return polarities

    @staticmethod
    def polarities_of_protein(protein: Protein,
                              assign_attribute: bool = False,
                              key: str = "polarity") -> List[List[int]]:
        """
        Returns the polarities of the residues in the protein.

        Args:
            protein (Protein): The protein for which to determine the polarities.
            assign_attribute (bool): Whether to assign the polarities to the protein.
            key (str): The key to use for the attribute.

        Returns:
            List[List[int]]: The polarities of the residues in the protein.
        """
        polarities = [Polarity.polarities_of_chain(chain, assign_attribute=assign_attribute, key=key) for chain in protein.chains]
        if assign_attribute:
            protein.set_attribute(key, polarities)
        return polarities

    @staticmethod
    def polarities_of_sequence(sequence: Sequence,
                                 assign_attribute: bool = False,
                                 key: str = "polarity") -> List[int]:
        """
        Returns the polarities of the residues in the sequence.

        Args:
            sequence (Sequence): The sequence for which to determine the polarities.
            assign_attribute (bool): Whether to assign the polarities to the sequence.
            key (str): The key to use for the attribute.

        Returns:
            List[int]: The polarities of the residues in the sequence.
        """
        polarities = [Polarity.POLARITY.get(residue, Polarity.UNDEFINED) for residue in sequence]
        if assign_attribute:
            sequence.set_attribute(key, polarities)
        return polarities
