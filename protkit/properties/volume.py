#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `Volume` to calculate the volume of a protein,
chain, residue or sequence.

Volume is calculated by summing the volumes of residues in a chain, protein
or sequence. For residues, current volume assignments are fixed and based on
values as provided by the IMGT. In future versions,
volume calculations could be based on side-chain conformations of the residues.

See https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
"""

from typing import List

from protkit.structure.residue import Residue
from protkit.structure.chain import Chain
from protkit.structure.protein import Protein
from protkit.seq.sequence import Sequence


class Volume:
    # VOLUME values for residues as provided by IMGT. See
    # https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
    # Also see:
    # Pontius, J., Richelle, J., Wodak, S.J. (1996) Deviations from standard
    # atomic volumes as a quality measure for protein crystal structures.
    # J. Mol. Biol. 264, 121-136.
    # Volumes can be assessed by one or three-letter codes.
    VOLUME = {
        "ALA": 88.6,
        "ARG": 173.4,
        "ASN": 114.1,
        "ASP": 111.1,
        "CYS": 108.5,
        "GLN": 143.8,
        "GLU": 138.4,
        "GLY": 60.1,
        "HIS": 153.2,
        "ILE": 166.7,
        "LEU": 166.7,
        "LYS": 168.6,
        "MET": 162.9,
        "PHE": 189.9,
        "PRO": 112.7,
        "SER": 89.0,
        "THR": 116.1,
        "TRP": 227.8,
        "TYR": 193.6,
        "VAL": 140.0,

        "A": 88.6,
        "R": 173.4,
        "N": 114.1,
        "D": 111.1,
        "C": 108.5,
        "Q": 143.8,
        "E": 138.4,
        "G": 60.1,
        "H": 153.2,
        "I": 166.7,
        "L": 166.7,
        "K": 168.6,
        "M": 162.9,
        "F": 189.9,
        "P": 112.7,
        "S": 89.0,
        "T": 116.1,
        "W": 227.8,
        "Y": 193.6,
        "V": 140.0
    }

    # Volume class as assigned by IMGT. See
    # https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html?
    # The classes are:
    # VS: Very Small [60-90]: G, A, S
    # S: Small [108-117]: C, D, N, P, T
    # M: Medium [138-154]: E, V, Q, H
    # L: Large [162-174]: M, I, L, K, R
    # VL: Very Large [189-228]: F, Y, W
    UNDEFINED = 0
    VS = 1
    S = 2
    M = 3
    L = 4
    VL = 5

    VOLUME_CLASS_STRING = [
        "Undefined",
        "Very Small",
        "Small",
        "Medium",
        "Large",
        "Very Large"
    ]

    VOLUME_CLASS = {
        "ALA": VS,
        "ARG": L,
        "ASN": S,
        "ASP": S,
        "CYS": S,
        "GLN": M,
        "GLU": M,
        "GLY": VS,
        "HIS": M,
        "ILE": L,
        "LEU": L,
        "LYS": L,
        "MET": L,
        "PHE": VL,
        "PRO": S,
        "SER": VS,
        "THR": S,
        "TRP": VL,
        "TYR": VL,
        "VAL": M,

        "A": VS,
        "R": L,
        "N": S,
        "D": S,
        "C": S,
        "Q": M,
        "E": M,
        "G": VS,
        "H": M,
        "I": L,
        "L": L,
        "K": L,
        "M": L,
        "F": VL,
        "P": S,
        "S": VS,
        "T": S,
        "W": VL,
        "Y": VL,
        "V": M
    }

    @staticmethod
    def volume_of_residue(residue: Residue,
                          assign_attribute: bool = False,
                          key: str = "volume") -> float:
        """
        Get the volume of a residue.

        Args:
            residue (Residue): The residue for which to calculate the volume.
            assign_attribute (bool): If True, the volume will be assigned as an attribute to the residue.
            key (str): The key to use for the attribute.

        Returns:
            float: A float representing the volume of the residue.
        """
        volume = Volume.VOLUME.get(residue.residue_type, 0.0)
        if assign_attribute:
            residue.set_attribute(key, volume)
        return volume

    @staticmethod
    def volume_class_of_residue(residue: Residue,
                                assign_attribute: bool = False,
                                key: str = "volume_class") -> int:
        """
        Get the volume class of a residue.

        Args:
            residue (Residue): The residue for which to calculate the volume class.
            assign_attribute (bool): If True, the volume class will be assigned as an attribute to the residue.
            key (str): The key to use for the attribute.

        Returns:
            int: A string representing the volume class of the residue.
        """
        volume_class = Volume.VOLUME_CLASS.get(residue.residue_type, Volume.UNDEFINED)
        if assign_attribute:
            residue.set_attribute(key, volume_class)
        return volume_class

    @staticmethod
    def volume_of_chain(chain: Chain,
                        assign_attribute: bool = False,
                        key: str = "volume") -> float:
        """
        Get the volume of a chain.

        Args:
            chain (Chain): The chain for which to calculate the volume.
            assign_attribute (bool): If True, the volume will be assigned as an attribute to the chain.
            key (str): The key to use for the attribute.

        Returns:
            float: A float representing the volume of the chain.
        """
        volume = sum(Volume.volume_of_residue(residue, assign_attribute=assign_attribute, key=key) for residue in chain.residues)
        if assign_attribute:
            chain.set_attribute(key, volume)
        return volume

    @staticmethod
    def volume_classes_of_chain(chain: Chain,
                                assign_attribute: bool = False,
                                key: str = "volume_class") -> List[int]:
        """
        Get the volume classes of a chain.

        Args:
            chain (Chain): The chain for which to calculate the volume classes.
            assign_attribute (bool): If True, the volume classes will be assigned as an attribute to the chain.
            key (str): The key to use for the attribute.

        Returns:
            list: A list of volume classes of the chain.
        """
        volume_classes = [Volume.volume_class_of_residue(residue, assign_attribute=assign_attribute, key=key) for residue in chain.residues]
        # if assign_attribute:
        #     chain.set_attribute(key, volume_classes)
        return volume_classes

    @staticmethod
    def volume_of_protein(protein: Protein,
                          assign_attribute: bool = False,
                          key: str = "volume") -> float:
        """
        Get the volume of a protein.

        Args:
            protein (Protein): The protein for which to calculate the volume.
            assign_attribute (bool): If True, the volume will be assigned as an attribute to the protein.
            key (str): The key to use for the attribute.

        Returns:
            float: A float representing the volume of the protein.
        """
        volume = sum(Volume.volume_of_chain(chain, assign_attribute=assign_attribute, key=key) for chain in protein.chains)
        if assign_attribute:
            protein.set_attribute(key, volume)
        return volume

    @staticmethod
    def volume_classes_of_protein(protein: Protein,
                                  assign_attribute: bool = False,
                                  key: str = "volume_class") -> List[List[int]]:
        """
        Get the volume classes of a protein.

        Args:
            protein (Protein): The protein for which to calculate the volume classes.
            assign_attribute (bool): If True, the volume classes will be assigned as an attribute to the protein.
            key (str): The key to use for the attribute.

        Returns:
            list: A list of volume classes of the protein.
        """
        volume_classes = [Volume.volume_classes_of_chain(chain, assign_attribute=assign_attribute, key=key) for chain in protein.chains]
        # if assign_attribute:
        #     protein.set_attribute(key, volume_classes)
        return volume_classes

    @staticmethod
    def volume_of_sequence(sequence: Sequence,
                           assign_attribute: bool = False,
                           key: str = "volume") -> float:
        """
        Get the volume of a sequence.

        Args:
            sequence (Sequence): The sequence for which to calculate the volume.
            assign_attribute (bool): If True, the volume will be assigned as an attribute to the sequence.
            key (str): The key to use for the attribute.

        Returns:
            float: A float representing the volume of the sequence.
        """
        volume = sum(Volume.VOLUME.get(residue, 0.0) for residue in sequence)
        if assign_attribute:
            sequence.set_attribute(key, volume)
        return volume

    @staticmethod
    def volume_classes_of_sequence(sequence: Sequence,
                                   assign_attribute: bool = False,
                                   key: str = "volume_class") -> List[int]:
        """
        Get the volume classes of a sequence.

        Args:
            sequence (Sequence): The sequence for which to calculate the volume classes.
            assign_attribute (bool): If True, the volume classes will be assigned as an attribute to the sequence.
            key (str): The key to use for the attribute.

        Returns:
            list: A list of volume classes of the sequence.
        """
        volume_classes = [Volume.VOLUME_CLASS.get(residue, Volume.UNDEFINED) for residue in sequence]
        if assign_attribute:
            sequence.set_attribute(key, volume_classes)
        return volume_classes
