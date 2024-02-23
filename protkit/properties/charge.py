#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `Charge` to calculate the charge of a residue,
chain, protein or sequence.

Charge values related to residues are defined at:
https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html#charge
"""

from protkit.structure.residue import Residue
from protkit.structure.chain import Chain
from protkit.structure.protein import Protein
from protkit.seq.sequence import Sequence


class Charge:
    # Charges of residues as defined by IMGT.
    # See: https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html#charge
    # Positively charged: R, H, K
    # Negatively charged: D, E
    # Neutral: A, C, F, G, I, L, M, N, P, Q, S, T, V, W, Y
    # As noted by IMGT: The sulfhydryl group of cystein and
    # phenolic hydroxyl group of tyrosine show some degree of pH-dependent ionization.

    NEUTRAL = 0
    POSITIVE = 1
    NEGATIVE = -1

    RESIDUE_CHARGE = {
        "ALA": 0,
        "ARG": 1,
        "ASN": 0,
        "ASP": -1,
        "CYS": 0,
        "GLN": 0,
        "GLU": -1,
        "GLY": 0,
        "HIS": 1,  # Note: The charge of histidine is pH-dependent. A value of 1 is used here, but 0.1 is also common.
        "ILE": 0,
        "LEU": 0,
        "LYS": 1,
        "MET": 0,
        "PHE": 0,
        "PRO": 0,
        "SER": 0,
        "THR": 0,
        "TRP": 0,
        "TYR": 0,
        "VAL": 0,

        "A": 0,
        "R": 1,
        "N": 0,
        "D": -1,
        "C": 0,
        "Q": 0,
        "E": -1,
        "G": 0,
        "H": 1,
        "I": 0,
        "L": 0,
        "K": 1,
        "M": 0,
        "F": 0,
        "P": 0,
        "S": 0,
        "T": 0,
        "W": 0,
        "Y": 0
    }

    @staticmethod
    def charge_of_residue(residue: Residue,
                          assign_attribute: bool = False,
                          key: str = "charge") -> int:
        """
        Returns the charge of a residue.

        Args:
            residue (Residue): The residue for which to determine the charge.
            assign_attribute (bool): Whether to assign the charge to the residue.
            key (str): The key to use for the attribute.

        Returns:
            int: The charge of the residue.
        """
        charge = Charge.RESIDUE_CHARGE.get(residue.residue_type, Charge.NEUTRAL)
        if assign_attribute:
            residue.set_attribute(key, charge)
        return charge

    @staticmethod
    def charge_of_chain(chain: Chain,
                        assign_attribute: bool = False,
                        key: str = "charge") -> int:
        """
        Returns the charge of a chain.

        Args:
            chain (Chain): The chain for which to determine the charge.
            assign_attribute (bool): Whether to assign the charge to the chain.
            key (str): The key to use for the attribute.

        Returns:
            int: The charge of the chain.
        """
        charge = sum([Charge.charge_of_residue(residue, assign_attribute=assign_attribute, key=key) for residue in chain.residues])
        if assign_attribute:
            chain.set_attribute(key, charge)
        return charge

    @staticmethod
    def charge_of_protein(protein: Protein,
                          assign_attribute: bool = False,
                          key: str = "charge") -> int:
        """
        Returns the charge of a protein.

        Args:
            protein (Protein): The protein for which to determine the charge.
            assign_attribute (bool): Whether to assign the charge to the protein.
            key (str): The key to use for the attribute.

        Returns:
            int: The charge of the protein.
        """
        charge = sum([Charge.charge_of_chain(chain, assign_attribute=assign_attribute, key=key) for chain in protein.chains])
        if assign_attribute:
            protein.set_attribute(key, charge)
        return charge

    @staticmethod
    def charge_of_sequence(sequence: Sequence,
                           assign_attribute: bool = False,
                           key: str = "charge") -> int:
        """
        Returns the charge of a sequence.

        Args:
            sequence (Sequence): The sequence for which to determine the charge.
            assign_attribute (bool): Whether to assign the charge to the sequence.
            key (str): The key to use for the attribute.

        Returns:
            int: The charge of the sequence.
        """
        charge = sum([Charge.RESIDUE_CHARGE.get(residue, Charge.NEUTRAL) for residue in sequence])
        if assign_attribute:
            sequence.set_attribute(key, charge)
        return charge
