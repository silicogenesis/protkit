#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `ProteinSequence` to represent a protein sequence.
"""

from typing import Optional, List, Union
from protkit.seq.sequence import Sequence


class ProteinSequence(Sequence):
    THREE_TO_ONE = {
        # Standard Amino Acids
        "ALA": "A",
        "CYS": "C",
        "ASP": "D",
        "GLU": "E",
        "PHE": "F",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LYS": "K",
        "LEU": "L",
        "MET": "M",
        "ASN": "N",
        "PRO": "P",
        "GLN": "Q",
        "ARG": "R",
        "SER": "S",
        "THR": "T",
        "VAL": "V",
        "TRP": "W",
        "TYR": "Y",

        # Non-standard Amino Acids
        "SEC": "U",     # Selenocysteine
        "PYL": "O",     # Pyrrolysine
        "ASX": "B",     # Asparagine or Aspartic Acid
        "GLX": "Z",     # Glutamine or Glutamic Acid
        "XLE": "J",     # Leucine or Isoleucine
        "UNK": "X",     # Unknown
        "XAA": "X"
    }

    ONE_TO_THREE = {
        # Standard Amino Acids
        "A": "ALA",
        "C": "CYS",
        "D": "ASP",
        "E": "GLU",
        "F": "PHE",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "K": "LYS",
        "L": "LEU",
        "M": "MET",
        "N": "ASN",
        "P": "PRO",
        "Q": "GLN",
        "R": "ARG",
        "S": "SER",
        "T": "THR",
        "V": "VAL",
        "W": "TRP",
        "Y": "TYR",

        # Non-standard Amino Acids
        "U": "SEC",
        "O": "PYL",
        "B": "ASX",
        "Z": "GLX",
        "J": "XLE",
        "X": "UNK"
    }

    def __init__(self,
                 sequence: Union[str, List[str]],
                 description: Optional[str] = None,
                 chain_id: Optional[str] = None
                 ):
        """
        Constructor.

        Args:
            sequence (Union[str, List[str]]): The sequence of residue names.
            description (Optional[str]): Optional description associated with the sequence.
            chain_id (Optional[str]): Optional chain ID associated with the sequence.

        Returns:
            None
        """
        super().__init__(sequence, description, chain_id)

    @staticmethod
    def from_sequence(seq: Sequence):
        """
        Creates a ProteinSequence object from a Sequence object.

        Args:
            seq (Sequence): The Sequence object.

        Returns:
            ProteinSequence: The ProteinSequence object.
        """
        sequence = ProteinSequence(seq.sequence, seq.description, seq.chain_id)
        return sequence

    def to_single_letter(self):
        """
        Converts the sequence to single-letter residue names.

        Args:
            None

        Returns:
            None
        """

        for i in range(len(self._sequence)):
            res = self._sequence[i].upper()
            if res in ProteinSequence.THREE_TO_ONE:
                self._sequence[i] = ProteinSequence.THREE_TO_ONE[res]
            else:
                self._sequence[i] = "X"

    def to_triple_letter(self):
        """
        Converts the sequence to three-letter residue names.

        Args:
            None

        Returns:
            None
        """

        for i in range(len(self._sequence)):
            res = self._sequence[i].upper()
            if res in ProteinSequence.ONE_TO_THREE:
                self._sequence[i] = ProteinSequence.ONE_TO_THREE[res]
            else:
                self._sequence[i] = "UNK"
