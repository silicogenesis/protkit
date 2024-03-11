#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `Sequence` to represent a sequence of residues.

This class provides a representation of a sequence. The class
provides core functionality for sequence alignment and comparison,
that is applicable to both nucleotide and protein sequences.

Specific functionality for nucleotide and protein sequences are
provided by the NucleotideSequence and ProteinSequence classes.

This class should be able to handle sequences of any length, and
of any type (nucleotide or protein). The underlying representation
of the sequence should be able to handle any type of sequence.
"""

from typing import Optional, List, Union
from protkit.core.extend_attributes import ExtendedAttributes


class Sequence(ExtendedAttributes):
    def __init__(self,
                 sequence: Union[str, List[str]],
                 description: Optional[str] = None,
                 chain_id: Optional[str] = None):
        """
        Constructor.

        Args:
            sequence (Union[str, List[str]]): The sequence of residue names.
            description (Optional[str]): Optional description associated with the sequence.
            chain_id (Optional[str]): Optional chain ID associated with the sequence.

        Returns:
            None

        Notes:
            The sequence can be provided as a string or a list of strings.
            If provided as a string, the string will be converted to a list of strings.

            A sequence can be respresented by a string of single letters, such as "AGILE",
            or a list of three-letter codes such as ["ALA", "GLY", "ILE", "LEU", "GLU"]. For
            consistency, the sequence is always represented as a list of codes.
        """
        super().__init__()

        if type(sequence) is str:
            # Convert string representations to a list representation.
            sequence = list(sequence)

        self._sequence = sequence
        self._chain_id = chain_id
        self._description = description

    def __str__(self):
        """
        Returns a string representation of the sequence.

        Returns:
            str: The string representation of the sequence.
        """
        if self._sequence is None:
            return ""
        if len(self.sequence) == 0:
            return ""
        if len(self._sequence[0]) == 1:
            spacer = ""
        else:
            spacer = " "
        return spacer.join(self._sequence)

    def to_string(self, start_index: int, end_index: int) -> str:
        """
        Returns a string representation of a subsequence.

        Args:
            start_index (int): The start index of the subsequence.
            end_index (int): The end index of the subsequence.

        Returns:
            str: The string representation of the subsequence.
        """
        if self._sequence is None or len(self.sequence) == 0:
            return ""

        sub_sequence = [self._sequence[i] for i in range(start_index, min(end_index + 1, self.length))]

        if len(self._sequence[0]) == 1:
            return "".join(sub_sequence)
        else:
            return " ".join(sub_sequence)


    def __len__(self):
        """
        Returns the length of the sequence.

        Returns:
            int: The length of the sequence.
        """
        return len(self._sequence)

    def __getitem__(self, index):
        """
        Returns the residue at the specified index.

        Args:
            index (int): The index of the residue to return.

        Returns:
            str: The residue at the specified index.
        """
        return self._sequence[index]

    def __setitem__(self, index, value):
        """
        Sets the residue at the specified index.

        Args:
            index (int): The index of the residue to set.
            value (str): The residue to set.

        Returns:
            None
        """
        self._sequence[index] = value

    @property
    def sequence(self):
        """
        Returns the sequence.

        Returns:
            str: The sequence.
        """
        return self._sequence

    @sequence.setter
    def sequence(self, sequence: Union[str, List[str]]):
        """
        Sets the sequence.

        Args:
            sequence (Union[str, List[str]]): The sequence.

        Returns:
            None
        """
        if type(sequence) is str:
            # Convert string representations to a list representation.
            self._sequence = list(sequence)
        else:
            self._sequence = sequence

    @property
    def description(self):
        """
        Returns the description.

        Returns:
            str: The description.
        """
        return self._description

    @description.setter
    def description(self, description):
        """
        Sets the description.

        Args:
            description (str): The description.

        Returns:
            None
        """
        self._description = description

    @property
    def chain_id(self):
        """
        Returns the chain ID.

        Returns:
            str: The chain ID.
        """
        return self._chain_id

    @chain_id.setter
    def chain_id(self, chain_id):
        """
        Sets the chain ID.

        Args:
            chain_id (str): The chain ID.

        Returns:
            None
        """
        self._chain_id = chain_id

    @property
    def length(self):
        """
        Returns the length of the sequence.

        Returns:
            int: The length of the sequence.
        """
        return len(self._sequence)
