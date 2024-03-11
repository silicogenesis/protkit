#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `FastaIO` to read and write
FASTA files. FASTA files contain one or more sequences
(protein or nucleotide) of biological data with their
associated metadata.

See https://en.wikipedia.org/wiki/FASTA_format for more information.
See https://www.rcsb.org/ for examples of FASTA files.

Methods are static and can be called without instantiating the class.
The main functions exposed by the class are:

- `load()` to load a protein from a Fasta file.
- `save()` to save a protein to a Fasta file.
"""

from typing import List
from protkit.seq.sequence import Sequence


class FastaIO:
    @staticmethod
    def load(file_path: str) -> List[Sequence]:
        """
        Loads a FASTA file and returns a list of sequences.

        Args:
            file_path (str): The path to the FASTA file.

        Returns:
            List[Sequence]: A list of sequences.
        """

        sequences = []
        with open(file_path, "rt") as file:
            description = None
            seq = ""

            for line in file:
                if line.startswith(">") or line.startswith(";"):
                    if seq != "":
                        # The previous sequence is added.
                        # There were no empty lines between the
                        # previous sequence and the current sequence.
                        sequences.append(Sequence(seq, description))
                        seq = ""
                        description = None

                    # Comments/description lines start in ";" or ">".
                    if description is None:
                        description = line[1:].strip()
                    else:
                        description += line[1:].strip()
                elif line.strip() == "":
                    # Empty lines are ignored.
                    # Empty lines are also used to separate sequences.
                    if description is not None or seq != "":
                        sequences.append(Sequence(seq, description))
                        seq = ""
                        description = None
                else:
                    # Sequence lines start with a letter.
                    # Sequence lines can end with a "*" character.
                    seq += line.strip()
                    if seq[-1] == "*":
                        seq = seq[:-1]

            # The last sequence is added assuming there is no empty line at the end.
            if description is not None or seq != "":
                sequences.append(Sequence(seq, description))

        return sequences

    @staticmethod
    def save(sequence: [Sequence, List[Sequence]], file_path: str, line_length: [int, None] = 80) -> None:
        """
        Saves a FASTA file.

        Args:
            sequence (Union[Sequence, List[Sequence]]): The sequence(s) to save.
            file_path (str): The path to the FASTA file.
            line_length (int): The length of the lines in the FASTA file.
                If None, the sequence will be saved on one line.

        Returns:
            None
        """
        if type(sequence) is Sequence:
            sequence = [sequence]

        with open(file_path, "w") as file:
            for seq in sequence:
                if seq.description is not None:
                    file.write(">" + seq.description + "\n")
                else:
                    file.write(">\n")

                if line_length is None:
                    file.write(seq.sequence + "\n")
                else:
                    for i in range((len(seq.sequence) - 1) // line_length + 1):
                        file.write(seq.to_string(i * line_length, (i + 1) * line_length - 1) + "\n")
                file.write("\n")
