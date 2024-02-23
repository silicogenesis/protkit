#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `SequenceEval` for evaluating protein sequences.

Various metrics, such as sequence identity, similarity, and coverage, can be calculated
to evaluate the quality of a protein sequence. These scores can be used to compare two
sequences or to evaluate a single sequence against a reference sequence.
"""

from protkit.seq import Sequence
from protkit.metrics.scoring_matrix import ScoringMatrix


class SequenceEval:
    @staticmethod
    def sequence_identity(seq1: Sequence,
                          seq2: Sequence) -> float:
        """
        Calculate the sequence identity between two sequences.

        Sequence identity is a measure of the similarity between two sequences. It is
        defined as the number of identical residues divided by the total number of residues.

        Args:
            seq1 (Sequence): The first sequence.
            seq2 (Sequence): The second sequence.

        Returns:
            float: A float representing the sequence identity between the two sequences.

        Raises:
            ValueError: If the sequences are not of equal length.
        """
        if seq1.length != seq2.length:
            raise ValueError("Sequences must be of equal length.")
        if seq1.length == 0:
            raise ValueError("Sequences must have a length greater than 0.")

        identity = 0
        for i in range(seq1.length):
            if seq1[i] == seq2[i]:
                identity += 1

        return identity / seq1.length

    @staticmethod
    def sequence_similarity(seq1: Sequence,
                            seq2: Sequence,
                            match_score: int = 2,
                            mismatch_score: int = -1,
                            scoring_matrix: ScoringMatrix = None) -> float:
        """
        Calculate the sequence similarity between two sequences.

        Sequence similarity is a measure of the similarity between two sequences. It is
        defined as the sum of the scores for matching residues divided by the total number
        of residues.

        Args:
            seq1 (str): The first sequence.
            seq2 (str): The second sequence.
            match_score (int): The score to assign to matching residues.
            mismatch_score (int): The score to assign to mismatching residues.
            scoring_matrix (ScoringMatrix): The scoring matrix to use for scoring residue pairs.

        Returns:
            float: A float representing the sequence similarity between the two sequences.

        Raises:
            ValueError: If the sequences are not of equal length.
        """
        if seq1.length != seq2.length:
            raise ValueError("Sequences must be of equal length.")
        if seq1.length == 0:
            raise ValueError("Sequences must have a length greater than 0.")

        similarity = 0
        if scoring_matrix is not None:
            for i in range(seq1.length):
                similarity += scoring_matrix.score(seq1[i], seq2[i])
        else:
            for i in range(seq1.length):
                if seq1[i] == seq2[i]:
                    similarity += match_score
                else:
                    similarity += mismatch_score

        return similarity / seq1.length

    @staticmethod
    def alignment_coverage(seq: Sequence, gap_symbol: str = "-") -> float:
        """
        Calculate the alignment coverage between two sequences.

        Alignment coverage is a measure of the proportion of residues in one sequence that
        are aligned with residues in another sequence. It is defined as the number of aligned
        residues divided by the total number of residues.

        Args:
            seq (Sequence): The sequence.
            gap_symbol (str): The symbol used to represent gaps in the alignment.

        Returns:
            float: A float representing the alignment coverage of the sequence.
        """
        if seq.length == 0:
            raise ValueError("Sequence must have a length greater than 0.")

        aligned_residues = 0
        for residue in seq:
            if residue != gap_symbol:
                aligned_residues += 1

        return aligned_residues / seq.length

    @staticmethod
    def edit_distance(seq1: Sequence, seq2: Sequence) -> float:
        """
        Calculate the edit distance (Levenshtein distance) between two sequences.

        The edit distance is a measure of the similarity between two sequences. It is defined
        as the minimum number of single-character edits (insertions, deletions, or substitutions)
        required to change one sequence into the other.

        Args:
            seq1 (str): The first sequence.
            seq2 (str): The second sequence.

        Returns:
            float: A float representing the edit distance between the two sequences.

        Raises:
            ValueError: If the sequences are not of equal length.
        """
        if seq1.length == 0 or seq2.length == 0:
            raise ValueError("Sequences must have a length greater than 0.")

        # Initialize the matrix.
        rows = seq1.length + 1
        cols = seq2.length + 1
        matrix = [[0] * cols for _ in range(rows)]

        # Fill the matrix.
        for i in range(rows):
            matrix[i][0] = i
        for j in range(cols):
            matrix[0][j] = j

        for i in range(1, rows):
            for j in range(1, cols):
                if seq1[i - 1] == seq2[j - 1]:
                    cost = 0
                else:
                    cost = 1
                matrix[i][j] = min(matrix[i - 1][j] + 1,    # Deletion
                                   matrix[i][j - 1] + 1,    # Insertion
                                   matrix[i - 1][j - 1] + cost)  # Substitution

        # The edit distance is the value in the bottom-right corner of the matrix.
        edit_distance = matrix[rows - 1][cols - 1]

        return edit_distance
