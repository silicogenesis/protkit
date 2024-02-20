#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `SequenceAlignment` to align two sequences.

- `local_align()` to align two sequences using the Smith-Waterman algorithm.
- `global_align()` to align two sequences using the Needleman-Wunsch algorithm.
"""

from typing import List
from protkit.seq.sequence import Sequence


class SequenceAlignment:
    @staticmethod
    def local_align(
            sequence1: Sequence,
            sequence2: Sequence,
            match_value: int = 2,
            mismatch_penalty: int = -1,
            gap_penalty: int = -1) -> (List[str], List[str]):
        """
        Aligns two sequences using the Smith-Waterman algorithm.

        Args:
            sequence1 (Sequence): The first sequence.
            sequence2 (Sequence): The second sequence.
            match_value (int): The score for a match.
            mismatch_penalty (int): The penalty for a mismatch.
            gap_penalty (int): The penalty for a gap.

        Returns:
            List[str], List[str]: The aligned sequences.

        Raises:
            None
        """
        sequence1 = sequence1.sequence
        sequence2 = sequence2.sequence

        # Initialize the matrix.
        rows = len(sequence1) + 1
        cols = len(sequence2) + 1
        matrix = [[0] * cols for _ in range(rows)]

        # Fill the matrix.
        for i in range(1, rows):
            for j in range(1, cols):
                match_score = matrix[i - 1][j - 1] + (match_value if sequence1[i - 1] == sequence2[j - 1] else mismatch_penalty)
                delete_score = matrix[i - 1][j] + gap_penalty
                insert_score = matrix[i][j - 1] + gap_penalty
                matrix[i][j] = max(0, match_score, delete_score, insert_score)

        # Find the maximum score.
        max_i = 0
        max_j = 0
        max_score = -1
        alignment1 = []
        alignment2 = []
        for i in range(1, rows):
            for j in range(1, cols):
                if matrix[i][j] > max_score:
                    max_score = matrix[i][j]
                    max_i = i
                    max_j = j

        # Traceback
        i = max_i
        j = max_j
        while i > 0 and j > 0 and matrix[i][j] > 0:
            if matrix[i][j] == matrix[i - 1][j - 1] + (match_value if sequence1[i - 1] == sequence2[j - 1] else mismatch_penalty):
                alignment1.append(sequence1[i - 1])
                alignment2.append(sequence2[j - 1])
                i -= 1
                j -= 1
            elif matrix[i][j] == matrix[i - 1][j] + gap_penalty:
                alignment1.append(sequence1[i - 1])
                alignment2.append("-")
                i -= 1
            else:
                alignment1.append("-")
                alignment2.append(sequence2[j - 1])
                j -= 1

        alignment1.reverse()
        alignment2.reverse()

        return alignment1, alignment2

    @staticmethod
    def global_align(
            sequence1: Sequence,
            sequence2: Sequence,
            match_value: int = 2,
            mismatch_penalty: int = -1,
            gap_penalty: int = -1
    ) -> (List[str], List[str]):
        """
        Aligns two sequences using the Needleman-Wunsch algorithm.

        Args:
            sequence1 (Sequence): The first sequence.
            sequence2 (Sequence): The second sequence.
            match_value (int): The score for a match.
            mismatch_penalty (int): The penalty for a mismatch.
            gap_penalty (int): The penalty for a gap.
        Returns:
            List[str], List[str]: The aligned sequences.

        Raises:
            None
        """
        sequence1 = sequence1.sequence
        sequence2 = sequence2.sequence

        # Initialize the scoring matrix.
        rows = len(sequence1) + 1
        cols = len(sequence2) + 1
        matrix = [[0] * cols for _ in range(rows)]

        # Assign gap penalties.
        for i in range(1, rows):
            matrix[i][0] = i * gap_penalty
        for j in range(1, cols):
            matrix[0][j] = j * gap_penalty

        # Fill the scoring matrix.
        for i in range(1, rows):
            for j in range(1, cols):
                match_score = matrix[i - 1][j - 1] + (match_value if sequence1[i - 1] == sequence2[j - 1] else mismatch_penalty)
                delete_score = matrix[i - 1][j] + gap_penalty
                insert_score = matrix[i][j - 1] + gap_penalty
                matrix[i][j] = max(match_score, delete_score, insert_score)

        # Traceback.
        i = rows - 1
        j = cols - 1
        alignment1 = []
        alignment2 = []

        while i > 0 or j > 0:
            if i > 0 and j > 0 and matrix[i][j] == matrix[i - 1][j - 1] + (match_value if sequence1[i - 1] == sequence2[j - 1] else mismatch_penalty):
                alignment1.append(sequence1[i - 1])
                alignment2.append(sequence2[j - 1])
                i -= 1
                j -= 1
            elif i > 0 and matrix[i][j] == matrix[i - 1][j] + gap_penalty:
                alignment1.append(sequence1[i - 1])
                alignment2.append("-")
                i -= 1
            else:
                alignment1.append("-")
                alignment2.append(sequence2[j - 1])
                j -= 1

        alignment1.reverse()
        alignment2.reverse()

        return alignment1, alignment2
