#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `StructureEval` for evaluating protein structures.

Various scores, such as RMSD, TMScore, GDT TS and Fnat, can be calculated to evaluate
the quality of a protein structure. These scores can be used to compare two structures
or to evaluate a single structure against a reference structure.
"""

from typing import List


class StructureEval:
    @staticmethod
    def RMSD(coordinates: List[List[float]], reference: List[List[float]]):
        """
        Calculates the root-mean-square deviation (RMSD) between two sets of coordinates.

        It provides a measure of the average distance between atoms of two superimposed structures.
        See https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions.

        Args:
            coordinates (List[List[float]]): A list of lists containing the coordinates to be evaluated.
            reference (List[List[float]]): A list of lists containing the reference coordinates.
        Returns:
            float: A float representing the RMSD between the two sets of coordinates.
        """
        rmsd = 0.0
        for i in range(len(coordinates)):
            rmsd += (coordinates[i][0] - reference[i][0])**2 + (coordinates[i][1] - reference[i][1])**2 + (coordinates[i][2] - reference[i][2])**2
        return (rmsd / len(coordinates))**0.5

    @staticmethod
    def TMscore(coordinates: List[List[float]], reference: List[List[float]]):
        """
        Calculates the template modelling (TM-score) between two sets of coordinates.

        The TM-score is a measure of the structural similarity between two protein structures.
        It is intended to be more robust to variations in protein length than the RMSD. It
        provides a score in the range (0, 1] where 1 indicates a perfect match. See
        https://en.wikipedia.org/wiki/Template_modeling_score

        Args:
            coordinates (List[List[float]]): A list of lists containing the coordinates to be evaluated.
            reference (List[List[float]]): A list of lists containing the reference coordinates.
        Returns:
            float: A float representing the TM-score between the two sets of coordinates.
        """
        raise NotImplementedError("TMscore() not implemented yet.")

    @staticmethod
    def GDT_TS(coordinates: List[List[float]], reference: List[List[float]]):
        """
        Calculates the GDT_TS score between two sets of coordinates.

        Args:
            coordinates (List[List[float]]): A list of lists containing the coordinates to be evaluated.
            reference (List[List[float]]): A list of lists containing the reference coordinates.
        Returns:
            float: A float representing the GDT_TS score between the two sets of coordinates.
        """
        raise NotImplementedError("GDS_TS() not implemented yet.")

    @staticmethod
    def Fnat():
        """
        Calculates the fraction of native contacts (Fnat).
        """
        raise NotImplementedError("Fnat() score not implemented yet.")