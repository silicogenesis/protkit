#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `CircularVariance` to calculate the circular variance of a protein.

Circular variance is a measure of uniformity of a set of points on a sphere. It is
calculated as 1 - |sum(unit_vectors)| / n, where n is the number of points and
unit_vectors are the vectors from the center to the points.

The circular variance of a protein can be calculated for each atom or residue.
For atoms, the circular variance is calculated using the coordinates of the atom.
For residues, the circular variance is calculated using the coordinates of the alpha carbon
atom of the residue.
"""

from typing import List
import numpy as np

from protkit.structure.protein import Protein
from protkit.geometry.space_query import SpaceQuery


class CircularVariance:
    @staticmethod
    def circular_variance(coordinates: List[List[float]], radius: float = 5.0) -> List[float]:
        # Set up the space query.
        space_query = SpaceQuery(coordinates)

        circular_variance = []

        coordinates = space_query.coordinates
        for i, coordinate in enumerate(coordinates):
            # Determine distances and neighbours to all the neighbouring atoms.
            neighbours, distances = space_query.query_distance_np([coordinate], radius)
            neighbours = neighbours[0]
            distances = distances[0]

            # Calculate circular variance.
            non_zero_indices = np.where(distances > 0.01)
            neighbour_coordinates = coordinates[neighbours[non_zero_indices]]
            coordinate = np.array([coordinate])
            unit_vectors = neighbour_coordinates - coordinate
            unit_vectors /= distances[non_zero_indices].reshape(-1, 1)
            cv = 1.0 - np.linalg.norm(np.sum(unit_vectors, axis=0)) / neighbours[non_zero_indices].shape[0]

            circular_variance.append(float(cv))

        return circular_variance

    @staticmethod
    def circular_variance_by_atom(protein: Protein,
                                  radius: float = 5.0,
                                  assign_attribute: bool = False,
                                  key: str = "cv_atom") -> List[float]:
        """
        Calculate the circular variance of a protein for each atom.

        Args:
            protein (Protein): The protein for which to calculate the circular variance.
            radius (float): The radius to use for the calculation.
            assign_attribute (bool): Whether to assign the circular variance to the atoms.
            key (str): The key to use for the attribute.

        Returns:
            List[float]: The circular variance of each atom in the protein.
        """
        # Prepare coordinates
        atom_coordinates = [(atom.x, atom.y, atom.z) for atom in protein.atoms]

        # Calculate circular variance
        circular_variance_by_atom = CircularVariance.circular_variance(atom_coordinates, radius)

        # Assign circular variance to atoms
        if assign_attribute:
            for atom, cv in zip(protein.atoms, circular_variance_by_atom):
                atom.set_attribute(key, cv)

        return circular_variance_by_atom

    @staticmethod
    def circular_variance_by_residue(protein: Protein,
                                     radius: float = 5.0,
                                     assign_attribute: bool = False,
                                     key: str = "cv_residue") -> List[float]:
        """
        Calculate the circular variance of a protein for each residue.

        Args:
            protein (Protein): The protein for which to calculate the circular variance.
            radius (float): The radius to use for the calculation.
            assign_attribute (bool): Whether to assign the circular variance to the residues.
            key (str): The key to use for the attribute.

        Returns:
            List[float]: The circular variance of each residue in the protein.
        """
        # Prepare coordinates
        residues = protein.filter_residues()
        residue_coordinates = [(residue.get_atom("CA").x, residue.get_atom("CA").y, residue.get_atom("CA").z) for residue in residues]

        # Assign circular variance to residues
        circular_variance_by_residue = CircularVariance.circular_variance(residue_coordinates, radius)
        if assign_attribute:
            residues = protein.filter_residues()
            for residue, cv in zip(residues, circular_variance_by_residue):
                residue.set_attribute(key, cv)

        return circular_variance_by_residue
