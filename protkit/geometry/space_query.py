#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `SpaceQuery` to do fast queries of atomic coordinates.

It relies on the KDTree implementation provided by scikit-learn. It
abstracts computations of distances and nearest neighbours using
numpy arrays to hide these implementation details from the user.
"""

from typing import List, Tuple
from sklearn.neighbors import KDTree
import numpy as np


class SpaceQuery:
    def __init__(self, coordinates: List[Tuple[float, float, float]], leaf_size: int = 50):
        """
        Constructor.  Constructs a KDTree from the atomic coordinates.

        Args:
            coordinates (List[Tuple[float, float, float]]): List of atomic coordinates.
                The coordinates should be in the form of a list of lists,
                where each sublist contains the x, y and z coordinates of
                an atom.
            leaf_size (int): Leaf size of the KDTree. Default is 50.

        Returns:
            None
        """
        self._coordinates = np.array(coordinates)
        self._tree = KDTree(self._coordinates, leaf_size=leaf_size)

    @property
    def coordinates(self):
        """
        Getter for the coordinates.

        Returns:
            np.ndarray: The coordinates.
        """
        return self._coordinates

    @property
    def tree(self):
        """
        Getter for the KDTree.

        Returns:
            KDTree: The KDTree.
        """
        return self._tree

    def query_nearest(self, coordinates: List[List[float]], k: int = 1) -> Tuple[List[int], List[float]]:
        """
        Query the KDTree for the k nearest neighbours of a given coordinates.

        Args:
            coordinates (List[List[float]]): List of coordinates to query.
            k (int): Number of nearest neighbours to return. Default is 1.

        Returns:
            Tuple[List[int], List[float]]: List of indices of the k nearest neighbours and their distances.
        """
        coordinates = np.array(coordinates)
        distances, indices = self._tree.query(coordinates, k=k)
        return indices, distances

    def query_partners(self, coordinates: List[Tuple[float, float, float]], radius: float) -> (List[Tuple], List[Tuple]):
        """
        Query the KDTree for all atoms within a certain radius of a given coordinates.

        Args:
            coordinates (List[Tuple[float, float, float]]): List of coordinates to query.
            radius (float): Radius within which to search for atoms.

        Returns:

        """
        coordinates = np.array(coordinates)
        indices = self._tree.query_radius(coordinates, r=radius)
        indices1 = set()
        for (i, index) in enumerate(indices):
            if len(index) > 0:
                for j in index:
                    indices1.add(j)
        indices1 = sorted(list(indices1))
        indices2 = [i for i in range(len(coordinates)) if len(indices[i]) > 0]
        return indices1, indices2

    def query_indices(self, coordinates: List[List[float]], radius: float) -> List[List[int]]:
        """
        Query the KDTree for all atoms within a certain radius of a given coordinates.

        Args:
            coordinates (List[List[float]]): List of coordinates to query.
            radius (float): Radius within which to search for atoms.

        Returns:
            List[List[int]]: List of lists of indices of atoms within the radius of the query coordinates.
        """
        coordinates = np.array(coordinates)
        indices = self._tree.query_radius(coordinates, r=radius)
        return indices

    def query_distance(self, coordinates: List[Tuple[float, float, float]], radius: float) -> Tuple[List[List[int]], List[List[float]]]:
        """
        Query the KDTree for all atoms within a certain radius of a given coordinates.

        Args:
            coordinates (List[List[float]]): List of coordinates to query.
            radius (float): Radius within which to search for atoms.

        Returns:
            List[List[float]]: List of lists of distances of atoms within the radius of the query coordinates.
        """
        coordinates = np.array(coordinates)
        indices, distances = self._tree.query_radius(coordinates, r=radius, return_distance=True)
        return indices, distances

    def query_distance_np(self, coordinates: List[List[float]], radius: float) -> Tuple[List[List[int]], List[List[float]]]:
        """
        Query the KDTree for all atoms within a certain radius of a given coordinates.

        Args:
            coordinates (List[List[float]]): List of coordinates to query.
            radius (float): Radius within which to search for atoms.

        Returns:
            List[List[float]]: List of lists of distances of atoms within the radius of the query coordinates.
        """
        coordinates = np.array(coordinates)
        indices, distances = self._tree.query_radius(coordinates, r=radius, return_distance=True)
        return indices, distances

