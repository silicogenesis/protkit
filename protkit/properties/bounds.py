#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `Bounds` to calculate the bounds and center of a residue, chain or protein.

Bounds are defined as the minimum and maximum x, y and z coordinates of the atoms in a residue, chain or protein.
Center is defined as the center of the bounds.
"""

from typing import Tuple, List
from protkit.structure.protein import Protein
from protkit.structure.chain import Chain
from protkit.structure.residue import Residue


class Bounds:
    @staticmethod
    def _bounds(xyz: List[Tuple[float, float, float]]):
        """
        Returns the bounds of a set of coordinates.

        Args:
            xyz (Tuple[float, float, float]): The coordinates.

        Returns:
            Tuple[float, float, float, float, float, float]: The bounds of the coordinates.
        """
        min_x = min_y = min_z = float("inf")
        max_x = max_y = max_z = float("-inf")

        for x, y, z in xyz:
            min_x = min(min_x, x)
            min_y = min(min_y, y)
            min_z = min(min_z, z)
            max_x = max(max_x, x)
            max_y = max(max_y, y)
            max_z = max(max_z, z)

        return min_x, min_y, min_z, max_x, max_y, max_z

    @staticmethod
    def _center(xyz: List[Tuple[float, float, float]]):
        """
        Returns the center of a set of coordinates.

        Args:
            xyz (Tuple[float, float, float]): The coordinates.

        Returns:
            Tuple[float, float, float]: The center of the coordinates.
        """
        min_x, min_y, min_z, max_x, max_y, max_z = Bounds._bounds(xyz)
        return (min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2

    @staticmethod
    def _bounds_and_center(xyz: List[Tuple[float, float, float]]):
        """
        Returns the bounds and center of a set of coordinates.

        Args:
            xyz (Tuple[float, float, float]): The coordinates.

        Returns:
            Tuple[float, float, float, float, float, float, float, float, float]: The bounds and center of the coordinates.
        """
        min_x, min_y, min_z, max_x, max_y, max_z = Bounds._bounds(xyz)
        center_x, center_y, center_z = (min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2
        return (min_x, min_y, min_z, max_x, max_y, max_z), (center_x, center_y, center_z)

    @staticmethod
    def _bounds_and_center_of_residue(residue: Residue,
                                      assign_bounds: bool = False,
                                      bounds_key: str = "bounds",
                                      assign_center: bool = False,
                                      center_key: str = "center") -> Tuple[Tuple[float, float, float, float, float, float], Tuple[float, float, float]]:
        """
        Returns the bounds and center of a residue.

        Args:
            residue (Residue): The residue.

        Returns:
            Tuple[Tuple[float, float, float, float, float, float], Tuple[float, float, float]]: The bounds and center of the residue.
        """
        xyz = [(atom.x, atom.y, atom.z) for atom in residue.atoms]
        bounds, center = Bounds._bounds_and_center(xyz)
        if assign_bounds:
            residue.set_attribute(bounds_key, bounds)
        if assign_center:
            residue.set_attribute(center_key, center)

        return bounds, center

    @staticmethod
    def bounds_of_residue(residue: Residue,
                          assign_attribute: bool = False,
                          key: str = "bounds") -> Tuple[float, float, float, float, float, float]:
        """
        Returns the bounds of a residue.

        Args:
            residue (Residue): The residue.
            assign_attribute (bool): Whether to assign the bounds of the residue as an attribute.
            key (str): The key to use for the attribute.

        Returns:
            Tuple[float, float, float, float, float, float]: The bounds of the residue.

        Raises:
            None
        """
        bounds, center = Bounds._bounds_and_center_of_residue(residue, assign_bounds=assign_attribute, bounds_key=key)
        return bounds

    @staticmethod
    def center_of_residue(residue: Residue,
                          assign_attribute: bool = False,
                          key: str = "center") -> Tuple[float, float, float]:
        """
        Returns the center of a residue.

        Args:
            residue (Residue): The residue.
            assign_attribute (bool): Whether to assign the center of the residue as an attribute.
            key (str): The key to use for the attribute.

        Returns:
            Tuple[float, float, float]: The center of the residue.

        Raises:
            None
        """
        bounds, center = Bounds._bounds_and_center_of_residue(residue, assign_center=assign_attribute, center_key=key)
        return center

    @staticmethod
    def _bounds_and_center_of_chain(chain: Chain,
                                    assign_bounds: bool = False,
                                    bounds_key: str = "bounds",
                                    assign_center: bool = False,
                                    center_key: str = "center") -> Tuple[Tuple[float, float, float, float, float, float], Tuple[float, float, float]]:
        """
        Returns the bounds and center of a chain.

        Args:
            chain (Chain): The chain.

        Returns:
            Tuple[float, float, float, float, float, float, float, float, float]: The bounds and center of the chain.
        """
        min_x = min_y = min_z = float("inf")
        max_x = max_y = max_z = float("-inf")

        for residue in chain.residues:
            bounds, center = Bounds._bounds_and_center_of_residue(residue, assign_bounds=assign_bounds, bounds_key=bounds_key, assign_center=assign_center, center_key=center_key)
            min_x = min(min_x, bounds[0])
            min_y = min(min_y, bounds[1])
            min_z = min(min_z, bounds[2])
            max_x = max(max_x, bounds[3])
            max_y = max(max_y, bounds[4])
            max_z = max(max_z, bounds[5])

        cx, cy, cz = (min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2

        if assign_bounds:
            chain.set_attribute(bounds_key, (min_x, min_y, min_z, max_x, max_y, max_z))
        if assign_center:
            chain.set_attribute(center_key, (cx, cy, cz))

        return (min_x, min_y, min_z, max_x, max_y, max_z), (cx, cy, cz)

    @staticmethod
    def bounds_of_chain(chain: Chain,
                        assign_attribute: bool = False,
                        key: str = "bounds") -> Tuple[float, float, float, float, float, float]:
        """
        Returns the bounds of a chain.

        Args:
            chain (Chain): The chain.

        Returns:
            Tuple[float, float, float, float, float, float]: The bounds of the chain.

        Raises:
            None
        """
        bounds, center = Bounds._bounds_and_center_of_chain(chain, assign_bounds=assign_attribute, bounds_key=key)
        return bounds

    @staticmethod
    def center_of_chain(chain: Chain,
                        assign_attribute: bool = False,
                        key: str = "center") -> Tuple[float, float, float]:
        """
        Returns the center of a chain.

        Args:
            chain (Chain): The chain.
            assign_attribute (bool): Whether to assign the center of the chain as an attribute.
            key (str): The key to use for the attribute.

        Returns:
            Tuple[float, float, float]: The center of the chain.1
        """
        bounds, center = Bounds._bounds_and_center_of_chain(chain, assign_center=assign_attribute, center_key=key)
        return center

    @staticmethod
    def _bounds_and_center_of_protein(protein: Protein,
                                      assign_bounds: bool = False,
                                      bounds_key: str = "bounds",
                                      assign_center: bool = False,
                                      center_key: str = "center") -> Tuple[Tuple[float, float, float, float, float, float], Tuple[float, float, float]]:
        """
        Returns the bounds and center of a protein.

        Args:
            protein (Protein): The protein.

        Returns:
            Tuple[Tuple[float, float, float, float, float, float], Tuple[float, float, float]]: The bounds and center of the protein.
        """
        min_x = min_y = min_z = float("inf")
        max_x = max_y = max_z = float("-inf")

        for chain in protein.chains:
            bounds, center = Bounds._bounds_and_center_of_chain(chain, assign_bounds=assign_bounds, bounds_key=bounds_key, assign_center=assign_center, center_key=center_key)
            min_x = min(min_x, bounds[0])
            min_y = min(min_y, bounds[1])
            min_z = min(min_z, bounds[2])
            max_x = max(max_x, bounds[3])
            max_y = max(max_y, bounds[4])
            max_z = max(max_z, bounds[5])

        cx, cy, cz = (min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2

        if assign_bounds:
            protein.set_attribute(bounds_key, (min_x, min_y, min_z, max_x, max_y, max_z))
        if assign_center:
            protein.set_attribute(center_key, (cx, cy, cz))

        return (min_x, min_y, min_z, max_x, max_y, max_z), (cx, cy, cz)

    @staticmethod
    def bounds_of_protein(protein: Protein,
                          assign_attribute: bool = False,
                          key: str = "bounds") -> Tuple[float, float, float, float, float, float]:
        """
        Returns the bounds of a protein.

        Args:
            protein (Protein): The protein.

        Returns:
            Tuple[float, float, float, float, float, float]: The bounds of the protein.
        """
        bounds, center = Bounds._bounds_and_center_of_protein(protein, assign_bounds=assign_attribute, bounds_key=key)
        return bounds

    @staticmethod
    def center_of_protein(protein: Protein,
                          assign_attribute: bool = False,
                          key: str = "center") -> Tuple[float, float, float]:
        """
        Returns the center of a protein.

        Args:
            protein (Protein): The protein.
            assign_attribute (bool): Whether to assign the center of the protein as an attribute.
            key (str): The key to use for the attribute.

        Returns:
            Tuple[float, float, float]: The center of the protein.
        """
        bounds, center = Bounds._bounds_and_center_of_protein(protein, assign_center=assign_attribute, center_key=key)
        return center
