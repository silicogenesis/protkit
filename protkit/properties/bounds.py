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

from typing import Tuple
from protkit.structure.protein import Protein
from protkit.structure.chain import Chain
from protkit.structure.residue import Residue


class Bounds:
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
        min_x = min_y = min_z = float("inf")
        max_x = max_y = max_z = float("-inf")

        for atom in residue.atoms:
            min_x = min(min_x, atom.x)
            min_y = min(min_y, atom.y)
            min_z = min(min_z, atom.z)
            max_x = max(max_x, atom.x)
            max_y = max(max_y, atom.y)
            max_z = max(max_z, atom.z)

        if assign_attribute:
            residue.set_attribute(key, (min_x, max_x, min_y, max_y, min_z, max_z))

        return min_x, max_x, min_y, max_y, min_z, max_z

    @staticmethod
    def center_of_residue(residue: Residue,
                          assign_attribute: bool = False,
                          key: str = "center",
                          assign_bounds: bool = False,
                          bounds_key: str = "bounds") -> Tuple[float, float, float]:
        """
        Returns the center of a residue.

        Args:
            residue (Residue): The residue.
            assign_attribute (bool): Whether to assign the center of the residue as an attribute.
            key (str): The key to use for the attribute.
            assign_bounds (bool): Whether to assign the bounds of the residue as an attribute.
            bounds_key (str): The key to use for the attribute.

        Returns:
            Tuple[float, float, float]: The center of the residue.

        Raises:
            None
        """
        min_x, max_x, min_y, max_y, min_z, max_z = Bounds.bounds_of_residue(residue, assign_attribute=assign_bounds, key=bounds_key)
        center_x, center_y, center_z = (min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2
        if assign_attribute:
            residue.set_attribute(key, (center_x, center_y, center_z))
        return center_x, center_y, center_z

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
        min_x = min_y = min_z = float("inf")
        max_x = max_y = max_z = float("-inf")

        for residue in chain.residues:
            min_x, max_x, min_y, max_y, min_z, max_z = Bounds.bounds_of_residue(residue, assign_attribute=assign_bounds, key=key)
            min_x = min(min_x, min_x)
            min_y = min(min_y, min_y)
            min_z = min(min_z, min_z)
            max_x = max(max_x, max_x)
            max_y = max(max_y, max_y)
            max_z = max(max_z, max_z)

        if assign_attribute:
            chain.set_attribute(key, (min_x, max_x, min_y, max_y, min_z, max_z))

        return min_x, max_x, min_y, max_y, min_z, max_z

    @staticmethod
    def center_of_chain(chain: Chain,
                        assign_attribute: bool = False,
                        key: str = "center",
                        assign_bounds: bool = False,
                        bounds_key: str = "bounds") -> Tuple[float, float, float]:
        """
        Returns the center of a chain.

        Args:
            chain (Chain): The chain.
            assign_attribute (bool): Whether to assign the center of the chain as an attribute.
            key (str): The key to use for the attribute.
            assign_bounds (bool): Whether to assign the bounds of the chain as an attribute.
            bounds_key (str): The key to use for the attribute.

        Returns:
            Tuple[float, float, float]: The center of the chain.1
        """

        min_x, max_x, min_y, max_y, min_z, max_z = Bounds.bounds_of_chain(chain, assign_attribute=assign_bounds, key=bounds_key)
        center_x, center_y, center_z = (min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2
        if assign_attribute:
            chain.set_attribute(key, (center_x, center_y, center_z))
        return center_x, center_y, center_z

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
        min_x = min_y = min_z = float("inf")
        max_x = max_y = max_z = float("-inf")

        for chain in protein.chains:
            min_x, max_x, min_y, max_y, min_z, max_z = Bounds.bounds_of_chain(chain, assign_attribute=assign_attribute, key=key)
            min_x = min(min_x, min_x)
            min_y = min(min_y, min_y)
            min_z = min(min_z, min_z)
            max_x = max(max_x, max_x)
            max_y = max(max_y, max_y)
            max_z = max(max_z, max_z)

        if assign_attribute:
            protein.set_attribute(key, (min_x, max_x, min_y, max_y, min_z, max_z))

        return min_x, max_x, min_y, max_y, min_z, max_z

    @staticmethod
    def center_of_protein(protein: Protein,
                          assign_attribute: bool = False,
                          key: str = "center",
                          assign_bounds: bool = False,
                          bounds_key: str = "bounds") -> Tuple[float, float, float]:
        """
        Returns the center of a protein.

        Args:
            protein (Protein): The protein.
            assign_attribute (bool): Whether to assign the center of the protein as an attribute.
            key (str): The key to use for the attribute.
            assign_bounds (bool): Whether to assign the bounds of the protein as an attribute.
            bounds_key (str): The key to use for the attribute.

        Returns:
            Tuple[float, float, float]: The center of the protein.
        """

        min_x, max_x, min_y, max_y, min_z, max_z = Bounds.bounds_of_protein(protein, assign_attribute=assign_bounds, key=bounds_key)
        center_x, center_y, center_z = (min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2
        if assign_attribute:
            protein.set_attribute(key, (center_x, center_y, center_z))
        return center_x, center_y, center_z
