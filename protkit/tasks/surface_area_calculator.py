#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS), Claudio Jardim (CJ)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `SurfaceAreaCalculator` to specify the interface for classes
that calculate the accessible surface area (ASA) of proteins.
"""

from abc import ABC, abstractmethod
from typing import List

from protkit.structure import Atom
from protkit.structure import Protein


class SurfaceAreaCalculator(ABC):
    """
    Abstract base class for calculating molecular surface properties.

    This class defines the interface for calculating the accessible surface area (ASA)
    of proteins, which is essential for understanding protein interactions and properties.
    """

    @abstractmethod
    def calculate_surface_area(self, atoms: List[Atom]) -> List[float]:
        """
        Calculate and update the accessible surface area of a protein.

        Args:
            atoms (List[Atom]): A list of Atom objects whose surface area needs to be calculated.

        Returns:
            List[float]: A list of surface areas for each atom in the input list.
        """
        pass

    def calculate_surface_area_for_protein(self, protein: Protein) -> List[float]:
        """
        Calculate and update the accessible surface area of a protein.

        Args:
            protein (Protein): A Protein object whose surface area needs to be calculated.

        Returns:
            Protein: An updated Protein object with new surface area properties.
        """
        atoms = list(protein.atoms)
        return self.calculate_surface_area(atoms)
