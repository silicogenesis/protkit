#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS), Claudio Jardim (CJ)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `Protonator`.  The Protonator class serves as a template
for classes that are responsible for adding hydrogen atoms to proteins in
order to simulate the effects of protonation at specific pH levels.
"""

from abc import ABC, abstractmethod

from protkit.structure.protein import Protein


class Protonator(ABC):
    """
    Abstract base class for protonating proteins.

    The Protonator class serves as a template for classes that are responsible
    for adding hydrogen atoms to proteins in order to simulate the effects of
    protonation at specific pH levels.
    """

    def __init__(self):
        """
        Constructs all the necessary attributes for the Protonator object.
        """
        pass

    @abstractmethod
    def protonate(self, protein: Protein) -> Protein:
        """
        Abstract method for protonating a protein.

        This method should be implemented by subclasses to add hydrogens to the
        amino acid residues of a protein structure, typically according to a
        specified pH level.

        Parameters:
        protein (Protein): The protein to be protonated.

        Returns:
        Protein: The protonated protein.
        """
        pass