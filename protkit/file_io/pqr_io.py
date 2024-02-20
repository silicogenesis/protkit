#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `PQRIO` to read and write protein structures in the PQR format.
"""

from typing import List, Optional
from protkit.file_io.pdb_io import PDBIO
from protkit.structure.protein import Protein


class PQRIO:
    @staticmethod
    def load(file_path: str,
             pdb_id: Optional[str] = None) -> List[Protein]:
        """
        Loads a protein from a PQR file.

        Args:
            file_path (str): The path to the PQR file.
            pdb_id (str): The PDB ID of the protein.

        Returns:
            Protein: The protein.
        """
        return PDBIO.load(file_path, is_pqr_format=True, pdb_id=pdb_id)

    @staticmethod
    def save(protein: Protein,
             file_path: str) -> None:
        """
        Saves a protein to a PQR file.

        Args:
            protein (Protein): The protein.
            file_path (str): The path to the PQR file.

        Returns:
            None
        """
        PDBIO.save(protein, file_path, is_pqr_format=True)