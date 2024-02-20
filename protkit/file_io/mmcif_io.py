#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `MMCIFIO` to read and write protein structures in the mmCIF format.

The mmCIF format is a text-based format that is used to store data from macromolecular
crystallography experiments. It is the successor to the PDB format and is the preferred
format for the PDB archive.

See
https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/beginner%E2%80%99s-guide-to-pdb-structures-and-the-pdbx-mmcif-format
for more information.

<font color="red">Note that the current implementation relies on converting the mmCIF
file to PDB format using BioPython. Metadata that are not supported by the PDB format
will be lost in the process. Future implementations of this class will correctly
handle additional metadata.</font>
"""
import tempfile
from typing import List

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBIO as IO

from protkit.file_io.prot_io import ProtIO
from protkit.file_io.pdb_io import PDBIO
from protkit.structure.protein import Protein


class MMCIFIO:
    @staticmethod
    def convert_to_pdb(mmcif_file_path: str, pdb_file_path: str) -> None:
        """
        Converts a mmCIF file to a PDB file.

        Args:
            mmcif_file_path (str): The path to the mmCIF file.
            pdb_file_path (str): The path to the PDB file.

        Returns:
            None

        Raises:
            FileNotFoundError: If the mmCIF file could not be found.
        """
        try:
            parser = MMCIFParser()
            structure = parser.get_structure("structure", mmcif_file_path)
            io = IO()
            io.set_structure(structure)
            io.save(pdb_file_path)
        except FileNotFoundError as e:
            raise e

    @staticmethod
    def convert_to_prot(mmcif_file_path: str, prot_file_path: str) -> None:
        """
        Converts a mmCIF file to a Prot file.

        Args:
            mmcif_file_path (str): The path to the mmCIF file.
            prot_file_path (str): The path to the Prot file.

        Returns:
            None

        Raises:
            FileNotFoundError: If the mmCIF file could not be found.
        """
        try:
            parser = MMCIFParser()
            structure = parser.get_structure("structure", mmcif_file_path)
            io = IO()
            io.set_structure(structure)

            with tempfile.TemporaryDirectory() as tmp:
                pdb_file_path = tmp + "/tmp.pdb"
                io.save(pdb_file_path)
                ProtIO.convert(pdb_file_path, prot_file_path)

        except FileNotFoundError as e:
            raise e

    @staticmethod
    def load(mmcif_file_path: str) -> List[Protein]:
        """
        Loads a protein from a mmCIF file.

        Args:
            mmcif_file_path (str): The path to the mmCIF file.

        Returns:
            List[Protein]: The protein(s) in the mmCIF file.

        Raises:
            FileNotFoundError: If the mmCIF file could not be found.
        """
        try:
            parser = MMCIFParser()
            structure = parser.get_structure("structure", mmcif_file_path)
            io = IO()
            io.set_structure(structure)

            with tempfile.TemporaryDirectory() as tmp:
                pdb_file_path = tmp + "/tmp.pdb"
                io.save(pdb_file_path)
                return PDBIO.load(pdb_file_path)

        except FileNotFoundError as e:
            raise e