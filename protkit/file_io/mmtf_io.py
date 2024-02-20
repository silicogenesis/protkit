#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `MMTFIO` to read and write protein structures in the MMTF format.

The MMTF format is a binary format that is used to store data from macromolecular
crystallography experiments. It is the successor to the PDB format and is the preferred
format for the PDB archive.

See https://mmtf.rcsb.org/ for more information.

<font color="red">Note that the current implementation relies on converting the MMTF
file to PDB format using BioPython (which in turn uses mmtf-python).
Metadata that are not supported by the PDB format
will be lost in the process. Future implementations of this class will correctly
handle additional metadata.</font>
"""
import tempfile
from typing import List

from Bio.PDB.mmtf import MMTFParser
from Bio.PDB import PDBIO as IO

from protkit.file_io.prot_io import ProtIO
from protkit.file_io.pdb_io import PDBIO
from protkit.structure.protein import Protein

class MMTFIO:
    @staticmethod
    def convert_to_pdb(mmtf_file_path: str, pdb_file_path: str) -> None:
        """
        Converts a MMTF file to a PDB file.

        Args:
            mmtf_file_path (str): The path to the MMTF file.
            pdb_file_path (str): The path to the PDB file.

        Returns:
            None

        Raises:
            FileNotFoundError: If the MMTF file could not be found.
        """
        try:
            parser = MMTFParser()
            structure = parser.get_structure(mmtf_file_path)
            io = IO()
            io.set_structure(structure)
            io.save(pdb_file_path)
        except FileNotFoundError as e:
            raise e

    @staticmethod
    def convert_to_prot(mmtf_file_path: str, prot_file_path: str) -> None:
        """
        Converts a MMTF file to a Prot file.

        Args:
            mmtf_file_path (str): The path to the MMTF file.
            prot_file_path (str): The path to the Prot file.

        Returns:
            None

        Raises:
            FileNotFoundError: If the MMTF file could not be found.
        """
        try:
            parser = MMTFParser()
            structure = parser.get_structure("structure", mmtf_file_path)
            io = IO()
            io.set_structure(structure)

            with tempfile.TemporaryDirectory() as tmp:
                pdb_file_path = tmp + "/tmp.pdb"
                io.save(pdb_file_path)
                ProtIO.convert(pdb_file_path, prot_file_path)

        except FileNotFoundError as e:
            raise e

    @staticmethod
    def load(mmtf_file_path: str) -> List[Protein]:
        """
        Loads a MMTF file and returns a protein.

        Args:
            mmtf_file_path (str): The path to the MMTF file.

        Returns:
            List[Protein]: A list of protein(s) in the MMTF file.

        Raises:
            FileNotFoundError: If the MMTF file could not be found.
        """
        try:
            parser = MMTFParser()
            structure = parser.get_structure(mmtf_file_path)
            io = IO()
            io.set_structure(structure)

            with tempfile.TemporaryDirectory() as tmp:
                pdb_file_path = tmp + "/tmp.pdb"
                io.save(pdb_file_path)
                return PDBIO.load(pdb_file_path)

        except FileNotFoundError as e:
            raise e
