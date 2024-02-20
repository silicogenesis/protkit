#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `ProtIO` to read and write protein structures in the Prot format.
The Prot format enables storage of attributes of the protein structure that are not
part of the PDB format. It closely resembles the internal data structure of the
Protein class.

The Prot format is a JSON format that is compressed using zlib to
reduce the file size. The JSON format is used because it is human-readable and
can be easily parsed by other programs.

Methods are static and can be called without instantiating the class. The main methods
are:

- `convert()` to convert a PDB file to a Prot file.
- `load()` to load a Prot file.
- `save()` to save a Prot file.
"""

import json
import zlib
from typing import List, Dict, Union, Optional

from protkit.file_io import PDBIO
from protkit.structure.protein import Protein
from protkit.structure.chain import Chain
from protkit.structure.residue import Residue
from protkit.structure.atom import Atom


class ProtIO:
    @staticmethod
    def convert(pdb_file_path: str,
                prot_file_path: str,
                is_pqr_format: bool = False,
                pdb_id: Optional[str] = None,
                compress: bool = True,
                indent: int = 2) -> Protein:
        """
        Converts a PDB file to a Prot file.

        Args:
            pdb_file_path (str): The path to the PDB file.
            prot_file_path (str): The path to the Prot file.
            is_pqr_format (bool): Whether the PDB file is in PQR format.
            pdb_id (str): The ID of the protein.
            compress (bool): Whether to compress the Prot file.
            indent (int): The indentation level of the Prot file.

        Returns:
            Protein: The protein.

        Raises:
            FileNotFoundError: If the PDB file could not be found.

        Examples:
            >> ProtIO.convert("1a00.pdb", "1a00.prot")
            >> ProtIO.convert("1a00.pdb", "1a00.prot", is_pqr_format=True)
            >> ProtIO.convert("1a00.pdb", "1a00.prot", pdb_id="1A00")
        """
        protein = PDBIO.load(pdb_file_path, is_pqr_format=is_pqr_format, pdb_id=pdb_id)[0]
        ProtIO.save(protein, prot_file_path, compress=compress, indent=indent)

        return protein

    @staticmethod
    def load(file_path: str,
             decompress: bool = True) -> List[Protein]:
        """
        Loads a protein from a Prot file.

        Args:
            file_path (str): The path to the Prot file.
            decompress (bool): Whether to decompress the Prot file.

        Returns:
            List[Protein]: The protein(s) in the Prot file.

        Raises:
            FileNotFoundError: If the Prot file could not be found.
        """
        if decompress:
            with open(file_path, "rb") as file:
                text = file.read()
            text = zlib.decompress(text)
        else:
            with open(file_path, "rt") as file:
                text = file.read()

        property_list = json.loads(text)

        proteins = [ProtIO.protein_from_dict(property_dict) for property_dict in property_list]
        return proteins

    @staticmethod
    def save(protein: Union[Protein, List[Protein]],
             file_path: str,
             compress: bool = True,
             indent: int = 2) -> None:
        """
        Saves a protein to a Prot file.

        Args:
            protein (Union[Protein, List[Protein]]): The protein(s) to save.
            file_path (str): The path to the Prot file.
            compress (bool): Whether to compress the Prot file.
            indent (int): The indentation level of the Prot file.

        Returns:
            None
        """
        if type(protein) is Protein:
            protein = [protein]

        if compress:
            text = json.dumps([ProtIO.protein_to_dict(prot) for prot in protein], indent=0)
            text = zlib.compress(text.encode('utf-8'))
            with open(file_path, "wb") as file:
                file.write(text)
        else:
            text = json.dumps([ProtIO.protein_to_dict(prot) for prot in protein], indent=indent)
            with open(file_path, "wt") as file:
                file.write(text)

    # ------------------------------------------------------------
    # Methods to convert between Prot and Python data structures
    # - atom_to_dict()
    # - atom_from_dict()
    # - residue_to_dict()
    # - residue_from_dict()
    # - chain_to_dict()
    # - chain_from_dict()
    # - protein_to_dict()
    # - protein_from_dict()
    # ------------------------------------------------------------

    @staticmethod
    def atom_to_dict(atom: Atom) -> Dict:
        """
        Converts an Atom object to a dictionary.

        Args:
            atom (Atom): The Atom object.

        Returns:
            Dict: The dictionary representation of the Atom object.

        Raises:
            None
        """
        props = dict()
        var = vars(atom)
        exclude = ["_residue"]
        for key in var:
            if key not in exclude:
                props[key[1:]] = var[key]
        return props

    @staticmethod
    def atom_from_dict(props: dict) -> Atom:
        """
        Converts a dictionary to an Atom object.

        Args:
            props (dict): The dictionary representation of the Atom object.

        Returns:
            Atom: The Atom object created from the dictionary.

        Raises:
            None
        """
        atom = Atom("", "", 0, 0, 0)
        for key in props:
            atom.set_attribute(key, props[key])
        return atom

    @staticmethod
    def residue_to_dict(residue: Residue) -> Dict:
        """
        Converts a Residue object to a dictionary.

        Args:
            residue (Residue): The Residue object.

        Returns:
            Dict: The dictionary representation of the Residue object.

        Raises:
            None
        """
        props = dict()
        var = vars(residue)
        exclude = ["_chain", "_atoms"]
        for key in var:
            if key not in exclude:
                props[key[1:]] = var[key]

        # Add atoms
        props["atoms"] = {}
        for atom in residue.atoms:
            props["atoms"][atom.atom_type] = ProtIO.atom_to_dict(atom)

        return props

    @staticmethod
    def residue_from_dict(props: dict) -> Residue:
        """
        Converts a dictionary to a Residue object.

        Args:
            props (dict): The dictionary representation of the Residue object.

        Returns:
            Residue: The Residue object created from the dictionary.

        Raises:
            None
        """
        residue = Residue("")
        exclude = ["chain", "atoms"]
        for key in props:
            if key not in exclude:
                residue.set_attribute(key, props[key])

        # Add atoms
        for key in props["atoms"]:
            atom = ProtIO.atom_from_dict(props["atoms"][key])
            residue.add_atom(key, atom)

        return residue

    @staticmethod
    def chain_to_dict(chain: Chain) -> Dict:
        """
        Converts a Chain object to a dictionary.

        Args:
            chain (Chain): The Chain object.

        Returns:
            Dict: The dictionary representation of the Chain object.

        Raises:
            None
        """
        props = dict()
        var = vars(chain)
        exclude = ["_protein", "_residues"]
        for key in var:
            if key not in exclude:
                props[key[1:]] = var[key]

        # Add residues
        props["residues"] = []
        for residue in chain.residues:
            props["residues"].append(ProtIO.residue_to_dict(residue))

        return props

    @staticmethod
    def chain_from_dict(props: dict) -> Chain:
        """
        Converts a dictionary to a Chain object.

        Args:
            props (dict): The dictionary representation of the Chain object.

        Returns:
            Chain: The Chain object created from the dictionary.

        Raises:
            None
        """
        chain = Chain("")
        exclude = ["protein", "residues"]
        for key in props:
            if key not in exclude:
                chain.set_attribute(key, props[key])

        # Add residues
        for res in props["residues"]:
            residue = ProtIO.residue_from_dict(res)
            chain.add_residue(residue)

        return chain

    @staticmethod
    def protein_to_dict(protein: Protein) -> Dict:
        """
        Converts a Protein object to a dictionary.

        Args:
            protein (Protein): The Protein object.

        Returns:
            Dict: The dictionary representation of the Protein object.

        Raises:
            None
        """
        props = dict()
        var = vars(protein)
        exclude = ["_chains"]
        for key in var:
            if key not in exclude:
                props[key[1:]] = var[key]

        # Add chains
        props["chains"] = {}
        for chain in protein.chains:
            props["chains"][chain.chain_id] = ProtIO.chain_to_dict(chain)

        return props

    @staticmethod
    def protein_from_dict(props: dict) -> Protein:
        """
        Converts a dictionary to a Protein object.

        Args:
            props (dict): The dictionary representation of the Protein object.

        Returns:
            Protein: The Protein object created from the dictionary.

        Raises:
            None
        """
        protein = Protein(None)
        exclude = ["chains"]
        for key in props:
            if key not in exclude:
                protein.set_attribute(key, props[key])

        # Add chains
        for chain in props["chains"]:
            chain = ProtIO.chain_from_dict(props["chains"][chain])
            protein.add_chain(chain.chain_id, chain)

        return protein
