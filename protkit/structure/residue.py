#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `Residue` to represent a residue in a protein.
"""

from __future__ import annotations
from typing import Dict, List, Set, TYPE_CHECKING, Optional, Any, Iterator
from copy import deepcopy

from protkit.structure.atom import Atom

if TYPE_CHECKING:
    from protkit.structure.chain import Chain


class Residue:
    HEAVY_ATOMS = {
        "ALA": {"N", "CA", "C", "O", "CB"},
        "ARG": {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
        "ASN": {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"},
        "ASP": {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"},
        "CYS": {"N", "CA", "C", "O", "CB", "SG"},
        "GLU": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"},
        "GLN": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"},
        "GLY": {"N", "CA", "C", "O"},
        "HIS": {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"},
        "ILE": {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"},
        "LEU": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"},
        "LYS": {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"},
        "MET": {"N", "CA", "C", "O", "CB", "CG", "SD", "CE"},
        "PHE": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
        "PRO": {"N", "CA", "C", "O", "CB", "CG", "CD"},
        "SER": {"N", "CA", "C", "O", "CB", "OG"},
        "THR": {"N", "CA", "C", "O", "CB", "OG1", "CG2"},
        "TRP": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
        "TYR": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
        "VAL": {"N", "CA", "C", "O", "CB", "CG1", "CG2"}
    }

    HYDROGEN_ATOMS = {
        "ALA": {"H", "HA", "HB1", "HB2", "HB3"},
        "ARG": {"H", "HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "HE", "HH11", "HH12", "HH21", "HH22"},
        "ASN": {"H", "HA", "HB2", "HB3", "HD21", "HD22"},
        "ASP": {"H", "HA", "HB2", "HB3", "HD2"},
        "CYS": {"H", "HA", "HB2", "HB3", "HG"},
        "GLU": {"H", "HA", "HB2", "HB3", "HG2", "HG3", "HE2"},
        "GLN": {"H", "HA", "HB2", "HB3", "HG2", "HG3", "HE21", "HE22"},
        "GLY": {"H", "HA", "HA2", "HA3"},
        "HIS": {"H", "HA", "HB2", "HB3", "HD1", "HD2", "HE1", "HE2"},
        "ILE": {"H", "HA", "HB", "HG12", "HG13", "HG21", "HG22", "HG23", "HD11", "HD12", "HD13"},
        "LEU": {"H", "HA", "HB2", "HB3", "HG", "HD11", "HD12", "HD13", "HD21", "HD22", "HD23"},
        "LYS": {"H", "HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "HE2", "HE3", "HZ1", "HZ2", "HZ3"},
        "MET": {"H", "HA", "HB2", "HB3", "HG2", "HG3", "HE1", "HE2", "HE3"},
        "PHE": {"H", "HA", "HB2", "HB3", "HD1", "HD2", "HE1", "HE2", "HZ"},
        "PRO": {"H", "HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3"},
        "SER": {"H", "HA", "HB2", "HB3", "HG"},
        "THR": {"H", "HA", "HB", "HG1", "HG21", "HG22", "HG23"},
        "TRP": {"H", "HA", "HB2", "HB3", "HD1", "HE1", "HE3", "HZ2", "HZ3", "HH2"},
        "TYR": {"H", "HA", "HB2", "HB3", "HD1", "HD2", "HE1", "HE2", "HH"},
        "VAL": {"H", "HA", "HB", "HG11", "HG12", "HG13", "HG21", "HG22", "HG23"}
    }

    SHORT_CODE = {
        "ALA": "A",
        "CYS": "C",
        "ASP": "D",
        "GLU": "E",
        "PHE": "F",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LYS": "K",
        "LEU": "L",
        "MET": "M",
        "ASN": "N",
        "PRO": "P",
        "GLN": "Q",
        "ARG": "R",
        "SER": "S",
        "THR": "T",
        "VAL": "V",
        "TRP": "W",
        "TYR": "Y"
    }

    # ------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------
    def __init__(self,
                 # Core
                 residue_type: str,

                 # PDB properties
                 sequence_no: Optional[int] = None,
                 insertion_code: Optional[str] = None,

                 # Optional
                 chain: Optional[Chain] = None):
        """
        Constructor for the Residue class.

        Args:
            residue_type (str): The residue type.
            sequence_no (Optional[int]): The residue's sequence number.
            insertion_code (Optional[str]): The residue's insertion code.
            chain (Optional[Chain]): The residue's chain.

        Returns:
            None
        """

        self._residue_type: str = residue_type
        if sequence_no is not None:
            self._sequence_no = sequence_no
        if insertion_code is not None:
            self._insertion_code = insertion_code

        self._atoms: Dict[str, Atom] = {}
        self._chain: Chain = chain

    def copy(self):
        """
        Returns a deep copy of the residue, excluding the chain.

        Returns:
            Residue: A deep copy of the residue, excluding the chain.
        """
        chain = self._chain
        self._chain = None
        residue = deepcopy(self)
        self._chain = chain
        return residue

    # ------------------------------------------------------------
    # Methods for managing the residue's core properties.
    # - id (property)
    # - residue_type (property)
    # - residue_type (setter)
    # - residue_code (property)
    # - sequence_code (property)
    # - sequence_no (property)
    # - sequence_no (setter)
    # - insertion_code (property)
    # - insertion_code (setter)
    # - is_disordered (property)
    # - is_hetero (property)
    # - short_code (property)
    # ------------------------------------------------------------

    @property
    def id(self) -> str:
        """
        Returns the residue's ID, which is the chain ID + sequence number.

        Returns:
            str: The residue's ID.
        """
        if self._chain is None:
            return self.residue_code
        else:
            return self._chain.id + ":" + self.residue_code

    @property
    def residue_type(self) -> str:
        """
        Returns the residue type.

        Returns:
            str: The residue type.
        """
        return self._residue_type

    @residue_type.setter
    def residue_type(self, residue_type: str):
        """
        Sets the residue type.

        Args:
            residue_type (str): The residue type.

        Returns:
            None
        """
        self._residue_type = residue_type

    @property
    def residue_code(self) -> str:
        """
        Returns the residue's residue code.

        Returns:
            str: The residue's residue code.
        """
        return str(self._residue_type) + self.sequence_code

    @property
    def sequence_code(self) -> str:
        """
        Returns the residue's sequence code.

        Returns:
            str: The residue's sequence code.
        """
        return str(self._sequence_no) + self._insertion_code

    @property
    def sequence_no(self) -> int:
        """
        Returns the residue's sequence number.

        Returns:
            int: The residue's sequence number.
        """
        return self._sequence_no

    @sequence_no.setter
    def sequence_no(self, sequence_no: int) -> None:
        """
        Sets the residue's sequence number.

        Args:
            sequence_no (int): The residue's sequence number.

        Returns:
            None
        """
        self._sequence_no = sequence_no

    @property
    def insertion_code(self) -> str:
        """
        Returns the residue's insertion code.

        Returns:
            str: The residue's insertion code.
        """
        return self._insertion_code

    @insertion_code.setter
    def insertion_code(self, insertion_code: str) -> None:
        """
        Sets the residue's insertion code.

        Args:
            insertion_code (str): The residue's insertion code.

        Returns:
            None
        """
        self._insertion_code = insertion_code

    @property
    def is_disordered(self) -> bool:
        """
        Returns True if the residue is a disordered residue.

        Returns:
            bool: True if the residue is a disordered residue.
        """
        for atom in self._atoms.values():
            if atom.is_disordered:
                return True
        return False

    @property
    def is_hetero(self) -> bool:
        """
        Returns True if the residue is a hetero residue.

        Returns:
            bool: True if the residue is a hetero residue.
        """
        for atom in self._atoms.values():
            if atom.is_hetero:
                return True
        return False

    @property
    def short_code(self) -> str:
        """
        Returns the residue's short code.

        Returns:
            str: The residue's short code.

        Raises:
            Exception: If the residue type is not found in Residue.SHORT_CODE.
        """
        if self._residue_type not in Residue.SHORT_CODE:
            raise Exception("Residue type {self._residue_type} not found in Residue.SHORT_CODE.")

        return Residue.SHORT_CODE[self._residue_type]

    # ------------------------------------------------------------
    # Methods for managing the residue's chain.
    # - chain (property)
    # - chain (setter)
    # ------------------------------------------------------------

    @property
    def chain(self) -> Chain:
        """
        Returns the residue's chain.

        Returns:
            Chain: The residue's chain.
        """
        return self._chain

    @chain.setter
    def chain(self, chain: Chain) -> None:
        """
        Sets the residue's chain.

        Args:
            chain (Chain): The residue's chain.

        Returns:
            None
        """
        self._chain = chain

    # ------------------------------------------------------------
    # Methods for managing the residue's atoms.
    # Create
    # - add_atom
    # Read
    # - atoms (property, iterator)
    # - num_atoms (property)
    # - num_heavy_atoms (property)
    # - num_hydrogen_atoms (property)
    # - num_disordered_atoms (property)
    # - num_hetero_atoms (property)
    # - get_atom
    # - filter_atoms (iterator)
    # - heavy_atom_types
    # - missing_heavy_atom_types
    # - extra_heavy_atom_types
    # - hydrogen_atom_types
    # - missing_hydrogen_atom_types
    # - extra_hydrogen_atom_types
    # Update
    # - fix_disordered_atoms
    # Delete
    # - remove_atoms
    # - remove_hydrogen_atoms
    # - keep_atoms
    # - keep_backbone_atoms
    # ------------------------------------------------------------

    def add_atom(self, atom_id: str, atom: Atom) -> Atom:
        """
        Adds an atom to the residue.

        Args:
            atom_id (str): The atom's ID.
            atom (Atom): The atom.

        Returns:
            Atom: The atom.
        """
        if atom_id in self._atoms:
            current_atom = self._atoms[atom_id]
            if current_atom.is_disordered and atom.is_disordered:
                current_atom.merge_disordered_atom(atom)
                return current_atom
            else:
                raise Exception("Atom ID {atom_id} already exists in residue.")
        else:
            self._atoms[atom_id] = atom
            atom.residue = self
            return atom

    @property
    def atoms(self) -> Iterator[Atom]:
        """
        Returns an iterator of the residue's atoms.

        Yields:
            Atom: The next atom in the residue.
        """
        for atom in self._atoms.values():
            yield atom

    @property
    def num_atoms(self) -> int:
        """
        Returns the number of atoms in the residue.

        Returns:
            int: The number of atoms in the residue.
        """
        return len(self._atoms)

    @property
    def num_heavy_atoms(self) -> int:
        """
        Returns the number of heavy atoms in the residue.

        Returns:
            int: The number of heavy atoms in the residue.
        """
        count = 0
        for atom in self._atoms.values():
            if atom.element != "H":
                count += 1
        return count

    @property
    def num_hydrogen_atoms(self) -> int:
        """
        Returns the number of hydrogen atoms in the residue.

        Returns:
            int: The number of hydrogen atoms in the residue.
        """
        count = 0
        for atom in self._atoms.values():
            if atom.element == "H":
                count += 1
        return count

    @property
    def num_disordered_atoms(self) -> int:
        """
        Returns the number of disordered atoms in the residue.

        Returns:
            int: The number of disordered atoms in the residue.
        """
        count = 0
        for atom in self._atoms.values():
            if atom.is_disordered:
                count += 1
        return count

    @property
    def num_hetero_atoms(self) -> int:
        """
        Returns the number of hetero atoms in the residue.

        Returns:
            int: The number of hetero atoms in the residue.
        """
        count = 0
        for atom in self._atoms.values():
            if atom.is_hetero:
                count += 1
        return count

    def get_atom(self, atom_type: str) -> Optional[Atom]:
        """
        Returns the atom with the specified atom type.

        Args:
            atom_type (str): The atom type.

        Returns:
            Atom: The atom with the specified atom type.
        """
        if atom_type in self._atoms:
            return self._atoms[atom_type]
        else:
            return None

    def filter_atoms(self, atom_criteria: Optional[List] = None) -> Iterator[Atom]:
        """
        Returns an iterator of atoms that match the specified criteria.

        Args:
            atom_criteria (Optional[List]): The criteria for selecting atoms.

        Yields:
            Atom: The next atom that matches the criteria.
        """
        if atom_criteria is None:
            atom_criteria = []
        for atom in self._atoms.values():
            match = True
            for key, value in atom_criteria:
                if type(value) is list:
                    if atom.get_attribute(key) not in value:
                        match = False
                        break
                elif atom.get_attribute(key) != value:
                    match = False
                    break
            if match:
                yield atom

    def heavy_atom_types(self) -> Set[str]:
        """
        Returns a set of heavy atom types.

        Returns:
            Set[str]: A set of heavy atom types.
        """
        atom_types = set()
        for atom in self._atoms.values():
            if atom.element != "H":
                atom_types.add(atom.atom_type)
        return atom_types

    def missing_heavy_atom_types(self) -> Set[str]:
        """
        Returns a set of missing heavy atom types.

        Returns:
            Set[str]: A set of missing heavy atom types.
        """
        heavy_atom_types = self.heavy_atom_types()
        return Residue.HEAVY_ATOMS[self._residue_type] - heavy_atom_types

    def extra_heavy_atom_types(self) -> Set[str]:
        """
        Returns a set of extra heavy atom types.

        Returns:
            Set[str]: A set of extra heavy atom types.
        """
        heavy_atom_types = self.heavy_atom_types()
        return heavy_atom_types - Residue.HEAVY_ATOMS[self._residue_type]

    def hydrogen_atom_types(self) -> Set[str]:
        """
        Returns a set of hydrogen atom types.

        Returns:
            Set[str]: A set of hydrogen atom types.
        """
        atom_types = set()
        for atom in self._atoms.values():
            if atom.element == "H":
                atom_types.add(atom.atom_type)
        return atom_types

    def missing_hydrogen_atom_types(self) -> Set[str]:
        """
        Returns a set of missing hydrogen atom types.

        Returns:
            Set[str]: A set of missing hydrogen atom types.
        """
        hydrogen_atom_types = self.hydrogen_atom_types()
        return Residue.HYDROGEN_ATOMS[self._residue_type] - hydrogen_atom_types

    def extra_hydrogen_atom_types(self) -> Set[str]:
        """
        Returns a set of extra hydrogen atom types.

        Returns:
            Set[str]: A set of extra hydrogen atom types.
        """
        hydrogen_atom_types = self.hydrogen_atom_types()
        return hydrogen_atom_types - Residue.HYDROGEN_ATOMS[self._residue_type]

    def fix_disordered_atoms(self):
        """
        Fixes disordered atoms by keeping the atom with the
        highest occupancy.

        Returns:
            None
        """
        for atom in self._atoms.values():
            atom.fix_disordered_atom()

    def remove_atoms(self, atom_types: List[str]):
        """
        Removes the specified atoms from the residue.

        Args:
            atom_types (List[str]): The atom types to remove.

        Returns:
            None
        """
        if atom_types is None:
            atom_types = []
        for atom_type in atom_types:
            if atom_type in self._atoms:
                del self._atoms[atom_type]

    def remove_hydrogen_atoms(self):
        """
        Removes hydrogen atoms from the residue.

        Returns:
            None
        """
        keys = list(self._atoms.keys())
        for key in keys:
            if self._atoms[key].element == "H":
                del self._atoms[key]

        # The implementation below fails because the dictionary
        # changes size during iteration.
        # for atom in self._atoms.values():
        #     if atom.element == "H":
        #         del self._atoms[atom.atom_type]

    def keep_atoms(self, atom_types: List[str]) -> None:
        """
        Removes all atoms except the specified atoms from the residue.

        Args:
            atom_types (List[str]): The atom types to keep.

        Returns:
            None
        """
        if atom_types is None:
            atom_types = []
        atoms = {}
        for atom_type in atom_types:
            if atom_type in self._atoms:
                atoms[atom_type] = self._atoms[atom_type]
        self._atoms = atoms

    def keep_backbone_atoms(self) -> None:
        """
        Removes all atoms except the backbone atoms from the residue.
        The backbone atoms are N, CA, C, and O.

        Returns:
            None
        """
        self.keep_atoms(["N", "CA", "C", "O"])

    # ------------------------------------------------------------
    # Methods for managing the residue's attributes.
    # - has_attribute
    # - get_attribute
    # - set_attribute
    # - delete_attribute
    # - sum_attribute
    # ------------------------------------------------------------

    def has_attribute(self, key: str) -> bool:
        """
        Returns True if the residue has the specified attribute.

        Args:
            key (str): The name of the attribute.

        Returns:
            bool: True if the residue has the specified attribute.
        """
        return hasattr(self, "_" + key)

    def get_attribute(self, key: str) -> Any:
        """
        Returns the value of the specified attribute.

        Args:
            key (str): The name of the attribute.

        Returns:
            Any: The value of the specified attribute.

        Raises:
            AttributeError: If the attribute does not exist.
        """
        if hasattr(self, "_" + key):
            return getattr(self, "_" + key)

    def set_attribute(self, key: str, value: Any) -> None:
        """
        Sets the value of the specified attribute.

        Args:
            key (str): The name of the attribute.
            value (Any): The value of the attribute.

        Returns:
            None
        """
        if value is not None:
            setattr(self, "_" + key, value)

    def delete_attribute(self, key: str) -> None:
        """
        Deletes the specified attribute.

        Args:
            key (str): The name of the attribute.

        Returns:
            None
        """
        if self.has_attribute(key):
            delattr(self, "_" + key)

    def sum_attribute(self, key: str) -> None:
        """
        Sums the values of the specified attribute.
        """
        atrribute_sum = 0.0
        for atom in self._atoms.values():
            atrribute_sum += atom.get_attribute(key)
        self.set_attribute(key, atrribute_sum)
