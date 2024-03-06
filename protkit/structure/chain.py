#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `Chain` to represent a chain in a protein.
"""

from __future__ import annotations
from typing import List, Dict, Set, Optional, TYPE_CHECKING, Any, Iterator, Union
from copy import deepcopy
from collections import defaultdict

from protkit.structure.residue import Residue
from protkit.structure.atom import Atom
from protkit.seq.sequence import Sequence

if TYPE_CHECKING:
    from protkit.structure.protein import Protein


class Chain:
    # ------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------
    def __init__(self,
                 # Core
                 chain_id: str,
                 # Optional
                 protein: Optional[Protein] = None
                 ) -> None:
        """
        Constructor for the Chain class.

        Args:
            chain_id (str): The chain's ID.
            protein (Optional[Protein]): The chain's protein.

        Returns:
            None
        """

        self._chain_id: str = chain_id
        self._residues: List[Residue] = []
        self._protein: Optional[Protein] = protein

    def copy(self):
        """
        Returns a deep copy of the chain, excluding the protein.

        Returns:
            Chain: A deep copy of the chain.
        """
        protein = self._protein
        self._protein = None
        chain = deepcopy(self)
        self._protein = protein
        return chain

    # ------------------------------------------------------------
    # Methods for managing the chain's ID.
    # - id (property)
    # - chain_id (property)
    # - chain_id (setter)
    # - sequence (property)
    # ------------------------------------------------------------

    @property
    def id(self) -> str:
        """
        Returns the chain's ID, which is the chain ID if the chain
        is not part of a protein, or the protein ID + chain ID if
        the chain is part of a protein.

        Returns:
            str: The chain's ID.
        """
        if self._protein is None:
            return self._chain_id
        elif self._protein.id is None or self._protein.id == "":
            return self._chain_id
        else:
            return self._protein.id + ":" + self._chain_id

    @property
    def chain_id(self) -> str:
        """
        Returns the chain's ID.

        Returns:
            str: The chain's ID.
        """
        return self._chain_id

    @chain_id.setter
    def chain_id(self, chain_id: str) -> None:
        """
        Sets the chain's ID.

        Args:
            chain_id (str): The chain's ID.

        Returns:
            None
        """
        self._chain_id = chain_id

    @property
    def sequence(self) -> str:
        """
        Returns the chain's sequence.

        Returns:
            str: The chain's sequence.
        """
        sequence = [residue.short_code for residue in self._residues]
        sequence = "".join(sequence)
        return sequence

    # ------------------------------------------------------------
    # Methods for managing the chain's protein.
    # - protein (property)
    # - protein (setter)
    # ------------------------------------------------------------

    @property
    def protein(self) -> Protein:
        """
        Returns the chain's protein.

        Returns:
            Protein: The chain's protein.
        """
        return self._protein

    @protein.setter
    def protein(self, protein: Protein) -> None:
        """
        Sets the chain's protein.

        Args:
            protein (Protein): The chain's protein.

        Returns:
            None
        """
        self._protein = protein

    # ------------------------------------------------------------
    # Methods for managing the chain's residues.
    # Create
    # - add_residue
    # Read
    # - residues (property, iterator)
    # - num_residues (property)
    # - num_hetero_residues (property)
    # - num_water_residues (property)
    # - num_disordered_residues (property)
    # - num_residues_by_type (property)
    # - get_residue
    # - filter_residues (iterator)
    # - seqres_analysis
    # - assign_segments
    # Update
    # - renumber_residues
    # Delete
    # - remove_hetero_residues
    # ------------------------------------------------------------

    def add_residue(self, residue: Residue) -> Residue:
        """
        Adds a residue to the end of the chain.

        Args:
            residue (Residue): The residue to be added.

        Returns:
            Residue: The residue that was added.
        """
        self._residues.append(residue)
        residue.chain = self
        return residue

    @property
    def residues(self) -> Iterator[Residue]:
        """
        Returns an iterator over the chain's residues.

        Yields:
            Residue: The next residue in the chain.
        """
        for residue in self._residues:
            yield residue

    @property
    def num_residues(self):
        """
        Returns the number of residues in the chain.

        Returns:
            int: The number of residues in the chain.
        """
        return len(self._residues)

    @property
    def num_disordered_residues(self) -> int:
        """
        Returns the number of disordered residues in the chain.

        Returns:
            int: The number of disordered residues in the chain.
        """
        count = 0
        for residue in self._residues:
            if residue.is_disordered:
                count += 1
        return count

    @property
    def num_hetero_residues(self) -> int:
        """
        Returns the number of hetero residues in the chain.

        Returns:
            int: The number of hetero residues in the chain.
        """
        count = 0
        for residue in self._residues:
            if residue.is_hetero:
                count += 1
        return count

    @property
    def num_water_residues(self) -> int:
        """
        Returns the number of water residues in the chain.

        Returns:
            int: The number of water residues in the chain.
        """
        count = 0
        for residue in self._residues:
            if residue.residue_type == "HOH":
                count += 1
        return count

    @property
    def num_residues_by_type(self) -> Dict[str, int]:
        """
        Returns the number of residues in the chain by residue type.

        Returns:
            Dict[str, int]: The number of residues in the chain by residue type.
        """
        num_resides = defaultdict(int)
        for residue in self._residues:
            num_resides[residue.residue_type] += 1
        return num_resides

    def get_residue(self, index) -> Residue:
        """
        Returns the residue at the specified index (0-based).

        Args:
            index (int): The index of the residue.

        Returns:
            Residue: The residue at the specified index.
        """
        if (index > len(self._residues) - 1) or (index < 0):
            raise IndexError("Residue index out of range.")

        return self._residues[index]

    def filter_residues(self, residue_criteria: Optional[List] = None) -> List[Residue]:
        """
        Yields residues that match the specified criteria.

        Args:
            residue_criteria (Optional[List]): List of residue criteria.

        Yields:
            Residue: The next residue that matches the specified criteria.
        """
        if residue_criteria is None:
            residue_criteria = []
        for residue in self._residues:
            match = True
            for key, value in residue_criteria:
                if type(value) is list:
                    if residue.get_attribute(key) not in value:
                        match = False
                        break
                elif residue.get_attribute(key) != value:
                    match = False
                    break
            if match:
                yield residue

    def seqres_analysis(self):
        """
        Compares the sequence of the chain to the SEQRES record.
        """
        sequence_residues = self.get_attribute("seqres")
        structure_residues = [residue.residue_type for residue in self._residues if not residue.is_hetero]
        alignment1, alignment2 = Sequence.align(sequence_residues, structure_residues)
        return alignment1, alignment2

    def renumber_residues(self, new_numbers: Optional[List[int]] = None) -> None:
        """
        Renumbers the residues in the chain. Residues are numbered
        sequentially, starting at 1.  If new_numbers is not specified,
        the residues are numbered sequentially, starting at 1.
        Insertion codes are removed.

        Args:
            new_numbers (Optional[List[int]]): List of new residue numbers.
                If not specified, the residues are numbered sequentially,
                starting at 1.

        Returns:
            None

        Raises:
            None
        """
        if new_numbers is None:
            new_numbers = list(range(1, len(self._residues) + 1))

        for i in range(len(self._residues)):
            self._residues[i].sequence_no = new_numbers[i]
            self._residues[i].insertion_code = ""

    def remove_hetero_residues(self, hetero_residue_types: Optional[Union[str, List[str]]] = None) -> None:
        """
        Removes hetero residues from the chain.

        Args:
            hetero_residue_types (Optional[Union[str, List[str]]]): The residue types to be removed.
                If not specified, all hetero residues are removed.

        Returns:
            None
        """
        if hetero_residue_types is None:
            self._residues = [residue for residue in self._residues if not residue.is_hetero]
        else:
            if type(hetero_residue_types) is str:
                hetero_residue_types = [hetero_residue_types]
            self._residues = [residue for residue in self._residues if residue.residue_type not in hetero_residue_types]

    def remove_water_residues(self) -> None:
        """
        Removes water residues from the chain.

        Returns:
            None
        """
        self._residues = [residue for residue in self._residues if not residue.residue_type == "HOH"]

    # ------------------------------------------------------------
    # Methods for managing the chain's atoms.
    # Create
    # Read
    # - atoms (property, iterator)
    # - num_atoms (property)
    # - num_heavy_atoms (property)
    # - num_hydrogen_atoms (property)
    # - num_disordered_atoms (property)
    # - num_hetero_atoms (property)
    # - filter_atoms (iterator)
    # - missing_heavy_atom_types
    # - extra_heavy_atom_types
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

    @property
    def atoms(self) -> Iterator[Atom]:
        """
        Returns an iterator over the chain's atoms.

        Yields:
            Atom: The next atom in the chain.
        """
        for residue in self._residues:
            for atom in residue.atoms:
                yield atom

    @property
    def num_atoms(self) -> int:
        """
        Returns the number of atoms in the chain.

        Returns:
            int: The number of atoms in the chain.
        """
        count = 0
        for residue in self._residues:
            count += residue.num_atoms
        return count

    @property
    def num_heavy_atoms(self) -> int:
        """
        Returns the number of heavy atoms in the chain.

        Returns:
            int: The number of heavy atoms in the chain.
        """
        count = 0
        for residue in self._residues:
            count += residue.num_heavy_atoms
        return count

    @property
    def num_hydrogen_atoms(self) -> int:
        """
        Returns the number of hydrogen atoms in the chain.

        Returns:
            int: The number of hydrogen atoms in the chain.
        """
        count = 0
        for residue in self._residues:
            count += residue.num_hydrogen_atoms
        return count

    @property
    def num_disordered_atoms(self) -> int:
        """
        Returns the number of disordered atoms in the chain.

        Returns:
            int: The number of disordered atoms in the chain.
        """
        count = 0
        for residue in self._residues:
            count += residue.num_disordered_atoms
        return count

    @property
    def num_hetero_atoms(self) -> int:
        """
        Returns the number of hetero atoms in the chain.

        Returns:
            int: The number of hetero atoms in the chain.
        """
        count = 0
        for residue in self._residues:
            count += residue.num_hetero_atoms
        return count

    def filter_atoms(self,
                     atom_criteria: Optional[List] = None,
                     residue_criteria: Optional[List] = None) -> Iterator[Atom]:
        """
        Yields atoms that match the specified criteria.

        Args:
            atom_criteria (Optional[List]): List of atom criteria.
            residue_criteria (Optional[List]): List of residue criteria.

        Yields:
            Atom: The next atom that matches the specified criteria.
        """
        for residue in self.filter_residues(residue_criteria):
            for atom in residue.filter_atoms(atom_criteria):
                yield atom

    def missing_heavy_atom_types(self) -> Dict[str, Set[str]]:
        """
        Returns a list of missing heavy atom types.

        Returns:
            Dict[str, Set[str]]: A list of missing heavy atom types.
        """
        missing_atom_types = {}
        for residue in self._residues:
            missing = residue.missing_heavy_atom_types()
            if len(missing) > 0:
                missing_atom_types[residue.residue_code] = missing
        return missing_atom_types

    def extra_heavy_atom_types(self) -> Dict[str, Set[str]]:
        """
        Returns a list of extra heavy atom types.

        Returns:
            Dict[str, Set[str]]: A list of extra heavy atom types.
        """
        extra_atom_types = {}
        for residue in self._residues:
            extra = residue.extra_heavy_atom_types()
            if len(extra) > 0:
                extra_atom_types[residue.residue_code] = extra
        return extra_atom_types

    def missing_hydrogen_atom_types(self) -> Dict[str, Set[str]]:
        """
        Returns a list of missing hydrogen atom types.

        Returns:
            Dict[str, Set[str]]: A list of missing hydrogen atom types.
        """
        missing_atom_types = {}
        for residue in self._residues:
            missing = residue.missing_hydrogen_atom_types()
            if len(missing) > 0:
                missing_atom_types[residue.residue_code] = missing
        return missing_atom_types

    def extra_hydrogen_atom_types(self) -> Dict[str, Set[str]]:
        """
        Returns a list of extra hydrogen atom types.

        Returns:
            Dict[str, Set[str]]: A list of extra hydrogen atom types.
        """
        extra_atom_types = {}
        for residue in self._residues:
            extra = residue.extra_hydrogen_atom_types()
            if len(extra) > 0:
                extra_atom_types[residue.residue_code] = extra
        return extra_atom_types

    def fix_disordered_atoms(self) -> None:
        """
        Fixes disordered atoms in the chain.

        Returns:
            None
        """
        for residue in self._residues:
            residue.fix_disordered_atoms()

    def remove_atoms(self, atom_types: List[str]) -> None:
        """
        Removes atoms from the chain.

        Args:
            atom_types (List[str]): The atom types to be removed.

        Returns:
            None
        """
        for residue in self._residues:
            residue.remove_atoms(atom_types)

    def remove_hydrogen_atoms(self):
        """
        Removes hydrogen atoms from the chain.

        Returns:
            None
        """
        for residue in self._residues:
            residue.remove_hydrogen_atoms()

    def keep_atoms(self, atom_types: List[str]) -> None:
        """
        Removes atoms from the chain.

        Args:
            atom_types (List[str]): The atom types to be kept.

        Returns:
            None
        """
        for residue in self._residues:
            residue.keep_atoms(atom_types)

    def keep_backbone_atoms(self) -> None:
        """
        Removes atoms from the chain.

        Returns:
            None
        """
        for residue in self._residues:
            residue.keep_backbone_atoms()

    # ------------------------------------------------------------
    # Methods for managing the chain's attributes.
    # - has_attribute
    # - get_attribute
    # - set_attribute
    # - delete_attribute
    # - sum_attribute
    # ------------------------------------------------------------

    def has_attribute(self, key: str) -> bool:
        """
        Returns True if the chain has the specified attribute.

        Args:
            key (str): The name of the attribute.

        Returns:
            bool: True if the chain has the specified attribute.
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
        if hasattr(self, "_" + key):
            delattr(self, "_" + key)

    def sum_attribute(self, key: str) -> None:
        """
        Sums the value of the specified attribute.
        """
        atrribute_sum = 0.0
        for residue in self._residues:
            if not hasattr(residue, key):
                residue.sum_attribute(key)
            atrribute_sum += residue.get_attribute(key)
        self.set_attribute(key, atrribute_sum)

    # ------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------

    def summary(self) -> str:
        text = [
            f"Chain ID: {self.id}",
            f"Number of residues: {self.num_residues}",
            f"Number of atoms: {self.num_atoms}"
        ]

        return "\n".join(text)
