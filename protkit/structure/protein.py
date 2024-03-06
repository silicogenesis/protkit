#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `Protein` to represent a protein's structural data.
"""

from typing import Dict, List, Set, Optional, Any, Iterator, Union
from copy import deepcopy
from collections import defaultdict

from protkit.structure.chain import Chain
from protkit.structure.residue import Residue
from protkit.structure.atom import Atom


class Protein:
    # ------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------

    def __init__(self,
                 pdb_id: Optional[str] = None) -> None:
        """
        Constructor for the Protein class.

        Args:
            pdb_id (str): The PDB ID of the protein.

        Returns:
            None
        """

        #: The PDB ID of the protein.
        self._pdb_id: Optional[str] = pdb_id

        #: The chains in the protein.
        self._chains: Dict[str, Chain] = dict()

    def copy(self,
             keep_chain_ids: Optional[Union[List[str], str]] = None,
             remove_chain_ids: Optional[Union[List[str], str]] = None):
        """
        Returns a deep copy of the protein.

        Args:
            keep_chain_ids (list of str): The chain IDs to keep.
            remove_chain_ids (list of str): The chain IDs to remove.

        Returns:
            Protein: A deep copy of the protein.
        """
        protein = deepcopy(self)
        if keep_chain_ids is not None:
            protein.keep_chains(keep_chain_ids)
        if remove_chain_ids is not None:
            protein.remove_chains(remove_chain_ids)

        return protein

    # ------------------------------------------------------------
    # Methods for managing the protein's ID.
    # - id (property)
    # - pdb_id (property)
    # - pdb_id (setter)
    # ------------------------------------------------------------

    @property
    def id(self) -> str:
        """
        Returns the protein's ID.

        Returns:
            str: The protein's ID.
        """
        return "" if self._pdb_id is None else self._pdb_id

    @property
    def pdb_id(self) -> Optional[str]:
        """
        Returns the protein's PDB ID.

        Returns:
            str: The protein's PDB ID.
        """
        return self._pdb_id

    @pdb_id.setter
    def pdb_id(self, pdb_id: str) -> None:
        """
        Sets the protein's PDB ID.

        Args:
            pdb_id (str): The protein's PDB ID.

        Returns:
            None
        """
        self._pdb_id = pdb_id

    # ------------------------------------------------------------
    # Methods for managing the protein's chains.
    # Create
    # - create_chain
    # - add_chain
    # Read
    # - chains (property, iterator)
    # - chain_ids (property)
    # - num_chains (property)
    # - has_chain
    # - get_chain
    # - filter_chains (iterator)
    # Update
    # - rename_chain
    # Delete
    # - remove_chains
    # - keep_chains
    # ------------------------------------------------------------

    def create_chain(self, chain_id: str) -> Chain:
        """
        Creates a new empty chain in the protein
        with the specified chain ID.

        Args:
            chain_id (str): The ID of the chain.

        Returns:
            Chain: The new chain.
        """
        if chain_id in self._chains:
            raise Exception(f"Chain {chain_id} already exists")
        else:
            self._chains[chain_id] = Chain(chain_id, self)
            return self._chains[chain_id]

    def add_chain(self, chain_id: str, chain: Chain) -> Chain:
        """
        Adds an existing chain to the protein.

        Args:
            chain_id (str): The ID of the chain.
            chain (Chain): The chain to add.

        Returns:
            Chain: The added chain.
        """
        if chain_id in self._chains:
            raise Exception(f"Chain {chain_id} already exists")

        self._chains[chain_id] = chain
        chain.chain_id = chain_id
        chain.protein = self
        return self._chains[chain_id]

    @property
    def chains(self) -> Iterator[Chain]:
        """
        Yields all chains in the protein.

        Yields:
            Chain: The next chain in the protein.
        """
        for chain in self._chains.values():
            yield chain

    @property
    def chain_ids(self) -> List[str]:
        """
        Returns a sorted list of chain IDs in the protein.

        Returns:
            list of str: The sorted list of chain IDs in the protein.
        """
        return sorted(list(self._chains.keys()))

    @property
    def num_chains(self) -> int:
        """
        Returns the number of chains in the protein.

        Returns:
            int: The number of chains in the protein.
        """
        return len(self._chains)

    def has_chain(self, chain_id: str) -> bool:
        """
        Checks if a chain exists in the protein.

        Args:
            chain_id (str): The ID of the chain.

        Returns:
            bool: True if the chain exists in the protein.
        """
        return chain_id in self._chains

    def get_chain(self, chain_id: str) -> Optional[Chain]:
        """
        Returns a chain from the protein if it exists.

        Args:
            chain_id (str): The ID of the chain.

        Returns:
            Chain: The chain if it exists, otherwise None.
        """
        if chain_id in self._chains:
            return self._chains[chain_id]
        else:
            return None

    def filter_chains(self, chain_criteria: Optional[List] = None) -> Iterator[Chain]:
        """
        Filters chains in the protein based on the specified criteria.
        Yielded chains must match all criteria.

        Args:
            chain_criteria (list of tuple): The criteria to filter the chains.

        Yields:
            Chain: The next chain in the protein that matches the criteria.
        """
        if chain_criteria is None:
            chain_criteria = []
        for chain in self._chains.values():
            match = True
            for key, value in chain_criteria:
                if type(value) is list:
                    if chain.get_attribute(key) not in value:
                        match = False
                        break
                elif chain.get_attribute(key) != value:
                    match = False
                    break
            if match:
                yield chain

    def rename_chain(self, chain_id_from: str, chain_id_to: str) -> None:
        """
        Renames a chain in the protein.

        Args:
            chain_id_from (str): The ID of the chain to rename.
            chain_id_to (str): The new ID of the chain.

        Returns:
            None
        """
        if chain_id_from not in self._chains:
            raise Exception(f"Chain {chain_id_from} does not exist")
        elif chain_id_to in self._chains:
            raise Exception(f"Chain {chain_id_to} already exists")
        else:
            self._chains[chain_id_to] = self._chains[chain_id_from]
            self._chains[chain_id_to].chain_id = chain_id_to
            self._chains.pop(chain_id_from)

    def remove_chains(self, chain_ids: Union[List[str], str]) -> None:
        """
        Removes chains from the protein.

        Args:
            chain_ids (list of str): The IDs of the chains to remove.

        Returns:
            None
        """
        if type(chain_ids) is str:
            chain_ids = [chain_ids]
        for chain_id in chain_ids:
            if chain_id not in self._chains:
                raise Exception(f"Chain {chain_id} does not exist")
            else:
                self._chains.pop(chain_id)

    def keep_chains(self, chain_ids: Union[List[str], str]) -> None:
        """
        Keeps only the selected chains in the protein.

        Args:
            chain_ids (list of str): The IDs of the chains to keep.

        Returns:
            None
        """
        if type(chain_ids) is str:
            chain_ids = [chain_ids]
        remove_chain_ids = list(set(self._chains.keys()) - set(chain_ids))
        self.remove_chains(remove_chain_ids)

    # ------------------------------------------------------------
    # Methods for managing the protein's residues.
    # Create
    # Read
    # - residues (property, iterator)
    # - num_residues (property)
    # - num_disordered_residues (property)
    # - num_hetero_residues (property)
    # - num_water_residues (property)
    # - num_residues_by_type (property)
    # - get_residue()
    # - filter_residues (iterator)
    # Update
    # - renumber_residues
    # Delete(keep_chain_ids=["A"])
    # - remove_hetero_residues
    # - remove_water_residues
    # ------------------------------------------------------------

    @property
    def residues(self) -> Iterator[Residue]:
        """
        Yields all residues in the protein.

        Yields:
            Residue: The next residue in the protein.
        """
        for chain in self.chains:
            for residue in chain.residues:
                yield residue

    @property
    def num_residues(self) -> int:
        """
        Returns the number of residues in the protein.

        Returns:
            int: The number of residues in the protein.
        """
        count = 0
        for chain in self._chains.values():
            count += chain.num_residues
        return count

    @property
    def num_disordered_residues(self) -> int:
        """
        Returns the number of disordered residues in the protein.

        Returns:
            int: The number of disordered residues in the protein.
        """
        count = 0
        for chain in self._chains.values():
            count += chain.num_disordered_residues
        return count

    @property
    def num_hetero_residues(self) -> int:
        """
        Returns the number of hetero residues in the protein.

        Returns:
            int: The number of hetero residues in the protein.
        """
        count = 0
        for chain in self._chains.values():
            count += chain.num_hetero_residues
        return count

    @property
    def num_water_residues(self) -> int:
        """
        Returns the number of water residues in the protein.

        Returns:
            int: The number of water residues in the protein.
        """
        count = 0
        for chain in self._chains.values():
            count += chain.num_water_residues
        return count

    @property
    def num_residues_by_type(self) -> Dict[str, int]:
        """
        Returns the number of residues in the protein by residue type.

        Returns:
            dict of str, int: The number of residues in the protein by residue type.
        """
        num_residues = defaultdict(int)
        for residue in self.residues:
            num_residues[residue.residue_type] += 1
        return num_residues

    def get_residue(self, chain_id: str, residue_index: int) -> Optional[Residue]:
        """
        Returns a residue from the protein if it exists.

        Args:
            chain_id (str): The ID of the chain.
            residue_index (int): The index of the residue.

        Returns:
            Residue: The residue if it exists, otherwise None.
        """
        if chain_id in self._chains:
            return self._chains[chain_id].get_residue(residue_index)
        else:
            return None

    def filter_residues(self,
                        chain_criteria: Optional[List] = None,
                        residue_criteria: Optional[List] = None) -> Iterator[Residue]:
        """
        Filters residues in the protein based on the specified criteria.
        Yielded residues must match all criteria.

        Args:
            chain_criteria (list of tuple): The criteria to filter the chains.
            residue_criteria (list of tuple): The criteria to filter the residues.

        Yields:
            Residue: The next residue in the protein that matches the criteria.
        """
        for chain in self.filter_chains(chain_criteria):
            for residue in chain.filter_residues(residue_criteria):
                yield residue

    def renumber_residues(self) -> None:
        """
        Renumbers residues in the protein sequentially.

        Returns:
            None
        """
        for chain in self.chains:
            chain.renumber_residues()

    def remove_hetero_residues(self, hetero_residue_types: Optional[Union[str, List[str]]] = None):
        """
        Removes all hetero residues from the protein.

        Args:
            hetero_residue_types (list of str): The residue types to remove.

        Returns:
            None
        """
        for chain in self.chains:
            chain.remove_hetero_residues(hetero_residue_types)

    def remove_water_residues(self):
        """
        Removes all water residues from the protein.

        Returns:
            None
        """
        for chain in self.chains:
            chain.remove_water_residues()

    # ------------------------------------------------------------
    # Methods for managing the protein's atoms.
    # Create
    # Read
    # - atoms (property, iterator)
    # - num_atoms (property)
    # - num_heavy_atoms (property)
    # - num_hydrogen_atoms (property)
    # - num_disordered_atoms (property)
    # - num_hetero_atoms (property)
    # - get_atom()
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
        Yields all atoms in the protein.

        Yields:
            Atom: The next atom in the protein.
        """
        for residue in self.residues:
            for atom in residue.atoms:
                yield atom

    @property
    def num_atoms(self) -> int:
        """
        Returns the number of atoms in the protein.

        Returns:
            int: The number of atoms in the protein.
        """
        count = 0
        for chain in self._chains.values():
            count += chain.num_atoms
        return count

    @property
    def num_heavy_atoms(self) -> int:
        """
        Returns the number of heavy atoms in the protein.

        Returns:
            int: The number of heavy atoms in the protein.
        """
        count = 0
        for chain in self._chains.values():
            count += chain.num_heavy_atoms
        return count

    @property
    def num_hydrogen_atoms(self) -> int:
        """
        Returns the number of hydrogen atoms in the protein.

        Returns:
            int: The number of hydrogen atoms in the protein.
        """
        count = 0
        for chain in self._chains.values():
            count += chain.num_hydrogen_atoms
        return count

    @property
    def num_disordered_atoms(self) -> int:
        """
        Returns the number of disordered atoms in the protein.

        Returns:
            int: The number of disordered atoms in the protein.
        """
        count = 0
        for chain in self._chains.values():
            count += chain.num_disordered_atoms
        return count

    @property
    def num_hetero_atoms(self) -> int:
        """
        Returns the number of hetero atoms in the protein.

        Returns:
            int: The number of hetero atoms in the protein.
        """
        count = 0
        for chain in self._chains.values():
            count += chain.num_hetero_atoms
        return count

    def get_atom(self, chain_id: str, residue_index: int, atom_name: str) -> Optional[Atom]:
        """
        Returns an atom from the protein if it exists.

        Args:
            chain_id (str): The ID of the chain.
            residue_index (int): The index of the residue.
            atom_name (str): The name of the atom.

        Returns:
            Atom: The atom if it exists, otherwise None.
        """
        if chain_id in self._chains:
            return self._chains[chain_id].get_residue(residue_index).get_atom(atom_name)
        else:
            return None

    def filter_atoms(self,
                     chain_criteria: Optional[List] = None,
                     atom_criteria: Optional[List] = None,
                     residue_criteria: Optional[List] = None) -> Iterator[Atom]:
        """
        Filters atoms in the protein based on the specified criteria.
        Yielded atoms must match all criteria.

        Args:
            chain_criteria (list of tuple): The criteria to filter the chains.
            atom_criteria (list of tuple): The criteria to filter the atoms.
            residue_criteria (list of tuple): The criteria to filter the residues.

        Yields:
            Atom: The next atom in the protein that matches the criteria.
        """
        for chain in self.filter_chains(chain_criteria):
            for residue in chain.filter_residues(residue_criteria):
                for atom in residue.filter_atoms(atom_criteria):
                    yield atom

    def missing_heavy_atom_types(self) -> Dict[str, Set[str]]:
        """
        Returns a dictionary of missing heavy atom types in the protein.

        Returns:
            dict of str, set of str: The missing heavy atom types in the protein.
        """
        missing_atom_types = {}
        for chain in self.chains:
            missing = chain.missing_heavy_atom_types()
            if len(missing) > 0:
                missing_atom_types[chain.chain_id] = missing
        return missing_atom_types

    def extra_heavy_atom_types(self) -> Dict[str, Set[str]]:
        """
        Returns a dictionary of extra heavy atom types in the protein.

        Returns:
            dict of str, set of str: The extra heavy atom types in the protein.
        """
        extra_atom_types = {}
        for chain in self.chains:
            extra = chain.extra_heavy_atom_types()
            if len(extra) > 0:
                extra_atom_types[chain.chain_id] = extra
        return extra_atom_types

    def missing_hydrogen_atom_types(self) -> Dict[str, Set[str]]:
        """
        Returns a dictionary of missing hydrogen atom types in the protein.

        Returns:
            dict of str, set of str: The missing hydrogen atom types in the protein.
        """
        missing_atom_types = {}
        for chain in self.chains:
            missing = chain.missing_hydrogen_atom_types()
            if len(missing) > 0:
                missing_atom_types[chain.chain_id] = missing
        return missing_atom_types

    def extra_hydrogen_atom_types(self) -> Dict[str, Set[str]]:
        """
        Returns a dictionary of extra hydrogen atom types in the protein.

        Returns:
            dict of str, set of str: The extra hydrogen atom types in the protein.
        """
        extra_atom_types = {}
        for chain in self.chains:
            extra = chain.extra_hydrogen_atom_types()
            if len(extra) > 0:
                extra_atom_types[chain.chain_id] = extra
        return extra_atom_types

    def fix_disordered_atoms(self):
        """
        Fixes disordered atoms in the protein.

        Returns:
            None
        """
        for chain in self.chains:
            chain.fix_disordered_atoms()

    def remove_atoms(self, atom_types: List[str]) -> None:
        """
        Removes all atoms of the specified type from the protein.

        Args:
            atom_types (list of str): The atom types to remove.

        Returns:
            None
        """
        for chain in self.chains:
            chain.remove_atoms(atom_types)

    def remove_hydrogen_atoms(self):
        """
        Removes all hydrogen atoms from the protein.

        Returns:
            None
        """
        for residue in self.residues:
            residue.remove_hydrogen_atoms()

    def keep_atoms(self, atom_types: List[str]) -> None:
        """
        Keeps only the specified atoms in the protein.

        Args:
            atom_types (list of str): The atom types to keep.

        Returns:
            None
        """
        for chain in self.chains:
            chain.keep_atoms(atom_types)

    def keep_backbone_atoms(self) -> None:
        """
        Keeps only the backbone atoms in the protein.

        Returns:
            None
        """
        for chain in self.chains:
            chain.keep_backbone_atoms()

    # ------------------------------------------------------------
    # Methods for managing the protein's attributes.
    # - has_attribute
    # - get_attribute
    # - set_attribute
    # - delete_attribute
    # - sum_attribute
    # - assign_list
    # ------------------------------------------------------------

    def has_attribute(self, key: str) -> bool:
        """
        Checks if the protein has the specified attribute.

        Args:
            key (str): The name of the attribute.

        Returns:
            bool: True if the protein has the specified attribute.
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
        if self.has_attribute(key):
            delattr(self, "_" + key)

    def sum_attribute(self, key: str) -> None:
        """
        Sums the value of the specified attribute.
        """
        attribute_sum = 0.0
        for chain in self._chains.values():
            if not hasattr(chain, key):
                chain.sum_attribute(key)
                attribute_sum += chain.get_attribute(key)
        self.set_attribute(key, attribute_sum)

    @staticmethod
    def assign_list(object_list: List[Any], value_list: List[Any], key: str) -> None:
        if len(object_list) != len(value_list):
            raise Exception("List lengths do not match")
        for i in range(len(object_list)):
            object_list[i].set_attribute(key, value_list[i])

    # ------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------

    def summary(self) -> str:
        """
        Returns a summary of the protein.
        """
        text = [
            f"Protein ID: {self.id}",
            f"# Chains: {self.num_chains}",
            f"# Residues: {self.num_residues}",
            f"# Disordered residues: {self.num_disordered_residues}",
            f"# Hetero residues: {self.num_hetero_residues}",
            f"# Water residues: {self.num_water_residues}",
            f"# Atoms: {self.num_atoms}",
            f"# Heavy atoms: {self.num_heavy_atoms}",
            f"# Hydrogen atoms: {self.num_hydrogen_atoms}",
            f"# Disordered atoms: {self.num_disordered_atoms}",
            f"# Hetero atoms: {self.num_hetero_atoms}"
        ]

        return "\n".join(text)
