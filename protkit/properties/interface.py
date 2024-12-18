#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `Interface` to represent an interface between two proteins.

An interface is defined as the set of atoms or residues within one chain or protein
that are within a specified distance of atoms or residues in another chain or protein.

The class is meant to be used as a utility class and is not meant to be instantiated.

It uses a SpaceQuery data structure to calculate interacting atoms and residues fast.
"""

from typing import List, Tuple


from protkit.structure.residue import Residue
from protkit.structure.atom import Atom

from protkit.geometry.space_query import SpaceQuery


class Interface:
    @staticmethod
    def interface_atoms(atoms1: List[Atom],
                        atoms2: List[Atom],
                        cutoff: float = 5.0,
                        assign_attribute: bool = False,
                        key: str = "in_interface") -> (List[Atom], List[Atom]):
        """
        Returns a list of atoms that are within a specified distance of each other.

        Args:
            atoms1 (List[Atom]): A list of atoms.
            atoms2 (List[Atom]): A list of atoms.
            cutoff (float): The cutoff distance.
            assign_attribute (bool): If True, the atoms will be assigned as an attribute to the atoms.
            key (str): The key to use for the attribute.

        Returns:
            List[Atom]: A list of atoms that are within the cutoff distance of each other.
        """

        # Get the coordinates of the atoms
        coordinates1 = [(atom.x, atom.y, atom.z) for atom in atoms1]
        coordinates2 = [(atom.x, atom.y, atom.z) for atom in atoms2]

        # Create a KDTree from the coordinates
        # Get the atoms that are within the cutoff distance of each other
        tree = SpaceQuery(coordinates1)
        indices1, indices2 = tree.query_partners(coordinates2, cutoff)
        interface_atoms1 = [atoms1[i] for i in indices1]
        interface_atoms2 = [atoms2[i] for i in indices2]

        # Assign the atoms as an attribute
        if assign_attribute:
            for atom in interface_atoms1:
                atom.set_attribute(key, True)
            for atom in interface_atoms2:
                atom.set_attribute(key, True)

        return interface_atoms1, interface_atoms2

    @staticmethod
    def interface_residues(residues1: List[Residue],
                           residues2: List[Residue],
                           cutoff: float = 5.0,
                           assign_attribute: bool = False,
                           key: str = "in_interface") -> (List[Residue], List[Residue]):
        """
        Returns a list of residues that are within a specified distance of each other. Two residues
        are considered to be within the cutoff distance of each other if any of their atoms are
        within the cutoff distance of each other.

        Args:
            residues1 (List[Residue]): A list of residues.
            residues2 (List[Residue]): A list of residues.
            cutoff (float): The cutoff distance.
            assign_attribute (bool): If True, the residues will be assigned as an attribute to the residues.
            key (str): The key to use for the attribute.

        Returns:
            List[Residue]: A list of residues that are within the cutoff distance of each other.
        """

        # Get all atoms in the residues
        atoms1 = [atom for residue in residues1 for atom in residue.atoms]
        atoms2 = [atom for residue in residues2 for atom in residue.atoms]

        # Get the atoms that are within the cutoff distance of each other
        interface_atoms1, interface_atoms2 = Interface.interface_atoms(atoms1, atoms2, cutoff, assign_attribute=assign_attribute, key=key)

        # Get the residues that contain the interface atoms
        interface_residues1 = []
        interface_residues2 = []
        last_residue = None
        for atom1 in interface_atoms1:
            if atom1.residue != last_residue:
                interface_residues1.append(atom1.residue)
                last_residue = atom1.residue
        last_residue = None
        for atom2 in interface_atoms2:
            if atom2.residue != last_residue:
                interface_residues2.append(atom2.residue)
                last_residue = atom2.residue

        # Get the residues that contain the interface atoms
        if assign_attribute:
            for residue in interface_residues1:
                residue.set_attribute(key, True)
            for residue in interface_residues2:
                residue.set_attribute(key, True)

        return interface_residues1, interface_residues2

    @staticmethod
    def interface_residues_from_alpha_carbon(residues1: List[Residue],
                                             residues2: List[Residue],
                                             cutoff: float = 5.0,
                                             assign_attribute: bool = False,
                                             key: str = "in_interface") -> (List[Residue], List[Residue]):
        """
        Returns a list of residues that are within a specified distance of each other. Two residues
        are considered to be within the cutoff distance of each other if any of their alpha carbon
        atoms are within the cutoff distance of each other.

        Args:
            residues1 (List[Residue]): A list of residues.
            residues2 (List[Residue]): A list of residues.
            cutoff (float): The cutoff distance.
            assign_attribute (bool): If True, the residues will be assigned as an attribute to the residues.
            key (str): The key to use for the attribute.

        Returns:
            List[Residue]: A list of residues that are within the cutoff distance of each other.
        """

        # Get all alpha carbon atoms in the residues
        atoms1 = [residue.get_atom("CA") for residue in residues1]
        atoms2 = [residue.get_atom("CA") for residue in residues2]

        # Get the coordinates of the atoms
        coordinates1 = [(atom.x, atom.y, atom.z) for atom in atoms1]
        coordinates2 = [(atom.x, atom.y, atom.z) for atom in atoms2]

        # Create a KDTree from the coordinates
        # Get the atoms that are within the cutoff distance of each other
        tree = SpaceQuery(coordinates1)
        indices1, indices2 = tree.query_partners(coordinates2, cutoff)
        interface_residues1 = [residues1[i] for i in indices1]
        interface_residues2 = [residues2[i] for i in indices2]

        # Assign the interface to residues as an attribute
        if assign_attribute:
            for residue in interface_residues1:
                residue.set_attribute(key, True)
            for residue in interface_residues2:
                residue.set_attribute(key, True)

        return interface_residues1, interface_residues2

    @staticmethod
    def interface_atom_pairs(atoms1: List[Atom],
                             atoms2: List[Atom],
                             cutoff: float = 5.0,
                             assign_attribute: bool = False,
                             key: str = "interacting_atoms") -> List[Tuple[Atom, Atom]]:
        """
        Returns a list of atom pairs (a1, a2) that are within a specified distance of each other.

        Two atoms are considered to be interacting if the distance between them
        is less than the specified cutoff.

        If assign_attribute is True, each atom will have the `key` attribute
        updated/created as a list of interacting partner atoms.

        Args:
            atoms1 (List[Atom]): A list of atoms (e.g., from one protein or selection).
            atoms2 (List[Atom]): Another list of atoms (e.g., from another protein or selection).
            cutoff (float): The distance cutoff for considering atoms to be interacting.
            assign_attribute (bool): If True, adds the interacting partner atoms as attributes to each atom.
            key (str): Attribute name to store the list of partner atoms.

        Returns:
            List[(Atom, Atom)]: A list of tuples, each containing an atom from atoms1 and an atom
                                 from atoms2 that are interacting.
        """

        # Extract coordinates from atoms2 for KD-tree construction
        atoms2_coords = [(a.x, a.y, a.z) for a in atoms2]

        # Create KD-tree from atoms2
        tree = SpaceQuery(atoms2_coords)

        interacting_pairs = []

        # For each atom in atoms1, query the KD-tree
        for a1 in atoms1:
            coords1 = [(a1.x, a1.y, a1.z)]
            idx1, idx2 = tree.query_partners(coords1, cutoff)
            # idx1 corresponds to matches in atoms2_coords
            for partner_index in idx1:
                a2 = atoms2[partner_index]
                interacting_pairs.append((a1, a2))

        # Assign attributes if requested
        if assign_attribute:
            for (atom_a, atom_b) in interacting_pairs:
                # Update atom_a's attribute
                if atom_a.has_attribute(key):
                    partners = atom_a.get_attribute(key)
                else:
                    partners = []
                partners.append(atom_b)
                atom_a.set_attribute(key, partners)

                # Update atom_b's attribute
                if atom_b.has_attribute(key):
                    partners = atom_b.get_attribute(key)
                else:
                    partners = []
                partners.append(atom_a)
                atom_b.set_attribute(key, partners)

        return interacting_pairs

    @staticmethod
    def interface_residue_pairs(residues1: List[Residue],
                                residues2: List[Residue],
                                cutoff: float = 5.0,
                                assign_attribute: bool = False,
                                key: str = "interacting_residues") -> List[Tuple[Residue, Residue]]:
        """
        Returns a list of residue pairs (r1, r2) that are within a specified distance of each other.

        Two residues are considered interacting if any atom of a residue in `residues1` is
        within `cutoff` distance of any atom of a residue in `residues2`.

        If assign_attribute is True, each residue will have the `key` attribute
        updated/created as a list of interacting partner residues.

        Args:
            residues1 (List[Residue]): A list of residues (e.g., from one protein or selection).
            residues2 (List[Residue]): Another list of residues (e.g., from another protein or selection).
            cutoff (float): The distance cutoff for considering residues to be interacting.
            assign_attribute (bool): If True, adds the interacting partner residues as attributes to each residue.
            key (str): Attribute name to store the list of partner residues.

        Returns:
            List[(Residue, Residue)]: A list of tuples, each containing a residue from residues1 and a residue
                                      from residues2 that are interacting.
        """

        # Flatten atoms from residues2 and keep a mapping to their residues
        atoms2_coords = []
        atoms2_res_map = []
        for r2 in residues2:
            for a2 in r2.atoms:
                atoms2_coords.append((a2.x, a2.y, a2.z))
                atoms2_res_map.append(r2)

        # Create a KD-tree for the second set of residues
        tree = SpaceQuery(atoms2_coords)

        interacting_pairs = []

        # For each residue in residues1, find interacting residues in residues2
        for r1 in residues1:
            found_partner_residues = set()

            for a1 in r1.atoms:
                # Query the KD-tree with the coordinates of a single atom from residue1
                idx1, idx2 = tree.query_partners([(a1.x, a1.y, a1.z)], cutoff)
                # idx1 are indices in atoms2_coords
                for partner_index in idx1:
                    found_partner_residues.add(atoms2_res_map[partner_index])

            for r2 in found_partner_residues:
                interacting_pairs.append((r1, r2))

        # Assign attributes if requested
        if assign_attribute:
            for (res_a, res_b) in interacting_pairs:
                # Update res_a's attribute
                if res_a.has_attribute(key):
                    partners = res_a.get_attribute(key)
                else:
                    partners = []
                partners.append(res_b)
                res_a.set_attribute(key, partners)

                # Update res_b's attribute
                if res_b.has_attribute(key):
                    partners = res_b.get_attribute(key)
                else:
                    partners = []
                partners.append(res_a)
                res_b.set_attribute(key, partners)

        return interacting_pairs