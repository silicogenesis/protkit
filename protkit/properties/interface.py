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

from typing import List

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
