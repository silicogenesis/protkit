#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `Atom` to represent an atom in a protein.
"""

from __future__ import annotations
from typing import TYPE_CHECKING, Optional, Any
from copy import deepcopy

if TYPE_CHECKING:
    from protkit.structure.residue import Residue


class Atom:
    # ------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------
    def __init__(self,
                 # Core
                 element: str,
                 atom_type: str,
                 x: float,
                 y: float,
                 z: float,

                 is_hetero: bool = False,
                 is_disordered: bool = False,

                 # Optional
                 alt_loc: Optional[str] = None,
                 occupancy: Optional[float] = None,
                 temp_factor: Optional[float] = None,
                 assigned_charge: Optional[str] = None,

                 residue: Optional[Residue] = None):
        """
        Constructor for the Atom class.

        Args:
            element (str): The element associated with the atom.
            atom_type (str): The atom type.
            x (float): The x coordinate of the atom.
            y (float): The y coordinate of the atom.
            z (float): The z coordinate of the atom.
            is_hetero (bool, optional): The hetero flag. Defaults to False.
            is_disordered (bool, optional): The disordered flag. Defaults to False.
            alt_loc (Optional[str], optional): The alt_loc identifier. Defaults to None.
            occupancy (Optional[float], optional): The occupancy. Defaults to None.
            temp_factor (Optional[float], optional): The temperature factor. Defaults to None.
            assigned_charge (Optional[str], optional): The assigned charge. Defaults to None.
            residue (Optional[Residue], optional): The residue the atom forms part of. Defaults to None.

        Returns:
            None
        """

        # Core properties
        self._x: float = x
        self._y: float = y
        self._z: float = z
        self._element: str = element
        self._atom_type: str = atom_type

        # Derived properties
        self._is_hetero: bool = is_hetero
        self._is_disordered: bool = is_disordered

        # Optional attributes
        if alt_loc is not None and alt_loc != "":
            self.set_attribute("alt_loc", alt_loc)
        if occupancy is not None and occupancy != 1.0:
            self.set_attribute("occupancy", occupancy)
        if temp_factor is not None:
            self.set_attribute("temp_factor", temp_factor)
        if assigned_charge is not None and assigned_charge != "":
            self.set_attribute("assigned_charge", assigned_charge)

        # If alt_loc is specified, the atom is disordered
        # Keep track of all disordered states.
        if self.get_attribute("alt_loc") is not None:
            self._is_disordered = True
            self._disordered_states = [{
                "alt_loc": self.get_attribute("alt_loc"),
                "occupancy": self.get_attribute("occupancy"),
                "temp_factor": self.get_attribute("temp_factor"),
                "x": self._x,
                "y": self._y,
                "z": self._z
            }]

        self._residue: Optional[Residue] = residue

    def copy(self):
        """
        Returns a deep copy of the atom, excluding the residue.

        Returns:
            Atom: A deep copy of the atom, excluding the residue.
        """
        residue = self._residue
        self._residue = None
        atom = deepcopy(self)
        self._residue = residue
        return atom

    # ------------------------------------------------------------
    # Methods for managing the atom's core properties.
    # - id (property)
    # - element (property)
    # - element (setter)
    # - atom_type (property)
    # - atom_type (setter)
    # - x (property)
    # - x (setter)
    # - y (property)
    # - y (setter)
    # - z (property)
    # - z (setter)
    # - is_hetero (property)
    # - is_disordered (property)
    # ------------------------------------------------------------

    @property
    def id(self) -> str:
        """
        Returns the atom's id, which is a combination of the
        residue id and the atom type.

        Returns:
            str: The atom's id.
        """
        if self._residue is None:
            return self._atom_type
        else:
            return self._residue.id + ":" + self._atom_type

    @property
    def element(self) -> str:
        """
        Returns the element associated with the atom.

        Returns:
            str: The element associated with the atom.
        """
        return self._element

    @element.setter
    def element(self, element: str) -> None:
        """
        Sets the element associated with the atom.

        Args:
            element (str): The element associated with the atom.

        Returns:
            None
        """
        self._element = element

    @property
    def atom_type(self) -> str:
        """
        Returns the atom type.

        Returns:
            str: The atom type.
        """
        return self._atom_type

    @atom_type.setter
    def atom_type(self, atom_type: str) -> None:
        """
        Sets the atom type.

        Args:
            atom_type (str): The atom type.

        Returns:
            None
        """
        self._atom_type = atom_type

    @property
    def x(self) -> float:
        """
        Returns the x coordinate of the atom.

        Returns:
            float: The x coordinate of the atom.
        """
        return self._x

    @x.setter
    def x(self, x: float) -> None:
        """
        Sets the x coordinate of the atom.

        Args:
            x (float): The x coordinate of the atom.

        Returns:
            None
        """
        self._x = x

    @property
    def y(self) -> float:
        """
        Returns the y coordinate of the atom.

        Returns:
            float: The y coordinate of the atom.
        """
        return self._y

    @y.setter
    def y(self, y: float) -> None:
        """
        Sets the y coordinate of the atom.

        Args:
            y (float): The y coordinate of the atom.

        Returns:
            None
        """
        self._y = y

    @property
    def z(self) -> float:
        """
        Returns the z coordinate of the atom.

        Returns:
            float: The z coordinate of the atom.
        """
        return self._z

    @z.setter
    def z(self, z: float) -> None:
        """
        Sets the z coordinate of the atom.

        Args:
            z (float): The z coordinate of the atom.

        Returns:
            None
        """
        self._z = z

    @property
    def is_hetero(self) -> bool:
        """
        Returns the hetero flag.

        Returns:
            bool: The hetero flag.
        """
        return self._is_hetero

    @property
    def is_disordered(self) -> bool:
        """
        Returns the disordered flag.

        Returns:
            bool: The disordered flag.

        """
        return self._is_disordered

    # ------------------------------------------------------------
    # Methods for managing the atom's residue
    # - residue (property)
    # - residue (setter)
    # ------------------------------------------------------------

    @property
    def residue(self) -> Optional[Residue]:
        """
        Returns the residue the atom forms part of.

        Returns:
            Optional[Residue]: The residue the atom forms part of.
        """
        return self._residue

    @residue.setter
    def residue(self, residue: Residue) -> None:
        """
        Sets the residue the atom forms part of.

        Args:
            residue (Residue): The residue the atom forms part of.

        Returns:
            None
        """
        self._residue = residue

    # ------------------------------------------------------------
    # Disordered state
    # - merge_disordered_atom
    # - fix_disordered_atom
    # ------------------------------------------------------------

    def merge_disordered_atom(self, atom: Atom) -> None:
        """
        Merges the disordered atom with the current atom.
        Keep track of the additional disordered information
        from the atom that the merge is done with.

        Args:
            atom (Atom): The atom to merge with.

        Returns:
            None
        """
        self._disordered_states.append({
            "alt_loc": atom.get_attribute("alt_loc"),
            "occupancy": atom.get_attribute("occupancy"),
            "temp_factor": atom.get_attribute("temp_factor"),
            "x": atom.x,
            "y": atom.y,
            "z": atom.z
        })
        self._is_disordered = True

    def fix_disordered_atom(self) -> None:
        """
        Fixes the disordered atom by selecting the disordered
        state with the highest occupancy.

        Returns:
            None
        """
        if self._is_disordered:
            # Sort disordered states by occupancy
            self._disordered_states.sort(key=lambda dis_state: dis_state["occupancy"], reverse=True)

            # Select the disordered state with the highest occupancy
            state = self._disordered_states[0]

            # Set the atom's properties to the disordered state with the highest occupancy
            self.set_attribute("alt_loc", state["alt_loc"])
            self.set_attribute("occupancy", state["occupancy"])
            self.set_attribute("temp_factor", state["temp_factor"])
            self._x = state["x"]
            self._y = state["y"]
            self._z = state["z"]

            # Remove the disordered states
            self._disordered_states = None
            self._is_disordered = False

    # ------------------------------------------------------------
    # Methods for managing the atom's attributes.
    # - has_attribute
    # - get_attribute
    # - set_attribute
    # - delete_attribute
    # ------------------------------------------------------------

    def has_attribute(self, key: str) -> bool:
        """
        Returns True if the atom has the specified attribute.

        Args:
            key (str): The name of the attribute.

        Returns:
            bool: True if the atom has the specified attribute.
        """
        return hasattr(self, "_" + key)

    def get_attribute(self, key: str) -> Any:
        """
        current_residueReturns the value of the specified attribute.

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
        """"
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
