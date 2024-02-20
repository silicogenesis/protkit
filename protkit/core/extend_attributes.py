#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Class `ExtendedAttributes` defines a set of methods that is used
to manage the properties for an object. These methods are available
to the object by inheritance.

Methods include:

- has_attribute(key): Checks if the object has the specified attribute.
- get_attribute(key): Returns the value of the specified attribute.
- set_attribute(key, value): Sets the value of the specified attribute.
- delete_attribute(key): Deletes the specified attribute.

Classes such as Protein, Chain, Residue, Atom, Sequence, etc.
inherit these methods to provide a consistent interface for
managing object properties.
"""

from typing import Any, List


class ExtendedAttributes:
    def __init__(self):
        """
        Constructor.
        """
        self._extended_attributes = set()

    def has_attribute(self, key: str) -> bool:
        """
        Checks if the object has the specified attribute.

        Args:
            key (str): The name of the attribute.

        Returns:
            bool: True if the object has the specified attribute.
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
            self._extended_attributes.add(key)

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
            self._extended_attributes.remove(key)

    def list_attributes(self) -> List[str]:
        """
        Lists the attributes that have been defined for the object.

        Returns:
            List[str]: A list of attribute names.
        """
        return list(self._extended_attributes)