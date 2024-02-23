#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `Math` to perform various mathematical operations. These mathematical
operations include vector operations, such as calculating the magnitude, dot product,
cross product, Euclidean distance, angle and dihedral angle between vectors.

The class is meant to be used as a utility class and is not meant to be instantiated.
It is unaware of the context in which it is used and does not have any dependencies
on other classes or modules in the protkit package.
"""

import math
from typing import Tuple


class Math:
    @staticmethod
    def magnitude(x: float, y: float, z: float) -> float:
        """
        Calculates the magnitude of a vector.

        Args:
            x (float): The x-coordinate of the vector.
            y (float): The y-coordinate of the vector.
            z (float): The z-coordinate of the vector.

        Returns:
            float: The magnitude of the vector.
        """
        return math.sqrt(x * x + y * y + z * z)

    @staticmethod
    def dot_product(x1: float, y1: float, z1: float, x2: float, y2: float, z2: float) -> float:
        """
        Calculates the dot product of two vectors.

        Args:
            x1 (float): The x-coordinate of the first vector.
            y1 (float): The y-coordinate of the first vector.
            z1 (float): The z-coordinate of the first vector.
            x2 (float): The x-coordinate of the second vector.
            y2 (float): The y-coordinate of the second vector.
            z2 (float): The z-coordinate of the second vector.

        Returns:
            float: The dot product of the two vectors.
        """
        return x1 * x2 + y1 * y2 + z1 * z2

    @staticmethod
    def cross_product(x1: float, y1: float, z1: float, x2: float, y2: float, z2: float) -> Tuple[float, float, float]:
        """
        Calculates the cross product of two vectors.

        Args:
            x1 (float): The x-coordinate of the first vector.
            y1 (float): The y-coordinate of the first vector.
            z1 (float): The z-coordinate of the first vector.
            x2 (float): The x-coordinate of the second vector.
            y2 (float): The y-coordinate of the second vector.
            z2 (float): The z-coordinate of the second vector.

        Returns:
            Tuple[float, float, float]: The cross product of the two vectors.
        """
        return y1 * z2 - z1 * y2, z1 * x2 - x1 * z2, x1 * y2 - y1 * x2

    @staticmethod
    def euclidean_distance(x1: float, y1: float, z1: float, x2: float, y2: float, z2: float) -> float:
        """
        Calculates the Euclidean distance between two points.

        Args:
            x1 (float): The x-coordinate of the first point.
            y1 (float): The y-coordinate of the first point.
            z1 (float): The z-coordinate of the first point.
            x2 (float): The x-coordinate of the second point.
            y2 (float): The y-coordinate of the second point.
            z2 (float): The z-coordinate of the second point.

        Returns:
            float: The distance between the two points.
        """
        dx = x1 - x2
        dy = y1 - y2
        dz = z1 - z2
        return math.sqrt(dx * dx + dy * dy + dz * dz)

    @staticmethod
    def vector_angle(x1: float, y1: float, z1: float,
                     x2: float, y2: float, z2: float) -> float:
        """
        Calculates the angle between two vectors.

        Args:
            x1 (float): The x-coordinate of the first vector.
            y1 (float): The y-coordinate of the first vector.
            z1 (float): The z-coordinate of the first vector.
            x2 (float): The x-coordinate of the second vector.
            y2 (float): The y-coordinate of the second vector.
            z2 (float): The z-coordinate of the second vector.

        Returns:
            float: The angle in degrees between the two vectors.
        """
        # Calculate the dot product of the vectors
        dot_product = Math.dot_product(x1, y1, z1, x2, y2, z2)

        # Calculate the magnitude of the vectors
        mag_u = Math.magnitude(x1, y1, z1)
        mag_v = Math.magnitude(x2, y2, z2)

        # Calculate the cosine of the angle
        cos_angle = dot_product / (mag_u * mag_v)

        # Edge cases
        if cos_angle > 1:
            cos_angle = 1
        elif cos_angle < -1:
            cos_angle = -1

        angle_radians = math.acos(cos_angle)
        angle = math.degrees(angle_radians)

        return angle

    @staticmethod
    def angle(x1: float, y1: float, z1: float,
              x2: float, y2: float, z2: float,
              x3: float, y3: float, z3: float) -> float:
        """
        Calculates the angle between three points.

        Args:
            x1 (float): The x-coordinate of the first point.
            y1 (float): The y-coordinate of the first point.
            z1 (float): The z-coordinate of the first point.
            x2 (float): The x-coordinate of the second point.
            y2 (float): The y-coordinate of the second point.
            z2 (float): The z-coordinate of the second point.
            x3 (float): The x-coordinate of the third point.
            y3 (float): The y-coordinate of the third point.
            z3 (float): The z-coordinate of the third point.

        Returns:
            float: The angle in degrees between the three points.
        """

        # Calculate the vectors between the points
        ux = x1 - x2
        uy = y1 - y2
        uz = z1 - z2
        vx = x3 - x2
        vy = y3 - y2
        vz = z3 - z2

        return Math.vector_angle(ux, uy, uz, vx, vy, vz)

    @staticmethod
    def dihedral_angle(x1: float, y1: float, z1: float,
                       x2: float, y2: float, z2: float,
                       x3: float, y3: float, z3: float,
                       x4: float, y4: float, z4: float) -> float:
        """
        Calculates the dihedral angle between four points.

        Args:
            x1 (float): The x-coordinate of the first point.
            y1 (float): The y-coordinate of the first point.
            z1 (float): The z-coordinate of the first point.
            x2 (float): The x-coordinate of the second point.
            y2 (float): The y-coordinate of the second point.
            z2 (float): The z-coordinate of the second point.
            x3 (float): The x-coordinate of the third point.
            y3 (float): The y-coordinate of the third point.
            z3 (float): The z-coordinate of the third point.
            x4 (float): The x-coordinate of the fourth point.
            y4 (float): The y-coordinate of the fourth point.
            z4 (float): The z-coordinate of the fourth point.

        Returns:
            float: The dihedral angle in degrees between the four points.
        """

        # Calculate the vectors between the points
        ux, uy, uz = x2 - x1, y2 - y1, z2 - z1
        vx, vy, vz = x3 - x2, y3 - y2, z3 - z2
        wx, wy, wz = x4 - x3, y4 - y3, z4 - z3

        # Calculate the normal vectors
        n1x, n1y, n1z = Math.cross_product(ux, uy, uz, vx, vy, vz)
        n2x, n2y, n2z = Math.cross_product(vx, vy, vz, wx, wy, wz)

        # Calculate the dot products
        angle_deg = Math.vector_angle(n1x, n1y, n1z, n2x, n2y, n2z)

        # Calculate the sign of the dihedral angle
        v_cross_x, vcross_y, vcross_z = Math.cross_product(n1x, n1y, n1z, n2x, n2y, n2z)
        dot_uvwx = Math.dot_product(v_cross_x, vcross_y, vcross_z, vx, vy, vz)
        if dot_uvwx < 0:
            angle_deg = -angle_deg

        return angle_deg
