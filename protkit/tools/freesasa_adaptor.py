#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS), Claudio Jardim (CJ)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `FreeSASAAdaptor` to calculate the accessible surface area (ASA)
of proteins using the FreeSASA library.
"""

import freesasa
from typing import List
import numpy as np

from protkit.structure import Atom
from protkit.tasks.surface_area_calculator import SurfaceAreaCalculator
from protkit.properties.vdw_radius import VDWRadius


class FreeSASAAdaptor(SurfaceAreaCalculator):
    """
    Adapter class for the FreeSASA tool.

    This class provides concrete implementations of surface area calculations
    by interfacing with the FreeSASA library and updating Protein objects with
    the calculated surface properties.
    """
    LEE_RICHARDS = freesasa.LeeRichards
    SHRAKE_RUPLEY = freesasa.ShrakeRupley

    DEFAULT_PROBE_RADIUS = 1.4
    DEFAULT_LEE_RICHARDS_SLICES = 20
    DEFAULT_SHRAKE_RUPLEY_POINTS = 100
    DEFAULT_N_THREADS = 1

    def __init__(self,
                 algorithm: str = LEE_RICHARDS,
                 probe_radius: float = DEFAULT_PROBE_RADIUS,
                 n_points: int = DEFAULT_SHRAKE_RUPLEY_POINTS,
                 n_slices: int = DEFAULT_LEE_RICHARDS_SLICES,
                 n_threads: int = DEFAULT_N_THREADS):
        """
        Initialize the FreeSASA adaptor with the desired settings.

        Args:
            algorithm (str): The algorithm to use for surface area calculations.
                Options: LEE_RICHARDS, SHRAKE_RUPLEY
            probe_radius (float): The radius of the probe to use for surface area calculations.
            n_points (int): The number of points to use for the surface area calculations.
            n_slices (int): The number of slices to use for the surface area calculations.
            n_threads (int): The number of threads to use for the surface area calculations.
        """
        super().__init__()

        self._algorithm = algorithm
        self._probe_radius = probe_radius
        self._n_points = n_points
        self._n_slices = n_slices
        self._n_threads = n_threads

    def calculate_surface_area(self, atoms: List[Atom]) -> List[float]:
        """
        Calculate the accessible surface area of a protein using FreeSASA.

        Args:
            atoms (List[Atom]): A list of Atom objects whose surface area needs to be calculated.

        Returns:
            List[float]: A list of surface areas for each atom in the input list.
        """

        # Extract the atom coordinates and van der Waals radii
        num_atoms = len(atoms)
        atom_coordinates = np.array([(atom.x, atom.y, atom.z) for atom in atoms])
        vdw_radii = [VDWRadius.SURFACE_RADIUS.get(atom.element, 2.0) for atom in atoms]

        # Calculate the surface area of the protein
        parameters = freesasa.Parameters({
            'algorithm': self._algorithm,
            'probe-radius': self._probe_radius,
            'n-points': self._n_points,
            'n-slices': self._n_slices,
            'n-threads': self._n_threads
        })

        # Calculate the solvent accessible surface area
        result = freesasa.calcCoord(atom_coordinates.flatten(), vdw_radii, parameters=parameters)
        atom_sasa = [result.atomArea(i) for i in range(num_atoms)]

        # assert (sum(atom_sasa) == result.totalArea())
        # for i in range(len(vdw_radii)):
        #     atom_sasa[i] = result.atomArea(i)
        # atom_sasa_total = result.totalArea()

        # if self._assign_to_protein:
        #     protein.assign_list(list(protein.atoms), atom_sasa, "sasa")
        #     protein.set_attribute("sasa_total", atom_sasa_total)
        #
        #     for chain in protein.chains:
        #         for residue in chain.residues:
        #             residue.sum_attribute("sasa")
        #         chain.sum_attribute("sasa")

        return atom_sasa




