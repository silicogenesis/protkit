#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `SurfaceArea` to calculate the surface area of proteins.

The ASA or SASA was first described by Lee and Richards in 1971. Shrake and Rupley
developed the rolling ball algorithm in 1973.

RASA can be calculated by dividing the ASA by the MaxASA. MasASA values were calculated
by Miller et al. (1987).

Levy (2010) describes a method for defining interior, surface, support, core and rim
residues of a protein based on ASA calculations.

Abbreviations:
--------------

ASA    - Accessible Surface Area
BASA   - Buried Surface Area
MaxASA - Maximum Accessible Surface Area
RASA   - Relative Accessible Surface Area = ASA / MaxASA
SASA   - Solvent Accessible Surface Area (same as ASA)

Papers referenced in this code:
-------------------------------

Lee B., Richards F.M. (1971)
The interpretation of protein structures: estimation of static accessibility.
J. Mol. Biol. Vol. 55, pp. 279-400.

Levy E.D. (2010)
A simple definition of structural regions in proteins and its use in analyzing
interface evolution.
J. Mol. Biol. Vol. 403, pp. 660-670.

Miller S., Janin J., Lesk A.M., Chothia C. (1987).
Interior and surface of monomeric proteins.
J. Mol. Biol. Vol. 196, pp. 641-656.

Mitternacht S. (2016)
FreeSASA: An open source C library for solvent accessible surface area calculations.
F1000 Research, 5:189.

Shrake A., Rupley J.A. (1973)
Environment and exposure to solvent of protein atoms. Lysozyme and insulin.
J. Mol. Biol. Vol. 79, pp. 351-371.

External packages used in this code:
------------------------------------

FreeSASA - https://freesasa.github.io/
Provides C implementations of the Lee-Richards and Shrake-Rupley algorithms
with implementation hooks for Python.

"""

from typing import TYPE_CHECKING
from protkit.structure import Protein
from protkit.tools.freesasa_adaptor import FreeSASAAdaptor


class SurfaceArea:
    MAX_ASA_MILLER = {
        "ALA": 113.0,
        "ARG": 241.0,
        "ASN": 158.0,
        "ASP": 151.0,
        "CYS": 140.0,
        "GLN": 189.0,
        "GLU": 183.0,
        "GLY": 85.0,
        "HIS": 194.0,
        "ILE": 182.0,
        "LEU": 180.0,
        "LYS": 211.0,
        "MET": 204.0,
        "PHE": 218.0,
        "PRO": 143.0,
        "SER": 122.0,
        "THR": 146.0,
        "TRP": 259.0,
        "TYR": 229.0,
        "VAL": 160.0
    }

    @staticmethod
    def surface_area_of_protein(protein: Protein,
                                assign_attribute: bool = False,
                                key: str = "surface_area") -> float:
        """
        Calculate the surface area of a protein.

        Args:
            protein (Protein): The protein object for which to calculate the surface area.
            assign_attribute (bool): Whether to assign the surface area to the protein object.
            key (str): The attribute key to use when assigning the surface area to the protein object.

        Returns:
            float: The surface area of the protein.
        """

        # Calculate the solvent accessible surface area using Lee Richards
        freesasa_adaptor = FreeSASAAdaptor(algorithm=FreeSASAAdaptor.LEE_RICHARDS)
        sasa = freesasa_adaptor.calculate_surface_area(list(protein.atoms))

        if assign_attribute:
            index = 0
            protein_surface_area = 0
            for chain in protein.chains:
                chain_surface_area = 0
                for residue in chain.residues:
                    residue_surface_area = 0
                    for atom in residue.atoms:
                        atom.set_attribute(key, sasa[index])
                        residue_surface_area += sasa[index]
                        index += 1
                    chain_surface_area += residue_surface_area
                    residue.set_attribute(key, residue_surface_area)
                protein_surface_area += chain_surface_area
                chain.set_attribute(key, chain_surface_area)
            protein.set_attribute(key, protein_surface_area)
            return protein_surface_area
        else:
            return sum(sasa)
