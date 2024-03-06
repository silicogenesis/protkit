#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3


"""
Implements class `Mass` to calculate the mass of an atom, residue, chain, protein or sequence.

Masses can be calculated in three ways:

- Atomic weight. In the first method, a mass is assigned to each atom.
  It is then propagated to the residue, chain and protein levels.
  The first method is applicable if the complete atomic structure is available.
- Residue weight. In the second method, a mass is assigned to each residue.
  It is then propagated to the chain and protein. It can also be assigned to a sequence,
  since it is dependent on knowledge of the residue sequence only.
- Molecular weight. In the third method, the molecular mass is assigned to each residue.
  This is the mass of the unbound amino acid. It would not make sense to
  assign this mass to the chain or protein. It can be assigned to a sequence.

Masses are in Daltons. Calculated values can be added as attributes to the respective objects.

For more information, see:
https://en.wikipedia.org/wiki/Atomic_mass_unit

The Dalton values were obtained from:
https://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html
"""

from protkit.structure.atom import Atom
from protkit.structure.residue import Residue
from protkit.structure.chain import Chain
from protkit.structure.protein import Protein
from protkit.seq.sequence import Sequence


class Mass:
    # Molecular weights of residues in Daltons.
    # Note that these are the weights of amino acids,
    # not the residues when bound in a chain.
    MOLECULAR_MASS = {
        "ALA": 89.10,
        "ARG": 174.20,
        "ASN": 132.12,
        "ASP": 133.11,
        "CYS": 121.16,
        "GLN": 146.15,
        "GLU": 147.13,
        "GLY": 75.07,
        "HIS": 155.16,
        "ILE": 131.18,
        "LEU": 131.18,
        "LYS": 146.19,
        "MET": 149.21,
        "PHE": 165.19,
        "PRO": 115.13,
        "SER": 105.09,
        "THR": 119.12,
        "TRP": 204.23,
        "TYR": 181.19,
        "VAL": 117.15,

        "A": 89.10,
        "R": 174.20,
        "N": 132.12,
        "D": 133.11,
        "C": 121.16,
        "Q": 146.15,
        "E": 147.13,
        "G": 75.07,
        "H": 155.16,
        "I": 131.18,
        "L": 131.18,
        "K": 146.19,
        "M": 149.21,
        "F": 165.19,
        "P": 115.13,
        "S": 105.09,
        "T": 119.12,
        "W": 204.23,
        "Y": 181.19,
        "V": 117.15
    }

    # Mass values of residues in Daltons.
    # These are the weights of residues when bound in a chain.
    # It excludes the mass of water that is released during
    # the formation of the peptide bond.
    MASS_RESIDUE = {
        "ALA": 71.0788,
        "ARG": 156.1875,
        "ASN": 114.1038,
        "ASP": 115.0886,
        "CYS": 103.1388,
        "GLN": 128.1307,
        "GLU": 129.1155,
        "GLY": 57.0519,
        "HIS": 137.1411,
        "ILE": 113.1594,
        "LEU": 113.1594,
        "LYS": 128.1741,
        "MET": 131.1926,
        "PHE": 147.1766,
        "PRO": 97.1167,
        "SER": 87.0782,
        "THR": 101.1051,
        "TRP": 186.2132,
        "TYR": 163.1760,
        "VAL": 99.1326,

        "A": 71.0788,
        "R": 156.1875,
        "N": 114.1038,
        "D": 115.0886,
        "C": 103.1388,
        "Q": 128.1307,
        "E": 129.1155,
        "G": 57.0519,
        "H": 137.1411,
        "I": 113.1594,
        "L": 113.1594,
        "K": 128.1741,
        "M": 131.1926,
        "F": 147.1766,
        "P": 97.1167,
        "S": 87.0782,
        "T": 101.1051,
        "W": 186.2132,
        "Y": 163.1760,
        "V": 99.1326
    }

    # Mass values of elements in Daltons.
    MASS_ATOM = {
        "H": 1.00794,
        "C": 12.0107,
        "N": 14.0067,
        "O": 15.9994,
        "S": 32.065,
        "Se": 78.96
    }

    @staticmethod
    def residue_mass_of_residue(residue: Residue,
                                assign_attribute: bool = False,
                                key: str = "residue_mass") -> float:
        """
        Returns the mass of the residue.

        Args:
            residue (Residue): The residue for which to determine the mass.
            assign_attribute (bool): If True, the mass will be assigned as an attribute to the residue.
            key (str): The key to use for the attribute.

        Returns:
            float: The mass of the residue.

        """
        mass = Mass.MASS_RESIDUE.get(residue.residue_type, 0.0)
        if assign_attribute:
            residue.set_attribute(key, mass)
        return mass

    @staticmethod
    def residue_mass_of_chain(chain: Chain,
                              assign_attribute: bool = False,
                              key: str = "residue_mass") -> float:
        """
        Returns the mass of the chain.

        Args:
            chain (Chain): The chain for which the mass will be returned.
            assign_attribute (bool): If True, the mass will be assigned as an attribute to the chain.
            key (str): The key to use for the attribute.

        Returns:
            float: The mass of the chain.
        """
        mass = sum([Mass.residue_mass_of_residue(residue, assign_attribute=assign_attribute, key=key) for residue in chain.residues])
        if assign_attribute:
            chain.set_attribute(key, mass)
        return mass

    @staticmethod
    def residue_mass_of_protein(protein: Protein,
                                assign_attribute: bool = False,
                                key: str = "residue_mass") -> float:
        """
        Returns the mass of the protein.

        Args:
            protein (Protein): The protein for which the mass will be returned.
            assign_attribute (bool): If True, the mass will be assigned as an attribute to the protein.
            key (str): The key to use for the attribute.

        Returns:
            float: The mass of the protein.
        """
        mass = sum([Mass.residue_mass_of_chain(chain, assign_attribute=assign_attribute, key=key) for chain in protein.chains])
        if assign_attribute:
            protein.set_attribute(key, mass)
        return mass

    @staticmethod
    def residue_mass_of_sequence(sequence: Sequence,
                                 assign_attribute: bool = False,
                                 key: str = "residue_mass") -> float:
        """
        Returns the mass of the sequence.

        Args:
            sequence (Sequence): The sequence for which the mass will be returned.
            assign_attribute (bool): If True, the mass will be assigned as an attribute to the sequence.
            key (str): The key to use for the attribute.

        Returns:
            float: The mass of the sequence.
        """
        mass = sum([Mass.MASS_RESIDUE.get(residue, 0.0) for residue in sequence])
        if assign_attribute:
            sequence.set_attribute(key, mass)
        return mass

    @staticmethod
    def atomic_mass_of_atom(atom: Atom,
                            assign_attribute: bool = False,
                            key: str = "atomic_mass") -> float:
        """
        Returns the mass of the atom.

        Args:
            atom (Atom): The atom for which to determine the mass.
            assign_attribute (bool): If True, the mass will be assigned as an attribute to the atom.
            key (str): The key to use for the attribute.

        Returns:
            float: The mass of the atom
        """
        mass = Mass.MASS_ATOM.get(atom.element, 0.0)
        if assign_attribute:
            atom.set_attribute(key, mass)
        return mass

    @staticmethod
    def atomic_mass_of_residue(residue: Residue,
                               assign_attribute: bool = False,
                               key: str = "atomic_mass") -> float:
        """
        Returns the mass of the residue, calculated as the sum of atomic masses.

        Args:
            residue (Residue): The residue for which to determine the mass.
            assign_attribute (bool): If True, the mass will be assigned as an attribute to the residue.
            key (str): The key to use for the attribute.

        Returns:
            float: The mass of the residue.

        """
        mass = sum([Mass.atomic_mass_of_atom(atom, assign_attribute=assign_attribute, key=key) for atom in residue.atoms])
        if assign_attribute:
            residue.set_attribute(key, mass)
        return mass

    @staticmethod
    def atomic_mass_of_chain(chain: Chain,
                             assign_attribute: bool = False,
                             key: str = "atomic_mass") -> float:
        """
        Returns the mass of the chain, calculated as the sum of atomic masses.

        Args:
            chain (Chain): The chain for which to determine the mass.
            assign_attribute (bool): If True, the mass will be assigned as an attribute to the chain.
            key (str): The key to use for the attribute.

        Returns:
            float: The mass of the chain.

        """
        mass = sum([Mass.atomic_mass_of_residue(residue, assign_attribute=assign_attribute, key=key) for residue in chain.residues])
        if assign_attribute:
            chain.set_attribute(key, mass)
        return mass

    @staticmethod
    def atomic_mass_of_protein(protein: Protein,
                               assign_attribute: bool = False,
                               key: str = "atomic_mass") -> float:
        """
        Returns the mass of the protein, calculated as the sum of atomic masses.

        Args:
            protein (Protein): The protein for which to determine the mass.
            assign_attribute (bool): If True, the mass will be assigned as an attribute to the protein.
            key (str): The key to use for the attribute.

        Returns:
            float: The mass of the protein.

        """
        mass = sum([Mass.atomic_mass_of_chain(chain, assign_attribute=assign_attribute, key=key) for chain in protein.chains])
        if assign_attribute:
            protein.set_attribute(key, mass)
        return mass

    @staticmethod
    def molecular_mass_of_residue(residue: Residue,
                                  assign_attribute: bool = False,
                                  key: str = "molecular_mass") -> float:
        """
        Returns the molecular mass of the residue.

        Args:
            residue (Residue): The residue for which to determine the mass.
            assign_attribute (bool): If True, the mass will be assigned as an attribute to the residue.
            key (str): The key to use for the attribute.

        Returns:
            float: The mass of the residue.

        """
        mass = Mass.MOLECULAR_MASS.get(residue.residue_type, 0.0)
        if assign_attribute:
            residue.set_attribute(key, mass)
        return mass

    @staticmethod
    def molecular_mass_of_chain(chain: Chain,
                                assign_attribute: bool = False,
                                key: str = "molecular_mass") -> float:
        """
        Returns the molecular mass of the chain.

        Args:
            chain (Chain): The chain for which to determine the mass.
            assign_attribute (bool): If True, the mass will be assigned as attributes to the residues in the chain.
            key (str): The key to use for the attribute.

        Returns:
            float: The mass of the chain.

        """
        mass = sum([Mass.molecular_mass_of_residue(residue, assign_attribute=assign_attribute, key=key) for residue in chain.residues])
        return mass

    @staticmethod
    def molecular_mass_of_protein(protein: Protein,
                                  assign_attribute: bool = False,
                                  key: str = "molecular_mass") -> float:
        """
        Returns the molecular mass of the protein.

        Args:
            protein (Protein): The protein for which to determine the mass.
            assign_attribute (bool): If True, the mass will be assigned as attributes to the residues in the protein.
            key (str): The key to use for the attribute.

        Returns:
            float: The mass of the protein.

        """
        mass = sum([Mass.molecular_mass_of_chain(chain, assign_attribute=assign_attribute, key=key) for chain in protein.chains])
        return mass