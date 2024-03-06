#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `DonorsAcceptors` to represent the donors and acceptors of hydrogen bonds.
"""

from typing import List

from protkit.structure.atom import Atom
from protkit.structure.residue import Residue
from protkit.structure.chain import Chain
from protkit.structure.protein import Protein
from protkit.seq.sequence import Sequence


class DonorsAcceptors:
    # IMGT classifies the amino acids into donors and acceptors of hydrogen bonds.
    # The following dictionary is based on the IMGT classification.
    # Donors: R, K, W
    # Acceptors: D, E
    # Both: N, Q, H, S, T, Y
    # Neither: A, C, F, G, I, L, M, P, V
    DONOR_RESIDUES = {
        "ALA": False,
        "ARG": True,
        "ASN": True,
        "ASP": False,
        "CYS": False,
        "GLN": True,
        "GLU": False,
        "GLY": False,
        "HIS": True,
        "ILE": False,
        "LEU": False,
        "LYS": True,
        "MET": False,
        "PHE": False,
        "PRO": False,
        "SER": True,
        "THR": True,
        "TRP": True,
        "TYR": True,
        "VAL": False,

        "A": False,
        "R": True,
        "N": True,
        "D": False,
        "C": False,
        "Q": True,
        "E": False,
        "G": False,
        "H": True,
        "I": False,
        "L": False,
        "K": True,
        "M": False,
        "F": False,
        "P": False,
        "S": True,
        "T": True,
        "W": True,
        "Y": True,
        "V": False
    }

    ACCEPTOR_RESIDUES = {
        "ALA": False,
        "ARG": False,
        "ASN": True,
        "ASP": True,
        "CYS": False,
        "GLN": True,
        "GLU": True,
        "GLY": False,
        "HIS": True,
        "ILE": False,
        "LEU": False,
        "LYS": False,
        "MET": False,
        "PHE": False,
        "PRO": False,
        "SER": True,
        "THR": True,
        "TRP": False,
        "TYR": True,
        "VAL": False,

        "A": False,
        "R": False,
        "N": True,
        "D": True,
        "C": False,
        "Q": True,
        "E": True,
        "G": False,
        "H": True,
        "I": False,
        "L": False,
        "K": False,
        "M": False,
        "F": False,
        "P": False,
        "S": True,
        "T": True,
        "W": False,
        "Y": True,
        "V": False
    }

    # Donor and acceptor atoms as defined in
    # https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/charge/index.html#hydrogen
    # with the addition of Nitrogen for donor atoms
    # and Oxygen and Terminal Oxygen for acceptor atoms.

    DONOR_ATOMS = {
        "ALA": {"N"},
        "ARG": {"N", "NE", "NH1", "NH2"},
        "ASN": {"N", "ND2"},  # OD1?
        "ASP": {"N"},
        "CYS": {"N"},  # SG?
        "GLU": {"N"},
        "GLN": {"N", "NE2"},  # OE1?
        "GLY": {"N"},
        "HIS": {"N", "ND1", "NE2"},  # CD2?, CE1?
        "ILE": {"N"},
        "LEU": {"N"},
        "LYS": {"N", "NZ"},
        "MET": {"N"},
        "PHE": {"N"},
        "PRO": {},
        "SER": {"N", "OG"},
        "THR": {"N", "OG1"},
        "TRP": {"N", "NE1"},
        "TYR": {"N", "OH"},
        "VAL": {"N"},

        "A": {"N"},
        "R": {"N", "NE", "NH1", "NH2"},
        "N": {"N", "ND2"},  # OD1?
        "D": {"N"},
        "C": {"N"},  # SG?
        "Q": {"N"},
        "E": {"N", "NE2"},  # OE1?
        "G": {"N"},
        "H": {"N", "ND1", "NE2"},  # CD2?, CE1?
        "I": {"N"},
        "L": {"N"},
        "K": {"N", "NZ"},
        "M": {"N"},
        "F": {"N"},
        "P": {},
        "S": {"N", "OG"},
        "T": {"N", "OG1"},
        "W": {"N", "NE1"},
        "Y": {"N", "OH"},
        "V": {"N"},
    }

    ACCEPTOR_ATOMS = {
        "ALA": {"O", "OXT"},
        "ARG": {"O", "OXT"},
        "ASN": {"O", "OD1", "OXT"},  # ND2?
        "ASP": {"O", "OD1", "OD2", "OXT"},
        "CYS": {"O", "OXT"},  # SG?
        "GLU": {"O", "OE1", "OE2", "OXT"},
        "GLN": {"O", "OE1", "OXT"},  # NE2?
        "GLY": {"O", "OXT"},
        "HIS": {"O", "ND1", "NE2", "OXT"},  # CD2?, CE1?
        "ILE": {"O", "OXT"},
        "LEU": {"O", "OXT"},
        "LYS": {"O", "OXT"},
        "MET": {"O", "OXT"},  # SD?
        "PHE": {"O", "OXT"},
        "PRO": {"O", "OXT"},
        "SER": {"O", "OG", "OXT"},
        "THR": {"O", "OG1", "OXT"},
        "TRP": {"O", "OXT"},
        "TYR": {"O", "OH", "OXT"},
        "VAL": {"O", "OXT"},

        "A": {"O", "OXT"},
        "R": {"O", "OXT"},
        "N": {"O", "OD1", "OXT"},  # ND2?
        "D": {"O", "OD1", "OD2", "OXT"},
        "C": {"O", "OXT"},  # SG?
        "Q": {"O", "OE1", "OE2", "OXT"},
        "E": {"O", "OE1", "OXT"},  # NE2?
        "G": {"O", "OXT"},
        "H": {"O", "ND1", "NE2", "OXT"},  # CD2?, CE1?
        "I": {"O", "OXT"},
        "L": {"O", "OXT"},
        "K": {"O", "OXT"},
        "M": {"O", "OXT"},  # SD?
        "F": {"O", "OXT"},
        "P": {"O", "OXT"},
        "S": {"O", "OG", "OXT"},
        "T": {"O", "OG1", "OXT"},
        "W": {"O", "OXT"},
        "Y": {"O", "OH", "OXT"},
        "V": {"O", "OXT"}
    }

    @staticmethod
    def donor_residue(residue: Residue,
                      assign_attribute: bool = False,
                      key: str = "is_donor_residue") -> bool:
        """
        Returns whether the residue is a donor of hydrogen bonds.

        Args:
            residue (str): The residue for which to check if it is a donor.
            assign_attribute (bool): Whether to assign the donor property to the residue.
            key (str): The key to use for the attribute.

        Returns:
            bool: True if the residue is a donor, False otherwise.
        """
        is_donor = DonorsAcceptors.DONOR_RESIDUES.get(residue.residue_type, False)
        if assign_attribute:
            residue.set_attribute(key, is_donor)
        return is_donor

    @staticmethod
    def acceptor_residue(residue: Residue,
                         assign_attribute: bool = False,
                         key: str = "is_acceptor_residue") -> bool:
        """
        Returns whether the residue is an acceptor of hydrogen bonds.

        Args:
            residue (str): The residue for which to check if it is an acceptor.
            assign_attribute (bool): Whether to assign the acceptor property to the residue.
            key (str): The key to use for the attribute.

        Returns:
            bool: True if the residue is an acceptor, False otherwise.
        """
        is_acceptor = not DonorsAcceptors.ACCEPTOR_RESIDUES.get(residue.residue_type, False)
        if assign_attribute:
            residue.set_attribute(key, is_acceptor)
        return is_acceptor

    @staticmethod
    def donor_residues_of_chain(chain: Chain,
                                assign_attribute: bool = False,
                                key: str = "is_donor_residue") -> list:
        """
        Returns the donor residues of a chain.

        Args:
            chain (Chain): The chain for which to determine the donor residues.
            assign_attribute (bool): Whether to assign the donor residues to the chain.
            key (str): The key to use for the attribute.

        Returns:
            list: A list of donor residues.
        """
        donor_residues = [DonorsAcceptors.donor_residue(residue, assign_attribute=assign_attribute, key=key) for residue in chain.residues]
        if assign_attribute:
            chain.set_attribute(key, donor_residues)
        return donor_residues

    @staticmethod
    def acceptor_residues_of_chain(chain: Chain,
                                   assign_attribute: bool = False,
                                   key: str = "is_acceptor_residue") -> List[bool]:
        """
        Returns the acceptor residues of a chain.

        Args:
            chain (Chain): The chain for which to determine the acceptor residues.
            assign_attribute (bool): Whether to assign the acceptor residues to the chain.
            key (str): The key to use for the attribute.

        Returns:
            list: A list of acceptor residues.
        """
        acceptor_residues = [DonorsAcceptors.acceptor_residue(residue, assign_attribute=assign_attribute, key=key) for residue in chain.residues]
        if assign_attribute:
            chain.set_attribute(key, acceptor_residues)
        return acceptor_residues

    @staticmethod
    def donor_residues_of_protein(protein: Protein,
                                  assign_attribute: bool = False,
                                  key: str = "is_donor_residue") -> List[List[bool]]:
        """
        Returns the donor residues of a protein.

        Args:
            protein (Protein): The protein for which to determine the donor residues.
            assign_attribute (bool): Whether to assign the donor residues to the protein.
            key (str): The key to use for the attribute.

        Returns:
            list: A list of donor residues.
        """
        donor_residues = [DonorsAcceptors.donor_residues_of_chain(chain, assign_attribute=assign_attribute, key=key) for chain in protein.chains]
        if assign_attribute:
            protein.set_attribute(key, donor_residues)
        return donor_residues

    @staticmethod
    def acceptor_residues_of_protein(protein: Protein,
                                     assign_attribute: bool = False,
                                     key: str = "is_acceptor_residue") -> List[List[bool]]:
        """
        Returns the acceptor residues of a protein.

        Args:
            protein (Protein): The protein for which to determine the acceptor residues.
            assign_attribute (bool): Whether to assign the acceptor residues to the protein.
            key (str): The key to use for the attribute.

        Returns:
            list: A list of acceptor residues.
        """
        acceptor_residues = [DonorsAcceptors.acceptor_residues_of_chain(chain, assign_attribute=assign_attribute, key=key) for chain in protein.chains]
        if assign_attribute:
            protein.set_attribute(key, acceptor_residues)
        return acceptor_residues

    @staticmethod
    def donor_residues_of_sequence(sequence: Sequence,
                                   assign_attribute: bool = False,
                                   key: str = "is_donor_residue") -> List[bool]:
        """
        Returns the donor residues of a sequence.

        Args:
            sequence (Sequence): The sequence for which to determine the donor residues.
            assign_attribute (bool): Whether to assign the donor residues to the sequence.
            key (str): The key to use for the attribute.

        Returns:
            list: A list of donor residues.
        """
        donor_residues = [DonorsAcceptors.DONOR_RESIDUES.get(residue, False) for residue in sequence]
        if assign_attribute:
            sequence.set_attribute(key, donor_residues)
        return donor_residues

    @staticmethod
    def acceptor_residues_of_sequence(sequence: Sequence,
                                      assign_attribute: bool = False,
                                      key: str = "is_acceptor_residue") -> List[bool]:
        """
        Returns the acceptor residues of a sequence.

        Args:
            sequence (Sequence): The sequence for which to determine the acceptor residues.
            assign_attribute (bool): Whether to assign the acceptor residues to the sequence.
            key (str): The key to use for the attribute.

        Returns:
            list: A list of acceptor residues.
        """
        acceptor_residues = [not DonorsAcceptors.ACCEPTOR_RESIDUES.get(residue, False) for residue in sequence]
        if assign_attribute:
            sequence.set_attribute(key, acceptor_residues)
        return acceptor_residues

    # ----------------- Donor and Acceptor Atoms -----------------

    @staticmethod
    def donor_atom(atom: Atom,
                   assign_attribute: bool = False,
                   key: str = "is_donor_atom"):
        """
        Returns whether the atom is a donor atom.

        Args:
            atom (Atom): The atom for which to check if it is a donor atom.
            assign_attribute (bool): Whether to assign the donor property to the atom.
            key (str): The key to use for the attribute.

        Returns:
            bool: True if the atom is a donor, False otherwise.
        """
        is_donor = atom.atom_type in DonorsAcceptors.DONOR_ATOMS.get(atom.residue.residue_type, set())
        if assign_attribute:
            atom.set_attribute(key, is_donor)
        return is_donor

    @staticmethod
    def acceptor_atom(atom: Atom,
                      assign_attribute: bool = False,
                      key: str = "is_acceptor_atom"):
        """
        Returns whether the atom is an acceptor atom.

        Args:
            atom (Atom): The atom for which to check if it is an acceptor atom.
            assign_attribute (bool): Whether to assign the acceptor property to the atom.
            key (str): The key to use for the attribute.

        Returns:
            bool: True if the atom is an acceptor, False otherwise.
        """
        is_acceptor = atom.atom_type in DonorsAcceptors.ACCEPTOR_ATOMS.get(atom.residue.residue_type, set())
        if assign_attribute:
            atom.set_attribute(key, is_acceptor)
        return is_acceptor

    @staticmethod
    def donor_atoms_of_residue(residue: Residue,
                               assign_attribute: bool = False,
                               key: str = "is_donor_atom") -> List[bool]:
        """
        Returns the donor atoms of a residue.

        Args:
            residue (Residue): The residue for which to determine the donor atoms.
            assign_attribute (bool): Whether to assign the donor atoms to the residue.
            key (str): The key to use for the attribute.

        Returns:
            set: A set of donor atoms.
        """
        donor_atoms = [DonorsAcceptors.donor_atom(atom, assign_attribute=assign_attribute, key=key) for atom in residue.atoms]
        if assign_attribute:
            residue.set_attribute(key, donor_atoms)
        return donor_atoms

    @staticmethod
    def acceptor_atoms_of_residue(residue: Residue,
                                  assign_attribute: bool = False,
                                  key: str = "is_acceptor_atom") -> List[bool]:
        """
        Returns the acceptor atoms of a residue.

        Args:
            residue (Residue): The residue for which to determine the acceptor atoms.
            assign_attribute (bool): Whether to assign the acceptor atoms to the residue.
            key (str): The key to use for the attribute.

        Returns:
            set: A set of acceptor atoms.
        """
        acceptor_atoms = [DonorsAcceptors.acceptor_atom(atom, assign_attribute=assign_attribute, key=key) for atom in residue.atoms]
        if assign_attribute:
            residue.set_attribute(key, acceptor_atoms)
        return acceptor_atoms

    @staticmethod
    def donor_atoms_of_chain(chain: Chain,
                             assign_attribute: bool = False,
                             key: str = "is_donor_atom") -> List[List[bool]]:
        """
        Returns the donor atoms of a chain.

        Args:
            chain (Chain): The chain for which to determine the donor atoms.
            assign_attribute (bool): Whether to assign the donor atoms to the chain.
            key (str): The key to use for the attribute.

        Returns:
            list: A list of donor atoms.
        """
        donor_atoms = [DonorsAcceptors.donor_atoms_of_residue(residue, assign_attribute=assign_attribute, key=key) for residue in chain.residues]
        if assign_attribute:
            chain.set_attribute(key, donor_atoms)
        return donor_atoms

    @staticmethod
    def acceptor_atoms_of_chain(chain: Chain,
                                assign_attribute: bool = False,
                                key: str = "is_acceptor_atom") -> List[List[bool]]:
        """
        Returns the acceptor atoms of a chain.

        Args:
            chain (Chain): The chain for which to determine the acceptor atoms.
            assign_attribute (bool): Whether to assign the acceptor atoms to the chain.
            key (str): The key to use for the attribute.

        Returns:
            list: A list of acceptor atoms.
        """
        acceptor_atoms = [DonorsAcceptors.acceptor_atoms_of_residue(residue, assign_attribute=assign_attribute, key=key) for residue in chain.residues]
        if assign_attribute:
            chain.set_attribute(key, acceptor_atoms)
        return acceptor_atoms

    @staticmethod
    def donor_atoms_of_protein(protein: Protein,
                               assign_attribute: bool = False,
                               key: str = "is_donor_atom") -> List[List[List[bool]]]:
        """
        Returns the donor atoms of a protein.

        Args:
            protein (Protein): The protein for which to determine the donor atoms.
            assign_attribute (bool): Whether to assign the donor atoms to the protein.
            key (str): The key to use for the attribute.

        Returns:
            list: A list of donor atoms.
        """
        donor_atoms = [DonorsAcceptors.donor_atoms_of_chain(chain, assign_attribute=assign_attribute, key=key) for chain in protein.chains]
        if assign_attribute:
            protein.set_attribute(key, donor_atoms)
        return donor_atoms

    @staticmethod
    def acceptor_atoms_of_protein(protein: Protein,
                                  assign_attribute: bool = False,
                                  key: str = "is_acceptor_atom") -> List[List[List[bool]]]:
        """
        Returns the acceptor atoms of a protein.

        Args:
            protein (Protein): The protein for which to determine the acceptor atoms.
            assign_attribute (bool): Whether to assign the acceptor atoms to the protein.
            key (str): The key to use for the attribute.

        Returns:
            list: A list of acceptor atoms.
        """
        acceptor_atoms = [DonorsAcceptors.acceptor_atoms_of_chain(chain, assign_attribute=assign_attribute, key=key) for chain in protein.chains]
        if assign_attribute:
            protein.set_attribute(key, acceptor_atoms)
        return acceptor_atoms
