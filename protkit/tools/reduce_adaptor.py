#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS), Claudio Jardim (CJ)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class ReduceAdaptor. The ReduceAdaptor class provides an interface
to the Reduce software, which is used to add hydrogens to protein structures.

The Reduce software is a widely used tool for adding hydrogen atoms to protein
structures. It is used to prepare protein structures for molecular dynamics
simulations, docking studies, and other computational analyses.

To use the ReduceAdaptor class, you need to have the Reduce software installed.
You can download Reduce from the following website:
https://github.com/rlabduke/reduce

Reduce was developed by the Richardson Lab at Duke University. If you use
Reduce in your research, please cite the following paper:

Word, J. M., Lovell, S. C., Richardson, J. S., & Richardson, D. C. (1999).
Asparagine and glutamine: using hydrogen atom contacts in the choice of side-chain
amide orientation. Journal of Molecular Biology, 285(4), 1735-1747.
https://doi.org/10.1006/jmbi.1998.2401
"""

import subprocess

from protkit.file_io import PDBIO
from protkit.structure import Protein
from protkit.tasks.protonator import Protonator


class ReduceProtonationError(Exception):
    pass


class ReduceDeprotonationError(Exception):
    pass


class ReduceAdaptor(Protonator):
    """
    Adapter class for the Reduce software.

    This class provides a concrete implementation of the Protonator class,
    by interfacing with the Reduce software and updating Protein objects
    with the protonated structure.
    """

    def __init__(self,
                 reduce_bin_path: str,
                 quiet: bool = False):
        """
        Initialize the Reduce adaptor.

        Args:
            reduce_bin_path (str): The path to the Reduce binary.
            quiet (bool): Whether to suppress Reduce output.
        """
        super().__init__()

        self._reduce_bin_path = reduce_bin_path
        self._quiet = quiet

    @staticmethod
    def merge_copy(source_protein: Protein, updated_protein: Protein) -> Protein:
        """
        Merge the updated protein into the source protein. The source protein
        is copied before the merge operation.

        Args:
            source_protein (Protein): The source protein.
            updated_protein (Protein): The updated protein.

        Returns:
            Protein: The merged protein.
        """
        merged_protein = source_protein.copy()
        ReduceAdaptor.merge_into(merged_protein, updated_protein)
        return merged_protein

    @staticmethod
    def merge_into(source_protein: Protein, updated_protein: Protein):
        """
        Merge the updated protein into the source protein. The source
        protein is updated in place.

        Args:
            source_protein (Protein): The source protein.
            updated_protein (Protein): The updated protein.

        Returns:
            None
        """
        for chain in updated_protein.chains:
            source_chain = source_protein.get_chain(chain.chain_id)
            if source_chain is None:
                source_protein.add_chain(chain.chain_id, chain)
            else:
                for index, residue in enumerate(chain.residues):
                    source_residue = source_chain.get_residue(index)
                    if source_residue is None:
                        source_chain.add_residue(residue)
                    else:
                        for atom in residue.atoms:
                            source_atom = source_residue.get_atom(atom.atom_type)
                            if source_atom is None:
                                source_residue.add_atom(atom.atom_type, atom)
                            else:
                                source_atom.x = atom.x
                                source_atom.y = atom.y
                                source_atom.z = atom.z
                                if atom.has_attribute("temp_factor"):
                                    source_atom.set_attribute("temp_factor", atom.get_attribute("temp_factor"))

    def _reduce_protonate(self, pdb_str: str) -> str:
        """
        Protonate a protein using the Reduce software.

        Args:
            pdb_str (str): The PDB file as a string.

        Returns:
            str: The protonated PDB file as a string.
        """
        if self._quiet:
            args = [self._reduce_bin_path, "-quiet", "-build", "-"]
        else:
            args = [self._reduce_bin_path, "-build", "-"]
        protonate_process = subprocess.Popen(args,
                                             stdin=subprocess.PIPE,
                                             stdout=subprocess.PIPE)
        protonate_process.stdin.write(pdb_str.encode())
        stdout, stderr = protonate_process.communicate()
        protonate_process_out = stdout.decode("utf8")
        protonate_process.stdin.close()
        if protonate_process.stderr is not None:
            raise ReduceProtonationError(
                "Error while protonating the protein."
            )
        return protonate_process_out

    def _reduce_deprotonate(self, pdb_str: str) -> str:
        """
        Deprotonate a protein using the Reduce software.

        Args:
            pdb_str (str): The PDB file as a string.

        Returns:
            str: The deprotonated PDB file as a string.
        """
        if self._quiet:
            args = [self._reduce_bin_path, "-quiet", "-Trim", "-"]
        else:
            args = [self._reduce_bin_path, "-Trim", "-"]
        deprotonate_process = subprocess.Popen(args,
                                               stdin=subprocess.PIPE,
                                               stdout=subprocess.PIPE)
        deprotonate_process.stdin.write(pdb_str.encode())
        stdout, stderr = deprotonate_process.communicate()
        deprotonate_process_output = stdout.decode("utf8")
        deprotonate_process.stdin.close()
        if deprotonate_process.stderr is not None:
            raise ReduceDeprotonationError(
                "Error while deprotonating the protein before protonation."
            )
        return deprotonate_process_output

    def protonate(self, protein: Protein, reduce_output_pdb_file_path: str=None) -> Protein:
        """
        Protonate a protein using the Reduce software.

        Args:
            protein (Protein): The protein to be protonated.
            reduce_output_pdb_file_path (str): If set, the protonated PDB file
                as produced by Reduce will be saved to this path.

        Returns:
            Protein: The protonated protein.
        """

        # Create an in memory PDB file
        input_pdb_str = PDBIO.save_to_string(protein)

        # Remove protons first, in case the structure is already protonated
        deprotonated_pdb_str = self._reduce_deprotonate(input_pdb_str)

        # Add hydrogens to the protein.
        protonated_pdb_str = self._reduce_protonate(deprotonated_pdb_str)
        if reduce_output_pdb_file_path is not None:
            with open(reduce_output_pdb_file_path, "w") as file:
                file.write(protonated_pdb_str)

        # Merge the protonated protein with the original protein
        protein_out = PDBIO.load_from_string(protonated_pdb_str)[0]
        protein_merged = ReduceAdaptor.merge_copy(protein, protein_out)

        return protein_merged

    def deprotonate(self, protein: Protein, reduce_output_pdb_file_path: str=None) -> Protein:
        """
        Deprotonate a protein using the Reduce software.

        Args:
            protein (Protein): The protein to be deprotonated.
            reduce_output_pdb_file_path (str): If set, the deprotonated PDB file
                as produced by Reduce will be saved to this path.

        Returns:
            Protein: The deprotonated protein.
        """

        # Create an in memory PDB file
        input_pdb_str = PDBIO.save_to_string(protein)

        # Remove hydrogen atoms from the structure.
        deprotonated_pdb_str = self._reduce_deprotonate(input_pdb_str)
        if reduce_output_pdb_file_path is not None:
            with open(reduce_output_pdb_file_path, "w") as file:
                file.write(deprotonated_pdb_str)

        # Merge the deprotonated protein with the original protein
        protein_out = PDBIO.load_from_string(deprotonated_pdb_str)[0]
        protein_merged = ReduceAdaptor.merge_copy(protein, protein_out)

        return protein_merged
