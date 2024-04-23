#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS), Claudio Jardim (CJ)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `PropkaAdaptor`. The PropkaAdaptor class provides an interface
to the propka software, which is used to calculate pKa values of ionizable
residues and protein-ligand complexes based on the 3D structure.

To use the PropkaAdaptor class, you need to have the propka software installed.
The propka software is automatically installed as a PyPi package when you install
Protkit.

The propka software was developed by the Jensen Lab at the University of Copenhagen.

The source code is available at:
https://github.com/jensengroup/propka

Documentation is available at:
https://propka.readthedocs.io/en/latest/

If you use propka in your research, please cite the following papers:

Olsson, M. H. M., Søndergaard, C. R., Rostkowski, M., & Jensen, J. H. (2011).
PROPKA3: Consistent Treatment of Internal and Surface Residues in Empirical pKa Predictions.
Journal of Chemical Theory and Computation, 7(2), 525-537.
https://doi.org/10.1021/ct100578z

Søndergaard, C. R., Olsson, M. H. M., Rostkowski, M., & Jensen, J. H. (2011).
Improved Treatment of Ligands and Coupling Effects in Empirical Calculation and Rationalization of pKa Values.
Journal of Chemical Theory and Computation, 7(7), 2284-2295.
https://doi.org/10.1021/ct200133y

"""

import tempfile
from typing import Dict

from protkit.file_io import PDBIO
from protkit.structure import Protein
from propka.run import single


class PropkaAdaptorError(Exception):
    pass


class PropkaAdaptor():
    """
    Adapter class for the propka software.
    """

    def __init__(self,
                 ph: float = 7.0):
        """
        Constructs all the necessary attributes for the PropkaAdaptor object.

        Args:
            ph (float): The pH value.
        """
        super().__init__()

        self.ph = ph

    def _get_pka(self, molecule) -> Dict[str, float]:
        """
        Extracts the pKa values from the propka output.

        Args:
            molecule (MolecularContainer): The propka molecule.

        Returns:
            dict: The pKa values.

        Notes:
            The implementation was adapted from
            write_pka(), get_summary_string(), get_summary_section() in Propka.
        """

        pka_values = {}
        conformation = "AVR"

        for group in molecule.conformations[conformation].groups:
            if not group.coupled_titrating_group:
                # Propka computes pKa values for additional group types as well,
                # for example N+. Further work is needed if we want to include
                # these as well.
                # if group.residue_type in molecule.version.parameters.write_out_order:
                if group.residue_type in ["ASP", "GLU", "HIS", "CYS", "TYR", "LYS", "ARG", "SER"]:
                    pka_values[group.label] = group.pka_value

        return pka_values

    def calculate_pka(self,
                      protein: Protein,
                      output_pdb_file_path: str = None) -> Protein:
        """
        Calculates the pKa values of ionizable residues in a protein structure.

        Args:
            protein (Protein): The protein.
            output_pdb_file_path (str): Optional. The path to the output PDB file.

        Returns:
            Protein: The protein with pKa values.
        """
        if output_pdb_file_path is None:
            temp_file_pdb = tempfile.NamedTemporaryFile(delete=True, suffix=".pdb")
            output_pdb_file_path = temp_file_pdb.name
        PDBIO.save(protein, output_pdb_file_path)

        # TODO: Propka offers the ability to read data from a string stream.
        # In case that the PDB file is not written to disk, we can use this feature.
        molecule = single(output_pdb_file_path, optargs=[f"--pH={self.ph}"], write_pka=False)
        pkas = self._get_pka(molecule)

        for label in pkas:
            residue_type, residue_no, chain_id = label.split()
            pka = pkas[label]
            chain = protein.get_chain(chain_id)
            # TODO: Need to verify how Propka handles residues with insertion codes.
            residue = chain.get_residue_by_sequence_code(residue_no)
            assert residue.residue_type == residue_type
            residue.set_attribute("pka", pka)

        protein_out = protein.copy()

        return protein_out
