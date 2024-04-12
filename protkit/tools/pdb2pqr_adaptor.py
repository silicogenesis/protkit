#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS), Claudio Jardim (CJ)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `PDB2PQRAdaptor`. The PDB2PQRAdaptor class provides an interface
to the pdb2pqr software, which is used to convert PDB files to PQR files.

The pdb2pqr software is a widely used tool for adding charges and radii to protein
structures. It is used to prepare protein structures for molecular dynamics
simulations, docking studies, and other computational analyses.  pdb2pqr can add a
limited number of missing heavy atoms to structures, as well as
hydrogen atoms. It may change the coordinates of some atoms in the structure.

To use the PDB2PQRAdaptor class, you need to have the pdb2pqr software installed.
The pdb2pqr software is automatically installed as a PyPi package when you install
Protkit.

PDB2PQR was developed by the Baker Lab at the University of Washington. If you use
pdb2pqr in your research, please cite the following papers:

Dolinsky, T. J., Nielsen, J. E., McCammon, J. A., & Baker, N. A. (2004).
PDB2PQR: an automated pipeline for the setup of Poisson-Boltzmann electrostatics
calculations. Nucleic Acids Research, 32, W665-W667.
https://doi.org/10.1093/nar/gkh381

Dolinsky, T. J., Czodrowski, P., Li, H., Nielsen, J. E., Jensen, J. H., Klebe, G., & Baker, N. A. (2007).
PDB2PQR: expanding and upgrading automated preparation of biomolecular structures for molecular simulations.
Nucleic Acids Research, 35, W522-W525.
https://doi.org/10.1093/nar/gkm276

Unni, S., Huang, Y., Hanson, R. M., Tobias, M., Krishnan, S., Li, W. W., Nielsen, J. E., & Baker, N. A. (2011).
Web servers and services for electrostatics calculations with APBS and PDB2PQR.
Journal of Computational Chemistry, 32(7), 1488-1491.
https://doi.org/10.1002/jcc.21720

Jurrus, E., Engel, D., Star, K., Monson, K., Brandi, J., Felberg, L. E., Brookes, D. H., Wilson, L., Chen, J., Liles, K., Chun, M., Li, P., Gohara, D. W., Dolinsky, T., Konecny, R., Koes, D. R., Nielsen, J. E., Head-Gordon, T., Geng, W., Krasny, R., Wei, G. W., Holst, M. J., McCammon, J. A., & Baker, N. A. (2018).
Improvements to the APBS biomolecular solvation software suite.
Protein Science, 27(1), 112-128.
https://doi.org/10.1002/pro.3280
"""

import subprocess
import tempfile
import os

from protkit.file_io import PDBIO
from protkit.file_io import PQRIO
from protkit.structure import Protein


class PDB2PQRAdaptorError(Exception):
    pass


class PDB2PQRAdaptor():
    """
    Adapter class for the pdb2pqr software.
    """

    AMBER = "AMBER"
    CHARMM = "CHARMM"
    PARSE = "PARSE"
    TYL06 = "TYL06"
    PEOEPB = "PEOEPB"
    SWANSON = "SWANSON"

    def __init__(self,
                 pdb2pqr_bin_path: str = None,
                 force_field: str = PARSE):
        """
        Initialize the PDB2PQR adaptor.

        Args:
            pdb2pqr_bin_path (str): The path to the pdb2pqr binary.  If
                None, the adaptor will use the pdb2pqr binary in the PATH.
            force_field (str): The force field to use for the PQR file.
                Options: AMBER, CHARMM, PARSE, TYL06, PEOEPB, SWANSON

        Returns:
            None
        """
        super().__init__()

        self._pdb2pqr_bin_path = pdb2pqr_bin_path
        self._force_field = force_field

    def run(self,
            protein: Protein,
            output_pqr_file_path: str=None,
            output_pdb_file_path: str=None) -> Protein:
        """
        Run pdb2pqr on the specified PDB file.

        Args:
            protein (Protein): The protein to be converted to a PQR file.
            output_pqr_file_path (str): If set, the PQR file will be written to this file.
            output_pdb_file_path (str): If set, the PDB file that is used as input
                to pdb2pqr will be written to this file.

        Returns:
            Protein: The protein with the PQR file information.
        """

        # Save PDB file to serve as input to PDB2PQR
        if output_pdb_file_path is None:
            temp_file_pdb = tempfile.NamedTemporaryFile(delete=True)
            output_pdb_file_path = temp_file_pdb.name
        PDBIO.save(protein, output_pdb_file_path)

        # Allocate a temporary file for the PQR file
        if output_pqr_file_path is None:
            temp_file_pqr = tempfile.NamedTemporaryFile(delete=True)
            temp_file_pqr.close()
            output_pqr_file_path = temp_file_pqr.name

        # Ensure the directory for the log file exists
        log_directory = os.path.dirname(output_pqr_file_path)
        if not os.path.exists(log_directory):
            os.makedirs(log_directory)

        # Setup PDB2PQR
        if self._pdb2pqr_bin_path is None:
            args = ["pdb2pqr", output_pdb_file_path, output_pqr_file_path, "--ff=" + self._force_field, "--keep-chain"]
        else:
            args = [self._pdb2pqr_bin_path, output_pdb_file_path, output_pqr_file_path, "--ff=" + self._force_field, "--keep-chain"]
        pdb2pqr_process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Run PDB2PQR
        stdout, stderr = pdb2pqr_process.communicate()
        exit_code = pdb2pqr_process.returncode
        if exit_code != 0:
            raise PDB2PQRAdaptorError(
                "Error while running pdb2pqr." + stderr.decode("utf8")
            )

        # Load PQR file
        protein = PQRIO.load(output_pqr_file_path)[0]

        return protein
