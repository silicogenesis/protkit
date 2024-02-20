"""
Package `file_io` contains classes to read and write data from and to files
containing biological data.

Structural Data:

- `PDBIO` class in `pdb_io`: Reads and writes data from and to a PDB file.
- `PQRIO` class in `pqr_io`: Reads and writes data from and to a PQR file.
- `MMCIFIO` class in `mmcif_io`: Reads and writes data from and to a MMCIF file.
- `MMTFIO` class in `mmtf_io`: Reads and writes data from and to a MMTF file.
- `ProtIO` class in `prot_io`: Reads and writes data from and to a Prot file.

Sequence Data:

- `FastaIO` class in `fasta_io`: Reads and writes data from and to a FASTA file.

"""

from .pdb_io import PDBIO
from .pqr_io import PQRIO
from .mmcif_io import MMCIFIO
from .mmtf_io import MMTFIO
from .prot_io import ProtIO
from .fasta_io import FastaIO
