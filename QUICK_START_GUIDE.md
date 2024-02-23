# Protkit Quick Start Guide

Protkit is an open-source toolkit that makes it easy to work with data related to molecules and proteins and to handle various tasks in computational biology.

Protkit is designed to be intuitive and easy to use.  This Quick Start Guide will help you get started with Protkit and show you how to use some of its features.

## 1. Installation

### Installation from PyPI

The easiest way to get started with Protkit is to install it from the Python Package Index (PyPI).

Protkit requires Python 3.6 or higher.  It can be installed using `pip`:

```bash
pip install protkit
```
A number of dependencies will be installed automatically, such as `numpy`, `joblib`, `requests` and others.

See [Protkit on PyPI](https://pypi.org/project/protkit/) for more details.

### Cloning the Repository

You can also clone the repository and install it from source:

```bash
git clone https://github.com/silicogenesis/protkit.git
```

You can install the project requirements using `pip`:

```bash
pip install -r requirements.txt
```

## 2. Quick Start Example

Here is a simple example to get you started.  It illustrates how powerful computation can be done with Protkit in just a few lines of code.   

In the example, we download a PDB file from the RCSB, extract the A and B chains and do some cleanup like removing hetero atoms and fixing disordered atoms.  We then compute dihedral angles and surface areas for the protein and save it to a file.  We then load the protein from the file and print the surface area and a note that we added to the protein.

```python
from protkit.download import Download
from protkit.file_io.pdb_io import PDBIO
from protkit.file_io.prot_io import ProtIO
from protkit.properties.dihedral_angles import DihedralAngles
from protkit.properties.surface_area import SurfaceArea

# Download a PDB file from the RCSB PDB database and save it to a file.
Download.download_pdb_file_from_rcsb("1ahw", "1ahw.pdb")

# Load a PDB file into a Protein object.
protein = PDBIO.load("1ahw.pdb")[0]

# Print the number of chains in the protein.
print(protein.num_chains)

# Keep only the A and B chains
protein.keep_chains(["A", "B"])
print(protein.get_chain('A').sequence)

# Do a bit of cleanup, by removing any hetero atoms and fixing disordered atoms.
protein.remove_hetero_residues()
protein.fix_disordered_atoms()

# Compute dihedral angles for the protein, and assign them as extended attributes to residues.
DihedralAngles.dihedral_angles_of_protein(protein, assign_attribute=True)
print(protein.get_chain('A').get_residue(1).get_attribute('dihedral_angles')['PHI'])

# Compute surface areas for the protein. Surface areas are automatically computed and assigned
# at the residue, chain and protein level.
SurfaceArea.surface_area_of_protein(protein, assign_attribute=True)
print(protein.get_attribute('surface_area'))

# Save the protein to a protkit (.prot) file.  All attributes, such as the
# computed dihedral angles and surface areas, will be saved as well and
# is available for later retrieval!
protein.set_attribute("note", "Experimenting with Protkit")
ProtIO.save(protein, "1ahw.prot")
protein2 = ProtIO.load("1ahw.prot")[0]
print(protein2.get_attribute('surface_area'))
print(protein2.get_attribute('note'))
```

## 3. Downloading Files

One of the first steps in working with proteins is to download structural data
or sequence data.  Protkit provides simple methods for downloading such data.

The `Download` class in the `protkit.Download` package provides a number of methods to download molecular data
from data sources such as RCSB PDB, Uniprot and SAbDab.

Single files can be downloaded, or multiple files can be efficiently downloaded in parallel.

Most of the methods in the `Download` class take a unique identifier for a file, such
as a PDB id or Uniprot id as input. The methods also take a `file_name` argument, which is the name 
of the file on the disk it will be saved to.

### 3.1 Downloading PDB Files

PDB files can be downloaded from the RCSB PDB database or the SAbDab database,
as illustrated by the example below.

The PDB id of the file is used as an identifier of the file. The files are downloaded 
with a specified filename or directory. In the case of a directory, the directory
will automatically be created if it does not exist and the file will be saved with the
appropriate extension (.pdb).

In the case of downloading files in parallel, the `n_jobs` argument specifies the number
of parallel jobs to run. In the case that `n_jobs` is set to -1, the number of jobs will
be set to the number of cores on the machine.

```python
from protkit.download import Download

# Download a PDB file from the RCSB PDB database and save it to a file.
Download.download_pdb_file_from_rcsb("1ahw", "data/pdb_files/rcsb/1ahw.pdb")

# Download a PDB file from the SAbDab database and save it to a file.
Download.download_pdb_file_from_sabdab("1ahw", "data/pdb_files/sabdab/1ahw.pdb")

# Download multiple PDB files from the RCSB PDB in parallel and save them to a directory.
Download.download_pdb_files_from_rcsb(["1ahw", "1a4y", "1a6m"], "data/pdb_files/rcsb/", n_jobs=3)
```

### 3.2 Downloading Fasta Files

Uniprot files can be downloaded from the Uniprot database, as illustrated by the example below.

The Uniprot id of the file is used as an identifier of the file. The files are downloaded
with a specified filename or directory, as in the previous example.

```python
from protkit.download import Download

# Download a Uniprot file from the Uniprot database and save it to a file.
Download.download_uniprot_file("P12345", "data/uniprot_files/P12345.xml")

# Download multiple Uniprot files from the Uniprot database in parallel and save them to a directory.
Download.download_uniprot_files(["P12345", "P12346", "P12347"], "data/uniprot_files/", n_jobs=3)
```

### 3.3 Downloading CIF and Binary CIF Files

In a similar way to PDB files, CIF files and Binary CIF files can be downloaded from the RCSB PDB database.

```python
from protkit.download import Download

# Download a CIF file from the RCSB PDB database and save it to a file.
Download.download_cif_file_from_rcsb("1ahw", "data/cif_files/rcsb/1ahw.cif")

# Download multiple CIF files from the RCSB PDB database in parallel and save them to a directory.
Download.download_cif_files_from_rcsb(["1ahw", "1a4y", "1a6m"], "data/cif_files/rcsb/")

# Download a binary CIF file from the RCSB PDB database and save it to a file.
Download.download_binary_cif_file_from_rcsb("1ahw", "data/cif_files/rcsb/1ahw.bcif")

# Download multiple binary CIF files from the RCSB PDB database in parallel and save them to a directory.
Download.download_binary_cif_files_from_rcsb(["1ahw", "1a4y", "1a6m"], "data/cif_files/rcsb/")
```

