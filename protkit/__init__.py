"""
# Protkit

Protkit is an open source Python library that can be used for a variety of tasks in computational biology
and bioinformatics, focusing on structural bioinformatics, protein engineering and machine learning.

It is designed to support the broad community of computational biologists,
bioinformaticians and machine learning researchers in academia, industry
and government labs.

Protkit can be used for a variety of computational biology tasks across the computational biology pipeline, such as:

- **Reading and writing data** from popular structure file formats, such as
    PDB, PQR, MMTF, mmCIF; and sequence file formats, such as FASTA.
- **Downloading** data from popular databases of protein structures, such as the PDB RCSB, UniProt and SAbDab.
- **Data structures** for representing proteins, protein complexes, chains,
    residues, atoms and sequences. These data structures provide capabilities to extract data
    in both hierarchical and linear formats. It is extensible and easy to add
    new properties to the data structure. It has a rich set of methods for extracting
    and filtering data from the data structure.
- **Detecting and fixing anomalies** in protein structures, such as missing atoms,
    missing residues, detecting sequence gaps, detecting atomic clashes, removing
    hetero residues or water molecules, and removing alternate conformations.
- **Calculating properties** of proteins, such as hydrophobicity, charge, surface areas,
    secondary structures, dihedral angles, interface residues and more.
- **Geometric operations** on proteins, such as aligning and superimposing
    structures.
- **Metrics** for comparing proteins, such as RMSD and Sequence Similarity.
- **Featurization** of proteins and their properties enabling preparation of datasets
    for **machine learning** applications.
- Performing and enabling a large variety of **computational tasks** on proteins,
    such as protein folding, protein docking, protein-protein binding affinity prediction,
    humanisation of antibodies, prediction of developability characteristics etc. Care is taken
    that the various tools are interoperable and can be used together in a seamless manner.

Protkit is an open source library that is free to use and modify.  We welcome
contributions from the community.

---

## Installation

### Installation from PyPI

`protkit` requires Python 3.6 or higher.  It can be installed using `pip`:

```bash
pip install protkit
```

A number of dependencies will be installed automatically, such as `numpy`, `joblib`, `requests` and others.

See [Protkit](https://pypi.org/project/protkit/) on PyPI for more details.

### Cloning the Repository

You can also clone the repository and install it from source:

```bash
git clone https://github.com/silicogenesis/protkit.git
```

You can install the project requirements using `pip`:

```bash
pip install -r requirements.txt
```

---

## Quick Start Example

Protkit is designed to be intuitive and easy to use.  An extensive set of examples can be found in the [Quick Start Guide](QUICK_START_GUIDE.md).

Here is a simple example to get you started.  It illustrates how powerful computation can be done with Protkit in just a few lines of code.

In the example, we download a PDB file from the RCSB, extract the A and B chains and do some cleanup like removing hetero atoms and fixing disordered atoms.  We then compute dihedral angles and surface areas for the protein and save it to a file.  We then load the protein from the file and print the surface area and a note that we added to the protein.

```python
from protkit.download import Download
from protkit.file_io import PDBIO, ProtIO
from protkit.properties import DihedralAngles, SurfaceArea

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

Please consult the [Quick Start Guide](QUICK_START_GUIDE.md) for more examples.

---
"""