# Protkit Quick Start Guide

Protkit is an open-source toolkit that makes it easy to work with data related to molecules and proteins and to handle various tasks in computational biology.

Protkit is designed to be intuitive and easy to use.  This Quick Start Guide will help you get started with Protkit and show you how to use some of its features.

---

## Table of Contents

1. [Installation](#1-installation)
2. [Quick Start Example](#2-quick-start-example)
3. [Downloading Files](#3-downloading-files)
    1. [Downloading PDB Files](#31-downloading-pdb-files)
    2. [Downloading Fasta Files](#32-downloading-fasta-files)
    3. [Downloading CIF and Binary CIF Files](#33-downloading-cif-and-binary-cif-files)
4. [File I/O](#4-file-io)
    1. [Working with PDB Files](#41-working-with-pdb-files)
    2. [Working with Prot Files](#42-working-with-prot-files)
    3. [Working with Fasta Files](#43-working-with-fasta-files)
    4. [Working with PQR Files](#44-working-with-pqr-files)
    5. [Working with mmCIF Files](#45-working-with-mmcif-files)
    6. [Working with MMTF Files](#46-working-with-mmtf-files)
5. [Representing Protein Structure](#5-representing-protein-structure)
    1. [Aims with the Molecular Data Structures](#51-aims-with-the-molecular-data-structures)
    2. [Protein, Chain, Residue and Atom Class Hierarchy](#52-protein-chain-residue-and-atom-class-hierarchy)
    3. [Protein Class](#53-protein-class)
    4. [Chain Class](#54-chain-class)
    5. [Residue Class](#55-residue-class)
    6. [Atom Class](#56-atom-class)
    7. [Extending the Structure Classes with Properties](#57-extending-the-structure-classes-with-properties)
6. [Ensuring Data Quality](#6-ensuring-data-quality)
    1. [Creating a Copy of a Protein](#61-creating-a-copy-of-a-protein)
    2. [Removing Water Molecules](#62-removing-water-molecules)
    3. [Removing Hetero Residues](#63-removing-hetero-residues)
    4. [Fixing Disordered Atoms](#64-fixing-disordered-atoms)
    5. [Removing Hydrogen Atoms (Deprotonation)](#65-removing-hydrogen-atoms-deprotonation) 

---

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

---

## 2. Quick Start Example

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
---

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

---

## 4. File I/O

Protkit provides a number of classes for reading and writing molecular data files.  These classes are located in the `protkit.file_io` package.  For example, the `PDBIO` class provides methods for reading and writing PDB files, while the `ProtIO` class provides methods for reading and writing Prot files, etc.

Generally, the classes provide a method to read one or more proteins from a file, using the `load()` method, and a method to write one or more proteins to a file, using the `save()` method. Conversions between different file formats are also possible, through the use of the `convert()` method.

The `load()` method returns a list of `Protein` objects, while the `save()` method takes a list of `Protein` objects as input. We will discuss the `Protein` class in more detail in the next section.

### 4.1 Working with PDB Files

#### Reading PDB Files

The PDB file format is a standard format for storing structural information about proteins. OpenSG provides a simple interface for reading and writing PDB files.

```python
from protkit.file_io import PDBIO

protein = PDBIO.load("1ahw.pdb")[0]
```

The `load()` method returns a list of proteins. Generally, PDB files generated through X-ray crystallography contain only one protein, while NMR files may contain multiple proteins. In Protkit, the `load()` method always returns a list of proteins, even if the file contains only one protein.   

The above code loads the protein 1AHW from the file 1AHW.pdb. In this case, the list contains only one protein, so we take the first element of the list, using the `[0]` index in the array.

The `Protein` object contains all the information about the protein. The `Protein` class is discussed in more detail in the next section.

#### Saving PDB Files

Saving a protein to a PDB file is accomplished using the `save()` method. Let's modify the 
protein by extracting the A and B chains and saving it to a new file.

```python
protein.keep_chains(["A", "B"])

PDBIO.save(protein, "1ahw_AB.pdb")
```

#### PDB Metadata

PDB files contain metadata such as the resolution of the structure, the method used to solve the structure, the date the structure was deposited, etc.

Protkit does not currently read such metadata. Since Protkit is primarily aimed towards structural biology, it reads 
data such as the sequence, the coordinates of the atoms, the occupancy and temperature factors, etc. contained
in records such as ATOM, HETATM, MODEL, SEQRES, etc.

Future versions of Protkit may include methods to read and write such metadata.

### 4.2 Working with Prot Files

Although the PDB file format is the dominant format for storing structural information about proteins, it is an old format that has many limitations.  For example, it cannot store additional information about the protein, such as surface areas, molecular masses or hydrophobicity values.

To overcome the limitations, we have developed the Prot file format, which is a simple (JSON) text-based format that can store all the information about a protein in a single file. The Prot file format is designed to be easy to read and write, and to be easily extensible.  The Prot file format closely resembles the `Protein` data structure (discussed later), and can be efficiently serialised and deserialised.

When working with protein data, it is recommended to use the Prot file format.

#### Advantages of Prot Files over PDB files

- Prot files can store additional information about the protein, such as surface areas, molecular masses, hydrophobicity values, etc.
- Prot files are based on JSON, which is a text-based format that is easy to read and write. Many serialization libraries are available for JSON, and it is supported by many programming languages.
- Prot files can be compressed using the gzip algorithm, which reduces the file size. Prot files in compressed JSON format is typically smaller than similar PDB files. Compression is enabled by default.
- Prot files can store multiple (different) proteins in a single file, which facilitates the creation of databases of proteins that can easily be managed through a single file. 

#### Converting PDB Files to Prot Files

The easiest way to create a Prot file is to convert a PDB file to a Prot file.
The conversion returns a Protein object. The following code converts the
protein 1AHW from a PDB file to a Prot file.

```python
from protkit.file_io import ProtIO

protein = ProtIO.convert("1AHW.pdb", "1AHW.prot") 
```

The Prot file could also have been generated by first loading the protein from
the PDB file and then saving it to a Prot file.

#### Reading Prot Files

The ProtIO class provides a `load()` method for reading Prot files, in the same 
way as the PDBIO class provides a load method for reading PDB files.

```python
protein = ProtIO.load("1AHW.prot")[0]
```

#### Saving Prot files

The ProtIO class provides a `save()` method for saving Prot files, as in the example:

```python
ProtIO.save(protein, "1AHW.prot")
```

#### Saving Multiple Proteins in a Single Prot File

Note that multiple proteins can be saved to a single Prot file. This facilitates
the creation of databases of proteins that can easily be managed through a
single file.

```python
protein1 = PDBIO.load("1AHW.pdb")[0]
protein2 = PDBIO.load("4NKQ.pdb")[0]

ProtIO.save([protein1, protein2], "database.prot")
``` 

#### Compression and Decompression of Prot Files

Prot files are based on the JSON format, which is a text-based format. This
could lead to large file sizes. To reduce the file size, Prot files are
compressed using the gzip algorithm. Compression can be enabled or disabled
by parameters passed to the save and load functions.  Compression is enabled
by default.

```python
protein = ProtIO.load("1AHW.prot", decompress=True)[0]
ProtIO.save(protein, "1AHW.prot.json", compress=False)
```

### 4.3 Working with Fasta Files

The Fasta file format is a standard format for storing sequence information about proteins. Protkit provides a simple interface for reading and writing Fasta files.

#### Reading Fasta Files

The FastaIO class provides a `load()` method for reading Fasta files. The `load()` method returns 
a list of `Sequence` objects, describing the sequences in the file. The sequence objects contain the
sequence and the header of the sequence. The `Sequence` class is discussed in more detail later in this guide.

```python

```python
from protkit.file_io import FastaIO

sequences = FastaIO.load("1AHW.fasta")
```

#### Saving Fasta Files

The FastaIO class provides a `save()` method for saving Fasta files. The `save()` method takes a single `Sequence` object or a list of `Sequence` objects as input.

```python
from protkit.file_io import FastaIO

sequences = FastaIO.load("1AHW.fasta")
FastaIO.save(sequences, "1AHW_copy.fasta")
FastaIO.save(sequences[0], "1AHW_A.fasta")
```

### 4.4 Working with PQR Files

Protkit supports the PQR file format, which is a variant of the PDB file format that includes atomic charges and radii. The PQR file format is used in molecular dynamics simulations and in the calculation of electrostatic potentials.

It is typically created by programs such as APBS, PDB2PQR, etc.

#### Reading PQR Files

The PQRIO class provides a `load()` method for reading PQR files. The `load()` method returns a list of `Protein` objects, describing the proteins in the file.

The difference between the `PQRIO` class and the `PDBIO` class is that the `PQRIO` class reads the atomic charges and radii from the PQR file and assigns them to the atoms in the protein.

```python
from protkit.file_io import PQRIO
protein = PQRIO.load("1AHW.pqr")[0]
```

### Saving PQR Files

The PQRIO class provides a `save()` method for saving PQR files. The `save()` method takes a single `Protein` object or a list of `Protein` objects as input.

```python
PQRIO.save(protein, "1AHW.pqr")
```

### 4.5 Working with mmCIF Files

The mmCIF format is a text-based format that is used to store data from macromolecular crystallography experiments. It is the successor to the PDB format and is the preferred format for the PDB archive.

#### Reading mmCIF Files

The mmCIFIO class provides a `load()` method for reading mmCIF files. The `load()` method returns a list of `Protein` objects, describing the proteins in the file.

```python
from protkit.file_io import MMTFIO

protein = MMTFIO.load("1AHW.cif")[0]
```

#### Saving mmCIF Files

The mmCIFIO class provides a `save()` method for saving mmCIF files. The `save()` method takes a single `Protein` object or a list of `Protein` objects as input.

```python
MMTFIO.save(protein, "1AHW.cif")
```

#### Implementation Note

Note that the current implementation relies on converting the mmCIF file to PDB format using BioPython. Metadata that are not supported by the PDB format will be lost in the process. Future implementations of this class will correctly handle additional metadata.

### 4.6 Working with MMTF Files

The MMTF format is a binary format that is used to store data from macromolecular crystallography experiments. It is the successor to the PDB format and is the preferred format for the PDB archive.

#### Reading MMTF Files

The MMTFIO class provides a `load()` method for reading MMTF files. The `load()` method returns a list of `Protein` objects, describing the proteins in the file.

```python
from protkit.file_io import MMTFIO

protein = MMTFIO.load("1AHW.mmtf")[0]
```

#### Implementation Note

Saving MMTF files is not currently supported in Protkit. Future versions of Protkit may include support for saving MMTF files.

Note that the current implementation relies on converting the mmCIF file to PDB format using BioPython. Metadata that are not supported by the PDB format will be lost in the process. Future implementations of this class will correctly handle additional metadata.

---

## 5. Representing Protein Structure

Protkit provides a number of classes for representing protein molecular data in structural format.

Classes for representing structural data include the `Protein`, `Chain`, `Residue` and `Atom` classes. 
These classes are available through the `protkit.structure` package.

### 5.1 Aims with Protkit Data Structure Representations

In designing these classes, we have aimed to make them intuitive and easy to use. We also aimed
to give them a set of desirable properties that make them suitable for use in computational biology (and beyond what is currently available in other libraries).

The following design considerations are taken into account:

- The hierarchy present in proteins should be maintained, i.e. a protein may contain chains, chains may contain residues, and residues may contain atoms.
- It should be easy to filter these structures based on certain criteria, to present a linear view of the data (for example to be used in machine learning applications). For example, it should be possible to easily return the coordinates of all the atoms in the protein.
- The attributes associated with each structure should be extensible, beyond the core attributes associated with each structure. For example, it should be possible to add surface areas or hydrophobicity values to residues. This should be accomplished in a way that does not bloat the class with unnecessary attributes.
- The structures should be easily serializable and deserializable. For example, it should be possible to save and load protein files easily.

### 5.2 Protein, Chain, Residue and Atom Class Hierarchy

The Protein class is the top-level class for storing information about a protein.
It contains a set of Chain objects. Each chain object contains a list of Residue
objects. In turn, each residue object contains a set of Atom objects. The Protein,
Chain, Residue and Atom classes are part of the core data representation of 
molecules and is accessible via the `protkit.structure` package.

At each level of the hierarchy, there are methods for accessing, updating or 
removing the objects at lower levels of the hierarchy. Thus, the Protein class has
methods to operate on Chains, Residues and Atoms; the Chain class has methods to
operate on Residues and Atoms and the Residue class has methods to operate on Atoms.
Objects at lower levels of the hierarchy can also access their parent objects.

### 5.3 Protein Class

The Protein class is the top-level class for storing information about a protein.

#### Loading a Protein

A protein object is typically constructed by loading a protein from a PDB or Prot file.
When a protein is loaded from such a file, the protein object is automatically populated
with chains, residues and atoms.

```python
from protkit.file_io import ProtIO

protein = ProtIO.load("1AHW.prot")[0]
```

#### Accessing Chains in a Protein

The protein class exposes many methods for accessing information about its chains,
residues and atoms. For example, ```.chains``` provides an iterator through all 
the chains of the protein. ```get_chain()``` returns a specific chain. 
```rename_chain()``` allows you to rename a chain, etc.

The following example shows how to iterate through all the chains in a protein, 
returning the number of residues and the sequence of each chain.

```python
for chain in protein.chains:
    print(f"{chain.chain_id}: {chain.num_residues} residues")
    print(chain.sequence)
```

Running the code above would produce the following output:

```commandline
A: 214 residues
DIKMTQSPSSMYASLGERVTITCKASQDIRKYLNWYQQKPWKSPKTLIYYATSLADGVPSRFSGSGSGQDYSLTISSLE
SDDTATYYCLQHGESPYTFGGGTKLEINRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNG
VLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC
B: 214 residues
EIQLQQSGAELVRPGALVKLSCKASGFNIKDYYMHWVKQRPEQGLEWIGLIDPENGNTIYDPKFQGKASITADTSSNTA
YLQLSSLTSEDTAVYYCARDNSYYFDYWGQGTTLTVSSAKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTW
NSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSTWPSETVTCNVAHPASSTKVDKKI
```

To access a specific chain by name, use the `get_chain()` method.

```python
chain_a = protein.get_chain("A")
print(chain_a.sequence)
```

#### Isolating Chains

The `keep_chains()` method allows you to keep only the chains you are interested in.
You might want to make a copy of the protein before doing this, as the method modifies the protein in place.
This is accomplished by calling the `copy()` method on the protein.

```python
protein2 = protein.copy()
protein2.keep_chains(["A", "B"])
print(protein.num_chains)
print(protein2.num_chains)
```

Running this example will illustrate that the original protein had 6 chains, while the new protein has only 2 chains.
A useful and similar method is `remove_chains()`, which removes the chains you are not interested in.

```python
protein3 = protein.copy()
protein3.remove_chains(["C", "D"])
```

#### Renaming Chains

The `rename_chain()` method allows you to rename a chain. Although the new chain
name can be any string, it is recommended to use a single uppercase letter.

```python
protein.rename_chain("A", "Z")
```

#### Accessing Residues in a Protein

The residues that are part of the protein can be accessed in a similar way as the chains. 
For example, the `residues` attribute provides an iterator through all the residues of the protein.
The iterator steps through all the chains in the protein, and then through all the residues in each chain.

```python
for residue in protein.residues:
    print(f"{residue.residue_id}: {residue.name}")
```

#### Accessing Atoms in a Protein

Similary, the atoms that are part of the protein can be accesses through the `atoms` attribute.
For example, to access the coordinates of all the atoms in the protein, you can use the following code:

```python
coordinates = [(atom.x, atom.y, atom.z) for atom in protein.atoms]
print(coordinates)
```

#### Powerful Filters

Protkit provides powerful filters for filtering the atoms, residues and chains of a protein.
The functionality is provided by the `filter_atoms()`, `filter_residues()` and `filter_chains()` methods.
These methods provide iterators through the atoms, residues and chains of the protein, respectively.

For example, suppose we are interested only in the coordinates of the backbone atoms (N, CA, C) of Proline (PRO) 
residues in the A chain of protein. We can do this by applying a filter at the protein level.

In the example, the filter takes in chain, residue and atom criteria, although any criteria that is not applicable can also be skipped. The criteria are specified as a list of tuples. Each tuple contains a field name and a list of values. The filter returns an iterator to all entities that match the criteria.

```python
atoms = protein.filter_atoms(
    chain_criteria=[("chain_id": "A")], 
    residue_criteria=[("residue_type", "PRO")], 
    atom_criteria=[("atom_type", ["C", "CA", "N"])])
for atom in atoms:
    print(f"{atom.residue.residue_code}, {atom.atom_type}, {atom.x}, {atom.y}, {atom.z}")
```

Running this example would produce output similar to the following:

```commandline
PRO8, N, 12.446, 8.496, 43.176
PRO8, CA, 12.683, 9.89, 42.836
PRO8, C, 12.013, 10.129, 41.497
PRO40, N, 8.979, 5.367, 25.606
PRO40, CA, 8.307, 6.502, 24.936
PRO40, C, 7.392, 6.029, 23.842
...
```

Note that the filters provide a convenient way to access the data in a linear fashion, which is useful for machine learning applications.

#### Protein Statistics

Proteins maintain a number of statistics about the protein, such as the number of chains, residues and atoms, etc. Consult the API documentation for more details.

```python
print(protein.num_chains)

print(protein.num_residues)
print(protein.num_disordered_residues)
print(protein.num_hetero_residues)
print(protein.num_water_residues)

print(protein.num_atoms)
print(protein.num_disordered_atoms)
print(protein.num_hetero_atoms)
print(protein.num_heavy_atoms)
print(protein.num_hydrogen_atoms)
```

Chains and residues maintain similar statistics.

### 5.4 Chain Class

The `Chain` class represents a chain in a protein. It contains a list of `Residue` objects.

#### Identifying a Chain

A chain object can be identified by its chain id, which (usually) is a single uppercase letter. The chain id is a unique identifier for the chain in the protein, i.e. no other chains in the protein can have the same identifier.

To access the chain id of a chain, use the `chain_id` attribute.

```python
chain = protein.get_chain("A")
print(chain.chain_id)
```

If you would like to change the chain id of a chain that is part of a protein, 
call the `rename_chain()` method of the protein.

```python
protein.rename_chain("A", "Z")
```

Chains can also exist independently of a protein, although this is not common. In this case,
the chain id can be changed at the chain level:

```python
chain.chain_id = "Z"
```

#### Chain sequence

The sequence of a chain can be accessed through the `sequence` attribute. The sequence is a string of the one-letter codes of the residues in the chain.
Note that the sequence is automatically generated from the residues in the chain, and is not stored as an attribute of the chain. In some cases, such as 
when some residue information is missing, the sequence may not be the true sequence of the chain. In 
such cases, it is preferable to work with `Sequence` objects instead.

```python
print(chain.sequence)
```

#### Referencing the Chain's Parent Protein

A chain object maintains a reference to its parent protein if is part of a protein. This is useful for accessing the protein from the chain.

```python
protein = chain.protein
```

#### Accessing Residues in a Chain

The residues that are part of the chain can be accessed in a similar way as for the protein, through the `residues` attribute, which returns an iterator through all the residues of the chain.


```python
for residue in chain.residues:
    print(f"{residue.residue_id}: {residue.name}")
```

To get a residue by its residue index, use the `get_residue()` method. Note that the method uses 0-based indexing, i.e. to get the first residue, use `get_residue(0)`.

```python
residue = chain.get_residue(0)
```

#### Accessing Atoms in a Chain

Similarly, the atoms that are part of the chain can be accessed through the `atoms` attribute, which returns an iterator through all the atoms of the chain.

```python
for atom in chain.atoms:
    print(f"{atom.atom_id}: {atom.atom_type}")
```

#### Powerful Filters

The `filter_residues()` and `filter_atoms()` methods provide powerful filters for filtering the residues and atoms of a chain. The functionality is similar to the `filter_atoms()` and `filter_residues()` methods of the protein class.

### 5.5 Residue Class

The `Residue` class represents a residue in a protein. It contains a list of `Atom` objects.

#### Identifying a Residue

Within a chain, a residue is identified by its 0-based index. To get the first residue in a chain, use `get_residue(0)`.

```python
residue = chain.get_residue(0)
```

#### Core Residue Properties

Residues have a set of core properties that are maintained as part of the residue.

`residue_type` is the type of the residue, such as "ALA", "GLY", "PRO", etc.  Since the residue types can be residues besides the 20 standard amino acids, the three-letter code is used to identify the residue type. To access the one-letter code of the residue, use the `short_code` attribute.

`sequence_no` and `insertion_code` refers to the position of the residue in the chain, as provided in the source data, or assigned at a later time. 

`is_disordered` is a boolean that indicates whether the residue contains disordered atoms.  `is_hetero` is a boolean that indicates whether the residue is a hetero residue. These properties are derived from the atoms in the residue.

```python
print(residue.residue_type)    # eg. GLY
print(residue.short_code)      # eg. G
print(residue.sequence_no)     # eg. 100
print(residue.insertion_code)  # eg. A
print(residue.is_disordered)   # eg. True
print(residue.is_hetero)       # eg. False
```

Residues also maintain additional properties that make it easy to identify them. For example, `sequence_code` is a combination of the sequence number and the insertion code.  `residue_code` provides more context about the residue, by combining the `residue_type` and `sequence_code`.  `sequence_code` and `residue_code` thus uniquely identifies the residue in the chain.  

`residue_id` combines the chain id of the residue's chain, together with the `residue_code`. The `residue_id` is thus unique within a protein.

```python
print(residue.sequence_code)    # eg. 100A
print(residue.residue_code)     # eg. GLY100A
print(residue.residue_id)       # eg. A:GLY100A
```

#### Referencing the Residue's Parent Chain

A residue can maintain reference to its parent chain. This is useful for accessing the chain from the residue.

```python
chain = residue.chain
```

#### Accessing Atoms in a Residue

The atoms that are part of the residue can be accessed through the `atoms` attribute, which returns an iterator through all the atoms of the residue.

```python
for atom in residue.atoms:
    print(f"{atom.atom_id}: {atom.atom_type}")
```

To get a specific atom in a residue, use the `get_atom()` method. The method takes the atom type as input.  For example, to get the CA atom in a residue, use `get_atom("CA")`.

```python
atom = residue.get_atom("CA")
```

#### Modifying the Atoms in a Residue

The atoms in a residue can be modified in place. For example, to remove all the atoms except backbone atoms, call the `keep_backbone_atoms()` method. This will remove all the atoms in the residue except the backbone atoms (N, CA, C, O).

```python
residue.keep_backbone_atoms()
```

To keep only specific types of atoms, call the `keep_atoms()` method. This method takes a list of atom types as input, and removes all atoms that are not in the list.

```python
residue.keep_atoms(["N", "CA", "C"])
```

To remove only specific types of atoms, call the `remove_atoms()` method. This method takes a list of atom types as input, and removes all atoms that are in the list.

```python
residue.remove_atoms(["OXT"])
```

To remove all hydrogen atoms (deprotonate the residue), call the `remove_hydrogen_atoms()` method.

```python
residue.remove_hydrogen_atoms()
```

Note that most of these methods are available at the chain and protein level as well. When these methods are called at the chain or protein level, they are applied to all the residues in the chain or protein.

#### Powerful Filters

As is the case for the protein and chain classes, the `filter_atoms()` method provides a powerful filter for filtering the atoms of a residue. The functionality is similar to the `filter_atoms()` method of the protein and chain classes.

### 5.6 Atom Class

The `Atom` class represents an atom in a protein. It is the lowest level of the hierarchy of the protein data structure. 
It contains the coordinates of the atom, the atom type, the residue it belongs to, etc.

#### Identifying an Atom

An atom is identified by its atom type, which is a string that identifies the type of the atom, such as "N", "CA", "C", "O", etc. The atom type is unique within a residue, i.e. no two atoms in the same residue can have the same atom type.

To access the atom type of the atom, use the `atom_type` attribute.

```python
atom = residue.get_atom("CA")
print(atom.atom_type)
```

#### Core Atom Properties

Atoms have a set of core properties that are maintained as part of the atom.

`element` is the element of the atom, such as "C", "N", "O", etc.

`x`, `y` and `z` are the coordinates of the atom.

`is_disordered` is a boolean that indicates whether the atom is disordered.  `is_hetero` is a boolean that indicates whether the atom is a hetero atom.

```python
print(atom.element)             # eg. C
print(atom.x, atom.y, atom.z)   # eg. 12.446, 8.496, 43.176
print(atom.is_disordered)       # eg. False
print(atom.is_hetero)           # eg. False
```

#### Fixing Disordered Atoms



### 5.7 Extending the Structure Classes with Properties

The Protein, Chain, Residue and Atom classes can be extended with additional properties. For example, it is possible to add surface areas, hydrophobicity values, etc. to residues. This is accomplished by adding attributes to the objects.

```python
protein = ProtIO.load("1AHW.prot")[0]
protein.get_chain("A").get_residue(1).set_attribute("surface_area", 100.0)
```
---

## 6. Ensuring Data Quality

Working with data from the PDB is challenging because the data is often of poor
quality. For example, data about residues may be missing or incomplete. Atoms 
within residues may be missing or incomplete. The coordinates of atoms may be
uncertain and flagged as disordered. There may be hetero atoms or water molecules
in the structure which may need to be removed before experimentation. The list 
goes on.

Protkit provides a set of methods to ensure that the data is of high quality and
to remove of fix data where possible.

### 6.1 Creating a Copy of a Protein

Fixes or changes to a protein are made in place. To preserve the original protein,
it is recommended that a copy of the protein be made before any changes.

The following code creates a deep copy of the protein 3i40. The copy is identical
to the original protein. Changes made to the copy will not affect the original
protein.

```python
from protkit.file_io import ProtIO

protein = ProtIO.load("3i40.prot")[0]
protein2 = protein.copy()
```

Copies of other structures such as chains, residues and atoms can be made in the 
same way.

### 6.2 Removing Water Molecules

Water molecules are often present in PDB files. They are not part of the protein
structure and are often removed before experimentation.  Waters are usually 
identified by their residue code, which is "HOH". They could be present in any
chain or be within chains of their own.  

6BOM is the crystal structure of mutant 2-methylcitrate synthase mcsAG396A from 
Aspergillus furmigatus. It contains a large number of additional water molecules.
The following code removes all water molecules from the protein.

```python
from protkit.file_io import ProtIO
protein = ProtIO.load("6BOM.prot")[0]
print(f"{protein.num_water_residues}")
protein.remove_water_residues()
print(f"{protein.num_water_residues}")
```

Water molecules can also be removed at the chain level: 

```python
protein.get_chain("A").remove_water_residues()
```

### 6.3 Removing Hetero Residues

Hetero atoms are atoms that are not part of the protein structure. They are
often removed before experimentation. Hetero atoms are usually identified by
their residue code, which is not one of the 20 standard amino acids. They could
be present in any chain or be within chains of their own.

6BOM contains a TRS (2-AMINO-2-HYDROXYMETHYL-PROPANE-1,3-DIOL) molecule in the
A chain. The following code removes all TRS residues from the protein.

```python
from protkit.file_io import ProtIO

protein = ProtIO.load("6BOM.prot")[0]
print(f"{protein.num_hetero_residues}")
protein.remove_hetero_residues("TRS")
print(f"{protein.num_hetero_residues}")
```

The ```remove_hetero_residues()``` method is overloaded to take a string or a list of hetero 
residue codes. Thus, the following code removes all TRS and SO4 residues from the protein.

```python
protein.remove_hetero_residues(["TRS", "SO4"])
```

If no argument is specified, all hetero residues (including waters) are removed.

```remove_hetero_residues()``` can also be called at the chain level.

### 6.4 Fixing Disordered Atoms

Disordered atoms are atoms whose coordinates are uncertain. They are usually 
flagged as disordered in the PDB file. Disordered atoms are usually identified
by their "alt_loc" code, such "A", "B", "C", etc. They typically have an occupancy
value, indicating the probability that the atom is present at the given location.

For each atom, Protkit stores information about the alt_loc code, occupancy value,
coordinates and temperature factor as part of the `Atom` object. 

A crystal structure for human insulin is provided in the structure 3i40. It
contains a number of disordered atoms associated with the TYR14 residue in the
A chain. The following code illustrates that there are 8 disordered atoms
in the structure. The CB atom of the TYR14 residue is disordered, but the CA
atom is not. The CB atom of the TYR14 residue is disordered because it has
two alt_loc codes, "A" and "B". The occupancy value for both alt_loc codes is
0.5, indicating that the probability of the atom being at either location is
equal.

```python
from protkit.file_io import ProtIO

protein = ProtIO.load("3i40.prot")[0]
print(f"{protein.num_disordered_atoms}")
print(protein.get_chain("A").get_residue(13).is_disordered)
print(protein.get_chain("A").get_residue(13).get_atom("CA").is_disordered)
print(protein.get_chain("A").get_residue(13).get_atom("CB").is_disordered)
```

Running the above code would produce the following output:

```commandline
True
False
True
```

Protkit provides a method to fix disordered atoms by choosing the atom with the
highest occupancy value. The following code fixes all disordered atoms in the
protein.

```python
protein.fix_disordered_atoms()
print(protein.get_chain("A").get_residue(13).get_atom("CB").is_disordered)
print(f"{protein.num_disordered_atoms}")
```

Running the above code would produce the following output:

```commandline
False
0
```

```fix_disordered_atoms()``` can also be called at the chain or residue level.

### 6.5 Removing Hydrogen Atoms (Deprotonation)

Hydrogen atoms are often removed before experimentation. The following code
removes all hydrogen atoms from the protein.

```python
from protkit.file_io import ProtIO

protein = ProtIO.load("1A4Y_A_B.prot")[0]

print(f"{protein.num_atoms} atoms")
print(f"{protein.num_heavy_atoms} heavy atoms")
print(f"{protein.num_hydrogen_atoms} hydrogen atoms")

protein.remove_hydrogen_atoms()

print(f"{protein.num_hydrogen_atoms} hydrogen atoms after removal")
```

Running the above code would produce the following output:

```commandline 
8849 atoms
4484 heavy atoms
4365 hydrogen atoms
0 hydrogen atoms after removal
```

```remove_hydrogen_atoms()``` can also be called at the chain or residue level.

***This guide is in active development. More sections will be added soon.***