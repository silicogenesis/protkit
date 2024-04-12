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
7. [Calculating Protein Properties](#7-calculating-protein-properties)
    1. [Methodology behind Calculating Properties](#71-methodology-behind-calculating-properties)
    2. [Hydrophobicity](#72-hydrophobicity)
    3. [Hydrophobicity Class](#73-hydrophobicity-class)
    4. [Mass](#74-mass)
    5. [Chemical Class](#75-chemical-class)
    6. [Residue Charge](#76-residue-charge)
    7. [Residue Polarity](#77-residue-polarity)
    8. [Donors and Acceptors](#78-donors-and-acceptors)
    9. [Surface Area](#79-surface-area)
    10. [Volume](#710-volume)
    11. [Volume Class](#711-volume-class)
    12. [Object Bounds and Center](#712-object-bounds-and-center)
    13. [Bond Lengths](#713-bond-lengths)
    14. [Bond Angles](#714-bond-angles)
    15. [Torsion / Dihedral Angles](#715-torsion-dihedral-angles)
    16. [Circular Variance](#716-circular-variance)
    17. [Interfaces](#717-interfaces)
8. [Tasks, Tools and Pipelines](#8-tasks-tools-and-pipelines)
    1. [Tasks](#81-tasks)
    2. [Tools](#82-tools)
    3. [Pipelines](#83-pipelines)
    4. [Current Tool Integrations](#84-current-tool-integrations)
    
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

In the example, we download a PDB file from the RCSB, extract the A and B chains and do some cleanup like
removing hetero atoms and fixing disordered atoms.  We then compute dihedral angles and surface areas for the
protein and save it to a file.  We then load the protein from the file and print the surface area and a note
that we added to the protein.

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
of the file it will be saved to on disk.

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
 Download.download_pdb_file_from_rcsb("1ahw", "data/pdb/rcsb/1ahw.pdb")

 # Download a PDB file from the SAbDab database and save it to a file.
 Download.download_pdb_file_from_sabdab("1ahw", "data/pdb/sabdab/1ahw.pdb")

 # Download multiple PDB files from the RCSB PDB in parallel and save them to a directory.
 Download.download_pdb_files_from_rcsb(["1ahw", "1a4y", "1a6m"], "data/pdb/rcsb/", n_jobs=3)
```

### 3.2 Downloading Fasta Files

Fasta files can be downloaded from the Uniprot database, as illustrated by the example below. 
The Uniprot id of the file is used as an identifier of the file. The files are downloaded
with a specified filename or directory, as in the previous example.

Fasta files can similarly be downloaded from the RCSB database. In this case,
the PDB id of the file is used as an identifier of the file.

```python
from protkit.download import Download

    # Download a Fasta file from the Uniprot database and save it to a file.
    Download.download_fasta_file_from_uniprot("P01308", "data/fasta/uniprot/P01308.fasta")

    # Download multiple Fasta files from the Uniprot database in parallel and save them to a directory.
    Download.download_fasta_files_from_uniprot(["P01308", "P68871", "P00533"], "data/fasta/uniprot", n_jobs=3)

    # Download a Fasta file from the RCSB PDB database and save it to a file.
    Download.download_fasta_file_from_rcsb("1ahw", "data/fasta/rcsb/1ahw.fasta")

    # Download multiple Fasta files from the RCSB PDB database in parallel and save them to a directory.
    Download.download_fasta_files_from_rcsb(["1ahw", "1a4y", "1a6m"], "data/fasta/rcsb/", n_jobs=3)
```

### 3.3 Downloading CIF and Binary CIF Files

In a similar way to PDB files, CIF files and Binary CIF files can be downloaded from the RCSB PDB database.

```python
from protkit.download import Download

    # Download a CIF file from the RCSB PDB database and save it to a file.
    Download.download_cif_file_from_rcsb("1ahw", "data/cif/rcsb/1ahw.cif")

    # Download multiple CIF files from the RCSB PDB database in parallel and save them to a directory.
    Download.download_cif_files_from_rcsb(["1ahw", "1a4y", "1a6m"], "data/cif/rcsb/")

    # Download a binary CIF file from the RCSB PDB database and save it to a file.
    Download.download_binary_cif_file_from_rcsb("1ahw", "data/cif/rcsb/1ahw.bcif")

    # Download multiple binary CIF files from the RCSB PDB database in parallel and save them to a directory.
    Download.download_binary_cif_files_from_rcsb(["1ahw", "1a4y", "1a6m"], "data/cif/rcsb/")
```

---

## 4. File I/O

Protkit provides a number of classes for reading and writing molecular data files.  These classes are located in 
the `protkit.file_io` package.  For example, the `PDBIO` class provides methods for reading and writing PDB files, 
while the `ProtIO` class provides methods for reading and writing Prot files, etc.

Generally, the classes provide a method to read one or more proteins from a file, using the `load()` method,
and a method to write one or more proteins to a file, using the `save()` method. Conversions between different
file formats are also possible, through the use of the `convert()` method.

The `load()` method returns a list of `Protein` objects, while the `save()` method takes a list of `Protein` objects
as input. We will discuss the `Protein` class in more detail in the next section.

### 4.1 Working with PDB Files

#### Reading PDB Files

The PDB file format is a standard format for storing structural information about proteins. OpenSG provides a simple
interface for reading and writing PDB files.

```python
from protkit.file_io import PDBIO

protein = PDBIO.load("1ahw.pdb")[0]
```

The `load()` method returns a list of proteins. Generally, PDB files generated through X-ray crystallography contain
only one protein, while NMR files may contain multiple proteins. In Protkit, the `load()` method always returns a 
list of proteins, even if the file contains only one protein.   

The above code loads the protein 1AHW from the file 1AHW.pdb. In this case, the list contains only one protein, 
so we take the first element of the list, using the `[0]` index in the array.

The `Protein` object contains all the information about the protein. The `Protein` class is discussed in more detail 
in the next section.

#### Saving PDB Files

Saving a protein to a PDB file is accomplished using the `save()` method. Let's modify the 
protein by extracting the A and B chains and saving it to a new file.

```python
protein.keep_chains(["A", "B"])

PDBIO.save(protein, "1ahw_AB.pdb")
```

#### PDB Metadata

PDB files contain metadata such as the resolution of the structure, the method used to solve the structure,
the date the structure was deposited etc.

Protkit does not currently read such metadata. Since Protkit is primarily aimed towards structural biology, it reads 
data such as the sequence, the coordinates of the atoms, the occupancy and temperature factors, etc. contained
in records such as ATOM, HETATM, MODEL, SEQRES, etc.
   
Future versions of Protkit may include methods to read and write such metadata.

### 4.2 Working with Prot Files

Although the PDB file format is the dominant format for storing structural information about proteins, it is an
old format that has many limitations.  For example, it cannot store additional information about the
protein, such as surface areas, molecular masses or hydrophobicity values.

To overcome the limitations, we have developed the Prot file format, which is a simple (JSON) text-based format that 
can store all the information about a protein in a single file. The Prot file format is designed to be easy to read
and write, and to be easily extended.  The Prot file format closely resembles the `Protein` data structure 
(discussed later), and can be efficiently serialised and deserialised.

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

protein = ProtIO.convert("1ahw.pdb", "1ahw.prot")
```

The Prot file could also have been generated by first loading the protein from
the PDB file and then saving it to a Prot file.

#### Reading Prot Files

The ProtIO class provides a `load()` method for reading Prot files, in the same 
way as the PDBIO class provides a load method for reading PDB files.

```python
protein = ProtIO.load("1ahw.prot")[0]
```

#### Saving Prot files

The ProtIO class provides a `save()` method for saving Prot files, as in the example:

```python
ProtIO.save(protein, "1ahw.prot")
```

#### Saving Multiple Proteins in a Single Prot File

Note that multiple proteins can be saved to a single Prot file. This facilitates
the creation of databases of proteins that can easily be managed through a
single file.

```python
from protkit.file_io import ProtIO, PDBIO

protein1 = PDBIO.load("1ahw.pdb")[0]
protein2 = PDBIO.load("4nkq.pdb")[0]

ProtIO.save([protein1, protein2], "database.prot")
``` 

#### Compression and Decompression of Prot Files

Prot files are based on the JSON format, which is a text-based format. This
could lead to large file sizes. To reduce the file size, Prot files are
compressed using the gzip algorithm. Compression can be enabled or disabled
by parameters passed to the save and load functions.  Compression is enabled
by default.

```python
from protkit.file_io import ProtIO

protein = ProtIO.load("1ahw.prot", decompress=True)[0]
ProtIO.save(protein, "1ahw.prot.json", compress=False)
```

### 4.3 Working with Fasta Files

The Fasta file format is a standard format for storing sequence information about proteins. Protkit provides 
a simple interface for reading and writing Fasta files.

#### Reading Fasta Files

The FastaIO class provides a `load()` method for reading Fasta files. The `load()` method returns 
a list of `Sequence` objects, describing the sequences in the file. The sequence objects contain the
sequence and the header of the sequence. The `Sequence` class is discussed in more detail later in this guide.

```python
from protkit.file_io import FastaIO

sequences = FastaIO.load("1ahw.fasta")
for sequence in sequences:
    print(f"{sequence.chain_id}: {sequence}")
```

#### Saving Fasta Files

The FastaIO class provides a `save()` method for saving Fasta files. The `save()` method takes a single `Sequence` 
object or a list of `Sequence` objects as input.

```python
from protkit.file_io import FastaIO

sequences = FastaIO.load("1AHW.fasta")
FastaIO.save(sequences, "1AHW_copy.fasta")
FastaIO.save(sequences[0], "1AHW_A.fasta")
```

### 4.4 Working with PQR Files

Protkit supports the PQR file format, which is a variant of the PDB file format that includes atomic charges and 
radii. The PQR file format is used in molecular dynamics simulations and in the calculation of electrostatic 
potentials.

It is typically created by programs such as APBS, PDB2PQR, etc.

#### Reading PQR Files

The PQRIO class provides a `load()` method for reading PQR files. The `load()` method returns a list of `Protein` objects,
describing the proteins in the file.

The difference between the `PQRIO` class and the `PDBIO` class is that the `PQRIO` class reads the atomic charges 
and radii from the PQR file and assigns them to the atoms in the protein.

```python
from protkit.file_io import PQRIO
    
protein = PQRIO.load("1ahw.pqr")[0]
```

### Saving PQR Files

The PQRIO class provides a `save()` method for saving PQR files. The `save()` method takes a single `Protein` object 
or a list of `Protein` objects as input.

```python
PQRIO.save(protein, "1ahw_copy.pqr")
```

### 4.5 Working with mmCIF Files

The mmCIF format is a text-based format that is used to store data from macromolecular crystallography experiments. 
It is the successor to the PDB format and is the preferred format for the PDB archive.

#### Reading mmCIF Files

The mmCIFIO class provides a `load()` method for reading mmCIF files. The `load()` method returns a list of `Protein` 
objects, describing the proteins in the file.

```python
from protkit.file_io import MMTFIO

protein = MMTFIO.load("1AHW.cif")[0]
```

#### Saving mmCIF Files

The mmCIFIO class provides a `save()` method for saving mmCIF files. The `save()` method takes a 
single `Protein` object or a list of `Protein` objects as input.

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

protein = ProtIO.load("1ahw.prot")[0]
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
print(protein3.num_chains)
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
from protkit.file_io import ProtIO

protein = ProtIO.load("1ahw.prot")[0]

for residue in protein.residues:
    print(f"{residue.id}: {residue.residue_type}")
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
    chain_criteria=[("chain_id", "A")], 
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
from protkit.file_io import ProtIO

protein = ProtIO.load("1ahw.prot")[0]
chain = protein.get_chain("A")
print(chain.chain_id)
```

If you would like to change the chain id of a chain that is part of a protein, 
call the `rename_chain()` method of the protein.

```python
protein.rename_chain("A", "Z")
print(chain.chain_id)
```

Chains can also exist independently of a protein, although this is not common. In this case,
the chain id can be changed at the chain level:

```python
chain.chain_id = "X"
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
    print(f"{residue.id}: {residue.residue_type}")
```

To get a residue by its residue index, use the `get_residue()` method. Note that the method uses 0-based indexing, i.e. to get the first residue, use `get_residue(0)`.

```python
residue = chain.get_residue(0)
print(residue.id)
```

#### Accessing Atoms in a Chain

Similarly, the atoms that are part of the chain can be accessed through the `atoms` attribute, which returns an iterator through all the atoms of the chain.

```python
for atom in chain.atoms:
    print(f"{atom.id}")
```

#### Powerful Filters

The `filter_residues()` and `filter_atoms()` methods provide powerful filters for filtering the residues and atoms of a chain. The functionality is similar to the `filter_atoms()` and `filter_residues()` methods of the protein class.

### 5.5 Residue Class

The `Residue` class represents a residue in a protein. It contains a list of `Atom` objects.

#### Identifying a Residue

Within a chain, a residue is identified by its 0-based index. To get the first residue in a chain, use `get_residue(0)`.

```python
from protkit.file_io import ProtIO

protein = ProtIO.load("1ahw.prot")[0]
chain = protein.get_chain("A")
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

`id` combines the chain id of the residue's chain, together with the `residue_code`. The `id` is thus unique within a protein.

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
    print(f"{atom.id}")
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

To access the atom type of the atom, use the `atom_type` attribute.  To identify the full atom, use the `id` attribute.

```python
from protkit.file_io import ProtIO

protein = ProtIO.load("1ahw.prot")[0]
residue = protein.get_chain("A").get_residue(0)
atom = residue.get_atom("CA")

print(atom.id)
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

Sometimes, additional properties such as charge and temperature factor are maintained as part of the atom. 
It may or may not be present in PDB files.  These properties will be stored internally
within the atom class, only if they are present in the PDB file. The can be accessed 
through the `get_attribute()` method.

```python
print(atom.get_attribute("assigned_charge"))
print(atom.get_attribute("temp_factor"))
```

#### Referencing the Atom's Parent Residue

An atom maintains a reference to its parent residue.

```python
residue = atom.residue
```

#### Fixing Disordered Atoms

Sometimes, the crystal structure of the protein may contain disordered atoms. 
Disordered atoms are atoms whose coordinates are uncertain. They are usually 
flagged as disordered in the PDB file. Disordered atoms are identified by their 
"alt_loc" code, such as "A", "B", "C", etc. They typically have an occupancy 
value, indicating the probability that the atom is present at the given location.

Atoms, maintain information about the alt_loc code, occupancy value, coordinates 
and temperature factor as part of the `Atom` object. Multiple such records may be
present in the PDB file for a single atom.  In Protkit, each atom maintains a list of
such records. Note that this is different to other libraries, which may create
an atom for each record.

To create only a single entry for each atom, call the `fix_disordered_atoms()` method.
This could be called at the protein, chain or residue level.

```python
from protkit.file_io import ProtIO

protein = ProtIO.load("1ahw.prot")[0]
protein.fix_disordered_atoms()
```

### 5.7 Extending the Structure Classes with Properties

The Protein, Chain, Residue and Atom classes can be extended with additional 
properties. For example, it is possible to add surface areas, 
hydrophobicity values, etc. to proteins, chains, residues or atoms.

The ```set_attribute()``` and ```get_attribute()``` methods can be used to set
and get attributes of proteins, chains, residues and atoms.  The attributes are
stored as key-value pairs.  The key is a string and the value can be any object.
The ```has_attribute()``` method can be used to check if an attribute exists.  
The ```delete_attribute()``` method can be used to remove an attribute.

The rationale behind this design is that it allows for the extension of the
Protein, Chain, Residue and Atom classes with additional properties, without
bloating the classes with unnecessary attributes. It also means that the 
properties can be computed once and is then conveniently stored with the
structure.

Various property calculators can add different calculated properties to 
the structure classes in a convenient way. How this is achieved will be discussed
in upcoming sections.

The following example shows how to manage attributes of a protein.

```python
from protkit.file_io import ProtIO

protein = ProtIO.load("3i40.prot")[0]

protein.set_attribute("note", "The file describes Human Insulin")
protein.set_attribute("organism", "Homo Sapiens")
protein.set_attribute("resolution", 1.85)
protein.get_chain("A").set_attribute("name", "Insulin A chain")
protein.get_chain("A").get_residue(0).set_attribute('first', True)

print(protein.get_attribute("note"))
print(protein.get_attribute("organism"))
print(protein.get_attribute("resolution"))
print(protein.get_chain("A").has_attribute("name"))
print(protein.get_chain("A").get_attribute("name"))
print(protein.get_chain("B").has_attribute("name"))
print(protein.get_chain("A").get_residue(0).get_attribute('first'))
```

Running the above example produces the following output:

```python
The file describes Human Insulin
Homo Sapiens
1.85
True
Insulin A chain
False
True
```

Attributes added to objects in this way will be persisted when the object is
saved to a Prot file. Of course, files such as PDB files do not support such
attributes, and they will be lost when the object is saved to a PDB file.

Note that values that conform to elementary types such as int, float, str, bool, 
or simple data structures such as lists, tuples, sets, dictionaries, etc. can be
stored and will be persisted correctly.  More complex objects that are instantiated
as classes that are set as attributes will not be persisted correctly. 
This will be fixed in a future release.

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

protein = ProtIO.load("1a4y_A_B.prot")[0]

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

## 7. Calculating Protein Properties

It is often useful to compute properties, attributes or features of a protein, chain,
residue or atom. For example, we may want to compute the surface area 
of a protein, the hydrophobicity of a residue or the mass of an atom.

Protkit provides an extensive set of methods for computing such properties.
These methods are available through the `protkit.properties` package.

### 7.1 Methodology behind Calculating Properties

The design decision was made to compute the properties externally to the
protein, chain, residue or atom objects. This was done to avoid bloating the
code with a large number of methods inside these classes.

Instead, properties are calculated using classes from the `protkit.properties`
package. These classes are designed to be easily extendable, and to be easily
used in a variety of contexts.  Each class within the `protkit.properties` package
is used to compute one or a few related properties.  For example, the `SurfaceArea`
class is used to compute the surface area of a protein, while the `Hydrophobicity`
class is used to compute the hydrophobicity of a residue.

These classes are typically called with a protein, chain, residue or atom as input,
and return the computed property as output. The function names of the methods 
are typically descriptive of the property being computed. For example, the 
`surface_area_of_protein()` method of the `SurfaceArea` class computes the surface
area of a protein, while the `hydrophobicity_of_chain()` method of the `Hydrophobicity`
class computes the hydrophobicity of a chain.

The computed property is not stored 
as part of the corresponding object, but is returned as a separate object.  This
is done to avoid bloating the objects with unnecessary attributes, in cases they are 
not needed downstream. The computed property can be stored as part of the object,
by setting the `assign_attribute` parameter to `True` when calling the method. An 
optional `key` parameter can be used to specify the name of the attribute.
Properties computed in this way are stored as part of the object and are available for 
future use.  It is persisted in the Prot file when the object is saved.

Note that some properties are computed at several levels of the hierarchy. For example,
when the mass of a protein is calculated, the mass of the protein is computed as the sum
of the masses of the chains in the protein, which in turn is the sum of the masses of the
residues in the chain. When the `assign_attribute` parameter is set to `True`, the mass
values are assigned across all levels of the hierarchy.

### 7.2 Hydrophobicity

- Applicable to: `Protein`, `Chain`, `Residue` and `Sequence`
- Default key: `hydrophobicity`
- Property Type: `Additive`

Hydrophobicity is computed using the ```Hydrophobicity``` class in the
```src.properties``` package. 

Hydrophobicity is computed using
the Kyte-Doolittle scale, which ranges from -4.5 to 4.5.  The hydrophobicity
values are only defined for the 20 standard amino acids.  The values for 
non-standard amino acids are set to 0.0.

The methods `hydrophobicity_of_protein()`, `hydrophobicity_of_chain()`, `hydrophobicity_of_residue()` and
`hydrophobicity_of_sequence()` are used to compute the hydrophobicity. Hydrophobicity is not defined
at the atom level.

The example below illustrates how hydrophobicity is computed at the protein level
and assigned as an attribute of the protein (which in turn assigns the hydrophobicity to
the chains and residues in the protein).

```python
from protkit.file_io import ProtIO
from protkit.properties import Hydrophobicity

protein = ProtIO.load("1ahw.prot")[0]
Hydrophobicity.hydrophobicity_of_protein(protein, assign_attribute=True)

print(protein.get_attribute("hydrophobicity"))
print(protein.get_chain("A").get_attribute("hydrophobicity"))
print(protein.get_chain("B").get_attribute("hydrophobicity"))
print(protein.get_chain("A").get_residue(0).get_attribute("hydrophobicity"))
```

### 7.3 Hydrophobicity Class

- Applicable to: `Protein`, `Chain`, `Residue` and `Sequence`
- Default key: `hydrophobicity_class`
- Property type: `Categorical`

Hydrophobicity class is computed using the ```Hydrophobicity``` class in the
```src.properties.hydrophobicity``` module. 

The hydrophobicity class is a categorical property that assigns a residue to one of
three classes: hydrophobic, hydrophilic or neutral. The hydrophobicity classes are
assigned based on the definitions by the IMGT.

```python
from protkit.file_io import ProtIO
from protkit.properties import Hydrophobicity

protein = ProtIO.load("1ahw.prot")[0]
Hydrophobicity.hydrophobicity_classes_of_protein(protein, assign_attribute=True)
print(Hydrophobicity.HYDROPHOBICITY_CLASS_STRING[protein.get_chain("A").get_residue(0).get_attribute("hydrophobicity_class")])
```

### 7.4 Mass

Mass can be computed in three ways:

- Atomic mass: A mass is assigned to each atom.
  It is then propagated to the residue, chain and protein levels.
  The first method is applicable if the complete atomic structure is available.
  - Applicable to: `Protein`, `Chain`, `Residue` and `Atom`
  - Default key: `atomic_mass`
  - Property type: `Additive`
- Residue mass. A mass is assigned to each residue.
  It is then propagated to the chain and protein. It can also be assigned to a sequence,
  since it is dependent on knowledge of the residue sequence only.
  - Applicable to: `Protein`, `Chain`, `Residue` and `Sequence`
  - Default key: `residue_mass`
  - Property type: `Additive`
- Molecular mass. In the third method, the molecular mass is assigned to each residue.
  This is the mass of the unbound amino acid. It would not make sense to
  assign this mass to the chain or protein. It can be assigned to a sequence.
  - Applicable to: `Residue` and `Sequence` (can be called at the protein or chain level, but no assignments are made)
  - Default key: `molecular_mass`
  - Property type: `Normal`

```python
from protkit.file_io import ProtIO
from protkit.properties import Mass

protein = ProtIO.load("1ahw.prot")[0]
Mass.residue_mass_of_protein(protein, assign_attribute=True)
Mass.atomic_mass_of_protein(protein, assign_attribute=True)
Mass.molecular_mass_of_protein(protein, assign_attribute=True)

residue = protein.get_chain("A").get_residue(0)
print(residue.get_attribute("residue_mass"))
print(residue.get_attribute("atomic_mass"))
print(residue.get_attribute("molecular_mass"))
```

### 7.5 Chemical Class

- Applicable to: `Protein`, `Chain`, `Residue` and `Sequence`
- Default key: `chemical_class`
- Property type: `Categorical`

The chemical class of a residue is computed using the ```ChemicalClass``` class in the
```src.properties``` module.

The chemical class is a categorical property that assigns a residue to one of
seven classes: aliphatic, aromatic, acidic, amide, basic, hydroxyl and sulfur. 
The chemical classes follow the definitions by the IMGT.

```python
from protkit.file_io import ProtIO
from protkit.properties import ChemicalClass

protein = ProtIO.load("1ahw.prot")[0]
ChemicalClass.chemical_classes_of_protein(protein, assign_attribute=True)
print(ChemicalClass.CHEMICAL_CLASS_STRING[protein.get_chain("A").get_residue(0).get_attribute("chemical_class")])
```

### 7.6 Residue Charge

- Applicable to: `Protein`, `Chain`, `Residue` and `Sequence`
- Default key: `charge`
- Property type: `Additive`

The charge of a residue, chain or protein is computed using the ```Charge``` class in the
```src.properties``` module. 

Note that these charges are not calculated using the pKa values of the residues. 
Instead, they are based on residue charges as defined by the IMGT.
The charges are then propagated from residues to chains and proteins.

```python
from protkit.file_io import ProtIO
from protkit.properties import Charge

protein = ProtIO.load("1ahw.prot")[0]
Charge.charge_of_protein(protein, assign_attribute=True)
chain = protein.get_chain("A")
residue = chain.get_residue(0)
print(protein.get_attribute("charge"))
print(chain.get_attribute("charge"))
print(residue.get_attribute("charge"))
```

### 7.7 Residue Polarity

- Applicable to: `Protein`, `Chain`, `Residue` and `Sequence`
- Default key: `polarity`
- Property type: `Categorical`

The polarity of a residue is computed using the ```Polarity``` class in the
```src.properties``` module.

The polarity is a categorical property that assigns a residue to one of
two classes: polar or non-polar. The polarity classes are assigned based on the
definitions by the IMGT.

```python
from protkit.file_io import ProtIO
from protkit.properties import Polarity

protein = ProtIO.load("1ahw.prot")[0]
Polarity.polarities_of_protein(protein, assign_attribute=True)
residue = protein.get_chain("A").get_residue(0)
print(Polarity.POLARITY_STRING[residue.get_attribute("polarity")])
```

### 7.8 Donors and Acceptors

Donors and Acceptors are calculated using the ```DonorsAcceptors``` class in the
```src.properties``` module.

Donors and Acceptors can refer to atoms or residues. Specific atoms are assigned
as donors or acceptors based on their atom type. Residues are assigned as donors
or acceptors based on the presence of specific atoms in the residue. 

Donor and acceptor values are assigned based on the definitions by the IMGT. 

- Donor and Acceptor Residues
  - Applicable to: `Protein`, `Chain`, `Residue` and `Sequence`
  - Default keys: `is_donor_residue` and `is_acceptor_residue`
  - Property type: `Boolean`
- Donor and Acceptor Atoms
  - Applicable to: `Protein`, `Chain`, `Residue` and `Atom`
  - Default keys: `is_donor_atom` and `is_acceptor_atom`
  - Property type: `Boolean`

```python
from protkit.file_io import ProtIO
from protkit.properties import DonorsAcceptors

protein = ProtIO.load("1ahw.prot")[0]

DonorsAcceptors.donor_residues_of_protein(protein, assign_attribute=True)
DonorsAcceptors.acceptor_residues_of_protein(protein, assign_attribute=True)
residue = protein.get_chain("A").get_residue(0)
print(residue.get_attribute("is_donor_residue"))
print(residue.get_attribute("is_acceptor_residue"))

DonorsAcceptors.donor_atoms_of_protein(protein, assign_attribute=True)
DonorsAcceptors.acceptor_atoms_of_protein(protein, assign_attribute=True)
atom = protein.get_chain("A").get_residue(0).get_atom("N")
print(atom.get_attribute("is_donor_atom"))
print(atom.get_attribute("is_acceptor_atom"))
```

### 7.9 Surface Area

- Applicable to: `Protein`, `Chain`, `Residue` and `Atom`
- Default key: `surface_area`
- Property type: `Additive`

The surface area of a protein, chain, residue or atom is computed using the ```SurfaceArea``` class in the
```src.properties``` module.

The surface area refers to the solvent accessible surface area (SASA) of the protein, chain, residue or atom. 
The surface area can be computed using the Lee-Richards or the Shrake-Rupley algorithm. 
The surface area is then propagated from atoms to residues, chains and proteins.

The current implementation uses FreeSASA for surface computations. FreeSASA is installed
as one of the dependencies when Protkit is installed.


```python
from protkit.file_io import ProtIO
from protkit.properties import SurfaceArea

protein = ProtIO.load("1ahw.prot")[0]
chain = protein.get_chain("A")
residue = chain.get_residue(0)

SurfaceArea.surface_area_of_protein(protein, assign_attribute=True)
print(protein.get_attribute("surface_area"))
print(chain.get_attribute("surface_area"))
print(residue.get_attribute("surface_area"))
```

### 7.10 Volume

- Applicable to: `Protein`, `Chain`, `Residue` and `Sequence`.
- Default key: `volume`
- Property type: `Additive`

The volume of a protein, chain, residue or sequence is computed using the ```Volume``` class in the
```src.properties``` module.

Note that volumes are based on average estimated values for the residues. The volume of a protein is the sum of the volumes of the chains in the protein, which in turn is the sum of the volumes of the residues in the chain. 
The volume of a sequence is the sum of the volumes of the residues in the sequence. The volume
values are defined by the IMGT for the 20 standard amino acids.

```python
from protkit.file_io import ProtIO, FastaIO
from protkit.properties import Volume

protein = ProtIO.load("1ahw.prot")[0]
chain = protein.get_chain("A")
residue = chain.get_residue(0)

Volume.volume_of_protein(protein, assign_attribute=True)
print(protein.get_attribute("volume"))
print(chain.get_attribute("volume"))
print(residue.get_attribute("volume"))

sequence = FastaIO.load("1ahw.fasta")[0]
Volume.volume_of_sequence(sequence, assign_attribute=True)
print(sequence.get_attribute("volume"))
```

### 7.11 Volume Class

- Applicable to: `Protein`, `Chain`, `Residue` and `Sequence`
- Default key: `volume_class`
- Property type: `Categorical`

The volume class of a residue is computed using the ```Volume``` class in the
```src.properties``` module.

The volume class is a categorical property that assigns a residue to one of five classes: 
very small, small, medium, large or very large. The volume classes are assigned based on the
definitions by the IMGT.

```python
from protkit.file_io import ProtIO, FastaIO
from protkit.properties import Volume

protein = ProtIO.load("1ahw.prot")[0]
residue = protein.get_chain("A").get_residue(0)

Volume.volume_classes_of_protein(protein, assign_attribute=True)
print(Volume.VOLUME_CLASS_STRING[residue.get_attribute("volume_class")])

sequence = FastaIO.load("1ahw.fasta")[0]
Volume.volume_classes_of_sequence(sequence, assign_attribute=True)
for volume_class in sequence.get_attribute("volume_class"):
    print(Volume.VOLUME_CLASS_STRING[volume_class])
```

### 7.12 Object Bounds and Center

- Applicable to: `Protein`, `Chain` and `Residue`
- Default keys: `bounds` and `center`
- Property type: `Normal`

The bounds and center of a protein, chain or residue are computed using the 
```Bounds``` class in the ```src.properties``` package.

The bounds of an object are the minimum and maximum coordinates of the object
in the x, y and z dimensions, based on the coordinates of the atoms that
comprise the object. The center of an object is the average of the minimum and
maximum coordinates of the object in the x, y and z dimensions.

```python
from protkit.file_io import ProtIO
from protkit.properties import Bounds

protein = ProtIO.load("3i40.prot")[0]

Bounds.bounds_of_protein(protein, assign_attribute=True)
Bounds.center_of_protein(protein, assign_attribute=True)

print(protein.get_attribute("bounds"))
print(protein.get_chain("A").get_attribute("bounds"))
print(protein.get_attribute("center"))
print(protein.get_chain("A").get_residue(0).get_attribute("center"))
```

### 7.13 Bond Lengths

- Applicable to: `Protein`, `Chain` and `Residue`
- Default key: `bond_lengths`
- Property type: `Normal`

Bond lengths for residues can be computed using the
```BondLengths``` class in the ```src.properties``` package.  

Each different type of residue  has a unique set expected covalent bonds.  
This collection of bond lengths is defined by the `HEAVY_ATOM_BONDS` dictionary
in the `BondLengths` class.  For example, an Alanine (ALA) residue has 4 covalent
bonds, namely N-CA, CA-C, C-O and CA-CB.

Bond lengths are computed for each residue and can be stored as an attribute of the
residue.  The bond lengths are stored as a dictionary, where the keys are tuples of
atom types and the values are the bond lengths. In the case that one of the
expected atoms are not present in the residue, the bond length is set to `None`.

```python
from protkit.file_io import ProtIO
from protkit.properties import BondLengths

protein = ProtIO.load("1ahw.prot")[0]
residue = protein.get_chain("A").get_residue(0)

BondLengths.bond_lengths_of_protein(protein, assign_attribute=True)
for (atom1, atom2), length in residue.get_attribute("bond_lengths").items():
    print(f"{atom1}-{atom2}: {length}")
```

To calculate the bond lengths (peptide lengths) between residues in a chain, use the
`bond_lengths_of_chain()` or `bond_lengths_of_protein()` method.  This will calculate
the bond lengths between the C atom of one residue and the N atom of the next residue.

If the `assign_attribute` parameter is set to `True, the bond lengths are stored as 
a list of floating point values in the `bond_lengths` attribute of the corresponding
chain(s).

```python
from protkit.file_io import ProtIO
from protkit.properties import BondLengths

protein = ProtIO.load("1ahw.prot")[0]
chain = protein.get_chain("A")

BondLengths.peptide_bond_lengths_of_protein(protein, assign_attribute=True)
print(chain.get_attribute("peptide_bond_lengths"))
```

### 7.14 Bond Angles

- Applicable to: `Protein`, `Chain` and `Residue`
- Default key: `bond_angles`
- Property type: `Normal`

Bond angles for residues can be computed using the
```BondAngles``` class in the ```src.properties``` package.

Bond angles are defined as the angle between three atoms in a residue.  Each 
residue has a unique set of expected bond angles.  This collection of bond angles
is defined by the `HEAVY_ATOM_ANGLES` dictionary in the `BondAngles` class.  For
example, an Alanine (ALA) residue has 4 bond angles, namely N-CA-C, CA-C-O, 
N-CA-CB and C-CA-CB.

To compute the bond angles of a residue, use the `bond_angles_of_residue()`,
`bond_angles_of_chain()` or `bond_angles_of_protein()` method.  This will calculate
the bond angles for each residue in the chain, or for each residue in the protein.

```python
from protkit.file_io import ProtIO
from protkit.properties import BondAngles

protein = ProtIO.load("1ahw.prot")[0]
residue = protein.get_chain("A").get_residue(0)

BondAngles.bond_angles_of_protein(protein, assign_attribute=True)
for (atom1, atom2, atom3), angle in residue.get_attribute("bond_angles").items():
    print(f"{atom1}-{atom2}-{atom3}: {angle}")
```

### 7.15 Torsion / Dihedral Angles

- Applicable to: `Protein`, `Chain` and `Residue`
- Default key: `dihedral_angles`
- Attribute type: `Normal`

Dihedral angles can be computed using the
```DihedralAngles``` class in the ```src.properties``` package.

Dihedral angles are defined as the angle between four atoms in a residue.  Each
residue has a unique set of expected dihedral angles.  This collection of dihedral
angles is defined by the `DIHEDRAL_ANGLES` dictionary in the `DihedralAngles`
class.  For example, an Aspargine (ASN) residue has 2 dihedral angles associated
with its side chain, namely N-CA-CB-CG and CA-CB-CG-OD1. These are known as 
the chi1 and chi2 angles.  

In addition, the phi, psi and omega angles are also
calculated - these are defined as the dihedral angles between two neighbouring
residues. As is the usual convention, the phi, psi and omega angles are 
associated with the second residue.

To compute the dihedral angles of a residue, use the `dihedral_angles_of_residue()`,
`dihedral_angles_of_chain()` or `dihedral_angles_of_protein()` method.  This will calculate
the dihedral angles for each residue in the chain, or for each residue in the protein.

```python
from protkit.file_io import ProtIO
from protkit.properties import DihedralAngles

protein = ProtIO.load("1ahw.prot")[0]
residue = protein.get_chain("A").get_residue(1)

DihedralAngles.dihedral_angles_of_protein(protein, assign_attribute=True)
for angle_name, angle in residue.get_attribute("dihedral_angles").items():
    print(f"{angle_name}: {angle}")
```

### 7.16 Circular Variance

- Applicable to: `Residue` and `Atom`
- Default key: `cv_residue` and `cv_atom`
- Property type: `Normal`

The circular variance of a residue or atom is computed using the ```CircularVariance``` class in the
```src.properties``` package.

Circular variance (CV) is a measure of the spread of atoms in the local environment of a residue or atom.
It can be calculated in two ways: 

```python
from protkit.file_io import ProtIO
from protkit.properties import CircularVariance

protein = ProtIO.load("1ahw.prot")[0]
residue = protein.get_chain("A").get_residue(1)
atom = residue.get_atom("CA")

CircularVariance.circular_variance_by_residue(protein, assign_attribute=True)
print(residue.get_attribute("cv_residue"))

CircularVariance.circular_variance_by_atom(protein, assign_attribute=True)
print(atom.get_attribute("cv_atom"))
```

### 7.17 Interfaces

- Applicable to: `Residue` and `Atom`
- Default key: `in_interface`
- Property type: `Boolean`

The interface of a residue or atom is computed using the ```Interface``` class in the
```src.properties``` package.

The interface is a boolean property that indicates whether the residue or atom 
is part of an interface. The interface is defined as the region of the protein
that is in contact with another protein or ligand. The interface is computed
using the distance between atoms in the residue or atom and atoms in other
proteins or ligands.

Interfaces can be computed at the residue or atom level.  For each, two sets of 
objects are provided as input: the objects (atoms or residues) of interest and a set of other objects
that are in contact with the object of interest.  

For atoms, two atoms are considered to be interacting if the distance between them is 
less than a specified `cutoff` value. 

```python
from protkit.file_io import ProtIO
from protkit.properties import Interface

protein = ProtIO.load("1ahw.prot")[0]
atoms1 = list(protein.filter_atoms(chain_criteria=[("chain_id", ["A", "B"])]))
atoms2 = list(protein.filter_atoms(chain_criteria=[("chain_id", ["C"])]))

Interface.interface_atoms(atoms1, atoms2, cutoff=5.0, assign_attribute=True)
for atom in atoms1:
    if atom.get_attribute("in_interface"):
        print(atom.id)
```

For residues, two residues are considered to be
interacting if the distance between any two atoms in the two residues is less than the
specified `cutoff` value. A second method for computing interface residues is to
consider two residues to be interacting if the distance between the C-alpha atoms of
the two residues is less than the specified `cutoff` value.

```python
from protkit.file_io import ProtIO
from protkit.properties import Interface

protein = ProtIO.load("1ahw.prot")[0]
residues1 = list(protein.filter_residues(chain_criteria=[("chain_id", ["A", "B"])]))
residues2 = list(protein.filter_residues(chain_criteria=[("chain_id", ["C"])]))

Interface.interface_residues(residues1, residues2, cutoff=5.0, assign_attribute=True)
Interface.interface_residues_from_alpha_carbon(residues1, residues2, cutoff=6.0, assign_attribute=True, key="ca_in_interface")

for residue in residues1:
    if residue.get_attribute("in_interface"):
        print(residue.id)

for residue in residues1:
    if residue.get_attribute("ca_in_interface"):
        print(residue.id)
```

## 8. Tasks, Tools and Pipelines

Protkit distinguishes between Tasks, Tools and Pipelines.

### 8.1 Tasks

In computational biology, we often need to perform various tasks on a protein
structure or sequence. These tasks can be fairly diverse.  

A ***Task*** is a collection of steps that perform a specific function. For example, a 
***structure-prediction task*** may take a protein sequence as input and predict the
structure of the protein as the output of the task.

Examples of tasks could be calculating the surface area of a protein, docking two
proteins together, predicting the structure upon mutation of a residue, predicting
the structure of a protein from sequence alone, etc.

Note that in Protkit, a task does not directly perform the computation (computation is performed
by tools - discussed next). Instead, a task provides a definition of the computation
to be performed. For example, in a structure-prediction task, a transformation of
a protein sequences (```Sequence``` object) to a protein structure (```Protein``` object)
is defined. The actual computation is performed by a tool.

Tasks are defined in the `protkit.tasks` package. Each task is specified as a abstract
class that defines the input and output of the task. Tools that implement the task
must adhere to the input and output specifications of the task.

### 8.2 Tools

***Tools*** implement the functionality of tasks. For example, a ***structure-prediction task***
may be implemented by a ***structure-prediction tool*** such as AlphaFold or RoseTTAFold.
Modelling the effect of a mutation on the structure of a protein (task) may be implemented by
a ***mutation-modelling tools*** such as FoldX or EvoEF (tools).

The computation performed by a tool can be achieved in a variety of ways. For example, the
computation may be performed by calling an external executable program, 
by calling a web service, by calling a library function installed via PyPI, etc.

Tools are defined in the `protkit.tools` package. Each tool is specified as a class that
implements the functionality of a task. The tool must adhere to the input and output
specifications of the task. We often call the classes associated with tools adaptors, to
indicate that they provide a "bridge" between the task specification and the actual computation.
For example, the `AlphaFoldAdaptor` class is an adaptor for the `AlphaFold` tool, which
implements the `StructurePrediction` task.

### 8.3 Pipelines

***Pipelines*** are a sequence of tasks that are executed in a specific order and 
automates the execution of a series of tasks. The output of one task may be fed to the input
of another task in the pipeline. For example, a ***structure-prediction pipeline*** may
consist of a ***sequence-alignment task***, a ***structure-prediction task*** and a ***structure-refinement
task***.

Pipelines are defined in the `protkit.pipelines` package. Each pipeline is specified as a class
that defines the sequence of tasks to be executed. The pipeline also defines the input and output
of the pipeline.

Pipelines provide the ability to automate complex tasks. Pipelines can be executed in a single
command, and the output of one task is automatically fed to the input of the next task.  The tools
used in the pipeline are interchangeable, as long as they adhere to the input and output specifications
defined by the tasks.

### 8.4 Current Tool Integrations

#### 8.4.1 Reduce

The Reduce software is a widely used tool for adding hydrogen atoms to protein 
structures. It is used to prepare protein structures for molecular dynamics 
simulations, docking studies, and other computational analyses.

The Reduce software is not included in the Protkit package. It must be installed
separately. You can download Reduce and follow installation instructions from the following website: 
https://github.com/rlabduke/reduce

The `ReduceAdaptor` class in the `protkit.tools.reduce_adaptor` module provides an
interface to the Reduce software.  The following example illustrates how easy 
it is to use the `ReduceAdaptor` class to protonate or deprotonate a protein structure.

```python
from protkit.tools.reduce_adaptor import ReduceAdaptor
from protkit.file_io import ProtIO

reduce_bin_path = "/usr/local/bin/reduce"
reduce = ReduceAdaptor(reduce_bin_path)

protein = ProtIO.load("1ahw.prot")[0]
protein_protonated = reduce.protonate(protein)
protein_deprotonated = reduce.deprotonate(protein)
```

The `RecuceAdaptor` implements the `Protonation` task, specified in the `protkit.tasks.protonation`
module. The `Protonation` task specifies that the input is a `Protein` object and the output is a `Protein`
object. The `ReduceAdaptor` class adheres to these specifications.

Note that a protonated or deprotonated `Protein` object is returned by the `protonate()` and `deprotonate()`
methods. The advantage for researchers is that they do not need to worry about the details of how the
protonation or deprotonation is performed.

#### 8.4.2 FreeSASA

FreeSASA is a widely used tool for calculating the solvent accessible surface area (SASA) of proteins.
It is used to calculate the surface area of proteins, which is useful for understanding protein-protein
interactions, protein-ligand interactions, protein folding, etc.

FreeSASA is available as a Python library. By default, Protkit will install FreeSASA as a dependency
when it is installed. 

The `FreeSASAAdaptor` class in the `protkit.tools.freesasa_adaptor` module provides 
an interface to the FreeSASA software. The following example illustrates how the surface
area for each atom in a protein can be calculated using the `FreeSASAAdaptor` class.

```python
from protkit.tools.freesasa_adaptor import FreeSASAAdaptor
from protkit.file_io import ProtIO

freesasa = FreeSASAAdaptor()
protein = ProtIO.load("1ahw.prot")[0]

atoms = list(protein.atoms)
atom_areas = freesasa.calculate_surface_area(atoms)
protein_area = sum(atom_areas)
```

The `FreeSASAAdaptor` can be initialised with a number of parameters, such as the probe radius (default 1.4A),
the algorithm used to calculate the surface area (default Lee-Richards, but Shrake-Rupley can also be specified)
and various other parameters that control the calculation of the surface area.

Internally, Protkit also uses FreeSASA to calculate the surface area of proteins. The
`SurfaceArea` class in the `protkit.properties.surface_area` is build on top of the `FreeSASAAdaptor`
class and provides a different interface for assigning the surface area to the protein, chain, residue
or atom objects.

### 8.4.3 PDB2PQR

The pdb2pqr software is a widely used tool for adding charges and radii to protein
structures. It is used to prepare protein structures for molecular dynamics
simulations, docking studies, and other computational analyses.  pdb2pqr can add a
limited number of missing heavy atoms to structures, as well as
hydrogen atoms. It may change the coordinates of some atoms in the structure.

pdb2pqr is available as a Python library. By default, Protkit will install pdb2pqr 
as a dependency when it is installed.

The `PDB2PQRAdaptor` class in the `protkit.tools.pdb2pqr_adaptor` module provides an
interface to the pdb2pqr software. The use of the `PDB2PQRAdaptor` class is similar
to the `ReduceAdaptor` class. The following example illustrates how to use the
`PDB2PQRAdaptor` class to add charges and radii to a protein structure.

```python
from protkit.tools.pdb2pqr_adaptor import PDB2PQRAdaptor
from protkit.file_io import ProtIO

pdb2pqr = PDB2PQRAdaptor(force_field=PDB2PQRAdaptor.AMBER)

protein = ProtIO.load("1ahw.prot")[0]
protein_out = pdb2pqr.run(protein)
```

The `PDB2PQRAdaptor` can be initialised with a number of parameters, such as the force field.
The pdb2pqr software implements a number of force fields, such as AMBER, CHARMM, PARSE, etc.
If no force field is specified, PARSE is used.

The `run` command is used to run the pdb2pqr software on a `Protein` object. It 
returns a new `Protein` object with charges and radii added to the atoms. The run
command takes additional parameters, such the paths for storing pdb or pqr files
that are used as input and output for the pdb2pqr software.

***This guide is in active development. More sections will be added soon.***