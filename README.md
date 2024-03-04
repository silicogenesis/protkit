<p align="center">
    <!--img src="media/logo/protkit800px.png" width="400"-->
    <img src="https://raw.githubusercontent.com/silicogenesis/protkit/main/media/logo/protkit800px.png" width="400px">
</p>

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

You can clone the repository and install it from source:

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
# are available for later retrieval.
protein.set_attribute("note", "Experimenting with Protkit")
ProtIO.save(protein, "1ahw.prot")
protein2 = ProtIO.load("1ahw.prot")[0]
print(protein2.get_attribute('surface_area'))
print(protein2.get_attribute('note'))
```

Please consult the [Quick Start Guide](QUICK_START_GUIDE.md) for more examples.

---

## A Unified Approach to Computational Biology

### Why Protkit?

Protkit is designed to be easy to use, modular, scalable and extensible.
Here are some of the reasons why you might want to use Protkit:

- **Open Source**: Protkit is open source.  It is free to use and modify.  It is not
    a proprietary system where you are locked in and unable to make changes.
- **Community**: Protkit is supported by a community of computational biologists,
    bioinformaticians and machine learning researchers. We actively encourage
    contributions from the community.
- **Easy to use**: Protkit is designed to be easy to use.  It is intuitive and has a
    consistent API.  It is easy to learn and use.
- **Modular**: Protkit is modular.  It is designed to be used as a library, where you
    can import only the modules that you need.  It is not a monolithic framework.
- **Extensible**: Protkit is extensible.  It is easy to add new functionality to the
    library.
- **Interoperable**: Protkit modules are interoperable.  We have designed the library
    so that modules can be used together in a seamless manner. For example, you can
    easily add docking engines or protein folding engines to the library. While adhering
    to the same API, these engines can be used with other modules in the library.
- **Scalable**: Protkit is scalable.  It can be used to create large datasets and
    perform large-scale analyses on biological data. We designed the library to run
    efficiently on both single machines and clusters.

### Challenges in Structural Biology (and how Protkit can help)

**Data**

Computational biology is a rapidly growing field.  It is a multidisciplinary field
that combines biology, computer science, statistics, mathematics and engineering.
It is a field that is driven by data.  The amount of data that is being generated
in the field is growing exponentially.  This is due to the rapid advances in
experimental techniques, such as cryo-electron microscopy, that are used to
generate data, as well as in-silico data generation such as sequence to structure prediction.

- **Access to data**.  Protkit makes it easy to  downloading protein data from
    popular protein structures and sequence databases, such as PDB RCSB
    and UnitProt. It supports reading and writing of data stored in a variety of file formats
    such as PDB and Fasta.

- **Ensuring data quality**. Protein structures are often incomplete or contain
    missing atoms, missing residues, sequence gaps, atomic clashes, alternate
    conformations, hetero residues or water molecules. Unfortunately, we often
    encounter these anomalies in protein structures in published datasets, which
    is then used for training machine learning models or as starting points for
    computational tasks.  Protkit provides methods to detect and fix these anomalies.

**Data Representations**

Unfortunately, there is no unified data representation for protein structures across
various research groups and industry. New tools and applications are often built from
scratch and are based on different data representations which are often incompatible
with each other. This makes it difficult to share data between different tools and applications
and many tools revert to using PDB files to share data between applications. PDB files are very
limited in the data they can represent and are not suitable for many applications.

- **Unified data representations**. Protkit provides a unified data representation for protein structures.  This data
    representation is based on a hierarchical data structure that can represent the
    structure of proteins, protein complexes, chains, residues, atoms and sequences. It
    provides capabilities to extract data in both hierarchical and linear formats.

    Our goal is that this data representation will be adopted by the community and will
    become a standard for representing protein structures.  This will enable interoperability
    between different tools and applications and will make it easier to share data between
    different tools and applications. Protkit surpasses the capabilities of frameworks such
    as Biopython (a limited hierarchical structure) or BioPandas (a linear view of the data) and
    provides a rich set of methods for extracting and filtering data from the data structure.

- **Extensible data representations and metadata management**.  Protkit provides an extensible data representation framework
    for protein structures.  It is easy to add new properties to the data structure.  This
    makes it easy to add new functionality to the library. For example, you can easily add
    new properties to the data structure such as hydrophobicity, charge, surface areas etc.

    Protkit serializes the data structure when data is stored to disk (in .prot files), meaning that the data
    structure can be easily stored and shared between different tools and applications, preserving
    properties in the process.

- **Property computation**.  We often find that different tools and applications compute the same
    properties in different ways. For example, protein-protein-interfaces are often computed
    using cutoff distances that are different between different tools and applications.  In other
    cases, the change in surface accessible area is used to compute the interface.  This makes it
    difficult to compare results between different tools and applications.  Protkit provides a
    unified way to compute properties, which makes it easier to compare results between different
    tools and applications.

    Protkit provides a rich set of methods for computing properties
    of proteins, such as hydrophobicity, charge, surface areas, secondary structures, dihedral
    angles, interface residues and more.  These properties can be added to the data structure
    and can be used to filter and extract data from the data structure.

    Protkit was designed in a modular way, so that it is easy to add new modules to compute
    new properties.  Over time, as the community grows, modules will be
    added to the library to compute new properties.

**Algorithms and Methods**

The field is also driven by the development of new algorithms and methods, especially
with the rapid advances in machine learning. These
methods are used to analyse the data and extract meaningful insights from it.
These methods are also used to develop new drugs and therapies for diseases.

Unfortunately, many of these methods are not interoperable with each other.  For example,
one docking tool may require a specific format to specify residues that take part or do not
take part in the docking process.  Another docking tool may require a different format. This
makes it difficult to use these tools together in a seamless manner.  Similarly, there may
be difficulties in using the output of one tool as input to another tool.

- **Interoperability of tools**. Protkit was designed to allow tools to be used together in a
    seamless manner. Common tasks that may be performed on proteins, such as folding, docking,
    binding affinity prediction, humanisation of antibodies, prediction of developability, etc.
    all have task definitions in the form of an API that is to be adhered to by the community. Any tool that
    adheres to the API can seamlessly be used with other tools in the library.  This makes it easy to
    combine tools together to perform complex tasks. We are working on adaptors for various
    tools that will allow them to be used within the Protkit ecosystem.

- **Modular design**. Protkit is modular in design. We often see code repositories associated
    with publications that are monolithic in nature.  Unfortunately, this makes it difficult
    to work with those tools as they make inherent assumptions about the data and the
    computational infrastructure that is available.  Protkit is designed to be modular and
    interoperable.  New modules can be added to the library with ease.  As researchers adopt
    the framework it allows them to focus on the development of new algorithms and methods,
    rather than having to worry about how these tools will be combined.

**Machine Learning**

We are seeing rapid advances in machine learning and deep learning applied to computational
biology. Unfortunately, the way in which datasets are prepared for machine learning applications
often leaves a lot to be desired. For example, datasets are often prepared in an ad-hoc manner
and are not reproducible. In other cases, the datasets are not balanced and are biased towards
particular families of proteins that are overrepresented in the population.

- **Dataset creation**. Protkit provides a rich set of methods for creating datasets for machine
    learning applications.  These methods are designed to be reproducible and
    extensible.  We are taking care to ensure that tools are built into the process to ensure
    that datasets are balanced and are not biased towards particular families of proteins.

- **Support for machine learning frameworks**. Protkit provides support for different
    machine learning frameworks such as PyTorch and Tensorflow.  We are working on
    building dataloaders for these frameworks that can be used to load datasets into
    these frameworks across a wide variety of machine learning models.

- **Metrics and evaluation**. Protkit provides a rich set of metrics for evaluating machine
    learning models.  We often see that metrics are not reported in publications, or that
    metrics are computed slightly differently between publications.  Protkit provides a
    consistent set of metrics that can be used to evaluate machine learning models.

- **Published datasets**. We are working on creating a repository of datasets that can be used by the
    community.  These datasets will be created using the methods in Protkit and will be
    reproducible and extensible.

**Computational Infrastructure**

Datasets are often large and require computational infrastructure to process them.  Some
computational processes can be very computationally intensive and require high-performance
computing clusters to process them.  Unfortunately, many researchers do not have access to
such infrastructure.

- **Flexible compute architectures**. Protkit is designed to run on a variety of compute
    architectures, including multi-core CPUs, GPUs, and clusters.

- **Scalable infrastructure**. Protkit is designed to be scalable.  It can be used to
    create large datasets and perform large-scale analyses on biological data. We designed
    the library to run efficiently on both single machines and clusters.

- **Cloud-based computational infrastructure**
    We are working on providing cloud-based access to computational infrastructure
    to the community.

---

## Contributing to Protkit

Protkit is an open source project and we welcome contributions from the community.  If you would like to contribute, please see the [Contibuting](CONTRIBUTING.md) file for details. Please also adhere to the [Code of Conduct](CODE_OF_CONDUCT.md).

---

## Acknowledgements

Protkit was conceived and developed by the scientists and engineers at [Silicogenesis](https://www.silicogenesis.com). SilicoGenesis is a company that is dedicated to the development of computational tools for the life sciences.  We are grateful to the community for their support and contributions to the project. 

We would like to extend our thanks to the following people:

- Fred Senekal
- Lionel Bisschoff
- Mechiel Nieuwoudt
- Claudio Jardim
- Dean Sherry

If you use Protkit in a scientific publication we would appreciate using the following citation:

F. Senekal, L. Bisschoff, M. Nieuwoudt, C. Jardim, D. Sherry, Protkit: A
unified toolkit for protein engineering., 2024. URL: https://protkit.silicogenesis.com/.

Bibtex entry:
```latex
@misc{protkit,
  author = {Senekal, Fred and Bisschoff, Lionel and Nieuwoudt, Mechiel and Jardim, Claudio and Sherry, Dean},
  title = {Protkit: A Unified Toolkit for Protein Engineering.},
  url ={https://protkit.silicogenesis.com/},
  year = {2024}
}
```

---

## Licence

Protkit is licensed under the GPL v3.0 license.  See the [License](LICENSE) file for details.



