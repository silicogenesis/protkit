# Protkit

Protkit is a Python library for that can be used for a variety of tasks
in computational biology and bioinformatics, focusing and structural bioinformations
and protein engineering.

It is used to process the structure of proteins and protein complexes as well
as their sequences.

It is designed to support the broad community of computational biologists,
bioinformaticians, and machine learning researchers in academia, industry,
and government labs.

Protkit is an open source library that is free to use and modify.  We welcome
contributions from the community.

---

## A Unified Approach to Computational Biology

### End-to-end Computational Biology

Protkit can be used for a variety of computational biology tasks across the computational biology pipeline such as:

- **Reading and writing data** from popular structure file formats, such as
    PDB and MMTF; and sequence file formats, such as FASTA.
- **Downloading** data from popular databases of protein structures, such as the PDB RCSB.
- **Data structures** for representing proteins, protein complexes, chains,
    residues, atoms and sequences. These data structure provides capabilities to extract data
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
- **Metrics** for comparing proteins, such as RMSD and TM-score.
- **Featurization** of proteins and their properties and preparation of datasets
    for **machine learning** applications.
- Performing and enabling a large variety of **computational tasks** on proteins,
    such as protein folding, protein docking, protein-protein binding affinity prediction,
    humanisation of antibodies, prediction of developability characteristics. Care is taken
    that the various tools are interoperable and can be used together in a seamless manner.

### Why Protkit?

Protkit is designed to be easy to use, modular, scalable and extensible.
Here are some of the reasons why you might want to use Protkit:

- **Open Source**: Protkit is open source.  It is free to use and modify.  It is not
    a proprietary system where you are locked in and unable to make changes.
- **Community**: Protkit is supported by a community of computational biologists,
    bioinformaticians, and machine learning researchers. We actively encourage
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
that combines biology, computer science, statistics, mathematics, and engineering.
It is a field that is driven by data.  The amount of data that is being generated
in the field is growing exponentially.  This is due to the rapid advances in
experimental techniques, such as cryo-electron microscopy, that are used to
generate data.

- **Access to data**.  Protkit provides access and downloading of data from
    popular databases of protein structures and sequences, such as the PDB RCSB
    and UnitProt. It supports reading and writing of data stored in a variety of file formats,
    such as PDB and Fasta.

- **Ensuring data quality**. Protein structures are often incomplete or contain
    missing atoms, missing residues, sequence gaps, atomic clashes, alternate
    conformations, hetero residues or water molecules. Unfortunately, we often
    encounter these anomalies in protein structures in published datasets, which
    is then used for training machine learning models or as starting points for
    computational tasks.  Protkit provides methods to detect and fix these anomalies.

**Data Representations**

Unfortunately, there is no unified data representation for protein structures across
various research groups and industry. New tools and applications are often build from
scratch and are based on different data representations, which are often incompatible
with each other. This makes it difficult to share data between different tools and applications
and many tools revert to using PDB files to share data between applications. PDB files are very
limited in the data they can represent and are not suitable for many applications.

- **Unified data representations**. Protkit provides a unified data representation for protein structures.  This data
    representation is based on a hierarchical data structure that can represent the
    structure of proteins, protein complexes, chains, residues, atoms and sequences. It
    provides capabilities to extract data in both hierarchical and linear formats.

    Our hope is that this data representation will be adopted by the community and will
    become a standard for representing protein structures.  This will enable interoperability
    between different tools and applications and will make it easier to share data between
    different tools and applications. Protkit surpasses the capabilities of frameworks such
    as Biopython (a limited hierarchical structure) or BioPandas (a linear view of the data) and
    provides a rich set of methods for extracting and filtering data from the data structure.

- **Extensible data representations and metadata management**.  Protkit provides an extensible data representation
    for protein structures.  It is easy to add new properties to the data structure.  This
    makes it easy to add new functionality to the library. For example, you can easily add
    new properties to the data structure such as hydrophobicity, charge, surface areas, etc.

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
    new properties.  Over time, as the community grows, we hope that various modules will be
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
    all have task definitions in the form of an API that need to be adhered to. Any tool that
    adheres to the API can be used with other tools in the library.  This makes it easy to
    combine tools together to perform complex tasks. We are working on adaptors for various
    tools that will allow them to be used within the Protkit ecosystem.

- **Modular design**. Protkit is modular in design. We often see code repositories associated
    with publications that are monolithic in nature.  Unfortunately, this makes it difficult
    to work with those tools as they make inherent assumptions about the data and the
    computational infrastructure that is available.  Protkit is designed to be modular and
    interoperable.  New modules can be added to the library with ease.  As researchers adopts
    the framework, it allows them to focus on the development of new algorithms and methods,
    rather than having to worry about how these tools will be combined together.

**Machine Learning**

We are seeing rapid advances in machine learning and deep learning applied to computational
biology. Unfortunately, the way in which datasets are prepared for machine learning applications
often leaves a lot to be desired. For example, datasets are often prepared in an ad-hoc manner
and are not reproducible. In other cases, the datasets are not balanced and are biased towards a
particular families of proteins that are overrepresented in the population.

- **Dataset creation**. Protkit provides a rich set of methods for creating datasets for machine
    learning applications.  These methods are designed to be reproducible and
    extensible.  We are taking care to ensure that tools are build into the process to ensure
    that datasets are balanced and are not biased towards a particular families of proteins.

- **Support for machine learning frameworks**. Protkit provides support for different
    machine learning frameworks, such as PyTorch and Tensorflow.  We are working on
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