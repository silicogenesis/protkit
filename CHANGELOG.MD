# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] - 2024-05-08

### Added

- Added the ability to export proteins, chains, residues, atoms and sequences as Pandas Dataframes, through the ProteinToDataFrame class.
- Added additional properties to proteins, chains, residues and atoms, to make querying and filtering easier.

## [0.2.4] - 2024-04-23

### Added
- Support for Propka tool (https://propka.readthedocs.io/en/latest/).
- Allows lookup of residues sequence number and insertion code.

### Fixes.
- Fixed an issue where PDB files with HETATM records were not being saved correctly.

## [0.2.3] - 2024-04-11

### Added
- Support for PDB2PQR tool.
- HETATM support in PDB files.
- Atom class supports initialization with data from PQR files.

### Changed
- TER lines terminate chains in PDB and PQR files.

### Fixes
- PQRIO correctly reads and writes PQR files.
- Fixed an error that can occur when fixing disordered atoms in PDB files contain alt_loc information with occupancy set to 1.00.
- Fixes a circular import issue that can occur when calculating surface properties.

## [0.2.2] - 2024-04-01

### Added
- Support for Reduce tool.
- Support for FreeSASA tool.

### Fixes
- Bug fixes in downloading.
- Bug fixes in Fasta parsing.

## [0.2.1] - 2024-03-08

### Added
- Quick Start Guide examples.

### Fixes
- Slight fixes in bounds and center computation.

## [0.2.0] = 2024-02-26

### Added
- Property computation (Surface Area, Volume, Hydrophobicity, Dihedral Angles, etc.)
- Metrics
- Task and Tool definitions.

## [0.1.0] - 2024-02-20 

### Added
- Initial release of the project. 
- Data representations (Protein, Chain, Residue, Atom, Sequence, etc)
- File IO (PDB, Fasta, Prot, etc.)
- File Download (RCSB, Sabdab, etc.)




