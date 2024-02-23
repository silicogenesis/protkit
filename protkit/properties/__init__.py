"""
Package `properties` contains classes to calculate properties of proteins.

Properties are divided into various categories:

- Chemical Identity: Properties that describe the chemical composition of molecules.
- Physiochemical: Properties that describe the physical and chemical characteristics of molecules.
- Pharmacophore: Properties that relate to the recognition and binding of a molecule to its biological target.
- Structural: Properties related to the three-dimensional structure of a molecule.
"""

from .bond_angles import BondAngles
from .bond_lengths import BondLengths
from .bounds import Bounds
from .charge import Charge
from .chemical_class import ChemicalClass
from .charge import Charge
from .circular_variance import CircularVariance
from .dihedral_angles import DihedralAngles
from .donors_acceptors import DonorsAcceptors
from .hydrophobicity import Hydrophobicity
from .interface import Interface
from .mass import Mass
from .polarity import Polarity
from .surface_area import SurfaceArea
from .vdw_radius import VDWRadius
from .volume import Volume
