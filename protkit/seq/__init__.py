"""
Package `seq` contains classes to represent sequences in computational biology.

Sequence Representations:

- `Sequence` class in `sequence`: Represents sequence data.
- `ProteinSequence` class in `protein_sequence`: Represents sequence data of a protein.
- `AntibodySequence` class in `antibody_sequence`: Represents sequence data of an antibody.
- `NucleotideSequence` class in `nucleotide_sequence`: Represents sequence data of a nucleotide.
"""

from .sequence import Sequence
from .protein_sequence import ProteinSequence
from .antibody_sequence import AntibodySequence
from .nucleotide_sequence import NucleotideSequence

from .sequence_alignment import SequenceAlignment
