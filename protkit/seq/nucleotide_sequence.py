#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `NucleotideSequence` to represent a nucleotide sequence.

IUPAC base symbols for nucleotide sequences:
    A   A               Adenine
    C   C               Cytosine
    G   G               Guanine
    T   T               Thymine
    U   U               Uracil

    W   A or T          Weak
    S   G or C          Strong
    R   A or G          Purine
    Y   C or T          Pyrimidine
    K   G or T          Ketone
    M   A or C          Amino

    B   C or G or T     Not A
    D   A or G or T     Not C
    H   A or C or T     Not G
    V   A or C or G     Not T

    N   A, C, G or T    Any base
    -                   Gap

See: https://en.wikipedia.org/wiki/Nucleic_acid_notation
"""

from typing import Optional, List, Union
from protkit.seq.sequence import Sequence


class NucleotideSequence(Sequence):
    def __init__(self,
                 sequence: Union[str, List[str]],
                 description: Optional[str] = None,
                 chain_id: Optional[str] = None
                 ):
        super().__init__(sequence, description, chain_id)
