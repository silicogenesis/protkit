#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `AntibodySequence` to represent an antibody sequence.
"""

from typing import Optional, List, Union
from protkit.seq.protein_sequence import ProteinSequence


class AntibodySequence(ProteinSequence):
    def __init__(self,
                 sequence: Union[str, List[str]],
                 description: Optional[str] = None,
                 chain_id: Optional[str] = None
                 ):
        """
        Constructor.

        Args:
            sequence (Union[str, List[str]]): The sequence of residue names.
            description (Optional[str]): Optional description associated with the sequence.
            chain_id (Optional[str]): Optional chain ID associated with the sequence.

        Returns:
            None
        """
        super().__init__(sequence, description, chain_id)
