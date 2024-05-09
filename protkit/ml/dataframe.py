#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `ProteinToDataframe` to convert a protein structure into a pandas DataFrame.
"""

import pandas as pd
from typing import List

from protkit.structure import Protein
from protkit.seq import Sequence


class ProteinToDataframe:
    DEFAULT_PROTEIN_FIELDS = ['pdb_id', 'num_chains', 'num_residues', 'num_atoms']
    EXCLUDED_PROTEIN_FIELDS = ['chains']

    DEFAULT_CHAIN_FIELDS = ['pdb_id', 'chain_id', 'num_residues', 'num_atoms']
    EXCLUDED_CHAIN_FIELDS = ['residues', 'protein', 'seqres']

    DEFAULT_RESIDUE_FIELDS = ['pdb_id', 'chain_id', 'sequence_no', 'insertion_code', 'residue_type',
                              'is_disordered', 'is_hetero', 'num_atoms', 'num_heavy_atoms', 'num_hydrogen_atoms']
    EXCLUDED_RESIDUE_FIELDS = ['atoms', 'chain']

    DEFAULT_ATOM_FIELDS = ['pdb_id', 'chain_id', 'sequence_no', 'insertion_code', 'residue_type',
                             'atom_type', 'element', 'x', 'y', 'z', 'is_disordered', 'is_hetero',
                            'alt_loc', 'occupancy', 'temp_factor', 'assigned_charge', 'calculated_charge', 'radius']
    EXCLUDED_ATOM_FIELDS = ['residue']

    DEFAULT_SEQUENCE_FIELDS = ['chain_id', 'description', 'length', 'sequence']
    EXCLUDED_SEQUENCE_FIELDS = ['extended_attributes']

    def __init__(self):
        pass

    def proteins_dataframe(self,
                           protein: [Protein, List[Protein]],
                           additional_fields: List[str] = None) -> pd.DataFrame:
        """
        Featurize the protein(s).

        Args:
            protein: Protein object or list of Protein objects.
            additional_fields: List of additional fields to include in the dataframe.
                If None, all default fields as well as any additional fields will be included.
                If List, all default fields as well as the specified fields will be included.

        Returns:
            pd.DataFrame: Dataframe containing protein features.
        """
        if not isinstance(protein, list):
            protein = [protein]

        protein_data = []

        for p in protein:
            protein_attributes = {
                'pdb_id': p.pdb_id,
                'num_chains': p.num_chains,
                'num_residues': p.num_residues,
                'num_atoms': p.num_atoms
            }
            if additional_fields is not None:
                for field in additional_fields:
                    if p.has_attribute(field):
                        protein_attributes[field] = p.get_attribute(field)
                    else:
                        protein_attributes[field] = None
            else:
                for key, value in p.__dict__.items():
                    key = key[1:] if key.startswith('_') else key
                    if key not in ProteinToDataframe.DEFAULT_PROTEIN_FIELDS:
                        if key not in ProteinToDataframe.EXCLUDED_PROTEIN_FIELDS:
                            protein_attributes[key] = value

            protein_data.append(protein_attributes)

        df = pd.DataFrame(protein_data)

        return df

    def chains_dataframe(self,
                         protein: [Protein, List[Protein]],
                         additional_fields: List[str] = None) -> pd.DataFrame:
        """
        Featurize the chain(s) of a protein or list of proteins.

        Args:
            protein: Protein object or list of Protein objects.
            additional_fields: List of additional fields to include in the dataframe.
                If None, all default fields as well as any additional fields will be included.
                If List, all default fields as well as the specified fields will be included.

        Returns:
            pd.DataFrame: Dataframe containing chain features.
        """
        if not isinstance(protein, list):
            protein = [protein]

        chain_data = []

        for p in protein:
            for chain in p.chains:
                chain_attributes = {
                    'pdb_id': chain.pdb_id,
                    'chain_id': chain.chain_id,
                    'num_residues': chain.num_residues,
                    'num_atoms': chain.num_atoms
                }
                if additional_fields is not None:
                    for field in additional_fields:
                        if chain.has_attribute(field):
                            chain_attributes[field] = chain.get_attribute(field)
                        else:
                            chain_attributes[field] = None
                else:
                    for key, value in chain.__dict__.items():
                        key = key[1:] if key.startswith('_') else key
                        if key not in ProteinToDataframe.DEFAULT_CHAIN_FIELDS:
                            if key not in ProteinToDataframe.EXCLUDED_CHAIN_FIELDS:
                                chain_attributes[key] = value

                chain_data.append(chain_attributes)

        df = pd.DataFrame(chain_data)

        return df

    def residues_dataframe(self,
                           protein: [Protein, List[Protein]],
                           additional_fields: List[str] = None) -> pd.DataFrame:
        """
        Featurize the residues of a protein or list of proteins.

        Args:
            protein: Protein object or list of Protein objects.
            additional_fields: List of additional fields to include in the dataframe.
                If None, all default fields as well as any additional fields will be included.
                If List, all default fields as well as the specified fields will be included.

        Returns:
            pd.DataFrame: Dataframe containing residue features.
        """
        residue_data = []

        if not isinstance(protein, list):
            protein = [protein]

        for p in protein:
            for residue in p.residues:
                residue_attributes = {
                    'pdb_id': residue.pdb_id,
                    'chain_id': residue.chain_id,
                    'sequence_no': residue.sequence_no,
                    'insertion_code': residue.insertion_code,
                    'residue_type': residue.residue_type,

                    'is_disordered': residue.is_disordered,
                    'is_hetero': residue.is_hetero,

                    'num_atoms': residue.num_atoms,
                    'num_heavy_atoms': residue.num_heavy_atoms,
                    'num_hydrogen_atoms': residue.num_hydrogen_atoms
                }
                if additional_fields is not None:
                    for field in additional_fields:
                        if residue.has_attribute(field):
                            residue_attributes[field] = residue.get_attribute(field)
                        else:
                            residue_attributes[field] = None
                else:
                    for key, value in residue.__dict__.items():
                        key = key[1:] if key.startswith('_') else key
                        if key not in ProteinToDataframe.DEFAULT_RESIDUE_FIELDS:
                            if key not in ProteinToDataframe.EXCLUDED_RESIDUE_FIELDS:
                                residue_attributes[key] = value

                residue_data.append(residue_attributes)


        df = pd.DataFrame(residue_data)

        return df

    def atoms_dataframe(self,
                        protein: [Protein, List[Protein]],
                        additional_fields: List[str] = None) -> pd.DataFrame:
        """
        Featurize the atoms of a protein or list of proteins.

        Args:
            protein: Protein object or list of Protein objects.
            additional_fields: List of additional fields to include in the dataframe.
                If None, all default fields as well as any additional fields will be included.
                If List, all default fields as well as the specified fields will be included.

        Returns:
            pd.DataFrame: Dataframe containing atom features.
        """
        if not isinstance(protein, list):
            protein = [protein]

        atom_data = []

        for p in protein:
            for atom in p.atoms:
                atom_attributes = {
                    'pdb_id': atom.pdb_id,
                    'chain_id': atom.chain_id,
                    'sequence_no': atom.residue.sequence_no,
                    'insertion_code': atom.residue.insertion_code,
                    'residue_type': atom.residue.residue_type,

                    'atom_type': atom.atom_type,
                    'element': atom.element,

                    'x': atom.x,
                    'y': atom.y,
                    'z': atom.z,

                    'is_disordered': atom.is_disordered,
                    'is_hetero': atom.residue.is_hetero,

                    'alt_loc': atom.get_attribute('alt_loc'),
                    'occupancy': atom.get_attribute('occupancy'),
                    'temp_factor': atom.get_attribute('temp_factor'),
                    'assigned_charge': atom.get_attribute('assigned_charge'),
                    'calculated_charge': atom.get_attribute('calculated_charge'),
                    'radius': atom.get_attribute('radius'),
                }
                if additional_fields is not None:
                    for field in additional_fields:
                        if atom.has_attribute(field):
                            atom_attributes[field] = atom.get_attribute(field)
                        else:
                            atom_attributes[field] = None
                else:
                    for key, value in atom.__dict__.items():
                        key = key[1:] if key.startswith('_') else key
                        if key not in ProteinToDataframe.DEFAULT_ATOM_FIELDS:
                            if key not in ProteinToDataframe.EXCLUDED_ATOM_FIELDS:
                                atom_attributes[key] = value

                atom_data.append(atom_attributes)

        df = pd.DataFrame(atom_data)

        return df

    def sequences_dataframe(self,
                            sequences: List[Sequence],
                            additional_fields: List[str] = None) -> pd.DataFrame:
        """
        Featurize the sequences of a protein or list of proteins.

        Args:
            sequences: List of Sequence objects.
            additional_fields: List of additional fields to include in the dataframe.

        Returns:
            pd.DataFrame: Dataframe containing sequence features.
        """
        sequence_data = []

        for seq in sequences:
            sequence_attributes = {
                'chain_id': seq.chain_id,
                'description': seq.description,
                'length': seq.length,
                'sequence': str(seq)
            }
            if additional_fields is not None:
                for field in additional_fields:
                    if seq.has_attribute(field):
                        sequence_attributes[field] = seq.get_attribute(field)
                    else:
                        sequence_attributes[field] = None
            else:
                for key, value in seq.__dict__.items():
                    key = key[1:] if key.startswith('_') else key
                    if key not in ProteinToDataframe.DEFAULT_SEQUENCE_FIELDS:
                        if key not in ProteinToDataframe.EXCLUDED_SEQUENCE_FIELDS:
                            sequence_attributes[key] = value

            sequence_data.append(sequence_attributes)

        df = pd.DataFrame(sequence_data)

        return df
