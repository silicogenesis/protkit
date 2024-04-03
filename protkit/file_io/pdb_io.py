#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements class `PDBIO` to read and write data
from and to PDB files.  PDB files contain protein structural
information.

The PDB File Format Specification is available at:
https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html

Currently, PDB records related to structural and sequence data is parsed
and saved.  Other fields are ignored. The following fields are processed:

- MODEL (model start)
- ENDMDL (model end)
- ATOM (structural data)
- HETATOM (structural data)
- SEQRES (sequence data)
- TER (chain termination)
- MASTER (record keeping)

Methods are static and can be called without instantiating the class.
The main functions exposed by the class are:

- `load()` to load a protein from a PDB file.
- `save()` to save a protein to a PDB file.

<font color="red">Saving a file in PQR format is not supported yet.</font>

"""

from typing import List, Optional
from protkit.structure.protein import Protein
from protkit.structure.chain import Chain
from protkit.structure.residue import Residue
from protkit.structure.atom import Atom


class PDBIO():
    # --------------------------------------------------------------------------------
    # MODEL
    # --------------------------------------------------------------------------------

    @staticmethod
    def parse_pdb_model_line(line: str):
        """
        Parse a MODEL (start model) line from a PDB file.

        COLUMNS        DATA  TYPE    FIELD          DEFINITION
        ---------------------------------------------------------------------------------------
         1 -  6        Record name   "MODEL "
        11 - 14        Integer       serial         Model serial number.
        """
        pdb_model_dict = {
            "record_name": "MODEL",
            "serial": line[10:14]
        }

        return pdb_model_dict

    @staticmethod
    def create_pdb_model_line(serial: int):
        """
        Creates a MODEL (start model) line for a PDB file.
        """
        return "MODEL".ljust(10) + str(serial).rjust(4) + "".ljust(66)

    # --------------------------------------------------------------------------------
    # ENDMDL
    # --------------------------------------------------------------------------------

    @staticmethod
    def parse_pdb_end_model_line(line: str):
        """
        Parse a ENDMDL (end model) line from a PDB file.

        COLUMNS       DATA  TYPE     FIELD        DEFINITION
        ------------------------------------------------------------------
        1 - 6         Record name   "ENDMDL"

        """
        pdb_end_model_dict = {
            "record_name": "ENDMDL"
        }

        return pdb_end_model_dict

    @staticmethod
    def create_pdb_end_model_line():
        """
        Creates a ENDMDL (end model) line for a PDB file.
        """
        return "ENDMDL".ljust(80)

    # --------------------------------------------------------------------------------
    # MASTER
    # --------------------------------------------------------------------------------

    @staticmethod
    def create_pdb_master_line(num_atoms, num_hetatoms, num_ter, num_remark=0, num_het=0, num_helix=0, num_sheet=0,
                               num_site=0, num_trans=0, num_conect=0, num_seqres=0):
        """
        Creates a MASTER (recordkeeping) line for a PDB file.

        COLUMNS         DATA TYPE     FIELD          DEFINITION
        ----------------------------------------------------------------------------------
         1 -  6         Record name   "MASTER"
        11 - 15         Integer       numRemark      Number of REMARK records
        16 - 20         Integer       "0"
        21 - 25         Integer       numHet         Number of HET records
        26 - 30         Integer       numHelix       Number of HELIX records
        31 - 35         Integer       numSheet       Number of SHEET records
        36 - 40         Integer       numTurn        deprecated
        41 - 45         Integer       numSite        Number of SITE records
        46 - 50         Integer       numXform       Number of coordinate transformation
                                                     records  (ORIGX+SCALE+MTRIX)
        51 - 55         Integer       numCoord       Number of atomic coordinate records
                                                     records (ATOM+HETATM)
        56 - 60         Integer       numTer         Number of TER records
        61 - 65         Integer       numConect      Number of CONECT records
        66 - 70         Integer       numSeq         Number of SEQRES records
        """

        entries = [
            "MASTER".ljust(10),
            str(num_remark).rjust(5),
            "0".rjust(5),
            str(num_het).rjust(5),
            str(num_helix).rjust(5),
            str(num_sheet).rjust(5),
            "0".rjust(5),
            str(num_site).rjust(5),
            str(num_trans).rjust(5),
            str(num_atoms + num_hetatoms).rjust(5),
            str(num_ter).rjust(5),
            str(num_conect).rjust(5),
            str(num_seqres).rjust(5),
            "".ljust(10)
        ]

        return "".join(entries)

    # --------------------------------------------------------------------------------
    # SEQRES
    # --------------------------------------------------------------------------------

    @staticmethod
    def parse_pdb_seqres_line(line: str):
        """
        Parse a SEQRES (sequence) line from a PDB file.
        """
        pdb_seqres_dict = {
            "record_name": "SEQRES",
            "chain_id": line[11],
            "residues": line[19:].strip().split()
        }

        return pdb_seqres_dict

    @staticmethod
    def create_seqres_line(line_no, chain_id, num_residues, residues):
        """
        Creates a SEQRES (sequence) line for a PDB file.

         1 -  6        Record name    "SEQRES"
         8 - 10        Integer        serNum       Serial number of the SEQRES record for  the
                                                   current  chain. Starts at 1 and increments
                                                   by one  each line. Reset to 1 for each chain.
        12             Character      chainID      Chain identifier. This may be any single
                                                   legal  character, including a blank which
                                                   is used if there is only one chain.
        14 - 17        Integer        numRes       Number of residues in the chain.
                                                   This  value is repeated on every record.
        20 - 22        Residue name   resName      Residue name.
        24 - 26        Residue name   resName      Residue name.
        28 - 30        Residue name   resName      Residue name.
        32 - 34        Residue name   resName      Residue name.
        36 - 38        Residue name   resName      Residue name.
        40 - 42        Residue name   resName      Residue name.
        44 - 46        Residue name   resName      Residue name.
        48 - 50        Residue name   resName      Residue name.
        52 - 54        Residue name   resName      Residue name.
        56 - 58        Residue name   resName      Residue name.
        60 - 62        Residue name   resName      Residue name.
        64 - 66        Residue name   resName      Residue name.
        68 - 70        Residue name   resName      Residue name.
        """
        entries = [
            "SEQRES ",
            str(line_no).rjust(3) + " ",
            chain_id + " ",
            str(num_residues).rjust(4) + "  "
        ]
        for residue in residues:
            entries.append(residue["residue"].rjust(3) + " ")
        return "".join(entries).ljust(80)

    # --------------------------------------------------------------------------------
    # TER
    # --------------------------------------------------------------------------------

    @staticmethod
    def parse_pdb_ter_line(line: str):
        """
        Parse a TER (chain termination) line from a PDB file.

        COLUMNS        DATA  TYPE    FIELD           DEFINITION
        -------------------------------------------------------------------------
         1 -  6        Record name   "TER   "
         7 - 11        Integer       serial          Serial number.
        18 - 20        Residue name  resName         Residue name.
        22             Character     chainID         Chain identifier.
        23 - 26        Integer       resSeq          Residue sequence number.
        27             AChar         iCode           Insertion code.

        """
        if line == "TER":
            pdb_ter_dict = {
                "record_name": "TER",
                "serial": "",
                "res_name": "",
                "chain_id": "",
                "res_seq": "",
                "icode": ""
            }
        else:
            line = line.ljust(80)
            pdb_ter_dict = {
                "record_name": "TER",
                "serial": line[6:11],
                "res_name": line[17:21],
                "chain_id": line[21],
                "res_seq": line[22:26],
                "icode": line[26]
            }
        return pdb_ter_dict

    @staticmethod
    def create_pdb_ter_line(serial, res_name, chain_id, res_seq, icode):
        """
        Creates a TER (chain termination) line for a PDB file.
        """
        text = "TER".ljust(6)
        text += str(serial).rjust(5) + "      "
        text += res_name.rjust(3) + " "
        text += chain_id
        text += res_seq.rjust(4)
        text += icode.ljust(1)

        return text.ljust(80)

    # --------------------------------------------------------------------------------
    # ATOM | HETATM
    # --------------------------------------------------------------------------------

    @staticmethod
    def parse_pdb_atom_line(line, is_pqr_format: bool = False):
        """
        Parse a ATOM (individual atom) or HETATM (hetero atom) line from a PDB file.

        COLUMNS        DATA  TYPE    FIELD        DEFINITION
        -------------------------------------------------------------------------------------
         1 -  6        Record name   "ATOM  "
         7 - 11        Integer       serial       Atom  serial number.
        13 - 16        Atom          name         Atom name.
        17             Character     altLoc       Alternate location indicator.
        18 - 20        Residue name  resName      Residue name.
        22             Character     chainID      Chain identifier.
        23 - 26        Integer       resSeq       Residue sequence number.
        27             AChar         iCode        Code for insertion of residues.
        31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
        39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
        47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
        55 - 60        Real(6.2)     occupancy    Occupancy.
        61 - 66        Real(6.2)     tempFactor   Temperature  factor.
        77 - 78        LString(2)    element      Element symbol, right-justified.
        79 - 80        LString(2)    charge       Charge  on the atom.`
        """
        if not is_pqr_format:
            pdb_atom_dict = {
                "record_name": line[0:6].strip().upper(),
                "serial": line[6:11].strip(),
                # Open space
                "atom_name": line[12:16].strip(),
                "alt_loc": line[16].strip(),
                "res_name": line[17:20].strip(),
                # Open space
                "chain_id": line[21].strip(),
                "res_seq": line[22:26].strip(),
                "icode": line[26].strip(),
                # 3x Open space
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "occupancy": float(line[54:60]),
                # "temp_factor": None,
                "temp_factor": float(line[60:66]),
                # 66 - 75?
                "element": line[76:78].strip(),
                "assigned_charge": line[78:80].strip(),
                "calculated_charge": None,
                "radius": None
            }
        else:
            pdb_atom_dict = {
                "record_name": line[0:6].strip().upper(),
                "serial": line[6:11].strip(),
                # Open space
                "atom_name": line[12:16].strip(),
                "alt_loc": line[16].strip(),
                "res_name": line[17:20].strip(),
                # Open space
                "chain_id": line[21].strip(),
                "res_seq": line[22:26].strip(),
                "icode": line[26].strip(),
                # 3x Open space
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),

                "occupancy": 1.0,
                "temp_factor": None,

                "element": line[76:78].strip(),
                "assigned_charge": None,
                "calculated_charge": float(line[54:62]),
                "radius": float(line[62:70])
            }
        return pdb_atom_dict

    @staticmethod
    def create_pdb_atom_line(atom, serial: int, record_name="ATOM"):
        """
        Creates a ATOM line for a PDB file.

        The record_name is expected to be ATOM or HETATM
        """
        # text = "ATOM".ljust(6)
        # text += atom["serial"].rjust(5) + " "
        text = record_name.ljust(6)
        text += str(serial).rjust(5) + " "

        # The PDB spec says that if the element has two characters, e.g. Fe (iron)
        # it should start at the position 13, otherwise it should start at position 14.
        # This seems a little problematic, as some atoms could use 4 characters
        # for example, HG11. The implementation here is consistent with observed PDBs
        # from the RCSB.
        if len(atom["element"]) == 2 or len(atom["atom_name"]) == 4:
            text += atom["atom_name"].ljust(4)
        else:
            text += " " + atom["atom_name"].ljust(3)
        text += ("" if atom["alt_loc"] is None else atom["alt_loc"]).ljust(1)
        text += atom["res_name"].rjust(3) + " "
        text += atom["chain_id"]
        text += atom["res_seq"].rjust(4)
        text += atom["icode"].ljust(1) + "   "
        text += f'{atom["x"]:8.3f}'
        text += f'{atom["y"]:8.3f}'
        text += f'{atom["z"]:8.3f}'
        text += f'{1.0 if atom["occupancy"] is None else atom["occupancy"]:6.2f}'
        text += f'{0.0 if atom["temp_factor"] is None else atom["temp_factor"]:6.2f}'
        text += "".ljust(10)
        text += atom["element"].rjust(2)
        text += ("" if atom["charge"] is None else atom["charge"]).rjust(2)
        return text

    # --------------------------------------------------------------------------------
    # Main parsing and saving functions
    # --------------------------------------------------------------------------------

    @staticmethod
    def parse_pdb_entries_from_string(pdb_string: str, is_pqr_format: bool = False):
        """
        Parse a PDB file in string format into a list of entries.

        Args:
            pdb_string (str): The PDB file as a string.
            is_pqr_format (bool): If the file is in PQR format.

        Returns:
            List[Dict]: A list of parsed entries.
        """

        # Split the string into lines
        lines = pdb_string.split("\n")

        # Parse every line into dictionary information
        entries = []
        for line in lines:
            record_name = line[0:6].strip().upper()
            if record_name == "ATOM":
                pdb_atom_dict = PDBIO.parse_pdb_atom_line(line)
                entries.append(pdb_atom_dict)
            elif record_name == "HETATM":
                pdb_atom_dict = PDBIO.parse_pdb_atom_line(line, is_pqr_format=is_pqr_format)
                entries.append(pdb_atom_dict)
            elif record_name == "TER":
                pdb_ter_dict = PDBIO.parse_pdb_ter_line(line)
                entries.append(pdb_ter_dict)
            elif record_name == "MODEL":
                pdb_model_dict = PDBIO.parse_pdb_model_line(line)
                entries.append(pdb_model_dict)
            elif record_name == "ENDMDL":
                pdb_end_model_dict = PDBIO.parse_pdb_end_model_line(line)
                entries.append(pdb_end_model_dict)
            elif record_name == "SEQRES":
                pdb_seqres_dict = PDBIO.parse_pdb_seqres_line(line)
                entries.append(pdb_seqres_dict)
        return entries


    @staticmethod
    def parse_pdb_entries(file_path: str, is_pqr_format: bool = False):
        """
        Parse PDB file according to the spec provided at:
        https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        """
        # Open the file and extract the individual lines.
        with open(file_path) as file:
            pdb_string = file.read()

        return PDBIO.parse_pdb_entries_from_string(pdb_string, is_pqr_format=is_pqr_format)

    @staticmethod
    def load_from_string(pdb_string: str,
                         is_pqr_format: bool = False,
                         pdb_id: Optional[str] = None) -> List[Protein]:
        protein = Protein(pdb_id=pdb_id)
        models = []
        seqres = {}
        current_model = None
        current_chain: Optional[Chain] = None
        current_chain_id = None
        current_residue: Optional[Residue] = None
        current_residue_id = None

        entries = PDBIO.parse_pdb_entries_from_string(pdb_string, is_pqr_format=is_pqr_format)

        for entry in entries:
            record_name = entry["record_name"]
            if record_name == "SEQRES":
                chain_id = entry["chain_id"]
                if chain_id not in seqres:
                    seqres[chain_id] = []
                seqres[chain_id].extend(entry["residues"])
            elif record_name == "MODEL":
                # Usually, the MODEL field will only be found in structures based on NMR.
                # We only use the first model in a file if provided.
                # There is nothing to do, so we
                pass
            elif record_name == "ENDMDL":
                # Usually, the ENDMDL field will only be found in structures based on NMR.
                # it closes a corresponding MODEL field.
                # Since the end of the model is reached, no more entries should be parsed.
                break
            elif record_name == "TER":
                # The TER field terminates the current chain.
                current_chain = None
                current_chain_id = None
            elif record_name == "ATOM" or record_name == "HETATM":
                # If we switched chains, set the chain.
                if entry["chain_id"] != current_chain_id:
                    current_chain_id = entry["chain_id"]
                    if protein.has_chain(current_chain_id):
                        current_chain = protein.get_chain(current_chain_id)
                    else:
                        current_chain = protein.create_chain(current_chain_id)

                # Check if the residue exists
                residue_id = current_chain_id + entry["res_seq"] + entry["icode"]
                if residue_id != current_residue_id:
                    current_residue_id = residue_id
                    if record_name == "ATOM":
                        # if entry["res_name"] in ResidueTemplate.VALID_DNA_RESIDUES:
                        #     current_residue = current_chain.add_dna_residue(
                        #         residue_type=entry["res_name"],
                        #         sequence_no=int(entry["res_seq"]),
                        #         insertion_code=entry["icode"]
                        #     )
                        # else:
                        current_residue = Residue(
                            residue_type=entry["res_name"],
                            sequence_no=int(entry["res_seq"]),
                            insertion_code=entry["icode"],
                            chain=current_chain)
                        current_chain.add_residue(current_residue)
                    elif record_name == "HETATM":
                        current_residue = Residue(
                            residue_type=entry["res_name"],
                            sequence_no=int(entry["res_seq"]),
                            insertion_code=entry["icode"],
                            chain=current_chain)
                        current_chain.add_residue(current_residue)

                        # current_residue = current_chain.add_het_residue(
                        #     residue_type=entry["res_name"],
                        #     sequence_no=int(entry["res_seq"]),
                        #     insertion_code=entry["icode"]
                        # )

                # Temporary fix need to see if anything breaks (For Evo eg. 2HH1 -> HH12)
                # Fix atom naming
                if entry["atom_name"][0].isnumeric():
                    entry["atom_name"] = entry["atom_name"][1:] + entry["atom_name"][0]
                    if entry["element"] != "H":
                        print(0)
                    if entry["atom_name"][0] == "H":
                        entry["element"] = "H"

                # Add the atom to the residue
                # if entry["alt_loc"] == "":
                #     atom_id = entry["atom_name"]
                # else:
                #     atom_id = entry["atom_name"] + "_" + entry["alt_loc"]
                atom = Atom(
                    element=entry["element"],
                    atom_type=entry["atom_name"],
                    x=entry["x"],
                    y=entry["y"],
                    z=entry["z"],

                    # Serial is dependent on the atom ordering,
                    # so it is not used in the Atom class
                    # pdb_serial=entry["serial"],
                    alt_loc=entry["alt_loc"],
                    occupancy=entry["occupancy"],
                    temp_factor=entry["temp_factor"],
                    assigned_charge=entry["assigned_charge"],
                    # pqr_calculated_charge=entry["calculated_charge"],
                    # pqr_radius=entry["radius"],
                    #
                    is_hetero=(record_name == "HETATM"),

                    residue=current_residue
                )

                # current_residue.add_atom(atom_id, atom)
                current_residue.add_atom(atom.atom_type, atom)

        if seqres:
            for chain_id in seqres:
                if seqres[chain_id]:
                    protein.get_chain(chain_id).set_attribute("seqres", seqres[chain_id])

        # Housekeeping
        # protein.set_residue_indexing()

        return [protein]


    @staticmethod
    def load(file_path: str,
             is_pqr_format: bool = False,
             pdb_id: Optional[str] = None) -> List[Protein]:
        """
        Load a protein from a PDB file.

        Args:
            file_path (str): The path to the PDB file.
            is_pqr_format (bool): If the file is in PQR format.
            pdb_id (str): The PDB ID of the protein.

        Returns:
            List[Protein]: A list of proteins.
        """

        # Open the file and extract the individual lines.
        with open(file_path) as file:
            pdb_string = file.read()

        return PDBIO.load_from_string(pdb_string, is_pqr_format=is_pqr_format, pdb_id=pdb_id)

    @staticmethod
    def save_to_string(protein: Protein,
                       is_pqr_format=False) -> str:
        text_entries = []
        num_seqres = 0
        num_atoms = 0
        num_hetatoms = 0
        num_ter = 0
        serial = 1

        for chain in protein.chains:
            for residue in chain.residues:
                for atom in residue.atoms:
                    if not atom.is_hetero:
                        atom_entry = {
                            "element": atom.element,
                            "atom_name": atom.atom_type,
                            "alt_loc": atom.get_attribute("pdb_alt_loc"),
                            "res_name": residue.residue_type,
                            "chain_id": chain.chain_id,
                            "res_seq": str(residue.sequence_no),
                            "icode": residue.insertion_code,
                            "x": atom.x,
                            "y": atom.y,
                            "z": atom.z,
                            "occupancy": atom.get_attribute("occupancy"),
                            # "temp_factor": atom._pdb_temp_factor if atom._pdb_temp_factor is not None else 0.0,
                            "temp_factor": atom.get_attribute("temp_factor"),
                            "charge": atom.get_attribute("assigned_charge")
                        }
                        text_entries.append(PDBIO.create_pdb_atom_line(atom_entry, serial))
                        serial += 1
                        num_atoms += 1

        # Bookkeeping records
        text_entries.append(PDBIO.create_pdb_master_line(num_atoms, num_hetatoms, num_ter, num_seqres=num_seqres))
        text_entries.append("END".ljust(80))

        return "\n".join(text_entries)

    @staticmethod
    def save(protein: Protein,
             file_path: str,
             is_pqr_format=False) -> None:

        pdb_string = PDBIO.save_to_string(protein, is_pqr_format=is_pqr_format)

        with open(file_path, "w") as file:
            file.write(pdb_string)

    # TODO: from_pdb_text and to_pdb_text was added by Claudio
    # The intent is that one should parse proteins from text
    # without working with files. However, a lot of code copying
    # was made and it need to be reworked significantly.

    # @staticmethod
    # def from_pdb_text(lines) -> List[Protein]:
    #     lines = lines.split("\n")
    #
    #     entries = []
    #     for line in lines:
    #         record_name = line[0:6].strip().upper()
    #         if record_name == "ATOM":
    #             pdb_atom_dict = PDBIO.parse_pdb_atom_line(line)
    #             entries.append(pdb_atom_dict)
    #         elif record_name == "HETATM":
    #             pdb_atom_dict = PDBIO.parse_pdb_atom_line(line)
    #             entries.append(pdb_atom_dict)
    #         elif record_name == "TER":
    #             pdb_ter_dict = PDBIO.parse_pdb_ter_line(line)
    #             entries.append(pdb_ter_dict)
    #         elif record_name == "MODEL":
    #             pdb_model_dict = PDBIO.parse_pdb_model_line(line)
    #             entries.append(pdb_model_dict)
    #         elif record_name == "ENDMDL":
    #             pdb_end_model_dict = PDBIO.parse_pdb_end_model_line(line)
    #             entries.append(pdb_end_model_dict)
    #         elif record_name == "SEQRES":
    #             pdb_seqres_dict = PDBIO.parse_pdb_seqres_line(line)
    #             entries.append(pdb_seqres_dict)
    #
    #     protein = Protein()
    #     models = []
    #     seqres = {}
    #     current_model = None
    #     current_chain: Optional[Chain] = None
    #     current_chain_id = None
    #     current_residue: Optional[Residue] = None
    #     current_residue_id = None
    #
    #     for entry in entries:
    #         record_name = entry["record_name"]
    #         if record_name == "SEQRES":
    #             chain_id = entry["chain_id"]
    #             if chain_id not in seqres:
    #                 seqres[chain_id] = []
    #             seqres[chain_id].extend(entry["residues"])
    #         elif record_name == "MODEL":
    #             # Usually, the MODEL field will only be found in structures based on NMR.
    #             # We only use the first model in a file if provided.
    #             # There is nothing to do, so we
    #             pass
    #         elif record_name == "ENDMDL":
    #             # Usually, the ENDMDL field will only be found in structures based on NMR.
    #             # it closes a corresponding MODEL field.
    #             # Since the end of the model is reached, no more entries should be parsed.
    #             break
    #         elif record_name == "TER":
    #             # The TER field terminates the current chain.
    #             current_chain = None
    #             current_chain_id = None
    #         elif record_name == "ATOM" or record_name == "HETATM":
    #             # If we switched chains, set the chain.
    #             if entry["chain_id"] != current_chain_id:
    #                 current_chain_id = entry["chain_id"]
    #                 if protein.has_chain(current_chain_id):
    #                     current_chain = protein.get_chain(current_chain_id)
    #                 else:
    #                     current_chain = protein.create_chain(current_chain_id)
    #
    #             # Check if the residue exists
    #             residue_id = current_chain_id + entry["res_seq"] + entry["icode"]
    #             if residue_id != current_residue_id:
    #                 current_residue_id = residue_id
    #                 if record_name == "ATOM":
    #                     # if entry["res_name"] in ResidueTemplate.VALID_DNA_RESIDUES:
    #                     #     current_residue = current_chain.add_dna_residue(
    #                     #         residue_type=entry["res_name"],
    #                     #         sequence_no=int(entry["res_seq"]),
    #                     #         insertion_code=entry["icode"]
    #                     #     )
    #                     # else:
    #                     current_residue = Residue(
    #                         residue_type=entry["res_name"],
    #                         sequence_no=int(entry["res_seq"]),
    #                         insertion_code=entry["icode"],
    #                         chain=current_chain)
    #                     current_chain.add_residue(current_residue)
    #                 elif record_name == "HETATM":
    #                     current_residue = Residue(
    #                         residue_type=entry["res_name"],
    #                         sequence_no=int(entry["res_seq"]),
    #                         insertion_code=entry["icode"],
    #                         chain=current_chain)
    #                     current_chain.add_residue(current_residue)
    #
    #                     # current_residue = current_chain.add_het_residue(
    #                     #     residue_type=entry["res_name"],
    #                     #     sequence_no=int(entry["res_seq"]),
    #                     #     insertion_code=entry["icode"]
    #                     # )
    #
    #             # Temporary fix need to see if anything breaks (For Evo eg. 2HH1 -> HH12)
    #             # Fix atom naming
    #             if entry["atom_name"][0].isnumeric():
    #                 entry["atom_name"] = entry["atom_name"][1:] + entry["atom_name"][0]
    #                 if entry["element"] != "H":
    #                     print(0)
    #                 if entry["atom_name"][0] == "H":
    #                     entry["element"] = "H"
    #
    #             # Add the atom to the residue
    #             if entry["alt_loc"] == "":
    #                 atom_id = entry["atom_name"]
    #             else:
    #                 atom_id = entry["atom_name"] + "_" + entry["alt_loc"]
    #             atom = Atom(
    #                 element=entry["element"],
    #                 atom_type=entry["atom_name"],
    #                 x=entry["x"],
    #                 y=entry["y"],
    #                 z=entry["z"],
    #
    #                 # pdb_serial=entry["serial"],
    #                 # pdb_alt_loc=entry["alt_loc"],
    #                 # pdb_occupancy=entry["occupancy"],
    #                 # pdb_temp_factor=entry["temp_factor"],
    #                 # pdb_assigned_charge=entry["assigned_charge"],
    #                 # pqr_calculated_charge=entry["calculated_charge"],
    #                 # pqr_radius=entry["radius"],
    #                 #
    #                 is_hetero=(record_name == "HETATM"),
    #
    #                 residue=current_residue
    #             )
    #             current_residue.add_atom(atom_id, atom)
    #
    #     # if seqres:
    #     #     for chain_id in seqres:
    #     #         if seqres[chain_id]:
    #     #             protein.chains[chain_id].sequence = seqres[chain_id]
    #
    #     # Housekeeping
    #     # protein.set_residue_indexing()
    #
    #     return [protein]

    # @staticmethod
    # def to_pdb_text(protein: Protein):
    #     text_entries = ""
    #     num_seqres = 0
    #     num_atoms = 0
    #     num_hetatoms = 0
    #     num_ter = 0
    #     serial = 1
    #
    #     for chain in protein.chains:
    #         for residue in chain.residues:
    #             for atom in residue.atoms:
    #                 if not atom.is_hetero:
    #                     atom_entry = {
    #                         "element": atom.element,
    #                         "atom_name": atom.atom_type,
    #                         "alt_loc": atom.get_attribute("pdb_alt_loc"),
    #                         "res_name": residue.residue_type,
    #                         "chain_id": chain.chain_id,
    #                         "res_seq": str(residue.sequence_no),
    #                         "icode": residue.insertion_code,
    #                         "x": atom.x,
    #                         "y": atom.y,
    #                         "z": atom.z,
    #                         "occupancy": atom.get_attribute("occupancy"),
    #                         # "temp_factor": atom._pdb_temp_factor if atom._pdb_temp_factor is not None else 0.0,
    #                         "temp_factor": atom.get_attribute("temp_factor"),
    #                         "charge": atom.get_attribute("assigned_charge")
    #                     }
    #                     text_entries += PDBIO.create_pdb_atom_line(atom_entry, serial)
    #                     text_entries += "\n"
    #                     serial += 1
    #                     num_atoms += 1
    #
    #     # Bookkeeping records
    #     text_entries += PDBIO.create_pdb_master_line(num_atoms, num_hetatoms, num_ter, num_seqres=num_seqres)
    #     text_entries += "END".ljust(80)
    #     text_entries += "\n"
    #
    #     return text_entries

    # @staticmethod
    # def save_to_pdb(file_path: str, protein):
    #     entries = []
    #
    #     num_seqres = 0
    #     num_atoms = 0
    #     num_hetatoms = 0
    #     num_ter = 0
    #     serial = 1
    #
    #     for chain in protein.chains:
    #         for residue in chain.residues:
    #             # ATOM entries
    #             for atom in residue.atoms:
    #                 if not atom.is_hetero:
    #                     atom_entry = {
    #                         "element": atom.element,
    #                         "atom_name": atom.atom_type,
    #                         "alt_loc": atom._pdb_alt_loc,
    #                         "res_name": residue.residue_type,
    #                         "chain_id": chain.chain_id,
    #                         "res_seq": str(residue._pdb_sequence_no),
    #                         "icode": residue._pdb_insertion_code,
    #                         "x": atom.x,
    #                         "y": atom.y,
    #                         "z": atom.z,
    #                         "occupancy": atom._pdb_occupancy,
    #                         "temp_factor": atom._pdb_temp_factor if atom._pdb_temp_factor is not None else 0.0,
    #                         "charge": atom._pdb_assigned_charge
    #                     }
    #                     entries.append(PDBParser.create_pdb_atom_line(atom_entry, serial))
    #                     serial += 1
    #                     num_atoms += 1
    #
    #         # TER entry
    #         if chain.num_residues - chain.num_hetero_residues > 0:
    #             entries.append(
    #                 PDBParser.create_pdb_ter_line(serial, atom_entry["res_name"], chain.chain_id, atom_entry["res_seq"],
    #                                               atom_entry["icode"]))
    #             serial = serial + 1
    #             num_ter += 1
    #
    #     for chain in protein.chains:
    #         for residue in chain.residues:
    #             # HETATM entries
    #             for atom in residue.atoms:
    #                 if atom.is_hetero:
    #                     atom_entry = {
    #                         "element": atom.element,
    #                         "atom_name": atom.atom_type,
    #                         "alt_loc": atom._pdb_alt_loc,
    #                         "res_name": residue.residue_type,
    #                         "chain_id": chain.chain_id,
    #                         "res_seq": str(residue._pdb_sequence_no),
    #                         "icode": residue._pdb_insertion_code,
    #                         "x": atom.x,
    #                         "y": atom.y,
    #                         "z": atom.z,
    #                         "occupancy": atom._pdb_occupancy,
    #                         "temp_factor": atom._pdb_temp_factor if atom._pdb_temp_factor is not None else 0.0,
    #                         "charge": atom._pdb_assigned_charge
    #                     }
    #                     entries.append(PDBParser.create_pdb_atom_line(atom_entry, serial, "HETATM"))
    #                     serial += 1
    #                     num_hetatoms += 1
    #
    #         # Hetero residues have no TER records, so no need to write a termination code.
    #
    #     # Bookkeeping records
    #     entries.append(PDBParser.create_pdb_master_line(num_atoms, num_hetatoms, num_ter, num_seqres=num_seqres))
    #     entries.append("END".ljust(80))
    #
    #     # Write file
    #     with open(file_path, "w") as file:
    #         file.write("\n".join(entries))

    # @staticmethod
    # def load_from_pdb(file_path: str):
    #     entries = PDBParser.parse_pdb_entries(file_path)
    #
    #     protein = Protein()
    #     models = []
    #     seqres = {}
    #     current_model = None
    #     current_chain: Chain = None
    #     current_chain_id = None
    #     current_residue: Residue = None
    #     current_residue_id = None
    #
    #     for entry in entries:
    #         record_name = entry["record_name"]
    #         if record_name == "SEQRES":
    #             chain_id = entry["chain_id"]
    #             if chain_id not in seqres:
    #                 seqres[chain_id] = []
    #             seqres[chain_id].extend(entry["residues"])
    #         elif record_name == "MODEL":
    #             # Usually, the MODEL field will only be found in structures based on NMR.
    #             # We only use the first model in a file if provided.
    #             # There is nothing to do, so we
    #             pass
    #         elif record_name == "ENDMDL":
    #             # Usually, the ENDMDL field will only be found in structures based on NMR.
    #             # it closes a corresponding MODEL field.
    #             # Since the end of the model is reached, no more entries should be parsed.
    #             break
    #         elif record_name == "TER":
    #             # The TER field terminates the current chain.
    #             current_chain = None
    #             current_chain_id = None
    #         elif record_name == "ATOM" or record_name == "HETATM":
    #             # If we switched chains, set the chain.
    #             if entry["chain_id"] != current_chain_id:
    #                 current_chain_id = entry["chain_id"]
    #                 if current_chain_id in protein._chains:
    #                     current_chain = protein._chains[current_chain_id]
    #                 else:
    #                     current_chain = protein.create_chain(current_chain_id)
    #
    #             # Check if the residue exists
    #             residue_id = current_chain_id + entry["res_seq"] + entry["icode"]
    #             if residue_id != current_residue_id:
    #                 current_residue_id = residue_id
    #                 if record_name == "ATOM":
    #                     # if entry["res_name"] in ResidueTemplate.VALID_DNA_RESIDUES:
    #                     #     current_residue = current_chain.add_dna_residue(
    #                     #         residue_type=entry["res_name"],
    #                     #         sequence_no=int(entry["res_seq"]),
    #                     #         insertion_code=entry["icode"]
    #                     #     )
    #                     # else:
    #                     current_residue = Residue(
    #                         residue_type=entry["res_name"],
    #                         pdb_sequence_no=int(entry["res_seq"]),
    #                         pdb_insertion_code=entry["icode"],
    #                         chain=current_chain)
    #                     current_chain.add_residue(current_residue)
    #                 elif record_name == "HETATM":
    #                     current_residue = Residue(
    #                         residue_type=entry["res_name"],
    #                         pdb_sequence_no=int(entry["res_seq"]),
    #                         pdb_insertion_code=entry["icode"],
    #                         chain=current_chain)
    #                     current_chain.add_residue(current_residue)
    #
    #                     # current_residue = current_chain.add_het_residue(
    #                     #     residue_type=entry["res_name"],
    #                     #     sequence_no=int(entry["res_seq"]),
    #                     #     insertion_code=entry["icode"]
    #                     # )
    #
    #             # Temporary fix need to see if anything breaks (For Evo eg. 2HH1 -> HH12)
    #             # Fix atom naming
    #             if entry["atom_name"][0].isnumeric():
    #                 entry["atom_name"] = entry["atom_name"][1:] + entry["atom_name"][0]
    #                 if entry["element"] != "H":
    #                     print(0)
    #                 if entry["atom_name"][0] == "H":
    #                     entry["element"] = "H"
    #
    #             # Add the atom to the residue
    #             if entry["alt_loc"] == "":
    #                 atom_id = entry["atom_name"]
    #             else:
    #                 atom_id = entry["atom_name"] + "_" + entry["alt_loc"]
    #             atom = Atom(
    #                 element=entry["element"],
    #                 atom_type=entry["atom_name"],
    #                 x=entry["x"],
    #                 y=entry["y"],
    #                 z=entry["z"],
    #
    #                 pdb_serial=entry["serial"],
    #                 pdb_alt_loc=entry["alt_loc"],
    #                 pdb_occupancy=entry["occupancy"],
    #                 pdb_temp_factor=entry["temp_factor"],
    #                 pdb_assigned_charge=entry["assigned_charge"],
    #                 pqr_calculated_charge=entry["calculated_charge"],
    #                 pqr_radius=entry["radius"],
    #
    #                 is_hetero=(record_name == "HETATM"),
    #
    #                 residue=current_residue
    #             )
    #             current_residue.add_atom(atom_id, atom)
    #
    #     # if seqres:
    #     #     for chain_id in seqres:
    #     #         if seqres[chain_id]:
    #     #             protein.chains[chain_id].sequence = seqres[chain_id]
    #
    #     # Housekeeping
    #     protein.set_residue_indexing()
    #
    #     return protein
