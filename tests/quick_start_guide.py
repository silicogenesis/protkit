def quick_start_example():
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


def prep_data():
    from protkit.download import Download
    from protkit.file_io import ProtIO

    # Download a PDB file from the RCSB PDB database and save it to a file.
    Download.download_pdb_file_from_rcsb("1ahw", "1ahw.pdb")
    Download.download_pdb_file_from_rcsb("3i40", "3i40.pdb")
    Download.download_pdb_file_from_rcsb("4nkq", "4nkq.pdb")
    Download.download_pdb_file_from_rcsb("6bom", "6bom.pdb")

    # Load a PDB file into a Protein object.
    ProtIO.convert("1ahw.pdb", "1ahw.prot")
    ProtIO.convert("3i40.pdb", "3i40.prot")
    ProtIO.convert("4nkq.pdb", "4nkq.prot")
    ProtIO.convert("6bom.pdb", "6bom.prot")

    # Download a FASTA file from RCSB PDB and save it to file.
    Download.download_fasta_file_from_rcsb("1ahw", "1ahw.fasta")
    Download.download_fasta_file_from_rcsb("3i40", "1ahw.3i40")
    Download.download_fasta_file_from_rcsb("4nkq", "4nkq.fasta")
    Download.download_fasta_file_from_rcsb("6bom", "6bom.fasta")


def download_pdb_example():
    from protkit.download import Download

    # Download a PDB file from the RCSB PDB database and save it to a file.
    Download.download_pdb_file_from_rcsb("1ahw", "data/pdb/rcsb/1ahw.pdb")

    # Download a PDB file from the SAbDab database and save it to a file.
    Download.download_pdb_file_from_sabdab("1ahw", "data/pdb/sabdab/1ahw.pdb")

    # Download multiple PDB files from the RCSB PDB in parallel and save them to a directory.
    Download.download_pdb_files_from_rcsb(["1ahw", "1a4y", "1a6m"], "data/pdb/rcsb/", n_jobs=3)


def download_fasta_example():
    from protkit.download import Download

    # Download a Fasta file from the Uniprot database and save it to a file.
    Download.download_fasta_file_from_uniprot("P01308", "data/fasta/uniprot/P01308.fasta")

    # Download multiple Fasta files from the Uniprot database in parallel and save them to a directory.
    Download.download_fasta_files_from_uniprot(["P01308", "P68871", "P00533"], "data/fasta/uniprot", n_jobs=3)

    # Download a Fasta file from the RCSB PDB database and save it to a file.
    Download.download_fasta_file_from_rcsb("1ahw", "data/fasta/rcsb/1ahw.fasta")

    # Download multiple Fasta files from the RCSB PDB database in parallel and save them to a directory.
    Download.download_fasta_files_from_rcsb(["1ahw", "1a4y", "1a6m"], "data/fasta/rcsb/", n_jobs=3)


def download_cif_example():
    from protkit.download import Download

    # Download a CIF file from the RCSB PDB database and save it to a file.
    Download.download_cif_file_from_rcsb("1ahw", "data/cif/rcsb/1ahw.cif")

    # Download multiple CIF files from the RCSB PDB database in parallel and save them to a directory.
    Download.download_cif_files_from_rcsb(["1ahw", "1a4y", "1a6m"], "data/cif/rcsb/")

    # Download a binary CIF file from the RCSB PDB database and save it to a file.
    Download.download_binary_cif_file_from_rcsb("1ahw", "data/cif/rcsb/1ahw.bcif")

    # Download multiple binary CIF files from the RCSB PDB database in parallel and save them to a directory.
    Download.download_binary_cif_files_from_rcsb(["1ahw", "1a4y", "1a6m"], "data/cif/rcsb/")


def file_io_pdb():
    from protkit.file_io import PDBIO

    protein = PDBIO.load("1ahw.pdb")[0]

    protein.keep_chains(["A", "B"])

    PDBIO.save(protein, "1ahw_AB.pdb")


def file_io_prot():
    from protkit.file_io import ProtIO

    protein = ProtIO.convert("1ahw.pdb", "1ahw.prot")

    protein = ProtIO.load("1ahw.prot")[0]

    ProtIO.save(protein, "1ahw.prot")


def file_io_prot_multisave():
    from protkit.file_io import ProtIO, PDBIO

    protein1 = PDBIO.load("1ahw.pdb")[0]
    protein2 = PDBIO.load("4nkq.pdb")[0]

    ProtIO.save([protein1, protein2], "database.prot")


def file_io_prot_compression():
    from protkit.file_io import ProtIO

    protein = ProtIO.load("1ahw.prot", decompress=True)[0]
    ProtIO.save(protein, "1ahw.prot.json", compress=False)


def file_io_fasta_load():
    from protkit.file_io import FastaIO

    sequences = FastaIO.load("1ahw.fasta")
    for sequence in sequences:
        print(f"{sequence.chain_id}: {sequence}")


def file_io_fasta_save():
    from protkit.file_io import FastaIO

    sequences = FastaIO.load("1ahw.fasta")
    FastaIO.save(sequences, "1ahw_copy.fasta")
    FastaIO.save(sequences[0], "1ahw_A.fasta")


def file_io_pqr():
    from protkit.file_io import PQRIO

    protein = PQRIO.load("1ahw.pqr")[0]

    PQRIO.save(protein, "1ahw_copy.pqr")


def file_io_cif():
    from protkit.file_io import MMTFIO

    protein = MMTFIO.load("1ahw.cif")[0]

    MMTFIO.save(protein, "1ahw_copy.cif")


def file_io_mmtf():
    from protkit.file_io import MMTFIO

    protein = MMTFIO.load("1ahw.mmtf")[0]


def representation_protein():
    # Loading a Protein
    from protkit.file_io import ProtIO

    protein = ProtIO.load("1ahw.prot")[0]

    # Accessing Chains in a Protein
    for chain in protein.chains:
        print(f"{chain.chain_id}: {chain.num_residues} residues")
        print(chain.sequence)

    chain_a = protein.get_chain("A")
    print(chain_a.sequence)

    protein2 = protein.copy()
    protein2.keep_chains(["A", "B"])
    print(protein.num_chains)
    print(protein2.num_chains)

    protein3 = protein.copy()
    protein3.remove_chains(["C", "D"])
    print(protein3.num_chains)

    # Isolating Chains
    protein2 = protein.copy()
    protein2.keep_chains(["A", "B"])
    print(protein.num_chains)
    print(protein2.num_chains)

    # Renaming Chains
    protein.rename_chain("A", "Z")


def representation_protein_access():
    from protkit.file_io import ProtIO

    protein = ProtIO.load("1ahw.prot")[0]

    # Accessing Residues in a Protein
    for residue in protein.residues:
        print(f"{residue.id}: {residue.residue_type}")

    # Accessing Atoms in a Protein
    coordinates = [(atom.x, atom.y, atom.z) for atom in protein.atoms]
    print(coordinates)

    # Powerful filters
    atoms = protein.filter_atoms(
        chain_criteria=[("chain_id", "A")],
        residue_criteria=[("residue_type", "PRO")],
        atom_criteria=[("atom_type", ["C", "CA", "N"])])
    for atom in atoms:
        print(f"{atom.residue.residue_code}, {atom.atom_type}, {atom.x}, {atom.y}, {atom.z}")

    # Protein statistics
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


def representation_chain():
    from protkit.file_io import ProtIO

    # Identifying a Chain
    protein = ProtIO.load("1ahw.prot")[0]
    chain = protein.get_chain("A")
    print(chain.chain_id)

    protein.rename_chain("A", "Z")
    print(chain.chain_id)
    chain.chain_id = "X"

    # Chain Sequence
    print(chain.sequence)

    # Accessing Residues in a Chain
    for residue in chain.residues:
        print(f"{residue.id}: {residue.residue_type}")

    residue = chain.get_residue(0)
    print(residue.id)

    # Accessing Atoms in a Chain
    for atom in chain.atoms:
        print(f"{atom.id}")


def representation_residue():
    # Identifying a Residue
    from protkit.file_io import ProtIO

    protein = ProtIO.load("1ahw.prot")[0]
    chain = protein.get_chain("A")
    residue = chain.get_residue(0)

    # Core Residue Properties
    print(residue.residue_type)  # eg. GLY
    print(residue.short_code)  # eg. G
    print(residue.sequence_no)  # eg. 100
    print(residue.insertion_code)  # eg. A
    print(residue.is_disordered)  # eg. True
    print(residue.is_hetero)  # eg. False

    print(residue.sequence_code)  # eg. 100A
    print(residue.residue_code)  # eg. GLY100A
    print(residue.id)  # eg. A:GLY100A

    # Accessing Atoms in a Residue
    for atom in residue.atoms:
        print(f"{atom.id}")

    # Modifying the Atoms in a Residue
    residue.keep_backbone_atoms()

    residue.keep_atoms(["N", "CA", "C"])

    residue.remove_atoms(["OXT"])

    residue.remove_hydrogen_atoms()


def representation_atom():
    # Identifying an Atom
    from protkit.file_io import ProtIO

    protein = ProtIO.load("1ahw.prot")[0]
    residue = protein.get_chain("A").get_residue(0)
    atom = residue.get_atom("CA")

    print(atom.id)
    print(atom.atom_type)

    print(atom.element)  # eg. C
    print(atom.x, atom.y, atom.z)  # eg. 12.446, 8.496, 43.176
    print(atom.is_disordered)  # eg. False
    print(atom.is_hetero)  # eg. False

def representation_atom_fix_disordered():
    from protkit.file_io import ProtIO

    protein = ProtIO.load("1ahw.prot")[0]
    protein.fix_disordered_atoms()


def representation_attributes():
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


def quality_create_copy():
    from protkit.file_io import ProtIO

    protein = ProtIO.load("3i40.prot")[0]
    protein2 = protein.copy()


def quality_remove_water_molecules():
    from protkit.file_io import ProtIO

    protein = ProtIO.load("6bom.prot")[0]
    print(f"{protein.num_water_residues}")
    protein.remove_water_residues()
    print(f"{protein.num_water_residues}")

    protein.get_chain("A").remove_water_residues()


def quality_remove_hetero_residues():
    from protkit.file_io import ProtIO

    protein = ProtIO.load("6bom.prot")[0]
    print(f"{protein.num_hetero_residues}")
    protein.remove_hetero_residues("TRS")
    print(f"{protein.num_hetero_residues}")

    protein.remove_hetero_residues(["TRS", "SO4"])


def quality_fixing_disordered_atoms():
    from protkit.file_io import ProtIO

    protein = ProtIO.load("3i40.prot")[0]
    print(f"{protein.num_disordered_atoms}")
    print(protein.get_chain("A").get_residue(13).is_disordered)
    print(protein.get_chain("A").get_residue(13).get_atom("CA").is_disordered)
    print(protein.get_chain("A").get_residue(13).get_atom("CB").is_disordered)

    protein.fix_disordered_atoms()
    print(protein.get_chain("A").get_residue(13).get_atom("CB").is_disordered)
    print(f"{protein.num_disordered_atoms}")


def quality_removing_hydrogen_atoms():
    from protkit.file_io import ProtIO

    protein = ProtIO.load("1a4y_A_B.prot")[0]

    print(f"{protein.num_atoms} atoms")
    print(f"{protein.num_heavy_atoms} heavy atoms")
    print(f"{protein.num_hydrogen_atoms} hydrogen atoms")

    protein.remove_hydrogen_atoms()

    print(f"{protein.num_hydrogen_atoms} hydrogen atoms after removal")


def properties_hydrophobicity():
    from protkit.file_io import ProtIO
    from protkit.properties import Hydrophobicity

    protein = ProtIO.load("1ahw.prot")[0]
    Hydrophobicity.hydrophobicity_of_protein(protein, assign_attribute=True)

    print(protein.get_attribute("hydrophobicity"))
    print(protein.get_chain("A").get_attribute("hydrophobicity"))
    print(protein.get_chain("B").get_attribute("hydrophobicity"))
    print(protein.get_chain("A").get_residue(0).get_attribute("hydrophobicity"))


def properties_hydrophobicity_class():
    from protkit.file_io import ProtIO
    from protkit.properties import Hydrophobicity

    protein = ProtIO.load("1ahw.prot")[0]
    Hydrophobicity.hydrophobicity_classes_of_protein(protein, assign_attribute=True)
    print(Hydrophobicity.HYDROPHOBICITY_CLASS_STRING[protein.get_chain("A").get_residue(0).get_attribute("hydrophobicity_class")])


def properties_mass():
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


def properties_chemical_class():
    from protkit.file_io import ProtIO
    from protkit.properties import ChemicalClass

    protein = ProtIO.load("1ahw.prot")[0]
    ChemicalClass.chemical_classes_of_protein(protein, assign_attribute=True)
    print(ChemicalClass.CHEMICAL_CLASS_STRING[protein.get_chain("A").get_residue(0).get_attribute("chemical_class")])


def properties_charge():
    from protkit.file_io import ProtIO
    from protkit.properties import Charge

    protein = ProtIO.load("1ahw.prot")[0]
    Charge.charge_of_protein(protein, assign_attribute=True)
    chain = protein.get_chain("A")
    residue = chain.get_residue(0)
    print(protein.get_attribute("charge"))
    print(chain.get_attribute("charge"))
    print(residue.get_attribute("charge"))


def properties_polarity():
    from protkit.file_io import ProtIO
    from protkit.properties import Polarity

    protein = ProtIO.load("1ahw.prot")[0]
    Polarity.polarities_of_protein(protein, assign_attribute=True)
    residue = protein.get_chain("A").get_residue(0)
    print(Polarity.POLARITY_STRING[residue.get_attribute("polarity")])


def properties_donors_acceptors():
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


def properties_surface_area():
    from protkit.file_io import ProtIO
    from protkit.properties import SurfaceArea

    protein = ProtIO.load("1ahw.prot")[0]
    chain = protein.get_chain("A")
    residue = chain.get_residue(0)

    SurfaceArea.surface_area_of_protein(protein, assign_attribute=True)
    print(protein.get_attribute("surface_area"))
    print(chain.get_attribute("surface_area"))
    print(residue.get_attribute("surface_area"))


def properties_volume():
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


def properties_volume_class():
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


def properties_bounds_and_center():
    from protkit.file_io import ProtIO
    from protkit.properties import Bounds

    protein = ProtIO.load("1ahw.prot")[0]
    chain = protein.get_chain("A")
    residue = chain.get_residue(0)

    Bounds.bounds_of_protein(protein, assign_attribute=True)
    Bounds.center_of_protein(protein, assign_attribute=True)

    print(protein.get_attribute("bounds"))
    print(chain.get_attribute("bounds"))
    print(residue.get_attribute("bounds"))

    print(protein.get_attribute("center"))
    print(chain.get_attribute("center"))
    print(residue.get_attribute("center"))


def properties_bond_lengths():
    from protkit.file_io import ProtIO
    from protkit.properties import BondLengths

    protein = ProtIO.load("1ahw.prot")[0]
    residue = protein.get_chain("A").get_residue(0)

    BondLengths.bond_lengths_of_protein(protein, assign_attribute=True)
    for (atom1, atom2), length in residue.get_attribute("bond_lengths").items():
        print(f"{atom1}-{atom2}: {length}")


def properties_peptide_lengths():
    from protkit.file_io import ProtIO
    from protkit.properties import BondLengths

    protein = ProtIO.load("1ahw.prot")[0]
    chain = protein.get_chain("A")

    BondLengths.peptide_bond_lengths_of_protein(protein, assign_attribute=True)
    print(chain.get_attribute("peptide_bond_lengths"))


def properties_bond_angles():
    from protkit.file_io import ProtIO
    from protkit.properties import BondAngles

    protein = ProtIO.load("1ahw.prot")[0]
    residue = protein.get_chain("A").get_residue(0)

    BondAngles.bond_angles_of_protein(protein, assign_attribute=True)
    for (atom1, atom2, atom3), angle in residue.get_attribute("bond_angles").items():
        print(f"{atom1}-{atom2}-{atom3}: {angle}")


def properties_dihedral_angles():
    from protkit.file_io import ProtIO
    from protkit.properties import DihedralAngles

    protein = ProtIO.load("1ahw.prot")[0]
    residue = protein.get_chain("A").get_residue(1)

    DihedralAngles.dihedral_angles_of_protein(protein, assign_attribute=True)
    for angle_name, angle in residue.get_attribute("dihedral_angles").items():
        print(f"{angle_name}: {angle}")


def properties_circular_variance():
    from protkit.file_io import ProtIO
    from protkit.properties import CircularVariance

    protein = ProtIO.load("1ahw.prot")[0]
    residue = protein.get_chain("A").get_residue(1)
    atom = residue.get_atom("CA")

    CircularVariance.circular_variance_by_residue(protein, assign_attribute=True)
    print(residue.get_attribute("cv_residue"))

    CircularVariance.circular_variance_by_atom(protein, assign_attribute=True)
    print(atom.get_attribute("cv_atom"))


def properties_interface_atoms():
    from protkit.file_io import ProtIO
    from protkit.properties import Interface

    protein = ProtIO.load("1ahw.prot")[0]
    atoms1 = list(protein.filter_atoms(chain_criteria=[("chain_id", ["A", "B"])]))
    atoms2 = list(protein.filter_atoms(chain_criteria=[("chain_id", ["C"])]))

    Interface.interface_atoms(atoms1, atoms2, cutoff=5.0, assign_attribute=True)
    for atom in atoms1:
        if atom.get_attribute("in_interface"):
            print(atom.id)


def properties_interface_residues():
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

def properties_interface_atom_pairs():
    from protkit.file_io import PDBIO
    from protkit.properties import Interface

    protein = PDBIO.load("1ahw.pdb")[0]
    atoms1 = list(protein.filter_atoms(chain_criteria=[("chain_id", ["A", "B"])]))
    atoms2 = list(protein.filter_atoms(chain_criteria=[("chain_id", ["C"])]))

    pairs = Interface.interface_atom_pairs(atoms1, atoms2, cutoff=5.0, assign_attribute=True)

    # pairs is now a list of (atom_from_A, atom_from_B) that interact.
    for aA, aB in pairs:
        print(f"{aA.id} interacts with {aB.id}")

def properties_interface_residue_pairs():
    from protkit.file_io import PDBIO
    from protkit.properties import Interface

    protein = PDBIO.load("1ahw.pdb")[0]
    residues1 = list(protein.filter_residues(chain_criteria=[("chain_id", ["A", "B"])]))
    residues2 = list(protein.filter_residues(chain_criteria=[("chain_id", ["C"])]))

    pairs = Interface.interface_residue_pairs(residues1, residues2, cutoff=5.0, assign_attribute=True)

    for rA, rB in pairs:
        print(f"{rA.id} interacts with {rB.id}")

def tools_reduce():
    from protkit.tools.reduce_adaptor import ReduceAdaptor
    from protkit.file_io import ProtIO

    reduce_bin_path = "/usr/local/bin/reduce"
    reduce = ReduceAdaptor(reduce_bin_path)

    protein = ProtIO.load("1ahw.prot")[0]
    protein_protonated = reduce.protonate(protein)
    protein_deprotonated = reduce.deprotonate(protein)

def tools_freesasa():
    from protkit.tools.freesasa_adaptor import FreeSASAAdaptor
    from protkit.file_io import ProtIO

    freesasa = FreeSASAAdaptor()
    protein = ProtIO.load("1ahw.prot")[0]

    atoms = list(protein.atoms)
    atom_areas = freesasa.calculate_surface_area(atoms)
    protein_area = sum(atom_areas)
    print(protein_area)

def ml_dataframe():
    from protkit.file_io import ProtIO
    from protkit.ml import ProteinToDataframe

    featurizer = ProteinToDataframe()
    protein = ProtIO.load("1ahw.prot")[0]
    protein.pdb_id = "1ahw"

    df_chains = featurizer.chains_dataframe(protein)
    df_residues = featurizer.residues_dataframe(protein)
    df_atoms = featurizer.atoms_dataframe(protein)

    print(df_atoms.head())

def ml_dataframe2():
    from protkit.file_io import ProtIO
    from protkit.ml import ProteinToDataframe

    featurizer = ProteinToDataframe()
    protein1 = ProtIO.load("1ahw.prot")[0]
    protein1.pdb_id = "1ahw"
    protein2 = ProtIO.load("3i40.prot")[0]
    protein2.pdb_id = "3i40"

    df_proteins = featurizer.proteins_dataframe([protein1, protein2])

    print(df_proteins.head())

def ml_dataframe3():
    from protkit.file_io import ProtIO
    from protkit.ml import ProteinToDataframe
    from protkit.properties import SurfaceArea, Hydrophobicity, Mass

    featurizer = ProteinToDataframe()
    protein = ProtIO.load("1ahw.prot")[0]
    protein.pdb_id = "1ahw"

    # Assign additional properties to the protein
    SurfaceArea.surface_area_of_protein(protein, assign_attribute=True)
    Hydrophobicity.hydrophobicity_of_protein(protein, assign_attribute=True)
    Mass.residue_mass_of_protein(protein, assign_attribute=True)

    # Create a DataFrame for the residues of the protein
    df_residues_all = featurizer.residues_dataframe(protein)
    df_residues_none = featurizer.residues_dataframe(protein, additional_fields=[])
    df_residues_specific = featurizer.residues_dataframe(protein, additional_fields=["surface_area", "hydrophobicity"])

    print(df_residues_all.head())
    print(df_residues_none.head())
    print(df_residues_specific.head())

# # prep_data()
#
# quick_start_example()
#
# download_pdb_example()
# download_fasta_example()
# download_cif_example()
#
# properties_hydrophobicity()
# properties_hydrophobicity_class()
# properties_mass()
# properties_chemical_class()
# properties_charge()
# properties_polarity()
# properties_donors_acceptors()
# properties_surface_area()
# properties_volume()
# properties_volume_class()
# properties_bounds_and_center()
# properties_bond_lengths()
# properties_peptide_lengths()
# properties_bond_angles()
# properties_dihedral_angles()
# properties_circular_variance()  # -> double check first residue
# # properties_interface_atoms()
# # properties_interface_residues()
properties_interface_atom_pairs()
properties_interface_residue_pairs()

# tools_reduce()
# tools_freesasa()

# ml_dataframe()
# ml_dataframe2()
# ml_dataframe3()