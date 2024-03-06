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

def download_pdb_example():
    from protkit.download import Download

    # Download a PDB file from the RCSB PDB database and save it to a file.
    Download.download_pdb_file_from_rcsb("1ahw", "data/pdb_files/rcsb/1ahw.pdb")

    # Download a PDB file from the SAbDab database and save it to a file.
    Download.download_pdb_file_from_sabdab("1ahw", "data/pdb_files/sabdab/1ahw.pdb")

    # Download multiple PDB files from the RCSB PDB in parallel and save them to a directory.
    Download.download_pdb_files_from_rcsb(["1ahw", "1a4y", "1a6m"], "data/pdb_files/rcsb/", n_jobs=3)

def download_fasta_example():
    from protkit.download import Download

    # Download a Fasta file from the Uniprot database and save it to a file.
    Download.download_fasta_file_from_uniprot("P12345", "data/fasta_files/uniprot/P12345.xml")

    # Download multiple Fasta files from the RCSB database in parallel and save them to a directory.
    Download.download_fasta_files_from_rcsb(["P12345", "P12346", "P12347"], "data/fasta_files/rcsb/")

    # Download multiple Fasta files from the Uniprot database in parallel and save them to a directory.
    Download.download_fasta_files_from_uniprot(["P12345", "P12346", "P12347"], "data/fasta_files/rcsb/")


def download_cif_example():
    from protkit.download import Download

    # Download a CIF file from the RCSB PDB database and save it to a file.
    Download.download_cif_file_from_rcsb("1ahw", "data/cif_files/rcsb/1ahw.cif")

    # Download multiple CIF files from the RCSB PDB database in parallel and save them to a directory.
    Download.download_cif_files_from_rcsb(["1ahw", "1a4y", "1a6m"], "data/cif_files/rcsb/")

    # Download a binary CIF file from the RCSB PDB database and save it to a file.
    Download.download_binary_cif_file_from_rcsb("1ahw", "data/cif_files/rcsb/1ahw.bcif")

    # Download multiple binary CIF files from the RCSB PDB database in parallel and save them to a directory.
    Download.download_binary_cif_files_from_rcsb(["1ahw", "1a4y", "1a6m"], "data/cif_files/rcsb/")


quick_start_example()
# download_pdb_example()
# download_fasta_example()
# download_cif_example()

