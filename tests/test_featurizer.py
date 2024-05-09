def test_featurizer():
    from protkit.ml.dataframe import ProteinToDataframe
    from protkit.file_io import PDBIO, FastaIO
    from protkit.properties.surface_area import SurfaceArea

    file_unprotonated = "data/pdb/rcsb/1ahw.pdb"
    protein = PDBIO.load(file_unprotonated)[0]
    protein.pdb_id = "1ahw"
    protein.set_attribute("note", "This is a test protein")
    SurfaceArea.surface_area_of_protein(protein, assign_attribute=True)

    file_unprotonated2 = "data/pdb/rcsb/1a4y.pdb"
    protein2 = PDBIO.load(file_unprotonated2)[0]
    protein2.pdb_id = "1a4y"
    protein2.set_attribute("note2", "This is a another protein")

    featurizer = ProteinToDataframe()

    df_protein1 = featurizer.proteins_dataframe([protein, protein2])
    print(df_protein1.head().to_string())

    df_protein = featurizer.proteins_dataframe([protein, protein2], additional_fields=["note", "note2", "surface_area"])
    print(df_protein.head().to_string())

    df_chain1 = featurizer.chains_dataframe([protein, protein2])
    print(df_chain1.to_string())
    df_chain = featurizer.chains_dataframe([protein, protein2], additional_fields=["note", "surface_area"])
    # print(df_chain.head().to_string())
    print(df_chain.to_string())

    df_residues1 = featurizer.residues_dataframe([protein, protein2])
    print(df_residues1.head().to_string())

    df_residues = featurizer.residues_dataframe(protein, additional_fields=["note", "surface_area"])
    print(df_residues.head().to_string())

    df_atoms1 = featurizer.atoms_dataframe([protein, protein2])
    print(df_atoms1.head().to_string())

    df_atoms = featurizer.atoms_dataframe(protein, additional_fields=["note", "surface_area"])
    print(df_atoms.head().to_string())
    # print(df_atoms.to_string())

    fasta_file_path = "data/fasta/rcsb/1ahw.fasta"
    sequences = FastaIO.load(fasta_file_path)
    df_sequence = featurizer.sequences_dataframe(sequences)
    print(df_sequence.head().to_string())

test_featurizer()