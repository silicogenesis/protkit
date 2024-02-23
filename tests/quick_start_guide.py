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


# download_pdb_example()
# download_fasta_example()
download_cif_example()