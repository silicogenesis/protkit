#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Authors:  Fred Senekal (FS)
# Contact:  fred@silicogenesis.com
# License:  GPLv3

"""
Implements the `Download` class to download biological data from the internet.

Currently, downloading of the following data is supported:

- PDB files from the RCSB and Sabdab
- CIF files from the RCSB
- Binary CIF files from the RCSB
- FASTA files from the RCSB and Uniprot

For more information about the various data sources, see the following URLs:

- RCSB: https://www.rcsb.org/
- Uniprot: https://www.uniprot.org/
- Sabdab: https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/

A list of download services provided by RCSB is available at:

- https://www.rcsb.org/docs/programmatic-access/file-download-services

"""
import os.path

import requests
from joblib import Parallel, delayed
from typing import List, Union


class Download:
    """
    Class `Download` is a base class for downloading
    biological data from the internet.
    """

    @staticmethod
    def download_file(url: str, file_path: str) -> None:
        """
        Downloads a single file from the internet.

        Uses the `requests` library to download the file.

        Args:
            url (str): The URL of the file.
            file_path (str): The path to the file.

        Returns:
            None

        Raises:
            Exception: If the file could not be downloaded or saved.
        """
        try:
            response = requests.get(url)
            if response.status_code == 200:
                # Create the directory if it does not exist.
                directory = os.path.dirname(file_path)
                if directory != "":
                    if not os.path.exists(directory):
                        os.makedirs(directory)

                # Write the file to disk.
                with open(file_path, "wb") as file:
                    file.write(response.content)
        except Exception as e:
            raise e

    @staticmethod
    def _download_file(url: str, file_path: str) -> None:
        """
        Downloads a single file from the internet.

        Uses the `urllib` library to download the file.

        Args:
            url (str): The URL of the file.
            file_path (str): The path to the file.

        Returns:
            None

        Raises:
            Exception: If the file could not be downloaded or saved.
        """
        import urllib.request
        try:
            with urllib.request.urlopen(url) as response:
                content = response.read()

                with open(file_path, 'wb') as file:
                    file.write(content)
        except Exception as e:
            raise e

    @staticmethod
    def parallel_download(
            urls: List[str],
            file_paths: List[str],
            n_jobs: int = -1) -> None:
        """
        Downloads multiple files from the internet in parallel.

        Args:
            urls (List[str]): The URLs of the files.
            file_paths (List[str]): The paths to the files.
            n_jobs (int): The number of jobs to run in parallel.
                If -1, the number of jobs is set to the number of CPU cores.

        Returns:
            None

        Raises:
            Exception: If the files could not be downloaded or saved.
        """
        Parallel(n_jobs=n_jobs)(delayed(Download.download_file)(url, file_path)
                                for url, file_path in zip(urls, file_paths))

    @staticmethod
    def download_fasta_file_from_rcsb(
            pdb_id: str,
            file_path_or_directory: str) -> None:
        """
        Downloads a single FASTA file from the RCSB.

        Args:
            pdb_id (str): The ID of the PDB file.
            file_path_or_directory (str): The path where the FASTA file should be saved.
                If the path is a directory, the file is saved in the directory
        """
        url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
        if os.path.isdir(file_path_or_directory):
            file_path_or_directory = os.path.join(file_path_or_directory, f"{pdb_id}.fasta")

        Download.download_file(url, file_path_or_directory)

    @staticmethod
    def download_fasta_files_from_rcsb(
            pdb_ids: List[str],
            directory: str,
            n_jobs: int = -1) -> None:
        """
        Downloads multiple FASTA files from the RCSB.

        Args:
            pdb_ids (List[str]): The IDs of the PDB files.
            directory (str): The directory where the FASTA files should be saved.
            n_jobs (int): The number of jobs to run in parallel.
                If -1, the number of jobs is set to the number of CPU cores.

        Returns:
            None
        """
        urls = [f"https://www.rcsb.org/fasta/entry/{pdb_id}" for pdb_id in pdb_ids]
        file_paths = [os.path.join(directory, f"{pdb_id}.fasta") for pdb_id in pdb_ids]
        Download.parallel_download(urls, file_paths, n_jobs=n_jobs)

    @staticmethod
    def download_fasta_file_from_uniprot(
            uniprot_id: str,
            file_path_or_directory: str) -> None:
        """
        Downloads a single FASTA file from Uniprot.

        Args:
            uniprot_id (str): The ID of the UniProt file.
            file_path_or_directory (str): The path where the FASTA file should be saved.
                If the path is a directory, the file is saved in the directory
        """
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
        if os.path.isdir(file_path_or_directory):
            file_path_or_directory = os.path.join(file_path_or_directory, f"{uniprot_id}.fasta")

        Download.download_file(url, file_path_or_directory)

    @staticmethod
    def download_fasta_files_from_uniprot(
            uniprot_ids: List[str],
            directory: str,
            n_jobs: int = -1) -> None:
        """
        Downloads multiple FASTA files from Uniprot.

        Args:
            uniprot_ids (List[str]): The IDs of the UniProt files.
            directory (str): The directory where the FASTA files should be saved.
            n_jobs (int): The number of jobs to run in parallel.
                If -1, the number of jobs is set to the number of CPU cores.
        Returns:
            None
        """
        urls = [f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta" for uniprot_id in uniprot_ids]
        file_paths = [os.path.join(directory, f"{uniprot_id}.fasta") for uniprot_id in uniprot_ids]
        Download.parallel_download(urls, file_paths, n_jobs=n_jobs)

    @staticmethod
    def download_pdb_file_from_rcsb(
            pdb_id: str,
            file_path_or_directory: str) -> None:
        """
        Downloads a single PDB file from the RCSB.

        Args:
            pdb_id (str): The ID of the PDB file.
            file_path_or_directory (str): The path where the PDB file should be saved.
                If the path is a directory, the file is saved in the directory
                with the name <pdb_id>.pdb.
        """
        download_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        if os.path.isdir(file_path_or_directory):
            file_path_or_directory = os.path.join(file_path_or_directory, f"{pdb_id}.pdb")

        Download.download_file(download_url, file_path_or_directory)

    @staticmethod
    def download_pdb_files_from_rcsb(
            pdb_ids: List[str],
            directory: str,
            n_jobs: int = -1) -> None:
        """
        Downloads multiple PDB files from the RCSB.

        Args:
            pdb_ids (List[str]): The IDs of the PDB files.
            directory (str): The directory where the PDB files should be saved.
            n_jobs (int): The number of jobs to run in parallel.
                If -1, the number of jobs is set to the number of CPU cores.

        Returns:
            None
        """
        urls = [f"https://files.rcsb.org/download/{pdb_id}.pdb" for pdb_id in pdb_ids]
        file_paths = [os.path.join(directory, f"{pdb_id}.pdb") for pdb_id in pdb_ids]
        Download.parallel_download(urls, file_paths, n_jobs=n_jobs)

    @staticmethod
    def download_pdb_file_from_sabdab(
            pdb_id: str,
            file_path_or_directory: str) -> None:
        """
        Downloads a single PDB file from Sabdab.

        Args:
            pdb_id (str): The ID of the PDB file.
            file_path_or_directory (str): The path where the PDB file should be saved.
                If the path is a directory, the file is saved in the directory
                with the name <pdb_id>.pdb.
        """
        download_url = f"https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/pdb/{pdb_id.lower()}/?raw=true"
        if os.path.isdir(file_path_or_directory):
            file_path_or_directory = os.path.join(file_path_or_directory, f"{pdb_id}.pdb")

        Download.download_file(download_url, file_path_or_directory)

    @staticmethod
    def download_pdb_files_from_sabdab(
            pdb_ids: List[str],
            directory: str,
            n_jobs: int = -1) -> None:
        """
        Downloads multiple PDB files from Sabdab.

        Args:
            pdb_ids (List[str]): The IDs of the PDB files.
            directory (str): The directory where the PDB files should be saved.
            n_jobs (int): The number of jobs to run in parallel.
                If -1, the number of jobs is set to the number of CPU cores.

        Returns:
            None
        """
        urls = [f"https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/pdb/{pdb_id.lower()}/?raw=true" for pdb_id in pdb_ids]
        file_paths = [os.path.join(directory, f"{pdb_id}.pdb") for pdb_id in pdb_ids]
        Download.parallel_download(urls, file_paths, n_jobs=n_jobs)

    @staticmethod
    def download_cif_file_from_rcsb(
            pdb_id: str,
            file_path_or_directory: str) -> None:
        """
        Downloads a single CIF file from the RCSB.

        Args:
            pdb_id (str): The ID of the CIF file.
            file_path_or_directory (str): The path where the CIF file should be saved.
                If the path is a directory, the file is saved in the directory
                with the name <pdb_id>.cif.
        """
        download_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        if os.path.isdir(file_path_or_directory):
            file_path_or_directory = os.path.join(file_path_or_directory, f"{pdb_id}.cif")

        Download.download_file(download_url, file_path_or_directory)

    @staticmethod
    def download_cif_files_from_rcsb(
            pdb_ids: List[str],
            directory: str,
            n_jobs: int = -1) -> None:
        """
        Downloads multiple CIF files from the RCSB.

        Args:
            pdb_ids (List[str]): The IDs of the CIF files.
            directory (str): The directory where the CIF files should be saved.
            n_jobs (int): The number of jobs to run in parallel.
                If -1, the number of jobs is set to the number of CPU cores.

        Returns:
            None
        """
        urls = [f"https://files.rcsb.org/download/{pdb_id}.cif" for pdb_id in pdb_ids]
        file_paths = [os.path.join(directory, f"{pdb_id}.cif") for pdb_id in pdb_ids]
        Download.parallel_download(urls, file_paths, n_jobs=n_jobs)


    @staticmethod
    def download_binary_cif_file_from_rcsb(
            pdb_id: str,
            file_path_or_directory: str) -> None:
        """
        Downloads a single Binary CIF file from the RCSB.

        Args:
            pdb_id (str): The ID of the Binary CIF file.
            file_path_or_directory (str): The path where the Binary CIF file should be saved.
                If the path is a directory, the file is saved in the directory
                with the name <pdb_id>.bcif.
        """
        download_url = f"https://models.rcsb.org/{pdb_id}.bcif"
        if os.path.isdir(file_path_or_directory):
            file_path_or_directory = os.path.join(file_path_or_directory, f"{pdb_id}.bcif")

        Download.download_file(download_url, file_path_or_directory)

    @staticmethod
    def download_binary_cif_files_from_rcsb(
            pdb_ids: List[str],
            directory: str,
            n_jobs: int = -1) -> None:
        """
        Downloads multiple Bianry CIF files from the RCSB.

        Args:
            pdb_ids (List[str]): The IDs of the Binary CIF files.
            directory (str): The directory where the CIF files should be saved.
            n_jobs (int): The number of jobs to run in parallel.
                If -1, the number of jobs is set to the number of CPU cores.

        Returns:
            None
        """
        urls = [f"https://models.rcsb.org/{pdb_id}.bcif" for pdb_id in pdb_ids]
        file_paths = [os.path.join(directory, f"{pdb_id}.bcif") for pdb_id in pdb_ids]
        Download.parallel_download(urls, file_paths, n_jobs=n_jobs)


