"""
Package `download` contains classes to download
biological data from the internet.

- Class `Download` in module `download` is used to download
    files such as PDB and Fasta file from the internet. It
    can be used to download a single file or multiple files
    in parallel. It downloads files from sources such as
    RCSB PDB, UniProt, and SAbDab.
"""

from .download import Download
