"""Data, file, and design patterns"""

import os
import time
from typing import List, Dict


def get_aa_identifiers() -> Dict[str, List[str]]:
    """
    Amino acid naming conventions and basic features.
    Contains the full name, the one letter code, three-letter code,
    molecular formula, and molecular weight.
    Returns
    -------
    aa : dict[str, list[str, str, str]]
        Dictionary containing the name, one letter code, three-letter code, and molecular formula.
    """

    aa = {
        "alanine": ["A", "ALA"],
        "arginine": ["R", "ARG"],
        "asparagine": ["N", "ASN"],
        "aspartate": ["D", "ASP"],
        "cysteine": ["C", "CYS"],
        "glutamate": ["E", "GLU"],
        "glutamine": ["Q", "GLN"],
        "glycine": ["G", "GLY"],
        "histidine": ["H", "HIS"],
        "isoleucine": ["I", "ILE"],
        "leucine": ["L", "LEU"],
        "lysine": ["K", "LYS"],
        "methionine": ["M", "MET"],
        "phenylalanine": ["F", "PHE"],
        "proline": ["P", "PRO"],
        "serine": ["S", "SER"],
        "threonine": ["T", "THR"],
        "tryptophan": ["W", "TRP"],
        "tyrosine": ["Y", "TYR"],
        "valine": ["V", "VAL"],
        "2-aminoisobutyrate": ["U", "AIB"],
        "acetyl": ["na", "ACE"],
    }

    return aa


def sequence(protein) -> List[str]:
    """
    Returns the sequence of one of the mimochromes.
    Add additonal sequences to this function as lists.
    Parameters
    ----------
    protein : str
        The name of the protei/peptide whose sequence has been requested.
    Returns
    -------
    seq : list[str]
        A list of the amino acids in the requested protein/peptide.
    """

    # Sequence of mimochrome MC6
    if protein == "mc6":
        seq = [
            "ASP1",
            "GLU2",
            "GLN3",
            "GLN4",
            "LEU5",
            "HIS6",
            "SER7",
            "GLN8",
            "LYS9",
            "ARG10",
            "LYS11",
            "ILE12",
            "THR13",
            "LEU14",
            "ASP15",
            "GLU16",
            "GLN17",
            "GLN18",
            "LEU19",
            "SER20",
            "SER21",
            "GLN22",
            "LYS23",
            "ARG24",
        ]

    # Sequence of mimochrome MC6*
    if protein == "mc6s":
        seq = [
            "ASP1",
            "LEU2",
            "GLN3",
            "GLN4",
            "LEU5",
            "HIS6",
            "SER7",
            "GLN8",
            "LYS9",
            "ARG10",
            "LYS11",
            "ILE12",
            "THR13",
            "LEU14",
            "ASP15",
            "GLU16",
            "GLN17",
            "GLN18",
            "LEU19",
            "SER20",
            "SER21",
            "GLN22",
            "LYS23",
            "ARG24",
        ]

    # Sequence of mimochrome MC6*a
    if protein == "mc6sa":
        seq = [
            "ASP1",
            "LEU2",
            "GLN3",
            "GLN4",
            "LEU5",
            "HIS6",
            "SER7",
            "GLN8",
            "LYS9",
            "ARG10",
            "LYS11",
            "ILE12",
            "THR13",
            "LEU14",
            "ASP15",
            "GLU16",
            "AIB17",
            "GLN18",
            "LEU19",
            "AIB20",
            "SER21",
            "GLN22",
            "LYS23",
            "ARG24",
        ]

    # All amino acids in the sequence of MC6
    if protein == "mc6_all":
        seq = [
            "ACE1",
            "ASP2",
            "GLU3",
            "GLN4",
            "GLN5",
            "LEU6",
            "HIS7",
            "SER8",
            "GLN9",
            "LYS10",
            "ARG11",
            "LYS12",
            "ILE13",
            "THR14",
            "LEU15",
            "NHE16",
            "ACE17",
            "ASP18",
            "GLU19",
            "GLN20",
            "GLN21",
            "LEU22",
            "SER23",
            "SER24",
            "GLN25",
            "LYS26",
            "ARG27",
            "NHE28",
            "HEME",
            "FE",
        ]

    if protein == "2jof":
        seq = [
            "ASP1",
            "ALA2",
            "TYR3",
            "ALA4",
            "GLN5",
            "TRP6",
            "LEU7",
            "LYS8",
            "ASP9",
            "GLY10",
            "GLY11",
            "PRO12",
            "SER13",
            "SER14",
            "GLY15",
            "ARG16",
            "PRO17",
            "PRO18",
            "PRO19",
            "SER20",
        ]
    return seq
