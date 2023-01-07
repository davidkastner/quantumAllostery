"""Library of data patterns."""


def get_aa_identifiers() -> dict[str, list[str, str, str, float]]:
    """
    Amino acid naming conventions and basic features.

    Contains the full name, the one letter code, three-letter code,
    molecular formula, and molecular weight.

    Returns
    -------
    aa : dict[str, list[str, str, str, float]]
        Dictionary containing the name, one letter code, three-letter code, molecular formula, and molecular weight.

    """

    aa = {
    "alanine" : ["A", "ALA", "C3H5NO", 71.08],
    "arginine" : ["R", "ARG", "C6H12N4O", 156.19],
    "asparagine" : ["N", "ASN", "C4H6N2O2", 114.11],
    "aspartate" : ["D", "ASP", "C4H5NO3", 115.09],
    "cysteine" : ["C", "CYS", "C3H5NOS", 103.15],
    "glutamate" : ["E", "GLU", "C5H7NO3", 129.12],
    "glutamine" : ["Q", "GLN", "C5H8N2O2", 128.13],
    "glycine" : ["G", "GLY", "C2H3NO", 57.05],
    "histidine" : ["H", "HIS", "C6H7N3O", 137.14],
    "isoleucine" : ["I", "ILE", "C6H11NO", 113.16],
    "leucine" : ["L", "LEU", "C6H11NO", 113.16],
    "lysine" : ["K", "LYS", "C6H12N2O", 128.18],
    "methionine" : ["M", "MET", "C5H9NOS", 131.20],
    "phenylalanine" : ["F", "PHE", "C9H9NO", 147.18],
    "proline" : ["P", "PRO", "C5H7NO", 97.12],
    "serine" : ["S", "SER", "C3H5NO2", 87.08],
    "threonine" : ["T", "THR", "C4H7NO2", 101.11],
    "tryptophan" : ["W", "TRP", "C11H10N2O	", 186.22],
    "tyrosine" : ["Y", "TYR", "C9H9NO2", 163.18],
    "valine" : ["V", "VAL", "C5H9NO", 99.13],
    "2-aminoisobutyrate" : ["U", "AIB", "C4H7NO", 85.12],
    "acetyl" : ["na", "ACE", "C2H3O", 43]
    }

    return aa

def sequence(protein) -> list[str]:
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

    return seq
