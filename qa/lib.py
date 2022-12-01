"""Library of data patterns."""


def sequence(protein) -> list[str]:
    """
    Generates mutual information and corss-correlation matrices.

    Add additonal sequences to this function as lists.

    """

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
