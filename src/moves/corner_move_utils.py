"""
Utility functions for handling corner movements in a lattice-based system.
"""

from typing import Tuple, Dict

def find_potential_corner(lattice, chain_index: int) -> Tuple[int, int]:
    """
    Find potential new corner coordinates for a residue in the sequence.

    Args:
        lattice: The lattice object containing sequence and grid information.
        chain_index (int): The index of the amino acid in the chain.

    Returns:
        Tuple[int, int]: The new (x, y) coordinates for the potential corner position.
    """
    if lattice.sequence.aa_coord[chain_index - 1]["x"] == lattice.sequence.aa_coord[chain_index - 2]["x"]:
        y_new = lattice.sequence.aa_coord[chain_index - 2]["y"]
        x_new = lattice.sequence.aa_coord[chain_index]["x"]
    else:
        x_new = lattice.sequence.aa_coord[chain_index - 2]["x"]
        y_new = lattice.sequence.aa_coord[chain_index]["y"]
    return x_new, y_new


def get_corner_aa_to_check(lattice, chain_index: int) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, int]]:
    """
    Retrieve residue coordinates needed to check for potential corner moves.

    Args:
        lattice: The lattice object containing sequence and grid information.
        chain_index (int): The index of the amino acid in the chain.

    Returns:
        Tuple[Dict[str, int], Dict[str, int], Dict[str, int]]: The coordinates of the amino acid at `chain_index - 1`, 
                                                               `chain_index - 2`, and `chain_index` respectively.
    """
    res_i = lattice.sequence.aa_coord[chain_index - 1]
    res_iminus1 = lattice.sequence.aa_coord[chain_index - 2]
    res_iplus1 = lattice.sequence.aa_coord[chain_index]
    return res_i, res_iminus1, res_iplus1


def execute_corner_move(lattice, chain_index: int, res_i: Dict[str, int], x_new: int, y_new: int) -> None:
    """
    Execute the corner move by updating the lattice grid and amino acid coordinates.

    Args:
        lattice: The lattice object containing sequence and grid information.
        chain_index (int): The index of the amino acid to be moved.
        res_i (Dict[str, int]): The current coordinates of the amino acid.
        x_new (int): The new x-coordinate for the amino acid.
        y_new (int): The new y-coordinate for the amino acid.

    Returns:
        None
    """
    lattice.lattice[x_new, y_new] = chain_index
    lattice.lattice[res_i["x"], res_i["y"]] = 0
    lattice.sequence.aa_coord_update(chain_index - 1, x_new, y_new)  # Update the coordinates of the amino acid in the sequence
