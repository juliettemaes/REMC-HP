"""
Utility functions for handling U-moves in a lattice-based amino acid chain system.
"""

from typing import Dict, Tuple, Optional
from moves.move_utils import relative_position, check_occupancy


def get_cks_aa_to_check(lattice, chain_index: int) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, int], Dict[str, int]]:
    """
    Retrieve amino acid coordinates needed to check U-moves based on the chain index.

    Args:
        lattice: The lattice object containing sequence and grid information.
        chain_index (int): The index of the amino acid in the chain.

    Returns:
        Tuple[Dict[str, int], Dict[str, int], Dict[str, int], Dict[str, int]]: 
        The coordinates of the amino acid at `chain_index - 1`, `chain_index - 2`, `chain_index`, 
        and `chain_index + 1` respectively.
    """
    if chain_index < lattice.sequence.length - 1:
        res_i = lattice.sequence.aa_coord[chain_index - 1]
        res_iminus1 = lattice.sequence.aa_coord[chain_index - 2]
        res_iplus1 = lattice.sequence.aa_coord[chain_index]
        res_next2 = lattice.sequence.aa_coord[chain_index + 1]
    elif chain_index == lattice.sequence.length - 1:
        res_i = lattice.sequence.aa_coord[chain_index - 2]
        res_iminus1 = lattice.sequence.aa_coord[chain_index - 3]
        res_iplus1 = lattice.sequence.aa_coord[chain_index - 1]
        res_next2 = lattice.sequence.aa_coord[chain_index]
    return res_i, res_iminus1, res_iplus1, res_next2


def get_alternative_positions(lattice, chain_index: int) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, int], Dict[str, int]]:
    """
    Retrieve alternative amino acid positions for U-move validation.

    Args:
        lattice: The lattice object containing sequence and grid information.
        chain_index (int): The index of the amino acid in the chain.

    Returns:
        Tuple[Dict[str, int], Dict[str, int], Dict[str, int], Dict[str, int]]: 
        The coordinates of `chain_index - 3`, `chain_index - 2`, `chain_index - 1`, and `chain_index`.
    """
    res_prev2 = lattice.sequence.aa_coord[chain_index - 3]
    res_iminus1 = lattice.sequence.aa_coord[chain_index - 2]
    res_i = lattice.sequence.aa_coord[chain_index - 1]
    res_iplus1 = lattice.sequence.aa_coord[chain_index]
    return res_prev2, res_iminus1, res_i, res_iplus1


def are_U(res1: Dict[str, int], res2: Dict[str, int], res3: Dict[str, int], res4: Dict[str, int]) -> bool:
    """
    Check if four amino acids form a U-structure on the lattice.

    Args:
        res1 (Dict[str, int]): Coordinates of the first amino acid.
        res2 (Dict[str, int]): Coordinates of the second amino acid.
        res3 (Dict[str, int]): Coordinates of the third amino acid.
        res4 (Dict[str, int]): Coordinates of the fourth amino acid.

    Returns:
        bool: True if the amino acids form a U-structure, False otherwise.
    """
    cdt1 = abs(res1['x'] - res2['x']) + abs(res1['y'] - res2['y']) == 1
    cdt2 = abs(res2['x'] - res3['x']) + abs(res2['y'] - res3['y']) == 1
    cdt3 = abs(res3['x'] - res4['x']) + abs(res3['y'] - res4['y']) == 1
    cdt4 = abs(res4['x'] - res1['x']) + abs(res4['y'] - res1['y']) == 1
    return (cdt1 and cdt2 and cdt3 and cdt4)


def find_potential_U(res1: Dict[str, int], res2: Dict[str, int], res3: Dict[str, int], res4: Dict[str, int]) -> Tuple[int, int, int, int]:
    """
    Calculate the potential new coordinates for residues after a U-move.

    Args:
        res1 (Dict[str, int]): Coordinates of the first amino acid.
        res2 (Dict[str, int]): Coordinates of the second amino acid, to be moved.
        res3 (Dict[str, int]): Coordinates of the third amino acid, to be moved.
        res4 (Dict[str, int]): Coordinates of the fourth amino acid.

    Returns:
        Tuple[int, int, int, int]: The new (x, y) coordinates for the amino acids after a potential U-move.
    """
    x_new1 = res1['x'] - relative_position(res2, res1)[0]
    y_new1 = res1['y'] - relative_position(res2, res1)[1]

    x_new2 = res4['x'] - relative_position(res3, res4)[0]
    y_new2 = res4['y'] - relative_position(res3, res4)[1]

    return x_new1, y_new1, x_new2, y_new2


def validate_U(lattice, res_iminus1: Dict[str, int], res_i: Dict[str, int], res_iplus1: Dict[str, int], res_next2: Dict[str, int], chain_index: int) -> Tuple[bool, Optional[int], Optional[int], Optional[int], Optional[int]]:
    """
    Validate whether the amino acids form a U-structure and if the U-move is possible.

    Args:
        lattice: The lattice object containing sequence and grid information.
        res_iminus1 (Dict[str, int]): Coordinates of the amino acid at `chain_index - 2`.
        res_i (Dict[str, int]): Coordinates of the amino acid at `chain_index - 1`.
        res_iplus1 (Dict[str, int]): Coordinates of the amino acid at `chain_index`.
        res_next2 (Dict[str, int]): Coordinates of the amino acid at `chain_index + 1`.
        chain_index (int): The index of the amino acid in the chain.

    Returns:
        Tuple[bool, Optional[int], Optional[int], Optional[int], Optional[int]]: 
        A tuple with a boolean indicating if the U-move is possible and the new (x, y) coordinates for the move if valid.
    """
    if are_U(res_iminus1, res_i, res_iplus1, res_next2) and chain_index < lattice.sequence.length - 1:
        x_new1, y_new1, x_new2, y_new2 = find_potential_U(res_iminus1, res_i, res_iplus1, res_next2)
        if check_occupancy(lattice, x_new1, y_new1) and check_occupancy(lattice, x_new2, y_new2):
            return True, x_new1, y_new1, x_new2, y_new2
        else:
            return False, x_new1, y_new1, x_new2, y_new2
    return False, None, None, None, None


def execute_u_move(lattice, chain_index: int, res_i: Dict[str, int], res_iplus1: Dict[str, int], x_new1: int, y_new1: int, x_new2: int, y_new2: int) -> None:
    """
    Execute a U-move by updating the lattice grid and amino acid coordinates.
    The residues at `chain_index - 1` and `chain_index` are moved to new positions.

    Args:
        lattice: The lattice object containing sequence and grid information.
        chain_index (int): The index of the amino acid to be moved.
        res_i (Dict[str, int]): The current coordinates of the amino acid at `chain_index - 1`.
        res_iplus1 (Dict[str, int]): The current coordinates of the amino acid at `chain_index`.
        x_new1 (int): The new x-coordinate for the amino acid at `chain_index - 1`.
        y_new1 (int): The new y-coordinate for the amino acid at `chain_index - 1`.
        x_new2 (int): The new x-coordinate for the amino acid at `chain_index`.
        y_new2 (int): The new y-coordinate for the amino acid at `chain_index`.

    Returns:
        None
    """
    lattice.lattice[x_new1, y_new1] = chain_index
    lattice.lattice[x_new2, y_new2] = chain_index + 1
    lattice.lattice[res_i["x"], res_i["y"]] = 0
    lattice.lattice[res_iplus1["x"], res_iplus1["y"]] = 0
    lattice.sequence.aa_coord_update(chain_index - 1, x_new1, y_new1)
    lattice.sequence.aa_coord_update(chain_index, x_new2, y_new2)


def execute_alternative_u_move(
    lattice, 
    chain_index: int, 
    res_i: Dict[str, int], 
    res_iminus1: Dict[str, int], 
    x_new1: int, 
    y_new1: int, 
    x_new2: int, 
    y_new2: int
) -> None:
    """
    Execute an alternative U-move, relocating two amino acids in the lattice.
    The residues at `chain_index - 2` and `chain_index - 1` are moved to new positions.

    Args:
        lattice: The lattice object containing sequence and grid information.
        chain_index (int): The index of the amino acid in the chain.
        res_i (Dict[str, int]): Coordinates of the amino acid at `chain_index - 1`.
        res_iminus1 (Dict[str, int]): Coordinates of the amino acid at `chain_index - 2`.
        x_new1 (int): The new x-coordinate for the amino acid at `chain_index - 2`.
        y_new1 (int): The new y-coordinate for the amino acid at `chain_index - 2`.
        x_new2 (int): The new x-coordinate for the amino acid at `chain_index - 1`.
        y_new2 (int): The new y-coordinate for the amino acid at `chain_index - 1`.

    Returns:
        None
    """
    # Move the amino acids to the new positions
    lattice.lattice[x_new1, y_new1] = chain_index - 1
    lattice.lattice[x_new2, y_new2] = chain_index

    # Free the previous positions of the amino acids
    lattice.lattice[res_iminus1["x"], res_iminus1["y"]] = 0
    lattice.lattice[res_i["x"], res_i["y"]] = 0

    # Update the coordinates of the amino acids in the sequence object
    lattice.sequence.aa_coord_update(chain_index - 2, x_new1, y_new1)
    lattice.sequence.aa_coord_update(chain_index - 1, x_new2, y_new2)
