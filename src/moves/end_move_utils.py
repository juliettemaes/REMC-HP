"""
Functions for handling end moves in a lattice-based amino acid chain system. 
"""

from typing import List, Tuple
import random as rd

def is_end_move_possible(lattice, chain_index: int) -> int:
    """
    Check if the amino acid is at the start or end of the chain, 
    and return the index of the residue next to it. 

    Args:
        lattice: The lattice object containing sequence and grid information.
        chain_index (int): The index of the amino acid in the chain.

    Returns:
        int: The index of the chain next to the amino acid, to latter evaluate 
        the move possibilities.
             Returns 2 if at the beginning, or the sequence length minus one if at the end.

    Raises:
        ValueError: If the amino acid is not at the start or end of the chain.
    """
    if chain_index == 1:  # The amino acid is at the beginning of the chain
        return 2
    elif chain_index == lattice.sequence.length:  # The amino acid is at the end of the chain
        return lattice.sequence.length - 1
    else:
        raise ValueError("Impossible move. The amino acid is not at the end of the chain")


def execute_end_move(lattice, possible_pos: List[Tuple[int, int]], chain_index: int, old_x: int, old_y: int) -> None:
    """
    Execute a move for an amino acid at the end of the chain to a new position in the lattice.

    Args:
        lattice: The lattice object containing sequence and grid information.
        possible_pos (List[Tuple[int, int]]): A list of possible (x, y) coordinates to move to.
        chain_index (int): The index of the amino acid to be moved.
        old_x (int): The current x-coordinate of the amino acid.
        old_y (int): The current y-coordinate of the amino acid.

    Returns:
        None
    """
    new_pos = rd.choice(possible_pos)
    x_new, y_new = new_pos[0], new_pos[1]
    lattice.lattice[old_x, old_y] = 0
    lattice.lattice[x_new, y_new] = chain_index
    lattice.sequence.aa_coord_update(chain_index - 1, x_new, y_new)  # Update the coordinates of the amino acid in the sequence