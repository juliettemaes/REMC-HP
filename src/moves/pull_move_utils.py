"""
Utility functions for performing pull moves in a lattice-based protein structure simulation.
"""

from typing import Dict, List, Tuple
from moves.move_utils import find_empty_diagonal, find_empty_neighbors, are_topological_neighbors, move, free_position


def get_pull_aa_to_check(lattice, chain_index: int) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, int], Dict[str, int]]:
    """
    Retrieve the amino acids that need to be checked for a pull move.

    Args:
        lattice: The lattice object containing sequence and grid information.
        chain_index (int): The index of the amino acid in the chain.

    Returns:
        Tuple containing dictionaries for the amino acids at indices (chain_index - 3),
        (chain_index - 2), (chain_index - 1), and (chain_index).
    """
    res_i = lattice.sequence.aa_coord[chain_index - 1]  # equals to CHAIN INDEX
    res_iplus1 = lattice.sequence.aa_coord[chain_index]
    res_iminus1 = lattice.sequence.aa_coord[chain_index - 2]  # equals to CHAIN INDEX -1
    res_iminus2 = lattice.sequence.aa_coord[chain_index - 3]  # equals to CHAIN INDEX -2
    return res_iminus2, res_iminus1, res_i, res_iplus1

def find_L_positions(lattice, res_i: Dict[str, int], res_i1: Dict[str, int]) -> Tuple[Tuple[int, int], bool]:
    """
    Find possible positions for the 'L' amino acid in a pull move.

    Args:
        lattice: The lattice object containing sequence and grid information.
        res_i (Dict[str, int]): Coordinates of the amino acid at `chain_index - 1`.
        res_i1 (Dict[str, int]): Coordinates of the amino acid at `chain_index`.

    Returns:
        A tuple where the first element is the position of 'L' and the second element is a boolean
        indicating if a valid position was found.
    """
    diag_pos = find_empty_diagonal(lattice, res_i["x"], res_i["y"])
    neighbors_pos = find_empty_neighbors(lattice, res_i1["x"], res_i1["y"])
    intersection = list(set(diag_pos) & set(neighbors_pos))
    if len(intersection) > 0:
        return intersection[0], True  # Keep the first position for now
    else:
        return (None, None), False

def find_C_positions(lattice, L: Tuple[int, int], res_i: Dict[str, int], res_iminus1: Dict[str, int]) -> Tuple[Tuple[int, int], str, bool]:
    """
    Find possible positions for the 'C' amino acid in a pull move.

    Args:
        lattice: The lattice object containing sequence and grid information.
        L (Tuple[int, int]): Coordinates of the 'L' amino acid.
        res_i (Dict[str, int]): Coordinates of the amino acid at `chain_index - 1`.
        res_iminus1 (Dict[str, int]): Coordinates of the amino acid at `chain_index - 2`.

    Returns:
        A tuple where the first element is the position of 'C', the second element is a string
        describing the occupancy status of 'C', and the third element is a boolean indicating if a valid position was found.
    """
    neighbors_pos_L = find_empty_neighbors(lattice, L[0], L[1])
    neighbors_pos_i = find_empty_neighbors(lattice, res_i["x"], res_i["y"])
    possible_pos = list(set(neighbors_pos_L) & set(neighbors_pos_i))  # Find the common empty neighbors
    if len(possible_pos) > 0:
        return possible_pos[0], "empty", True  # 'C' is empty
    elif are_topological_neighbors({"x": L[0], "y": L[1]}, res_iminus1):  # Check if 'C' is occupied by res_iminus1
        return [(res_iminus1["x"], res_iminus1["y"])], "C is res_iminus1", True
    else:
        return (None, None), "C is occupied", False

def execute_pull_move(lattice, chain_index: int, res_i: Dict[str, int], res_iminus1: Dict[str, int], pos_L: Tuple[int, int], pos_C: Tuple[int, int]) -> None:
    """
    Execute a pull move by moving two amino acids to new positions and updating their coordinates.

    Args:
        lattice: The lattice object containing sequence and grid information.
        chain_index (int): The index of the amino acid in the chain.
        res_i (Dict[str, int]): Coordinates of the amino acid at `chain_index - 1`.
        res_iminus1 (Dict[str, int]): Coordinates of the amino acid at `chain_index - 2`.
        pos_L (Tuple[int, int]): New position for the 'L' amino acid.
        pos_C (Tuple[int, int]): New position for the 'C' amino acid.

    Returns:
        None
    """
    # Move `res_i` and `res_iminus1` to positions `pos_L` and `pos_C`, and free the previous positions
    move(lattice, res_i, pos_L)
    move(lattice, res_iminus1, pos_C)

    free_position(lattice, res_i)
    free_position(lattice, res_iminus1)

    # Update the coordinates of the amino acids in the sequence object
    lattice.sequence.aa_coord_update(chain_index - 1, pos_L[0], pos_L[1])
    lattice.sequence.aa_coord_update(chain_index - 2, pos_C[0], pos_C[1])

def propagate_pull(
    lattice, 
    loop_index: int, 
    empty_pos: List[Tuple[int, int]]
) -> None:
    """
    Propagate the pull move through the chain, updating positions as necessary.

    Args:
        lattice: The lattice object containing sequence and grid information.
        loop_index (int): The index in the chain from which to start propagating.
        empty_pos (List[Tuple[int, int]]): List of empty positions available for moves.

    Returns:
        None
    """
    while loop_index >= 2:
        res_loop_iminus_1 = lattice.sequence.aa_coord[loop_index - 1]  # Get res i-1
        res_loop_iminus_2 = lattice.sequence.aa_coord[loop_index - 2]  # Get res i-2
        is_neighbor = are_topological_neighbors(res_loop_iminus_1, res_loop_iminus_2)

        if not is_neighbor:
            temp_pos = (res_loop_iminus_2["x"], res_loop_iminus_2["y"])
            free_position(lattice, res_loop_iminus_2)
            new_position = empty_pos.pop(0)
            move(lattice, lattice.sequence.aa_coord[loop_index - 2], new_position)
            lattice.sequence.aa_coord_update(loop_index - 2, new_position[0], new_position[1])
            empty_pos.append(temp_pos)
        if is_neighbor:
            break
        loop_index -= 1
