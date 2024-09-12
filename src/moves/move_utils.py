"""
General utility functions for handling movements in a lattice-based system.
"""

from typing import Dict, List, Tuple
import numpy as np

def is_hydrophobic(res: Dict[str, str]) -> bool:
    """
    Check if a given amino acid is hydrophobic.

    Args:
        res (Dict[str, str]): A dictionary containing information about an amino acid, 
                              with 'type' being a key indicating its hydrophobicity.

    Returns:
        bool: True if the amino acid is hydrophobic (indicated by 'H'), False otherwise.
    """
    return res['type'] == 'H'


def are_topological_neighbors(res1: Dict[str, int], res2: Dict[str, int]) -> bool:
    """
    Check if two amino acids are topological neighbors, meaning they are adjacent 
    on a grid (either horizontally or vertically).

    Args:
        res1 (Dict[str, int]): A dictionary representing the first amino acid.
        res2 (Dict[str, int]): A dictionary representing the second amino acid.

    Returns:
        bool: True if the amino acids are topological neighbors, False otherwise.
    """
    return abs(res1['x'] - res2['x']) + abs(res1['y'] - res2['y']) == 1


def are_not_connected_neighbors(res1: Dict[str, int], res2: Dict[str, int]) -> bool:
    """
    Check if two amino acids are directly connected or not.

    Args:
        res1 (Dict[str, int]): A dictionary representing the first amino acid.
        res2 (Dict[str, int]): A dictionary representing the second amino acid.

    Returns:
        bool: True if the amino acids are connected, False otherwise.
    """
    return abs(res1['index'] - res2['index']) > 1


def are_corner(res1: Dict[str, int], res2: Dict[str, int]) -> bool:
    """
    Check if two amino acids are positioned at a potential corner.

    Args:
        res1 (Dict[str, int]): A dictionary representing the first amino acid.
        res2 (Dict[str, int]): A dictionary representing the second amino acid.

    Returns:
        bool: True if the amino acids are in a corner, False otherwise.
    """
    return (abs(res1['x'] - res2['x']) + abs(res1['y'] - res2['y']) == 2) and (res1['x'] != res2['x']) and (res1['y'] != res2['y'])


def check_occupancy(lattice, i: int, j: int) -> bool:
    """
    Check if a given position in the lattice is occupied by an amino acid.

    Args:
        lattice: The lattice object containing the grid structure.
        i (int): The x-coordinate of the position.
        j (int): The y-coordinate of the position.

    Returns:
        bool: True if the position is unoccupied, False otherwise.
    """
    return lattice.lattice[i, j] == 0


def find_empty_diagonal(lattice, i: int, j: int) -> List[Tuple[int, int]]:
    """
    Find the diagonally adjacent empty positions for a given location in the lattice.

    Args:
        lattice: The lattice object containing the grid structure.
        i (int): The x-coordinate of the current position.
        j (int): The y-coordinate of the current position.

    Returns:
        List[Tuple[int, int]]: A list of tuples representing the coordinates of the empty diagonal positions.
    """
    possible_pos = []
    if lattice.lattice[i-1, j-1] == 0:
        possible_pos.append((i-1, j-1))
    if lattice.lattice[i+1, j+1] == 0:
        possible_pos.append((i+1, j+1))
    if lattice.lattice[i-1, j+1] == 0:
        possible_pos.append((i-1, j+1))
    if lattice.lattice[i+1, j-1] == 0:
        possible_pos.append((i+1, j-1))
    return possible_pos


def find_empty_neighbors(lattice, i: int, j: int) -> List[Tuple[int, int]]:
    """
    Find the empty neighboring positions (horizontally or vertically) for a given location in the lattice.

    Args:
        lattice: The lattice object containing the grid structure.
        i (int): The x-coordinate of the current position.
        j (int): The y-coordinate of the current position.

    Returns:
        List[Tuple[int, int]]: A list of tuples representing the coordinates of the empty neighboring positions.
    """
    possible_pos = []
    if lattice.lattice[i-1, j] == 0:
        possible_pos.append((i-1, j))
    if lattice.lattice[i+1, j] == 0:
        possible_pos.append((i+1, j))
    if lattice.lattice[i, j-1] == 0:
        possible_pos.append((i, j-1))
    if lattice.lattice[i, j+1] == 0:
        possible_pos.append((i, j+1))
    return possible_pos


def relative_position(res1: Dict[str, int], res2: Dict[str, int]) -> Tuple[int, int]:
    """
    Calculate the relative position of one amino acid compared to another.

    Args:
        res1 (Dict[str, int]): A dictionary representing the first amino acid.
        res2 (Dict[str, int]): A dictionary representing the second amino acid.

    Returns:
        Tuple[int, int]: A tuple representing the relative position as (dx, dy).
    """
    dx = res1['x'] - res2['x']
    dy = res1['y'] - res2['y']
    return dx, dy


def move(lattice, res_to_be_moved: Dict[str, int], place_to_move: Tuple[int, int]) -> None:
    """
    Move an amino acid to a new position in the lattice. 
    It only updates the lattice. We will update the amino acid's coordinates separately.

    Args:
        lattice: The lattice object containing the grid structure.
        res_to_be_moved (Dict[str, int]): A dictionary representing the amino acid to be moved.
        place_to_move (Tuple[int, int]): A tuple representing the coordinates of the new position.

    Returns:
        None
    """
    lattice.lattice[place_to_move[0], place_to_move[1]] = res_to_be_moved["chain_index"]


def free_position(lattice, res_to_free: Dict[str, int]) -> None:
    """
    Free the position occupied by a given amino acid in the lattice.

    Args:
        lattice: The lattice object containing the grid structure.
        res_to_free (Dict[str, int]): A dictionary representing the amino acid to be freed.

    Returns:
        None
    """
    lattice.lattice[res_to_free["x"], res_to_free["y"]] = 0


def translate_chain(lattice, x, y):
    # translate the lattice by x and y to recenter the conformation
    lattice.lattice = np.roll(lattice.lattice, x, axis = 0)
    lattice.lattice = np.roll(lattice.lattice, y, axis = 1)
    # update the coordinates of the amino acids in the sequence object
    for i in range(lattice.sequence.length):
        lattice.sequence.aa_coord_update(i, lattice.sequence.aa_coord[i]["x"] + x, lattice.sequence.aa_coord[i]["y"] + y)


def need_translation_x(lattice):
    list_x = [temp_dict["x"] for temp_dict in lattice.sequence.aa_coord]
    max_x = max(list_x)
    min_x = min(list_x)
    center = lattice.size // 2
    middle = (min_x + max_x) // 2
    dx = center - middle
    if max_x >= lattice.size - 2 or min_x < 1:
        return True, dx
    else:
        return False, 0

def need_translation_y(lattice):
    list_y = [temp_dict["y"] for temp_dict in lattice.sequence.aa_coord]
    max_y = max(list_y)
    min_y = min(list_y)
    center = lattice.size // 2
    middle = (min_y + max_y) // 2
    dy = center - middle
    if max_y >= lattice.size - 2 or min_y < 1:
        return True, dy
    else:
        return False, 0