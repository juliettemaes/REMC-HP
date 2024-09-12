
from typing import Dict, List, Tuple
import random as rd

def is_end_move_possible(lattice, chain_index: int) -> int:
    if chain_index == 1: # the amino acid is at the beginning of the chain
        return 2
    elif chain_index == lattice.sequence.length: # the amino acid is at the end of the chain
        return lattice.sequence.length - 1
    else:
        raise ValueError("Impossible move. The amino acid is not at the end of the chain")

def execute_end_move(lattice, possible_pos: List[Tuple[int,int]], chain_index: int, old_x: int, old_y: int) -> None:
    new_pos = rd.choice(possible_pos)
    x_new, y_new = new_pos[0], new_pos[1]
    lattice.lattice[old_x, old_y] = 0
    lattice.lattice[x_new, y_new] = chain_index
    lattice.sequence.aa_coord_update(chain_index - 1, x_new, y_new) # update the coordinates of the amino acid in the sequence object