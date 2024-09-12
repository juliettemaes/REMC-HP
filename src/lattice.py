from copy import deepcopy
import numpy as np
import random as rd
from variables import GRID_SIZE_FACTOR
from moves.move_utils import *
from moves.end_move_utils import *
from moves.corner_move_utils import *
from moves.cks_utils import *
from moves.pull_move_utils import *
import sequence as seq
from typing import Tuple, List

class Lattice:
    """
    Class representing the lattice for the protein, its sequence, along with methods to manipulate it.

    Attributes:
        sequence (seq.Sequence): The sequence of amino acids.
        size (int): The size of the lattice.
        lattice (np.ndarray): The 2D array representing the lattice.
        lattice_initial (np.ndarray): The initial configuration of the lattice.
        energy (Optional[float]): The energy of the lattice.
    """

    def __init__(self, sequence: seq.Sequence):
        """
        Initialize the lattice with a given sequence.

        Args:
            sequence (seq.Sequence): The sequence of amino acids.
        """
        self.sequence = sequence
        self.size = self.sequence.length * GRID_SIZE_FACTOR
        self.lattice = np.zeros((self.size, self.size), dtype=int)
        self._initialize_random()
        self.lattice_initial = deepcopy(self.lattice)
        self.energy = None

    def __str__(self) -> str:
        """
        Return the string representation of the lattice.

        Returns:
            str: The string representation of the lattice.
        """
        return str(self.lattice)

    def _initialize_random(self) -> None:
        """
        Initialize the lattice with a random valid conformation.
        
        This method places the amino acids on center of the lattice in a valid conformation.
        """
        while True:
            try:
                self.lattice.fill(0)
                start_pos_i = self.size // 2
                start_pos_j = self.size // 2
                i, j = start_pos_i, start_pos_j

                self.lattice[i, j] = 1
                self.sequence.aa_coord_update(0, i, j)

                for aa_index in range(1, self.sequence.length):
                    possible_pos = find_empty_neighbors(self, i, j)

                    if len(possible_pos) == 0:
                        raise ValueError("No possible position for the next amino acid, restarting initialization")

                    pos = rd.choice(possible_pos)
                    i, j = pos[0], pos[1]
                    self.lattice[i, j] = aa_index + 1
                    self.sequence.aa_coord_update(aa_index, i, j)

                break

            except ValueError:
                continue

    def calculate_energy(self) -> float:
        """
        Calculate the energy of the lattice based on the hydrophobic interactions.

        Returns:
            float: The calculated energy of the lattice.
        """
        energy = 0
        n = self.sequence.length
        for j in range(n - 1):
            if is_hydrophobic(self.sequence.aa_coord[j]):
                for k in range(j + 1, n):
                    if is_hydrophobic(self.sequence.aa_coord[k]):
                        if are_topological_neighbors(self.sequence.aa_coord[j], self.sequence.aa_coord[k]):
                            if are_not_connected_neighbors(self.sequence.aa_coord[j], self.sequence.aa_coord[k]):
                                energy -= 1
        return energy

    def end_move(self, chain_index: int) -> None:
        """
        Attempt to perform an end move for the amino acid at the specified chain index.
        The move is only performed if it is possible.
        If the move is not possible, the method does nothing.

        Args:
            chain_index (int): The chain index of the amino acid to move.
        """
        res_ref = is_end_move_possible(self, chain_index)
        x, y = self.sequence.aa_coord[chain_index - 1]["x"], self.sequence.aa_coord[chain_index - 1]["y"]
        possible_pos = find_empty_neighbors(self, self.sequence.aa_coord[res_ref - 1]["x"], self.sequence.aa_coord[res_ref - 1]["y"])

        if possible_pos:
            execute_end_move(self, possible_pos, chain_index, x, y)

    def corner_move(self, chain_index: int) -> None:
        """
        Attempt to perform a corner move for the amino acid at the specified chain index.
        The move is only performed if it is possible.
        If the move is not possible, the method does nothing.

        Args:
            chain_index (int): The chain index of the amino acid to move.

        Raises:
            ValueError: If the amino acid is at the end of the chain, where a corner move is not possible.
        """
        if chain_index not in (1, self.sequence.length):
            res_i, res_iminus1, res_iplus1 = get_corner_aa_to_check(self, chain_index)
            is_corner = are_corner(res_iminus1, res_iplus1)
            
            if is_corner:
                x_new, y_new = find_potential_corner(self, chain_index)
                if check_occupancy(self, x_new, y_new):
                    execute_corner_move(self, chain_index, res_i, x_new, y_new)
        else:
            raise ValueError("No Corner move possible. The amino acid is at the end of the chain")

    def cks_move(self, chain_index: int) -> None:
        """
        Attempt to perform a Crankshaft move for the amino acid at the specified chain index.
        The move is only performed if it is possible.
        If the move is not possible, the method does nothing.

        Args:
            chain_index (int): The chain index of the amino acid to move.

        Raises:
            ValueError: If the amino acid is at the end of the chain, where a CKS move is not possible.
        """
        if chain_index not in (1, self.sequence.length):
            res_i, res_iminus1, res_iplus1, res_next2 = get_cks_aa_to_check(self, chain_index)
            is_u_valid, x_new1, y_new1, x_new2, y_new2 = validate_U(self, res_iminus1, res_i, res_iplus1, res_next2, chain_index)

            if is_u_valid:
                execute_u_move(self, chain_index, res_i, res_iplus1, x_new1, y_new1, x_new2, y_new2)
            elif chain_index != 2:
                res_prev2, res_iminus1, res_i, res_iplus1 = get_alternative_positions(self, chain_index)
                is_alt_u_valid, x_new1, y_new1, x_new2, y_new2 = validate_U(self, res_prev2, res_iminus1, res_i, res_iplus1, chain_index)
                if is_alt_u_valid:
                    execute_alternative_u_move(self, chain_index, res_i, res_iminus1, x_new1, y_new1, x_new2, y_new2)
        else:
            raise ValueError("No CKS move possible. The amino acid is at the end of the chain")

    def pull_move(self, chain_index: int) -> None:
        """
        Attempt to perform a pull move for the amino acid at the specified chain index.
        The move is only performed if it is possible.
        If the move is not possible, the method does nothing.

        Args:
            chain_index (int): The chain index of the amino acid to move.
        """
        res_iminus2, res_iminus1, res_i, res_iplus1 = get_pull_aa_to_check(self, chain_index)
        pos_L, found_L = find_L_positions(self, res_i, res_iplus1)

        if found_L:
            pos_C, status_C, found_C = find_C_positions(self, pos_L, res_i, res_iminus1)

            if found_C:
                if status_C == "C is res_iminus1":
                    self.corner_move(chain_index)
                elif status_C == "empty":
                    empty_pos = [(res_i["x"], res_i["y"]), (res_iminus1["x"], res_iminus1["y"])]
                    execute_pull_move(self, chain_index, res_i, res_iminus1, pos_L, pos_C)
                    if res_i["chain_index"] >= 3:
                        loop_index = chain_index - 1
                        propagate_pull(self, loop_index, empty_pos)

    def translation_check_and_apply(self) -> None:
        """
        Check if the protein has reached the edge of the lattice and translate the chain to the center if needed.
        """
        check_x, dx = need_translation_x(self)
        check_y, dy = need_translation_y(self)
        if check_x or check_y:
            translate_chain(self, dx, dy)
