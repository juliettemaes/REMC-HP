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

# Classes:
class Lattice():
    def __init__(self, sequence : seq.Sequence):
        self.sequence = sequence # str of the HP model sequence
        self.size = self.sequence.length * GRID_SIZE_FACTOR # size of the lattice
        self.lattice = np.zeros((self.size, self.size), dtype = int) # creation of an empty lattice
        self._initialize_random() # initialize the lattice with a random conformation
        self.lattice_initial = deepcopy(self.lattice) # save the initial lattice as an attribute
        self.energy = None # energy of the system


    def __str__(self):
        return str(self.lattice)

    def _initialize_random(self):
        while True:  # Loop until a valid conformation is found
            try:
                # Clear the lattice and sequence coordinates before restarting
                self.lattice.fill(0)  # Clear the lattice 
                
                # Start placing the HP sequence in the lattice with a random conformation
                start_pos_i = self.size // 2
                start_pos_j = self.size // 2
                i = start_pos_i
                j = start_pos_j

                self.lattice[i, j] = 1
                self.sequence.aa_coord_update(0, i, j)  # Update the coordinates of the amino acid in the sequence object

                for aa_index in range(1, self.sequence.length):
                    possible_pos = find_empty_neighbors(self, i, j)

                    # If no possible positions are available, raise an exception to restart the initialization
                    if len(possible_pos) == 0:
                        raise ValueError("No possible position for the next amino acid, restarting initialization")

                    # Select a random position from the available neighbors
                    pos = rd.choice(possible_pos)
                    i = pos[0]
                    j = pos[1]

                    # Place the amino acid on the lattice
                    self.lattice[i, j] = aa_index + 1  # Code the amino acid number of the sequence in the lattice

                    # Update position for the next amino acid
                    self.sequence.aa_coord_update(aa_index, i, j)  # Update the coordinates of the amino acid in the sequence object

                # If the loop completes without errors, break out of the while loop
                break

            except ValueError:
                # If an exception is raised, restart the initialization (continue the loop)
                continue


    def calculate_energy(self):
        energy = 0
        n = self.sequence.length
        for j in range(n-1):
            if is_hydrophobic(self.sequence.aa_coord[j]):
                for k in range(j+1, n):
                    if is_hydrophobic(self.sequence.aa_coord[k]):
                        if are_topological_neighbors(self.sequence.aa_coord[j], self.sequence.aa_coord[k]):
                            if are_not_connected_neighbors(self.sequence.aa_coord[j], self.sequence.aa_coord[k]):
                                energy -= 1
        return energy

    
    def end_move(self, chain_index):
        #Check if the move is possible znd return the reference position
        res_ref = is_end_move_possible(self, chain_index)

        # Find the possible positions
        x, y = self.sequence.aa_coord[chain_index - 1]["x"], self.sequence.aa_coord[chain_index - 1]["y"]
        possible_pos = find_empty_neighbors(self, self.sequence.aa_coord[res_ref - 1]["x"], self.sequence.aa_coord[res_ref - 1]["y"])

        # Swap if possible
        if len(possible_pos) != 0: # if there are possible new position(s)
            execute_end_move(self, possible_pos, chain_index, x, y)
            

    def corner_move(self, chain_index):
        # Check if we are not dealing with the first or last amino acid that can't be corner moved.
        if chain_index not in (1, self.sequence.length):
            
            # Get the amino acids & check if they are corner
            res_i, res_iminus1, res_iplus1 = get_corner_aa_to_check(self, chain_index)
            is_corner = are_corner(res_iminus1, res_iplus1)
            
            if is_corner:  # If it is a corner, find the potential corner position
                x_new, y_new = find_potential_corner(self, chain_index)

                if check_occupancy(self, x_new, y_new): # If the corner is free, execute the move
                    execute_corner_move(self, chain_index, res_i, x_new, y_new)
        else:
            raise ValueError("No Corner move possible. The amino acid is at the end of the chain")

    
    def cks_move(self, chain_index):
        # can't perform the move if the amino acid is at the end of the chain
        if chain_index not in (1, self.sequence.length): 

            res_i,res_iminus1,res_iplus1,res_next2 = get_cks_aa_to_check(self, chain_index)

            # check if 4 connected neighbors for a U & if the move is possible
            is_u_valid, x_new1, y_new1, x_new2, y_new2 = validate_U(self, res_iminus1, res_i, res_iplus1, res_next2, chain_index)

            if is_u_valid:
                execute_u_move(self, chain_index, res_i, res_iplus1, x_new1, y_new1, x_new2, y_new2)

            # check for alternative cases
            elif chain_index !=2 :
                res_prev2, res_iminus1, res_i, res_iplus1 = get_alternative_positions(self, chain_index)
                is_alt_u_valid, x_new1, y_new1, x_new2, y_new2 = validate_U(self, res_prev2, res_iminus1, res_i, res_iplus1, chain_index)

                if is_alt_u_valid:
                    execute_alternative_u_move(self, chain_index, res_i, res_iminus1, x_new1, y_new1, x_new2, y_new2)
        else:
            raise ValueError("No CKS move possible. The amino acid is at the end of the chain")


    def pull_move(self, chain_index):
        # Get the amino acids to check
        res_iminus2, res_iminus1, res_i, res_iplus1 = get_pull_aa_to_check(self, chain_index)
        # Find the potential positions of L
        pos_L,found_L = find_L_positions(self, res_i, res_iplus1)

        if found_L:
            # If L is found, find the potential positions of C
            pos_C, status_C, found_C = find_C_positions(self, pos_L, res_i, res_iminus1)

            if found_C is True:
                if status_C == "C is res_iminus1":
                    # perform corner move
                    self.corner_move(chain_index)
                elif status_C == "empty":
                    # perform pull move & save positions of i and i-1
                    empty_pos = [(res_i["x"], res_i["y"]), (res_iminus1["x"], res_iminus1["y"])]
                    # move i and i+1 to L and C and free the prevous positions
                    execute_pull_move(self, chain_index, res_i, res_iminus1, pos_L, pos_C)
                    # propagate pull moves
                    if res_i["chain_index"] >= 3:
                        loop_index = chain_index -1 #init index to i-1
                        propagate_pull(self, loop_index, empty_pos)


    def translation_check_and_apply(self):
        check_x, dx = need_translation_x(self)
        check_y, dy = need_translation_y(self)
        if check_x or check_y:
            translate_chain(self, dx, dy)