from copy import deepcopy
import numpy as np
import random as rd
from variables import GRID_SIZE_FACTOR
import sequence as seq

# Classes:
class Lattice():
    def __init__(self, sequence : seq.Sequence):
        self.sequence = sequence # str of the HP model sequence
        self.size = self.sequence.length * GRID_SIZE_FACTOR # size of the lattice
        self.lattice = np.zeros((self.size, self.size), dtype = int) # creation of an empty lattice
        self._initialize_random() # initialize the lattice with a random conformation
        # self.lattice_hp = None # lattice with the HP sequence
        self.lattice_initial = deepcopy(self.lattice) # save the initial lattice as an attribute
        self.energy = None # energy of the system

        # add conditions to say that it is rejected if the sequence is not valid (HP model)

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
                    possible_pos = self._find_empty_neighbors(i, j)

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


    def calculate_energy(self): ### add the check as separate functions
        # calculate the energy of the conformation
        energy = 0
        for i in range(self.size - 1):
            for j in range(self.size - 1):
                if (self.lattice[i,j] != 0) and (self.lattice[i, j+1] != 0): # checks if the amino acid is present
                    if (abs(self.lattice[i,j] - self.lattice[i, j+1]) != 1): # checks if the amino acids are connected neighbors
                        if self.sequence.hp_sequence[self.lattice[i,j]-1] == "H" and self.sequence.hp_sequence[self.lattice[i,j+1]-1] == "H": # checks if the amino acids are hydrophobic
                            energy -= 1
                if (self.lattice[i,j] != 0) and (self.lattice[i+1, j] != 0):
                    if abs(self.lattice[i,j] - self.lattice[i+1, j]) != 1:
                        if self.sequence.hp_sequence[self.lattice[i,j]-1] == "H" and self.sequence.hp_sequence[self.lattice[i+1,j]-1] == "H":
                            energy -= 1
        # handle the bottom and right edges of the lattice
        j = self.size - 1
        for i in range(self.size - 1):
            if self.lattice[i,j] != 0 and self.lattice[i+1, j] != 0:
                if abs(self.lattice[i,j] - self.lattice[i+1, j]) != 1:
                    if self.sequence.hp_sequence[self.lattice[i,j]-1] == "H" and self.sequence.hp_sequence[self.lattice[i+1,j]-1] == "H":
                        energy -= 1
        i = self.size - 1
        for j in range(self.size - 1):
            if self.lattice[i,j] != 0 and self.lattice[i, j+1] != 0:
                if abs(self.lattice[i,j] - self.lattice[i, j+1]) != 1:
                    if self.sequence.hp_sequence[self.lattice[i,j]-1] == "H" and self.sequence.hp_sequence[self.lattice[i,j+1]-1] == "H":
                        energy -= 1
        self.energy = energy

        

    def _find_empty_neighbors(self, i, j):
        # find the possible positions for the next amino acid
        possible_pos = []
        if self.lattice[i-1,j] == 0:
            possible_pos.append((i-1, j))
        if self.lattice[i+1,j] == 0:
            possible_pos.append((i+1, j))
        if self.lattice[i,j-1] == 0:
            possible_pos.append((i, j-1))
        if self.lattice[i,j+1] == 0:
            possible_pos.append((i, j+1))
        return possible_pos


    def end_move(self, chain_index):
        if chain_index == 1: # the amino acid is at the beginning of the chain
            res_ref = 2
        elif chain_index == self.sequence.length: # the amino acid is at the end of the chain
            res_ref = self.sequence.length - 1
        else:
            raise ValueError("Impossible move. The amino acid is not at the end of the chain")

        x, y = self.sequence.aa_coord[chain_index - 1]["x"], self.sequence.aa_coord[chain_index - 1]["y"]
        possible_pos = self._find_empty_neighbors(self.sequence.aa_coord[res_ref - 1]["x"], self.sequence.aa_coord[res_ref - 1]["y"])
        print(possible_pos)

        if len(possible_pos) != 0: # if there are possible new position(s)
            new_pos = rd.choice(possible_pos)
            x_new, y_new = new_pos[0], new_pos[1]
            self.lattice[x, y] = 0
            self.lattice[x_new, y_new] = chain_index
            self.sequence.aa_coord_update(chain_index - 1, x_new, y_new) # update the coordinates of the amino acid in the sequence object
            

    def corner_move(self, chain_index):
        if chain_index not in (1, self.sequence.length):
            res = self.sequence.aa_coord[chain_index - 1]
            res_prev = self.sequence.aa_coord[chain_index - 2]
            res_next = self.sequence.aa_coord[chain_index]
            
            if self._are_corner(res_prev, res_next):  # check if 3 connected neighbors are corner
                print(f"residue {chain_index} is a corner")
                x_new, y_new = self._find_potential_corner(chain_index)
                print(f"{x_new}, {y_new} is a potential corner available for residue {chain_index}")

                if self._check_occupancy(x_new, y_new):
                    print(f"position ({x_new}, {y_new}) is free")
                    self.lattice[x_new, y_new] = chain_index
                    self.lattice[res["x"], res["y"]] = 0
                    self.sequence.aa_coord_update(chain_index - 1, x_new, y_new) # update the coordinates of the amino acid in the sequence object
                    print(f"Corner move for amino acid {chain_index} at position ({res['x']}, {res['y']}) to position ({x_new}, {y_new})")
        else:
            raise ValueError("No Corner move possible. The amino acid is at the end of the chain")


    def _are_adjacent(self, res1, res2):
        return abs(res1['x'] - res2['x']) + abs(res1['y'] - res2['y']) == 1


    def _are_corner(self, res1, res2):
        return ((abs(res1['x'] - res2['x']) + abs(res1['y'] - res2['y']) == 2) and (res1['x'] != res2['x']) and (res1['y'] != res2['y']))

    def _check_occupancy(self, i, j):
        # returns boolean, check if the positions are free by an amino acids
        return (self.lattice[i,j] == 0)

    def _find_potential_corner(self, chain_index):
        if self.sequence.aa_coord[chain_index - 1]["x"] == self.sequence.aa_coord[chain_index - 2]["x"]:
            y_new = self.sequence.aa_coord[chain_index - 2]["y"]
            x_new = self.sequence.aa_coord[chain_index ]["x"]
        else:
            x_new = self.sequence.aa_coord[chain_index - 2]["x"]
            y_new = self.sequence.aa_coord[chain_index]["y"]
        return x_new, y_new
    

    def cks_move(self, chain_index):
        if chain_index not in (1, self.sequence.length): # can't perform the move if the amino acid is at the end of the chain
            res = self.sequence.aa_coord[chain_index - 1]
            res_prev = self.sequence.aa_coord[chain_index - 2]
            res_next = self.sequence.aa_coord[chain_index]
            res_next2 = self.sequence.aa_coord[chain_index + 1]
            # check if 4 connected neighbors for a U
            if self._are_U(res_prev, res, res_next, res_next2):
                print(f"residue {chain_index} is a U")
                x_new1, y_new1, x_new2, y_new2 = self._find_potential_U(res_prev, res, res_next, res_next2)
                if self._check_occupancy(x_new1, y_new1) and self._check_occupancy(x_new2, y_new2):
                    # move the amino acids to the new positions
                    self.lattice[x_new1, y_new1] = chain_index
                    self.lattice[x_new2, y_new2] = chain_index + 1
                    # free the positions of the amino acids
                    self.lattice[res["x"], res["y"]] = 0
                    self.lattice[res_next["x"], res_next["y"]] = 0
                    # update the coordinates of the amino acid in the sequence object
                    self.sequence.aa_coord_update(chain_index - 1, x_new1, y_new1)
                    self.sequence.aa_coord_update(chain_index, x_new2, y_new2)
                    print(f"CKS move for amino acid {chain_index} at position ({res['x']}, {res['y']}) to position ({x_new1}, {y_new1})")

            else: # check for alternative U position
                res_prev2 = self.sequence.aa_coord[chain_index - 3]
                res_prev = self.sequence.aa_coord[chain_index - 2]
                res = self.sequence.aa_coord[chain_index - 1]
                res_next = self.sequence.aa_coord[chain_index]

                if self._are_U(res_prev2, res_prev, res, res_next):
                    print(f"residue {chain_index} is a U")
                    x_new1, y_new1, x_new2, y_new2 = self._find_potential_U(res_prev2, res_prev, res, res_next)
                    if self._check_occupancy(x_new1, y_new1) and self._check_occupancy(x_new2, y_new2):
                        # move the amino acids to the new positions
                        self.lattice[x_new1, y_new1] = chain_index - 1
                        self.lattice[x_new2, y_new2] = chain_index
                        # free the positions of the amino acids
                        self.lattice[res_prev["x"], res_prev["y"]] = 0
                        self.lattice[res["x"], res["y"]] = 0
                        # update the coordinates of the amino acid in the sequence object
                        self.sequence.aa_coord_update(chain_index - 2, x_new1, y_new1)
                        self.sequence.aa_coord_update(chain_index - 1, x_new2, y_new2)
                        print(f"CKS move for amino acid {chain_index} at position ({res['x']}, {res['y']}) to position ({x_new2}, {y_new2})")




    def _are_U(self, res1, res2, res3, res4):
        cdt1 = abs(res1['x'] - res2['x']) + abs(res1['y'] - res2['y']) == 1
        cdt2 = abs(res2['x'] - res3['x']) + abs(res2['y'] - res3['y']) == 1
        cdt3 = abs(res3['x'] - res4['x']) + abs(res3['y'] - res4['y']) == 1
        cdt4 = abs(res4['x'] - res1['x']) + abs(res4['y'] - res1['y']) == 1
        return (cdt1 and cdt2 and cdt3 and cdt4)

    def _find_potential_U(self, res1, res2, res3, res4):
        x_new1 = res1['x'] - self._relative_position(res2, res1)[0]
        y_new1 = res1['y'] - self._relative_position(res2, res1)[1]

        x_new2 = res4['x'] - self._relative_position(res3, res4)[0]
        y_new2 = res4['y'] - self._relative_position(res3, res4)[1]
        
        return x_new1, y_new1, x_new2, y_new2


    def _relative_position(self, res1, res2):
        # returns the relative position of res1 compared to res2
        dx = res1['x'] - res2['x']
        dy = res1['y'] - res2['y']
        return dx, dy


    def pull_move(self, chain_index):
        pass



































############## TO DEBUG ###########
    def lattice_hp(self):
        # Show the HP sequence on the lattice instead of the amino acid number
        print(type(self.lattice))
        self.lattice_hp = deepcopy(self.lattice)
        print(self.lattice_hp)
        self.lattice_hp.astype(str)
        print(self.lattice_hp)
        for i in range(self.size):
            for j in range(self.size):
                if self.lattice[i,j] != 0:
                    print(self.lattice[i,j])
                    self.lattice[i,j] = str(self.sequence.hp_sequence[self.lattice[i,j]-1])


######## TO DO ########
    def translate_to_center(self, x, y):
        # translate the lattice by x and y to recenter the conformation
        pass

    def rmsd(self):
        # calculate the root mean square deviation between the initial and final conformation
        pass






    def _check_connectivity(self, i, j):
        # check if the amino acid is connected to the chain
        return (abs(self.lattice[i,j] - self.lattice[i, j+1]) != 1)


    def _check_hydrophobicity(self, i, j):
        # check if the amino acid is hydrophobic
        pass

    # def _check_occupancy(self, i, j):
    #     # returns boolean, check if the positions are occupied by an amino acids
    #     return (self.lattice[i,j] != 0) and (self.lattice[i, j+1] != 0)
   
    #### TO KEEP OR NOT ? ####
    def initialize_extended(self):
        # place the HP sequence in the lattice with an extended conformation at the center of the lattice
        start_pos_i = self.size // 2
        start_pos_j = self.size // 2 - self.sequence.length // 2
        for j in range(start_pos_j, start_pos_j + self.sequence.length):
            self.lattice[start_pos_i,j] = self.lattice[start_pos_i,j-1] + 1


    def _are_topological_neighbors(self, res1, res2):
        # check if two amino acids are topological neighbors
        return abs(res1['x'] - res2['x']) + abs(res1['y'] - res2['y']) == 1

    def _are_not_connected_neighbors(self, res1, res2):
        # check if two amino acids are connected neighbors
        return abs(res1['index'] - res2['index']) > 1

    def _are_hydrophobic(self, res1, res2):
        # check if two amino acids are hydrophobic
        return res1['type'] == 'H' and res2['type'] == 'H'

    def calculate_energy2(self):
        energy = 0
        n = self.sequence.length
        for j in range(n-1):
            for k in range(j+1, n):
                if self._are_hydrophobic(self.sequence.aa_coord[j], self.sequence.aa_coord[k]):
                    if self._are_topological_neighbors(self.sequence.aa_coord[j], self.sequence.aa_coord[k]):
                        if self._are_not_connected_neighbors(self.sequence.aa_coord[j], self.sequence.aa_coord[k]):
                            energy -= 1
        return energy


    # # back up


    # def _initialize_random(self):
    # # place the HP sequence in the lattice with a random conformation
    # start_pos_i = self.size // 2
    # start_pos_j = self.size // 2
    # i = start_pos_i
    # j = start_pos_j
    # self.lattice[i,j] = 1
    # self.sequence.aa_coord_update(0, i, j) # update the coordinates of the amino acid in the sequence object
    # for aa in range(1,self.sequence.length):
    #     possible_pos = self._find_empty_neighbors(i, j)
    #     if len(possible_pos) == 0:
    #         raise ValueError("No possible position for the next amino acid, restart the initialization")
    #         ### add a condition to restart the initialization, test with try and except!!
    #     else:
    #         pos = rd.choice(possible_pos)
    #         i = pos[0]
    #         j = pos[1]
    #         self.lattice[i,j] = aa + 1 # code the amino acid nummber of the sequence in the lattice
    #         # update position for the next amino acid
    #         self.sequence.aa_coord_update(aa, i, j) # update the coordinates of the amino acid in the sequence object
    #         i = pos[0]
    #         j = pos[1]