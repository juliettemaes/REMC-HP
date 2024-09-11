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
        # print(possible_pos)

        if len(possible_pos) != 0: # if there are possible new position(s)
            new_pos = rd.choice(possible_pos)
            x_new, y_new = new_pos[0], new_pos[1]
            self.lattice[x, y] = 0
            self.lattice[x_new, y_new] = chain_index
            self.sequence.aa_coord_update(chain_index - 1, x_new, y_new) # update the coordinates of the amino acid in the sequence object
            

    def corner_move(self, chain_index):
        if chain_index not in (1, self.sequence.length):
            res_i = self.sequence.aa_coord[chain_index - 1]
            res_iminus1 = self.sequence.aa_coord[chain_index - 2]
            res_iplus1 = self.sequence.aa_coord[chain_index]
            
            if self._are_corner(res_iminus1, res_iplus1):  # check if 3 connected neighbors are corner
                # print(f"residue {chain_index} is a corner")
                x_new, y_new = self._find_potential_corner(chain_index)
                # print(f"{x_new}, {y_new} is a potential corner available for residue {chain_index}")

                if self._check_occupancy(x_new, y_new):
                    # print(f"position ({x_new}, {y_new}) is free")
                    self.lattice[x_new, y_new] = chain_index
                    self.lattice[res_i["x"], res_i["y"]] = 0
                    self.sequence.aa_coord_update(chain_index - 1, x_new, y_new) # update the coordinates of the amino acid in the sequence object
                    # print(f"Corner move for amino acid {chain_index} at position ({res_i['x']}, {res_i['y']}) to position ({x_new}, {y_new})")
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
    
    #### TO IMPROVE ####
    def cks_move(self, chain_index):
        if chain_index not in (1, self.sequence.length): # can't perform the move if the amino acid is at the end of the chain
            if chain_index < self.sequence.length - 1:
                res_i = self.sequence.aa_coord[chain_index - 1]
                res_iminus1 = self.sequence.aa_coord[chain_index - 2]
                res_iplus1 = self.sequence.aa_coord[chain_index]
                res_next2 = self.sequence.aa_coord[chain_index + 1]
            elif chain_index == self.sequence.length - 1:
                res_i = self.sequence.aa_coord[chain_index - 2]
                res_iminus1 = self.sequence.aa_coord[chain_index - 3]
                res_iplus1 = self.sequence.aa_coord[chain_index - 1]
                res_next2 = self.sequence.aa_coord[chain_index]

            # check if 4 connected neighbors for a U
            if self._are_U(res_iminus1, res_i, res_iplus1, res_next2) and chain_index < self.sequence.length - 1:
                # print(f"residue {chain_index} is a U")
                x_new1, y_new1, x_new2, y_new2 = self._find_potential_U(res_iminus1, res_i, res_iplus1, res_next2)
                if self._check_occupancy(x_new1, y_new1) and self._check_occupancy(x_new2, y_new2):
                    # move the amino acids to the new positions
                    self.lattice[x_new1, y_new1] = chain_index
                    self.lattice[x_new2, y_new2] = chain_index + 1
                    # free the positions of the amino acids
                    self.lattice[res_i["x"], res_i["y"]] = 0
                    self.lattice[res_iplus1["x"], res_iplus1["y"]] = 0
                    # update the coordinates of the amino acid in the sequence object
                    self.sequence.aa_coord_update(chain_index - 1, x_new1, y_new1)
                    self.sequence.aa_coord_update(chain_index, x_new2, y_new2)
                    # print(f"CKS move for amino acid {chain_index} at position ({res_i['x']}, {res_i['y']}) to position ({x_new1}, {y_new1})")

            elif chain_index !=2 : # check for alternative U position
                res_prev2 = self.sequence.aa_coord[chain_index - 3]
                res_iminus1 = self.sequence.aa_coord[chain_index - 2]
                res_i = self.sequence.aa_coord[chain_index - 1]
                res_iplus1 = self.sequence.aa_coord[chain_index]

                if self._are_U(res_prev2, res_iminus1, res_i, res_iplus1):
                    # print(f"residue {chain_index} is a U")
                    x_new1, y_new1, x_new2, y_new2 = self._find_potential_U(res_prev2, res_iminus1, res_i, res_iplus1)
                    if self._check_occupancy(x_new1, y_new1) and self._check_occupancy(x_new2, y_new2):
                        # move the amino acids to the new positions
                        self.lattice[x_new1, y_new1] = chain_index - 1
                        self.lattice[x_new2, y_new2] = chain_index
                        # free the positions of the amino acids
                        self.lattice[res_iminus1["x"], res_iminus1["y"]] = 0
                        self.lattice[res_i["x"], res_i["y"]] = 0
                        # update the coordinates of the amino acid in the sequence object
                        self.sequence.aa_coord_update(chain_index - 2, x_new1, y_new1)
                        self.sequence.aa_coord_update(chain_index - 1, x_new2, y_new2)
                        # print(f"CKS move for amino acid {chain_index} at position ({res_i['x']}, {res_i['y']}) to position ({x_new2}, {y_new2})")




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
        # print("-------------------------")
        # print("STARTING A NEW PULL MOVE")
        # print("-------------------------")
        # print("INFORMATIONS - CURRENT CHAIN_INDEX : ", chain_index, "CURRENT RESIDU : ", self.sequence.aa_coord[chain_index - 1])
        # print("-------------------------")
        # print("INFORMATIONS - I = ", chain_index, "I+1 : ", chain_index+1)
        # print("-------------------------")
        # print("INFORMATIONS - CURRENT LATTICE BEFORE MOVE : ", self.lattice)

        res_i = self.sequence.aa_coord[chain_index - 1] # equals to CHAIN INDEX
        res_iplus1 = self.sequence.aa_coord[chain_index]
        res_iminus1 = self.sequence.aa_coord[chain_index - 2] # equals to CHAIN INDEX -1
        res_iminus2 = self.sequence.aa_coord[chain_index - 3] # equals to CHAIN INDEX -2
        pos_L,found_L = self._find_L_positions(res_i, res_iplus1)
        if found_L:
            pos_C, status_C, found_C = self._find_C_positions(pos_L, res_i, res_iminus1)
            if found_C is False:
                print("No possible move as no C res_i available.")
            else:
                if status_C == "C is res_iminus1":
                    # perform corner move
                    self.corner_move(chain_index)
                    # print("SUCCESSFULL CORNER MOVE FOR RESIDUE", chain_index)
                elif status_C == "empty":
                    # perform pull move
                    # save positions of i and i-1
                    empty_pos = [(res_i["x"], res_i["y"]), (res_iminus1["x"], res_iminus1["y"])]
                    # print("POSITIONS TO BE EMPTIED : ", empty_pos)
                    # move i and i+1 to L and C and free the prevous positions
                    self._move(res_i, pos_L)
                    self._move(res_iminus1, pos_C)

                    # print("new position for", chain_index, "is", pos_L)
                    # print("coordinates of", chain_index, "are", self.sequence.aa_coord[chain_index - 1])

                    self._free_position(res_i)
                    self._free_position(res_iminus1)
                    # update the coordinates of the amino acid in the sequence object
                    self.sequence.aa_coord_update(chain_index - 1, pos_L[0], pos_L[1])
                    self.sequence.aa_coord_update(chain_index - 2, pos_C[0], pos_C[1])

                   #  print("SUCCESSFULL PULL MOVE FOR RESIDUE FOR", res_i["chain_index"], "AND", res_iminus1["chain_index"])
                   #  print(self.lattice)
                   # print("Current index", res_i["chain_index"])

                    if res_i["chain_index"] >= 3:
                        # print("Iterating")
                        # pull the chain
                        
                        loop_index = chain_index -1 #init index to i-1

                        while loop_index >= 2 :
                            # print("loop index", loop_index)
                            res_loop_iminus_1 = self.sequence.aa_coord[loop_index-1] # get res i-1
                            res_loop_iminus_2= self.sequence.aa_coord[loop_index-2] # get res i-2
                            # print("loop minus1", res_loop_iminus_1)
                            # print("loop minus2", res_loop_iminus_2)
                            is_neighbor = self._are_topological_neighbors(res_loop_iminus_1, res_loop_iminus_2)
                            # print("is neighbor", is_neighbor)
                            if is_neighbor is False:
                                # print("MOVING I MINUS 2")
                                temp_pos = (res_loop_iminus_2["x"], res_loop_iminus_2["y"])
                                self._free_position(res_loop_iminus_2)
                                new_position = empty_pos.pop(0)
                                # print("position to take : ",new_position)
                                self._move(self.sequence.aa_coord[loop_index-2], new_position)
                                self.sequence.aa_coord_update(loop_index-2, new_position[0], new_position[1])
                                empty_pos.append(temp_pos)
                            if is_neighbor is True:
                                # print("LOOP BREAK")
                                break
                            loop_index -= 1

        else:
            print("No possible move as no L res_i available.")

    
    def _move(self, res_to_be_moved, place_to_move):
        self.lattice[place_to_move[0], place_to_move[1]] = res_to_be_moved["chain_index"]
    
    def _free_position(self, res_to_free):
        self.lattice[res_to_free["x"], res_to_free["y"]] = 0


    def _are_diagonal(self, res1, res2):
        return (abs(res1['x'] - res2['x']) == 1 and abs(res1['y'] - res2['y']) == 1)

    def _find_empty_diagonal(self, i, j):
        possible_pos = []
        if self.lattice[i-1,j-1] == 0:
            possible_pos.append((i-1, j-1))
        if self.lattice[i+1,j+1] == 0:
            possible_pos.append((i+1, j+1))
        if self.lattice[i-1,j+1] == 0:
            possible_pos.append((i-1, j+1))
        if self.lattice[i+1,j-1] == 0:
            possible_pos.append((i+1, j-1))
        return possible_pos
    
    def _find_L_positions(self, res_i, res_i1):
        diag_pos = self._find_empty_diagonal(res_i["x"], res_i["y"])
        neighbors_pos = self._find_empty_neighbors(res_i1["x"], res_i1["y"])
        intersection = list(set(diag_pos) & set(neighbors_pos))
        if len(intersection) > 0:
            # print(f"{res_i["chain_index"]} ({res_i["x"]}, {res_i["y"]} and {res_i1["chain_index"]} have a L position at {intersection}")
            # print(f"we keep {intersection[0]}")
            return intersection[0], True #on garde le premier pour l'instant
        else:
            
            return (None,None), False 

    def _find_C_positions(self, L, res_i, res_iminus1):
        neighbors_pos_L = self._find_empty_neighbors(L[0], L[1])
        neighbors_pos_i = self._find_empty_neighbors(res_i["x"], res_i["y"])
        # print(f'neighbors_pos_L: {neighbors_pos_L}')
        # print(f'neighbors_pos_i: {neighbors_pos_i}')
        possible_pos = list(set(neighbors_pos_L) & set(neighbors_pos_i)) # find the common empty neighbors
        if len(possible_pos) > 0:
            # print(f"{res_i["chain_index"]} and {res_iminus1["chain_index"]} have a C empty position at {possible_pos}")
            return possible_pos[0], "empty", True # C is empty
        elif self._are_topological_neighbors({"x":L[0],"y":L[1]}, res_iminus1): # check if C is occupied by res_iminus1
            # print(f"there is res_iminus1 at C position {res_iminus1["chain_index"]} at coordinates {res_iminus1["x"], res_iminus1["y"]}")
            return [(res_iminus1["x"], res_iminus1["y"])], "C is res_iminus1", True
        else:
            # print(f"no possible C position for {res_i["chain_index"]}")
            return (None, None), "C is occupied", False

    def _translate_chain(self, x, y):
        # translate the lattice by x and y to recenter the conformation
        self.lattice = np.roll(self.lattice, x, axis = 0)
        self.lattice = np.roll(self.lattice, y, axis = 1)
        # update the coordinates of the amino acids in the sequence object
        for i in range(self.sequence.length):
            self.sequence.aa_coord_update(i, self.sequence.aa_coord[i]["x"] + x, self.sequence.aa_coord[i]["y"] + y)


    def _need_translation_x(self):
        list_x = [temp_dict["x"] for temp_dict in self.sequence.aa_coord]
        print(list_x)
        max_x = max(list_x)
        min_x = min(list_x)
        print(max_x)
        print(min_x)
        center = self.size // 2
        middle = (min_x + max_x) // 2
        dx = center - middle
        print("max_x >", self.size - 2 , "or min_x < 1:", max_x > self.size - 2 or min_x < 1)
        if max_x >= self.size - 2 or min_x < 1:
            return True, dx
        else:
            return False, 0

    def _need_translation_y(self):
        list_y = [temp_dict["y"] for temp_dict in self.sequence.aa_coord]
        print(list_y)
        max_y = max(list_y)
        min_y = min(list_y)
        print(max_y)
        print(min_y)
        center = self.size // 2
        middle = (min_y + max_y) // 2
        dy = center - middle
        print("max_y > ", self.size - 2 ,"or min_y < 1:", max_y > self.size - 2 or min_y < 1)
        if max_y >= self.size - 2 or min_y < 1:
            return True, dy
        else:
            return False, 0

    def translation_check_and_apply(self):
        check_x, dx = self._need_translation_x()
        check_y, dy = self._need_translation_y()
        if check_x or check_y:
            self._translate_chain(dx, dy)


































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