from moves.move_utils import find_empty_diagonal,find_empty_neighbors,are_topological_neighbors,move,free_position

def get_pull_aa_to_check(lattice, chain_index):

    res_i = lattice.sequence.aa_coord[chain_index - 1] # equals to CHAIN INDEX
    res_iplus1 = lattice.sequence.aa_coord[chain_index]
    res_iminus1 = lattice.sequence.aa_coord[chain_index - 2] # equals to CHAIN INDEX -1
    res_iminus2 = lattice.sequence.aa_coord[chain_index - 3] # equals to CHAIN INDEX -2
    return res_iminus2, res_iminus1, res_i, res_iplus1

def find_L_positions(lattice, res_i, res_i1):
    diag_pos = find_empty_diagonal(lattice,res_i["x"], res_i["y"])
    neighbors_pos = find_empty_neighbors(lattice, res_i1["x"], res_i1["y"])
    intersection = list(set(diag_pos) & set(neighbors_pos))
    if len(intersection) > 0:
        return intersection[0], True #on garde le premier pour l'instant
    else:
        return (None,None), False 

def find_C_positions(lattice, L, res_i, res_iminus1):
    neighbors_pos_L = find_empty_neighbors(lattice, L[0], L[1])
    neighbors_pos_i = find_empty_neighbors(lattice, res_i["x"], res_i["y"])
    possible_pos = list(set(neighbors_pos_L) & set(neighbors_pos_i)) # find the common empty neighbors
    if len(possible_pos) > 0:
        return possible_pos[0], "empty", True # C is empty
    elif are_topological_neighbors({"x":L[0],"y":L[1]}, res_iminus1): # check if C is occupied by res_iminus1
        return [(res_iminus1["x"], res_iminus1["y"])], "C is res_iminus1", True
    else:
        return (None, None), "C is occupied", False


def execute_pull_move(lattice, chain_index, res_i, res_iminus1, pos_L, pos_C):
    # move i and i+1 to L and C and free the prevous positions
    move(lattice,res_i, pos_L)
    move(lattice,res_iminus1, pos_C)

    free_position(lattice, res_i)
    free_position(lattice, res_iminus1)
    # update the coordinates of the amino acid in the sequence object
    lattice.sequence.aa_coord_update(chain_index - 1, pos_L[0], pos_L[1])
    lattice.sequence.aa_coord_update(chain_index - 2, pos_C[0], pos_C[1])

def propagate_pull(lattice, loop_index, empty_pos):

    while loop_index >= 2 :
        res_loop_iminus_1 = lattice.sequence.aa_coord[loop_index-1] # get res i-1
        res_loop_iminus_2= lattice.sequence.aa_coord[loop_index-2] # get res i-2
        is_neighbor = are_topological_neighbors(res_loop_iminus_1, res_loop_iminus_2)

        if is_neighbor is False:
            temp_pos = (res_loop_iminus_2["x"], res_loop_iminus_2["y"])
            free_position(lattice, res_loop_iminus_2)
            new_position = empty_pos.pop(0)
            move(lattice, lattice.sequence.aa_coord[loop_index-2], new_position)
            lattice.sequence.aa_coord_update(loop_index-2, new_position[0], new_position[1])
            empty_pos.append(temp_pos)
        if is_neighbor is True:
            break
        loop_index -= 1