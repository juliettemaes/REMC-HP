from moves.move_utils import relative_position,check_occupancy

def get_cks_aa_to_check(lattice, chain_index):

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

def get_alternative_positions(lattice, chain_index):
    res_prev2 = lattice.sequence.aa_coord[chain_index - 3]
    res_iminus1 = lattice.sequence.aa_coord[chain_index - 2]
    res_i = lattice.sequence.aa_coord[chain_index - 1]
    res_iplus1 = lattice.sequence.aa_coord[chain_index]
    return res_prev2, res_iminus1, res_i, res_iplus1


def are_U(res1, res2, res3, res4):
    cdt1 = abs(res1['x'] - res2['x']) + abs(res1['y'] - res2['y']) == 1
    cdt2 = abs(res2['x'] - res3['x']) + abs(res2['y'] - res3['y']) == 1
    cdt3 = abs(res3['x'] - res4['x']) + abs(res3['y'] - res4['y']) == 1
    cdt4 = abs(res4['x'] - res1['x']) + abs(res4['y'] - res1['y']) == 1
    return (cdt1 and cdt2 and cdt3 and cdt4)

def find_potential_U(res1, res2, res3, res4):
    x_new1 = res1['x'] - relative_position(res2, res1)[0]
    y_new1 = res1['y'] - relative_position(res2, res1)[1]

    x_new2 = res4['x'] - relative_position(res3, res4)[0]
    y_new2 = res4['y'] - relative_position(res3, res4)[1]
    
    return x_new1, y_new1, x_new2, y_new2


def validate_U(lattice, res_iminus1, res_i, res_iplus1, res_next2, chain_index):
    # check if the amino acids form a U
    if are_U(res_iminus1, res_i, res_iplus1, res_next2) and chain_index < lattice.sequence.length - 1:
        # find the potential positions of the amino acids after the move
        x_new1, y_new1, x_new2, y_new2 = find_potential_U(res_iminus1, res_i, res_iplus1, res_next2)
        if check_occupancy(lattice, x_new1, y_new1) and check_occupancy(lattice, x_new2, y_new2):
            return True, x_new1, y_new1, x_new2, y_new2
        else:
            return False, x_new1, y_new1, x_new2, y_new2
    return False, None, None, None, None

def execute_u_move(lattice, chain_index, res_i, res_iplus1, x_new1, y_new1, x_new2, y_new2):
    lattice.lattice[x_new1, y_new1] = chain_index
    lattice.lattice[x_new2, y_new2] = chain_index + 1
    # free the positions of the amino acids
    lattice.lattice[res_i["x"], res_i["y"]] = 0
    lattice.lattice[res_iplus1["x"], res_iplus1["y"]] = 0
    # update the coordinates of the amino acid in the sequence object
    lattice.sequence.aa_coord_update(chain_index - 1, x_new1, y_new1)
    lattice.sequence.aa_coord_update(chain_index, x_new2, y_new2)

def execute_alternative_u_move(lattice, chain_index, res_i, res_iminus1, x_new1, y_new1, x_new2, y_new2):
    # move the amino acids to the new positions
    lattice.lattice[x_new1, y_new1] = chain_index - 1
    lattice.lattice[x_new2, y_new2] = chain_index
    # free the positions of the amino acids
    lattice.lattice[res_iminus1["x"], res_iminus1["y"]] = 0
    lattice.lattice[res_i["x"], res_i["y"]] = 0
    # update the coordinates of the amino acid in the sequence object
    lattice.sequence.aa_coord_update(chain_index - 2, x_new1, y_new1)
    lattice.sequence.aa_coord_update(chain_index - 1, x_new2, y_new2)