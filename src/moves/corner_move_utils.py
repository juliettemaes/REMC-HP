
def find_potential_corner(lattice, chain_index):
    if lattice.sequence.aa_coord[chain_index - 1]["x"] == lattice.sequence.aa_coord[chain_index - 2]["x"]:
        y_new = lattice.sequence.aa_coord[chain_index - 2]["y"]
        x_new = lattice.sequence.aa_coord[chain_index ]["x"]
    else:
        x_new = lattice.sequence.aa_coord[chain_index - 2]["x"]
        y_new = lattice.sequence.aa_coord[chain_index]["y"]
    return x_new, y_new

def get_corner_aa_to_check(lattice, chain_index):
    res_i = lattice.sequence.aa_coord[chain_index - 1]
    res_iminus1 = lattice.sequence.aa_coord[chain_index - 2]
    res_iplus1 = lattice.sequence.aa_coord[chain_index]
    return res_i, res_iminus1, res_iplus1

def execute_corner_move(lattice, chain_index, res_i, x_new, y_new):
    lattice.lattice[x_new, y_new] = chain_index
    lattice.lattice[res_i["x"], res_i["y"]] = 0
    lattice.sequence.aa_coord_update(chain_index - 1, x_new, y_new) # update the coordinates of the amino acid in the sequence object