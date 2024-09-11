# Libraries
from copy import deepcopy
import numpy as np
import random as rd
import sequence as seq
import lattice as lat
import MC_search as mc_search
import amino_acid



# Main program
if __name__ == "__main__":
    # sequence = seq.Sequence("ARNDCEQGHILKMFGTKPSTWYV")
    # sequence.HP_convert()
    # print(type(sequence.hp_sequence))
    # print(sequence.length)
    # print(sequence.hp_sequence[1])
    # print(sequence.hp_sequence)
    # print(sequence.aa_coord)
    # print()
    # print(sequence.aa_coord[1]["type"])



    sequence = seq.Sequence("ARNDCEKMYV")
    print(sequence.HP_convert())
    # lattice2 = lat.Lattice(sequence = sequence)
    # #print(lattice2.size)
    # # lattice2.initialize_random()
    # #print(lattice2.lattice_initial)
    # print(lattice2)
    # print(lattice2.sequence.hp_sequence)
    #print(lattice2.size)

    # for i, aa in enumerate(lattice2.sequence.hp_sequence):
    #     print(sequence.aa_coord[i])
    # lattice2.calculate_energy()
    # print(lattice2.energy)
    # print(lattice2.calculate_energy2())
    #lattice2.end_move(1)
    # lattice2.end_move(10)
    # print(lattice2.lattice)
    # lattice2.end_move(1)
    # print("lattice original")
    # print(lattice2.lattice)
    # for aa in range(2, lattice2.sequence.length):
    #     print(f"lattice after {aa} move")
    #     lattice2.pull_move(aa)


    lat = lat.Lattice(sequence = sequence)
    print(lat.lattice)
    energy_initial = lat.calculate_energy2()
    mc_search = mc_search.MC_search(lattice = lat, temperature = 160, probability = 0.5, max_iteration = 5000)
    print("INITIAL LATTICE")
    print("ENERGY INITIAL", energy_initial)
    print(mc_search.lattice.lattice_initial)
    print("FINAL LATTICE")
    print(mc_search.lattice.lattice)
    print("FINAL ENERGY", mc_search.lattice.energy)
    print(len(mc_search.trajectory))
    




    
    # test = np.array([[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  8,  7,  6,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 10,  9,  0,  5,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  2,  3,  4,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,],
    # [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,]])

    # lattice_test = lat.Lattice(sequence = sequence)
    # lattice_test.lattice = np.transpose(test)
    # for i in range(len(lattice_test.sequence.hp_sequence)):
    #     lattice_test.sequence.aa_coord_update(
    #         i, 
    #         int(np.where(lattice_test.lattice == i+1)[0][0]), 
    #         int(np.where(lattice_test.lattice == i+1)[1][0])
    #         )

    # print(lattice_test.lattice)
    # for aa in lattice_test.sequence.aa_coord:
    #     print(aa)

    # print(lattice_test._need_translation_x())
    # print(lattice_test._need_translation_y())
    # print("BEFORE TRANSLATION")
    # print(lattice_test.lattice)
    # lattice_test.translation_check_and_apply()
    # print("AFTER TRANSLATION")
    # print(lattice_test.lattice)