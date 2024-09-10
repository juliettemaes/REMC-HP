# Libraries
from copy import deepcopy
import numpy as np
import random as rd
import sequence as seq
import lattice as lat
import amino_acid



# Main program
if __name__ == "__main__":
    sequence = seq.Sequence("ARNDCEQGHILKMFGTKPSTWYV")
    sequence.HP_convert()
    print(type(sequence.hp_sequence))
    print(sequence.length)
    print(sequence.hp_sequence[1])
    print(sequence.hp_sequence)
    print(sequence.aa_coord)
    print()
    print(sequence.aa_coord[1]["type"])



    sequence = seq.Sequence("ARNDCEKMYV")
    print(sequence.HP_convert())
    lattice2 = lat.Lattice(sequence = sequence)
    #print(lattice2.size)
    # lattice2.initialize_random()
    #print(lattice2.lattice_initial)
    print(lattice2)
    print(lattice2.sequence.hp_sequence)
    #print(lattice2.size)

    for i, aa in enumerate(lattice2.sequence.hp_sequence):
        print(sequence.aa_coord[i])
    lattice2.calculate_energy()
    print(lattice2.energy)
    print(lattice2.calculate_energy2())
    #lattice2.end_move(1)
    # lattice2.end_move(10)
    # print(lattice2.lattice)
    # lattice2.end_move(1)
    print("lattice original")
    print(lattice2.lattice)
    for aa in range(2, lattice2.sequence.length-1):
        print(f"lattice after {aa} move")
        lattice2.cks_move(aa)
        print(lattice2.lattice)