# Libraries
from copy import deepcopy
import numpy as np
import random as rd
import sequence as seq
import lattice as lat
import math
from variable import NB_ITER, TEMP, RHO

class MC_search(Lattice):
    def __init__(self, lattice):
        self.lattice = lattice

## NOT WORKING : NEED TO SOLVE HOW TO DEFINE THE CONFORMATION + separate the function which choose if we accept the confornmation or not
    def run(self, nb_iter = NB_ITER, temp = TEMP ):
        for i in range(nb_iter):
            new_lattice = deepcopy(self.lattice)
            # choose a random amino acid
            aa = rd.randint(0, self.sequence.length - 1) + 1
            # choose a move
            if rd.random() < RHO:
                # perform pull move
            else: # perform VSHD move
                

            # calculate the energy of the new lattice
            new_lattice.calculate_energy()

            # accept or reject the new conformation
            if new_lattice.energy < self.lattice.energy:
                self.lattice = deepcopy(new_lattice)
            else:
                if rd.random() > math.exp((self.lattice.energy - new_lattice.energy) / temp):
                    self.lattice = deepcopy(new_lattice)
        return self.lattice ## return or set an attribute? 

