# Libraries
from copy import deepcopy
import numpy as np
import random as rd
import sequence as seq
import lattice as lat
import math
from variables import NB_ITER, TEMP, RHO

class mc_search():
    
    def __init__(self, lattice, temperature: int, probability: float, max_iteration: int, target_energy: int = None):
        self.lattice = lattice
        self.trajectory = [lattice]
        self.temperature = temperature
        self.probability = probability
        self.max_iteration = max_iteration
        self.lattice.energy = self.lattice.calculate_energy()
        self.target_energy = target_energy
        

    def run(self):
        for i in range(self.max_iteration): # add the energy condition

            #print("Current energy : ",self.lattice.energy)
            #print("Target energy : ",self.target_energy)
            if self.target_energy is not None and self.lattice.energy == self.target_energy:
                break

            # choose a random amino acid
            aa = rd.randint(1, self.lattice.sequence.length)
            # make a move
            self.make_move(aa)
            # choose wether to accept the conformation or not
            if self.accept_conformation():
                # check if the conf should be translated
                self.lattice.translation_check_and_apply()
                self.trajectory.append(deepcopy(self.lattice)) # accept the conformation and add it to the trajectory
            else:
                self.lattice = deepcopy(self.trajectory[-1]) # reject the conformation


    def make_move(self, aa):
        proba = rd.random()
        # check if the amino acid is at the end of the sequence
        if aa == 1 or aa == self.lattice.sequence.length:
            self.lattice.end_move(aa)
        elif proba < self.probability:
            # perform pull move
            self.lattice.pull_move(aa)
        else: # perform VSHD move
            move = rd.choice(["corner_move", "cks_move"])
            if move == "corner_move":
                self.lattice.corner_move(aa)
            else:
                self.lattice.cks_move(aa)

    def accept_conformation(self):
        # calculate the energy of the new lattice
        energy_new_conf = self.lattice.calculate_energy()
        # accept or reject the new conformation
        # compare if a move has been made
        if np.array_equal(self.lattice.lattice, self.trajectory[-1].lattice):
            return False
        else:
            if energy_new_conf <= self.lattice.energy:
                self.lattice.energy = energy_new_conf
                return True
            else:
                rand_num = rd.random()
                if rand_num > math.exp((self.lattice.energy - energy_new_conf) / self.temperature):
                    self.lattice.energy = energy_new_conf
                    return True
                else:
                    return False

