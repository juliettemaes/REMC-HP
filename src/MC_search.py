# Libraries
from copy import deepcopy
import numpy as np
import random as rd
import sequence as seq
import lattice as lat
import math
from variables import NB_ITER, TEMP, RHO

class MC_search():
    
    def __init__(self, lattice, temperature: int, probability: float, max_iteration: int):
        self.lattice = lattice
        self.trajectory = [lattice]
        self.temperature = temperature
        self.probability = probability
        self.max_iteration = max_iteration
        self.lattice.energy = self.lattice.calculate_energy2()
        self.run()
        

    def run(self):
        for i in range(self.max_iteration): # add the energy condition
            # choose a random amino acid
            aa = rd.randint(1, self.lattice.sequence.length)
            print("TRY A MOVE ON AA", aa)
            # make a move
            self.make_move(aa)
            print("LATTICE AFTER MOVE")
            print(self.lattice.lattice)
            # choose wether to accept the conformation or not
            if self.accept_conformation():
                # check if the conf should be translated
                self.lattice.translation_check_and_apply()
                self.trajectory.append(deepcopy(self.lattice)) # accept the conformation and add it to the trajectory
            else:
                self.lattice = deepcopy(self.trajectory[-1]) # reject the conformation


    def make_move(self, aa):
        proba = rd.random()
        if aa == 1 or aa == self.lattice.sequence.length:
            print("try end move")
            self.lattice.end_move(aa)
        elif proba < self.probability:
            # perform pull move
            print("try pull move")
            self.lattice.pull_move(aa)
        else: # perform VSHD move
            move = rd.choice(["corner_move", "cks_move"])
            if move == "corner_move":
                print("try corner move")
                self.lattice.corner_move(aa)
            else:
                print("try cks move")
                self.lattice.cks_move(aa)

    def accept_conformation(self):
        # calculate the energy of the new lattice
        energy_new_conf = self.lattice.calculate_energy2()
        print("ENERGY NEW CONF", energy_new_conf)
        print("ENERGY OLD CONF", self.lattice.energy)
        # accept or reject the new conformation
        # compare if a move has been made
        print("NEW LATTICE")
        print(self.lattice.lattice)
        print("LAST OF TRAJECTORY")
        print(self.trajectory[-1].lattice)
        if np.array_equal(self.lattice.lattice, self.trajectory[-1].lattice):
            print("NO MOVE MADE")
            return False
        else:
            print("MOVE MADE")
            if energy_new_conf <= self.lattice.energy:
                self.lattice.energy = energy_new_conf
                print("ACCEPTED CONFORMATION")
                return True
            else:
                if rd.random() <= math.exp((self.lattice.energy - energy_new_conf) / self.temperature):
                    self.lattice.energy = energy_new_conf
                    print("ACCEPTED CONFORMATION")
                    return True
                else:
                    print("REJECTED CONFORMATION")
                    return False

