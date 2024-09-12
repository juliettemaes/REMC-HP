# Libraries
from copy import deepcopy
import math
import numpy as np
import random as rd
import lattice as lat
import MC_search as mc
import sequence as seq
from variables import *

def temp_range(start, end, nb_replicat):
    step = (end - start) // nb_replicat
    return list(range(start, end, step))

def exchange_replicates(replicate1, replicate2):
    temp1 = replicate1.temperature
    temp2 = replicate2.temperature
    energy1 = replicate1.lattice.energy
    energy2 = replicate2.lattice.energy
    # compute delta
    delta = compute_delta(energy1, energy2, temp1, temp2)
    # accept or reject the exchange
    if evaluate_exchange(delta):
        # exchange the lattice
        temporary_lattice = deepcopy(replicate1.lattice)
        replicate1.lattice = deepcopy(replicate2.lattice)
        replicate2.lattice = deepcopy(temporary_lattice)
        # update the trajectory
        replicate1.trajectory.append(deepcopy(replicate1.lattice))
        replicate2.trajectory.append(deepcopy(replicate2.lattice))
        # print("exchange accepted")
        # # update the coordinates
        # for i in range(replicate1.lattice.sequence.length):
        #     print(np.where(replicate1.lattice.lattice == i+1))
        #     replicate1.lattice.sequence.aa_coord_update(
        #         i, 
        #         int(np.where(replicate1.lattice.lattice == i+1)[0][0]),
        #         int(np.where(replicate1.lattice.lattice == i+1)[1][0])
        #         )
        #     print("i", i)
        #     print(np.where(replicate2.lattice.lattice == i+1))
        #     replicate2.lattice.sequence.aa_coord_update(
        #         i, 
        #         int(np.where(replicate2.lattice.lattice == i+1)[0][0]),
        #         int(np.where(replicate2.lattice.lattice == i+1)[1][0])
        #         )
        # # update the energy
        # print("previous energy lattice i", energy1, "previous energy lattice j", energy2)
        # replicate1.lattice.energy = replicate1.lattice.calculate_energy()
        # replicate2.lattice.energy = replicate2.lattice.calculate_energy()
        # print("new energy lattice i", replicate1.lattice.energy, "new energy lattice j", replicate2.lattice.energy)


def compute_delta(energy_i, energy_j, temp_i, temp_j):
    beta_i = 1 / temp_i
    beta_j = 1 / temp_j
    return (beta_j - beta_i) * (energy_i - energy_j)

def evaluate_exchange(delta):
    # print("delta", delta)
    if delta < 0:
        return True
    else:
        probability = rd.random()
        return probability > math.exp(-delta)



def REMC_search(sequence, Tmin, Tmax, nb_replica, energy_optimal, max_iteration = 500, probability = 0.5):
    temperatures = temp_range(Tmin, Tmax, nb_replica)
    energy_best = 0
    offset = 0
    # initialize the replicates
    replicates = []
    for temperature in temperatures:
        lattice =  lat.Lattice(sequence = sequence)
        mc_search = mc.mc_search(lattice = lattice, temperature = temperature, probability = probability, max_iteration = max_iteration)
        replicates.append(mc_search)
    # print(replicates)
    # print(energy_best, "energy best")
    # print(energy_optimal, "energy optimal")
    # print(energy_best < energy_optimal)
    # run the replica exchange
    # print("STARTING REPLICA EXCHANGE")
    while energy_best > energy_optimal:
        for replica in range(nb_replica):
            replicates[replica].run()
            if replicates[replica].lattice.energy < energy_best:
                energy_best = replicates[replica].lattice.energy
            
        i = offset + 1
        while i + 1 <= nb_replica:
            j = i + 1
            # print("tentative exchange between", i, "and", j)
            exchange_replicates(replicates[i-1], replicates[j-1])
            i += 2
        offset = 1 - offset
    # return the conformation associated with the best energy
    for i, replicate in enumerate(replicates):
        if replicate.lattice.energy == energy_best:
            print("energy best", energy_best, "is associated with lattice", i)
            return replicate.lattice

    





        

# Main program
if __name__ == "__main__":

    sequence1 = seq.Sequence(SI_I)
    print(sequence1.hp_sequence)
    print(REMC_search(sequence1, 160, 220, 10, energy_optimal = -5, max_iteration = 10000, probability = 0.5))
    print(sequence1.hp_sequence)
    # print(lat.sequence.hp_sequence)
    # lat = lat.Lattice(sequence = sequence1)
    # energy_initial = lat.calculate_energy()
    # mc_search = mc.mc_search(lattice = lat, temperature = 160, probability = 0.5, max_iteration = 10000)
    # mc_search.run()
    # print("INITIAL ENERGY", energy_initial)
    # print("ENERGY LAST LATTICE", lat.energy)
    # print("FINAL ENERGY", mc_search.lattice.energy)