from copy import deepcopy
import math
import numpy as np
import random as rd
import lattice as lat
import MC_search as mc

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
        return True
    else:
        return False


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
        mc_search = mc.mc_search(lattice = lattice, temperature = temperature, probability = probability, max_iteration = max_iteration, target_energy = energy_optimal)
        replicates.append(mc_search)
    # run the replica exchange
    while energy_best > energy_optimal:
        for replica in range(nb_replica):
            replicates[replica].run()
            if replicates[replica].lattice.energy < energy_best:
                energy_best = replicates[replica].lattice.energy
            print("replica", replica, "energy", replicates[replica].lattice.energy)
            print("energy best", energy_best, "target energy", energy_optimal)
            
        i = offset + 1
        while i + 1 <= nb_replica:
            j = i + 1
            test = exchange_replicates(replicates[i-1], replicates[j-1])
            if test:
                print("exchange between", i-1, "and", j-1, "successful")
            i += 2
        offset = 1 - offset
    # return the conformation associated with the best energy
    for i, replicate in enumerate(replicates):
        if replicate.lattice.energy == energy_best:
            print("energy best", energy_best, "is associated with lattice", i)
            return replicate.lattice
