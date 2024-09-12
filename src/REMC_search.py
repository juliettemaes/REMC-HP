from copy import deepcopy
import math
import numpy as np
import random as rd
import lattice as lat
import MC_search as mc
from typing import List, Optional, Union

def temp_range(start: int, end: int, nb_replicat: int) -> List[int]:
    """
    Generate a list of temperatures for replica exchange within a given range.

    Args:
        start (int): The min temperature.
        end (int): The max temperature.
        nb_replicat (int): The number of desired temperatures.
    Returns:
        List[int]: A list of temperatures.
    """
    step = (end - start) // nb_replicat
    return list(range(start, end, step))

def exchange_replicates(replicate1: mc.mc_search, replicate2: mc.mc_search) -> bool:
    """
    Attempt to exchange the states of two replica MC simulations based on the Metropolis like criterion.

    Args:
        replicate1 (mc.MCSearch): The first replicate to exchange.
        replicate2 (mc.MCSearch): The second replicate to exchange.

    Returns:
        bool: True if the exchange is accepted, False otherwise.
    """
    temp1 = replicate1.temperature
    temp2 = replicate2.temperature
    energy1 = replicate1.lattice.energy
    energy2 = replicate2.lattice.energy
    # Compute delta
    delta = compute_delta(energy1, energy2, temp1, temp2)
    # Accept or reject the exchange
    if evaluate_exchange(delta):
        # Exchange the lattice
        temporary_lattice = deepcopy(replicate1.lattice)
        replicate1.lattice = deepcopy(replicate2.lattice)
        replicate2.lattice = deepcopy(temporary_lattice)
        # Update the trajectory
        replicate1.trajectory.append(deepcopy(replicate1.lattice))
        replicate2.trajectory.append(deepcopy(replicate2.lattice))
        return True
    else:
        return False

def compute_delta(energy_i: float, energy_j: float, temp_i: float, temp_j: float) -> float:
    """
    Compute the energy difference (delta) for the exchange criterion.

    Args:
        energy_i (float): The energy of the first replica.
        energy_j (float): The energy of the second replica.
        temp_i (float): The temperature of the first replica.
        temp_j (float): The temperature of the second replica.

    Returns:
        float: The computed delta.
    """
    beta_i = 1 / temp_i
    beta_j = 1 / temp_j
    return (beta_j - beta_i) * (energy_i - energy_j)

def evaluate_exchange(delta: float) -> bool:
    """
    Evaluate whether to accept or reject the replicat exchange based on the computed delta.

    Args:
        delta (float): The computed delta from the energy difference.

    Returns:
        bool: True if the exchange is accepted, False otherwise.
    """
    if delta < 0:
        return True
    else:
        probability = rd.random()
        return probability > math.exp(-delta)

def REMC_search(
    sequence: str,
    Tmin: float,
    Tmax: float,
    nb_replica: int,
    energy_optimal: float,
    max_iteration: int = 500,
    probability: float = 0.5
) -> Optional[lat.Lattice]:
    """
    Perform a replica exchange Monte Carlo (REMC) search to find the lattice conformation
    with the lowest energy.

    Args:
        sequence (str): The amino acid sequence to be used in the lattice.
        Tmin (float): The minimum temperature for replica exchange range.
        Tmax (float): The maximum temperature for replica exchange range.
        nb_replica (int): The number of replicas (temperatures) to use.
        energy_optimal (float): The target energy to achieve.
        max_iteration (int, optional): The maximum number of iterations for each MC search. Defaults to 500.
        probability (float, optional): The probability of acceptance for each MC move. Defaults to 0.5.

    Returns:
        Optional[lat.Lattice]: The lattice with the best energy, or None if no lattice meets the criteria.
    """
    temperatures = temp_range(Tmin, Tmax, nb_replica)
    energy_best = float('inf')
    offset = 0
    # Initialize the replicates
    replicates: List[mc.mc_search] = []
    for temperature in temperatures:
        lattice = lat.Lattice(sequence=sequence)
        mc_search = mc.mc_search(
            lattice=lattice,
            temperature=temperature,
            probability=probability,
            max_iteration=max_iteration,
            target_energy=energy_optimal
        )
        replicates.append(mc_search)
    # Run the replica exchange
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
    # Return the conformation associated with the best energy
    for i, replicate in enumerate(replicates):
        if replicate.lattice.energy == energy_best:
            print("energy best", energy_best, "is associated with lattice", i)
            return replicate.lattice
