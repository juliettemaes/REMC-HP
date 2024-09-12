# Libraries
from copy import deepcopy
import numpy as np
import random as rd
import sequence as seq
import lattice as lat
import math
from variables import NB_ITER, TEMP, RHO
from typing import Optional, Union

class mc_search:
    """
    Monte Carlo Search for protein 2D structure optimization.

    Attributes:
        lattice (lat.Lattice): The lattice containing the protein to be optimized.
        trajectory (List[lat.Lattice]): The list of lattice during the search.
        temperature (int): The temperature for the Monte Carlo simulation.
        probability (float): The probability of performing a pull move.
        max_iteration (int): The maximum number of iterations for the search.
        target_energy (Optional[int]): The target energy to reach (if specified).
    """

    def __init__(self, lattice: lat.Lattice, temperature: int, probability: float, max_iteration: int, target_energy: Optional[int] = None):
        """
        Initialize the Monte Carlo search object.

        Args:
            lattice (lat.Lattice): The lattice to be optimized.
            temperature (int): The temperature for the Monte Carlo simulation.
            probability (float): The probability of performing a pull move.
            max_iteration (int): The maximum number of iterations for the search.
            target_energy (Optional[int], optional): The target energy to reach. Defaults to None.
        """
        self.lattice = lattice
        self.trajectory = [lattice]
        self.temperature = temperature
        self.probability = probability
        self.max_iteration = max_iteration
        self.lattice.energy = self.lattice.calculate_energy()
        self.target_energy = target_energy

    def run(self) -> None:
        """
        Run the Monte Carlo search.

        This method performs the Monte Carlo simulation, making moves and accepting or rejecting conformations
        based on the acceptance criterion (Metropolis Criterion). 
        It updates the trajectory with accepted conformations.
        """
        for _ in range(self.max_iteration):
            # Check if the target energy has been reached
            if self.target_energy is not None and self.lattice.energy == self.target_energy:
                break

            # Choose a random amino acid
            aa = rd.randint(1, self.lattice.sequence.length)
            # Make a move
            self.make_move(aa)
            # Choose whether to accept the conformation or not
            if self.accept_conformation():
                # Check if the conformation should be translated
                self.lattice.translation_check_and_apply()
                self.trajectory.append(deepcopy(self.lattice))  # Accept the conformation and add it to the trajectory
            else:
                self.lattice = deepcopy(self.trajectory[-1])  # Reject the conformation

    def make_move(self, aa: int) -> None:
        """
        Attend to make a move based on the specified amino acid.

        Args:
            aa (int): The index of the amino acid to be moved.
        """
        proba = rd.random()
        # Check if the amino acid is at the end of the sequence
        if aa == 1 or aa == self.lattice.sequence.length:
            self.lattice.end_move(aa)
        elif proba < self.probability:
            # Perform pull move
            self.lattice.pull_move(aa)
        else:
            # Perform VSHD move
            move = rd.choice(["corner_move", "cks_move"])
            if move == "corner_move":
                self.lattice.corner_move(aa)
            else:
                self.lattice.cks_move(aa)

    def accept_conformation(self) -> bool:
        """
        Evaluate whether to accept or reject the new conformation based on the energy.

        Returns:
            bool: True if the new conformation is accepted, False otherwise.
        """
        # Calculate the energy of the new lattice
        energy_new_conf = self.lattice.calculate_energy()
        # Accept or reject the new conformation
        # Compare if a move has been made
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
