# Library
import numpy as np

class Move():
    def __init__(self, lattice, residue):
        self.residue = residue
        self.lattice = lattice
        self.position = None

    def __str__(self):
        return f"{self.residue} at position {self.position}"

    def end_move(self):
        pass

    def corner_move(self):
        pass

    def cks_move(self):
        pass

    def pull_move(self):
        pass



def end_move(aa, x,y, lattice):
    # move the amino acid to the end of the chain
    if aa == 0 or aa == lattice.sequence.length-1: # the aa is at the end of the chain
        print(lattice.sequence.length)
        possible_pos = self._find_empty_neighbors(x,y)
        if len(possible_pos) != 0:
            pos = rd.choice(possible_pos)
        lattice.lattice[x,y] = 0
        x = pos[0]
        y = pos[1]
        lattice.lattice[x,y] = aa + 1

def corner_move(aa, x,y, lattice):
    pass


test = np.array([[1,2,3],[4,5,6],[7,8,9]])
print(test)
print(np.where(test == 4))
print(test[np.where(test == 4)])