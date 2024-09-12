# Variables used in this project

# sequence variables
# Dictionary mapping 1-letter amino acid codes to HP model classification
AA_DICT = {
    'A': 'H',  # Alanine (Hydrophobic)
    'R': 'P',  # Arginine (Polar)
    'N': 'P',  # Asparagine (Polar)
    'D': 'P',  # Aspartic acid (Polar)
    'C': 'H',  # Cysteine (Hydrophobic)
    'Q': 'P',  # Glutamine (Polar)
    'E': 'P',  # Glutamic acid (Polar)
    'G': 'H',  # Glycine (Hydrophobic)
    'H': 'P',  # Histidine (Polar)
    'I': 'H',  # Isoleucine (Hydrophobic)
    'L': 'H',  # Leucine (Hydrophobic)
    'K': 'P',  # Lysine (Polar)
    'M': 'H',  # Methionine (Hydrophobic)
    'F': 'H',  # Phenylalanine (Hydrophobic)
    'P': 'P',  # Proline (Polar)
    'S': 'P',  # Serine (Polar)
    'T': 'P',  # Threonine (Polar)
    'W': 'H',  # Tryptophan (Hydrophobic)
    'Y': 'P',  # Tyrosine (Polar)
    'V': 'H',  # Valine (Hydrophobic)
}

# Sequence examples
seq_ex1 = "GRAIDGLGIVKPGYPGVWKPGVW" # would translate to (HP)2PH2PHP2HPH2P2HPH
SI_1 = "HPHPPHHPHPPHPHHPPHPH"
SI_2 = "HHPPHPPHPPHPPHPPHPPHPPHH"
SI_4 = "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP"
SI_3 = "PPHPPHHPPPPHHPPPPHHPPPPHH"

# Lattice variables
GRID_SIZE_FACTOR = 2

# Monte Carlo search variables
NB_ITER = 500
TEMP = 160  # to adjust
RHO = 0.5


# Replicat exchange Monte Carlo search variables
T_MIN = 160
T_MAX = 220
