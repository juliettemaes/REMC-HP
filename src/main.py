# Libraries
import sequence as seq
import argparse
from variables import *
from REMC_search import *

# Main program
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='REMC search')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--hpsequence', type = str, help = 'HP sequence')
    group.add_argument('--aasequence', type = str, help = 'AA sequence')
    parser.add_argument('--optimal_energy', type = int, help = 'Target Energy', required=True)

    args = parser.parse_args()
    
    if args.aasequence:
        sequence = seq.Sequence(sequence = args.aasequence)
    else:
        sequence = seq.Sequence(hp_sequence = args.hpsequence)

    REMC_search(sequence, T_MIN, T_MAX, STEP, energy_optimal = args.optimal_energy, max_iteration = MAX_ITERATIONS, probability = PROBABILITY)