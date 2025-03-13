import math
import random
from copy import deepcopy

# INIT
with open("P12956.fasta", "r") as f:
    lines = f.readlines()
    o_sequence = list("".join([l.strip() for l in lines[1:]]))
    f.close()

with open("deep_mutational_scan.tsv", "r") as f:
    lines = f.readlines()
    dms_labels = lines[0].strip().split("\t")
    dms = [l.strip().split("\t") for l in lines[1:]]
    f.close()

DEEPMUTDATALABELS = dms_labels
DEEPMUTDATA = dms
PROTEIN = o_sequence
N = 10**5

def calc_identity(mutated_sequence):

    diff = 0

    for i in range(len(mutated_sequence)):
        if mutated_sequence[i] != PROTEIN[i]:
            diff += 1

    return diff/len(mutated_sequence)


# Algorithme évolutif

mutated_sequences = []

for _ in range(500):

    mutated_sequence = deepcopy(PROTEIN)

    while(calc_identity(mutated_sequence) < 0.50):

        s = "NA"
        k = 0

        # 1. Generate a random mutation
        i = random.randint(0, len(PROTEIN)-1)
        while (s == "NA"):
            k = random.randint(0, 19)
            # Taking s from DEEPMUTDATA matrix instead of creating a vector for each amino acid of the protein
            s = DEEPMUTDATA[ DEEPMUTDATALABELS.index(mutated_sequence[i]) ][k]

        # 2. Calculate the probability of fixation
        s = float(s.replace(",", "."))
        pfix = (1 - math.e**(-2*s) ) / (1 - math.e**(-4*N*s))

        # 3. Determine if mutation is successful
        r = random.random()
        if r < pfix:
            mutated_sequence[i] = DEEPMUTDATALABELS[k]

    mutated_sequences.append("".join(mutated_sequence))


# output sequences
outstream = open("mutated_sequences.txt", "w")
outstream.write("\n".join(mutated_sequences))
outstream.close()