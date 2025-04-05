import math
import random
from copy import deepcopy

# INIT
with open("P12956.fasta", "r") as f:
    lines = f.readlines()
    o_sequence = list("".join([l.strip() for l in lines[1:]]))
    f.close()

with open("tableau_suplementaire_2.csv", "r") as f:
    lines = f.readlines()
    # dms_labels = lines[0].strip().split("\t")
    dms = [l.split(";") for l in lines[1:]]
    f.close()

DEEPMUTDATALABELS = {
    "W" : 0, "F" : 1, "Y" : 2, "P" : 3, "M" : 4,
    "I" : 5, "L" : 6, "V" : 7, "A" : 8, "G" : 9,
    "C" : 10,"S" : 11, "T" : 12, "Q" : 13, "N" : 14,
    "D" : 15,"E" : 16, "H" : 17, "R" : 18, "K" : 19
}

# DEEPMUTDATALABELS = dms_labels
PROTEIN = [DEEPMUTDATALABELS[aa] for aa in o_sequence[557:]]
DEEPMUTDATA = [["NA"] * 20] * len(PROTEIN)


for line in dms:
    if line[4].isdigit() and line[5] in DEEPMUTDATALABELS.keys():
        DEEPMUTDATA[int(line[4]) - 558][DEEPMUTDATALABELS[line[5]]] = line[14]

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
            s = DEEPMUTDATA[i][k]

            if s == "NA":
                continue

        # 2. Calculate the probability of fixation
        s = float(s.replace(",", "."))
        pfix = (1 - math.e**(-2*s) ) / (1 - math.e**(-4*N*s))

        # 3. Determine if mutation is successful
        r = random.random()
        if r < pfix:
            mutated_sequence[i] = k

    mutated_sequences.append("".join([list(DEEPMUTDATALABELS.keys())[aa_index] for aa_index in mutated_sequence]))


# output sequences
outstream = open("mutated_sequences.txt", "w")
outstream.write("\n".join(mutated_sequences))
outstream.close()